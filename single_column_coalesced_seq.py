# (C) Copyright 2018- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from more_itertools import split_at

from loki.expression import symbols as sym
from loki.transform import resolve_associates
from loki import ir
from loki import (
    Transformation, FindNodes, FindScopes, FindVariables,
    FindExpressions, Transformer, NestedTransformer,
    SubstituteExpressions, SymbolAttributes, BasicType, DerivedType,
    pragmas_attached, CaseInsensitiveDict, as_tuple, flatten,
    demote_variables, SubroutineItem, CallStatement
)


__all__ = ['SingleColumnCoalescedTransformationSeq']


def get_integer_variable(routine, name):
    """
    Find a local variable in the routine, or create an integer-typed one.

    Parameters
    ----------
    routine : :any:`Subroutine`
        The subroutine in which to find the variable
    name : string
        Name of the variable to find the in the routine.
    """
    if name in routine.variable_map:
        v_index = routine.variable_map[name]
    else:
        dtype = SymbolAttributes(BasicType.INTEGER)
        v_index = sym.Variable(name=name, type=dtype, scope=routine)
    return v_index


def kernel_remove_vector_loops(routine, horizontal):
    """
    Remove all vector loops over the specified dimension.

    Parameters
    ----------
    routine : :any:`Subroutine`
        The subroutine in the vector loops should be removed.
    horizontal : :any:`Dimension`
        The dimension specifying the horizontal vector dimension
    """
    loop_map = {}
    for loop in FindNodes(ir.Loop).visit(routine.body):
        if loop.variable == horizontal.index:
            loop_map[loop] = loop.body
    routine.body = Transformer(loop_map).visit(routine.body)


def kernel_get_locals_to_demote(routine, horizontal, arguments_to_demote, successors):

    argument_names = [v.name for v in routine.arguments]

    def _is_constant(d):
        """Establish if a given dimensions symbol is a compile-time constant"""
        if isinstance(d, sym.IntLiteral):
            return True

        if isinstance(d, sym.RangeIndex):
            if d.lower:
                return _is_constant(d.lower) and _is_constant(d.upper)
            return _is_constant(d.upper)

        if isinstance(d, sym.Scalar) and isinstance(d.initial , sym.IntLiteral):
            return True

        return False

    def _get_local_arrays(section):
        """
        Filters out local argument arrays that solely buffer the
        horizontal vector dimension
        """
        arrays = FindVariables(unique=False).visit(section)
        # Only demote local arrays with the horizontal as fast dimension
        arrays = [v for v in arrays if isinstance(v, sym.Array)]
        arrays = [v for v in arrays if v.name not in argument_names]
        arrays = [v for v in arrays if v.shape and v.shape[0] == horizontal.size]

        # Also demote arrays whose remaning dimensions are known constants
        arrays = [v for v in arrays if all(_is_constant(d) for d in v.shape[1:])]
        return arrays

    # Create a list of all local horizontal temporary arrays
    to_demote = _get_local_arrays(routine.spec)

    # Create a list of arguments variable that are included in arguments_to_demote 
    arg_arrays = FindVariables().visit(routine.spec)
    arg_arrays = [a for a in arg_arrays if a.name in arguments_to_demote]

    for array in arg_arrays :
        to_demote.append(array)

    # Successors = scheduler items that are called from the current routine
    successor_map={successor.routine.name: successor for successor in successors if isinstance(successor, SubroutineItem)} 
        
    for call in FindNodes(CallStatement).visit(routine.body) :
        if call.name.name in successor_map :

            # Build the list of demoted arrays that are passed to the subroutine call
            # The list contains the name of the arguments in the subroutine
            args_list=[]
            arg_map = {arg[1].name : arg[0] for arg in call.arg_iter()}
            for array in to_demote:
                if array.name in arg_map:
                    args_list.append(arg_map[array.name])

            # Check if the subroutine has not already been treated previously
            if (successor_map[call.name.name].trafo_data != {} and args_list != []): 
                # Check if new list is contained in existing
                for arg in args_list :
                    if arg.name not in successor_map[call.name.name].trafo_data :
                        # If we encounter new non-optional arguments, we are facing multiple demotion patterns.
                        if not arg.type.optional :
                            raise RuntimeError(f'Argument demotion incompatible with previous call. \
                                Routine {call.name.name}, argument {arg.name}')

                        # If the new argument is optional, we assume it was not present
                        # in the previous calls. This assumption might not be correct in
                        # some edge cases.
                        else :
                            print("appended optional")
                            successor_map[call.name.name].trafo_data.append(arg.name)

            else :
                # Otherwise, we put the demoted arguments list in the trafo_data member of the item
                successor_map[call.name.name].trafo_data=[a.name for a in args_list]

    return to_demote


def kernel_annotate_sequential_loops_openacc(routine, horizontal):
    """
    Insert ``!$acc loop seq`` annotations around all loops that
    are not horizontal vector loops.

    Parameters
    ----------
    routine : :any:`Subroutine`
        The subroutine in which to annotate sequential loops
    horizontal: :any:`Dimension`
        The dimension object specifying the horizontal vector dimension
    """
    with pragmas_attached(routine, ir.Loop):

        for loop in FindNodes(ir.Loop).visit(routine.body):
            # Skip loops explicitly marked with `!$loki/claw nodep`
            if loop.pragma and any('nodep' in p.content.lower() for p in as_tuple(loop.pragma)):
                continue

            if loop.variable != horizontal.index:
                # Perform pragma addition in place to avoid nested loop replacements
                loop._update(pragma=ir.Pragma(keyword='acc', content='loop seq'))



def resolve_masked_stmts(routine, loop_variable):
    """
    Resolve :any:`MaskedStatement` (WHERE statement) objects to an
    explicit combination of :any:`Loop` and :any:`Conditional` combination.

    Parameters
    ----------
    routine : :any:`Subroutine`
        The subroutine in which to resolve masked statements
    loop_variable : :any:`Scalar`
        The induction variable for the created loops.
    """
    mapper = {}
    for masked in FindNodes(ir.MaskedStatement).visit(routine.body):
        # TODO: Currently limited to simple, single-clause WHERE stmts
        assert len(masked.conditions) == 1 and len(masked.bodies) == 1
        ranges = [e for e in FindExpressions().visit(masked.conditions[0]) if isinstance(e, sym.RangeIndex)]
        exprmap = {r: loop_variable for r in ranges}
        assert len(ranges) > 0
        assert all(r == ranges[0] for r in ranges)
        bounds = sym.LoopRange((ranges[0].start, ranges[0].stop, ranges[0].step))
        cond = ir.Conditional(condition=masked.conditions[0], body=masked.bodies[0], else_body=masked.default)
        loop = ir.Loop(variable=loop_variable, bounds=bounds, body=cond)
        # Substitute the loop ranges with the loop index and add to mapper
        mapper[masked] = SubstituteExpressions(exprmap).visit(loop)

    routine.body = Transformer(mapper).visit(routine.body)


def resolve_vector_dimension(routine, loop_variable, bounds):
    """
    Resolve vector notation for a given dimension only. The dimension
    is defined by a loop variable and the bounds of the given range.

    TODO: Consolidate this with the internal
    `loki.transform.transform_array_indexing.resolve_vector_notation`.

    Parameters
    ----------
    routine : :any:`Subroutine`
        The subroutine in which to resolve vector notation usage.
    loop_variable : :any:`Scalar`
        The induction variable for the created loops.
    bounds : tuple of :any:`Scalar`
        Tuple defining the iteration space of the inserted loops.
    """
    bounds_str = f'{bounds[0]}:{bounds[1]}'

    bounds_v = (sym.Variable(name=bounds[0]), sym.Variable(name=bounds[1]))

    mapper = {}
    for stmt in FindNodes(ir.Assignment).visit(routine.body):
        ranges = [e for e in FindExpressions().visit(stmt)
                  if isinstance(e, sym.RangeIndex) and e == bounds_str]
        if ranges:
            exprmap = {r: loop_variable for r in ranges}
            loop = ir.Loop(variable=loop_variable, bounds=sym.LoopRange(bounds_v),
                           body=(SubstituteExpressions(exprmap).visit(stmt),) )
            mapper[stmt] = loop

    routine.body = Transformer(mapper).visit(routine.body)


class SingleColumnCoalescedTransformationSeq(Transformation):
    """
    Single Column Coalesced: Direct CPU-to-GPU transformation for
    block-indexed gridpoint routines.

    This transformation will remove individiual CPU-style
    vectorization loops from "kernel" routines and either either
    re-insert the vector loop at the highest possible level (without
    interfering with subroutine calls), or completely strip it and
    promote the index variable to the driver if
    ``hoist_column_arrays`` is set.

    Unlike the CLAW-targetting SCA extraction, this will leave the
    block-based array passing structure in place, but pass a
    thread-local array index into any "kernel" routines. The
    block-based argument passing should map well to coalesced memory
    accesses on GPUs.

    Note, this requires preprocessing with the
    :class:`DerivedTypeArgumentsTransformation`.

    Parameters
    ----------
    horizontal : :any:`Dimension`
        :any:`Dimension` object describing the variable conventions used in code
        to define the horizontal data dimension and iteration space.
    vertical : :any:`Dimension`
        :any:`Dimension` object describing the variable conventions used in code
        to define the vertical dimension, as needed to decide array privatization.
    block_dim : :any:`Dimension`
        Optional ``Dimension`` object to define the blocking dimension
        to use for hoisted column arrays if hoisting is enabled.
    directive : string or None
        Directives flavour to use for parallelism annotations; either
        ``'openacc'`` or ``None``.
    hoist_column_arrays : bool
        Flag to trigger the more aggressive "column array hoisting"
        optimization.
    """

    def __init__(self, horizontal, vertical=None, block_dim=None, directive=None,
                 demote_local_arrays=True, hoist_column_arrays=True):
        self.horizontal = horizontal
        self.vertical = vertical
        self.block_dim = block_dim

        assert directive in [None, 'openacc']
        self.directive = directive

        self.demote_local_arrays = demote_local_arrays
        self.hoist_column_arrays = hoist_column_arrays

    def transform_subroutine(self, routine, **kwargs):
        """
        Apply transformation to convert a :any:`Subroutine` to SCC format.

        Parameters
        ----------
        routine : :any:`Subroutine`
            Subroutine to apply this transformation to.
        role : string
            Role of the subroutine in the call tree; either
            ``"driver"`` or ``"kernel"``
        targets : list of strings
            Names of all kernel routines that are to be considered "active"
            in this call tree and should thus be processed accordingly.
        """
        role = kwargs['role']
        item = kwargs.get('item', None)
        targets = kwargs.get('targets', None)
        successors = kwargs.get('successors', ())

        if role == 'driver':
            self.process_driver(routine, targets=targets)

        if role == 'kernel':
            demote_locals = self.demote_local_arrays
            if item:
                demote_locals = item.config.get('demote_locals', self.demote_local_arrays)
                print(f'item {item.name}  transfo_data : {item.trafo_data} ')
            self.process_kernel(routine, demote_locals=demote_locals, arguments_to_demote=item.trafo_data, successors=successors)

    def process_kernel(self, routine, demote_locals=True, arguments_to_demote=[], successors=[]):
        """
        Applies the SCC loop layout transformation to a "kernel"
        subroutine. This will primarily strip the innermost vector
        loops and either re-insert the vector loop at the highest
        possible level (without interfering with subroutine calls),
        or completely strip it and promote the index variable to the
        driver if ``hoist_column_arrays`` is set.

        In both cases argument arrays are left fully dimensioned,
        allowing us to use them in recursive subroutine invocations.

        Parameters
        ----------
        routine : :any:`Subroutine`
            Subroutine to apply this transformation to.
        """

        pragmas = FindNodes(ir.Pragma).visit(routine.body)
        routine_pragmas = [p for p in pragmas if p.keyword.lower() in ['loki', 'acc']]
        routine_pragmas = [p for p in routine_pragmas if 'routine' in p.content.lower()]

        seq_pragmas = [r for r in routine_pragmas if 'seq' in r.content.lower()]
        if seq_pragmas:
            if self.directive == 'openacc':
                # Mark routine as acc seq
                mapper = {seq_pragmas[0]: ir.Pragma(keyword='acc', content='routine seq')}
                routine.body = Transformer(mapper).visit(routine.body)

            # Bail and leave sequential routines unchanged
            return

        vec_pragmas = [r for r in routine_pragmas if 'vector' in r.content.lower()]
        if vec_pragmas:
            if self.directive == 'openacc':
                # Bail routines that have already been marked and this processed
                # TODO: This is a hack until we can avoid redundant re-application
                return

        if self.horizontal.bounds[0] not in routine.variable_map:
            raise RuntimeError(f'No horizontal start variable found in {routine.name}')
        if self.horizontal.bounds[1] not in routine.variable_map:
            raise RuntimeError(f'No horizontal end variable found in {routine.name}')

        # raise RuntimeError(f'horizontal start and end variables found in {routine.name} !!!!')


        # Find the iteration index variable for the specified horizontal
        v_index = get_integer_variable(routine, name=self.horizontal.index)

        # Associates at the highest level, so they don't interfere
        # with the sections we need to do for detecting subroutine calls
        resolve_associates(routine)

        # Resolve WHERE clauses
        resolve_masked_stmts(routine, loop_variable=v_index)

        # Resolve vector notation, eg. VARIABLE(KIDIA:KFDIA)
        resolve_vector_dimension(routine, loop_variable=v_index, bounds=self.horizontal.bounds)

        # Remove all vector loops over the specified dimension
        kernel_remove_vector_loops(routine, self.horizontal)

        # Demote all private local variables having only horizontal dimension
        if demote_locals:
            to_demote = kernel_get_locals_to_demote(routine, self.horizontal, arguments_to_demote, successors)
            variables = tuple(v.name for v in to_demote)
            if variables:
                demote_variables(routine, variable_names=variables, dimensions=self.horizontal.size)

        # Add loop index variable
        if v_index not in routine.arguments:
            new_v = v_index.clone(type=v_index.type.clone(intent='in'))
            # Remove original variable first, since we need to update declaration
            routine.variables = as_tuple(v for v in routine.variables if v != v_index)
            routine.arguments += as_tuple(new_v)

        call_map = {}
        for call in FindNodes(ir.CallStatement).visit(routine.body):
            # Append new loop variable to call signature
            new_call = call.clone(arguments=call.arguments)
            new_call._update(kwarguments=new_call.kwarguments + ((self.horizontal.index, v_index),))
            call_map[call] = new_call
        routine.body = Transformer(call_map).visit(routine.body)


        # Mark all non-parallel loops as `!$acc loop seq`
        kernel_annotate_sequential_loops_openacc(routine, self.horizontal)

        # Mark routine as `!$acc routine seq` to make it device-callable
        routine.spec.append(ir.Pragma(keyword='acc', content='routine seq'))


    def process_driver(self, routine, targets=None):
        """
        Process the "driver" routine by inserting the other level
        parallel loops, and optionally hoisting temporary column
        arrays.

        Note that if ``hoist_column_arrays`` is set, the driver needs
        to be processed before any kernels are trnasformed. This is
        due to the use of an interprocedural analysis forward pass
        needed to collect the list of "column arrays".

        Parameters
        ----------
        routine : :any:`Subroutine`
            Subroutine to apply this transformation to.
        targets : list or string
            List of subroutines that are to be considered as part of
            the transformation call tree.
        """

        # Resolve associates, since the PGI compiler cannot deal with
        # implicit derived type component offload by calling device
        # routines.
        resolve_associates(routine)

        with pragmas_attached(routine, ir.Loop, attach_pragma_post=True):

            # add_horizontal_loop_to_kernel_call applies a transformation to the routine body
            # This messes up the first loop in case of multiple driver calls
            # Putting the calls in a list then call add_horizontal_... after the first loop solves this
            calls_to_hoist=[]
            for call in FindNodes(ir.CallStatement).visit(routine.body):
                if not call.name in targets:
                    continue

                # Find the driver loop by checking the call's heritage
                ancestors = flatten(FindScopes(call).visit(routine.body))
                loops = [a for a in ancestors if isinstance(a, ir.Loop)]
                if not loops:
                    # Skip if there are no driver loops
                    continue
                loop = loops[0]
                # Mark driver loop as "gang parallel".
                if self.directive == 'openacc':
                    if loop.pragma is None:
                        loop._update(pragma=(ir.Pragma(keyword='acc', content='parallel loop gang vector_length(32)'), ))
                        loop._update(pragma_post=(ir.Pragma(keyword='acc', content='end parallel loop'), ))
                calls_to_hoist.append(call)        
            
            # Apply hoisting of temporary "column arrays"
            for call in calls_to_hoist:
                self.add_horizontal_loop_to_kernel_call(routine, call)

    def add_horizontal_loop_to_kernel_call(self, routine, call):
        """
        Add an horizontal loop around kernel call at the driver level
        with an  ``!$acc loop vector directive``. 
        Also passes the loop variable as an additional positional
        argument in the kernel call.

        Parameters
        ----------
        routine : :any:`Subroutine`
            Subroutine to apply this transformation to.
        call : :any:`CallStatement`
            Call to subroutine from which we hoist the column arrays.
        """
        if call.not_active or call.routine is BasicType.DEFERRED:
            raise RuntimeError(
                '[Loki] SingleColumnCoalescedTransform: Target kernel is not attached '
                'to call in driver routine.'
            )

        if not self.block_dim:
            raise RuntimeError(
                '[Loki] SingleColumnCoalescedTransform: No blocking dimension found '
                'for column hoisting.'
            )

        kernel = call.routine
        call_map = {}

        # Find the iteration index variable for the specified horizontal
        v_index = get_integer_variable(routine, name=self.horizontal.index)
        if v_index.name not in routine.variable_map:
            routine.variables += as_tuple(v_index)

        # Append new loop variable to call signature
        new_call = call.clone(arguments=call.arguments)
        new_call._update(kwarguments=new_call.kwarguments + ((self.horizontal.index, v_index),))

        # Create a vector loop around the kernel invocation
        pragma = None
        if self.directive == 'openacc':
            pragma = ir.Pragma(keyword='acc', content='loop vector')

        v_start =kernel.variable_map[self.horizontal.bounds[0]]
        v_end = kernel.variable_map[self.horizontal.bounds[1]]
        bounds = sym.LoopRange((v_start, v_end))
        vector_loop = ir.Loop(variable=v_index, bounds=bounds, body=[new_call], pragma=(pragma,))
        call_map[call] = vector_loop

        routine.body = Transformer(call_map).visit(routine.body)
