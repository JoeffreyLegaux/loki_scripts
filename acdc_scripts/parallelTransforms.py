from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType  )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

from loki.expression import FindTypedSymbols, FindVariables
from loki.expression.expr_visitors import SubstituteExpressions
from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList, LoopRange, IntLiteral
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from storable import retrieve

from syncTransforms import MakeSync, addFieldAPIPointers

from fieldAPITransforms import FieldAPIPtr
from commonTransforms import InlineMemberCalls, RemoveComments, RemovePragmas, RemovePragmaRegions, RemoveEmptyConditionals, \
                                AddSuffixToCalls, RemoveLoops, RemoveUnusedVariables, AddACCRoutineDirectives, RemoveUnusedImports

from arpege_parameters import params

# import errors should have already been caught in acdc_pragmas.py
import logical_lst
from openacc_transform import scc_transform_routine , alloc_temp

from termcolor import colored


class FindNodesOutsidePragmaRegion(FindNodes):
    """
     Find :any:`Node` instances that match a given criterion,
     but ignores nodes that are inside a PragmaRegion

    """
    def __init__(self, match, greedy=False):
        super().__init__(match, mode='type', greedy=greedy)


    def visit_tuple(self, o, **kwargs):
        """
        Visit all elements that are not PragmaRegion in the iterable and return the combined result.
        """
        ret = kwargs.pop('ret', self.default_retval())
        for i in o:
            if not isinstance(i, PragmaRegion):
                ret = self.visit(i, ret=ret, **kwargs)
        return ret or self.default_retval()

    def visit_Node(self, o, **kwargs):
        """
        Add the node to the returned list if it matches the criteria
        before visiting all children that are not PragmaRegion nodes.
        """

        ret = kwargs.pop('ret', self.default_retval())
        if self.rule(self.match, o):
            ret.append(o)
            if self.greedy:
                return ret
        for i in o.children:
            if not isinstance(i, PragmaRegion):
                ret = self.visit(i, ret=ret, **kwargs)
        return ret or self.default_retval()



class MakeParallel(Transformation):
    def __init__(self):        
        self.block_symbol = DeferredTypeSymbol(name=params.block_dimension, parent = DeferredTypeSymbol(name=params.cpg_opts_variable))
        self.block_counter = Variable(name='IBL', type=SymbolAttributes(BasicType.INTEGER, kind=Variable(name='JPIM')))
        self.outlined_routines = []
        self.has_column_loops = False

    def boundariesArgument(self, var):
        lower_bounds=[]
        upper_bounds=[]
        for dim in var.dimensions:
            if isinstance(dim, RangeIndex):
                lower_bounds.append(dim.lower if dim.lower else 1)
                upper_bounds.append(dim.upper)
            else:
                lower_bounds.append(1)
                upper_bounds.append(dim)

        lower_bounds.append(1)
        upper_bounds.append(self.block_symbol)

        has_lower = any(x != 1 for x in lower_bounds)

        bounds_kwargs = ( ('UBOUNDS', LiteralList(tuple(upper_bounds)) ), )
        if has_lower:
            bounds_kwargs += ( ('LBOUNDS', LiteralList(tuple(lower_bounds))), )
        return bounds_kwargs

    def getPrivateNameList(self, region):
        var_list=[]
        for assign in FindNodes(Assignment).visit(region.body):
            if isinstance(assign.lhs, Scalar):
                if (assign.lhs not in var_list):
                    var_list.append(assign.lhs)
        nameList = ""
        for var in var_list:
            nameList += ", " + var.name

        return nameList

    def decodeParallelPragma(self, pragma):
        #remove opening curly bracket
        trimmed_bracket = pragma.lower().split('{')[0]

        splitCommas = trimmed_bracket.split(',')
        clauses = {}
        for clause in splitCommas:
            if '=' in clause:
                clauses[clause.split('=')[0].strip()] = clause.split('=')[1].strip()

        target = None
        if (not clauses['target']):
            print("No target specified, defaulting to HOST")
            target = "HOST"
        elif (clauses['target'] == 'host' or clauses['target'] == 'device'):
            target = clauses['target'].upper()
        else :
            print(colored("Unkown target clause !!!!", "red"))
            exit(1)

        scc = False
        if ('vector' in clauses):
            if (clauses['vector'] == 'singlecolumn'):
                scc = True
            else:
                print(colored("Unkown vector clause !!!!", "red"))
                exit(1)


        outline = False
        if ('outline' in clauses):
            if (int(clauses['outline']) > 0):
                outline = True
            else:
                print(colored("Unkown outline clause !!!!", "red"))
                exit(1)


        directive = None
        if ('directive' in clauses):
            if (clauses['directive'] == 'openmp' or clauses['directive'] == 'openacc' ):
                directive = clauses['directive']
            else:
                print(colored("Unkown directive clause !!!!", "red"))
                exit(1)


        return(target, scc, outline, directive)





    def npromaToArray_nD(self, routine):
        """
        Search nproma-sized arrays (from params.nproma_aliases),
        transforms them into relevant ARRAY type,
        and add calls to its INIT/COPY/FINAL/WIPE methods.
        """
        arguments_map = {}
        nproma_arrays_map = {}
        nproma_names_map = {}
        nproma_arrays_dimensions = {}
        locals_dimensions = {}
        
        arrays_types_dimensions = set()

        init_copy_calls = ()
        wipe_final_calls = ()


        for var in FindVariables().visit(routine.spec):

            if isinstance(var, Array):

                to_transform = False

                firstdim = var.dimensions[0] 

                if isinstance(firstdim, RangeIndex):
                    firstdim = firstdim.upper

                # if firstdim in self.nproma_aliases:
                if firstdim in params.nproma_aliases:
                    isArgument = var in routine.arguments
                    type_string = "" 
                    if (var.type.dtype == BasicType.INTEGER):
                        type_string = "INT"
                    elif (var.type.dtype == BasicType.LOGICAL):
                        type_string = "LOG"

                    array_type_dim = f'{type_string}{len(var.dimensions)+1}'

                    arrays_types_dimensions.add(array_type_dim)
                        
                    new_var = Variable(     name=f'Y{("D" if isArgument else "L")}_{var.name}', 
                                            type=SymbolAttributes(  DerivedType(name=f'ARRAY_{array_type_dim}D'), 
                                                                    intent=var.type.intent, 
                                                                    scope = routine
                                                                    ) 
                                    )

                    nproma_arrays_map[var] = new_var
                    nproma_names_map[var.name] = new_var
                    nproma_arrays_dimensions[var.name] = var.dimensions

                    if not isArgument :

                        init_copy_calls += (CallStatement(name = DeferredTypeSymbol(name='INIT', parent = new_var), 
                                                            arguments=(),
                                                            kwarguments = self.boundariesArgument(var) + (('PERSISTENT', LogicLiteral(True)), ),
                                                            scope=routine) ,)

                        init_copy_calls += (CallStatement(name = DeferredTypeSymbol(name='COPY'), arguments=(new_var,), scope=routine), )
                        
                        wipe_final_calls += (CallStatement(name = DeferredTypeSymbol(name='WIPE'), arguments=(new_var,), scope=routine), )
                        
                        wipe_final_calls += (CallStatement(name = DeferredTypeSymbol(name='FINAL', parent = new_var), 
                                                                arguments=(), scope=routine), )


        routine.spec = SubstituteExpressions(nproma_arrays_map).visit(routine.spec)
        
        # Updates routine signature with new arguments through the _dummies property
        dummies_list = list(routine._dummies)
        for index, dummy in enumerate(dummies_list):
            if dummy.upper() in nproma_names_map:
                dummies_list[index] = nproma_names_map[dummy.upper()].name
        
        routine._dummies = tuple(dummies_list)

        body_vars_map={}
        for var in FindVariables().visit(routine.body):
            if isinstance(var, Array):
                if var.name in nproma_names_map:
                    new_var = nproma_names_map[var.name]
                    body_vars_map[var] = new_var.clone(dimensions = var.dimensions)

        routine.body = SubstituteExpressions(body_vars_map).visit(routine.body)


        routine.body.prepend(init_copy_calls)
        routine.body.append(wipe_final_calls)

        # Create imports for module ARRAY_MOD and UTIL_ARRAY_xD_MOD
        routine.spec.prepend(Import(module="ARRAY_MOD", 
                symbols=tuple([DeferredTypeSymbol(name=f'ARRAY_{varl}D') for varl in arrays_types_dimensions]) ))

        for type_dim in arrays_types_dimensions:
            routine.spec.prepend(Import(module=f'UTIL_ARRAY_{type_dim}D_MOD'))

        # We will very likely need to declare FIELD_BASIC variables
        routine.spec.prepend(Import(module="FIELD_MODULE", symbols=(DeferredTypeSymbol(name='FIELD_BASIC'),)))

        return ([nproma_names_map[var] for var in nproma_names_map], locals_dimensions)

    def makeOpenMPLoop(self, routine, boundary_variable, region):

        loop_pragma = Pragma(   keyword="OMP", 
                                content = "PARALLEL DO FIRSTPRIVATE ("+ boundary_variable.name +") PRIVATE (IBL" + 
                                self.getPrivateNameList(region) + ") "
                            )

        #We prepare the call to update_view for the block index
        loop_new_statements = CallStatement(name = 
                                                DeferredTypeSymbol(name="UPDATE_VIEW", 
                                                                    parent = boundary_variable), arguments=(), 
                                                                    kwarguments=(('BLOCK_INDEX', self.block_counter),
                                                                    ),    
                                            scope=None
                                            )
    
        return (loop_pragma, loop_new_statements)


    def makeOpenACCLoop(self, routine, boundary_variable):

        loop_pragma = Pragma(   keyword="ACC", 
                                    content = "PARALLEL LOOP GANG FIRSTPRIVATE ("+ boundary_variable.name +") PRIVATE (IBL)"
                            )

        # We compute the new indices
            #   YDCPG_BNDS%KBL    = IBL
            #   YDCPG_BNDS%KSTGLO = 1 + (IBL - 1) * YDCPG_BNDS%KLON
            #   YDCPG_BNDS%KFDIA  = MIN (YDCPG_BNDS%KLON, YDCPG_BNDS%KGPCOMP - YDCPG_BNDS%KSTGLO + 1)

        loop_new_statements = (Assignment(lhs = Variable(name='KBL', parent = boundary_variable), rhs = Variable(name = 'IBL')),)
        
        rhs_string = '1 + (IBL - 1) * YDCPG_BNDS%KLON'

        loop_new_statements += (Assignment( lhs = Variable(name='KSTGLO', parent = boundary_variable), 
                                            rhs = parse_fparser_expression(rhs_string, scope=routine)
                                          )
                                ,)

        rhs_string = 'MIN (YDCPG_BNDS%KLON, YDCPG_BNDS%KGPCOMP - YDCPG_BNDS%KSTGLO + 1)'

        loop_new_statements += (Assignment( lhs = Variable(name='KFDIA', parent = boundary_variable), 
                                            rhs = parse_fparser_expression(rhs_string, scope=routine)
                                          )
                                 ,)

        return (loop_pragma, loop_new_statements)

    def addJLONDefinition(self, routine):
        # jlon_found = False
        if "JLON" not in [v.name for v in routine.variables]:
            routine.variables += (Variable( name='JLON', 
                                            type=SymbolAttributes(BasicType.INTEGER,  kind=Variable(name='JPIM')))
                                ,)


    def makeOpenACCColumnsLoop(self, call_parallel, boundary_variable, local_boundary_variable, routine):
        # Inner loop on JLON

        # YLCPG_BNDS = YDCPG_BNDS
        assignments = (Assignment(lhs = local_boundary_variable, rhs = boundary_variable), ) 

        # YLCPG_BNDS%KIDIA = JLON
        assignments += (Assignment( lhs = Variable(name='KIDIA', parent = local_boundary_variable), 
                                    rhs = Variable(name='JLON') ), )
        # YLCPG_BNDS%KFDIA = JLON
        assignments += (Assignment( lhs = Variable(name='KFDIA', parent = local_boundary_variable), 
                                    rhs = Variable(name='JLON') ), )
        # Compute local stack variables boundaries
        assignments += (Assignment( lhs = Variable(name='YLSTACK'), rhs = Variable(name='YDSTACK')), )

        # YLSTACK%L = LOC (YDSTACK%F_P%DEVPTR (1,YDCPG_BNDS%KBL))
        L_string = 'LOC (YDSTACK%F_P%DEVPTR (1,IBL))'
        assignments += (Assignment( lhs = Variable(name='L', parent = Variable(name='YLSTACK')),
                                    rhs = parse_fparser_expression(L_string, scope=routine)),)
        # YLSTACK%U = YLSTACK%L + KIND (YDSTACK%F_P%DEVPTR) * SIZE (YDSTACK%F_P%DEVPTR (:,YDCPG_BNDS%KBL))
        U_string = 'YLSTACK%L + KIND (YDSTACK%F_P%DEVPTR) * SIZE (YDSTACK%F_P%DEVPTR (:,YDCPG_BNDS%KBL))'
        assignments += (Assignment( lhs = Variable(name='U', parent = Variable(name='YLSTACK')),
                                    rhs = parse_fparser_expression(U_string, scope=routine)  ),)                                                   

        columns_loop = Loop(    variable = Variable(name='JLON'), 
                                bounds = LoopRange( (Variable(name='KIDIA', parent = boundary_variable),
                                                    Variable(name='KFDIA', parent = boundary_variable))),
                                body = (assignments, call_parallel,),
                                pragma = (Pragma(keyword="ACC", content = f'LOOP VECTOR PRIVATE({local_boundary_variable},  YLSTACK)'),)

                            )


        if not self.has_column_loops:
            self.addJLONDefinition(routine)
            self.has_column_loops = True

        return columns_loop

    def transform_parallel_regions(self, routine):

        unmodified_spec = routine.spec.clone()

        # Add block counter declaration
        routine.variables += (self.block_counter,)

        # Use custom visitor to add PARALLEL suffix and stack variable to CallStatements outside PragmaRegions
        routine.apply(AddSuffixToCalls(suffix='_PARALLEL', custom_visitor = FindNodesOutsidePragmaRegion, additional_variables = ['YDSTACK']) )


        # Transform arrays of NPROMA size into ARRAY_XD types
        nproma_arrays, local_dimensions=self.npromaToArray_nD(routine)

        # Search the local variable that contains the boundaries
        boundary_variable = None
        for var in routine.variables :
            if (isinstance(var.type.dtype, DerivedType)):
                #                              'CPG_BNDS_TYPE'
                if (var.type.dtype.name == params.boundaries_type):
                    boundary_variable = var

        regions_map = {}
        region_num = -1

        total_FAPI_pointers = 0

        # Variable that tracks the need to add a secondary CPG_BNDS_TYPE variable 
        # for inner parallel loop in openacc
        local_boundary_variable=boundary_variable.clone(
                                                name='YL'+ boundary_variable.name[2:],
                                                type=boundary_variable.type.clone(intent=None))

        local_stack = Variable(name="YLSTACK", type=SymbolAttributes(DerivedType(name="STACK"), scope = routine))

        for region in FindNodes(PragmaRegion).visit(routine.body):
            print("region : ", region.pragma)
            region_num = region_num + 1

            (target, scc, outline, directive) = self.decodeParallelPragma(region.pragma.content)
    

            # Directives treatments :
            # TARGET = HOST|DEVICE
            # OUTLINE = 0 | 1+   =>  if >0 content of the section is placed in an separate subroutine 
            #                             (mandatory for declaring stack pointers)
            # VECTOR = SINGLECOLUMN  =>   SCC transformation
            # DIRECTIVE = OPENMP|OPENACC   
            # 
            # TARGET=HOST only => OpenMP directive, no outlining 
            # TARGET=DEVICE, VECTOR=SINGLECOLUMN, OUTLINE=1  => default GPU implementation  (DIRECTIVE=OPENACC implied)
            # OpenMP compatible with SINGLECOLUMN and OUTLINE, but not DEVICE !
            # OpenACC only compatible with DEVICE + SINGLECOLUMN + OUTLINE
            
            #Default values
            if not directive :
                directive = 'openmp'

            if directive == 'openacc' and not (target == 'DEVICE' and scc and outline):
                print(colored(f'Incompatible acdc clauses in pragma {region.pragma.content}', 'red'))
                exit(1)

            if directive == 'openmp' and target == 'DEVICE':
                print(colored(f'Incompatible acdc clauses (openmp and device) in pragma {region.pragma.content}', 'red'))
                exit(1)


            # Create synchronisation call 
            # building a subroutine that contains only the current region
            # Apply the Makesync transformation, then append it as a member routine
            sync_name = routine.name+"_PARALLEL_" + str(region_num) + "_SYNC_" + target
            sync_routine = routine.clone(body = region.body, spec = unmodified_spec, name = sync_name)

            sync_routine.apply(RemoveLoops())
            sync_routine.apply(AddSuffixToCalls(suffix='_SYNC_' + target))
            sync_routine.apply(RemoveUnusedImports())

            syncTransformation = MakeSync(pointerType=target.lower(), sections = (sync_routine.body,), nproma_arrays=nproma_arrays )
            sync_routine.apply(syncTransformation)
            # We keep the maximum amount of local pointers needed to pass FieldAPI variables to Sync subroutines
            total_FAPI_pointers = max (total_FAPI_pointers, syncTransformation.total_FAPI_pointers)

            sync_routine.apply(RemoveComments())
            sync_routine.apply(RemovePragmas())
            sync_routine.apply(RemoveEmptyConditionals())


            # We insert the sync routine as a member routine of its caller
            # Therefore we do not need any variable passing or declaration as we manipulate only variables of the caller
            if not routine.contains:
                routine.contains = Section(Intrinsic(text="CONTAINS"))

            includes = Section(tuple([n for n in FindNodes(Import).visit(sync_routine.spec) if n.c_import]))
            
            routine.contains.append(sync_routine.clone(body=sync_routine.body, args=None, spec=includes, contains=None))


            call_sync = CallStatement(name = DeferredTypeSymbol(name=sync_name), arguments=(), scope=None )


            if directive == 'openmp': 
                (loop_pragma, loop_new_statements) = self.makeOpenMPLoop(routine, boundary_variable, region)


            elif directive == 'openacc':
                (loop_pragma, loop_new_statements) = self.makeOpenACCLoop(routine, boundary_variable)



            if not outline:
                region = AddSuffixToCalls(  suffix=('_SINGLE_COLUMN' if scc else '') + '_FIELD_API_' + target, 
                                            node=region, 
                                            additional_variables=['YLSTACK'] if scc else []).transform_subroutine(routine = routine)

            if outline:

                new_subroutine = routine.clone( name = routine.name + "_PARALLEL_" + str(region_num), 
                                                body = region.body,
                                                spec = routine.spec.clone(),
                                                contains = None)

                new_subroutine.apply(AddSuffixToCalls(  suffix=('_SINGLE_COLUMN' if scc else '') + '_FIELD_API_' + target, 
                                            # node=region, 
                                            additional_variables=['YLSTACK'] if scc else []))

                new_subroutine.apply(RemoveUnusedImports())

                if scc:
                    true_symbols, false_symbols=logical_lst.symbols()
                    false_symbols.append('LHOOK')

                    scc_transform_routine(new_subroutine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols)

                new_subroutine.apply(FieldAPIPtr(pointerType=target.lower()))
                new_subroutine.apply(RemoveComments())
                new_subroutine.apply(RemovePragmas())
  
                # We do not want JLON and YSTACK in the arguments
                remove_transform = RemoveUnusedVariables(['JLON', 'YLSTACK'])
                new_subroutine.apply(remove_transform)

                new_subroutine.spec.append(Pragma(keyword='acc', content='routine seq'))
                if directive == 'openacc':
                    new_subroutine.apply(AddACCRoutineDirectives())

                self.outlined_routines.append(new_subroutine)


                call_arguments = remove_transform.used_symbols
                call_arguments = SubstituteExpressions({
                                        Scalar(name=boundary_variable.name):Scalar(name=local_boundary_variable.name),
                                        Scalar(name='YDSTACK'):Scalar(name='YLSTACK')
                                        }).visit(call_arguments)

                call_parallel = CallStatement(name=new_subroutine.procedure_symbol, arguments=call_arguments)



                if directive == 'openacc':
                    columns_loop = self.makeOpenACCColumnsLoop(call_parallel, boundary_variable, local_boundary_variable, routine)
                
            if not outline:
                loop_body = region.body
            else:
                if directive == 'openacc':
                    loop_body = columns_loop
                else:
                    loop_body = call_parallel

            blocks_loop = Loop( variable=self.block_counter, 
                                bounds=LoopRange(( IntLiteral(1), self.block_symbol )), 
                                body = (loop_new_statements, loop_body,), 
                                pragma = (loop_pragma,)   
                                )
            routine.body = Transformer({region:(call_sync,blocks_loop,)}).visit(routine.body)

            if not outline:
                # Default behaviour : apply FieldAPI transformation to the loop body
                routine.apply(FieldAPIPtr(pointerType=target.lower(), node=blocks_loop ))


        # We remove remaining imports only after having treated all Parallel blocks 
        routine.apply(RemoveUnusedImports())

        routine.variables += (local_boundary_variable, local_stack,)   

        routine.body = Transformer(regions_map).visit(routine.body)
        addFieldAPIPointers(routine, total_FAPI_pointers)

    def transform_subroutine(self, routine, **kwargs):

        
        self.transform_parallel_regions(routine)

        # In a routine with parallel blocks, we expect to have received a stack variable
        routine.spec.prepend(Import(module="STACK_MOD"))
        routine.arguments += (Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK"), intent='in')) ,)

        # get updatable fields and change their intent into INOUT
        FieldAPI_updatables = retrieve('../scripts/update_view.dat')

        #                      'CPG_BNDS_TYPE'
        FieldAPI_updatables[params.boundaries_type] = 1
        var_map = {}
        for var in routine.arguments :
            if (isinstance(var.type.dtype, DerivedType)):
                typename = var.type.dtype.name
                if typename in FieldAPI_updatables:
                    if FieldAPI_updatables[typename]:
                        new_type=var.type.clone(intent='inout')
                        var_map[var]= var.clone(type=new_type)

        routine.spec = SubstituteExpressions(var_map).visit(routine.spec)

