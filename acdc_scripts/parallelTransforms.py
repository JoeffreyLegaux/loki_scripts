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
from commonTransforms import InlineMemberCalls, RemoveComments, RemovePragmas, RemovePragmaRegions, RemoveEmptyConditionals, AddSuffixToCalls, RemoveLoops, RemoveUnusedVariables

from arpege_parameters import params

# import errors should have already been caught in acdc_pragmas.py
import logical_lst
from openacc_transform import scc_transform_routine , alloc_temp

from termcolor import colored

class MakeParallel(Transformation):
    def __init__(self):        
        self.block_symbol = DeferredTypeSymbol(name=params.block_dimension, parent = DeferredTypeSymbol(name=params.cpg_opts_variable))
        self.block_counter = Variable(name='IBL', type=SymbolAttributes(BasicType.INTEGER, kind=Variable(name='JPIM')))
        self.outlined_routines = []

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
        arguments_map = {}
        nproma_arrays_map = {}
        nproma_names_map = {}
        nproma_arrays_dimensions = {}
        arguments_dimensions = {}
        locals_dimensions = {}
        
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
                    new_var = Variable(     name=f'Y{("D" if isArgument else "L")}_{var.name}', 
                                            type=SymbolAttributes(  DerivedType(name=f'ARRAY_{type_string}{len(var.dimensions)+1}D'), 
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


        routine.variables = SubstituteExpressions(nproma_arrays_map).visit(routine.variables)

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
        array_lengths = set([len(arguments_dimensions[var])+1 for var in arguments_dimensions] )

        routine.spec.prepend(Import(module="ARRAY_MOD", 
                symbols=tuple([DeferredTypeSymbol(name=f'ARRAY_{varl}D') for varl in array_lengths]) ))

        for length in array_lengths:
            routine.spec.prepend(Import(module=f'UTIL_ARRAY_{length}D_MOD'))

        return ([nproma_names_map[var] for var in nproma_names_map], locals_dimensions, arguments_dimensions)


    def transform_parallel_regions(self, routine):

        unmodified_spec = routine.spec.clone()

        # Add block counter declaration
        routine.variables += (self.block_counter,)

        # Transform arrays of NPROMA size into ARRAY_XD types
        nproma_arrays, local_dimensions, arguments_dimensions=self.npromaToArray_nD(routine)


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
        for region in FindNodes(PragmaRegion).visit(routine.body):
            print("region : ", region.pragma)
            region_num = region_num + 1

            (target, scc, outline, directive) = self.decodeParallelPragma(region.pragma.content)
    

            # Directives treatments :
            # TARGET = HOST|DEVICE
            # OUTLINE = 0 | 1+   =>  if >0 content of the section is placed in an separate subroutine (mandatory for declaring stack pointers)
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
            syncTransformation = MakeSync(pointerType=target.lower(), sections = (sync_routine.body,), nproma_arrays=nproma_arrays )
            sync_routine.apply(syncTransformation)
            # We keep the maximum amount of local pointers needed to pass FieldAPI variables to Sync subroutines
            total_FAPI_pointers = max (total_FAPI_pointers, syncTransformation.total_FAPI_pointers)

            sync_routine.apply(RemoveComments())
            sync_routine.apply(RemovePragmas())
            sync_routine.apply(RemoveEmptyConditionals())


            # We insert the sync routine as a member routine of its caller
            # Therefore we do not need any variable passing or declaration as we manipulate only variables of the caller
            routine.contains.append(sync_routine.clone(body=sync_routine.body, args=None, spec=None, contains=None))

            call_sync = CallStatement(name = DeferredTypeSymbol(name=sync_name), arguments=(), scope=None )


            if directive == 'openmp': 
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
                inner_loop = None
    
            elif directive == 'openacc':
                loop_pragma = Pragma(   keyword="ACC", 
                                        content = "PARALLEL LOOP GANG FIRSTPRIVATE ("+ boundary_variable.name +") PRIVATE (IBL)"
                                    )

                # We compute the new indices
                    #   YDCPG_BNDS%KBL    = IBL
                    #   YDCPG_BNDS%KSTGLO = 1 + (IBL - 1) * YDCPG_BNDS%KLON
                    #   YDCPG_BNDS%KFDIA  = MIN (YDCPG_BNDS%KLON, YDCPG_BNDS%KGPCOMP - YDCPG_BNDS%KSTGLO + 1)
                loop_new_statements = (Assignment(lhs = Variable(name='KBL', parent = boundary_variable), rhs = Variable(name = 'IBL')),)
                
                lhs_string = '1 + (IBL - 1) * YDCPG_BNDS%KLON'
                loop_new_statements += (Assignment( lhs = Variable(name='KSTGLO', parent = boundary_variable), 
                                                    rhs = parse_fparser_expression(lhs_string, scope=routine)
                                                  )
                                        ,)

                lhs_string = 'MIN (YDCPG_BNDS%KLON, YDCPG_BNDS%KGPCOMP - YDCPG_BNDS%KSTGLO + 1)'
                loop_new_statements += (Assignment( lhs = Variable(name='KFDIA', parent = boundary_variable), 
                                                    rhs = parse_fparser_expression(lhs_string, scope=routine)
                                                  )
                                        ,)


                jlon_assign = Assignment( lhs = Variable(name='KIDIA', parent = Variable(name='YLCPG_BNDS')),
                            rhs = Variable(name='JLON'))
                                        

                # Inner loop on JLON
                inner_loop = Loop(  variable = Variable(name='JLON'), 
                                    bounds = LoopRange( (Variable(name='KIDIA', parent = boundary_variable),
                                                        Variable(name='KFDIA', parent = boundary_variable))),

                                    body = (jlon_assign, new_region.body,)
                                )



            new_region = AddSuffixToCalls( suffix=('_SINGLE_COLUMN' if scc else '') + '_FIELD_API_' + target, node=region).transform_subroutine(routine = routine)

            new_loop = Loop(variable=self.block_counter, 
                                bounds=LoopRange(( IntLiteral(1), self.block_symbol )), 
                                body = (loop_new_statements, inner_loop if inner_loop else new_region.body,), 
                                pragma = (loop_pragma,)   
                            )

            if outline:
                new_subroutine = routine.clone( name = routine.name + "_PARALLEL_" + str(region_num), 
                                                body = new_region.body,
                                                # SCC requires NPROMA arrays declaration, temporary use original spec
                                                spec = unmodified_spec if scc else routine.spec.clone(), 
                                                contains = None)


                if scc:
                    true_symbols, false_symbols=logical_lst.symbols()
                    false_symbols.append('LHOOK')
    
                    scc_transform_routine(new_subroutine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols)
                    # self.npromaToArray_nD(new_subroutine)
                    new_subroutine.spec = routine.spec.clone()

                new_subroutine.apply(FieldAPIPtr(pointerType='host'))
                new_subroutine.apply(RemoveComments())
  
                remove_transform = RemoveUnusedVariables()
                new_subroutine.apply(remove_transform)

                self.outlined_routines.append(new_subroutine)

                call_parallel = CallStatement(name=new_subroutine.procedure_symbol, arguments=remove_transform.used_symbols)
                
                new_loop = Loop(variable=self.block_counter, 
                                bounds=LoopRange(( IntLiteral(1), self.block_symbol )), 
                                body = (loop_new_statements, call_parallel,), 
                                pragma = (loop_pragma,)   
                            )

                routine.body = Transformer({new_region:(call_sync,new_loop,)}).visit(routine.body)
            else:
                new_loop = Loop(variable=self.block_counter, 
                                bounds=LoopRange(( IntLiteral(1), self.block_symbol )), 
                                body = (loop_new_statements, new_region.body,), 
                                pragma = (loop_pragma,)   
                            )
                routine.body = Transformer({new_region:(call_sync,new_loop,)}).visit(routine.body)

                # Default behaviour : apply FieldAPI transformation to the loop body
                routine.apply(FieldAPIPtr(pointerType='host', node=new_loop ))

        routine.body = Transformer(regions_map).visit(routine.body)
        addFieldAPIPointers(routine, total_FAPI_pointers)

    def transform_subroutine(self, routine, **kwargs):

        
        self.transform_parallel_regions(routine)


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

