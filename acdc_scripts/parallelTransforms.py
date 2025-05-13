from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType, FindTypedSymbols, FindVariables, SubstituteExpressions  )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

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
from openacc_transform import scc_transform_routine, alloc_temp

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
    def __init__(self,  FieldAPI_pointers={}):        
        self.FieldAPI_pointers = FieldAPI_pointers
        self.block_symbol = DeferredTypeSymbol( name=params.block_dimension, 
                                                parent = DeferredTypeSymbol(name=params.cpg_opts_variable))
        
        self.block_counter = Variable(  name=params.block_counter, 
                                        type=SymbolAttributes(BasicType.INTEGER, kind=Variable(name='JPIM')))
        
        self.block_min_symbol = DeferredTypeSymbol( name=self.block_counter.name+'MIN', 
                                                    parent = DeferredTypeSymbol(name=params.cpg_opts_variable))
        
        self.block_max_symbol = DeferredTypeSymbol( name=self.block_counter.name+'MAX',
                                                    parent = DeferredTypeSymbol(name=params.cpg_opts_variable))
        
        self.blocks_loop = Loop(variable=self.block_counter,
                                #bounds=LoopRange(( IntLiteral(1), self.block_symbol )),
                                bounds=LoopRange(( self.block_min_symbol, self.block_max_symbol )),
                                body = None,
                                pragma = None
                                )
         
        self.outlined_routines = []
        self.has_column_loops = False
        self.subroutines_to_transform = {}

    def boundariesArgument(self, var):
        lower_bounds=[]
        upper_bounds=[]
        for dim in var.dimensions:
            if isinstance(dim, RangeIndex):
                lower_bounds.append(dim.lower)
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

    def decodeParallelPragma2(self, pragma):
        #remove opening curly bracket
        trimmed_bracket = pragma.split('{')[0]

        splitCommas = trimmed_bracket.split(',')
        clauses = {}
        for clause in splitCommas:
            if '=' in clause:
                clauses[clause.split('=')[0].strip()] = clause.split('=')[1].strip()

        targets = None
        if ('TARGET' not in clauses):
            print("No target specified, defaulting to OPENMP")
            targets = ['OpenMP']
        else:
            targets =  clauses['TARGET'].split('/') 
         
            for target in targets:
                if (target not in ['OpenMP','OpenMPSingleColumn','OpenACCSingleColumn']):
                    print(colored("Unkown target clause !!!!", "red"))
                    exit(1)

        
        name = clauses['NAME'] if 'NAME' in clauses else None
        return(targets, name)

    def addTransform(self, subroutine, transform):
        if subroutine not in self.subroutines_to_transform:
            self.subroutines_to_transform[subroutine] = {transform}
        else:
            self.subroutines_to_transform[subroutine].add(transform)
 


    def npromaToFieldAPI(self, routine):
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

        init_calls = ()
        final_calls = ()


        for var in FindVariables().visit(routine.spec):

            if isinstance(var, Array):

                to_transform = False

                firstdim = var.dimensions[0] 

                if isinstance(firstdim, RangeIndex):
                    firstdim = firstdim.upper

                if firstdim in params.nproma_aliases:
                    isArgument = var in routine.arguments
                
                    #print("variable to FAPI : ", var, type(var))
                    #print("type : ", type(var.type.kind.name))

                    array_type_dim = f'{len(var.dimensions)+1}{var.type.kind.name[-2:]}'


                    arrays_types_dimensions.add(array_type_dim)
                        
                    new_var = Variable(     name=f'Y{("D" if isArgument else "L")}_{var.name}', 
                            type=SymbolAttributes(  DerivedType(name=f'FIELD_{array_type_dim}_ARRAY'), 
                                                                    intent=var.type.intent, 
                                                                    scope = routine
                                                                    ) 
                                    )

                    nproma_arrays_map[var] = new_var
                    nproma_names_map[var.name] = new_var
                    nproma_arrays_dimensions[var.name] = var.dimensions

                    if not isArgument :

                        init_calls += (CallStatement(name = DeferredTypeSymbol(name='INIT', parent = new_var), 
                                                            
                                                            kwarguments = self.boundariesArgument(var) + (('PERSISTENT', LogicLiteral(True)), ),
                                                            scope=routine) ,)

                        final_calls += (CallStatement(name = DeferredTypeSymbol(name='FINAL', parent = new_var), 
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


        routine.body.prepend(init_calls)
        routine.body.append(final_calls)

        # Create imports for module ARRAY_MOD and UTIL_ARRAY_xD_MOD
        #routine.spec.prepend(Import(module="ARRAY_MOD", 
        #        symbols=tuple([DeferredTypeSymbol(name=f'ARRAY_{varl}D') for varl in arrays_types_dimensions]) ))

        #for type_dim in arrays_types_dimensions:
        #    routine.spec.prepend(Import(module=f'UTIL_ARRAY_{type_dim}D_MOD'))

        # We will very likely need to declare FIELD_BASIC variables
        routine.spec.prepend(Import(module="FIELD_BASIC_MODULE", symbols=(DeferredTypeSymbol(name='FIELD_BASIC'),)))

        routine.spec.prepend(Import(module="YOMPARALLELMETHOD"))
        routine.spec.prepend(Import(module="FIELD_ARRAY_MODULE"))
        routine.spec.prepend(Import(module="STACK_MOD"))
        routine.spec.prepend(Import(module='stack.h', c_import=True))


        return ([nproma_names_map[var] for var in nproma_names_map], locals_dimensions)

    def makeOpenMPLoop(self, routine, boundary_variable, region):

        init_bounds = CallStatement(name = DeferredTypeSymbol( name='INIT',
                                                               parent = boundary_variable),
                                    arguments=(Variable(name=params.cpg_opts_variable)),
                                    scope = routine
                                    )
        loop_pragma = Pragma(   keyword="OMP", 
                                content = "PARALLEL DO FIRSTPRIVATE ("+ boundary_variable.name +") PRIVATE (" + 
                                self.block_counter.name + self.getPrivateNameList(region) + ") "
                            )

        #We prepare the call to update_view for the block index
        loop_new_statements = CallStatement(name = DeferredTypeSymbol(name='UPDATE',
                                                                      parent = boundary_variable), 
                                            arguments=(self.block_counter), 
                                            scope=routine
                                            )
   
        new_loop = self.blocks_loop.clone(body=(loop_new_statements, region.body,), pragma=(loop_pragma,) )

        return (init_bounds, new_loop,)

    def makeOpenMPSCCLoop(self, routine, boundary_variable, region):

        loop_pragma = Pragma(   keyword="OMP", 
                                content = "PARALLEL DO PRIVATE ("+ self.block_counter.name + ', JLON,'+
                                boundary_variable.name + ', YSTACK' + self.getPrivateNameList(region) + ") "
                            )

        columns_loop = self.makeColumnsLoop(routine, boundary_variable, region)              
    
        new_loop = self.blocks_loop.clone(body=(columns_loop,), pragma=(loop_pragma,) )

        return (new_loop,)


    def makeOpenACCSCCLoop(self, routine, boundary_variable, region):

        loop_pragma = Pragma(   keyword="ACC", 
                                    content = "PARALLEL LOOP GANG FIRSTPRIVATE ("+ boundary_variable.name +") PRIVATE (IBL)"
                            )

        columns_loop = self.makeColumnsLoop(routine, boundary_variable, region, 'OpenACC')              

        new_loop = self.blocks_loop.clone(body=(columns_loop,), pragma=(loop_pragma,) )

        return (new_loop,)

    def addJLONDefinition(self, routine):
        # jlon_found = False
        if "JLON" not in [v.name for v in routine.variables]:
            routine.variables += (Variable( name='JLON', 
                                            type=SymbolAttributes(BasicType.INTEGER,  kind=Variable(name='JPIM')))
                                ,)


    def makeColumnsLoop(self, routine, boundary_variable, region, directive = None):
        # Inner loop on JLON

        # YLCPG_BNDS = YDCPG_BNDS
        #assignments = (Assignment(lhs = local_boundary_variable, rhs = boundary_variable), ) 

        # YLCPG_BNDS%KIDIA = JLON
        assignments = (Assignment( lhs = Variable(name='KIDIA', parent = boundary_variable), 
                                    rhs = Variable(name='JLON') ), )
        # YLCPG_BNDS%KFDIA = JLON
        assignments += (Assignment( lhs = Variable(name='KFDIA', parent = boundary_variable), 
                                    rhs = Variable(name='JLON') ), )
        # Compute local stack variables boundaries
        #assignments += (Assignment( lhs = Variable(name='YLSTACK'), rhs = Variable(name='YDSTACK')), )

        # YLSTACK%L8 = stack_l8 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        l8_string = 'stack_l8 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='L8', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(l8_string, scope=routine)  ),)                                                   


        # YLSTACK%U8 = stack_u8 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        u8_string = 'stack_u8 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='U8', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(u8_string, scope=routine)  ),)                                                   


        # YLSTACK%L4 = stack_l4 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        l4_string = 'stack_l4 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='L4', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(l4_string, scope=routine)  ),)                                                   


        #YLSTACK%U4 = stack_u4 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        u4_string = 'stack_u4 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='U4', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(u4_string, scope=routine)  ),)                                                   

        
        


        upper_bound_string = 'MIN (YDCPG_OPTS%KLON, YDCPG_OPTS%KGPCOMP - (JBLK - 1) * YDCPG_OPTS%KLON)'

        # parse_fparser_expression(upper_bound_string, scope=routine)

        if (directive == 'OpenACC'):
            loop_pragma = (Pragma(keyword="ACC", content = f'LOOP VECTOR PRIVATE({boundary_variable},  YLSTACK)'),)
        else : 
            loop_pragma = None


        columns_loop = Loop(    variable = Variable(name='JLON'), 
                                bounds = LoopRange( (IntLiteral(1),
                                                    parse_fparser_expression(upper_bound_string, scope=routine))),
                                body = (assignments, region.body,),
                                pragma = loop_pragma
                            )


        if not self.has_column_loops:
            self.addJLONDefinition(routine)
            self.has_column_loops = True

        return columns_loop

    def makeRegionSyncCall(self, routine, target, region_num, body, spec, nproma_arrays):
        # Create synchronisation call 
        # building a subroutine that contains only the current region
        # Apply the Makesync transformation, then append it as a member routine
        sync_name = routine.name+"_PARALLEL_2_" + str(region_num) + "_SYNC_" + target
        sync_routine = routine.clone(body = body, spec = spec, name = sync_name)

        sync_routine.apply(RemoveLoops())
        add_suffix_transform = AddSuffixToCalls(suffix='_SYNC_' + target)
        sync_routine.apply(add_suffix_transform)
        for subroutine in add_suffix_transform.routines_called:
            self.addTransform(subroutine, 'SYNC_' + target)
        sync_routine.apply(RemoveUnusedImports())
           
        print("makesync")
        syncTransformation = MakeSync(pointerType=target.lower(), sections = (sync_routine.body,), nproma_arrays=nproma_arrays )
        sync_routine.apply(syncTransformation)
        print("after makesync, found ", syncTransformation.total_FAPI_pointers, " pointers")

        # We keep the maximum amount of local pointers needed to pass FieldAPI variables to Sync subroutines
        #total_FAPI_pointers = max (total_FAPI_pointers, syncTransformation.total_FAPI_pointers)

        sync_routine.apply(RemoveComments())
        sync_routine.apply(RemovePragmas())
        sync_routine.apply(RemoveEmptyConditionals())

        print("sync routine built")
        #return (sync_name, syncTransformation.total_FAPI_pointers)
        # We insert the sync routine as a member routine of its caller
        # Therefore we do not need any variable passing or declaration as we manipulate only variables of the caller
        if not routine.contains:
            routine.contains = Section(Intrinsic(text="CONTAINS"))

        includes = Section(tuple([n for n in FindNodes(Import).visit(sync_routine.spec) if n.c_import]))

        routine.contains.append(sync_routine.clone(body=sync_routine.body, args=None, spec=includes, contains=None))

        return (sync_name, syncTransformation.total_FAPI_pointers)


    def containsOnlyCalls(self, region):
        for n in FindNodes(Node).visit(region.body):
            if not (isinstance(n, CallStatement) 
                    or isinstance(n, Comment)
                    or isinstance(n, CommentBlock)):
                        return False 
        return True

    def makeParallelSubroutine(self, routine, region, region_num, scc, target='host'):
        print("make parallel subroutine : ", routine, target, region)
        new_subroutine = routine.clone( name = routine.name + "_PARALLEL_" + str(region_num), 
                                                body = region.body,
                                                spec = routine.spec.clone(),
                                                contains = None)

        add_suffix_transform = AddSuffixToCalls(  suffix=('_SCC'), additional_variables=['YDSTACK=YLSTACK'] ) 
        new_subroutine.apply(add_suffix_transform)
        for subroutine in add_suffix_transform.routines_called:
            self.addTransform(subroutine, 'SCC_'+target.upper())

        new_subroutine.apply(RemoveUnusedImports())
        if scc:
            true_symbols, false_symbols=logical_lst.symbols()
            false_symbols.append('LHOOK')

            scc_transform_routine(new_subroutine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols, FieldAPI_pointers=self.FieldAPI_pointers, is_node=True)

        new_subroutine.apply(RemoveComments())
        new_subroutine.apply(RemovePragmas())
        new_subroutine.apply(FieldAPIPtr(pointerType=target))  
        
        return new_subroutine.body
        # self.outlined_routines.append(new_subroutine)


    def transform_parallel_regions(self, routine):

        unmodified_spec = routine.spec.clone()

        # Add block counter declaration
        routine.variables += (self.block_counter,)

        # Use custom visitor to add PARALLEL suffix and stack variable to CallStatements outside PragmaRegions
        add_suffix_transform = AddSuffixToCalls(suffix='_PARALLEL_2', custom_visitor = FindNodesOutsidePragmaRegion, additional_variables = ['YDSTACK'])
        routine.apply(add_suffix_transform)


        for subroutine in add_suffix_transform.routines_called:
            self.addTransform(subroutine, 'PARALLEL')
       

        # Transform arrays of NPROMA size into ARRAY_XD types
        nproma_arrays, local_dimensions=self.npromaToFieldAPI(routine)

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
            print("region : ", region.pragma, region.pragma.content)
            region_num = region_num + 1

            #(target, scc, outline, directive) = self.decodeParallelPragma(region.pragma.content)
            (targets, name) = self.decodeParallelPragma2(region.pragma.content)
   
            if not name:
                name = f'REGION_{region_num}'

            print("targets and name : ", targets, name)
   
            archs = [False, False]
            for target in targets:  
                if 'OpenMP' in target:
                    archs[0] = True
                elif 'OpenACC' in target:
                    archs[1] = True

            #Generate sync subroutines for [HOST, DEVICE]
            sync_names = [None, None]
            if archs[0]:
                (sync_names[0], sync_FAPI_pointers) = self.makeRegionSyncCall(routine, 'HOST', region_num, region.body, unmodified_spec, nproma_arrays)
                # We keep the maximum amount of local pointers needed to pass FieldAPI variables to Sync subroutines
                total_FAPI_pointers = max (total_FAPI_pointers, sync_FAPI_pointers)
            if archs[1]:
                (sync_names[1], sync_FAPI_pointers) = self.makeRegionSyncCall(routine, 'DEVICE', region_num, region.body, unmodified_spec, nproma_arrays)
                total_FAPI_pointers = max (total_FAPI_pointers, sync_FAPI_pointers)

            print("sync_names ? ", sync_names)
             
            abort_calls=()
            for call in FindNodes(CallStatement).visit(region.body):
                print("call in region : ", call)
                self.addTransform(call.name.name, 'ABORT')
                abort_calls += (call.clone(name= DeferredTypeSymbol(name=f'{call.name}_ABORT')),)
            

            
            for index,target in enumerate(targets):

                call_sync = CallStatement(name = DeferredTypeSymbol(name=sync_names[1] if 'OpenACC' in target else sync_names[0]), 
                                                arguments=(), scope=None )
  
                #Create the replacement body for the region for current target

                if (target == 'OpenMP'):
                    new_region = FieldAPIPtr(pointerType='host').transform_node(region, routine)
                    new_loop = self.makeOpenMPLoop(routine, local_boundary_variable, new_region)
                    new_body = (call_sync, new_loop,)

                elif (target == 'OpenMPSingleColumn'):
                    # Check for outline necessity : if region only contains calls, directly call SCC variants
                    # Otherwise, create a contained subroutine with the transformed region
                    if self.containsOnlyCalls(region):
                        new_region = FieldAPIPtr(pointerType='host').transform_node(region, routine)
                        add_suffix_transform = AddSuffixToCalls(suffix=('_SCC_HOST'), additional_variables=['YDSTACK=YLSTACK'] )
                        new_region = add_suffix_transform.transform_node(new_region, routine)

                        for subroutine in add_suffix_transform.routines_called:
                            self.addTransform(subroutine, 'SCC_HOST')
                    else:
                        new_region = self.makeParallelSubroutine(routine, region, region_num, scc=True, target='host')
                        


                    new_loop = self.makeOpenMPSCCLoop(routine, local_boundary_variable, new_region)
                    new_body = (call_sync, new_loop,)                    
                
                elif (target == 'OpenACCSingleColumn'):
                    if self.containsOnlyCalls(region):
                        new_region = FieldAPIPtr(pointerType='device').transform_node(region, routine)

                        add_suffix_transform = AddSuffixToCalls(  suffix=('_SCC_DEVICE'), additional_variables=['YDSTACK=YLSTACK'] )
                        new_region = add_suffix_transform.transform_node(new_region, routine)
                        for subroutine in add_suffix_transform.routines_called:
                            self.addTransform(subroutine, 'SCC_DEVICE')
                    else:
                        new_region = self.makeParallelSubroutine(routine, region, region_num, scc=True, target='device')
                    
                    new_loop = self.makeOpenACCSCCLoop(routine, local_boundary_variable, new_region)
                    new_body = (call_sync, new_loop,)                    
                else :
                    print("Unknown target !!!! ", target)
                    new_body = (region.body,)


                #Create nested conditionals with call to LPARALLELMETHOD in conditions

                expr = InlineCall(function=DeferredTypeSymbol(name='LPARALLELMETHOD'), 
                        parameters=(StringLiteral(target.upper()),StringLiteral(f'{routine.name}:{name}'), ))
             
                if (index == 0) : 
                    outer_cond = Conditional(condition=expr, body=new_body)
                    current_cond = outer_cond

                elif (index < len(targets)):
                    new_cond = Conditional(condition=expr ,body=new_body )
                    current_cond._update(else_body = new_cond)
                    current_cond._update(has_elseif = True)
                    current_cond = new_cond

                #Final ELSE branch only contains a call to ABOR1
                current_cond._update(else_body = (CallStatement(name = DeferredTypeSymbol(name='ABOR1'),
                                                                arguments=(StringLiteral(f'{routine.name}_PARALLEL : METHOD WAS NOT FOUND')),
                                                                scope=routine
                                                                ),
                                                )
                                    )

            routine.body = Transformer({region:(abort_calls,outer_cond)}).visit(routine.body)

            addFieldAPIPointers(routine, total_FAPI_pointers)

            routine.variables += (local_boundary_variable, local_stack,)   

            # return



            # if not outline:
            #     region = AddSuffixToCalls(  suffix=('_SINGLE_COLUMN' if scc else '') + '_FIELD_API_' + target, 
            #                                 node=region, 
            #                                 additional_variables=['YLSTACK'] if scc else []).transform_subroutine(routine = routine)

            # if outline:

            #     new_subroutine = routine.clone( name = routine.name + "_PARALLEL_" + str(region_num), 
            #                                     body = region.body,
            #                                     spec = routine.spec.clone(),
            #                                     contains = None)

            #     new_subroutine.apply(AddSuffixToCalls(  suffix=('_SINGLE_COLUMN' if scc else '') + '_FIELD_API_' + target, 
            #                                 # node=region, 
            #                                 additional_variables=['YLSTACK'] if scc else []))

            #     new_subroutine.apply(RemoveUnusedImports())

            #     if scc:
            #         true_symbols, false_symbols=logical_lst.symbols()
            #         false_symbols.append('LHOOK')

            #         scc_transform_routine(new_subroutine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols)

            #     new_subroutine.apply(FieldAPIPtr(pointerType=target.lower()))
            #     new_subroutine.apply(RemoveComments())
            #     new_subroutine.apply(RemovePragmas())
  
                # We do not want JLON and YSTACK in the arguments
                # remove_transform = RemoveUnusedVariables(['JLON', 'YLSTACK'])
                # new_subroutine.apply(remove_transform)

                # new_subroutine.spec.append(Pragma(keyword='acc', content='routine seq'))
                # if directive == 'openacc':
                #     new_subroutine.apply(AddACCRoutineDirectives())

                # self.outlined_routines.append(new_subroutine)


                # call_arguments = remove_transform.used_symbols
                # call_arguments = SubstituteExpressions({
                #                         Scalar(name=boundary_variable.name):Scalar(name=local_boundary_variable.name),
                #                         Scalar(name='YDSTACK'):Scalar(name='YLSTACK')
                #                         }).visit(call_arguments)

                # call_parallel = CallStatement(name=new_subroutine.procedure_symbol, arguments=call_arguments)



                # if directive == 'openacc':
                #     columns_loop = self.makeOpenACCColumnsLoop(call_parallel, boundary_variable, local_boundary_variable, routine)
                
            # if not outline:
            #     loop_body = region.body
            # else:
            #     if directive == 'openacc':
            #         loop_body = columns_loop
            #     else:
            #         loop_body = call_parallel

            # blocks_loop = Loop( variable=self.block_counter, 
            #                     bounds=LoopRange(( IntLiteral(1), self.block_symbol )), 
            #                     body = (loop_new_statements, loop_body,), 
            #                     pragma = (loop_pragma,)   
            #                     )
            # routine.body = Transformer({region:(call_sync,blocks_loop,)}).visit(routine.body)

            # if not outline:
            #     # Default behaviour : apply FieldAPI transformation to the loop body
            #     routine.apply(FieldAPIPtr(pointerType=target.lower(), node=blocks_loop ))


        # We remove remaining imports only after having treated all Parallel blocks 
        # routine.apply(RemoveUnusedImports())

        # routine.variables += (local_boundary_variable, local_stack,)   

        # routine.body = Transformer(regions_map).visit(routine.body)
        # addFieldAPIPointers(routine, total_FAPI_pointers)

    def transform_subroutine(self, routine, **kwargs):


        self.transform_parallel_regions(routine)

        # In a routine with parallel blocks, we expect to have received a stack variable
        #routine.spec.prepend(Import(module="STACK_MOD"))
        #routine.arguments += (Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK"), intent='in')) ,)

        # get updatable fields and change their intent into INOUT
        FieldAPI_updatables = retrieve('../../update_view.dat')

        #                      'CPG_BNDS_TYPE'
        FieldAPI_updatables[params.boundaries_type] = 1
        var_map = {}
        #for var in routine.arguments :
        #    if (isinstance(var.type.dtype, DerivedType)):
        #        typename = var.type.dtype.name
        #        if typename in FieldAPI_updatables:
        #            if FieldAPI_updatables[typename]:
        #                new_type=var.type.clone(intent='inout')
        #                var_map[var]= var.clone(type=new_type)
        
        #routine.spec = SubstituteExpressions(var_map).visit(routine.spec)

