from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType, FindTypedSymbols, FindVariables, SubstituteExpressions  )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList, LoopRange, IntLiteral
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from storable import retrieve

from syncTransforms import MakeSync, addFieldAPIPointers

from fieldAPITransforms import FieldAPIPtr, get_fieldAPI_variables, get_fieldAPI_member, is_fieldAPI_ARRAY
from commonTransforms import InlineMemberCalls, RemoveComments, RemovePragmas, RemovePragmaRegions, RemoveEmptyConditionals, \
                                AddSuffixToCalls, RemoveLoops, RemoveUnusedVariables, AddACCRoutineDirectives, \
                                RemoveUnusedImports, FindNodesOutsidePragmaRegion, ReplaceArguments, AddFieldAPISuffixToCalls

from arpege_parameters import params

# import errors should have already been caught in acdc_pragmas.py
import logical_lst
from openacc_transform import scc_transform_routine, alloc_temp

from termcolor import colored

class MakeParallel(Transformation):
    
    def __init__(self,  FieldAPI_pointers={}, top_level_routine=False):        
        self.FieldAPI_pointers = FieldAPI_pointers
        self.top_level_routine = top_level_routine
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

        lower_bounds.append(self.block_min_symbol)
        upper_bounds.append(self.block_max_symbol)

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

        field_new_calls = ()
        field_delete_calls = ()

        for var in FindVariables().visit(routine.spec):

            if isinstance(var, Array):

                is_FAPI_pointer = var.name in self.FieldAPI_pointers

                is_nproma = False
                if not is_FAPI_pointer :
                    firstdim = var.dimensions[0] 
                
                    #Naively considering ranges with (xxx:nproma) are nproma-sized
                    if isinstance(firstdim, RangeIndex):
                        firstdim = firstdim.upper

                    is_nproma = firstdim in params.nproma_aliases
                
                if is_FAPI_pointer or is_nproma:

                    is_argument = var in routine.arguments
                
                    #print("variable to FAPI : ", var, type(var))
                    #print("type : ", type(var.type.kind.name))

                    array_type_dim = f'{len(var.dimensions)+1}{var.type.kind.name[-2:]}'


                    arrays_types_dimensions.add(array_type_dim)
                        
                    new_var = Variable(     name=f'Y{("D" if is_argument else "L")}_{var.name}', 
                                            type=var.type.clone(dtype=DerivedType(name=f'FIELD_{array_type_dim}'),
                                                                kind=None, 
                                                                # Abstract types require polymorphism
                                                                polymorphic=True, 
                                                                # Dummy arrays can stay target, local arrays have to become
                                                                # pointers which conflict with target attribute
                                                                #target = var.type.target if not is_argument else False,
                                                                target = False,
                                                                #pointer = False if (is_argument and var.type.target) else True
                                                                #pointer = True if is_argument else var.type.pointer
                                                                pointer = True
                                                                )
                                    )

                    nproma_arrays_map[var] = new_var
                    nproma_names_map[var.name] = new_var
                    nproma_arrays_dimensions[var.name] = var.dimensions

                    if not (is_argument or is_FAPI_pointer):
                        # Only create locally declared FieldAPI object
                        #if not is_argument:ddi    
                        field_new_calls += (CallStatement(name = DeferredTypeSymbol(name='FIELD_NEW'),
                                                                arguments = new_var,
                                                                kwarguments = self.boundariesArgument(var) + (('PERSISTENT', LogicLiteral(True)), ),
                                                                scope=routine
                                                            ),
        
                                               )

                        # For all FieldAPI object : call copy_object to have the structure offloaded 
                        # (device pointer is already present, but not the whole object)
                        #associated_call = InlineCall(function=DeferredTypeSymbol(name='ASSOCIATED'), parameters=(new_var,) ) 
                        #present_call = InlineCall(function=DeferredTypeSymbol(name='PRESENT'), parameters=(new_var,) )
                        copy_call = CallStatement(name = DeferredTypeSymbol(name='COPY_OBJECT', 
                                                                            parent = new_var
                                                                            ),
                                                    arguments=(), scope=routine
                                                    ) 
                        # Arguments might be unassociated pointers
                        #if is_argument and var.type.optional :
                        #    field_new_calls += (Conditional(condition=present_call, body = (copy_call,), inline = True), )
                        #else:
                        #if not (is_argument):
                        field_new_calls += (copy_call,)
                                             
                        # Symmetrically, call WIPE_OBJECT at the end to remove the structure from GPU memory
                        wipe_call = CallStatement(name = DeferredTypeSymbol(name='WIPE_OBJECT', 
                                                                            parent = new_var
                                                                           ),
                                                    arguments=(), scope=routine
                                                    )

                        #if is_argument and var.type.optional:

                        #    field_delete_calls += (Conditional(condition=present_call, body = (wipe_call,), inline = True), )
                        #else:
                        # Current dirty fix : do not wipe arguments to ensure presence in the caller
                        # We really should have this when they are constructed/destroyed !!!!
                        #if not is_argument:
                        field_delete_calls += (wipe_call,)
                         
                        # Finally, actually delete the local FieldAPI objects
                        #if not is_argument:
                        delete_call = CallStatement(name = DeferredTypeSymbol(name='FIELD_DELETE'), arguments = new_var)
                            #field_delete_calls += (Conditional(condition=associated_call, body = (delete_call,), inline = True), )                 
                        field_delete_calls += (delete_call,)
                        #init_calls += (CallStatement(name = DeferredTypeSymbol(name='INIT', parent = new_var), 
                        #                                    
                        #                                    kwarguments = self.boundariesArgument(var) + (('PERSISTENT', LogicLiteral(True)), ),
                        #                                    scope=routine) ,)

                        #final_calls += (CallStatement(name = DeferredTypeSymbol(name='FINAL', parent = new_var), 
                        #                                        arguments=(), scope=routine), )


                                       


        #print("nproma arrays map : ", nproma_arrays_map)

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


        routine.body.prepend(field_new_calls)
        routine.body.append(field_delete_calls)

        # Create imports for module ARRAY_MOD and UTIL_ARRAY_xD_MOD
        #routine.spec.prepend(Import(module="ARRAY_MOD", 
        #        symbols=tuple([DeferredTypeSymbol(name=f'ARRAY_{varl}D') for varl in arrays_types_dimensions]) ))

        #for type_dim in arrays_types_dimensions:
        #    routine.spec.prepend(Import(module=f'UTIL_ARRAY_{type_dim}D_MOD'))

        # We will very likely need to declare FIELD_BASIC variables
        routine.spec.prepend(Import(module="FIELD_BASIC_MODULE", symbols=(DeferredTypeSymbol(name='FIELD_BASIC'),)))

        routine.spec.prepend(Import(module="YOMPARALLELMETHOD"))
        #routine.spec.prepend(Import(module="FIELD_ARRAY_MODULE"))
        routine.spec.prepend(Import(module="FIELD_MODULE"))
        routine.spec.prepend(Import(module="FIELD_FACTORY_MODULE"))

        routine.spec.prepend(Import(module="STACK_MOD"))
        routine.spec.prepend(Import(module='stack.h', c_import=True))


        return ([nproma_names_map[var] for var in nproma_names_map], locals_dimensions)

    def makeOpenMPLoop(self, routine, boundary_variable, region):

       
        var_map={}
        for var in FindVariables().visit(region.body):
            if var.name == 'YDCPG_BNDS':
                var_map[var] = boundary_variable
            elif (var.name == 'YDCPG_BNDS%KIDIA' or var.name == 'YDCPG_BNDS%KFDIA'):
                var_map[var] = var.clone(parent=boundary_variable)

        new_body = SubstituteExpressions(var_map).visit(region.body)
        new_region = region.clone(body = new_body)
        
        new_region = FieldAPIPtr(pointerType='host').transform_node(new_region, routine)
        #new_region.apply( FieldAPIPtr(pointerType='host'))

        add_suffix_transform = AddFieldAPISuffixToCalls(argument=boundary_variable.name) 
        add_suffix_transform.transform_node(new_region, routine)
        for subroutine in add_suffix_transform.routines_called:
            self.addTransform(subroutine, 'FIELD_API')
        


        init_bounds = CallStatement(name = DeferredTypeSymbol( name='INIT',
                                                               parent = boundary_variable),
                                    arguments=(Variable(name=params.cpg_opts_variable)),
                                    scope = routine
                                    )
        loop_pragma = Pragma(   keyword="OMP", 
                                content = "PARALLEL DO FIRSTPRIVATE ("+ boundary_variable.name +") PRIVATE (" + 
                                self.block_counter.name + self.getPrivateNameList(new_region) + ") "
                            )

        #We prepare the call to update_view for the block index
        loop_new_statements = CallStatement(name = DeferredTypeSymbol(name='UPDATE',
                                                                      parent = boundary_variable), 
                                            arguments=(self.block_counter), 
                                            scope=routine
                                            )
   
        new_loop = self.blocks_loop.clone(body=(loop_new_statements, new_region.body,), pragma=(loop_pragma,) )

        return (init_bounds, new_loop,)

    def makeOpenMPSCCLoop(self, routine, boundary_variable, region):

        loop_pragma = Pragma(   keyword="OMP", 
                                content = "PARALLEL DO PRIVATE ("+ self.block_counter.name + ', JLON,'+
                                boundary_variable.name + ', YLSTACK' + self.getPrivateNameList(region) + ") "
                            )

        columns_loop = self.makeColumnsLoop(routine, boundary_variable, region)              
    
        new_loop = self.blocks_loop.clone(body=(columns_loop,), pragma=(loop_pragma,) )

        return (new_loop,)


    def makeOpenACCSCCLoop(self, routine, boundary_variable, region):

        loop_pragma = Pragma(   keyword="ACC", 
                                    content = "PARALLEL LOOP GANG DEFAULT(PRESENT) PRIVATE (" + params.block_counter + ")"
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
        l8_string = f'stack_l8 (YSTACK, {params.block_counter}-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='L8', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(l8_string, scope=routine)  ),)                                                   


        # YLSTACK%U8 = stack_u8 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        u8_string = f'stack_u8 (YSTACK, {params.block_counter}-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='U8', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(u8_string, scope=routine)  ),)                                                   


        # YLSTACK%L4 = stack_l4 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        l4_string = f'stack_l4 (YSTACK, {params.block_counter}-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'
         

        assignments += (Assignment( lhs = Variable(name='L4', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(l4_string, scope=routine)  ),)                                                   


        #YLSTACK%U4 = stack_u4 (YSTACK, JBLK-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)
        u4_string = f'stack_u4 (YSTACK, {params.block_counter}-YDCPG_OPTS%JBLKMIN+1, YDCPG_OPTS%KGPBLKS)'

        assignments += (Assignment( lhs = Variable(name='U4', parent = Variable(name='YLSTACK')),
                                 rhs = parse_fparser_expression(u4_string, scope=routine)  ),)                                                   

        
        


        upper_bound_string = f'MIN (YDCPG_OPTS%KLON, YDCPG_OPTS%KGPCOMP - ({params.block_counter} - 1) * YDCPG_OPTS%KLON)'

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

    def makeSyncRegion(self, routine, target, region_num, body, spec, nproma_arrays):
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
           
        syncTransformation = MakeSync(pointerType=target.lower(), sections = (sync_routine.body,), nproma_arrays=nproma_arrays )
        sync_routine.apply(syncTransformation)

        # We keep the maximum amount of local pointers needed to pass FieldAPI variables to Sync subroutines
        #total_FAPI_pointers = max (total_FAPI_pointers, syncTransformation.total_FAPI_pointers)

        sync_routine.apply(RemoveComments())
        sync_routine.apply(RemovePragmas())
        sync_routine.apply(RemoveEmptyConditionals())

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

    def makeParallelSubroutine(self, routine, region, region_num, local_boundary_variable, scc, target='host') :
        # Optimisation : when a block only contains calls, simply transform the calls
        if self.containsOnlyCalls(region):
            new_region = FieldAPIPtr(pointerType=target).transform_node(region, routine)
            add_suffix_transform = AddSuffixToCalls(suffix=('_SCC_'+target.upper()), additional_variables=['YDSTACK=YLSTACK'] )
            new_region = add_suffix_transform.transform_node(new_region, routine)
                
            for call in FindNodes(CallStatement).visit(new_region.body):
                ReplaceArguments(call, {'YDCPG_BNDS%KIDIA':Variable(name='KIDIA', parent = local_boundary_variable), \
                                        'YDCPG_BNDS%KFDIA':Variable(name='KFDIA', parent = local_boundary_variable)})
      
            for subroutine in add_suffix_transform.routines_called:
                self.addTransform(subroutine, 'SCC_'+target.upper())
            return new_region
        else:
            # We have to create a complete subroutine to apply the SCC transformation
            new_subroutine = routine.clone( name = routine.name + "_PARALLEL_" + str(region_num), 
                                            body = region.body,
                                            spec = routine.spec.clone(),
                                            contains = None)

            add_suffix_transform = AddSuffixToCalls(  suffix=('_SCC_'+target.upper()), additional_variables=['YDSTACK=YLSTACK'] ) 
            # Apply this transformation as a node, since we want to update imports in the complete routine
            new_subroutine = add_suffix_transform.transform_node(new_subroutine, routine)

            for call in FindNodes(CallStatement).visit(new_subroutine.body):
                ReplaceArguments(call, {'YDCPG_BNDS%KIDIA':Variable(name='KIDIA', parent = local_boundary_variable), \
                                        'YDCPG_BNDS%KFDIA':Variable(name='KFDIA', parent = local_boundary_variable)})

            for subroutine in add_suffix_transform.routines_called:
                self.addTransform(subroutine, 'SCC_'+target.upper())

            new_subroutine.apply(RemoveUnusedImports())
            if scc:
                true_symbols, false_symbols=logical_lst.symbols()
                false_symbols.append('LHOOK')

                prefixed_FieldAPI_pointers = {'YL_'+key:value for key,value in self.FieldAPI_pointers.items()}

                scc_transform_routine(new_subroutine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols, FieldAPI_pointers=prefixed_FieldAPI_pointers, is_node=True)

            new_subroutine.apply(RemoveComments())
            new_subroutine.apply(RemovePragmas())
            
            # Only apply node FieldAPIPtr transformation : we do not want to insert JBLK declaration and assignment
            new_subroutine.body = FieldAPIPtr(pointerType=target).transform_node(new_subroutine.body, new_subroutine)
        
            return new_subroutine.body
        # self.outlined_routines.append(new_subroutine)


    def transform_parallel_regions(self, routine):

        self.fieldAPI_variables = get_fieldAPI_variables(routine)
        unmodified_spec = routine.spec.clone()

        # Add block counter declaration
        routine.variables += (self.block_counter,)


        # Use custom visitor to add PARALLEL suffix and stack variable to CallStatements outside PragmaRegions
        add_suffix_transform = AddSuffixToCalls(suffix='_PARALLEL2', custom_visitor = FindNodesOutsidePragmaRegion)
        routine.apply(add_suffix_transform)


        for subroutine in add_suffix_transform.routines_called:
            self.addTransform(subroutine, 'PARALLEL')
       
        # Transform arrays of NPROMA size into ARRAY_XD types
        nproma_arrays, local_dimensions=self.npromaToFieldAPI(routine)
       
        # Transform arguments in the call
        for call in FindNodesOutsidePragmaRegion(CallStatement).visit(routine.body):
            # Change name of positional arguments receiving a fieldAPI array
            # (since their name changes also in the subroutine)
            new_kwargs = ()
            to_update = False
            for couple in call.kwarguments:
                if couple[1] in nproma_arrays:
                    #Left part of positional argument is a string (does not exist in current scope)
                    new_kwargs += ('YD_'+couple[0], couple[1]),
                    to_update = True
                else:
                    new_kwargs += (couple[0], couple[1]),

            call._update(kwarguments = new_kwargs)


            # If array view are passed as arguments, change them into their FieldAPI counterpart
            new_args=()
            updated_args = False
            for var in call.arguments:
                is_new_var = False
                if isinstance(var, DeferredTypeSymbol):
                    base_name = var.name_parts[0]
                    if base_name in self.fieldAPI_variables:
                        member_name = var.name_parts[1:]
                        member_name = [self.fieldAPI_variables[base_name]] + member_name
                        fAPI_member = get_fieldAPI_member(member_name)
                        if fAPI_member :
                            new_args += (Variable(name = fAPI_member[0], parent=var.parent, scope=routine),)
                            updated_args = True
                            is_new_var = True
                if not is_new_var:
                    new_args +=(var,)
            if updated_args:
                call._update(arguments = new_args)
                

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
            #print("region : ", region.pragma, region.pragma.content)
            region_num = region_num + 1

            #(target, scc, outline, directive) = self.decodeParallelPragma(region.pragma.content)
            (targets, name) = self.decodeParallelPragma2(region.pragma.content)
   
            if not name:
                name = f'REGION_{region_num}'

            #print("targets and name : ", targets, name)
   
            archs = [False, False]
            for target in targets:  
                if 'OpenMP' in target:
                    archs[0] = True
                elif 'OpenACC' in target:
                    archs[1] = True

            #Generate sync subroutines for [HOST, DEVICE]
            sync_names = [None, None]
            if archs[0]:
                (sync_names[0], sync_FAPI_pointers) = self.makeSyncRegion(routine, 'HOST', region_num, region.body, unmodified_spec, nproma_arrays)
                # We keep the maximum amount of local pointers needed to pass FieldAPI variables to Sync subroutines
                total_FAPI_pointers = max (total_FAPI_pointers, sync_FAPI_pointers)
            if archs[1]:
                (sync_names[1], sync_FAPI_pointers) = self.makeSyncRegion(routine, 'DEVICE', region_num, region.body, unmodified_spec, nproma_arrays)
                total_FAPI_pointers = max (total_FAPI_pointers, sync_FAPI_pointers)

             
            abort_calls=()
            for call in FindNodes(CallStatement).visit(region.body):
                #print("call in region : ", call)
                if call.name not in params.ignore_abort:
                    self.addTransform(call.name.name, 'ABORT')
                    # Adjust arguments : variables turned into fieldAPI should now pass their PTR
                    to_update = False
                    abort_call = call.clone() #name= DeferredTypeSymbol(name=f'{call.name}_ABORT'))
                    new_args = ()
                    for arg in abort_call.arguments:
                        if hasattr(arg, "type") and is_fieldAPI_ARRAY(arg.type.dtype.name):
                            new_args += (Variable(name='PTR', parent = arg, scope=routine),)
                            to_update = True
                        else:
                            new_args += (arg,)
                    new_kwargs = ()
                    for couple in abort_call.kwarguments:
                        if hasattr(couple[1], "type") and  is_fieldAPI_ARRAY(couple[1].type.dtype.name):
                            new_kwargs += ((couple[0], Variable(name='PTR', parent=couple[1], scope=routine)) ,)
                            to_update = True
                        else:
                            new_kwargs += ((couple[0], couple[1]),)

                    if to_update:
                        abort_call._update(arguments = new_args)
                        abort_call._update(kwarguments = new_kwargs)

                    abort_calls += (abort_call,)

            abort_calls = Section(abort_calls)
            transform = AddSuffixToCalls(suffix='_ABORT')
            abort_calls = transform.transform_node(abort_calls, routine)
            
            for index,target in enumerate(targets):

                call_sync = CallStatement(name = DeferredTypeSymbol(name=sync_names[1] if 'OpenACC' in target else sync_names[0]), 
                                                arguments=(), scope=None )
  
                #Create the replacement body for the region for current target

                if (target == 'OpenMP'):
                    new_loop = self.makeOpenMPLoop(routine, local_boundary_variable, region)
                    new_body = (call_sync, new_loop,)

                elif (target == 'OpenMPSingleColumn'):
                    new_region = self.makeParallelSubroutine(routine, region, region_num, local_boundary_variable, scc=True, target='host')
                    new_loop = self.makeOpenMPSCCLoop(routine, local_boundary_variable, new_region)
                    new_body = (call_sync, new_loop,)                    
                
                elif (target == 'OpenACCSingleColumn'):
                    new_region = self.makeParallelSubroutine(routine, region, region_num, local_boundary_variable, scc=True, target='device')
                    new_loop = self.makeOpenACCSCCLoop(routine, local_boundary_variable, new_region)
                    new_body = (call_sync, new_loop,)                    
                else :
                    print("Unknown target !!!! ", target)
                    new_body = (region.body,)


                #Create nested conditionals with call to LPARALLELMETHOD in conditions

                expr = InlineCall(function=DeferredTypeSymbol(name='LPARALLELMETHOD'), 
                        parameters=(StringLiteral(target.upper()),StringLiteral(f'{routine.name}_PARALLEL2:{name}'), ))
             
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
                                                                arguments=(StringLiteral(f'{routine.name}_PARALLEL2 : METHOD WAS NOT FOUND')),
                                                                scope=routine
                                                                ),
                                                )
                                    )

            routine.body = Transformer({region:(abort_calls,outer_cond)}).visit(routine.body)

            addFieldAPIPointers(routine, total_FAPI_pointers)

            routine.variables += (local_boundary_variable, local_stack,)   


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

