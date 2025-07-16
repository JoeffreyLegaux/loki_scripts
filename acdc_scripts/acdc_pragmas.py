from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, 
    PragmaRegion, MultiConditional )

from loki.ir import Section, Comment, VariableDeclaration, Pragma, FindVariables, SubstituteExpressions
from loki.ir.pragma_utils import PragmaRegionAttacher
from loki.transformations.build_system import DependencyTransformation
from loki.transformations.sanitise import do_resolve_associates

from loki.frontend.fparser import *


from termcolor import colored
from codetiming import Timer
from loki.logging import perf, info
import sys, traceback

#===================================================================================================
# Importing SCC transforms from Erwan 
#===================================================================================================
#sys.path.append('/home/legaux/meteo/models/build_scc_erwan/transformation/')
sys.path.append('/home/ext/cf/ccom/legauxj/gpupack/build_scc/transformation/')

try :
    from openacc_transform import scc_transform_routine, alloc_temp, assoc_alloc_pt
except ImportError as error:
    print (colored("build_scc : scc_transform_routine not imported, check sys.path !", "red"))
    print("Error message : ", error)
    traceback.print_exc()
    exit(1)
else:
    print (colored("build_scc : scc_transform_routine imported", "green"))

try :
   import logical_lst, logical
except ImportError as error:
   print (colored("build_scc : logical_lst or logical not imported, check sys.path !", "red"))
   print("Error message : ", error)
   traceback.print_exc()
   exit(1)
else:
   print (colored("build_scc : logical_lst imported", "green"))
#===================================================================================================


#Importing local transformation classes
from fieldAPITransforms import FieldAPIPtr, get_pointers_to_FieldAPI, nproma_to_FieldAPI
from arpege_parameters import params
from syncTransforms import MakeSync
from parallelTransforms import MakeParallel
from commonTransforms import (SplitArraysDeclarations, RemovePragmas, RemovePragmaRegions, 
        RemoveAssignments, ReplaceAbortRegions, AddSuffixToCalls, AddFieldAPISuffixToCalls, InlineMemberCalls, 
        RemoveComments, RemoveLoops, RemoveEmptyConditionals, AddACCRoutineDirectives, RepairMultiConditionals )




def add_to_transforms(routine, transformations):
    global treated_routines, routines_to_transform
    if routine in treated_routines:
        transformations = {trans for trans in transformations if trans not in treated_routines[routine]}

    if routine not in routines_to_transform:
        routines_to_transform[routine] = transformations
    else:
        routines_to_transform[routine] = routines_to_transform[routine].union(transformations)

def add_to_treated(routine, transformations):
    global treated_routines
    if routine not in treated_routines:
        treated_routines[routine] = transformations
    else:
        treated_routines[routine] = treated_routines[routine].union(transformations)

    

def attach_acdc_regions(routine):
    pairs = []
    for pragma in FindNodes(Pragma).visit(routine.body):
        #print("pragma found : ", pragma)
        if (pragma.keyword == 'ACDC'):
            if ('{' in pragma.content):
                start = pragma
            elif ('}' in pragma.content):
                if (not start):
                    print(colored('ACDC } without previous { !!!!', "red"))
                    exit(0)
                else:
                    print(colored(f'ACDC pragma region identified : {start.content}, {pragma.content}', "green"))
                    # extract_pragma_region(routine.body, start=start, end=pragma)
                    pairs.append((start,pragma))
                    start = None
        else:
            print(colored(f'Unknown pragma found : {pragma}', "red"))

    PragmaRegionAttacher(pragma_pairs=pairs, inplace=True).visit(routine.body)


true_symbols, false_symbols=logical_lst.symbols()

acdc_logger = info # perf pour log fichier

#source_path = '../../loki_WIP/src/local/'
source_path = '../../loki_WIP/src/'

#output_path = source_path + '/local/'
output_path = source_path + 'local/'


output_path_sources = output_path + 'ifsaux/loki_sources/'
output_path_interfaces = output_path + 'ifsaux/loki_interfaces/'


# Start at CPG_DYN_SLG with empty list of forced transformations

routines_to_transform = {}
routines_to_transform['CPG_DYN_SLG']={'PARALLEL'}
routines_to_transform['LACDYN'] ={'PARALLEL'}
#routines_to_transform['LACDYN'] ={'SYNC_HOST', 'SYNC_DEVICE'}
#routines_to_transform['LASSIE']={'ABORT','PARALLEL'}
#routines_to_transform['LAVENT']={'ABORT','PARALLEL'}
#routines_to_transform['LAVABO']={'ABORT','PARALLEL'}
#routines_to_transform['LASURE'] = {'SCC_HOST','SCC_DEVICE'}
#routines_to_transform['LATTES']={'ABORT','PARALLEL'}
#routines_to_transform['LATTEX']={'ABORT','PARALLEL'}
#routines_to_transform['LATTEX']={'PARALLEL'}
#routines_to_transform = {'VERDISINT':{'ABORT'}}
treated_routines = {}


#routines_to_transform['LAVABO_EXPL_LAITVSPCQM_PART1']={'SYNC_HOST'}
#routines_to_transform['VERDISINT']={'ABORT'}

#routines_to_transform['LATTEX_EXPL_VSPLTRANS']={'SYNC_HOST'}
# If set to False, only apply transformation listed in routines_to_transform
# If set to True, enqueue subroutines called during a transformation for further transformation
greedy_process = True #False #True

while (len(routines_to_transform) > 0):
    routine_name =  next(iter(routines_to_transform))

    print(colored( '\n' + '+' * 100, "light_blue"))
    print(colored(f'+++    Treating routine {routine_name}' + ' ' * (100-(27+len(routine_name))) + '+++', "light_blue"))
    print(colored( '+' * 100 + '\n', "light_blue"))
    transformations_to_generate = routines_to_transform[routine_name]

    # Remove routine from the list. If new transformations arise, treat them in the next pass
    routines_to_transform.pop(routine_name)
    print(colored(f'Transformations to generate : {transformations_to_generate} \n', "cyan")) 
    #print("treated routines : ", treated_routines)    
    if routine_name not in params.routines_to_files:
        print(f'ERROR : routine {routine} not found in files list !!!')
        exit(0)
    else:
        file = params.routines_to_files[routine_name]

    file_without_path = file.split('/')[-1:][0]

    source = Sourcefile.from_file(source_path+file, frontend=Frontend.FP)

    routine = source.subroutines[0]
   
    if routine_name != routine_name:
         print(f'ERROR : routine {routine} from file {file} does not correspond to name {routine_name} !!!')
         exit(0)
    
    attach_acdc_regions(routine)
    

    for transform in transformations_to_generate:

        logical.transform_subroutine(routine, true_symbols, false_symbols) 
        
        (FieldAPI_pointers, FieldAPI_pointers_names) = get_pointers_to_FieldAPI(routine, params.nproma_aliases)

        if ( len(FieldAPI_pointers) > 0) :
            print("FieldAPI_pointers found in routine : ", FieldAPI_pointers, type(FieldAPI_pointers))
       
        
        print(colored(f'\n  =====   Applying transformation {transform}  ===== \n', "light_magenta"))
        if (transform == 'FIELD_API'):
            #filename = '../loki_outputs/' + file[10:-4] + '_field_api_host.F90'
            filename = output_path_sources + file_without_path[:-4] + '_field_api.F90'
            f = open(filename, 'w')
            if True:
            #for routine in routines:

                # Temporary remove contained routines for transformation
                new_routine = routine.clone(contains = None)


                new_routine.apply(ReplaceAbortRegions(abort_call = False))
                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(InlineMemberCalls())

                new_routine.apply(RemovePragmas())
                new_routine.apply(RemoveComments())

                new_routine.apply(FieldAPIPtr(pointerType='host'))

                

              
                add_suffix_transform = AddFieldAPISuffixToCalls(argument='YDCPG_BNDS')
                new_routine.apply(add_suffix_transform) 
            
                # Add suffix and generate interface file
            
                if greedy_process:
                    for subroutine in add_suffix_transform.routines_called:
                        add_to_transforms(subroutine, {'FIELD_API'})


                
                transfodep = DependencyTransformation(suffix='_FIELD_API', include_path=output_path_interfaces)
                new_routine.apply(transfodep, role='kernel')
 
                # Get contained routines back
                new_routine.contains = routine.contains
            
                print(colored(f'writing to file : {filename}', "light_green"))
                f.write(new_routine.to_fortran())
            f.write('\n') #Add eol so vim doesn't complain
            f.close()

        elif (transform == 'SCC_HOST' or transform == 'SCC_DEVICE'):
            with Timer(logger=acdc_logger, text=f'[ACDC] {transform} complete transformation' + ' in {:.2f}s'):
                isHost = (transform == 'SCC_HOST')
                #filename = output_path + file[:-4] + '_scc' + ('_host.F90' if isHost else '_device.F90')
                filename = output_path_sources + file_without_path[:-4] + '_scc' + ('_host.F90' if isHost else '_device.F90')
                f = open(filename, 'w')
                if True:
                #for routine in routines:

                    new_routine = routine.clone()
                    new_routine.apply(ReplaceAbortRegions(abort_call = False))
                    new_routine.apply(RemovePragmaRegions())
                    new_routine.apply(RemovePragmas())
                    new_routine.apply(InlineMemberCalls())
                    new_routine.apply(RemoveComments())

                    true_symbols, false_symbols=logical_lst.symbols()
                    false_symbols.append('LHOOK')
                    
                    with Timer(logger=acdc_logger, text='[ACDC] scc_transform_routine in {:.2f}s'):
                        scc_transform_routine(new_routine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols, FieldAPI_pointers=FieldAPI_pointers)
                    
                    # Change called subroutines names, import their interface and add !$acc routine directives if relevant
                    add_suffix_transform =(AddSuffixToCalls(suffix='_'+transform, additional_kwvariables=[('YDSTACK','YLSTACK')]))
                    new_routine.apply(add_suffix_transform)
                    
                    if greedy_process:
                        for subroutine in add_suffix_transform.routines_called:
                            add_to_transforms(subroutine, {transform})

                    if not isHost:
                        new_routine.spec.append(Pragma(keyword='acc', content='routine seq'))
                        new_routine.apply(AddACCRoutineDirectives())

                    # Add suffix and generate interface file
                    transfodep = DependencyTransformation(suffix='_'+transform,  include_path=output_path_interfaces)
                    new_routine.apply(transfodep, role='kernel')

                    # Turn local arrays into cray pointers on the stack through "alloc" and "temp" macros
                    alloc_temp(new_routine)
                    assoc_alloc_pt(new_routine, FieldAPI_pointers_names)
                
                    print(colored(f'writing to file : {filename}', "light_green"))
                    f.write(new_routine.to_fortran())

                    
                    add_to_treated(routine_name, {transform})
                f.write('\n') #Add eol so vim doesn't complain
                f.close()


        elif (transform == 'SYNC_DEVICE' or transform == 'SYNC_HOST'):
            isHost = (transform == 'SYNC_HOST')
            #filename = output_path + file[:-4] + ('_sync_host.F90' if isHost else '_sync_device.F90')
            filename = output_path_sources + file_without_path[:-4] + ('_sync_host.F90' if isHost else '_sync_device.F90')
             
            f = open(filename, 'w')
            if True:
            #for routine in routines:
                new_routine = routine.clone()
                new_routine.apply(ReplaceAbortRegions(abort_call = False))
                
                # Cleanup
                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(RemovePragmas())
                new_routine.apply(InlineMemberCalls())
                new_routine.apply(RemoveComments())
                #new_routine.apply(RemoveLoops(['JLON', 'JLEV']))
                new_routine.apply(RemoveLoops(indices=params.nproma_loop_indices + params.vertical_loop_indices))
                add_suffix_transform = AddSuffixToCalls(suffix='_SYNC_HOST' if isHost else '_SYNC_DEVICE')
                new_routine.apply(add_suffix_transform)
                
                if greedy_process:
                    for subroutine in add_suffix_transform.routines_called:
                        add_to_transforms(subroutine, {transform})

                # Resolve associates because variables might be associated to a Field API member and we need to identify them
                do_resolve_associates(new_routine)
                # Main transformation : analyse dataflow and insert relevant sync calls
                new_routine.apply(MakeSync(pointerType='host' if isHost else 'device', nproma_pointers=FieldAPI_pointers))

                # new_routine.apply(RemoveEmptyConditionals())

                # Add suffix and generate interface file
                transfodep = DependencyTransformation(suffix='_SYNC_HOST' if isHost else '_SYNC_DEVICE',
                                                         # mode='strict' , 
                                                         include_path=output_path_interfaces)
                new_routine.apply(transfodep, role='kernel')
                

                print(colored(f'writing to file : {filename}', "light_green"))
                f.write(new_routine.to_fortran())
            
                add_to_treated(routine_name, {transform})

            f.write('\n') #Add eol so vim doesn't complain
            f.close()


        elif (transform == 'ABORT') : #and routine not in params.ignore_abort):
            #filename = output_path + file[:-4] + '_abort.F90'
            filename = output_path_sources + file_without_path[:-4] + '_abort.F90'
            #print("ouput file opened : ", filename)
            f = open(filename, 'w')
            if True:
            #for routine in routines:
                new_routine = routine.clone()
                do_resolve_associates(new_routine)
                new_routine.apply(RemoveAssignments())                
                new_routine.apply(RemoveLoops())

               
                new_routine.apply(ReplaceAbortRegions(abort_call=True))
             
                new_routine.apply(RemovePragmaRegions(empty = True))
                add_suffix_transform = AddSuffixToCalls(suffix='_ABORT')
                new_routine.apply(add_suffix_transform)
            
                new_routine.apply(RemovePragmas())
                new_routine.apply(RemoveComments())

                new_routine.apply(RemoveEmptyConditionals())
                # Multi conditionals with emptied bodies get empty value instead of tuple of empties
                # Probably a bug ?
                new_routine.apply(RepairMultiConditionals())

                if greedy_process:
                    for subroutine in add_suffix_transform.routines_called:
                        add_to_transforms(subroutine, {'ABORT'})


                transfodep = DependencyTransformation(suffix='_ABORT', include_path=output_path_interfaces)
                new_routine.apply(transfodep, role='kernel')

                print(colored(f'writing to file : {filename}', "light_green"))
                f.write(new_routine.to_fortran())

                add_to_treated(routine_name, {'ABORT'})

            f.write('\n') 
            f.close()


        elif (transform == 'PARALLEL'):
            #filename = output_path + file[:-4] + '_parallel.F90'
            filename = output_path_sources + file_without_path[:-4] + '_parallel2.F90'
            #print("ouput file opened : ", filename)
            f = open(filename, 'w')
            if True:
            #for routine in routines:
                new_routine = routine.clone()
                #new_routine.apply(RemoveComments())
                #new_routine.apply(InlineMemberCalls())
                do_resolve_associates(new_routine)

                new_routine.apply(ReplaceAbortRegions())

                parallel_transform = MakeParallel(FieldAPI_pointers) 
                new_routine.apply(parallel_transform)
                #print("routine parallele : ", parallel_transform.subroutines_to_transform)
                # We probably will have for new routines 'PARALLEL' 'SYNC' and 'ABORT' transforms to add
                #print("treated routines : ", treated_routines)
                
                if greedy_process:
                    for subroutine in parallel_transform.subroutines_to_transform:
                        add_to_transforms(subroutine, parallel_transform.subroutines_to_transform[subroutine])
                        #print("added transfo ", subroutine, parallel_transform.subroutines_to_transform[subroutine])

                #print("routines to trasnsform ? ", routines_to_transform)

                new_routine.apply(RemovePragmaRegions())

                transfodep = DependencyTransformation(suffix='_PARALLEL2', include_path=output_path_interfaces)
                new_routine.apply(transfodep, role='kernel')
                
                print(colored(f'writing to file : {filename}', "light_green"))
                f.write(new_routine.to_fortran())


                for subroutine in parallel_transform.outlined_routines:
                    f.write('\n')
                    f.write(subroutine.to_fortran())
                new_routine = None



                add_to_treated(routine_name, {'PARALLEL'})
            f.write('\n') 
            f.close()

exit(0)



