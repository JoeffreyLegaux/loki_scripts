from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, 
    PragmaRegion )

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
    from openacc_transform import scc_transform_routine , alloc_temp
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
from fieldAPITransforms import FieldAPIPtr, get_pointers_to_FieldAPI
from arpege_parameters import params
from syncTransforms import MakeSync
from parallelTransforms import MakeParallel
from commonTransforms import (RemovePragmas, RemovePragmaRegions, RemoveAssignments, ReplaceAbortRegions, 
    AddSuffixToCalls, InlineMemberCalls, RemoveComments, RemoveLoops, RemoveEmptyConditionals, AddACCRoutineDirectives )




def add_to_transforms(routine, transformations):
    if routine in treated_routines:
        transformations = {trans for trans in transformation if trans not in treated_routines[routine]}

    if routine not in routines_to_transform:
        routines_to_transform[routine] = transformations
    else:
        routines_to_transform[routine].union(transformations)

def add_to_treated(routine, transformations):
    if routine not in treated_routines:
        treated_routines[routine] = transformations
    else:
        treated_routines[routine].union(transformations)

    

def attach_acdc_regions(routine):
    pairs = []
    for pragma in FindNodes(Pragma).visit(routine.body):
        #print("pragma found : ", pragma)
        if (pragma.keyword == 'ACDC'):
            if ('{' in pragma.content):
                start = pragma
            elif ('}' in pragma.content):
                if (not start):
                    print("ACDC } without previous {")
                else:
                    print("ACDC pragma region identified")
                    # extract_pragma_region(routine.body, start=start, end=pragma)
                    pairs.append((start,pragma))
                    start = None
        else:
            print('Unknown pragma found : {pragma}')

    PragmaRegionAttacher(pragma_pairs=pairs, inplace=True).visit(routine.body)


true_symbols, false_symbols=logical_lst.symbols()

acdc_logger = info # perf pour log fichier

#source_path = '../../loki_WIP/src/local/'
source_path = '../../loki_WIP/src/'


output_path = source_path + '/local/'
output_path_scc = output_path
# Start at CPG_DYN_SLG with empty list of forced transformations

#routines_to_transform = {'CPG_DYN_SLG':{'PARALLEL'},}
#routines_to_transform = {'LACDYN':{'PARALLEL', },}  # 'ABORT'},}
routines_to_transform = {'VERINT':{'ABORT','SYNC_DEVICE','SCC_DEVICE'},}
#routines_to_transform = {'VERDISINT':{'SCC_DEVICE'},}
#routines_to_transform = {'SIGAM_GP':{'ABORT'},}
#routines_to_transform = {'LASSIE':{'ABORT','PARALLEL'},}
#routines_to_transform = {'GPRCP_EXPL':{'PARALLEL'},}

treated_routines = {}

while (len(routines_to_transform) > 0):
    routine =  next(iter(routines_to_transform))
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("treating ", routine)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    transformations_to_generate = routines_to_transform[routine]

    # Remove routine from the list. If new transformations arise, treat them in the next pass
    routines_to_transform.pop(routine)
    print("tranforms to generate : ", transformations_to_generate) 
    
    if routine not in params.routines_to_files:
        print(f'ERROR : routine {routine} not found in files list !!!')
        exit(0)
    else:
        file = params.routines_to_files[routine]
    
    print("file : ", file)

    source = Sourcefile.from_file(source_path+file, frontend=Frontend.FP)

    routines = source.subroutines

    for routine in routines:
        # start = None

        # pragmas_map = {}

        attach_acdc_regions(routine)
    

    for transform in transformations_to_generate:

        logical.transform_subroutine(routine, true_symbols, false_symbols) 
        
        FieldAPI_pointers = get_pointers_to_FieldAPI(routine, params.nproma_aliases)

        #print("FieldAPI_pointers found in acdc_pragmes : ", FieldAPI_pointers)
       

        # ========================================
        # Completely useless now ????
        # ========================================
        if (transform == 'FieldAPIHost'):
            filename = '../loki_outputs/' + file[10:-4] + '_field_api_host.F90'
            f = open(filename, 'w')
            for routine in routines:

                # Temporary remove contained routines for transformation
                new_routine = routine.clone(contains = None)

                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(RemovePragmas())
                new_routine.apply(RemoveComments())

                new_routine.apply(FieldAPIPtr(pointerType='host'))

                # Add suffix and generate interface file
                transfodep = DependencyTransformation(suffix='_FIELD_API_HOST',  include_path='../loki_outputs/')
                new_routine.apply(transfodep, role='kernel')
                
                # Get contained routines back
                new_routine.contains = routine.contains
                new_routine.apply(AddSuffixToCalls(suffix='_FIELD_API_HOST'))
            
                print("writing to file : ", filename)
                f.write(new_routine.to_fortran())
            f.write('\n') #Add eol so vim doesn't complain
            f.close()

        elif (transform == 'SCC_HOST' or transform == 'SCC_DEVICE'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            with Timer(logger=acdc_logger, text=f'[ACDC] {transform} complete transformation' + ' in {:.2f}s'):
                isHost = (transform == 'SCC_HOST')
                # suffix='_SINGLE_COLUMN_FIELD_API_' + ('HOST' if isHost else 'DEVICE')
                #filename = output_path + file[:-4] + '_scc' + ('_host.F90' if isHost else '_device.F90')
                filename = (source_path + file[:-4]).replace('main', 'local') + '_scc' + ('_host.F90' if isHost else '_device.F90')
                f = open(filename, 'w')
                for routine in routines:

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
                    
                    with Timer(logger=acdc_logger, text='[ACDC] FieldAPIPtr transform in {:.2f}s'):
                        new_routine.apply(FieldAPIPtr(pointerType='host' if isHost else 'device'))

                    # Change called subroutines names, import their interface and add !$acc routine directives if relevant
                    add_suffix_transform =(AddSuffixToCalls(suffix='_'+transform, additional_variables=['YLSTACK']))
                    new_routine.apply(add_suffix_transform)
                    for subroutine in add_suffix_transform.routines_called:
                        add_to_transforms(subroutine, {transform})

                    if not isHost:
                        new_routine.spec.append(Pragma(keyword='acc', content='routine seq'))
                        new_routine.apply(AddACCRoutineDirectives())

                    # Add suffix and generate interface file
                    transfodep = DependencyTransformation(suffix='_'+transform,  include_path='./')
                    new_routine.apply(transfodep, role='kernel')

                    # Turn local arrays into cray pointers on the stack through "alloc" and "temp" macros
                    alloc_temp(new_routine)

                
                    print("writing to file : ", filename)
                    f.write(new_routine.to_fortran())

                f.write('\n') #Add eol so vim doesn't complain
                f.close()


        elif (transform == 'SYNC_DEVICE' or transform == 'SYNC_HOST'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            isHost = (transform == 'SYNC_HOST')
            #filename = output_path + file[:-4] + ('_sync_host.F90' if isHost else '_sync_device.F90')
            filename = (source_path + file[:-4]).replace('main', 'local') + ('_sync_host.F90' if isHost else '_sync_device.F90')
             
            f = open(filename, 'w')
            for routine in routines:
                new_routine = routine.clone()
                new_routine.apply(ReplaceAbortRegions(abort_call = False))
                
                # Cleanup
                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(RemovePragmas())
                new_routine.apply(InlineMemberCalls())
                new_routine.apply(RemoveComments())
                #new_routine.apply(RemoveLoops(['JLON', 'JLEV']))
                new_routine.apply(RemoveLoops())
                add_suffix_transform = AddSuffixToCalls(suffix='_SYNC_HOST' if isHost else '_SYNC_DEVICE')
                new_routine.apply(add_suffix_transform)
                print("after add suffix")
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
                                                         include_path='./')
                new_routine.apply(transfodep, role='kernel')
                
                print("writing to file : ", filename)
                f.write(new_routine.to_fortran())
            
                add_to_treated(routine, {transform})

            f.write('\n') #Add eol so vim doesn't complain
            f.close()


        elif (transform == 'ABORT'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            #filename = output_path + file[:-4] + '_abort.F90'
            print("source path ", source_path)
            filename = (source_path + file[:-4]).replace('main', 'local') + '_abort.F90'
            print("ouput file opened : ", filename)
            f = open(filename, 'w')
            for routine in routines:
                new_routine = routine.clone()
                do_resolve_associates(new_routine)
                new_routine.apply(RemoveLoops())
                new_routine.apply(RemoveAssignments())                
                new_routine.apply(ReplaceAbortRegions())
                
                add_suffix_transform = AddSuffixToCalls(suffix='_ABORT')
                new_routine.apply(add_suffix_transform)
            

                new_routine.apply(RemovePragmas())
                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(RemoveComments())
                new_routine.apply(RemoveEmptyConditionals())

                for subroutine in add_suffix_transform.routines_called:
                    add_to_transforms(subroutine, {'ABORT'})


                transfodep = DependencyTransformation(suffix='_ABORT', include_path='./')
                new_routine.apply(transfodep, role='kernel')
                
                print("writing to file : ", filename)
                f.write(new_routine.to_fortran())

                add_to_treated(routine, {'ABORT'})

            f.write('\n') 
            f.close()


        elif (transform == 'PARALLEL'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            #filename = output_path + file[:-4] + '_parallel.F90'
            filename = (source_path + file[:-4]).replace('main', 'local') + '_parallel.F90'
            print("ouput file opened : ", filename)
            f = open(filename, 'w')
            for routine in routines:
                # new_routine = routine.clone(contains=None)
                new_routine = routine.clone()
                #new_routine.apply(RemoveComments())
                #new_routine.apply(InlineMemberCalls())
                print("after inline")
                do_resolve_associates(new_routine)

                new_routine.apply(ReplaceAbortRegions())

                parallel_transform = MakeParallel(FieldAPI_pointers) 
                new_routine.apply(parallel_transform)
                print("routine parallele : ", parallel_transform.subroutines_to_transform)
                # We probably will have for new routines 'PARALLEL' 'SYNC' and 'ABORT' transforms to add
                for subroutine in parallel_transform.subroutines_to_transform:
                    add_to_transforms(subroutine, parallel_transform.subroutines_to_transform[subroutine])
                    
                new_routine.apply(RemovePragmaRegions())

                transfodep = DependencyTransformation(suffix='_PARALLEL', include_path='./')
                new_routine.apply(transfodep, role='kernel')
                
                print("writing to file : ", filename)
                f.write(new_routine.to_fortran())


                for subroutine in parallel_transform.outlined_routines:
                    f.write('\n')
                    f.write(subroutine.to_fortran())
                new_routine = None



                add_to_treated(routine, {'PARALLEL'})

            f.write('\n') 
            f.close()

exit(0)



