from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, 
    PragmaRegion )

from loki.ir import Section, Comment, VariableDeclaration, Pragma, FindVariables, SubstituteExpressions
from loki.ir.pragma_utils import PragmaRegionAttacher
from loki.transformations.build_system import DependencyTransformation
from loki.transformations.sanitise import do_resolve_associates

from loki.frontend.fparser import *

from fieldAPITransforms import FieldAPIPtr

from arpege_parameters import params

from termcolor import colored
from codetiming import Timer
from loki.logging import perf, info
import sys, traceback

sys.path.append('/home/legaux/meteo/models/build_scc_erwan/transformation/')

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
    import logical_lst
except ImportError as error:
    print (colored("build_scc : logical_lst not imported, check sys.path !", "red"))
    print("Error message : ", error)
    traceback.print_exc()
    exit(1)
else:
    print (colored("build_scc : logical_lst imported", "green"))


from syncTransforms import MakeSync
from parallelTransforms import MakeParallel
from commonTransforms import RemovePragmas, RemovePragmaRegions, AddSuffixToCalls, InlineMemberCalls, RemoveComments, RemoveLoops, RemoveEmptyConditionals, AddACCRoutineDirectives


acdc_logger = info # perf pour log fichier

for file in ['../compute/cpg_dia_flu.F90', '../compute/cpcfu.F90', '../compute/cpxfu.F90', '../compute/dprecips_xfu.F90', '../compute/meanwind_xfu.F90' ]:
# for file in ['../compute/cpxfu_simple.F90', ]:
# for file in ['../compute/cpxfu.F90' ]:
# for file in ['../compute/cpg_dia_flu.F90', ]:
# for file in ['../compute/dprecips_xfu.F90', ]:

    source = Sourcefile.from_file(file, frontend=Frontend.FP)

    routines = source.subroutines

    for routine in routines:
        to_generate = []
        start = None

        pragmas_map = {}

        # Search spec for GENERATE pragmas
        for pragma in FindNodes(Pragma).visit(routine.spec):
            if (pragma.keyword == 'ACDC'):
                if ('GENERATE' in pragma.content ):
                    splitted = pragma.content.split('TARGET')

                    # Get value after TARGET keyword and strip '=' character
                    splitted = splitted[-1].split('=')[-1]
                    splitted = splitted.split('/')
                    for s in splitted:
                        if s not in to_generate:
                            to_generate.append(s)
                pragmas_map[pragma] = None

        print("transformations to generate : ", to_generate)

        if 'Parallel' in to_generate:
        #search and bind pragmas for parallel regions

            pairs = []
            for pragma in FindNodes(Pragma).visit(routine.body):
                print("pragma found : ", pragma)
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

    for transform in to_generate:
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

        elif (transform == 'SingleColumnFieldAPIHost' or transform == 'SingleColumnFieldAPIDevice'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            with Timer(logger=acdc_logger, text=f'[ACDC] {transform} complete transformation' + ' in {:.2f}s'):
                isHost = (transform == 'SingleColumnFieldAPIHost')
                suffix='_SINGLE_COLUMN_FIELD_API_' + ('HOST' if isHost else 'DEVICE')
                filename = '../loki_outputs/' + file[10:-4] + '_single_column_field_api' + ('_host.F90' if isHost else '_device.F90')
                f = open(filename, 'w')
                for routine in routines:

                    new_routine = routine.clone()
                    new_routine.apply(RemovePragmaRegions())
                    new_routine.apply(RemovePragmas())
                    new_routine.apply(InlineMemberCalls())
                    new_routine.apply(RemoveComments())

                    true_symbols, false_symbols=logical_lst.symbols()
                    false_symbols.append('LHOOK')
                    
                    with Timer(logger=acdc_logger, text='[ACDC] scc_transform_routine in {:.2f}s'):
                        scc_transform_routine(new_routine, params.nproma_aliases, params.nproma_loop_indices, params.nproma_bounds, true_symbols, false_symbols)
                    
                    with Timer(logger=acdc_logger, text='[ACDC] FieldAPIPtr transform in {:.2f}s'):
                        new_routine.apply(FieldAPIPtr(pointerType='host' if isHost else 'device'))

                    # Change called subroutines names, import their interface and add !$acc routine directives if relevant
                    new_routine.apply(AddSuffixToCalls(suffix=suffix, additional_variables=['YLSTACK']))
                    if not isHost:
                        new_routine.spec.append(Pragma(keyword='acc', content='routine seq'))
                        new_routine.apply(AddACCRoutineDirectives())

                    # Add suffix and generate interface file
                    transfodep = DependencyTransformation(suffix=suffix,  include_path='../loki_outputs/')
                    new_routine.apply(transfodep, role='kernel')

                    # Turn local arrays into cray pointers on the stack through "alloc" and "temp" macros
                    alloc_temp(new_routine)

                
                    print("writing to file : ", filename)
                    f.write(new_routine.to_fortran())
                f.write('\n') #Add eol so vim doesn't complain
                f.close()


        elif (transform == 'SyncDevice' or transform == 'SyncHost'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            isHost = (transform == 'SyncHost')
            filename = '../loki_outputs/' + file[10:-4] + ('_sync_host.F90' if isHost else '_sync_device.F90')
             
            f = open(filename, 'w')
            for routine in routines:
                new_routine = routine.clone()

                
                # Cleanup
                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(RemovePragmas())
                new_routine.apply(InlineMemberCalls())
                new_routine.apply(RemoveComments())
                #new_routine.apply(RemoveLoops(['JLON', 'JLEV']))
                new_routine.apply(RemoveLoops())
                new_routine.apply(AddSuffixToCalls(suffix='_SYNC_HOST' if isHost else '_SYNC_DEVICE'))
                print("after add suffix")

                # Resolve associates because variables might be associated to a Field API member and we need to identify them
                do_resolve_associates(new_routine)
                # Main transformation : analyse dataflow and insert relevant sync calls
                new_routine.apply(MakeSync(pointerType='host' if isHost else 'device'))

                # new_routine.apply(RemoveEmptyConditionals())

                # Add suffix and generate interface file
                transfodep = DependencyTransformation(suffix='_SYNC_HOST' if isHost else '_SYNC_DEVICE',
                                                         # mode='strict' , 
                                                         include_path='../loki_outputs/')
                new_routine.apply(transfodep, role='kernel')
                
                print("writing to file : ", filename)
                f.write(new_routine.to_fortran())
            
            f.write('\n') #Add eol so vim doesn't complain
            f.close()




        elif (transform == 'Parallel'):
            print(f'call {transform} ')
            print("====_____________================______________================______________=================")
            filename = '../loki_outputs/' + file[10:-4] + '_parallel.F90'
            f = open(filename, 'w')
            for routine in routines:
                # new_routine = routine.clone(contains=None)
                new_routine = routine.clone()
                new_routine.apply(RemoveComments())
                new_routine.apply(InlineMemberCalls())
                print("after inline")
                do_resolve_associates(new_routine)

                parallel_transform = MakeParallel() 
                new_routine.apply(parallel_transform)

                new_routine.apply(RemovePragmaRegions())
                new_routine.apply(RemovePragmas())


                transfodep = DependencyTransformation(suffix='_PARALLEL', 
                                                         include_path='../loki_outputs/')
                new_routine.apply(transfodep, role='kernel')
                
                print("writing to file : ", filename)
                f.write(new_routine.to_fortran())

                for subroutine in parallel_transform.outlined_routines:
                    f.write('\n')
                    f.write(subroutine.to_fortran())
                new_routine = None
            f.write('\n') 
            f.close()

exit(0)



