from dataclasses import dataclass

@dataclass
class params_class:
  nproma_aliases = ["KLON","YDCPG_OPTS%KLON","YDGEOMETRY%YRDIM%NPROMNH","YDGEOMETRY%YRDIM%NPROMA","KPROMA", "YDDIM%NPROMA", "NPROMA"]
  nproma_bounds = [["KIDIA", "YDCPG_BNDS%KIDIA","KST"],["KFDIA", "YDCPG_BNDS%KFDIA","KEND"]]
  nproma_loop_indices = ['JLON','JROF']
  cpg_opts_variable = 'YDCPG_OPTS'
  block_dimension = 'YDCPG_OPTS%KGPBLKS'
  block_counter = 'JBLK'
  boundaries_type = 'CPG_BNDS_TYPE'
  ignored_subroutines = ['ABOR1','GETENV','DGEMM']
  routines_to_files = { 'CPG_DYN_SLG':'main/arpifs/adiab/cpg_dyn_slg.F90', 
                        #'LACDYN':'arpifs/adiab/lacdyn.F90',
                        'LACDYN':'local/arpifs/adiab/lacdyn.F90',
                        'LAVABO':'main/arpifs/adiab/lavabo.F90',
                        'LATTEX':'main/arpifs/adiab/lattex.F90',
                        'LATTES':'main/arpifs/adiab/lattes.F90',
                        'LAVENT':'main/arpifs/adiab/lavent.F90',
                        'LASSIE':'main/arpifs/adiab/lassie.F90',
                        'LASURE':'main/arpifs/adiab/lasure.F90',
                        'GPRCP_EXPL':'main/.fypp/arpifs/adiab/gprcp_expl.F90',
                        'SIGAM_GP':'main/.fypp/arpifs/adiab/sigam_gp.F90',
                        'SITNU_GP':'main/.fypp/arpifs/adiab/sitnu_gp.F90',
                        'LAVABO_EXPL_LAITVSPCQM_PART1':'main/arpifs/adiab/lavabo_expl_laitvspcqm_part1.F90',
                        'LAVABO_EXPL_LAITVSPCQM_PART2':'main/arpifs/adiab/lavabo_expl_laitvspcqm_part2.F90',
                        'LATTEX_EXPL_2TL':'main/arpifs/adiab/lattex_expl_2tl.F90',
                        'LATTEX_EXPL_VSPLTRANS':'main/arpifs/adiab/lattex_expl_vspltrans.F90',
                        'VERDISINT':'main/arpifs/utility/verdisint.F90',
                        'VERINT':'main/arpifs/utility/verint.F90',
                        'VERINTS':'main/arpifs/utility/verints.F90',

                        'VERDER':'main/arpifs/utility/verder.F90'
                        }

params = params_class()
