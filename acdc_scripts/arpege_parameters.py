from dataclasses import dataclass

@dataclass
class params_class:
  nproma_aliases = ["KLON","YDCPG_OPTS%KLON","YDGEOMETRY%YRDIM%NPROMA","KPROMA", "YDDIM%NPROMA", "NPROMA"]
  nproma_bounds = [["KIDIA", "YDCPG_BNDS%KIDIA","KST"],["KFDIA", "YDCPG_BNDS%KFDIA","KEND"]]
  nproma_loop_indices = ['JLON','JROF']
  cpg_opts_variable = 'YDCPG_OPTS'
  block_dimension = 'YDCPG_OPTS%KGPBLKS'
  boundaries_type = 'CPG_BNDS_TYPE'

params = params_class()