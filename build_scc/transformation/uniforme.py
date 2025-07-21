# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.
class UniformizeLoops(Transformation):
  def __init__(self, horizontal, verbose=False):
    horizontal = horizontal
    verbose = verbose
  def transform_subroutine(self, routine, **kwargs):
    if verbose : print("transforming routine ", routine.name)
    if verbose : print(" ============================== ")
    total = 0
    if horizontal.index not in routine.variable_map:
      routine.variables += (Variable(name=horizontal.index,type=SymbolAttributes(BasicType.INTEGER), scope=routine), )
    loop_map = {}
    for loop in FindNodes(Loop).visit(routine.body):
      if loop.variable != horizontal.index:
        if (loop.bounds.upper == routine.variable_map[horizontal.bounds[1]] ) :
          to_replace = False
          if (loop.bounds.lower == routine.variable_map[horizontal.bounds[0]] ) :
            to_replace = True
          if ( isinstance(loop.bounds.lower, IntLiteral) ) :
            if (loop.bounds.lower.value == 1):
              to_replace = True
              print(colored(f'Assuming {loop.bounds.lower} is an alias for {routine.variable_map[horizontal.bounds[0]]} in {loop}', 'red' ))
          if (to_replace):
            if verbose : print (f'replacing loop counter : {loop.variable}')
            total += 1
            new_loop = loop.clone(variable=routine.variable_map[horizontal.index])
            var_map={}
            for var in FindVariables().visit(loop.body):
              if (var==loop.variable) :
                var_map[var] = routine.variable_map[horizontal.index]
            loop.body = SubstituteExpressions(var_map).visit(loop.body)
            loop.variable=routine.variable_map[horizontal.index]
    if (total > 0) : info(f'[Loki] {routine.name}:: {total} loops counters replaced by {horizontal.index}')

