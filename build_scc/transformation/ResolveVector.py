# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.

from loki.expression import symbols
from loki import ir
from loki import (
    Transformation, FindNodes, FindScopes, FindVariables,
    FindExpressions, Transformer, NestedTransformer,
    SubstituteExpressions, SymbolAttributes, BasicType, DerivedType,
    pragmas_attached, CaseInsensitiveDict, as_tuple, flatten,
    demote_variables
)



def resolve_vector_dimension(routine, loop_variable, bounds):
    """
    Resolve vector notation for a given dimension only. The dimension
    is defined by a loop variable and the bounds of the given range.

    TODO: Consolidate this with the internal
    `loki.transform.transform_array_indexing.resolve_vector_notation`.

    Parameters
    ----------
    routine : :any:`Subroutine`
        The subroutine in which to resolve vector notation usage.
    loop_variable : :any:`Scalar`
        The induction variable for the created loops.
    bounds : tuple of :any:`Scalar`
        Tuple defining the iteration space of the inserted loops.
    """
    bounds_str = f'{bounds[0]}:{bounds[1]}'

    bounds_v = (symbols.Variable(name=bounds[0].name), symbols.Variable(name=bounds[1].name))

    mapper = {}
    for stmt in FindNodes(ir.Assignment).visit(routine.body):
        ranges = [e for e in FindExpressions().visit(stmt)
                  if isinstance(e, symbols.RangeIndex) and e == bounds_str]
        if ranges:
            exprmap = {r: loop_variable for r in ranges}
            loop = ir.Loop(variable=loop_variable, bounds=symbols.LoopRange(bounds_v),
                           body=(SubstituteExpressions(exprmap).visit(stmt),) )
            mapper[stmt] = loop

    routine.body = Transformer(mapper).visit(routine.body)


