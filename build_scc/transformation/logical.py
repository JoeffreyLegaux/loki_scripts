# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.

from loki import (
    Sourcefile, FindNodes, CallStatement, 
    Transformer, Dimension, ir, 
    Scalar, Assignment, fgen,
    FindVariables, symbols, demote_variables,
    Intrinsic, Variable, SymbolAttributes,
    DerivedType, Conditional, FindTypedSymbols,
    DeferredTypeSymbol, LogicalNot, LogicalOr,
    LogicalAnd,
)

from loki.expression.symbols import Comparison
from loki.expression.symbols import IntLiteral

import os

def evaluate(comparison):
  """
  Evaluates comparison when the operators are integer values.
  :param comparison:.
  """
  operator=comparison.operator
  right=comparison.right
  left=comparison.left
  if operator=='==':
    return(left==right)
  elif operator=='<':
    return(left<right)
  elif operator=='>':
    return(left>right)
  elif operator=='!=':
    return(left!=right)
  elif operator=='>=':
    return(left>=right)
  elif operator=='<=':
    return(left<=right)
        
    
def evaluateCondition(my_condition, true_symbols, false_symbols):

  if isinstance(my_condition, DeferredTypeSymbol) or isinstance(my_condition, Scalar):
    # Simple terminal symboles : either True, False or not defined
    if (my_condition.name in true_symbols):
      return[True, True]
    elif (my_condition.name in false_symbols):
      return[True, False]
    else :
      return[False, False]

  elif isinstance(my_condition, Comparison):
    if isinstance(my_condition.left, IntLiteral) and isinstance(my_condition.right, IntLiteral):
      if evaluate(my_condition):
        return[True, True]
      else:
        return[True,False]
    else:
      return[False, False]
    
  elif isinstance(my_condition, LogicalNot):
    [evaluable, evaluation] = evaluateCondition(my_condition.child, true_symbols, false_symbols)
    return [evaluable, not evaluation]

  elif isinstance(my_condition, LogicalOr):
    # Evaluates to True if a child evaluates to True
    # Evaluates to False if both children evaluate to False
    # Cannot be evaluated otherwise
    maybe_false = True
    for c in my_condition.children :
      evaluable, evaluation = evaluateCondition(c, true_symbols, false_symbols)
      if (evaluable and evaluation ) :
        return[True, True]
      elif (not evaluable) :
        maybe_false = False

    # If all children evaluated and we did not return, they were false => evaluate OR as false
    if (maybe_false) :
      return [True, False]
    else :
      return [False, False]

  elif isinstance(my_condition, LogicalAnd):
    # Evaluates to True if both children evaluate to True
    # Evaluates to False if a child evaluates to False
    # Cannot be evaluated otherwise
    maybe_true = True
    for c in my_condition.children :
      evaluable, evaluation = evaluateCondition(c, true_symbols, false_symbols)
      if (evaluable and not evaluation) : 
        return [True, False]
      elif (not evaluable) :
        maybe_true = False

    # If all children evaluated and we did not return, they were true => evaluate AND as true
    if (maybe_true) :
      return [True, True]
    else :
      return [False, False]

  else :
    return[False, False]



def transform_subroutine(routine, true_symbols, false_symbols):
    callmap={}
    for cond in FindNodes(Conditional).visit(routine.body):

      # If at least one of the forced-value symbols is present in the condition, we can try to evaluate it
        evaluable = False
        for symbol in FindTypedSymbols().visit(cond.condition) :
            if (symbol in true_symbols or symbol in false_symbols):
                evaluable = True

        evaluable = True
        if evaluable :
          [evaluable, evaluation] = evaluateCondition(cond.condition, true_symbols, false_symbols)
          if evaluable :
            if evaluation :
              callmap[cond] = cond.body
            else:
              callmap[cond] = cond.else_body
            

    routine.body=Transformer(callmap).visit(routine.body)
