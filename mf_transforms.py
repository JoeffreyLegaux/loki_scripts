# from loki import FP, Sourcefile, Dimension, Subroutine
from loki import (Frontend, Sourcefile, Scheduler, FindNodes, Loop, Variable,
            Assignment, CallStatement, Transformation, Node, SymbolAttributes, 
            DerivedType, BasicType, Import, Transformer, Conditional, SchedulerConfig
            )
from loki.expression import LoopRange, FindVariables, FindTypedSymbols
from loki.expression.expr_visitors import SubstituteExpressions
from loki.expression.symbols import (
    Array, Scalar, InlineCall, TypedSymbol, FloatLiteral, IntLiteral, LogicLiteral,
    StringLiteral, IntrinsicLiteral, DeferredTypeSymbol, LogicalOr, LogicalAnd, LogicalNot, RangeIndex
)
from loki.ir import Section, Comment, VariableDeclaration, Import, Intrinsic
from pathlib import Path
from termcolor import colored
from loki.logging import info
import sys


class CleanSpec(Transformation):
  def transform_subroutine(self, routine, **kwargs):
    spec_map = {}
    for intr in FindNodes(Intrinsic).visit(routine.spec):
      if 'POINTER' in intr.text:
        spec_map[intr] = None # Remove Cray pointers declarations
    for incl in FindNodes(Import).visit(routine.spec):
      if (incl.c_import == True):
        spec_map[incl] = None
    routine.spec = Transformer(spec_map).visit(routine.spec)
    
class UniformizeLoops(Transformation):
  def __init__(self, horizontal, verbose=False):
    self.horizontal = horizontal
    self.verbose = verbose
  def transform_subroutine(self, routine, **kwargs):
    if self.verbose : print("transforming routine ", routine.name)
    if self.verbose : print(" ============================== ")
    total = 0
    if self.horizontal.index not in routine.variable_map:
      routine.variables += (Variable(name=self.horizontal.index,type=SymbolAttributes(BasicType.INTEGER), scope=routine), )
    loop_map = {}
    for loop in FindNodes(Loop).visit(routine.body):
      if loop.variable != self.horizontal.index:
        if (loop.bounds.upper == routine.variable_map[self.horizontal.bounds[1]] ) :
          to_replace = False
          if (loop.bounds.lower == routine.variable_map[self.horizontal.bounds[0]] ) :
            to_replace = True
          if ( isinstance(loop.bounds.lower, IntLiteral) ) :
            if (loop.bounds.lower.value == 1):
              to_replace = True
              print(colored(f'Assuming {loop.bounds.lower} is an alias for {routine.variable_map[self.horizontal.bounds[0]]} in {loop}', 'red' ))
          if (to_replace):
            if self.verbose : print (f'replacing loop counter : {loop.variable}')
            total += 1
            new_loop = loop.clone(variable=routine.variable_map[self.horizontal.index])
            var_map={}
            for var in FindVariables().visit(loop.body):
              if (var==loop.variable) :
                var_map[var] = routine.variable_map[self.horizontal.index]
            loop.body = SubstituteExpressions(var_map).visit(loop.body)
            loop.variable=routine.variable_map[self.horizontal.index]
    if (total > 0) : info(f'[Loki] {routine.name}:: {total} loops counters replaced by {self.horizontal.index}')

class DerivedTypeDimensions(Transformation):
  def __init__(self, horizontal, vertical, block_dim, verbose=False):
    self.routine_map = {}
    self.var_list = []
    self.verbose = verbose
    self.first_pass = True
    
    if '%' in horizontal.size      : self.var_list.append(horizontal.size)
    if len(horizontal.bounds) > 1 :
      if '%' in horizontal.bounds[0] : self.var_list.append(horizontal.bounds[0])
      if '%' in horizontal.bounds[1] : self.var_list.append(horizontal.bounds[1])
    
    if '%' in vertical.size      : self.var_list.append(vertical.size)
    if len(vertical.bounds) > 1 :
      if '%' in vertical.bounds[0] : self.var_list.append(vertical.bounds[0])
      if '%' in vertical.bounds[1] : self.var_list.append(vertical.bounds[1])
    
    if '%' in block_dim.size      : self.var_list.append(block_dim.size)
    if len(block_dim.bounds) > 1 :
      if '%' in block_dim.bounds[0] : self.var_list.append(block_dim.bounds[0])
      if '%' in block_dim.bounds[1] : self.var_list.append(block_dim.bounds[1])


  def transform_subroutine(self, routine, **kwargs):
    if self.verbose : print("transforming routine ", routine.name)
    if self.verbose : print(" ============================== ")
    
    if self.first_pass :
      self.routine_map[routine.name] = []
      to_add=[]
      for var in self.var_list :
        splitted = var.split('%')
        if splitted[0] in routine.variable_map :
          if self.verbose : print(f'Derived Type {splitted[0]} found in routine')
          var = Variable(name=splitted[1], parent=routine.variable_map[splitted[0]] )
          self.routine_map[routine.name].append(var)
          routine.variables += (var,)

        else :
          if self.verbose : print(f'Derived Type {splitted[0]} not found in routine !')
          if self.verbose : print(f'  member {var} not added')
    else:
      if self.verbose : print(f'Second pass, {len(routine_map)} variables to remove.')
      if (len(self.routine_map) > 0) :
        cleaned_variables = filter(lambda x : x not in self.routine_map[routine.name], routine.variables)
        routine.variables = tuple(cleaned_variables)




# returns [evaluable, evaluation] 
def evaluateCondition(my_condition, true_symbols, false_symbols):

  if isinstance(my_condition, DeferredTypeSymbol) or isinstance(my_condition, Scalar):
    # Simple terminal symboles : either True, False or not defined
    if (my_condition.name in true_symbols):
      return[True, True]
    elif (my_condition.name in false_symbols):
      return[True, False]
    else :
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
    if self.verbose : print("untreated class :", my_condition.__class__)
    return[False, False]



class LogicalsPreproc(Transformation):

  def __init__(self, true_symbols=[], false_symbols=[], verbose = False):
        self.true_symbols = true_symbols
        self.false_symbols = false_symbols
        self.verbose = verbose

  def transform_subroutine(self, routine, **kwargs):
    if self.verbose : print ("Treating routine ", routine.name)
    callmap={}
    for cond in FindNodes(Conditional).visit(routine.body):

      # If at least one of the forced-value symbols is present in the condition, we can try to evaluate it
      evaluable = False
      for symbol in FindTypedSymbols().visit(cond.condition) :
        if (symbol in self.true_symbols or symbol in self.false_symbols):
          evaluable = True

      if evaluable :
        [evaluable, evaluation] = evaluateCondition(cond.condition, self.true_symbols, self.false_symbols)
        if self.verbose : print ("condition : ", cond.condition)
        if self.verbose : print("evaluated : ", evaluable, "  ---  value : ", evaluation)
        if evaluable :
          if evaluation :
            callmap[cond] = cond.body
          else:
            callmap[cond] = cond.else_body
    routine.body=Transformer(callmap).visit(routine.body)


class ExplicitArraySyntaxes(Transformation):
  def __init__(self, horizontal, verbose = False):
    self.horizontal = horizontal
    self.verbose = verbose
  def transform_subroutine(self, routine, **kwargs):
    role = kwargs['role']
    if role == 'kernel':
      self.process_kernel(routine)

  def process_kernel(self, routine):
    if self.verbose : print("transforming routine ", routine.name)
    if self.verbose : print(" ============================== ")
    total = 0
    assign_map={}
    
    begin_index = None
    end_index = None

    if '%' in self.horizontal.bounds[0]:
      splitted = self.horizontal.bounds[0].split('%')
      for var in FindVariables().visit(routine.variables):
        if (var.name == splitted[0]):
          begin_index=Variable(name=f'{var.name}%{splitted[1]}', parent=var, scope=routine)
          if self.verbose : print(colored("derived type " + splitted[0] + " found", "green"))
    else:
      for var in FindVariables().visit(routine.variables):
        if (var.name == self.horizontal.bounds[0]):
          begin_index = var
          if self.verbose : print(colored("variable " + var.name + " found", "green"))

    if '%' in self.horizontal.bounds[1]:
      splitted = self.horizontal.bounds[1].split('%')
      for var in FindVariables().visit(routine.variables):
        if (var.name == splitted[0]):
          end_index=Variable(name=f'{var.name}%{splitted[1]}', parent=var, scope=routine)
          if self.verbose : print(colored("derived type " + splitted[0] + " found", "green"))
    else:
      for var in FindVariables().visit(routine.variables):
        if (var.name == self.horizontal.bounds[1]):
          end_index = var
          if self.verbose : print(colored("variable " + var.name + " found", "green"))


    if not (begin_index or end_index):
      raise RuntimeError(f'index variables not found in routine {routine.name}')


    for assign in FindNodes(Assignment).visit(routine.body):
      expression_map={}             
      for var in FindVariables().visit(assign):
        if isinstance(var, Array):
          if (len(var.dimensions) > 0):
            if (var.dimensions[0] == ':'):              
              # Check if variable seems to be of the correct size
              found_in_routine_vars = False
              if (var.name in routine.variable_map) :
                if self.verbose : print ("Variable", var.name, " found in routine variables")
                routine_var =routine.variable_map[var.name]
                dim_name = routine_var.dimensions[0].name
                if (dim_name == self.horizontal.size):
                  if self.verbose : print(colored("First dimension of array is "+dim_name, "green"))
                  newrange = RangeIndex( (begin_index, end_index, None) )
                  newdimensions = (newrange,) + var.dimensions[1:]
                  new_var = var.clone(dimensions=newdimensions)
                  expression_map[var] = new_var
                else:
                  print(colored("Unexpected first dimension of array : " + dim_name, "red"))

              else :
                print (colored("Variable not found in routine variables !", "red"))

      if expression_map:
        total += len(expression_map)
        explicited_assign = SubstituteExpressions(expression_map).visit(assign)
        assign_map[assign] = explicited_assign
        if self.verbose : print("Transforming ", assign)
        if self.verbose : print("        into ", explicited_assign,  "\n")
            
    routine.body=Transformer(assign_map).visit(routine.body)
    if (total > 0) : info(f'[Loki] {routine.name}:: {total} implicit array syntax expressions replaced with explicit boundaries')
