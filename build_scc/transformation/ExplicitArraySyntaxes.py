# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.
from loki import (Frontend, Sourcefile, Scheduler, FindNodes, Loop, Variable,
            Assignment, CallStatement, Transformation, Node, SymbolAttributes, 
            DerivedType, BasicType, Import, Transformer, Conditional, SchedulerConfig,
            FindVariables, FindTypedSymbols, SubstituteExpressions)
from loki.expression import LoopRange
from loki.expression.symbols import (
    Array, Scalar, InlineCall, TypedSymbol, FloatLiteral, IntLiteral, LogicLiteral,
    StringLiteral, IntrinsicLiteral, DeferredTypeSymbol, LogicalOr, LogicalAnd, LogicalNot, RangeIndex
)
from loki.ir import Section, Comment, VariableDeclaration, Import, Intrinsic
from pathlib import Path
from termcolor import colored
from loki.logging import info
import sys

import re


def is_FieldAPI_ARRAY(typename):
    #return (re.search("^FIELD_\d(IM|LM|RD|RM|RB)_ARRAY$", typename) != None)
    return (re.search("^FIELD_\d(IM|LM|RD|RM|RB)$", typename) != None)


# Check if variable is used in array-syntax
# Either full variable, either its whole first dimension is used
def check_dimensions(variable):
    if (len(variable.dimensions) == 0):
        return True
    else :
        return (variable.dimensions[0] == ':')

def ExplicitArraySyntaxes(routine, lst_horizontal_size, lst_horizontal_bounds, FieldAPI_pointers=None):
  """
  Remove array syntax for the horizontal dimension, and returns the upper and lower bounds of the horizontal loops.
  :param routine:.
  :param lst_horizontal_size: list of aliases of NPROMA
  :param lst_horizontal_bounds: list of diff names for the horizontal lower and upper bounds
  """
  total = 0
  assign_map={}
  
  begin_index = None
  end_index = None

  verbose=False
#  verbose=True

  define=False
  splitted1 = lst_horizontal_bounds[0][1].split('%')
  splitted2 = lst_horizontal_bounds[1][1].split('%')
  
  # Research of the horizontal bounds
  for var in FindVariables().visit(routine.variables):

    # Lower bound:
    if (var.name == splitted1[0]):
      begin_index=Variable(name=f'{var.name}%{splitted1[1]}', parent=var, scope=routine)
      if verbose : print(colored("derived type " + splitted1[0] + " found", "green"))


    if (var.name == lst_horizontal_bounds[0][0]):
      begin_index = var
      if verbose : print(colored("variable " + var.name + " found", "green"))

    if (var.name == lst_horizontal_bounds[0][2] or var.name == lst_horizontal_bounds[0][3]):
      begin_index = var
      if verbose : print(colored("variable " + var.name + " found", "green"))

    # Upper bound :
    if (var.name == splitted2[0]):
      end_index=Variable(name=f'{var.name}%{splitted2[1]}', parent=var, scope=routine)
      if verbose : print(colored("derived type " + splitted2[0] + " found", "green"))

    if (var.name == lst_horizontal_bounds[1][0]):
      end_index = var
      if verbose : print(colored("variable " + var.name + " found", "green"))

    if (var.name == lst_horizontal_bounds[1][2]):
      end_index = var
      if verbose : print(colored("variable " + var.name + " found", "green"))

  if verbose: print("begin_index=",begin_index) #lower bound
  if verbose: print("end_index=",end_index) #upper bound

  if not (begin_index or end_index):
    raise RuntimeError(f'index variables not found in routine {routine.name}')

  new_range = RangeIndex( (begin_index, end_index, None) )

  # Build list of pointers to nproma fields

  #print ("Field api pointers list : ", FieldAPI_pointers)
  verbose = False  
  #remove array syntax
  for assign in FindNodes(Assignment).visit(routine.body):
    if assign.ptr :
        if verbose : print("assign to pointer ignored : ", assign)
    else :
      is_pointer=False
      not_found=[]  
      expression_map={}             
      for var in FindVariables().visit(assign):
        if (var.type.pointer):
          is_pointer=True
          if verbose : print("Variable", var.name, " is a pointer")
        if isinstance(var, Array):
          if verbose : print("Variable", var.name, " is an array of dimensions : ", var.dimensions)

          # Only treat arrays whose first dimension is ':' OR have no dimensions at all
          if check_dimensions(var):

            if (var.name in routine.variable_map) :
              routine_var = routine.variable_map[var.name]
              if var.type.pointer or routine_var.type.pointer:
                  is_pointer = True
              
              if verbose : print ("Variable", var.name, " found in routine variables : ")
              #if not is_pointer:
              is_var_FieldAPI = is_FieldAPI_ARRAY(var.type.dtype.name)
              #else :
              if not is_var_FieldAPI and is_pointer:
                is_var_FieldAPI = var.name in FieldAPI_pointers.keys()
              if verbose : print ("is FieldAPI ? ", is_var_FieldAPI)

              if not (is_var_FieldAPI):
                #routine_var = routine.variable_map[var.name]
                if len(var.dimensions) > 0:
                  dim_name = routine_var.dimensions[0].name
                  if (dim_name in lst_horizontal_size):
                    #if verbose : 
                    print(colored("First dimension of array is "+dim_name, "green"))
                    is_var_FieldAPI = True
                  else:
                    print(colored("Unexpected first dimension of array : " + dim_name, "red"))


              if is_var_FieldAPI:
  
                new_dimensions = (new_range ,)

                if len(var.dimensions) == 0 :
                    if is_pointer :
                        dims = FieldAPI_pointers[var.name] - 1
                    else :
                        dims = len(routine.variable_map[var.name].dimensions) - 1
                    for i in range(dims):
                         new_dimensions += (RangeIndex((None, None, None)),)
                else:
                    new_dimensions += var.dimensions[1:]
                define=True
                new_var = var.clone(dimensions=new_dimensions)
                expression_map[var] = new_var

            else :
              not_found.append(var.name)

      if not is_pointer : 
        for var in not_found :
          print (colored("Variable not found in routine variables !", "red"))
          if verbose : print ("Variable", var, " not found in routine variables !")

      if expression_map:
        total += len(expression_map)
        explicited_assign = SubstituteExpressions(expression_map).visit(assign)
        assign_map[assign] = explicited_assign
        if verbose : print("Transforming ", assign)
        if verbose : print("        into ", explicited_assign,  "\n")
          
  routine.body=Transformer(assign_map).visit(routine.body)
  if verbose:
      if (total > 0) : info(f'[Loki] {routine.name}:: {total} implicit array syntax expressions replaced with explicit boundaries')
  if define==False :
    newrange=None
  return(end_index, begin_index, new_range)
