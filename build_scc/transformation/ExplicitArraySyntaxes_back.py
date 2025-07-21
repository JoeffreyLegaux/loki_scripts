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
    return (re.search("^FIELD_\d(IM|LM|RD|RM|RB)_ARRAY$", typename) != None)


def get_pointers_to_FieldAPI(routine):
    #print("nproma variables : ", nproma_variables)
    ptr_list = []
    for var in routine.variables :
        if var.type.pointer and var.type.dtype.name != 'FIELD_BASIC':
            ptr_list.append(var.name)
    print("ptr_list : ", ptr_list)
    FieldAPI_ptrs = set()
    for assign in FindNodes(Assignment).visit(routine.body):
        if assign.ptr :
            if assign.lhs.name in ptr_list:
                print("assign tavu : ", assign)
                if is_FieldAPI_ARRAY(assign.rhs.type.dtype.name):
               # if assign.rhs.name in nproma_variables:
                    FieldAPI_ptrs.add(assign.lhs.name)
    print("FAPIptrs : ", FieldAPI_ptrs)
    return FieldAPI_ptrs


def ExplicitArraySyntaxes(routine, lst_horizontal_size, lst_horizontal_bounds):
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
  #research of the horizontal bounds
  for var in FindVariables().visit(routine.variables):

#lower bound:
    if (var.name == splitted1[0]):
      begin_index=Variable(name=f'{var.name}%{splitted1[1]}', parent=var, scope=routine)
      if verbose : print(colored("derived type " + splitted1[0] + " found", "green"))


    if (var.name == lst_horizontal_bounds[0][0]):
      begin_index = var
      if verbose : print(colored("variable " + var.name + " found", "green"))

    if (var.name == lst_horizontal_bounds[0][2]):
      begin_index = var
      if verbose : print(colored("variable " + var.name + " found", "green"))

#upper bound :
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


  FieldAPI_ptrs = get_pointers_to_FieldAPI(routine)
  verbose = True
  #remove array syntax
  for assign in FindNodes(Assignment).visit(routine.body):
    if assign.ptr :
        print("assigne a pointeur ignoré : ", assign)
    else :
      is_pointer=False
      not_found=[]  
      expression_map={}             
      for var in FindVariables().visit(assign):
        print("varfound : ", var)
        if (var.type.pointer):
          is_pointer=True
          if verbose : print("Variable", var.name, " is a pointer")
        if isinstance(var, Array):
          if verbose : print("Variable", var.name, " is an array", var.dimensions)
          if (len(var.dimensions) > 0):
            if (var.dimensions[0] == ':'):              

              if not (begin_index or end_index):
                  raise RuntimeError(f'index variables not found in routine {routine.name}')
              # Check if variable seems to be of the correct size
              found_in_routine_vars = False
              if (var.name in routine.variable_map) :
                if verbose : print ("Variable", var.name, " found in routine variables : ")
                is_var_FielDAPI = is_FieldAPI_ARRAY(var.type.dtype.name)

                if not (is_var_FielDAPI):
                  routine_var = routine.variable_map[var.name]
                  dim_name = routine_var.dimensions[0].name
                  if (dim_name in lst_horizontal_size):
                    if verbose : print(colored("First dimension of array is "+dim_name, "green"))
                    is_var_FielDAPI = True
                  else:
                    print(colored("Unexpected first dimension of array : " + dim_name, "red"))

                if is_var_FielDAPI:
                  newrange = RangeIndex( (begin_index, end_index, None) )
                  define=True
                  newdimensions = (newrange,) + var.dimensions[1:]
                  new_var = var.clone(dimensions=newdimensions)
                  expression_map[var] = new_var

#            elif (var.type.pointer) :
#              is_pointer=True
#              if verbose : print("Variable", var.name, " is a pointer")
              
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
  return(end_index, begin_index, newrange)
