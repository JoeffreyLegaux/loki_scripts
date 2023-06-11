# from loki import FP, Sourcefile, Dimension, Subroutine
from loki import ( Frontend, Sourcefile, Scheduler, FindNodes, FindScopes, Loop, Variable, 
                Assignment, CallStatement, Transformation, Node, SymbolAttributes, DerivedType, 
                BasicType, Import, Transformer, flatten, pragmas_attached )
from loki.expression import LoopRange, FindVariables
from loki.expression.symbols import (
    Array, Scalar, InlineCall, TypedSymbol, FloatLiteral, IntLiteral, LogicLiteral,
    StringLiteral, IntrinsicLiteral, DeferredTypeSymbol, RangeIndex
)
from loki.ir import Section, Comment, VariableDeclaration, Pragma, Intrinsic
from pathlib import Path
from termcolor import colored
from loki.frontend.fparser import *
import sys
# Bootstrap the local transformations directory for custom transformations
sys.path.insert(0, str(Path(__file__).parent))
print("path  = ", str(Path(__file__).parent))
print("sys.path  = ", sys.path)
print ("dir =", dir)




class InsertCrayPointers(Transformation):
  def __init__(self, verbose = False):
    self.verbose = verbose
    self.stack_module_name="stack_mod"
    self.stack_type_name="STACK"
    self.stack_argument_name="YDSTACK"
    self.stack_local_name="YLSTACK"
    self.stack_global_name="YSTACK"
    self.stack_array_name="ZSTACK"

  def transform_subroutine(self, routine, **kwargs):
    role = kwargs['role']
    targets = kwargs.get('targets', None)
    if self.verbose : print("targets : ", targets)
    if role == 'driver':
      self.process_driver(routine, targets=targets)
    if role == 'kernel':
      self.process_kernel(routine, targets=targets)

  def process_kernel(self, routine, targets):
    if self.verbose :
      print(" ======================================= ")
      print("   Transforming routine ", routine.name)
      print(" ======================================= ")

    stack_argument = Variable(name=self.stack_argument_name, type=SymbolAttributes(DerivedType(name=self.stack_type_name), intent='in' ))
    stack_local = Variable(name=self.stack_local_name, type=SymbolAttributes(DerivedType(name=self.stack_type_name)), scope = routine)
    stack_member_L = Variable(name=f'{stack_local.name}%L', parent=stack_local, scope = routine)

    routine.arguments += (stack_argument,)
    routine.variables += (stack_argument, stack_local, )

    routine.spec.prepend(Import(module=self.stack_module_name))

    if self.verbose :print("Added stack module import statement")

    # Put arguments in a list to circumvent very slow direct accesses
    list_args = list(routine.arguments)

    local_array_variables = []
    for decl in FindNodes(VariableDeclaration).visit(routine.spec):
      for v in decl.symbols:
        if isinstance(v, Array):
          if ( v not in list_args):
            local_array_variables.append(v)

    if self.verbose : print("local arrays declarations found : ", len(local_array_variables))

    assignments = []
    declarations = []

    for n in FindNodes(Pragma).visit(routine.spec):
      if (n.keyword == "acc" and "routine" in n.content):
        main_pragma_node = n

    assignments.append(Assignment(lhs=stack_local, rhs=stack_argument))

    for v in local_array_variables:

      declarations.append(Intrinsic(text=f'POINTER(IP_{v.name}_, {v.name})' ))

      assignments.append(Assignment(lhs=Variable(name=f'IP_{v.name}_'), rhs=stack_member_L))      

      # Increase stack boundary by the size of the array
      # If variable does not have an explicit kind, use the first element instead
      if (v.type.kind == None):
        type_string=f'SIZEOF({v.name}(1'   
        for i in range(1, len(v.dimensions)):
          type_string+=',1'
        type_string+='))'
      else:
        type_string=str(v.type.kind) 

      temp = f'{stack_member_L.name} + {type_string} * SIZE ({v.name})'

      assignments.append(Assignment(lhs=stack_member_L, rhs=parse_fparser_expression(temp, scope=routine) ) )

    # replace pragma with declarations, then append the pragma afterwards
    map_declarations = {}
    map_declarations[main_pragma_node] = tuple(declarations)
    routine.spec=Transformer(map_declarations).visit(routine.spec)

    routine.spec.append(main_pragma_node)

    routine.body.prepend(assignments)


    call_mapper = {}
    for call in FindNodes(CallStatement).visit(routine.body):
      if (call.name in targets) :
        new_call = call.clone(arguments=call.arguments)
        new_call._update(kwarguments=new_call.kwarguments + ((stack_argument.name, stack_local),))
        call_mapper[call] = new_call

    routine.body = Transformer(call_mapper).visit(routine.body)

    if self.verbose :
      print("Stack argument added to Call statements ")
      for call in call_mapper:
        print (call.name)    

      print(" ======================================= ")
      print("   Routine ", routine.name, " transformed")
      print(" ======================================= \n ")


  def process_driver(self, routine, targets):
    if self.verbose :
      print(" ======================================= ")
      print("   Transforming main routine ", routine.name)
      print(" ======================================= ")

    stack_argument = Variable(name=self.stack_argument_name, type=SymbolAttributes(DerivedType(name=self.stack_type_name), intent='in' ))
    stack_global = Variable(name=self.stack_global_name, type=SymbolAttributes(DerivedType(name=self.stack_type_name) ), scope=routine)

    stack_member_L = Variable(name=f'{stack_global.name}%L', parent=stack_global, scope = routine)
    stack_member_U = Variable(name=f'{stack_global.name}%U', parent=stack_global, scope = routine)
    # Add stack module usage
    routine.spec.prepend(Import(module=self.stack_module_name))

    stack_array = Array(name=self.stack_array_name, \
      type=SymbolAttributes(BasicType.REAL, allocatable=True, kind=Variable(name='jprb')) , \
      dimensions= (RangeIndex((None, None, None)),RangeIndex((None, None, None)),) , scope = routine ) 

    routine.variables += (stack_global, stack_array,)


    with pragmas_attached(routine, Loop, attach_pragma_post=True):
      call_mapper={}
      first_call = True
      for call in FindNodes(CallStatement).visit(routine.body):
        if call.name in targets:
          if self.verbose : print ("call found : ", call)
          new_call = call.clone(arguments=call.arguments)
          new_call._update(kwarguments=new_call.kwarguments + ((stack_argument.name, stack_global),))

          ancestors = flatten(FindScopes(call).visit(routine.body))
          loops = [a for a in ancestors if isinstance(a, Loop)]
          
          if not loops:
            call_mapper[call] = new_call
          else : 
            dim = loops[0].variable
            
            if first_call:
              alloc_node = parse_fparser_expression(f'ALLOCATE({stack_array.name}(10000 * NPROMA, {loops[0].bounds.upper}))', scope=routine)
              call_mapper[loops[0]] = tuple([alloc_node, loops[0]])
              first_call = False
    
            temp = f'LOC ( {stack_array.name} (1, {dim} ))'
            assign_L = Assignment(lhs = stack_member_L, rhs = parse_fparser_expression(temp, scope=routine) )
            temp = f'{stack_member_L.name} + JPRB * NPROMA * 10000'
            assign_U = Assignment(lhs = stack_member_U, rhs = parse_fparser_expression(temp, scope=routine) )

            call_mapper[call] = tuple([assign_L, assign_U, new_call])

            if loops[1].pragma :
              prgm = loops[1].pragma[0]
              loops[1]._update(pragma=(Pragma(keyword=prgm.keyword, content=prgm.content + f' private({stack_global.name})'),) )

      routine.body = Transformer(call_mapper).visit(routine.body)
