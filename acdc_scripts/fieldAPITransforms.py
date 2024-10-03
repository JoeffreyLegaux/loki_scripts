from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType  )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

from loki.expression import FindTypedSymbols, FindVariables
from loki.expression.expr_visitors import SubstituteExpressions
from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from storable import retrieve

import re

from arpege_parameters import params

def is_FieldAPI_ARRAY(typename):
    return (re.search("^ARRAY_(|INT|LOG)\dD$", typename) != None)

def get_fieldAPI_member(var, types):
    head = var[0]

    if head not in types : return None
    if len(var) > 1:
        tail = var[1:]
        return get_fieldAPI_member(tail, types[head])
    else :
        if isinstance(types[head], dict):
            return None
        else :
            return types[head]

def get_FieldAPI_variables(routine, fieldAPI_types):
    fieldAPI_variables = {}
    for var in routine.variables :
        if (isinstance(var.type.dtype, DerivedType)):
            typename = var.type.dtype.name
            # Variables of type ARRAY_nD are automatically generated FieldAPI   
            if is_FieldAPI_ARRAY(typename) :
                fieldAPI_variables[var.name] = typename
            elif typename in fieldAPI_types:
                fieldAPI_variables[var.name] = typename
    return fieldAPI_variables


class FieldAPIPtr(Transformation):
    def __init__(self, pointerType='host', node=None):
        if (pointerType == 'host') :
            self.pointerSuffix = 'PTR'
        elif (pointerType == 'device') :
            self.pointerSuffix = 'DEVPTR'
        else : 
            error(f'Unrecognised pointer type : {pointerType}')
            exit(1)

        if node:
            assert hasattr(node, "body")
        self.node = node
        

    def transform_subroutine(self, routine, **kwargs):

        fieldAPI_types = retrieve('../scripts/types.dat')

        fieldAPI_variables = get_FieldAPI_variables(routine, fieldAPI_types)

        variables_map = {}

        boundary_variable = None

        for var in routine.variables :
            if (isinstance(var.type.dtype, DerivedType)):
                if (var.type.dtype.name == params.boundaries_type):
                    boundary_variable = var

        block_index = Variable(name='KBL', parent=boundary_variable, scope=routine)


        body = self.node.body if self.node else routine.body

        for var in FindVariables().visit(body) :
            
            base = var.name_parts[0]

            if base in fieldAPI_variables:

                fAPI_base = fieldAPI_variables[base]
                # Specific treatment for ARRAY_nD variables
                if is_FieldAPI_ARRAY(fAPI_base):
                    ndim = int(fAPI_base[-2])

                    if var.dimensions :
                        dimensions = var.dimensions
                    else : 
                        dimensions = ()
                        for i in range(ndim-1):    
                            dimensions += (RangeIndex( (None, None, None) ) ,)

                    dimensions = dimensions + (block_index,)

                    fAPI_var = Variable(name = "F_P", parent = Variable(name = base))                    
                    variables_map[var] = Array(name=self.pointerSuffix, parent=fAPI_var, dimensions=dimensions)

                else :
                    #Build the corresponding FieldAPI variable name
                    member = var.name_parts[1:]
                    member = [fAPI_base] + member
                    fAPI_member = get_fieldAPI_member(member, fieldAPI_types )
                    
                    if fAPI_member:                    
                        if isinstance(var, Array):
                            dimensions = var.dimensions
                        else:
                            dimensions = () 
                            for i in range(fAPI_member[1]):    
                                dimensions += (RangeIndex( (None, None, None) ) ,)

                        dimensions += (block_index,)
                        fAPI_var = Variable(name = fAPI_member[0], parent=var.parent, scope=routine)
                        variables_map[var] = Array(name=self.pointerSuffix, scope=routine, parent=fAPI_var, dimensions=dimensions)

        if self.node :
            # print("variables map  ", variables_map)
            new_node = self.node.clone(body = SubstituteExpressions(variables_map).visit(self.node.body) )

            routine.body = Transformer({self.node:new_node}).visit(routine.body)
        else:
            routine.body = SubstituteExpressions(variables_map).visit(routine.body)


