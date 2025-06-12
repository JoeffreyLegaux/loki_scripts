from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType  )

from loki.ir import (Section, Comment, CommentBlock, VariableDeclaration, Pragma, Import, Assignment, Conditional, LeafNode,
                     InternalNode, Associate, FindTypedSymbols, FindVariables, SubstituteExpressions)

from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from storable import retrieve

import re

from arpege_parameters import params


def is_fieldAPI_ARRAY(typename):

    #return (re.search("^FIELD_\d(IM|LM|RD|RM|RB)_ARRAY$", typename) != None)
    return (re.search("^FIELD_\d(IM|LM|RD|RM|RB)$", typename) != None)

# extract the fieldAPI member corresponding to a variable represented
# in the form of a list of its derived members
def get_fieldAPI_member(var, types = None):
    head = var[0]
    # Bootstrap recursion with complete FieldAPI types dict
    if not types:
        types = params.fieldAPI_types
    if head not in types : return None
    if len(var) > 1:
        tail = var[1:]
        return get_fieldAPI_member(tail, types[head])
    else :
        if isinstance(types[head], dict):
            return None
        else :
            return types[head]

def get_fieldAPI_variables(routine):
    fieldAPI_variables = {}
    for var in routine.variables :
        if (isinstance(var.type.dtype, DerivedType)):
            typename = var.type.dtype.name
            # Variables of type Field_nxx are automatically generated FieldAPI   
            if is_fieldAPI_ARRAY(typename) :
                fieldAPI_variables[var.name] = typename
            elif typename in params.fieldAPI_types:
                fieldAPI_variables[var.name] = typename
    return fieldAPI_variables


def get_pointers_to_FieldAPI(routine, nproma_variables):
    #print("nproma vairiables : ", nproma_variables)
    ptr_list = []
    for var in routine.variables :
        if var.type.pointer and var.type.dtype.name != 'FIELD_BASIC':
            ptr_list.append(var.name)
    #print("ptr_list : ", ptr_list)
    FieldAPI_ptrs_dims = {}
    FieldAPI_ptrs_vars = {}
    # probably won't encounter pointers to defferedtype so don't initlialize those
    fieldAPI_variables = None
    for assign in FindNodes(Assignment).visit(routine.body):
        if assign.ptr :
            # Only treat values once
            if assign.lhs.name in ptr_list and assign.lhs.name not in FieldAPI_ptrs_vars:
                # only treat direct assignment to arrays variables
                # Esp. avoid binding to NULL()
                #print("rhs ? ", assign.rhs, type(assign.rhs)) 
                if isinstance(assign.rhs, Array):
                    check = False
                    if is_fieldAPI_ARRAY(assign.rhs.type.dtype.name):
                        check = True
                        dims = int(assign.rhs.type.dtype.name[6])
                    elif assign.rhs.type.shape:
                        # Array with shape is plain old array, check if nproma size
                        if assign.rhs.type.shape[0].name in nproma_variables :
                            check = True
                            dims = len(assign.rhs.type.shape)
                    else:
                        # Array with no shape means deferredtype, it could be fieldapi !
                        base = assign.rhs.name_parts[0]                        
                        if not fieldAPI_variables:
                            fieldAPI_variables = get_fieldAPI_variables(routine)
                        #print("FieldAPI vars : ", fieldAPI_variables)
                        if base in fieldAPI_variables:
                            fAPI_base = fieldAPI_variables[base]
                            #Build the corresponding FieldAPI variable name
                            member = var.name_parts[1:]
                            member = [fAPI_base] + member
                            fAPI_member = get_fieldAPI_member(member)
                            if fAPI_member:
                                check = True
                                dims = fAPI_member[1]


                    if check :
                        FieldAPI_ptrs_dims[assign.lhs.name] = dims
                        FieldAPI_ptrs_vars[assign.lhs.name] = assign.rhs
    #print("FAPIptrs : ", FieldAPI_ptrs)
    return (FieldAPI_ptrs_dims, FieldAPI_ptrs_vars)


def nproma_to_FieldAPI(routine, FieldAPI_pointers):
        """
        Search nproma-sized arrays (from params.nproma_aliases),
        transforms them into relevant FIELD_NXX type,
        return the new arrays names along their new dimensions
        """
        arguments_map = {}
        nproma_arrays_map = {}
        nproma_names_map = {}
        nproma_arrays_dimensions = {}
        locals_dimensions = {}
        
        arrays_types_dimensions = set()

        field_new_calls = ()
        field_delete_calls = ()

        for var in FindVariables().visit(routine.spec):

            if isinstance(var, Array):

                is_FAPI_pointer = var.name in FieldAPI_pointers

                is_nproma = False
                if not is_FAPI_pointer :
                    firstdim = var.dimensions[0] 
                
                    #Naively considering ranges with (xxx:nproma) are nproma-sized
                    if isinstance(firstdim, RangeIndex):
                        firstdim = firstdim.upper

                    is_nproma = firstdim in params.nproma_aliases
                
                if is_FAPI_pointer or is_nproma:

                    is_argument = var in routine.arguments
                
                    #print("variable to FAPI : ", var, type(var))
                    #print("type : ", type(var.type.kind.name))

                    array_type_dim = f'{len(var.dimensions)+1}{var.type.kind.name[-2:]}'


                    arrays_types_dimensions.add(array_type_dim)
                        
                    new_var = Variable(     name=f'Y{("D" if is_argument else "L")}_{var.name}', 
                                            type=var.type.clone(dtype=DerivedType(name=f'FIELD_{array_type_dim}'),
                                                                kind=None, 
                                                                # Abstract types require polymorphism
                                                                polymorphic=True, 
                                                                # Dummy arrays can stay target, local arrays have to become
                                                                # pointers which conflict with target attribute
                                                                target = False,
                                                                pointer = True
                                                                )
                                    )

                    nproma_arrays_map[var] = new_var
                    nproma_names_map[var.name] = new_var
                    nproma_arrays_dimensions[var.name] = var.dimensions

        routine.spec = SubstituteExpressions(nproma_arrays_map).visit(routine.spec)
        
        # Updates routine signature with new arguments through the _dummies property
        dummies_list = list(routine._dummies)
        for index, dummy in enumerate(dummies_list):
            if dummy.upper() in nproma_names_map:
                dummies_list[index] = nproma_names_map[dummy.upper()].name
        
        routine._dummies = tuple(dummies_list)

        body_vars_map={}
        for var in FindVariables().visit(routine.body):
            if isinstance(var, Array):
                if var.name in nproma_names_map:
                    new_var = nproma_names_map[var.name]
                    body_vars_map[var] = new_var.clone(dimensions = var.dimensions)

        routine.body = SubstituteExpressions(body_vars_map).visit(routine.body)

        return None

class FieldAPIPtr(Transformation):
    def __init__(self, pointerType='host'):
        if (pointerType == 'host') :
            self.pointerSuffix = 'PTR'
        elif (pointerType == 'device') :
            self.pointerSuffix = 'DEVPTR'
        else : 
            error(f'Unrecognised pointer type : {pointerType}')
            exit(1)

       

    def transform_node(self, node, routine, inplace = False):
        if node:
            assert hasattr(node, "body")
        
        fieldAPI_variables = get_fieldAPI_variables(routine)

        variables_map = {}

        block_index = Variable(name=params.block_counter, scope=routine)

        for var in FindVariables().visit(node.body) :
            
            base = var.name_parts[0]

            if base in fieldAPI_variables:

                fAPI_base = fieldAPI_variables[base]
                # Specific treatment for ARRAY_nD variables
                if is_fieldAPI_ARRAY(fAPI_base):
                    ndim = int(fAPI_base[6])

                    if hasattr(var, "dimensions") and var.dimensions :
                        dimensions = var.dimensions
                    else : 
                        dimensions = ()
                        for i in range(ndim-1):    
                            dimensions += (RangeIndex( (None, None, None) ) ,)

                    dimensions = dimensions + (block_index,)

                    #fAPI_var = Variable(name = "F_P", parent = Variable(name = base))                    
                    fAPI_var = Variable(name = base)
                    variables_map[var] = Array(name=self.pointerSuffix, parent=fAPI_var, dimensions=dimensions)

                else :
                    #Build the corresponding FieldAPI variable name
                    member = var.name_parts[1:]
                    member = [fAPI_base] + member
                    fAPI_member = get_fieldAPI_member(member)
                    
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
                        #print( variables_map[var])

        if inplace:
            node.body = SubstituteExpressions(variables_map).visit(node.body)
            return None
        else:
            new_node = node.clone(body = SubstituteExpressions(variables_map).visit(node.body))
            return new_node


    def transform_subroutine(self, routine, **kwargs):
        self.transform_node(routine, routine, inplace=True)
