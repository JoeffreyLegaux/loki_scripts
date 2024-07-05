from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType  )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

from loki.expression import FindTypedSymbols, FindVariables
from loki.expression.expr_visitors import SubstituteExpressions
from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, VariableSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from storable import retrieve

from arpege_parameters import params

import re



# extract the fieldAPI member corresponding to a variable represented
# in the form of a list of its derived members
def get_fieldAPI_member_lst(var, types):
    head = var[0]

    if head not in types : return None
    if len(var) > 1:
        tail = var[1:]
        return get_fieldAPI_member_lst(tail, types[head])
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
            if typename in fieldAPI_types:
                fieldAPI_variables[var.name] = typename
    return fieldAPI_variables

def addFieldAPIPointers(routine, number_of_pointers):
    for i in range(number_of_pointers):
        routine.variables += (Variable( name=f'YLFLDPTR{i}', 
                                        type=SymbolAttributes(  DerivedType(name='FIELD_BASIC'),
                                                                pointer = True, 
                                                                polymorphic=True
                                                                ),
                                        scope = routine
                                        ) 
                              ,)            


class MakeSync(Transformation):


    def __init__(self, pointerType='host', sections=None, nproma_arrays=None):
        if (pointerType == 'host') :
            self.callSuffix = 'SYNC_HOST'
        elif (pointerType == 'device') :
            self.callSuffix = 'SYNC_DEVICE'
        else : 
            error(f'Unrecognised pointer type : {pointerType}')
            exit(1)
        self.map_reads = {}
        self.map_writes = {}
        self.map_nodes = {}
        self.nproma_vars_names = nproma_arrays if nproma_arrays else []
        self.sections = sections
        self.total_FAPI_pointers=0


    def get_fieldAPI_member(self, var, types):
        base_name = var.name_parts[0]
        if base_name in self.fieldAPI_variables:
            member_name = var.name_parts[1:]
            member_name = [self.fieldAPI_variables[base_name]] + member_name
            return get_fieldAPI_member_lst(member_name, self.fieldAPI_types )
            
        return None

    def isNpromaOrFieldAPI(self, var):
        base_name = var.name_parts[0]

        # Variables recasted as ARRAY_nD have become scalars ! 
        if isinstance(var, Array) or isinstance(var, Scalar):
            if base_name in self.nproma_vars_names:
                # If we have passed sections, we are using nproma arrays turned into ARRAY_nD arrays
                # Therefore we have to replace it by its FieldAPI pointer
                if self.sections:
                    return True, Variable(name = "F_P", parent = var) #, scope=routine)

                return True, var

        if isinstance(var, DeferredTypeSymbol):
            # if base_name in self.fieldAPI_variables:
                # member_name = var.name_parts[1:]
                # member_name = [self.fieldAPI_variables[base_name]] + member_name
            fAPI_member = self.get_fieldAPI_member(var, self.fieldAPI_types )
            if fAPI_member:
                return True, var.clone(name=fAPI_member[0])

        return False, None


    def npromaVariables(self, routine):
        var_list=[]
        for var in routine.variables:
            if isinstance(var, Array):
                # print("var found ", var)
                firstdim = var.dimensions[0] 
                if isinstance(firstdim, RangeIndex):
                    firstdim = firstdim.upper
                if firstdim :
                    if (firstdim.name == 'KLON' or firstdim.name == 'YDCPG_OPTS%KLON'):
                        var_list.append(var.name)
        return var_list
        
    # Transforms KLON arrays into POINTERS
    # return the list of names of transformed variables 
    def declarationsToPointers(self, routine, nproma_list):
        var_map = {}
        for var in routine.variables:
            if var.name in nproma_list:
                new_var = Variable(name=var.name, 
                            type=SymbolAttributes(
                                    DerivedType(name='FIELD_BASIC'), 
                                    intent=var.type.intent, pointer=True, polymorphic=True,
                                    # Local variables point to NULL()
                                    initial=None if var in routine.arguments else InlineCall(DeferredTypeSymbol('NULL'))
                                ),
                            scope = routine
                        )
                var_map[var] = new_var

        routine.variables = SubstituteExpressions(var_map).visit(routine.variables)

        return [v.name for v in var_map]


    # def findAssigns(self, node, balise, depth) : 
    def findAssigns(self, node) : 

        if isinstance(node, tuple):
            # sequence of actual nodes : get the list of unique reads and writes from each node

            writes_list = []
            reads_list = []

            # Empty nodes might appear from previous simplifications of the code, simply ignore them
            if (node != ()):
                for c in node:
                    (reads, writes) = self.findAssigns(c)

                    writes_list.extend([define for define in writes if define not in writes_list])
                    reads_list.extend([define for define in reads if define not in reads_list])

                # writes take precedence over reads
                reads_list = [r for r in reads_list if r not in writes_list]

                self.map_reads[node] = reads_list
                self.map_writes[node] = writes_list
            return(reads_list, writes_list)


        elif isinstance(node, Assignment):

            writes = [symbol
                        for test, symbol in [
                            self.isNpromaOrFieldAPI(s) 
                            for s in node.defines_symbols
                            ]
                        if test ]

            reads = [symbol
                        for test, symbol in [
                            self.isNpromaOrFieldAPI(s) 
                            for s in node.uses_symbols
                            # Writes take precedence over reads
                            if s not in node.defines_symbols
                            ]
                        if test ]                        

            return (reads, writes)

        elif isinstance(node, Section):
            # also covers Associate nodes
            # just recurse on the body

            (reads,writes) = self.findAssigns(node.body)

            return(reads,writes)


        elif isinstance(node, Conditional):

            (if_reads, if_writes) = self.findAssigns(node.body)
            (else_reads, else_writes) = self.findAssigns(node.else_body)

            all_writes = [w for w in if_writes if w in else_writes]
            all_reads = [r for r in if_reads if r in else_reads]

            # writes take precedence over reads
            all_reads = [r for r in all_reads if r not in all_writes]


            # print("all writes : ", all_writes)
            self.map_reads[node] = all_reads
            self.map_writes[node] = all_writes
            
            return(all_reads, all_writes)
        else :
            # print(f' noeud non traité !!!! {type(node)} {hasattr(node,"body")} ' )
            return([],[])


    def createSyncCallStatement(self, var, mode):

        if (mode == 'W'):
            call_name = self.callSuffix + '_RDWR'
        elif (mode == 'R'):
            call_name = self.callSuffix + '_RDONLY'
        else:
            print("no mode for call statement creation !!!")
            exit(1)

        associate_call = InlineCall(function=DeferredTypeSymbol(name='ASSOCIATED'), parameters=(var,) )
        sync_call = CallStatement(name = DeferredTypeSymbol(name=call_name, parent = var), arguments=(), scope=None)
        cond = Conditional(condition=associate_call, body = (sync_call,), inline = True)
        return cond

    def clearAssigns(self, node, upper_reads, upper_writes):

        if isinstance(node, tuple):
            if (node == ()):
                assert node not in self.map_reads, ('empty node in map !')
                assert node not in self.map_writes, ('empty node in map !')
                return True
            else:
                
                new_writes = [w for w in self.map_writes[node] if w not in upper_writes]
                new_reads = [r for r in self.map_reads[node] if r not in upper_reads and r not in upper_writes and r not in new_writes]

                self.map_nodes[node] = node
                empty = True

                for w in new_writes:
                    self.map_nodes[node] += (self.createSyncCallStatement(w, 'W'), ) 
                    empty = False

                for r in new_reads:
                    self.map_nodes[node] += (self.createSyncCallStatement(r, 'R'), ) 
                    empty = False

                for c in node :
                    empty_c = self.clearAssigns(c, upper_reads + new_reads, upper_writes + new_writes) 
                    empty = (empty and empty_c)

                if self.map_nodes[node] == ():
                    print("new node is None, node is : ", node)
                    if node == ():
                        del self.map_nodes[node]
                        print("deleted node from map, this is the way")

                if empty:
                    self.map_nodes[node] = None

                return empty


        elif isinstance(node, Assignment):
            self.map_nodes[node] = None
            return True
        
        elif isinstance(node, Section):
            return self.clearAssigns(node.body, upper_reads, upper_writes)
        
        elif isinstance(node, Conditional):

            new_writes = [w for w in self.map_writes[node] if w not in upper_writes]
            new_reads = [w for w in self.map_reads[node] if w not in upper_reads and w not in new_writes]

            new_nodes = ()
            if (new_writes != []):
                print("IF new_writes ", new_writes)
                new_nodes += Intrinsic(text=f'!IF writes {new_writes}')
            if (new_reads != []):
                print("IF new_reads ", new_reads)
                new_nodes += Intrinsic(text=f'!IF reads {new_reads}')
                

            empty_body = self.clearAssigns(node.body, upper_reads + new_reads, upper_writes + new_writes)
            empty_else = self.clearAssigns(node.else_body, upper_reads + new_reads, upper_writes + new_writes)

            if empty_else and empty_body :
                self.map_nodes[node] = None            

            return (empty_body and empty_else)


        elif isinstance(node, Comment):
            return True

        elif isinstance(node, CallStatement):
            return False

        else :
            print(f'noeud non traité !!!! {type(node)} {node} ' )
            return False



    def transform_subroutine(self, routine, **kwargs):

        # If there are no sections, we are treating a whole routine (GENERATE directive)
        # The main routine has turned its NPROMA local/arguments arrays into FieldAPI pointers
        # Therefore these array become pointers in the SYNC routine
        if not self.sections: 
            for var in FindVariables().visit(routine.spec):
                if isinstance(var, Array):
                    firstdim = var.dimensions[0] 
                    if isinstance(firstdim, RangeIndex):
                        firstdim = firstdim.upper
                    if firstdim in params.nproma_aliases:
                        self.nproma_vars_names.append(var.name)
            self.nproma_vars_names = self.declarationsToPointers(routine, self.nproma_vars_names)
        # Otherwise, we do nothing as the parallel transforms provided at construction the 
        # transformed names of nproma variables



        # Create the dict of FieldAPI variables used in this routine
        self.fieldAPI_types = retrieve('../scripts/types.dat')
        self.fieldAPI_variables = get_FieldAPI_variables(routine, self.fieldAPI_types)

        
        # print("nproma vars : ", self.nproma_vars_names)
        # print("FieldAPI variables : ", self.fieldAPI_variables)

        
        # Subroutines called inside a SYNC routine will only accept FIELD_BASIC arguments for FieldAPI variables.
        # We find those calls and swap FieldAPI arguments with a generic FIELD_BASIC pointer.
        # The replacement of variable usage with SYNC calls would erase the added pointers association,
        # so here we only build the map which will be applied later.

        calls_map = {}
        total_FAPI_pointers = 0
        for call in FindNodes(CallStatement).visit(routine.body):
            if call.name != "DR_HOOK":
                args_to_fAPI = {}
                print("call found : ", call)
                for arg in call.arguments:
                    if isinstance(arg,  Array):
                        if re.search("^ARRAY_.D$", arg.type.dtype.name):
                            args_to_fAPI[arg] = arg
                        elif arg.type.dtype.name != 'FIELD_BASIC' :
                            print('array not FIELD_BASIC or ARRAY_nD passed to subroutine !!!!!', arg, arg.type.dtype.name )

                    elif isinstance(arg, DeferredTypeSymbol):
                        fAPI_member = self.get_fieldAPI_member(arg, self.fieldAPI_types )
                        if fAPI_member :
                            args_to_fAPI[arg] = Variable(name = fAPI_member[0], parent=arg.parent, scope=routine)
                        
                if args_to_fAPI :
                    assignments = []
                    count = 0
                    args_to_pointers={}
                    for arg in args_to_fAPI:
                        args_to_pointers[arg] = Variable(name=f'YLFLDPTR{count}', scope=routine)
                        assignments.append(Assignment(lhs = args_to_pointers[arg], rhs = args_to_fAPI[arg], ptr = True))
                        count = count + 1
                    self.total_FAPI_pointers = max(count, self.total_FAPI_pointers)
                    new_call = call.clone(arguments = SubstituteExpressions(args_to_pointers).visit(call.arguments))
                    calls_map[call] = tuple(assignments) + (new_call,)

        # Add the required numbers of FIELD_BASIC pointers if there are no sections (full routine transformation)
        # If sections are present, the main transformation process will retrieve the # of pointers
        if not self.sections:
            addFieldAPIPointers(routine, self.total_FAPI_pointers)


        sect = routine.body

        if self.sections:
            if len(self.sections) > 1:
                print("error : multiple sections passed to Sync tranformation, aborting")
                return 0
            else :
                # print("une seule section")
                sect = self.sections[0]


        # We might end up with a condition with has_elseif attribute but without else_body.
        # This crashes reconstruction of the node, so we preventively eliminate all those attributes.
        for cond in FindNodes(Conditional).visit(routine.body):
            if cond.has_elseif:
                cond._update(has_elseif = False)


        # Initiate the data analysis process with loki built-in function
        map = {}
        with dataflow_analysis_attached(routine) :
            # First step : bottom-up traversal of the section
            # This step recursively propagates upwards variables read and writes when they
            # are common to multiple branches in the routine
            (top_reads, top_writes) = self.findAssigns(routine.body)

        top_reads = [r for r in top_reads if r not in top_writes]

        # Second step : top-down traversal of the section
        # This step inserts sync calls for variables at the highest level (identified in first step)
        self.clearAssigns(routine.body, top_reads, top_writes)


        # Search for position of first DR_HOOK call block (should be 0, but check anyway)
        idx = 0
        for stmt in routine.body.body:
            if isinstance(stmt, Conditional):
                if (stmt.condition == DeferredTypeSymbol(name='LHOOK')):
                    break
            idx = idx + 1


        for w in top_writes:
            routine.body.insert(idx+1, self.createSyncCallStatement(w, 'W')) 
        for r in top_reads:
            routine.body.insert(idx+1, self.createSyncCallStatement(r, 'R')) 

        routine.body = Transformer(self.map_nodes).visit(routine.body)

        # Remove inlining from conditionnals that contain more than one statement        
        cond_map = {}
        for cond in FindNodes(Conditional).visit(routine.body):
            if cond.inline: 
                if ( (cond.body and len(cond.body) > 1) or 
                    (cond.else_body and len(cond.else_body) > 1 ) ):
                    cond_map[cond] = cond.clone(inline=False)
        routine.body = Transformer(cond_map).visit(routine.body)


        # Transform the subroutine calls arguments into FIELD_BASIC pointers when relevant
        routine.body = Transformer(calls_map).visit(routine.body)
