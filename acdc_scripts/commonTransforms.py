from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType, FindTypedSymbols, FindVariables, SubstituteExpressions )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, PragmaRegion, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

from loki.transformations import inline_member_procedures

from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList
from loki.expression import symbolic
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from arpege_parameters import params

from codetiming import Timer
from storable import retrieve


# class RemoveImports(Transformation):
#     def transform_subroutine(self, routine, **kwargs):
#         imports_map = {}
#         for imp in FindNodes(Import).visit(routine.spec):
#             imp.c_import
#                 pragmas_map[pragma] = None
#         routine.spec = Transformer(pragmas_map).visit(routine.spec)

#class AddImports(Transformation):
#    def __init__(self, modules_names = []):
#        self.modules_names = modules_names
#    
#    def transform_subroutine(self, routine, **kwargs):
#        existing_modules = []
#        for imp in FindNodes(Import).visit(routine.spec):



class RemovePragmas(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        pragmas_map = {}
        for pragma in FindNodes(Pragma).visit(routine.body):
            #if (pragma.keyword.lower() == 'acdc'):
            pragmas_map[pragma] = None
        routine.body = Transformer(pragmas_map).visit(routine.body)
        pragmas_map = {}
        for pragma in FindNodes(Pragma).visit(routine.spec):
            #if (pragma.keyword.lower() == 'acdc'):
            pragmas_map[pragma] = None
        routine.spec = Transformer(pragmas_map).visit(routine.spec)


class RemovePragmaRegions(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        pragmas_map = {}
        for region in FindNodes(PragmaRegion).visit(routine.body):
            pragmas_map[region] = region.body
        routine.body = Transformer(pragmas_map).visit(routine.body)
        pragmas_map = {}
        for region in FindNodes(PragmaRegion).visit(routine.spec):
            pragmas_map[region] = region.spec
        routine.spec = Transformer(pragmas_map).visit(routine.spec)

class ReplaceAbortRegions(Transformation):
    def __init__(self, abort_call = True):
        self.abort_call = abort_call 
    def transform_subroutine(self, routine, **kwargs):
        regions_map = {}
        for region in FindNodes(PragmaRegion).visit(routine.body):
            if ('ABORT' in region.pragma.content):
                # If abort_call is set to false, we will blindly erase the region
                if not self.abort_call :
                    regions_map[region] = None
                # If the ABORT directive has a KEEPME clause, we keep itscontent
                elif ('KEEPME' not in region.pragma.content):
                    # Otherwise, we replace the region with an ABOR1 call
                    regions_map[region] = (CallStatement(name = DeferredTypeSymbol(name='ABOR1'), 
                                                        arguments=(StringLiteral('ERROR : WRONG SETTINGS')), 
                                                        scope=routine), 
                                           )

        routine.body = Transformer(regions_map).visit(routine.body)


class RemoveAssignments(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        assigns_map = {}
        for assign in FindNodes(Assignment).visit(routine.body):
            assigns_map[assign] = None
        routine.body = Transformer(assigns_map).visit(routine.body)
 

class AddSuffixToCalls(Transformation):
    def __init__(self, suffix, additional_variables = [], additional_kwvariables = [], custom_visitor = None):
        self.suffix = suffix
        self.additional_variables = additional_variables
        self.additional_kwvariables = additional_kwvariables
        self.custom_visitor = custom_visitor
        self.routines_called = set()
 
    def transform_node(self, node, routine, inplace = False):
        containedNames = []

        if routine.contains :
            for n in routine.contains.body:
                if isinstance(n, Subroutine):
                    containedNames.append(n.name)

        visitor = self.custom_visitor if self.custom_visitor else FindNodes

        calls_names = set()
        call_map = {}
        for call in visitor(CallStatement).visit(node.body):
            if (call.name not in params.ignored_subroutines):
                if call.name == 'DR_HOOK':
                    list_arguments = list(call.arguments)
                    list_arguments[0] = StringLiteral(list_arguments[0].value + self.suffix)
                    call_map[call] = call.clone(arguments = tuple(list_arguments) )
                elif call.name not in containedNames:
                    print("call foudn ", call.name)
                    new_args = call.arguments
                    new_kwargs = call.kwarguments

                    for var in self.additional_variables:
                        new_args += (Variable(name=var),)
                    for var_couple in self.additional_kwvariables:
                        new_kwargs += ((var_couple[0],Variable(name=var_couple[1])),)

                    #print("call kwargs : ", call.kwarguments)

                    new_call_name = call.name.name + self.suffix
                    self.routines_called.add(call.name.name) 
                    call_map[call] = call.clone(name=DeferredTypeSymbol(name=new_call_name), 
                                        arguments = new_args, kwarguments = new_kwargs)
                    calls_names.add(new_call_name.lower())

        # Update includes : check if they are not already imported
        treated_calls = []
        for imp in FindNodes(Import).visit(routine.spec):
            if imp.c_import:
                for call in calls_names :
                    if imp.module[:-8] == call:
                        treated_calls.append(call)

        for call in calls_names:
            if call not in treated_calls:
                routine.spec.append(Import(module=call + '.intfb.h', c_import = True))


        if inplace:
            node.body = Transformer(call_map).visit(node.body)
            return None
        else:
            new_node = node.clone(body = Transformer(call_map).visit(node.body) )
            return new_node


    def transform_subroutine(self, routine, **kwargs):
        self.transform_node(routine, routine, inplace=True)
            
class RemoveUnusedImports(Transformation):
    # Remove imports for subroutine that are not called 
    # generally needed after AddSuffixToCalls transformation
    def transform_subroutine(self, routine, **kwargs):
        calls_names = set()
        for call in FindNodes(CallStatement).visit(routine.body):
            calls_names.add(call.name.name.lower())
    
        imports_map={}
        for imp in FindNodes(Import).visit(routine.spec):
            if imp.c_import:
                if imp.module[:-8] not in calls_names:
                    imports_map[imp] = None
    
        routine.spec = Transformer(imports_map).visit(routine.spec)

class AddACCRoutineDirectives(Transformation):
    def __init__(self, nodes_list = None):
        self.nodes_list = nodes_list
 
    def transform_subroutine(self, routine,**kwargs):
        containedNames = []

        if routine.contains :
            for n in routine.contains.body:
                if isinstance(n, Subroutine):
                    containedNames.append(n.name)
        # print(f'contained names : {containedNames}')

        bodies_list = self.nodes_list if self.nodes_list else [routine.body]
        calls_names = set()
        call_map ={}
        for body in bodies_list:
            for call in FindNodes(CallStatement).visit(body):
                if call.name.name != 'DR_HOOK':
                    calls_names.add(call.name.name)

        for call in calls_names :
            routine.spec.append(Pragma(keyword='acc', content='routine ('+call+') seq'))
    




class InlineMemberCalls(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        # A call statement inside an inline conditional might become illegal after inlining of the call
        # if the inline body contains more than one line. Therefore we turn those condition into complete form
        for cond in FindNodes(Conditional).visit(routine.body):
            if cond.inline:
                if isinstance(cond.body[0], CallStatement):
                    cond._update(inline=False)
        routine.enrich(routine.members) 
        
        with Timer(logger=info, text='[ACDC] inline_member_procedures {:.2f}s'):
            inline_member_procedures(routine)
        

        with Timer(logger=info, text='[ACDC] search false optionals after inlining in {:.2f}s'):
            # After inlining, optional arguments might get illegally inserted in the routine.
            # However usage of such arguments are expected to happen inside a PRESENT(arg) condition.
            # These PRESENT conditions are transformed to .false. by the inlining procedure, we have to remove those blocks.
            cond_map={}
            for cond in FindNodes(Conditional).visit(routine.body):
                false_cond = False
                true_cond = False
                if (cond.condition == LogicLiteral(False) or 
                        symbolic.simplify(cond.condition, symbolic.Simplification.LogicEvaluation) == LogicLiteral(False)) :
                    
                    false_cond = True
                elif (cond.condition == LogicLiteral(True) or
                        symbolic.simplify(cond.condition, symbolic.Simplification.LogicEvaluation) == LogicLiteral(True)) :

                    true_cond = True
                    

                if false_cond:
                    cond_map[cond] = cond.else_body
                elif true_cond:
                    cond_map[cond] = cond.body

        
        with Timer(logger=info, text=f'[ACDC] apply transformer size {len(cond_map)}' + ' in {:.2f}s'):
            routine.body = Transformer(cond_map).visit(routine.body)


class RemoveComments(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        comments_map = {}
        for comment in FindNodes((Comment, CommentBlock)).visit(routine.body):
            comments_map[comment] = None
        routine.body = Transformer(comments_map).visit(routine.body)
        comments_map = {}
        for comment in FindNodes((Comment, CommentBlock)).visit(routine.spec):
            comments_map[comment] = None
        routine.spec = Transformer(comments_map).visit(routine.spec)
        # Remove comments starting before the first actual declaration
        routine.docstring = None



class RemoveLoops(Transformation):
    def __init__(self, indices=[]):
        for i in indices:
            if not isinstance(i, str):
                error(f'RemoveLoops should be initiliazed with a list of strings')
                exit(1)
        self.indices = indices

    def transform_subroutine(self, routine, **kwargs):
        loop_map = {}
        for loop in FindNodes(Loop).visit(routine.body):
            if (self.indices and loop.variable.name in self.indices) or (not self.indices) :
                loop_map[loop] = loop.body
        routine.body = Transformer(loop_map).visit(routine.body)

class RemoveEmptyConditionals(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        # Cleanup might leave ocnditional with empty ELSEIF which break reconstruction
        # Simply removing all elseif attributes prevents this
        for cond in FindNodes(Conditional).visit(routine.body):
            cond._update(has_elseif = False)
        cond_map = {}
        check = True
        while check:
            for cond in FindNodes(Conditional).visit(routine.body): 
                if not cond.body and not cond.else_body :
                    cond_map[cond] = None 
            if (cond_map == {}):
                check = False
            else :
                routine.body=NestedTransformer(cond_map).visit(routine.body)
                cond_map = {}


class RemoveUnusedVariables(Transformation):
    def __init__(self, symbols_to_remove=[]):
        self.used_symbols = ()
        self.symbols_to_remove = symbols_to_remove

    def transform_subroutine(self, routine, **kwargs):
        # Find called subroutine names as they end up as deferredtype symbols in FindVariables
        calls_names = [call.name.name for call in FindNodes(CallStatement).visit(routine.body)]

        body_variables = set()
        for var in FindVariables().visit(routine.body):
            if var.name not in calls_names:
                while var.parent:
                    var = var.parent
                body_variables.add(var.name)

        declarations_map = {}
        for decl in FindNodes(VariableDeclaration).visit(routine.spec):
            decl_used_symbols = ()

            for s in decl.symbols:
                if s.name in body_variables :
                    if s.name not in self.symbols_to_remove :
                        new_symbol = s.clone(type=s.type.clone(intent='inout'))
                        self.used_symbols += (new_symbol,)
                        decl_used_symbols += (new_symbol,)
                    else :
                        decl_used_symbols += (s,)

            if decl_used_symbols:
                declarations_map[decl] = decl.clone(symbols=decl_used_symbols)
            else:
                declarations_map[decl] = None

        routine.spec = Transformer(declarations_map).visit(routine.spec)
        routine.arguments = self.used_symbols        

