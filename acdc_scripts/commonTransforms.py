from loki import (Frontend, Sourcefile, FindNodes, Loop, Node, Intrinsic, Subroutine, Transformer, NestedTransformer, 
    PragmaRegion, DerivedType, Transformation, CallStatement, SymbolAttributes, BasicType  )

from loki.ir import Section, Comment, CommentBlock, VariableDeclaration, Pragma, PragmaRegion, Import, Assignment, Conditional, LeafNode, InternalNode, Associate

from loki.transformations import inline_member_procedures

from loki.expression import FindTypedSymbols, FindVariables
from loki.expression.expr_visitors import SubstituteExpressions
from loki.expression.symbols import DeferredTypeSymbol, TypedSymbol, Array, Scalar, RangeIndex, Variable, StringLiteral, InlineCall, LogicLiteral, LiteralList
from loki.expression import symbolic
from loki.frontend.fparser import *
from loki.logging import info, error
from loki.analyse import *

from codetiming import Timer

from storable import retrieve



class RemovePragmas(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        pragmas_map = {}
        for pragma in FindNodes(Pragma).visit(routine.body):
            if (pragma.keyword.lower() == 'acdc'):
                pragmas_map[pragma] = None
        routine.body = Transformer(pragmas_map).visit(routine.body)
        pragmas_map = {}
        for pragma in FindNodes(Pragma).visit(routine.spec):
            if (pragma.keyword.lower() == 'acdc'):
                pragmas_map[pragma] = None
        routine.spec = Transformer(pragmas_map).visit(routine.spec)


class RemovePragmaRegions(Transformation):
    def transform_subroutine(self, routine, **kwargs):
        pragmas_map = {}
        for region in FindNodes(PragmaRegion).visit(routine.body):
            pragmas_map[region] = region.body
        routine.body = Transformer(pragmas_map).visit(routine.body)
        # pragmas_map = {}
        # for region in FindNodes(PragmaRegion).visit(routine.spec):
        #     pragmas_map[region] = region.spec
        # routine.spec = Transformer(pragmas_map).visit(routine.spec)



class AddSuffixToCalls(Transformation):
    def __init__(self, suffix, node = None):
        self.suffix = suffix
        self.node = node
 
    def transform_subroutine(self, routine, **kwargs):
        containedNames = []

        if routine.contains :
            for n in routine.contains.body:
                if isinstance(n, Subroutine):
                    containedNames.append(n.name)

        body = self.node.body if self.node else routine.body

        calls_names = set()
        call_map ={}
        for call in FindNodes(CallStatement).visit(body):
            if call.name == 'DR_HOOK':
                list_arguments = list(call.arguments)
                list_arguments[0] = StringLiteral(list_arguments[0].value + self.suffix)
                call_map[call] = call.clone(arguments = tuple(list_arguments) )
            elif call.name not in containedNames:
                call_map[call] = call.clone(name=DeferredTypeSymbol(name=call.name.name + self.suffix))
                calls_names.add(call.name.name.lower())


        # update includes
        treated_calls = []
        imports_map = {}
        for imp in FindNodes(Import).visit(routine.spec):
            if imp.c_import:
                # print("module : ", imp.module[:-8])
                for call in calls_names :
                    if imp.module[:-8] == call:
                        imports_map[imp] = imp.clone(module=imp.module[:-8] + self.suffix.lower() + '.intfb.h')
                        treated_calls.append(call)
                    # Edge case : if import already changed in previous parallel region treatment
                    if imp.module[:-8] == call + self.suffix.lower():
                        treated_calls.append(call)


        for call in calls_names:
            if call not in treated_calls:
                routine.spec.append(Import(module=call + self.suffix.lower() + '.intfb.h', c_import = True))


        routine.spec = Transformer(imports_map).visit(routine.spec)


        if self.node:
            new_node = self.node.clone(body = Transformer(call_map).visit(self.node.body) )

            new_node2 = new_node.clone()
            routine.body = Transformer({self.node:new_node}).visit(routine.body)
            return new_node2
        else:
            routine.body = Transformer(call_map).visit(routine.body)
            return None
            
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
        cond_map = {}
        for cond in FindNodes(Conditional).visit(routine.body):
            if not cond.body and not cond.else_body :
                cond_map[cond] = None 
        routine.body=Transformer(cond_map).visit(routine.body)

class RemoveUnusedVariables(Transformation):
    def __init__(self):
        self.used_symbols = ()

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
                if s.name in body_variables:
                    new_symbol = s.clone(type=s.type.clone(intent='inout'))
                    self.used_symbols += (new_symbol,)
                    decl_used_symbols += (new_symbol,)

            if decl_used_symbols:
                declarations_map[decl] = decl.clone(symbols=decl_used_symbols)
            else:
                declarations_map[decl] = None

        routine.spec = Transformer(declarations_map).visit(routine.spec)
        routine.arguments = self.used_symbols        

