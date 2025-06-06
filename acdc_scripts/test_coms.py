from loki import (Frontend, backend, Sourcefile, FindNodes, FindExpressions, Loop, Node, Intrinsic, Subroutine, Transformer, PragmaRegion, symbols, Variable )

from loki.ir import Section, Comment,  VariableDeclaration, Pragma, FindVariables, SubstituteExpressions, PreprocessorDirective

source = Sourcefile.from_file('./test_comments.F90', frontend=Frontend.FP)



for routine in source.subroutines:

    print(routine.docstring)
    #exit()
   # for node in FindNodesi(Node).visit(routine.spec):
   #     print("node ? ", node)

   # exit()
    ast = routine.ir
    #routine._ast.children[0]=None

    my_map = {}
    my_map[routine.ir[0]] = tuple()

    new_ir = Transformer(my_map).visit(routine.ir)
    #routine.ir = Transformer(my_map).visit(routine.ir)


    for child in ast:
        print("==========")
        if isinstance(child, tuple):
            for stmt in child:
                print(stmt)
                stmt = None
        print(type(child))
        print(child)
        print("++++++++++")
    #print(type(ast))
    #print(dir(ast))
#    print(type(routine.ir))
    #ir_list = list(routine.ir)
    #ir_list.pop(0)
    new_routine = routine.clone()
    #new_routine = tuple(ir_list)
    new_routine.docstring = None


f = open('./test_mod.F90', 'w')
f.write(new_routine.to_fortran())
f.write('\n') #Add eol so vim doesn't complain
f.close()
