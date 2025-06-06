from loki import (Frontend, backend, Sourcefile, FindNodes, FindExpressions, Loop, Node, Intrinsic, Subroutine, Transformer, PragmaRegion, symbols, Variable )

from loki.ir import Section, Comment,  VariableDeclaration, Pragma, FindVariables, SubstituteExpressions, PreprocessorDirective

source = Sourcefile.from_file('./test.F90', frontend=Frontend.FP)



for routine in source.subroutines:
    zeu_map = {}
    for decls in FindNodes(VariableDeclaration).visit(routine.spec):
        for s in decls.symbols:
            #print(" s : ", s, type(s))
            if isinstance(s, symbols.Array) : 
                print(s)
                print(s.name)
                print(s.type) 
                #izeu_map[decls] = symbols.InlineCall(function=symbols.DeferredTypeSymbol(name='TUTU'))
                #zeu_map[decls] = decls.clone(symbols=decls.symbols.clone(type=None), kind=None, dimensions = ())

                routine.body.append(symbols.InlineCall(function=symbols.DeferredTypeSymbol(name='temp')))
                my_inlinecall = symbols.InlineCall(function=symbols.DeferredTypeSymbol(name='temp'),
                                        #parameters = (symbols.LogicalNot(child = Variable(name='TITI')),)    
                                        parameters = (s.type.kind, s.shape, s.name)
                                        )
                #routine.spec.prepend(VariableDeclaration(name='temp'))
                #print("shpae : ", s.shape, dir(s.shape))
                #zeu_map[decls] = Intrinsic(text='temp'+s.shape.to_str())
                params = ""
                for shape in s.shape:
                    params += backend.fgen(shape) + ","
                    print("shape : ", shape, shape.__str__)
                    #if isinstance(shape, symbols.InlineCall):
                        #print("params ", shape.parameters)
                        #el_not = shape.parameters[2]
                #zi_params = backend.fgen(s.shape) 
                print("zi params : ", params)
                #routine.spec.append(PreprocessorDirective(text='temp'+str(params)))
                #routine.spec.append(PreprocessorDirective(text=str(cohl)))
    
    routine.spec = Transformer(zeu_map).visit(routine.spec)

    zi_code = backend.fgen(my_inlinecall)
    print("zi code : ", zi_code)
    routine.spec.prepend(Intrinsic(text=zi_code))
    routine.spec.prepend(PreprocessorDirective(text=zi_code))


#print("el_not : ", el_not, type(el_not))
#exit()
expr_map = {}

for expr in FindExpressions().visit(routine.spec):
    if isinstance(expr, symbols.LogicalNot):
        print("expr found : ", expr, expr.child)
        expr_map[expr] = symbols.LogicalNot(child = Variable(name='TITI'))
        #expr_map[expr] = symbols.IntLiteral(42)

#routine.spec=SubstituteExpressions(expr_map).visit(routine.spec)

#print(routine.spec.children[0])

f = open('./test_mod.F90', 'w')
f.write(routine.to_fortran())
f.write('\n') #Add eol so vim doesn't complain
f.close()
