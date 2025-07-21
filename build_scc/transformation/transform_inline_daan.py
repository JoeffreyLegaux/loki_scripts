#!/usr/bin/env python3

# DE330 project

from loki import Sourcefile, fgen, FindNodes, analyse_dataflow, FindVariables, FindExpressions
from loki import Transformer, SubstituteExpressions
from loki import Assignment, CallStatement, Conditional, VariableDeclaration, DeferredTypeSymbol, Comment, Import
from loki import Scalar, Array, RangeIndex, Sum, Product, IntLiteral, SymbolAttributes, BasicType
from loki import do_resolve_associates
import copy
import sys

def rename_local_variables(parentroutine, subroutine):
    # prepend local variables of a subroutine with the subroutine name
    
    # note: first rename scalars, then rename arrays. Otherwise scalars used as indices won't get renamed.
    
    # find all local variables
    localVariables=[vv.name.lower() for vv in subroutine.variables if vv not in subroutine.arguments and vv.name.lower() not in ['jlon']]

    # rename scalar variables in spec and body
    rename_map={}
    for var in FindVariables().visit(subroutine.spec):
        if var.name.lower() in localVariables and isinstance(var,Scalar):
            rename_map[var]=var.clone(name=subroutine.name+'_'+var.name)
    subroutine.body=SubstituteExpressions(rename_map).visit(subroutine.body)
    subroutine.spec=SubstituteExpressions(rename_map).visit(subroutine.spec)
    
    # rename array variables in spec
    rename_map={}
    for var in FindVariables().visit(subroutine.spec):
        if var.name.lower() in localVariables and isinstance(var,Array):
            rename_map[var]=var.clone(name=subroutine.name+'_'+var.name)
    subroutine.spec=SubstituteExpressions(rename_map).visit(subroutine.spec)

    # rename array variables in body (different map because indices are part of variable)
    rename_map={}
    for var in FindVariables().visit(subroutine.body):
        if var.name.lower() in localVariables and isinstance(var,Array):
            rename_map[var]=var.clone(name=subroutine.name+'_'+var.name)
    subroutine.body=SubstituteExpressions(rename_map).visit(subroutine.body)
    
    # store declarations
    declarations=[]
    localVariables=[vv.name for vv in subroutine.variables if vv not in subroutine.arguments and vv.name.lower() not in ['jlon']]   # recreate list with modified names
    for localVariable in localVariables:
        declarations.append(VariableDeclaration((subroutine.variable_map[localVariable].clone(scope=parentroutine),)))
    return(declarations)

def substitute_arguments(callStmt,subroutine):
    # substitute variables in a subroutine by the arguments passed in a call
    
    # start from a deep copy to avoid modifying the actual subroutine
    subroutine_copy=subroutine.clone(body=copy.deepcopy(subroutine.body))
    
    # check for optional arguments
    for arg in subroutine_copy.arguments:
        if arg.type.optional:
            raise Exception("Optional arguments not supported. Error at "+str(callStmt.source))

    # loop over arguments
    for iArg in range(len(subroutine_copy.argnames)):

        # check if array dimensions match
        if isinstance(callStmt.arguments[iArg],Array):
            if callStmt.arguments[iArg].dimensions:
                nDimsPassed=0
                for dd in callStmt.arguments[iArg].dimensions:
                    if isinstance(dd,RangeIndex):
                        nDimsPassed=nDimsPassed+1
                        if dd.start or dd.stop:
                            raise Exception("Only full slices (:) or scalar indices are allowed when passing an argument. Error at "+str(callStmt.source))
            else:
                nDimsPassed=len(callStmt.arguments[iArg].scope.variable_map[callStmt.arguments[iArg].name].dimensions)
        else:
            nDimsPassed=0
        if isinstance(subroutine_copy.arguments[iArg],Array):
            nDimsExpected=len(subroutine_copy.arguments[iArg].dimensions)
        else:
            nDimsExpected=0
        if nDimsPassed != nDimsExpected:
            raise Exception('Number of dimensions of passed argument '+fgen(callStmt.arguments[iArg])+' does not match number of dimensions of dummy argument '
                            + fgen(subroutine_copy.arguments[iArg]) + '. Error at '+str(callStmt.source))

        # after these checks, fill rename map
        rename_map={}
        for var in FindVariables().visit(subroutine_copy.body):
            if var.name == subroutine_copy.argnames[iArg]:
                dims=()                
                # "merge dimensions" if they are in the call statement
                if isinstance(callStmt.arguments[iArg],Array):
                    if callStmt.arguments[iArg].dimensions:
                        dims=list(callStmt.arguments[iArg].dimensions)

                        if isinstance(var,Array):
                            if var.dimensions:
                                # replace range dimensions by indices provided in contained subroutine
                                jDim=0
                                for iDim in range(len(dims)):
                                    if isinstance(dims[iDim],RangeIndex):
                                        dims[iDim]=var.dimensions[jDim]
                                        jDim=jDim+1
                                    else:
                                        # scalar index
                                        pass
                            else:
                                # no dimensions in subroutine: just keep dimensions from call statement
                                pass
                        else:
                            # not supposed to happen, because passing an Array to a Scalar should have failed before
                            raise Exception("How did this pass the checks?")
                    else:
                        # no dimensions in call: just keep dimensions from subroutine
                        if isinstance(var,Array):
                            if var.dimensions:
                                dims=var.dimensions
                            else:
                                dims=None
                        else:
                            # not supposed to happen, because passing an Array to a Scalar should have failed before
                            raise Exception("How did this pass the checks?")

                    # that's all for an array argument
                    rename_map[var]=var.clone(scope=None,name=callStmt.arguments[iArg].name,dimensions=tuple(dims))
                    
                else:
                    # not passing an array: just substitute the argument (scalar or expression)
                    rename_map[var]=callStmt.arguments[iArg]

        print("rename_map=",rename_map)
        subroutine_copy.body=SubstituteExpressions(rename_map).visit(subroutine_copy.body)
               
    return(subroutine_copy)

def add_inlining_comments(callStmt,subroutine):            
    subroutine.body.prepend([
        Comment(''),
        Comment('!--- INLINING CALL TO '+subroutine.name+' ---'),
        Comment('!  original call statement: '+fgen(callStmt,linewidth=8*1024)),
        Comment('')                
    ])
    subroutine.body.append([
        Comment(''),
        Comment('!--- END OF INLINING CALL TO '+subroutine.name+' ---'),
        Comment('')
    ])

def inline_contained_subroutines(routine):
    routine_inlined=routine.clone(contains=None)  # contains section already removed
    analyse_dataflow.attach_dataflow_analysis(routine)

    # loop over contained subroutines
    for containedSubroutine in routine.subroutines:
        # move local module imports to calling routine
        routine_inlined.spec.insert(0,FindNodes(Import).visit(containedSubroutine.spec))
        
                
        # rename local variables to avoid conflicts with variables from calling routine
        # also move declarations to calling routine
        routine_inlined.spec.append(rename_local_variables(routine,containedSubroutine))

        # find calls to this subroutine
        containedSubroutineCalls=[cc for cc in FindNodes(CallStatement).visit(routine.body) if str(cc.name).lower() == containedSubroutine.name.lower()]
        inline_map={}
        for callStmt in containedSubroutineCalls:
            # substitute variables by passed arguments
            inlinedSubroutine=substitute_arguments(callStmt,containedSubroutine)

            # resolve associates
            do_resolve_associates(inlinedSubroutine)
            
            # add comments to indicate inlining region
            add_inlining_comments(callStmt,inlinedSubroutine)
            
            # add to mapper
            inline_map[callStmt]=inlinedSubroutine.body

        routine_inlined.body=Transformer(inline_map).visit(routine_inlined.body)        
    
    return(routine_inlined)

#if __name__ == '__main__':
#    if len(sys.argv) != 3:
#        print(sys.argv)
#        print("USAGE: transform_inline.py [file_in] [file_out]")
#    else:
#        file_in=sys.argv[1]
#        file_out=sys.argv[2]
#        source = Sourcefile.from_file(file_in)
#        routine_orig=source[source.subroutines[0].name]
#        routine=routine_orig.clone()
#        analyse_dataflow.attach_dataflow_analysis(routine)
#        routine_inlined=inline_contained_subroutines(routine)
#        with open(file_out,'w') as fid:
#            print("writing to "+file_out)
#            fid.write(fgen(routine_inlined))
#            print("done.")
#        
