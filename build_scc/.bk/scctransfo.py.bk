from loki import (
    Sourcefile, FindNodes, CallStatement, 
    Transformer, Dimension, ir, 
    Scalar, Assignment, fgen,
    FindVariables, symbols, demote_variables,
    Intrinsic, Variable, SymbolAttributes,
    DerivedType, VariableDeclaration, flatten,
    BasicType, FindInlineCalls, SubstituteExpressions,
    Nullify, analyse_dataflow, Comparison, RangeIndex, 
)

from loki.transformations.sanitise import resolve_associates

import os
import sys
import copy

from termcolor import colored
import logical
import ResolveVector
import ExplicitArraySyntaxes
import transform_inline_daan as inline_daan #daan inlining
import transform_inline_rolf as inline_rolf #rolf inlining

import re

from pathlib import Path

def load_subroutine(path, file, name):
    source=Sourcefile.from_file(path+file)
    return(source[name])

def save_subroutine(path, file):
    from pathlib import Path
    Sourcefile.to_file(fgen(routine), Path(path+name+'.F90'))


#*********************************************************
#*********************************************************
#*********************************************************
#       Some  routines  of  the  transformation
#*********************************************************
#*********************************************************
#*********************************************************




def add_openacc(routine):
    """
    Replace CALL CALLEE by CALL CALLEE_OPENACC, and CALL ABOR1 with CALL ABOR1_ACC
    :param routine:.
    """
    call_lst=[]
    suffix='_OPENACC'
    stack_argument=Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK"), intent='in'))
    stack_local=Variable(name="YLSTACK", type=SymbolAttributes(DerivedType(name="STACK")), scope=routine)
    for call in FindNodes(CallStatement).visit(routine.body):
        if call.name == "ABOR1" or call.name == "ABOR1_ACC":
            
            #add _ACC to CALL ABOR1 when not already added
            if call.name != "ABOR1_ACC":
                call.name.name=call.name.name+"_ACC"
#            call._update(kwarguments=call.kwarguments + ((stack_argument.name, stack_local),))

	#when the callee isn't ABOR1 or DR_HOOK : save call name, add _OPENACC to the call, and the stack to the args
        elif call.name != "DR_HOOK":
            call_lst.append(call.name.name) #save the call name in order to know which include to update with _openacc
            call.name.name=call.name.name+suffix
            call._update(kwarguments=call.kwarguments + ((stack_argument.name, stack_local),))

    #adding _openacc to the includes according to the lst created just above
    routine_importSPEC=FindNodes(ir.Import).visit(routine.spec)
    for imp in routine_importSPEC:
        name=imp.module.replace(".intfb.h","").upper()
        if imp.c_import==True and any(call==name for call in call_lst):
            new_name=imp.module.replace(".intfb.h","")+"_openacc"+".intfb.h"
            imp._update(module=f'{new_name}') #can also use map 

def remove_horizontal_loop(routine, lst_horizontal_idx):
    """
    Remove all the loops having their index in lst_horizontal_idx 
    :param routine:.
    :param lst_horizontal_idx: a list of given possible horizontal idx
    """
    loop_map={}
    for loop in FindNodes(ir.Loop).visit(routine.body):
        if loop.variable in lst_horizontal_idx:
            loop_map[loop]=loop.body
    routine.body=Transformer(loop_map).visit(routine.body)

def rename(routine):
    """
    Add _OPENACC to the name of the caller
    :param routine:.
    """
    routine.name=routine.name+'_OPENACC'

def acc_seq(routine):
    """
    Add acc seq to the routine
    :param routine:.
    """
    routine.spec.insert(0,ir.Pragma(keyword='acc', content='routine ('+routine.name+') seq'))
    routine.spec.insert(1,ir.Comment(text=''))

def jlon_kidia(routine, end_index, begin_index, new_range, horizontal_idx):
    """
    Add JLON=KIDIA to the routine
    :param routine:.
    :param end_idx: upper bond of the horizontal loop
    :param begin_idx: lower bond of the horizontal loop
    :param new_range: ???
    :param horizontal_idx: horizontal loop variable
    """
    routine.spec.append(Assignment(horizontal_idx, begin_index))

def stack_mod(routine):
    """
    Add before first variable declaration statement : USE STACK_MOD and #include stack.h
    :param routine:.
    """
    idx=-1 #look for idx(-1) of first variable declaration 
    for spec in routine.spec.body:
#        if type(spec)==ir.Intrinsic and spec.text=='IMPLICIT NONE':
#            break
        if type(spec)==ir.VariableDeclaration:
            break
        idx=idx+1
    routine.spec.insert(idx-1, ir.Import(module='STACK_MOD'))
    routine.spec.insert(idx, ir.Import(module='stack.h', c_import=True))
    routine.spec.insert(idx+1, ir.Comment(text=''))

def demote_horizontal(routine, horizontal_size):
    """
    Remove horizontal dimension from local (or not only loc according to demote_arg value) arrays that have the horizontal dim as only dimension and that aren't passed to other routines
    1) search variable
    2) call demote_variables
    :param routine:.
    :param horizontal_size: size of the horizontal dimension (NPROMA)
    """
    #demote_arg=True #rm KLON in function arg if True
    demote_arg=False #rm KLON in function arg if True
    routine_arg=[var.name for var in routine.arguments]
    variables=FindVariables(unique=True).visit(routine.spec)
    variables=[var for var in variables if isinstance(var, symbols.Array)]
    to_demote=[]
    for var in variables:
        not_cste=False
        if len(var.shape)==1:
            if var.shape[0] == horizontal_size:
                to_demote.append(var)
        else: 
            for shape in var.shape:
                if not isinstance(shape, symbols.IntLiteral):
                    if shape!=horizontal_size:
                        not_cste=True
                        break
            if not not_cste: #if the shape of the var are cste or horizontal
                to_demote.append(var) 
    
    if not demote_arg :
        to_demote = [var for var in to_demote if var.name not in routine_arg]


        calls = FindNodes(ir.CallStatement).visit(routine.body)
        call_args = flatten(call.arguments for call in calls)
        call_args += flatten(list(dict(call.kwarguments).values()) for call in calls)
        to_demote = [v for v in to_demote if v.name not in call_args]

    var_names=tuple(var.name for var in to_demote)
    #TODO demote over all horizontal dimensions if more than one
    if var_names:
        demote_variables(routine, var_names, dimensions=horizontal_size)

def alloc_temp(routine):
    """
    Replace local array declaration by alloc + temp macro : (the temporaries will be created in a piece of memory allocated outside of the kernel)
    :param routine:.
    """
    routine_arg=[var.name for var in routine.arguments]
    

    temp_map={}
    for decls in FindNodes(VariableDeclaration).visit(routine.spec):
        intrinsic_lst=[] #intrinsic lst for temp macro var decl
        var_lst=[] #var to keep in the decl.
        for s in decls.symbols: #a decl can have >1 var
            if isinstance(s, symbols.Array):
                if not s.type.pointer: #if pt, 
                    if s.name not in routine_arg:
                        not_cste=False
                        for shape in s.shape:

                            if (not isinstance(shape, symbols.IntLiteral)):
                                if isinstance(shape, symbols.Scalar) :
                                    if not shape.type.initial:
                                        not_cste=True
                                        break
                        if not_cste: 

                            if s.type.kind:
                                new_s='temp ('+s.type.dtype.name+' (KIND='+s.type.kind.name+'), '+s.name+', ('
                            else:
                                new_s='temp ('+s.type.dtype.name+', '+s.name+', ('
                            for shape in s.shape:
                                new_s=new_s+str(shape)+', '
                            new_s=new_s[:-2]
                            new_s=new_s+'))'
                            kind = symbols.DeferredTypeSymbol(name=f"KIND ({s.name})")
                            #kind = symbols.InlineCall(symbols.Variable(name='KIND'), parameters=(new_s,))
                            
                            cond8 = Comparison(
                            		left=kind,
                            		operator='==',
                            		right=symbols.IntLiteral(value=8))
                            cond4 = Comparison(
                            		left=kind,
                            		operator='==',
                            		right=symbols.IntLiteral(value=4))
                            alloc8='alloc8 ('+s.name+')'
                            alloc4='alloc4 ('+s.name+')'
                            stop='STOP 1'
                            else_cond = ir.Conditional(
                            		condition=cond4,
                            		body=(Intrinsic(alloc4),),
                            		else_body=(Intrinsic(stop),)
                                        )
                            alloc_block = ir.Conditional(
                            	condition=cond8,
                            	body=(Intrinsic(alloc8),),
                            	else_body=(else_cond,),
                            	hase_elseif=True
                            		)
                            routine.spec.append(alloc_block)
                            intrinsic_lst.append(Intrinsic(new_s))
                        else: 
                            var_lst=[s]
                            VAR=decls.clone(symbols=var_lst)
                            intrinsic_lst.append(VAR)
                    else: #if array in routine args
#                        var_lst.append(decls.clone(symbols=s))
                        var_lst=[s]
                        VAR=decls.clone(symbols=var_lst)
                        intrinsic_lst.append(VAR)
                else: #if s is a pointer
                    #print(colored("POINTER","red"))
                    var_lst=[s]
                    VAR=decls.clone(symbols=var_lst)
                    intrinsic_lst.append(VAR)
            else: #if not an array
#                var_lst.append(decls.clone(symbols=s))
                var_lst=[s]
                VAR=decls.clone(symbols=var_lst)
                intrinsic_lst.append(VAR)

        temp_map[decls]=tuple(intrinsic_lst)
#        VAR=decls.clone(symbols=var_lst)
#        intrinsic_lst.append(VAR)

    routine.spec=Transformer(temp_map).visit(routine.spec)
    

def get_horizontal_size(routine, lst_horizontal_size):
    """
    Look which alias of NPROMA is currently being used in the routine.
    :param routine:. 
    :param lst_horizontal_size: list of aliases of NPROMA, you may add new ones
    """
    #TODO: fix
    verbose=False
    #verbose=True
    #Looks for hor size in the routine variables, or if in derived type
    for name in lst_horizontal_size:
        if name in routine.variable_map:
            if verbose: print("horizontal size = ",name)
            return(name)
    for var in FindVariables().visit(routine.body):
        for vvar in var.name.split("%"):
            if vvar in lst_horizontal_size:
                if verbose: print(colored("Horizontal size in a derived type", "red"))
                if verbose: print("horizontal size = ",name)
                return(var.name)
    if verbose: print(colored("Horizontal size not found in routine args!", "red"))
#    if verbose: print("horizontal size = ", horizontal.size)
    
    
    #Looks for hor size if used in an array without being declared (--> TODO : fix/check this)
    hor_lst=[]
    var_lst=FindVariables(unique=True).visit(routine.spec)
    var_lst=[var for var in var_lst if isinstance(var, symbols.Array)]
    for var in var_lst:
        for shape in var.shape: 
            if shape in lst_horizontal_size: 
                if shape not in hor_lst:
                    hor_lst.append(shape)
    if len(hor_lst)>1:
        print(colored("diff horizontal size are used", "red"))
        return(horizontal.size)
    else:
        return(hor_lst[0])


def ystack1(routine):
    """
    Add STACK in the routine args.
    :param routine:.
    """
#    stack_argument=Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK")),scope=routine)
    stack_argument=Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK"), intent='in'))
    stack_local=Variable(name="YLSTACK", type=SymbolAttributes(DerivedType(name="STACK")), scope=routine)

    routine.arguments+=(stack_argument,)

def ystack2(routine):
    """
    Add STACK in the routine args.
    :param routine:.
    """
#    stack_argument=Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK")),scope=routine)
    stack_argument=Variable(name="YDSTACK", type=SymbolAttributes(DerivedType(name="STACK"), intent='in'))
    stack_local=Variable(name="YLSTACK", type=SymbolAttributes(DerivedType(name="STACK")), scope=routine)

    routine.variables+=(stack_argument, stack_local,)
    routine.spec.append(Assignment(stack_local, stack_argument))


def get_horizontal_idx(routine, lst_horizontal_idx):
    """
    Add and replace horizontal loop index by JLON if not already present.
    :param routine:.
    :param lst_horizontal_idx: a list of given possible horizontal idx 
    """

    is_present=False
    verbose=False
    var_lst=FindVariables(unique=True).visit(routine.spec)
    for var in var_lst:
        if var.name=="JLON":
            is_present=True
            loop_index=var
    
    if not is_present:
        jlon=Variable(name="JLON", type=SymbolAttributes(BasicType.INTEGER, kind=Variable(name='    JPIM')))
        loop_index=jlon
        routine.variables+=(jlon,)

    loop_map={}
#    for loop in FindNodes(ir.Loop).visit(routine.body):
#        if loop.variable in lst_horizontal_idx:
#            #new_loop=loop.clone(variable=routine.variable_map["JLON"])
#            if verbose: print("loop_idx=",loop_index)
#            new_loop=loop.clone(variable=routine.variable_map[loop_index.name])
#            var_map={}
#            for var in FindVariables().visit(loop.body):
#                if (var==loop.variable):
#                    var_map[var]=routine.variable_map[loop_index.name]
#            loop_map[loop]=SubstituteExpressions(var_map).visit(loop.body)
#            #loop.body=SubstituteExpressions(var_map).visit(loop.body)
#    #        loop.variable=routine.variable_map[loop_index.name]
#    
#    routine.body=Transformer(loop_map).visit(routine.body)
 
    #TODO : Change this part; if horizontal loop idx # in lst_horizontal_idx...
    #1) check bounds of the loop to try to get hor loop
    #2) change index everywhere cf "rename_hor" function.
    rename_map={}
    for var in FindVariables().visit(routine.body):
        if var.name in lst_horizontal_idx:
            rename_map[var]=loop_index
    routine.body=SubstituteExpressions(rename_map).visit(routine.body)

    return(loop_index)


def rm_sum(routine, horizontal_size):
    """
    Remove calls to SUM if the sum is performed on the horizontal dimension.

    TODO : add check to be sure the SUM is of this type : 
    ZTEST = SUM(A(:))
    IF (ZTEST > 0.0_JPRB) THEN
    ...

    :param routine:.
    """
    verbose = False
    #verbose = True
    call_map={}
    for assign in FindNodes(Assignment).visit(routine.body):
        for call in FindInlineCalls().visit(assign):
            if (call.name=="SUM"):
                if len(call.parameters)>1:
                    raise NotImplementedError("SUM should have only one arg")
                else:
                    arg = call.parameters[0]
                    if len(arg.shape)==1:
                        if arg.shape[0].name == horizontal_size:
                            if verbose : print(f"Removing an horizontal SUM : {fgen(assign)}")
                            call_map[call]=call.parameters[0]
                    else:
                        shape_idx = 0
                        for shape in arg.shape:
                            if arg.shape[0].name == horizontal_size:
                                break 
                            shape_idx += 1
                        if isinstance(arg.dimensions[shape_idx], RangeIndex):
                            if verbose : print(f"Removing an horizontal SUM : {fgen(assign)}")
                            call_map[call]=call.parameters[0]

    if call_map:
        routine.body=SubstituteExpressions(call_map).visit(routine.body)

def generate_interface(routine):
    """
    Generate the interface of the new routine.
    :param routine:. 
    """
    removal_map={}
    imports = FindNodes(ir.Import).visit(routine.spec)
    routine_new=routine.clone()
    #import must be removed 
    for im in imports:
        if im.c_import==True:
            removal_map[im]=None
    routine_new.spec = Transformer(removal_map).visit(routine_new.spec)
    #Sourcefile.to_file(fgen(routine_new.interface), Path(pathw+".intfb.h"))
    str_interface = fgen(routine_new.interface)
    return(str_interface)

def write_print(routine):
    """
    Change WRITE statements into PRINT statements.
    :param routine:. 
    """
    verbose=False
    intr_map={}
    for intr in FindNodes(Intrinsic).visit(routine.body):
        if "WRITE" in intr.text:
            if verbose : print(colored("WRITE found in routine","red"))
            pattern='WRITE\(NULERR, \*\)(.*)'
            match=re.search(pattern,intr.text)
            string=match.group(1)
            new_intr='PRINT *, '+string
            intr_map[intr]=Intrinsic(new_intr)
    routine.body=Transformer(intr_map).visit(routine.body)
            
    
#---------------------------------------------------------
#---------------------------------------------------------
#Pointers
#---------------------------------------------------------
#---------------------------------------------------------

def find_pt(routine):
    """
    Return lst of pointers and targets
    :param routine:.
    """
    #TODO : mv if target.dimensions[0]=='KPROMA' here...
    tmp_pt=[]
    tmp_target=[]
    #for var in FindVariables().visit(routine.variables):
    routine_args=[var.name for var in routine.arguments]
    for decls in FindNodes(VariableDeclaration).visit(routine.spec):
        for s in decls.symbols:
            if s.name not in routine_args: #we are looking only for loc pointers
                if s.type.pointer:
                    tmp_pt.append(s)
            if s.type.target:
                tmp_target.append(s)
    return(tmp_pt, tmp_target)
             
             
def get_dim_pt(routine, tmp_pt, tmp_target, horizontal_size):
    """
    Replace A=B by assoc(A,B). Where A : pointer and B : target.
    The first dimension of B must be the horizontal dim, and the size must be NPROMA. If there isn't any size, the size is assume to be # from NPROMA. 
    :param routine:. 
    :param tmp_pt: lst of temporary pointers declared in the subroutine
    :param tmp_pt: lst of temporary targets declared in the subroutine
    :param horizontal_size: size of the horizontal dimension (NPROMA)
    """
    #TODO : replace KPROMA by horizontal_size...
    assign_map={}
    tmp_pt_klon={} #dict associating pt and target 
    for assign in FindNodes(Assignment).visit(routine.body):
        assigned=False
        is_target=False
        is_pt=False
        rhs=FindVariables().visit(assign.rhs)
        rhs=list(rhs)
        lhs=FindVariables().visit(assign.lhs)
        lhs=list(lhs)
        if (len(rhs)==1 and len(lhs)==1):
            for pt in tmp_pt:
                if pt.name==lhs[0].name:
                   is_pt=True
                   break
            for target in tmp_target:
                if target.name==rhs[0].name:
                   is_target=True
                   break
            if is_pt and is_target:
#TODO: mv target.dimensions[0]=='KPROMA' in find_pt
                if target.dimensions[0]=='KPROMA': #replace KPROMA by hor_dim... here we assume that if they are no dim, then the dim isn't NPROMA.... 
                    tmp_pt_klon[lhs[0].name]=rhs[0]
                    new_assign='assoc ('+lhs[0].name+','+rhs[0].name+')'
                    assign_map[assign]=Intrinsic(new_assign)
                    assigned=True
        if not assigned:
            assign_map[assign]=assign
    routine.body=Transformer(assign_map).visit(routine.body) 
    return(tmp_pt_klon)

def nullify(routine, tmp_pt_klon):
    """
    Changes NULLIFY into nullptr macro 
    :param routine:.
    :param tmp_pt_klon: dict associating pt and target (only targets with horizontal dim of KPROMA)
    """
#    null=FindNodes(Nullify).visit(routine.body)
    null_map={}
    for null in FindNodes(Nullify).visit(routine.body):
        new_null=()
        for pt in null.variables:
            if null.variables[0].name in tmp_pt_klon:
                new_intr=Intrinsic("nullptr("+null.variables[0].name+")")
                new_null=new_null+(new_intr,)
        null_map[null]=new_null
    routine.body=Transformer(null_map).visit(routine.body)
#define target lst to avoid A.type => bug



def assoc_alloc_pt(routine, tmp_pt_klon):
    """
    Changes targets decl in tmp_pt_klon using the temp macro.
    :param routine:.
    :param tmp_pt_klon: dict associating pt and target (only targets with horizontal dim of KPROMA)
    """
    #TODO: create a func to replace allo by temp macro
    temp_map={}
    for decls in FindNodes(VariableDeclaration).visit(routine.spec):
        intrinsic_lst=[]
        var_lst=[]
        for s in decls.symbols:
            if s.name in tmp_pt_klon:
                var=tmp_pt_klon[s.name]
                if var.type.kind:
                    new_s='temp ('+var.type.dtype.name+' (KIND='+var.type.kind.name+'), '+s.name+', ('
                else:
                    new_s='temp ('+var.type.dtype.name+', '+s.name+', ('
                for shape in var.shape:
                    new_s=new_s+str(shape)+', '
                new_s=new_s[:-2]
                new_s=new_s+'))'
                #alloc='alloc ('+s.name+')'
                #routine.spec.append(Intrinsic(alloc))
                intrinsic_lst.append(Intrinsic(new_s))
            else:
                var_lst=[s]
                VAR=decls.clone(symbols=var_lst)
                intrinsic_lst.append(VAR)
        temp_map[decls]=tuple(intrinsic_lst)
    routine.spec=Transformer(temp_map).visit(routine.spec)


     
#*********************************************************
#*********************************************************
#*********************************************************
#       Defining     the       transformation
#*********************************************************
#*********************************************************
#*********************************************************

#def openacc_trans(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined):
def openacc_trans(routine, horizontal, lst_horizontal_idx, lst_horizontal_size, lst_horizontal_bounds, true_symbols, false_symbols, inline_match):
    """
    Meteo France SCC transformation

    """

    #----------------------------------------------
    #setup
    #----------------------------------------------
    #verbose=True
    verbose=False
    
    #----------------------------------------------
    #Inlining 
    #----------------------------------------------

    from loki.transformations import inline_member_procedures
    import transform_inline as tf_in
    
    #inline_method='daan'
    inline_method='ec'
    #inline_method='rolf'
    if inline_match:
        if inline_method=='ec':
            mv_include(routine) #ec trans pb with include
#            routine.enrich_calls(routine.members)
            tf_in.inline_member_procedures(routine)
            rename_hor(routine, lst_horizontal_idx)
        elif inline_method=='daan':
#            routine_orig=source.subroutines[0]
#            routine=routine_orig.clone()
            analyse_dataflow.attach_dataflow_analysis(routine)
            routine_inlined=inline_daan.inline_contained_subroutines(routine)
            routine=routine_inlined
        elif inline_method=='rolf':
            mv_include(routine) #ec trans pb with include
#            routine.enrich_calls(routine.members)
            inline_rolf.inline_member_procedures(routine)
            rename_hor(routine, lst_horizontal_idx)

    #----------------------------------------------
    #transformations:
    #----------------------------------------------
    horizontal_size=get_horizontal_size(routine, lst_horizontal_size)
    resolve_associates(routine)
    logical.transform_subroutine(routine, true_symbols, false_symbols)
    
    end_index, begin_index, new_range=ExplicitArraySyntaxes.ExplicitArraySyntaxes(routine, lst_horizontal_size, lst_horizontal_bounds)
    horizontal_idx=get_horizontal_idx(routine, lst_horizontal_idx)
    bounds=[begin_index, end_index]
    add_openacc(routine)
    
    rename(routine)
    acc_seq(routine)
    stack_mod(routine)
    rm_sum(routine, horizontal_size) #must be before demote_horizontal
    demote_horizontal(routine, horizontal_size)
    #TODO : change resolve_vector : changes dim for all possible idx and bounds..
    ResolveVector.resolve_vector_dimension(routine, horizontal_idx, bounds)
    
    remove_horizontal_loop(routine, lst_horizontal_idx)
    ###
    ystack1(routine)
    #rm_sum(routine) #must be before demote_horizontal
    str_interface = generate_interface(routine) #must be before temp allocation and y stack, or temp will be in interface
   
   ##----------------------------------------
   ##Pointers
   ##----------------------------------------
    tmp_pt, tmp_target=find_pt(routine)
    tmp_pt_klon=get_dim_pt(routine, tmp_pt, tmp_target, horizontal_size)
    assoc_alloc_pt(routine, tmp_pt_klon)
    nullify(routine, tmp_pt_klon)
    ##----------------------------------------
    ##----------------------------------------
    write_print(routine)
    
    ystack2(routine)
    alloc_temp(routine) #must be after demote_horizontal
    jlon_kidia(routine, end_index, begin_index, new_range, horizontal_idx) #at the end
    
    return(str_interface)
