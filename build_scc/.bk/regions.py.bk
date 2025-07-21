from loki import (
    Sourcefile, FindNodes, CallStatement, 
    Transformer, Dimension, ir, 
    Scalar, Assignment, fgen,
    FindVariables, symbols, demote_variables,
    Intrinsic, Variable, SymbolAttributes,
    DerivedType, VariableDeclaration, flatten,
    BasicType, FindInlineCalls, SubstituteExpressions,
    Nullify, analyse_dataflow
)

from loki.transformations import resolve_associates

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
    routine.body.prepend(Assignment(horizontal_idx, begin_index))

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
    demote_arg=False #rm KLON in function arg if True
    routine_arg=[var.name for var in routine.arguments]
    to_demote=FindVariables(unique=True).visit(routine.spec)
    to_demote=[var for var in to_demote if isinstance(var, symbols.Array)]
    to_demote=[var for var in to_demote if var.shape[-1] == horizontal_size]
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
    
    allocs_list=[]

    temp_map={}
    for decls in FindNodes(VariableDeclaration).visit(routine.spec):
        intrinsic_lst=[] #intrinsic lst for temp macro var decl
        var_lst=[] #var to keep in the decl.
        for s in decls.symbols: #a decl can have >1 var
            if isinstance(s, symbols.Array):
                if not s.type.pointer and not s.type.parameter : #if pt, 
                    if s.name not in routine_arg:
                        if s.type.kind:
                            new_s='temp ('+s.type.dtype.name+' (KIND='+s.type.kind.name+'), '+s.name+', ('
                        else:
                            new_s='temp ('+s.type.dtype.name+', '+s.name+', ('
                        for shape in s.shape:
                            new_s=new_s+str(shape)+', '
                        new_s=new_s[:-2]
                        new_s=new_s+'))'
                        alloc='alloc ('+s.name+')'
                        allocs_list.append(Intrinsic(alloc))
                        intrinsic_lst.append(Intrinsic(new_s))
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

    # Find the YLSTACK = YDSTACK assignement, insert "alloc" macros afterwards
    idx=0
    for node in FindNodes(ir.Node).visit(routine.body):
        if type(node)==ir.Assignment:
            if node.lhs.name == 'YLSTACK':
                break
        idx=idx+1
    routine.body.insert(idx, tuple(allocs_list))
    

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
    routine.body.prepend(Assignment(stack_local, stack_argument))


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


def rm_sum(routine):
    """
    Remove calls to SUM fortran function.
    :param routine:.
    """
    #TODO : Check if one wants to remove it everywhere? 
    call_map={}
    for assign in FindNodes(Assignment).visit(routine.body):
        for call in FindInlineCalls().visit(assign):
            if (call.name=="SUM"):
                call_map[call]=call.parameters[0]
    if call_map:
        routine.body=SubstituteExpressions(call_map).visit(routine.body)

def generate_interface(routine, pathw):
    """
    Generate the interface of the new routine.
    :param routine:. 
    :param pathw: absolute path (path + filename) where the routine_openacc is saved.
    """
    removal_map={}
    imports = FindNodes(ir.Import).visit(routine.spec)
    routine_new=routine.clone()
    #import must be removed 
    for im in imports:
        if im.c_import==True:
            removal_map[im]=None
    routine_new.spec = Transformer(removal_map).visit(routine_new.spec)
    Sourcefile.to_file(fgen(routine_new.interface), Path(pathw+".intfb.h"))

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
#      Functions    needed    for    inlining....
#*********************************************************
#*********************************************************
#*********************************************************

def add_contains(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined):
#def inline_calls(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined):
    """
    Add routine to inline as contained subroutine (loki's inliner can only inline contained subroutines)
    :param pathpack: absolute path to the pack
    :param pathview: path to the src/local/... or src/main/... folder
    :param pathfile: path (+name) to the file
    :param pathacc: path (+name) to the acc file 
    :param horizontal_opt: additional horizontal loop index
    :param inlined: lst of routines to inline
    """
    #verbose=False
    verbose=True
    pathr=pathpack+'/'+pathview+pathfile
    match_inline=False
    dict_callee_path={}
    if inlined != None:    
           
  #creation of a dict associating the name of each callee to inline to it's path; according to the path in the openacc.sh file 
  #TODO : create a file containing these path? 
        with open('/home/gmap/mrpm/cossevine/build_scc/openacc.sh', 'r', encoding='utf-8', errors='ignore') as file_lst_callee:
        #with open('/home/gmap/mrpm/cossevine/build_scc/openacc.sh', 'r') as file_lst_callee:
            lines=file_lst_callee.readlines()
            for callee_name in inlined:
                callee_path=None
                for line in lines:
                    callee_path=re.search('((\w*\/)+'+callee_name+')', line)
                    if callee_path:
                        
                        callee_path=pathpack+'/'+pathview+callee_path.group(0)
                        dict_callee_path[callee_name]=callee_path
        if verbose: print("dict_callee_path=",dict_callee_path)

   #open the current routine "caller", inserts each callee when a CALL CALLEE appears in the caller. 
        #with open(pathr, 'r') as file_caller:
        with open(pathr, 'r', encoding='utf-8', errors='ignore') as file_caller:
            caller = file_caller.read()
                
        for callee_name in inlined: #look for each callee sub  in the caller
            if caller.find("CALL "+callee_name.replace(".F90","").upper())!=-1:
                if verbose: print("callee_name = ", callee_name)
                #with open(dict_callee_path[callee_name], 'r') as file_callee:
                with open(dict_callee_path[callee_name], 'r', encoding='utf-8', errors='ignore') as file_callee:
                    callee = file_callee.read()
                if not match_inline: #add CONTAINS only for the first callee matching
                    match_inline=True
                    callee="CONTAINS\n\n"+callee
                    loc=caller.find("END SUBROUTINE")
                    if loc != -1 :
                        caller=caller[:loc]+callee+caller[loc:]
                else: #add the callee after the CONTAINS statement
                    loc=caller.find("CONTAINS")
                    if loc != -1 :
                         loc=loc+len("CONTAINS\n\n") #insert after CONTAINS
                         caller=caller[:loc]+callee+caller[loc:]

        if verbose: print(pathpack+"/tmp/"+os.path.basename(pathfile))
        if match_inline:
            with open(pathpack+"/tmp/"+os.path.basename(pathfile), "w", encoding='utf-8', errors='ignore') as file_caller:
           # with open(pathpack+"/tmp/"+os.path.basename(pathfile), "w") as file_caller:
                file_caller.write(caller)
    else:
        if verbose: print(colored("no routine to inline", "red"))


    return(match_inline)

#diff way to do : uniformize loop or re.search(_LOOPIDX) 
def rename_hor(routine, lst_horizontal_idx): 
    """
    When the inlining is done the loop index may be CALLEE_LOOPIDX, this routine changes it to LOOPIDX
    :param routine:.
    :param lst_horizontal_idx: list of aliases of NPROMA, you may add new ones
    """
    rename_map={}
    for var in FindVariables().visit(routine.body):
        for idx in lst_horizontal_idx:
            if '_'+idx in var.name:
                rename_map[var]=var.clone(name=idx)
    routine.body=SubstituteExpressions(rename_map).visit(routine.body)
                    
def mv_include(routine):
    """
    Ecmwf inlining function doesn't seem to move the include from the CALLEE to the CALLER : it's done in this function.
    :param routine:.
    """
    imp_lst=[imp.module for imp in FindNodes(ir.Import).visit(routine.spec)]
    for containedSubroutine in routine.subroutines:
        for imp in FindNodes(ir.Import).visit(containedSubroutine.spec):
            if not imp.module in imp_lst:
                imp_lst.append(imp.module)
                routine.spec.append(imp)
                #routine.spec.insert(0,imp)
    


     
#*********************************************************
#*********************************************************
#*********************************************************
#       Defining     the       transformation
#*********************************************************
#*********************************************************
#*********************************************************

def scc_transform_routine(routine, lst_horizontal_size, lst_horizontal_idx, lst_horizontal_bounds, true_symbols, false_symbols, pathw = None):
    #----------------------------------------------
    #transformations:
    #----------------------------------------------
    horizontal_size=get_horizontal_size(routine, lst_horizontal_size)
    resolve_associates(routine)
    logical.transform_subroutine(routine, true_symbols, false_symbols)
    
    end_index, begin_index, new_range=ExplicitArraySyntaxes.ExplicitArraySyntaxes(routine, lst_horizontal_size, lst_horizontal_bounds)
    horizontal_idx=get_horizontal_idx(routine, lst_horizontal_idx)
    bounds=[begin_index, end_index]
    
    
    
    
    stack_mod(routine)
    demote_horizontal(routine, horizontal_size)
    #TODO : change resolve_vector : changes dim for all possible idx and bounds..
    ResolveVector.resolve_vector_dimension(routine, horizontal_idx, bounds)
    
    remove_horizontal_loop(routine, lst_horizontal_idx)
    ###
    ystack1(routine)
    rm_sum(routine)

    if pathw:
        generate_interface(routine, pathw) #must be before temp allocation and y stack, or temp will be in interface

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
    
    jlon_kidia(routine, end_index, begin_index, new_range, horizontal_idx)
    ystack2(routine)
    # alloc_temp(routine)

def openacc_trans(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined):
    """
    Meteo France SCC transformation

    """

    #----------------------------------------------
    #setup
    #----------------------------------------------
    verbose=True
    #verbose=False
    pathr=pathpack+'/'+pathview+pathfile
    pathw=pathpack+pathacc+'/'+pathfile

    pathw=pathw.replace(".F90", "")+"_openacc"
    
    if verbose: print('pathr=', pathr)
    if verbose: print('pathw=', pathw)

    import logical_lst
    
    horizontal=Dimension(name='horizontal',size='KLON',index='JLON',bounds=['KIDIA','KFDIA'],aliases=['NPROMA','KDIM%KLON','D%INIT'])
    vertical=Dimension(name='vertical',size='KLEV',index='JLEV')
    
     #lst_horizontal_idx=['JLON','JROF','JL']
    lst_horizontal_idx=['JLON','JROF']
    #the JL idx have to be added only when it's used at an horizontal idx, because it's used as avertical idx in some places.... this should be fixed... The transformation script could check wether JL is a hor or vert idx instead of adding JL to the lst_horizontal_idx conditionally. 
    if horizontal_opt is not None:
        lst_horizontal_idx.append(horizontal_opt)
    if verbose: print("lst_horizontal_idx=",lst_horizontal_idx)
    
    lst_horizontal_size=["KLON","YDCPG_OPTS%KLON","YDGEOMETRY%YRDIM%NPROMA","KPROMA", "YDDIM%NPROMA", "NPROMA"]
    lst_horizontal_bounds=[["KIDIA", "YDCPG_BNDS%KIDIA","KST"],["KFDIA", "YDCPG_BNDS%KFDIA","KEND"]]
    
    true_symbols, false_symbols=logical_lst.symbols()
    false_symbols.append('LHOOK')
    
    #----------------------------------------------
    #inlining : 
    #----------------------------------------------
    inlined=list(inlined)
    if inlined:
        if verbose: print("****************INLINING****************")
        inline_match=add_contains(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined)
    else:
        inline_match=False

    if inline_match:
        filename=os.path.basename(pathfile)
        source=Sourcefile.from_file(pathpack+"/tmp/"+filename)
    else:
        source=Sourcefile.from_file(pathr)
    routine=source.subroutines[0]

    from loki.transform import inline_member_procedures
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



    add_openacc(routine)        
    rename(routine)
    acc_seq(routine)
    scc_transform_routine(routine, lst_horizontal_size, lst_horizontal_idx, lst_horizontal_bounds, true_symbols, false_symbols, pathw)    
    alloc_temp(routine)

    
    Sourcefile.to_file(fgen(routine), Path(pathw+".F90"))
    if inline_match:
        print("inline_match = ", inline_match)
        with open(pathw+".F90", 'r', encoding='utf-8', errors='ignore') as file_caller:
            caller = file_caller.read()
            new_caller=caller.replace("CONTAINS", "")

        with open(pathw+".F90", 'w', encoding='utf-8', errors='ignore') as file_caller:
            file_caller.write(new_caller)
#            file_caller.write(caller)