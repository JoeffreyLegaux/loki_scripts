# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.

my_dict_arpege = {
    "YDMODEL%YRML_DYN%YRDYN%NCURRENT_ITER": '0',
    "YDMODEL%YRML_DYN%YRDYN%RINTOPT" : '1._JPRB',
    "YDMODEL%YRML_DYN%YRDYN%NSPLTHOI" : '0',
    "YDMODEL%YRML_DYN%YRDYN%LADVF" : '.TRUE.',
    "YDMODEL%YRML_DYN%YRDYN%LRHS_CURV" : 'FALSE.',
    "YDMODEL%YRML_DYN%YRDYN%LSPLTHOIGFL" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LSLDIA" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LNHDYN" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LGWADV" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LPC_FULL" :  '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LPC_CHEAP" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LMIXETTLS" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LCOMAD" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LSLHD" : '.FALSE.',
    "YDMODEL%YRML_PHY_SLIN%YRPHLC%LSPHLC" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LTWOTL" : '.TRUE.',
    "YDMODEL%YRML_DYN%YRDYNA%LNHEE" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LNHQE" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LELTRA" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LSVTSM" : '.FALSE.',
    "YDGEOMETRY%YRCVER%LVERTFE" : '.TRUE.',
    "YDGEOMETRY%YRVERT_GEOM%YRCVER%LVERTFE" : '.TRUE.',
    "YDMODEL%YRML_DYN%YRDYNA%LSPRT" : '.TRUE.',
    "YDMODEL%YRML_DYN%YRDYNA%NVDVAR" : '3',
    "YDMODEL%YRML_DYN%YRDYNA%ND4SYS" : '2',
    "YDMODEL%YRML_DYN%YRDYNA%LSLINL" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LSACC" : '.FALSE.',
    "LLCT" : '.FALSE.',
    "LLCTC" : '.FALSE.',
    "YDMODEL%YRML_PHY_EC%YREPHY%LSLPHY" : '.FALSE.',
    "YDCPG_OPTS%LSFORC" : '.FALSE.',
    "LGWDIAGS_ON" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRSPNG%LNSPONGE" : '.FALSE.',
    "YDCPG_OPTS%LNUDG" : '.FALSE.',
    "YDMODEL%YRML_LBC%LTENC" : '.FALSE.',
    "YDMODEL%YRML_DIAG%YRLDDH%LRSLDDH" : '.FALSE.',
    "YDMODEL%YRML_DIAG%YRLDDH%LSDDH" : '.FALSE.',
    "YDCPG_OPTS%LSLPHY" : '.FALSE.',
    "YDMODEL%YRML_PHY_EC%YRECUCONVCA%LCUCONV_CA" : '.FALSE.',
    "YDMODEL%YRML_PHY_MF%YRSIMPHL%LTRAJPS" : '.FALSE.',
    "LTRAJSAVE" : '.FALSE.',
    "YDMODEL%YRML_DYN%YRDYNA%LSLAG" : '.TRUE.',
    "LLVERINT_ON_CPU" : '.FALSE.',
    "LLSIMPLE_DGEMM" : '.TRUE.',
}

my_dict = {
    "LMUSCLFA" : '.FALSE.',
    "YDLDDH%LFLEXDIA" : '.FALSE.',
    "YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA" : '.FALSE.',
    "YDSPP_CONFIG%LSPP" : '.FALSE.',
    "LMCAPEA" : '.FALSE.',
    "OCH1CONV" : '.FALSE.',
}


def symbols():
    false_symbols=[]
    true_symbols=[]
    for var in my_dict:
        if my_dict[var]==".FALSE." :
            false_symbols.append(var)
    
        if my_dict[var]==".TRUE." :
            true_symbols.append(var)

    for var in my_dict_arpege:
        if my_dict_arpege[var]==".FALSE." :
            false_symbols.append(var)
    
        if my_dict_arpege[var]==".TRUE." :
            true_symbols.append(var)

        #print("false_symbols=",false_symbols)
        #print("true_symbols=",true_symbols)
    return (true_symbols, false_symbols)
