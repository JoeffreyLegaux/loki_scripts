#!/usr/bin/bash

###~/loki/bin/activate
#source ~/loki_env_new/bin/activate
source ~/venv/bin/python
which python3

for f in \
   src2/actke.F90                       \
   src2/acbl89.F90                      \
   src2/acturb.F90                      \
   src2/acevolet.F90                    \
   src2/aclender.F90                    \
   src2/hl2fl.F90                       \
   src2/fl2hl.F90                       \
   src2/dprecips.F90                    \
   src2/ppwetpoint.F90                  \
   src2/acdrag.F90                      \
   src2/acpluis.F90                     \
   src2/acdnshf.F90                     \
   src2/radozcmf.F90                   \
   src2/suozon.F90                      \
   src2/cpphinp.F90                        \
   src2/aplpar_init.F90                 \
   src2/actqsat.F90                     \
   src2/acsol.F90                       \
   src2/acsolw.F90                      \
   src2/achmt.F90                       \
   src2/acntcls.F90                     \
   src2/achmtls.F90                     \
   src2/acclph.F90                      \
   src2/qngcor.F90                      \
   src2/acdrme.F90                      \
   src2/acevadcape.F90                  \
   src2/accldia.F90                     \
   src2/acvisih.F90                     \
   src2/mf_phys_precips.F90             \
   src2/mf_phys_bayrad.F90              \
   src2/acfluso.F90                     \
   src2/radaer.F90                     \
   src2/acnebcond.F90                   \
   src2/acnebsm.F90                     \
   src2/acnebn.F90                      \
   src2/acnebxrs.F90                    \
   src2/acnpart.F90                     \
   src2/acnpart_cloud_cover_wmo.F90     \
   src2/acnpart_cloud_cover.F90         \
   src2/acnpart_conv_base_top.F90       \
   src2/radheat.F90                    \
   src2/acdayd.F90                      \
   src2/acrso.F90                      \
   src2/acveg.F90                       \
   src2/acdifv1.F90                     \
   src2/arp_ground_param.F90            \
   src2/acdifv2.F90                     \
   src2/tridifv1.F90                    \
   src2/acuptq.F90                      \
   src2/acpluiz.F90                     \
   src2/acmicro.F90                     \
   src2/advprcs.F90                     \
   src2/acdrov.F90                      \
   src2/acdrov.F90                      \
   src2/cpwts.F90                          \
   src2/cptends.F90                        \
   src2/cptend_new.F90                     \
   src2/cpmvvps.F90                        
do 
python3 main.py $f
done

for f in \
   src2/cucalln_mf.F90                  \
   src2/cuadjtq.F90                      \
   src2/cuascn.F90                       \
   src2/cubasen.F90                      \
   src2/cuddrafn.F90                     \
   src2/cudlfsn.F90                      \
   src2/cududv.F90                       \
   src2/cuinin.F90                       \
   src2/cuadjtqs.F90                     \
   src2/cubasmcn.F90                     \
   src2/cubidiag.F90                     \
   src2/cuccdia.F90                      \
   src2/cuctracer.F90                    \
   src2/cudtdqn.F90                      \
   src2/cuentr.F90                       \
   src2/cuflxn.F90                       \
   src2/cumastrn.F90                     \
   src2/culight.F90                      \
   src2/satur.F90
do 
python3 main.py $f
done
