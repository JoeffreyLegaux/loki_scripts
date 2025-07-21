SUBROUTINE CUMASTRN_OPENACC (PPLDARE, PPLRG, YDTHF, YDCST, YDML_PHY_SLIN, YDML_PHY_EC, YGFL, YDCHEM, YDSPP_CONFIG, KIDIA, KFDIA,  &
& KLON, KLEV, PDX, LDTDKMF, LDMCAPEA, LDLAND, PTSPHY, PTEN, PQEN, PUEN, PVEN, PLITOT, PVERVEL, PQSEN, PQHFL, PAHFS, PAP, PAPH,  &
& PGEO, PGEOH, PGAW, PCUCONVCA, PGP2DSPP, PTENT, PTENQ, PTENU, PTENV, PTENTA, PTENQA, LDCUM, KTYPE, KCBOT, KCTOP, LDCUM_LIG,  &
& KCBOT_LIG, KCTOP_LIG, KBOTSC, LDSC, LDSHCV, PLCRIT_AER, PTU, PQU, PLU, PLUDE, PLUDELI, PSNDE, PENTH, PMFLXR, PMFLXS, PRAIN,  &
& PLRAIN, PRSUD, PMFU, PMFD, PMFUDE_RATE, PMFDDE_RATE, PCAPE, PWMEAN, PVDISCU, PDISS, KTRAC, PCEN, PTENC, PSCAV, PSCAV0, YDSTACK)
  
  !**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
  
  !     PURPOSE
  !     -------
  !          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
  !     PROGNOSTIC VARIABLES T,Q,U, V AND TRACERS DUE TO CONVECTIVE PROCESSES.
  !     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
  !     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
  !     SATURATED CUMULUS DOWNDRAFTS.
  
  !**   INTERFACE.
  !     ----------
  
  !          *CUMASTR* IS CALLED FROM *CUCALL*
  !     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
  !     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
  !     IT RETURNS ITS OUTPUT TO THE SAME SPACE
  !      1.MODIFIED TENDENCIES OF MODEL VARIABLES
  !      2.RATES OF CONVECTIVE PRECIPITATION
  !        (USED IN SUBROUTINE SURF)
  !      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
  !        (USED IN SUBROUTINE CLOUD)
  
  !     METHOD
  !     -------
  
  !     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
  !        (1) DEFINE CONSTANTS AND PARAMETERS
  !        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
  !            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
  !        (3) CALCULATE CLOUD BASE IN 'CUBASE'
  !            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
  !        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
  !        (5) DO DOWNDRAFT CALCULATIONS:
  !              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
  !              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
  !              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
  !                  EFFECT OF CU-DOWNDRAFTS
  !        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
  !        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
  !            DO EVAPORATION IN SUBCLOUD LAYER
  !        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
  !        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KTRAC*        NUMBER OF CHEMICAL TRACERS
  
  !     INPUT PARAMETERS (LOGICAL)
  
  !    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
  !    *LDTDKMF*      Arpege tuning (if TRUE)
  
  !     INPUT PARAMETERS (REAL)
  
  !    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
  !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
  !    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
  !    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
  !    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
  !    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATIONS KG/KG
  !    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT            KG/KG
  !    *PVERVEL*      VERTICAL VELOCITY                             PA/S
  !    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
  !    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
  !    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
  !    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
  !    *PGEO*         GEOPOTENTIAL                                  M2/S2
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
  !    *PGAW*       NORMALISED GAUSSIAN QUADRATURE WEIGHT / NUMBER OF LONGITUDE POINTS
  !                           LOCAL SUB-AREA == 4*RPI*RA**2 * PGAW
  !    *PLCRIT_AER*   CRITICAL LIQUID MMR FOR AUTOCONVERSION PROCESS KG/KG
  !    *PCUCONVCA*    COUPLING TO CELL AUTOMATON OR 2D ADV FIELD     ( )
  !    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)
  !    *PSCAV*        SCAVENGING COEFFICIENTS FOR CHEM TRANSPORT    UNITLESS
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PTENT*        TEMPERATURE TENDENCY                           K/S
  !    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
  !    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
  !    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
  !    *PTENC*        TENDENCY OF CHEMICAL TRACERS                   1/S
  !    *PTENTA*       TEMPERATURE TENDENCY DYNAMICS=TOT ADVECTION    K/S
  !    *PTENQA*       MOISTURE    TENDENCY DYNAMICS=TOT ADVECTION    1/S
  
  !    OUTPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  !    *LDCUM_LIG*    FLAG: .TRUE. FOR CONVECTIVE POINTS (FOR LIGHTNING PARAM)
  !    *LDSC*         FLAG: .TRUE. FOR SC-POINTS
  
  !    OUTPUT PARAMETERS (INTEGER):
  
  !    *KTYPE*        TYPE OF CONVECTION
  !                       1 = PENETRATIVE CONVECTION
  !                       2 = SHALLOW CONVECTION
  !                       3 = MIDLEVEL CONVECTION
  !    *KCBOT*        CLOUD BASE LEVEL
  !    *KCTOP*        CLOUD TOP LEVEL
  !    *KCBOT_LIG*    CLOUD BASE LEVEL (FOR LIGHTNING PARAM)
  !    *KCTOP_LIG*    CLOUD TOP LEVEL (FOR LIGHTNING PARAM)
  !    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PTU*          TEMPERATURE IN UPDRAFTS                         K
  !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
  !    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
  !    *PLUDE*        DETRAINED TOTAL CONDENSATE                    KG/(M2*S)
  !    *PLUDELI*      DETRAINED LIQUID, ICE, VAPOR, T               KG/(M2*S)
  !    *PSNDE*        DETRAINED SNOW/RAIN                           KG/(M2*S)
  !    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
  !    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
  !    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
  !    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
  !                   (NO EVAPORATION IN DOWNDRAFTS)
  !    *PLRAIN*       RAIN+SNOW CONTENT IN UPDRAFTS                 KG/KG
  !    *PRSUD*        CONVECT MEAN RAIN AND SNOW CONTENT IN COLUMN  KG/KG
  !    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
  !    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
  !    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                      KG/(M3*S)
  !    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M3*S)
  !    *PCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG
  !    *PWMEAN*       VERTICALLY AVERAGED UPDRAUGHT VELOCITY         M/S
  !    *PVDISCU*      VERT INTEGRATED CONVECTIVE DISSIPATION RATE    W/m2
  !    *PDISS*        DISSIPATION BY DEEP CONVECTION ESTIMATE       m2/s3^(1/3)
  
  !     EXTERNALS.
  !     ----------
  
  !       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
  !       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
  !       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
  !       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
  !       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
  !       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
  !       CUDQDT: UPDATES TENDENCIES FOR T AND Q
  !       CUDUDV: UPDATES TENDENCIES FOR U AND V
  !       CUCTRACER: TRACER TRANSPORT
  
  !     SWITCHES.
  !     --------
  
  !          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
  !          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
  !          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
  !          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
  !          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
  !          LMFTRAC=.false. TRACER TRANSPORT
  
  !     REFERENCE.
  !     ----------
  
  !          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
  !          DRAFT PAPER ON MASSFLUX SCHEME (NORDENG, 1995)
  !          Bechtold et al. (2008 QJRMS 134,1337-1351), Rooy et al. (2012 QJRMS)
  !          Bechtold et al. (2013 JAS)
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
  
  !     MODIFICATIONS.
  !     --------------
  !      03-08-29 : Clean-up, deep/shallow switches  P.Bechtold
  !      04-02-11 : Add tracer transport             P.Bechtold
  !      05-02-11 : Positive scaling of total Mflux  P.Bechtold
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      04-12-03 : Turn off shallow convection over stratocu. M.Ko"hler
  !      05-06-27 : Switch off ddraught if idtop<kctop
  !                 correction for detrainment rates P.Bechtold
  !      05-11-22 : Mods for coarser/finer physics D.Salmond + M.Hortal
  !      06-02-11 : Enable TQ implicit               P.Bechtold
  !      07-06-01 : Only single updraught call with  P.Bechtold
  !                 scaling, convective turnover time
  !                 scale, contain momentum computations in cumastrn
  !      07-10-09 : Added KE dissipation and convective scavenging   P. Bechtold
  !      12-03-02 : remove all entrainment stuff     P. Bechtold
  !      04-10-12 : Add RPLRG/RPLDARE for small planet  N.Semane+P.Bechtold
  !      13-02-23 : modif diurnal cycle CAPE closure P. Bechtold
  !      13-10-01 : modified option RCAPDCYCL=1 for diurnal cycle N. Semane
  !      16-01-27 : Introduced SPP scheme (LSPP)     M. Leutbecher & S.-J. Lock
  !      20180303 : Gabor: Just a comment line to force recompilation due to
  !                        compiler wrapper optimation exception liat change
  !      20-08-25 : Additional total moisture advection in CAPE closure P. Bechtold&T. Becker
  !      08-12-20 : Added separate fields for lightning parameterization  P. Lopez.
  !      20-10-12 : SPP abstraction and revision     M. Leutbecher & S. Lang
  !      20210913 : Modifications for Arpege Y.Bouteloup (LDTDKMF)
  !    2020-09-18 : Choose CAPE diagnostic (LMCAPEA) J.M. Piriou
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !      15-11-21 : Added PVDISCU, vert-integrated convective dissipation rate R.Forbes
  !------------------------------------------------------------------------------------------
  
!$acc routine( CUMASTRN_OPENACC ) seq
  
  USE MODEL_PHYSICS_ECMWF_MOD, ONLY: MODEL_PHYSICS_ECMWF_TYPE
  USE MODEL_PHYSICS_SIMPLINEAR_MOD, ONLY: MODEL_PHYSICS_SIMPLINEAR_TYPE
  USE YOM_YGFL, ONLY: TYPE_GFLD
  USE YOMCHEM, ONLY: TCHEM
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE SPP_MOD, ONLY: TSPP_CONFIG
  USE SPP_GEN_MOD, ONLY: SPP_PERT
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  REAL(KIND=JPRB), INTENT(IN) :: PPLDARE
  REAL(KIND=JPRB), INTENT(IN) :: PPLRG
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(MODEL_PHYSICS_ECMWF_TYPE), INTENT(IN) :: YDML_PHY_EC
  TYPE(TCHEM), INTENT(IN) :: YDCHEM
  TYPE(TSPP_CONFIG), INTENT(IN) :: YDSPP_CONFIG
  TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE), INTENT(IN) :: YDML_PHY_SLIN
  TYPE(TYPE_GFLD), INTENT(IN) :: YGFL
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KTRAC
  LOGICAL, INTENT(IN) :: LDMCAPEA
  LOGICAL, INTENT(IN) :: LDLAND(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
  REAL(KIND=JPRB), INTENT(IN) :: PLCRIT_AER(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLITOT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVERVEL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQHFL(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAHFS(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PGAW(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PCUCONVCA(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
  REAL(KIND=JPRB), INTENT(IN) :: PCEN(KLON, KLEV, KTRAC)
  REAL(KIND=JPRB), INTENT(IN) :: PSCAV(KTRAC)
  REAL(KIND=JPRB), INTENT(IN) :: PSCAV0(KTRAC)
  REAL(KIND=JPRB), INTENT(IN) :: PDX(KLON)
  LOGICAL, INTENT(IN) :: LDTDKMF
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENV(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTENTA(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTENQA(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENC(KLON, KLEV, KTRAC)
  LOGICAL, INTENT(OUT) :: LDCUM(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KTYPE(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCTOP(KLON)
  LOGICAL, INTENT(OUT) :: LDCUM_LIG(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCBOT_LIG(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCTOP_LIG(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KBOTSC(KLON)
  LOGICAL, INTENT(OUT) :: LDSC(KLON)
  LOGICAL, INTENT(IN) :: LDSHCV(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PTU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLUDE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLUDELI(KLON, KLEV, 4)
  REAL(KIND=JPRB), INTENT(OUT) :: PSNDE(KLON, KLEV, 2)
  REAL(KIND=JPRB), INTENT(OUT) :: PENTH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFLXR(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFLXS(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(OUT) :: PRAIN(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PLRAIN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PRSUD(KLON, KLEV, 2)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUDE_RATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFDDE_RATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PCAPE(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PWMEAN(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PVDISCU(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PDISS(KLON, KLEV)
  
  
  temp (REAL (KIND=JPRB), ZTENH, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZQENH, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZTENH2, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZQENH2, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZQSENH, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZTD, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZQD, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUQ, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDQ, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDMFUP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDMFDP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUL, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZRFL, (KLON))
  temp (REAL (KIND=JPRB), ZUU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZVU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZUD, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZVD, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZKINEU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZKINED, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUB, (KLON))
  REAL(KIND=JPRB) :: ZMFUB1
  REAL(KIND=JPRB) :: ZKHVFL
  REAL(KIND=JPRB) :: ZKHFL
  REAL(KIND=JPRB) :: ZDQCV
  REAL(KIND=JPRB) :: ZFACCA
  REAL(KIND=JPRB) :: ZMFCFL
  REAL(KIND=JPRB) :: ZSATFR
  temp (REAL (KIND=JPRB), ZDPMEL, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZLGLAC, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZDHPBL
  temp (REAL (KIND=JPRB), ZWUBASE, (KLON))
  
  temp (REAL (KIND=JPRB), ZDMFEN, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDMFDE, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZDX
  
  temp (INTEGER (KIND=JPIM), ILAB, (KLON, KLEV))
  temp (INTEGER (KIND=JPIM), IDTOP, (KLON))
  temp (INTEGER (KIND=JPIM), ICTOP0, (KLON))
  temp (INTEGER (KIND=JPIM), IDPL, (KLON))
  REAL(KIND=JPRB) :: ZCAPE
  REAL(KIND=JPRB) :: ZCAPE2
  REAL(KIND=JPRB) :: ZHEAT
  REAL(KIND=JPRB) :: ZCAPPBL
  REAL(KIND=JPRB) :: ZCAPDCYCL
  temp (LOGICAL, LLDDRAF, (KLON))
  temp (LOGICAL, LLDDRAF3, (KLON))
  temp (LOGICAL, LLDCUM, (KLON))
  temp (LOGICAL, LLSCVFLAG, (KLON))
  LOGICAL :: LLO1
  LOGICAL :: LLO2
  LOGICAL :: LLMIXS
  LOGICAL :: LLO3
  
  INTEGER(KIND=JPIM) :: IKB
  INTEGER(KIND=JPIM) :: IKD
  INTEGER(KIND=JPIM) :: ITOPM2
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZCONS2  !, ZALPHA
  REAL(KIND=JPRB) :: ZCONS  !, ZALPHA
  REAL(KIND=JPRB) :: ZDH  !, ZALPHA
  REAL(KIND=JPRB) :: ZDH2  !, ZALPHA
  REAL(KIND=JPRB) :: ZDQMIN  !, ZALPHA
  REAL(KIND=JPRB) :: ZDZ  !, ZALPHA
  REAL(KIND=JPRB) :: ZEPS  !, ZALPHA
  REAL(KIND=JPRB) :: ZFAC  !, ZALPHA
  REAL(KIND=JPRB) :: ZMFMAX  !, ZALPHA
  REAL(KIND=JPRB) :: ZPBMPT  !, ZALPHA
  REAL(KIND=JPRB) :: ZQUMQE  !, ZALPHA
  REAL(KIND=JPRB) :: ZRO  !, ZALPHA
  REAL(KIND=JPRB) :: ZMFA  !, ZALPHA
  REAL(KIND=JPRB) :: ZERATE  !, ZALPHA
  REAL(KIND=JPRB) :: ZDERATE  !, ZALPHA
  REAL(KIND=JPRB) :: ZORCPD  !, ZALPHA
  REAL(KIND=JPRB) :: ZRDOCPD  !, ZALPHA
  REAL(KIND=JPRB) :: ZDUTEN  !, ZALPHA
  REAL(KIND=JPRB) :: ZDVTEN  !, ZALPHA
  REAL(KIND=JPRB) :: ZTDIS  !, ZALPHA
  REAL(KIND=JPRB) :: ZTAURES  !, ZALPHA
  REAL(KIND=JPRB) :: ZAR  !, ZALPHA
  REAL(KIND=JPRB) :: ZAS  !, ZALPHA
  REAL(KIND=JPRB) :: ZBR  !, ZALPHA
  REAL(KIND=JPRB) :: ZBS  !, ZALPHA
  REAL(KIND=JPRB) :: ZAK  !, ZALPHA
  REAL(KIND=JPRB) :: ZRG  !, ZALPHA
  
  REAL(KIND=JPRB) :: ZTAU  ! adjustment time
  REAL(KIND=JPRB) :: ZTAUPBL  ! adjustment time
  REAL(KIND=JPRB) :: ZXTAU
  
  ! A bunch of SPP variables
  LOGICAL :: LLPERT_RTAU  ! SPP perturbation on?
  INTEGER(KIND=JPIM) :: IPRTAU  ! SPP random field pointer
  INTEGER(KIND=JPIM) :: IPN  !  SPP perturbation pointer
  TYPE(SPP_PERT) :: PN1  ! SPP pertn. configs. for RTAU
  
  ! scaling factor for momentum and tracer massflux
  REAL(KIND=JPRB) :: ZMFS
  temp (REAL (KIND=JPRB), ZMFUUS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDUS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUDR, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDDR, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZTENU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZTENV, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZMFUUB
  REAL(KIND=JPRB) :: ZMFUVB
  REAL(KIND=JPRB) :: ZMF_SHAL
  temp (REAL (KIND=JPRB), ZKMFL, (KLON))
  temp (REAL (KIND=JPRB), ZDMFUPC, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDMFDPC, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZUV2, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZSUM12
  REAL(KIND=JPRB) :: ZSUM22
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cuascn_openacc.intfb.h"
#include "cubasen_openacc.intfb.h"
#include "cuddrafn_openacc.intfb.h"
#include "cudlfsn_openacc.intfb.h"
#include "cudtdqn_openacc.intfb.h"
#include "cududv_openacc.intfb.h"
#include "cuflxn_openacc.intfb.h"
#include "cuinin_openacc.intfb.h"
#include "cuctracer_openacc.intfb.h"
  
#include "fcttre.func.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZTENH) == 8) THEN
    alloc8 (ZTENH)
  ELSE
    IF (KIND (ZTENH) == 4) THEN
      alloc4 (ZTENH)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZQENH) == 8) THEN
    alloc8 (ZQENH)
  ELSE
    IF (KIND (ZQENH) == 4) THEN
      alloc4 (ZQENH)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZTENH2) == 8) THEN
    alloc8 (ZTENH2)
  ELSE
    IF (KIND (ZTENH2) == 4) THEN
      alloc4 (ZTENH2)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZQENH2) == 8) THEN
    alloc8 (ZQENH2)
  ELSE
    IF (KIND (ZQENH2) == 4) THEN
      alloc4 (ZQENH2)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZQSENH) == 8) THEN
    alloc8 (ZQSENH)
  ELSE
    IF (KIND (ZQSENH) == 4) THEN
      alloc4 (ZQSENH)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZTD) == 8) THEN
    alloc8 (ZTD)
  ELSE
    IF (KIND (ZTD) == 4) THEN
      alloc4 (ZTD)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZQD) == 8) THEN
    alloc8 (ZQD)
  ELSE
    IF (KIND (ZQD) == 4) THEN
      alloc4 (ZQD)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUS) == 8) THEN
    alloc8 (ZMFUS)
  ELSE
    IF (KIND (ZMFUS) == 4) THEN
      alloc4 (ZMFUS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFDS) == 8) THEN
    alloc8 (ZMFDS)
  ELSE
    IF (KIND (ZMFDS) == 4) THEN
      alloc4 (ZMFDS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUQ) == 8) THEN
    alloc8 (ZMFUQ)
  ELSE
    IF (KIND (ZMFUQ) == 4) THEN
      alloc4 (ZMFUQ)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFDQ) == 8) THEN
    alloc8 (ZMFDQ)
  ELSE
    IF (KIND (ZMFDQ) == 4) THEN
      alloc4 (ZMFDQ)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDMFUP) == 8) THEN
    alloc8 (ZDMFUP)
  ELSE
    IF (KIND (ZDMFUP) == 4) THEN
      alloc4 (ZDMFUP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDMFDP) == 8) THEN
    alloc8 (ZDMFDP)
  ELSE
    IF (KIND (ZDMFDP) == 4) THEN
      alloc4 (ZDMFDP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUL) == 8) THEN
    alloc8 (ZMFUL)
  ELSE
    IF (KIND (ZMFUL) == 4) THEN
      alloc4 (ZMFUL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZRFL) == 8) THEN
    alloc8 (ZRFL)
  ELSE
    IF (KIND (ZRFL) == 4) THEN
      alloc4 (ZRFL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZUU) == 8) THEN
    alloc8 (ZUU)
  ELSE
    IF (KIND (ZUU) == 4) THEN
      alloc4 (ZUU)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZVU) == 8) THEN
    alloc8 (ZVU)
  ELSE
    IF (KIND (ZVU) == 4) THEN
      alloc4 (ZVU)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZUD) == 8) THEN
    alloc8 (ZUD)
  ELSE
    IF (KIND (ZUD) == 4) THEN
      alloc4 (ZUD)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZVD) == 8) THEN
    alloc8 (ZVD)
  ELSE
    IF (KIND (ZVD) == 4) THEN
      alloc4 (ZVD)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZKINEU) == 8) THEN
    alloc8 (ZKINEU)
  ELSE
    IF (KIND (ZKINEU) == 4) THEN
      alloc4 (ZKINEU)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZKINED) == 8) THEN
    alloc8 (ZKINED)
  ELSE
    IF (KIND (ZKINED) == 4) THEN
      alloc4 (ZKINED)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUB) == 8) THEN
    alloc8 (ZMFUB)
  ELSE
    IF (KIND (ZMFUB) == 4) THEN
      alloc4 (ZMFUB)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDPMEL) == 8) THEN
    alloc8 (ZDPMEL)
  ELSE
    IF (KIND (ZDPMEL) == 4) THEN
      alloc4 (ZDPMEL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZLGLAC) == 8) THEN
    alloc8 (ZLGLAC)
  ELSE
    IF (KIND (ZLGLAC) == 4) THEN
      alloc4 (ZLGLAC)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZWUBASE) == 8) THEN
    alloc8 (ZWUBASE)
  ELSE
    IF (KIND (ZWUBASE) == 4) THEN
      alloc4 (ZWUBASE)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDMFEN) == 8) THEN
    alloc8 (ZDMFEN)
  ELSE
    IF (KIND (ZDMFEN) == 4) THEN
      alloc4 (ZDMFEN)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDMFDE) == 8) THEN
    alloc8 (ZDMFDE)
  ELSE
    IF (KIND (ZDMFDE) == 4) THEN
      alloc4 (ZDMFDE)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ILAB) == 8) THEN
    alloc8 (ILAB)
  ELSE
    IF (KIND (ILAB) == 4) THEN
      alloc4 (ILAB)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (IDTOP) == 8) THEN
    alloc8 (IDTOP)
  ELSE
    IF (KIND (IDTOP) == 4) THEN
      alloc4 (IDTOP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ICTOP0) == 8) THEN
    alloc8 (ICTOP0)
  ELSE
    IF (KIND (ICTOP0) == 4) THEN
      alloc4 (ICTOP0)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (IDPL) == 8) THEN
    alloc8 (IDPL)
  ELSE
    IF (KIND (IDPL) == 4) THEN
      alloc4 (IDPL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (LLDDRAF) == 8) THEN
    alloc8 (LLDDRAF)
  ELSE
    IF (KIND (LLDDRAF) == 4) THEN
      alloc4 (LLDDRAF)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (LLDDRAF3) == 8) THEN
    alloc8 (LLDDRAF3)
  ELSE
    IF (KIND (LLDDRAF3) == 4) THEN
      alloc4 (LLDDRAF3)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (LLDCUM) == 8) THEN
    alloc8 (LLDCUM)
  ELSE
    IF (KIND (LLDCUM) == 4) THEN
      alloc4 (LLDCUM)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (LLSCVFLAG) == 8) THEN
    alloc8 (LLSCVFLAG)
  ELSE
    IF (KIND (LLSCVFLAG) == 4) THEN
      alloc4 (LLSCVFLAG)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUUS) == 8) THEN
    alloc8 (ZMFUUS)
  ELSE
    IF (KIND (ZMFUUS) == 4) THEN
      alloc4 (ZMFUUS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFDUS) == 8) THEN
    alloc8 (ZMFDUS)
  ELSE
    IF (KIND (ZMFDUS) == 4) THEN
      alloc4 (ZMFDUS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUDR) == 8) THEN
    alloc8 (ZMFUDR)
  ELSE
    IF (KIND (ZMFUDR) == 4) THEN
      alloc4 (ZMFUDR)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFDDR) == 8) THEN
    alloc8 (ZMFDDR)
  ELSE
    IF (KIND (ZMFDDR) == 4) THEN
      alloc4 (ZMFDDR)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZTENU) == 8) THEN
    alloc8 (ZTENU)
  ELSE
    IF (KIND (ZTENU) == 4) THEN
      alloc4 (ZTENU)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZTENV) == 8) THEN
    alloc8 (ZTENV)
  ELSE
    IF (KIND (ZTENV) == 4) THEN
      alloc4 (ZTENV)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZKMFL) == 8) THEN
    alloc8 (ZKMFL)
  ELSE
    IF (KIND (ZKMFL) == 4) THEN
      alloc4 (ZKMFL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDMFUPC) == 8) THEN
    alloc8 (ZDMFUPC)
  ELSE
    IF (KIND (ZDMFUPC) == 4) THEN
      alloc4 (ZDMFUPC)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDMFDPC) == 8) THEN
    alloc8 (ZDMFDPC)
  ELSE
    IF (KIND (ZDMFDPC) == 4) THEN
      alloc4 (ZDMFDPC)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZUV2) == 8) THEN
    alloc8 (ZUV2)
  ELSE
    IF (KIND (ZUV2) == 4) THEN
      alloc4 (ZUV2)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KIDIA
  
  !---------------------------------------------------------------------
  
  !     1.           SPECIFY CONSTANTS AND PARAMETERS
  !                  --------------------------------
  
  
  ZCONS2 = YDML_PHY_EC%YRECUMF%RMFCFL / (YDCST%RG*PTSPHY)
  ZCONS = 1.0_JPRB / (YDCST%RG*PTSPHY)
  ZORCPD = 1.0_JPRB / YDCST%RCPD
  ZRDOCPD = YDCST%RD*ZORCPD
  ZRG = 1.0_JPRB / YDCST%RG
  ZTAU = 0.0
  
  ! prepare SPP PERTURBATIONS
  LLPERT_RTAU = .false.
  
  !----------------------------------------------------------------------
  PCAPE(JLON) = 0.0_JPRB
  PVDISCU(JLON) = 0.0_JPRB
  LDCUM(JLON) = .false.
  
  !*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
  !                  ---------------------------------------------------
  
  CALL CUININ_OPENACC(YDCST, YDTHF, YDML_PHY_SLIN%YREPHLI, YDML_PHY_EC%YRECUMF, KIDIA, KFDIA, KLON, KLEV, PTEN, PQEN, PQSEN,  &
  & PUEN, PVEN, PGEO, PAPH, ILAB, ZTENH, ZQENH, ZQSENH, PGEOH, PTU, PQU, ZTD, ZQD, ZUU, ZVU, ZUD, ZVD, PLU, YDSTACK=YLSTACK)
  
  !---------------------------------------------------------------------
  
  !*    3.0          CLOUD BASE CALCULATIONS
  !                  -----------------------
  
  !*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
  !                  ---------------------------------------
  
  ZKMFL(JLON) = 0.0_JPRB
  LLMIXS = LDTDKMF
  CALL CUBASEN_OPENACC(PPLDARE, PPLRG, YDTHF, YDCST, YDML_PHY_SLIN%YREPHLI, YDML_PHY_EC%YRECLDP, YDML_PHY_EC%YRECUMF,  &
  & YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, YDML_PHY_EC%YRECUMF%NJKT1, LLMIXS, LDTDKMF, ZTENH, ZQENH, PGEOH, PAPH, PQHFL, PAHFS,  &
  & PGP2DSPP, ZKMFL, PTEN, PQEN, PQSEN, PGEO, PTU, PQU, PLU, ZKINEU, ZWUBASE, ILAB, LDCUM, LDSC, KCBOT, KBOTSC, ICTOP0, IDPL,  &
  & PCAPE, YDSTACK=YLSTACK)
  
  
  !*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
  !*                 DECIDE ON TYPE OF CUMULUS CONVECTION
  !*                 ONE THE BASIS OF THE DEPTH OF THE CONVECTION
  !*                 DEEP IF CLOUD DEPTH > 200MB
  !*                 SHALLOW IF CLOUD DEPTH <200MB
  !                  -----------------------------------------
  
  ! CALCULATE COLUMN AND SUB CLOUD LAYER MOISTURE CONVERGENCE
  ! AND SUB CLOUD LAYER MOIST STATIC ENERGY CONVERGENCE
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  ZDHPBL = 0.0_JPRB
  IDTOP(JLON) = 0
  ZCAPPBL = 0.
  ZKHVFL = (-PAHFS(JLON, KLEV + 1)*ZORCPD - YDCST%RETV*PTEN(JLON, KLEV)*PQHFL(JLON, KLEV + 1)) / (PPLRG*PPLDARE)
  ZKHFL = (-PAHFS(JLON, KLEV + 1) - YDCST%RLVTT*PQHFL(JLON, KLEV + 1)) / (PPLRG*PPLDARE)
  DO JK=YDML_PHY_EC%YRECUMF%NJKT2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LDCUM(JLON) .and. JK >= KCBOT(JLON)) THEN
      ZDZ = (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
      ZDHPBL = ZDHPBL + (YDCST%RLVTT*PTENQ(JLON, JK) + YDCST%RCPD*PTENT(JLON, JK))*ZDZ
      ZCAPPBL = ZCAPPBL + (PTENT(JLON, JK) + YDCST%RETV*PTEN(JLON, JK)*PTENQ(JLON, JK))*ZDZ
    END IF
  END DO
  
  !*                 ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
  !*                 CALCULATIONS IN CUASC AND INITIAL DETERMINATION OF
  !*                 CLOUD TYPE
  !*                 (MAX.POSSIBLE CLOUD HEIGHT
  !*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
  !                  -------------------------------------------------
  
  
  !*                 SPECIFY INITIAL CLOUD TYPE
  !*
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LDCUM(JLON)) THEN
    IKB = KCBOT(JLON)
    ITOPM2 = ICTOP0(JLON)
    ZPBMPT = PAPH(JLON, IKB) - PAPH(JLON, ITOPM2)
    IF (ZPBMPT >= YDML_PHY_EC%YRECUMF%RDEPTHS) THEN
      KTYPE(JLON) = 1
    ELSE
      KTYPE(JLON) = 2
    END IF
  ELSE
    KTYPE(JLON) = 0
  END IF
  
  !*             (C) calculate initial updraught mass flux
  !*                 and set lateral mixing rates
  !*
  !*                 for deep convection assume it is 10% of
  !*                 maximum value which is determined by the
  !*                 thickness of the layer and timestep
  !*
  !*                 for shallow convection calculated assuming
  !*                 a balance of moist static energy in the
  !*                 sub-cloud layer (ignores present of downdraughts)
  !                  ------------------------------------------
  
  IF (YDML_PHY_EC%YRECUMF%LMFWSTAR) THEN
    !DIR$ IVDEP
    !OCL NOVREC
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LDCUM(JLON)) THEN
      IKB = KCBOT(JLON)
      ZDZ = MAX(0.0_JPRB, MIN(1.5E3_JPRB, (PGEOH(JLON, IKB) - PGEOH(JLON, KLEV + 1)) / YDCST%RG))
      ZMF_SHAL = 0.07_JPRB*(YDCST%RG / PTEN(JLON, KLEV)*ZDZ*MAX(0.0_JPRB, ZKHVFL))**.3333
      ZMFMAX = (PAPH(JLON, IKB) - PAPH(JLON, IKB - 1))*ZCONS2
      ZMF_SHAL = MIN(ZMF_SHAL, ZMFMAX)
    END IF
  END IF
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LDCUM(JLON)) THEN
    IKB = KCBOT(JLON)
    ZMFMAX = (PAPH(JLON, IKB) - PAPH(JLON, IKB - 1))*ZCONS2
    
    ! deep convection
    
    IF (KTYPE(JLON) == 1) THEN
      ZMFUB(JLON) = ZMFMAX*0.1_JPRB
      
    ELSE IF (KTYPE(JLON) == 2) THEN
      
      ! shallow convection
      
      ZQUMQE = PQU(JLON, IKB) + PLU(JLON, IKB) - ZQENH(JLON, IKB)
      ZDQMIN = MAX(0.01_JPRB*ZQENH(JLON, IKB), 1.E-10_JPRB)
      ZDH = YDCST%RCPD*(PTU(JLON, IKB) - ZTENH(JLON, IKB)) + YDCST%RLVTT*ZQUMQE
      ZDH = YDCST%RG*MAX(ZDH, 1.E5_JPRB*ZDQMIN)
      IF (ZDHPBL > 0.0_JPRB) THEN
        ZMFUB(JLON) = ZDHPBL / ZDH
        ZMFUB(JLON) = MIN(ZMFUB(JLON), ZMFMAX)
      ELSE
        ZMFUB(JLON) = ZMFMAX*0.1_JPRB
        LDCUM(JLON) = .false.
      END IF
      IF (YDML_PHY_EC%YRECUMF%LMFWSTAR) ZMFUB(JLON) = ZMF_SHAL
    END IF
    
  ELSE
    
    ! no buoyancy cloud base from surface
    ! set cloud base mass flux and mixing rate
    ! to default value for safety
    
    ZMFUB(JLON) = 0.0_JPRB
  END IF
  
  !-----------------------------------------------------------------------
  
  !*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
  !                  -------------------------------------------
  
  !*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
  !*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
  !*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
  !                  -------------------------------------------------
  
  ! CALCULATIONS NOW DONE IS SECTION 3 ABOVE SO THAT
  ! INITIAL CLOUD DEPTH CAN BE USED TO SPECIFY
  ! THE TYPE OF CONVECTION
  
  !*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
  !                  --------------------------------------------
  
  CALL CUASCN_OPENACC(YDTHF, YDCST, YDML_PHY_SLIN%YREPHLI, YDML_PHY_EC%YRECLDP, YDML_PHY_EC%YRECUMF, YDSPP_CONFIG, YGFL, KIDIA,  &
  & KFDIA, KLON, KLEV, LDTDKMF, PTSPHY, ZTENH, ZQENH, PUEN, PVEN, PTEN, PQEN, PQSEN, PLITOT, PGEO, PGEOH, PAP, PAPH, PVERVEL,  &
  & ZWUBASE, PGP2DSPP, LDLAND, LDCUM, KTYPE, ILAB, LLSCVFLAG, PTU, PQU, PLU, PLRAIN, PMFU, ZMFUB, ZLGLAC, ZMFUS, ZMFUQ, ZMFUL,  &
  & PLUDE, PLUDELI, ZDMFUP, PLCRIT_AER, ZDMFEN, KCBOT, KCTOP, ICTOP0, IDPL, PMFUDE_RATE, ZKINEU, PWMEAN, YDSTACK=YLSTACK)
  
  !*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
  !              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
  !              -----------------------------------------------------
  
  !DIR$ IVDEP
  !OCL NOVREC
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LDCUM(JLON)) THEN
    IKB = KCBOT(JLON)
    ITOPM2 = KCTOP(JLON)
    ZPBMPT = PAPH(JLON, IKB) - PAPH(JLON, ITOPM2)
    IF (KTYPE(JLON) == 1 .and. ZPBMPT < YDML_PHY_EC%YRECUMF%RDEPTHS) KTYPE(JLON) = 2
    IF (KTYPE(JLON) == 2 .and. ZPBMPT >= YDML_PHY_EC%YRECUMF%RDEPTHS) KTYPE(JLON) = 1
    ICTOP0(JLON) = KCTOP(JLON)
  END IF
  ZRFL(JLON) = ZDMFUP(JLON, 1)
  DO JK=2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    ZRFL(JLON) = ZRFL(JLON) + ZDMFUP(JLON, JK)
  END DO
  DO JK=1,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    PMFD(JLON, JK) = 0.0_JPRB
    ZMFDS(JLON, JK) = 0.0_JPRB
    ZMFDQ(JLON, JK) = 0.0_JPRB
    ZDMFDP(JLON, JK) = 0.0_JPRB
    ZDPMEL(JLON, JK) = 0.0_JPRB
  END DO
  
  ! Needed for dissipation rate and if LMFUVDIS=T
  DO JK=1,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    ZTENU(JLON, JK) = PTENU(JLON, JK)
    ZTENV(JLON, JK) = PTENV(JLON, JK)
  END DO
  
  !-----------------------------------------------------------------------
  
  !*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
  !                  ------------------------------
  
  IF (YDML_PHY_EC%YRECUMF%LMFDD) THEN
    
    !*             (A) DETERMINE LFS IN 'CUDLFS'
    !                  -------------------------
    
    CALL CUDLFSN_OPENACC(YDTHF, YDCST, YDML_PHY_SLIN%YREPHLI, YDML_PHY_EC%YRECUMF, KIDIA, KFDIA, KLON, KLEV, KCBOT, KCTOP,  &
    & LDCUM, ZTENH, ZQENH, PTEN, PQSEN, PGEO, PGEOH, PAPH, PTU, PQU, ZMFUB, ZRFL, ZTD, ZQD, PMFD, ZMFDS, ZMFDQ, ZDMFDP, IDTOP,  &
    & LLDDRAF, YDSTACK=YLSTACK)
    
    !*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
    !                  -----------------------------------------------
    
    CALL CUDDRAFN_OPENACC(YDCST, YDTHF, YDML_PHY_SLIN%YREPHLI, YDML_PHY_EC%YRECUMF, KIDIA, KFDIA, KLON, KLEV, LLDDRAF, ZTENH,  &
    & ZQENH, PQSEN, PGEO, PGEOH, PAPH, ZRFL, ZTD, ZQD, PMFU, PMFD, ZMFDS, ZMFDQ, ZDMFDP, ZDMFDE, PMFDDE_RATE, ZKINED,  &
    & YDSTACK=YLSTACK)
    
  END IF
  
  !-----------------------------------------------------------------------
  
  !*    6.0          CLOSURE
  !                  ------
  
  !*                  RECALCULATE CLOUD BASE MASSFLUX FROM A
  !*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
  !*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
  !*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
  !                  --------------------------------------------
  
  !   DEEP CONVECTION
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  ZHEAT = 0.0_JPRB
  ZCAPE = 0.0_JPRB
  ZCAPE2 = 0.0_JPRB
  ZMFUB1 = ZMFUB(JLON)
  ZDQCV = 0.0_JPRB
  ZSATFR = 0.0_JPRB
  ZDX = 2*YDCST%RA*SQRT(YDCST%RPI*PGAW(JLON))
  ZDX = MAX(ZDX, 1.E2_JPRB)
  
  DO JK=YDML_PHY_EC%YRECUMF%NJKT2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    LLO1 = LDCUM(JLON) .and. KTYPE(JLON) == 1
    LLO3 = LLO1 .and. JK <= KCBOT(JLON) .and. JK > KCTOP(JLON)
    IF (LLO3) THEN
      ZDZ = (PGEO(JLON, JK - 1) - PGEO(JLON, JK))
      IF (LDTDKMF) THEN
        ZHEAT = ZHEAT + ((PTEN(JLON, JK - 1) - PTEN(JLON, JK) + ZDZ*ZORCPD) / ZTENH(JLON, JK) + YDCST%RETV*(PQEN(JLON, JK - 1) -  &
        & PQEN(JLON, JK)))*(YDCST%RG*(PMFU(JLON, JK) + PMFD(JLON, JK)))
      ELSE
        ZHEAT = ZHEAT + MAX(0.0_JPRB, ((PTEN(JLON, JK - 1) - PTEN(JLON, JK) + ZDZ*ZORCPD) / ZTENH(JLON, JK) +  &
        & YDCST%RETV*(PQEN(JLON, JK - 1) - PQEN(JLON, JK)))*(YDCST%RG*(PMFU(JLON, JK) + PMFD(JLON, JK))))
      END IF
      ZDZ = (PAP(JLON, JK) - PAP(JLON, JK - 1))
      ZCAPE = ZCAPE + ((PTU(JLON, JK) - ZTENH(JLON, JK)) / ZTENH(JLON, JK) + YDCST%RETV*(PQU(JLON, JK) - ZQENH(JLON, JK)) -  &
      & PLU(JLON, JK))*ZDZ        ! -PLRAIN(JL,JK) not added
    END IF
    IF (LLO3 .and. YDML_PHY_EC%YRECUMF%RCAPQADV > 0.0_JPRB) THEN
      ZTENH2(JLON, JK) = ZTENH(JLON, JK) - 0.5_JPRB*(PTENTA(JLON, JK) + PTENTA(JLON, JK - 1))*PTSPHY
      ZQENH2(JLON, JK) = ZQENH(JLON, JK) - 0.5_JPRB*(PTENQA(JLON, JK) + PTENQA(JLON, JK - 1))*PTSPHY
      ZCAPE2 = ZCAPE2 + ((PTU(JLON, JK) - ZTENH2(JLON, JK)) / ZTENH2(JLON, JK) + YDCST%RETV*(PQU(JLON, JK) - ZQENH2(JLON, JK)) -  &
      & PLU(JLON, JK))*ZDZ        ! -ZLRAIN(JL,JK) not added
    END IF
    ! IF(LLO1.AND.JK > KCTOP(JL)) THEN
    IF (LLO1) THEN
      ZDZ = (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
      ZDQCV = ZDQCV + PTENQA(JLON, JK)*ZDZ*(PQEN(JLON, JK) / PQSEN(JLON, JK))
      ! ZMFA=2.5_JPRB/MAX(2.5_JPRB,-PVERVEL(JL,JK))*PQSEN(JL,JK)/ZDX(JL)
      ! ZDQCV(JL)=ZDQCV(JL)+MIN(ZMFA,PTENQA(JL,JK))*ZDZ*(PQEN(JL,JK)/PQSEN(JL,JK))
    END IF
    IF (LLO1 .and. JK >= KCTOP(JLON)) THEN
      ZDZ = (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
      ZSATFR = ZSATFR + (PQEN(JLON, JK) / PQSEN(JLON, JK))*ZDZ
    END IF
  END DO
  ! time scale and subcloud contribution to CAPE to be substracted for better diurnal cycle over land
  !DIR$ LOOP_INFO EST_TRIPS(16)
  ZCAPDCYCL = 0.0_JPRB
  ZTAU = 0.0_JPRB
  IF (LDTDKMF) THEN
    ZDX = PDX(JLON)
    ZDX = MAX(ZDX, 1.E2_JPRB)
  END IF
  !! CFL limit as function of resolution used in DWD ICON model
  ! ZMFCFL(JL) = 2._JPRB*MIN(2._JPRB,1._JPRB + 2.5E-5_JPRB*ZDX(JL))
  IF (LDCUM(JLON) .and. KTYPE(JLON) == 1) THEN
    IF (YDML_PHY_EC%YRECUMF%LNEWTAU) THEN
      ZTAURES = YDML_PHY_EC%YRECUMF%RTAULIM*(1._JPRB + YDML_PHY_EC%YRECUMF%RTAUFAC / ((ZDX / 1000._JPRB -  &
      & YDML_PHY_EC%YRECUMF%RXMIN + 1._JPRB)**2))
    ELSE
      ZTAURES = 1.0_JPRB + 1.60_JPRB*ZDX / 125.E3_JPRB
      ! for the 10-1 km resolution range increase ZTAURES
      IF (ZDX < 8.E3_JPRB) THEN
        ZTAURES = 1.0_JPRB + LOG(8.E3_JPRB / ZDX)**2
      END IF
    END IF
    IF (ZDX > 125.E3_JPRB) ZTAURES = MIN(3.0_JPRB, ZTAURES)
    
    IKD = IDPL(JLON)
    IKB = KCBOT(JLON)
    IK = KCTOP(JLON)
    ZTAU = (PGEOH(JLON, IK) - PGEOH(JLON, IKB)) / ((2.0_JPRB + MIN(15.0_JPRB, PWMEAN(JLON)))*YDCST%RG) &
    & *ZTAURES*YDML_PHY_EC%YRECUMF%RTAUA
    ZXTAU = ZTAU
    
    LLO1 = PAPH(JLON, KLEV + 1) - PAPH(JLON, IKD) < 50.E2_JPRB
    IF (LLO1 .and. LDLAND(JLON) .and. YDML_PHY_EC%YRECUMF%RCAPDCYCL == 1.0_JPRB) THEN
      ZDZ = MIN(1.E4_JPRB, PGEOH(JLON, IKB) - PGEOH(JLON, KLEV + 1)) / YDCST%RG
      ZCAPDCYCL = ZXTAU*MAX(0.0_JPRB, ZKHVFL)*YDCST%RCPD / ZDZ
    END IF
    IF (LLO1 .and. YDML_PHY_EC%YRECUMF%RCAPDCYCL == 2.0_JPRB) THEN
      IF (LDLAND(JLON)) THEN
        ZCAPDCYCL = ZCAPPBL*ZTAU / ZTAURES
      ELSE
        ZDUTEN = 2.0_JPRB + SQRT(0.5*(PUEN(JLON, IKB)**2 + PVEN(JLON, IKB)**2 + PUEN(JLON, YDML_PHY_EC%YRECUMF%NJKT3)**2 +  &
        & PVEN(JLON, YDML_PHY_EC%YRECUMF%NJKT3)**2))
        ZTAUPBL = MIN(1.E4_JPRB, PGEOH(JLON, IKB) - PGEOH(JLON, KLEV + 1)) / (YDCST%RG*ZDUTEN)
        ZCAPDCYCL = ZCAPPBL*ZTAUPBL
      END IF
    END IF
    
    ZDQCV = ZDQCV*YDCST%RLVTT / PGEOH(JLON, IK)*ZXTAU / MAX(1.25_JPRB, ZTAURES)*YDML_PHY_EC%YRECUMF%RCAPQADV
    ZSATFR = ZSATFR / (PAPH(JLON, KLEV + 1) - PAPH(JLON, IK))
    
  ELSE
    ZTAU = 0.0_JPRB
  END IF
  
  IF (YDML_PHY_EC%YRECUMF%LMFCUCA) THEN
    !only allow cloud base mass flux to vary by certain amount
    !DIR$ LOOP_INFO EST_TRIPS(16)
    ZFACCA = MAX(0.3_JPRB, MIN(1.8_JPRB, PCUCONVCA(JLON)))
  ELSE
    ZFACCA = 1.0_JPRB
  END IF
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LDCUM(JLON) .and. KTYPE(JLON) == 1) THEN
    IKB = KCBOT(JLON)
    IF (LDTDKMF) THEN
      ZCAPE = MAX(YDML_PHY_EC%YRECUMF%RMINCAPE*ZCAPE, ZCAPE - ZCAPDCYCL)
    ELSE
      ZCAPE2 = YDML_PHY_EC%YRECUMF%RCAPQADV*ZCAPE2 + (1.0_JPRB - YDML_PHY_EC%YRECUMF%RCAPQADV)*ZCAPE
      ZCAPDCYCL = MAX(ZCAPDCYCL, -2*ZCAPE2)
      IF (ZSATFR <= 0.94 .or. PVERVEL(JLON, YDML_PHY_EC%YRECUMF%NJKT5) < -100._JPRB) THEN
        ZCAPE = MAX(YDML_PHY_EC%YRECUMF%RMINCAPE*ZCAPE, ZCAPE2 - ZCAPDCYCL + ZDQCV)
      ELSE
        ZCAPE = MAX(YDML_PHY_EC%YRECUMF%RMINCAPE*ZCAPE, ZCAPE - ZCAPDCYCL)
      END IF
    END IF
    ZCAPE = MIN(ZCAPE, 5000.0_JPRB)
    ZHEAT = MAX(1.E-4_JPRB, ZHEAT)
    !    ZXTAU(JL)=MAX(720._JPRB,ZXTAU(JL))
    ZXTAU = MAX(3.6E3_JPRB / (5.0_JPRB*PPLRG*PPLDARE), MIN(3.0_JPRB*3.6E3_JPRB / (PPLRG*PPLDARE), ZXTAU))
    ZMFUB1 = (ZCAPE*ZMFUB(JLON)) / (ZHEAT*ZXTAU)
    ZMFUB1 = MAX(ZMFUB1, 0.001_JPRB)*ZFACCA
    ZMFMAX = (PAPH(JLON, IKB) - PAPH(JLON, IKB - 1))*ZCONS2
    ZMFUB1 = MIN(ZMFUB1, ZMFMAX)
  END IF
  
  IF (LDMCAPEA) THEN
    ! Prefer ZCAPE computed from real CUASCN updraft rather than PCAPE estimated from CUBASEN.
    PCAPE(JLON) = YDCST%RG*ZCAPE
  END IF
  
  !  SHALLOW CONVECTION AND MID_LEVEL
  
  !DIR$ IVDEP
  !OCL NOVREC
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LDCUM(JLON) .and. KTYPE(JLON) >= 2) THEN
    IKB = KCBOT(JLON)
    IF (PMFD(JLON, IKB) < 0.0_JPRB) THEN
      ZEPS = -PMFD(JLON, IKB) / MAX(ZMFUB(JLON), 1.E-10_JPRB)
    ELSE
      ZEPS = 0.0_JPRB
    END IF
    ! maximum permissible value of ud base mass flux
    ZMFMAX = (PAPH(JLON, IKB) - PAPH(JLON, IKB - 1))*ZCONS2
    
    ! shallow convection
    
    IF (KTYPE(JLON) == 2) THEN
      IF (ZDHPBL > 0.0_JPRB) THEN
        ZQUMQE = PQU(JLON, IKB) + PLU(JLON, IKB) - ZEPS*ZQD(JLON, IKB) - (1.0_JPRB - ZEPS)*ZQENH(JLON, IKB)
        ZDH = YDCST%RCPD*(PTU(JLON, IKB) - ZEPS*ZTD(JLON, IKB) - (1.0_JPRB - ZEPS)*ZTENH(JLON, IKB)) + YDCST%RLVTT*ZQUMQE
        ZDH2 =  &
        & YDCST%RCPD*(PTU(JLON, IKB - 1) - ZEPS*ZTD(JLON, IKB - 1) - (1.0_JPRB - ZEPS)*ZTENH(JLON, IKB - 1)) + YDCST%RLVTT*ZQUMQE
        ZDH = 0.5_JPRB*(ZDH + ZDH2)
        ZDH = YDCST%RG*MAX(ZDH, 0.1_JPRB*YDCST%RCPD)
        ZMFUB1 = ZDHPBL / ZDH
        IF (.not.LDTDKMF) THEN
          IF (ZMFUB1 > 0.9_JPRB*ZMFMAX / YDML_PHY_EC%YRECUMF%RMFCFL .and. ZDH < 0.5_JPRB*YDCST%RG*YDCST%RCPD) THEN
            ZDH = 0.75_JPRB*YDCST%RCPD*YDCST%RG
            ZMFUB1 = ZDHPBL / ZDH
          END IF
        END IF
      ELSE
        ZMFUB1 = ZMFUB(JLON)
      END IF
      IF (YDML_PHY_EC%YRECUMF%LMFWSTAR) ZMFUB1 = ZMF_SHAL
      ZMFUB1 = MIN(ZMFUB1, ZMFMAX)
    END IF
    
    ! mid-level convection
    
    IF (KTYPE(JLON) == 3) THEN
      IF (LDTDKMF) THEN
        ZMFUB1 = MAX(ZMFUB(JLON), ZKHFL / PGEOH(JLON, IKB))*(1.0_JPRB + ZEPS)
      ELSE
        ZMFUB1 = MAX(ZMFUB(JLON), 5*ZKHFL / PGEOH(JLON, IKB))*(1.0_JPRB + ZEPS)
      END IF
      IF (PGEOH(JLON, IKB) < 1.5E3_JPRB) THEN
        ZMFUB1 = MIN(ZMFUB1, ZMFMAX / (1.05_JPRB*YDML_PHY_EC%YRECUMF%RMFCFL))
      ELSE
        ZMFUB1 = MIN(ZMFUB1, ZMFMAX)
      END IF
    END IF
    
  END IF
  
  ! rescale DD fluxes if deep and shallow convection
  
  DO JK=1,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LLDDRAF(JLON) .and. (KTYPE(JLON) == 1 .or. KTYPE(JLON) == 2)) THEN
      ZFAC = ZMFUB1 / MAX(ZMFUB(JLON), 1.E-10_JPRB)
      PMFD(JLON, JK) = PMFD(JLON, JK)*ZFAC
      ZMFDS(JLON, JK) = ZMFDS(JLON, JK)*ZFAC
      ZMFDQ(JLON, JK) = ZMFDQ(JLON, JK)*ZFAC
      ZDMFDP(JLON, JK) = ZDMFDP(JLON, JK)*ZFAC
      !  also rescale detrainment flux for ERA pp
      PMFDDE_RATE(JLON, JK) = PMFDDE_RATE(JLON, JK)*ZFAC
    END IF
  END DO
  
  !-----------------------------------------------------------------------
  
  !*    6.2          FINAL CLOSURE=SCALING
  !                  ---------------------
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LDCUM(JLON)) THEN
    ZMFS = ZMFUB1 / MAX(YDML_PHY_EC%YRECUMF%RMFCMIN, ZMFUB(JLON))
  END IF
  DO JK=2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
      IKB = KCBOT(JLON)
      IF (JK > IKB) THEN
        ZDZ = ((PAPH(JLON, KLEV + 1) - PAPH(JLON, JK)) / (PAPH(JLON, KLEV + 1) - PAPH(JLON, IKB)))
        PMFU(JLON, JK) = PMFU(JLON, IKB)*ZDZ
      END IF
      ZMFMAX = (PAPH(JLON, JK) - PAPH(JLON, JK - 1))*ZCONS2        !*zmfcfl(jl)/rmfcfl
      IF (.not.LDTDKMF) ZMFMAX = MIN(ZMFMAX, YDML_PHY_EC%YRECUMF%RMFLIA)
      IF (PMFU(JLON, JK)*ZMFS > ZMFMAX) THEN
        ZMFS = MIN(ZMFS, ZMFMAX / PMFU(JLON, JK))
        ZMFS = MAX(ZMFS, 1.E-10_JPRB)
      END IF
    END IF
  END DO
  
  !test additionally on divergence limit in case of direct Dynamics coupling
  IF (YDML_PHY_EC%YRECUMF%RMFADVW > 0.0_JPRB) THEN
    DO JK=2,KLEV - 1
      IF (KTYPE(JLON) == 1 .and. JK >= KCTOP(JLON) - 1) THEN
        ZFAC = ZMFS*YDCST%RG*YDML_PHY_EC%YRECUMF%RMFADVW*(PMFU(JLON, JK) - PMFU(JLON, JK + 1) +  &
        & YDML_PHY_EC%YRECUMF%RMFADVWDD*(PMFD(JLON, JK) - PMFD(JLON, JK + 1))) / (PAPH(JLON, JK) - PAPH(JLON, JK + 1))*PTSPHY
        IF (ABS(ZFAC) > 0.98_JPRB) ZMFS = ZMFS / MAX(1.02_JPRB, ABS(ZFAC))
      END IF
    END DO
    JK = KLEV
    IF (KTYPE(JLON) == 1) THEN
      ZFAC = ZMFS*YDCST%RG*YDML_PHY_EC%YRECUMF%RMFADVW*(PMFU(JLON, JK) + YDML_PHY_EC%YRECUMF%RMFADVWDD*PMFD(JLON, JK)) /  &
      & (PAPH(JLON, JK) - PAPH(JLON, JK + 1))*PTSPHY
      IF (ABS(ZFAC) > 0.98_JPRB) ZMFS = ZMFS / MAX(1.02_JPRB, ABS(ZFAC))
    END IF
  END IF
  
  DO JK=2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LDCUM(JLON) .and. JK <= KCBOT(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
      PMFU(JLON, JK) = PMFU(JLON, JK)*ZMFS
      ZMFUS(JLON, JK) = ZMFUS(JLON, JK)*ZMFS
      ZMFUQ(JLON, JK) = ZMFUQ(JLON, JK)*ZMFS
      ZMFUL(JLON, JK) = ZMFUL(JLON, JK)*ZMFS
      ZDMFUP(JLON, JK) = ZDMFUP(JLON, JK)*ZMFS
      ZDMFEN(JLON, JK) = ZDMFEN(JLON, JK)*ZMFS
      PLUDE(JLON, JK) = PLUDE(JLON, JK)*ZMFS
      PLUDELI(JLON, JK, 1) = PLUDELI(JLON, JK, 1)*ZMFS
      PLUDELI(JLON, JK, 2) = PLUDELI(JLON, JK, 2)*ZMFS
      PMFUDE_RATE(JLON, JK) = PMFUDE_RATE(JLON, JK)*ZMFS
    END IF
  END DO
  
  !-----------------------------------------------------------------------
  
  !*    6.5          IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
  !                  RESET LDCUM TO FALSE-> FLUXES SET TO ZERO IN CUFLXN
  !                  ---------------------------------------------------
  
  !                 exclude pathological KTYPE=2 KCBOT=KCTOP=KLEV-1
  
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (KTYPE(JLON) == 2 .and. KCBOT(JLON) == KCTOP(JLON) .and. KCBOT(JLON) >= KLEV - 1) THEN
    LDCUM(JLON) = .false.
    KTYPE(JLON) = 0
  END IF
  
  !                  turn off shallow convection if stratocumulus PBL type
  !DIR$ LOOP_INFO EST_TRIPS(16)
  LLO2 = .false.
  IF (.not.LDSHCV(JLON) .and. KTYPE(JLON) == 2) THEN
    LLO2 = .true.
    LDCUM(JLON) = .false.
  END IF
  
  
  IF (YDML_PHY_EC%YREPHY%LELIGHT .or. YDCHEM%LCHEM_LIGHT) THEN
    !DIR$ LOOP_INFO EST_TRIPS(16)
    ! The following is needed to allow lightning calculations in case deep convection
    ! is treated explicitly (i.e. when LMFPEN=false).
    ! Store convection flag.
    LDCUM_LIG(JLON) = LDCUM(JLON)
    ! Store convective cloud top and bottom levels.
    KCTOP_LIG(JLON) = KCTOP(JLON)
    KCBOT_LIG(JLON) = KCBOT(JLON)
  END IF
  
  IF (.not.YDML_PHY_EC%YRECUMF%LMFSCV .or. .not.YDML_PHY_EC%YRECUMF%LMFPEN) THEN
    !DIR$ LOOP_INFO EST_TRIPS(16)
    ! LLO2(JL)=.FALSE.
    IF (.not.YDML_PHY_EC%YRECUMF%LMFSCV .and. KTYPE(JLON) == 2 .or. .not.YDML_PHY_EC%YRECUMF%LMFPEN .and. KTYPE(JLON) == 1) THEN
      LLO2 = .true.
      LDCUM(JLON) = .false.
    END IF
  END IF
  
  !-----------------------------------------------------------------------
  
  !*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
  !                  ------------------------------------------
  
  !- set DD mass fluxes to zero above cloud top
  !  (because of inconsistency with second updraught)
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LLDDRAF(JLON) .and. IDTOP(JLON) <= KCTOP(JLON)) THEN
    IDTOP(JLON) = KCTOP(JLON) + 1
  END IF
  
  DO JK=2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LLDDRAF(JLON)) THEN
      IF (JK < IDTOP(JLON)) THEN
        PMFD(JLON, JK) = 0.0_JPRB
        ZMFDS(JLON, JK) = 0.0_JPRB
        ZMFDQ(JLON, JK) = 0.0_JPRB
        PMFDDE_RATE(JLON, JK) = 0.0_JPRB
        ZDMFDP(JLON, JK) = 0.0_JPRB
      ELSE IF (JK == IDTOP(JLON)) THEN
        PMFDDE_RATE(JLON, JK) = 0.0_JPRB
      END IF
    END IF
  END DO
  CALL CUFLXN_OPENACC(YDTHF, YDCST, YDML_PHY_SLIN%YREPHLI, YDML_PHY_EC%YRECUMF, KIDIA, KFDIA, KLON, KLEV, PTSPHY, PTEN, PQEN,  &
  & PQSEN, ZTENH, ZQENH, PAPH, PAP, PGEOH, LDLAND, LDCUM, LDTDKMF, KCBOT, KCTOP, IDTOP, ITOPM2, KTYPE, LLDDRAF, PMFU, PMFD,  &
  & ZMFUS, ZMFDS, ZMFUQ, ZMFDQ, ZMFUL, PLUDE, PLUDELI, PLRAIN, PSNDE, ZDMFUP, ZDMFDP, ZDPMEL, ZLGLAC, PMFLXR, PMFLXS, PRAIN,  &
  & PMFUDE_RATE, PMFDDE_RATE, YDSTACK=YLSTACK)
  
  !- correct DD detrainment rates if entrainment becomes negative
  !- correct UD detrainment rates if entrainment becomes negative
  !- conservation correction for precip
  
  DO JK=2,KLEV - 1
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LLDDRAF(JLON) .and. JK >= IDTOP(JLON) - 1) THEN
      ZERATE = -PMFD(JLON, JK) + PMFD(JLON, JK - 1) + PMFDDE_RATE(JLON, JK)
      IF (ZERATE < 0.0_JPRB) THEN
        PMFDDE_RATE(JLON, JK) = PMFDDE_RATE(JLON, JK) - ZERATE
      END IF
    END IF
    IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
      ZERATE = PMFU(JLON, JK) - PMFU(JLON, JK + 1) + PMFUDE_RATE(JLON, JK)
      IF (ZERATE < 0.0_JPRB) THEN
        PMFUDE_RATE(JLON, JK) = PMFUDE_RATE(JLON, JK) - ZERATE
      END IF
      ! ZDMFUP(JL,JK)=ZDMFUP(JL,JK)+ZDMFDP(JL,JK)
      ZDMFUP(JLON, JK) = PMFLXR(JLON, JK + 1) + PMFLXS(JLON, JK + 1) - PMFLXR(JLON, JK) - PMFLXS(JLON, JK)
      ZDMFDP(JLON, JK) = 0.0_JPRB
    END IF
  END DO
  
  ! avoid negative humidities at ddraught top
  !DIR$ LOOP_INFO EST_TRIPS(16)
  IF (LLDDRAF(JLON)) THEN
    JK = IDTOP(JLON)
    IK = MIN(JK + 1, KLEV)
    IF (ZMFDQ(JLON, JK) < 0.3_JPRB*ZMFDQ(JLON, IK)) THEN
      IF (YDML_PHY_EC%YRECUMF%RMFSOLTQ == 0.0_JPRB) THEN
        ZMFDQ(JLON, JK) = 0.3_JPRB*ZMFDQ(JLON, IK)
      ELSE
        PMFD(JLON, JK) = 0.3_JPRB*PMFD(JLON, IK)
      END IF
    END IF
  END IF
  
  ! avoid negative humidities near cloud top because gradient of precip flux
  ! and detrainment / liquid water flux too large
  DO JK=2,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1 .and. JK < KCBOT(JLON)) THEN
      ZDZ = PTSPHY*YDCST%RG / (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
      ZMFA = ZMFUQ(JLON, JK + 1) + ZMFDQ(JLON, JK + 1) - ZMFUQ(JLON, JK) - ZMFDQ(JLON, JK) + ZMFUL(JLON, JK + 1) - ZMFUL(JLON, JK &
      & ) + ZDMFUP(JLON, JK) - PSNDE(JLON, JK, 1) - PSNDE(JLON, JK, 2)
      ZMFA = (ZMFA - PLUDE(JLON, JK))*ZDZ
      IF (PQEN(JLON, JK) + ZMFA < 0.0_JPRB .and. PLUDE(JLON, JK) > 0.0_JPRB) THEN
        ZFAC = MIN(-(PQEN(JLON, JK) + ZMFA) / ZDZ, PLUDE(JLON, JK))
        PLUDE(JLON, JK) = PLUDE(JLON, JK) + ZFAC
        PLUDELI(JLON, JK, 1) = PLUDELI(JLON, JK, 1) + ZFAC*PLUDELI(JLON, JK, 1) / PLUDE(JLON, JK)
        PLUDELI(JLON, JK, 2) = PLUDELI(JLON, JK, 2) + ZFAC*PLUDELI(JLON, JK, 2) / PLUDE(JLON, JK)
      END IF
      IF (PLUDE(JLON, JK) < 0.0_JPRB) THEN
        PLUDE(JLON, JK) = 0.0_JPRB
        PLUDELI(JLON, JK, 1) = 0.0_JPRB
        PLUDELI(JLON, JK, 2) = 0.0_JPRB
      END IF
    END IF
    IF (.not.LDCUM(JLON)) THEN
      PMFUDE_RATE(JLON, JK) = 0.0_JPRB
    END IF
    IF (PMFD(JLON, JK - 1) == 0.0_JPRB) THEN
      PMFDDE_RATE(JLON, JK) = 0.0_JPRB
    END IF
  END DO
  !----------------------------------------------------------------------
  
  !*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
  !                  --------------------------------------------------
  
  IF (YDML_PHY_EC%YRECUMF%RMFSOLTQ > 0.0_JPRB) THEN
    ! derive draught properties for implicit
    
    DO JK=KLEV,2,-1
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LDCUM(JLON)) THEN
        IF (JK > KCBOT(JLON)) THEN
          ZMFA = 1.0_JPRB / MAX(1.E-15_JPRB, PMFU(JLON, JK))
          PQU(JLON, JK) = ZQENH(JLON, JK) + ZMFUQ(JLON, JK)*ZMFA
          PTU(JLON, JK) = ZTENH(JLON, JK) + ZMFUS(JLON, JK)*ZMFA*ZORCPD
          ZMFUS(JLON, JK) = PMFU(JLON, JK)*(YDCST%RCPD*PTU(JLON, JK) + PGEOH(JLON, JK))
          ZMFUQ(JLON, JK) = PMFU(JLON, JK)*PQU(JLON, JK)
          IF (LLDDRAF(JLON)) THEN
            ZMFA = 1.0_JPRB / MIN(-1.E-15_JPRB, PMFD(JLON, JK))
            ZQD(JLON, JK) = ZQENH(JLON, JK) + ZMFDQ(JLON, JK)*ZMFA
            ZTD(JLON, JK) = ZTENH(JLON, JK) + ZMFDS(JLON, JK)*ZMFA*ZORCPD
            ZMFDQ(JLON, JK) = PMFD(JLON, JK)*ZQD(JLON, JK)
            ZMFDS(JLON, JK) = PMFD(JLON, JK)*(YDCST%RCPD*ZTD(JLON, JK) + PGEOH(JLON, JK))
          END IF
        ELSE IF (JK <= KCBOT(JLON) .and. JK >= KCTOP(JLON)) THEN
          ZMFUS(JLON, JK) = PMFU(JLON, JK)*(YDCST%RCPD*PTU(JLON, JK) + PGEOH(JLON, JK))
          ZMFUQ(JLON, JK) = PMFU(JLON, JK)*PQU(JLON, JK)
          ZMFDS(JLON, JK) = PMFD(JLON, JK)*(YDCST%RCPD*ZTD(JLON, JK) + PGEOH(JLON, JK))
          ZMFDQ(JLON, JK) = PMFD(JLON, JK)*ZQD(JLON, JK)
        END IF
      END IF
    END DO
    
  END IF
  
  CALL CUDTDQN_OPENACC(YDTHF, YDCST, YDML_PHY_SLIN%YREPHLI, YDML_PHY_SLIN%YRPHNC, YDML_PHY_EC%YRECUMF, YDML_PHY_EC%YREPHY,  &
  & KIDIA, KFDIA, KLON, KLEV, ITOPM2, KTYPE, KCTOP, KCBOT, IDTOP, LDTDKMF, LDCUM, LLDDRAF, LLSCVFLAG, PTSPHY, PAPH, PAP, PGEOH,  &
  & PGEO, PTEN, ZTENH, PQEN, ZQENH, PQSEN, ZLGLAC, PLUDE, PLUDELI, PSNDE, PMFU, PMFD, ZMFUS, ZMFDS, ZMFUQ, ZMFDQ, ZMFUL, ZDMFUP,  &
  & ZDPMEL, PMFLXR, PMFLXS, PTENT, PTENQ, PENTH, YDSTACK=YLSTACK)
  
  !----------------------------------------------------------------------
  
  !*    9.0          COMPUTE MOMENTUM IN UPDRAUGHT AND DOWNDRAUGHT
  !                  ---------------------------------------------
  
  IF (YDML_PHY_EC%YRECUMF%LMFDUDV) THEN
    
    DO JK=KLEV - 1,2,-1
      IK = JK + 1
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LDCUM(JLON)) THEN
        IF (JK == KCBOT(JLON) .and. KTYPE(JLON) < 3) THEN
          IKB = IDPL(JLON)
          ZUU(JLON, JK) = PUEN(JLON, IKB - 1)
          ZVU(JLON, JK) = PVEN(JLON, IKB - 1)
        ELSE IF (JK == KCBOT(JLON) .and. KTYPE(JLON) == 3) THEN
          ZUU(JLON, JK) = PUEN(JLON, JK - 1)
          ZVU(JLON, JK) = PVEN(JLON, JK - 1)
        END IF
        IF (JK < KCBOT(JLON) .and. JK >= KCTOP(JLON)) THEN
          ZFAC = 0.0_JPRB
          IF (LDTDKMF .and. (KTYPE(JLON) == 1 .or. KTYPE(JLON) == 3)) ZFAC = 2.0_JPRB
          ! IF(KTYPE(JL)==2) ZFAC=1.0_JPRB
          ZERATE = PMFU(JLON, JK) - PMFU(JLON, IK) + (1.0_JPRB + ZFAC)*PMFUDE_RATE(JLON, JK)
          ZDERATE = (1.0_JPRB + ZFAC)*PMFUDE_RATE(JLON, JK)
          ZMFA = 1.0_JPRB / MAX(YDML_PHY_EC%YRECUMF%RMFCMIN, PMFU(JLON, JK))
          ZUU(JLON, JK) = (ZUU(JLON, IK)*PMFU(JLON, IK) + ZERATE*PUEN(JLON, JK) - ZDERATE*ZUU(JLON, IK))*ZMFA
          ZVU(JLON, JK) = (ZVU(JLON, IK)*PMFU(JLON, IK) + ZERATE*PVEN(JLON, JK) - ZDERATE*ZVU(JLON, IK))*ZMFA
        END IF
      END IF
    END DO
    DO JK=3,KLEV
      IK = JK - 1
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LDCUM(JLON)) THEN
        IF (JK == IDTOP(JLON)) THEN
          ZUD(JLON, JK) = 0.5_JPRB*(ZUU(JLON, JK) + PUEN(JLON, IK))
          ZVD(JLON, JK) = 0.5_JPRB*(ZVU(JLON, JK) + PVEN(JLON, IK))
        ELSE IF (JK > IDTOP(JLON)) THEN
          ZERATE = -PMFD(JLON, JK) + PMFD(JLON, IK) + PMFDDE_RATE(JLON, JK)
          ZMFA = 1.0_JPRB / MIN(-YDML_PHY_EC%YRECUMF%RMFCMIN, PMFD(JLON, JK))
          ZUD(JLON, JK) = (ZUD(JLON, IK)*PMFD(JLON, IK) - ZERATE*PUEN(JLON, IK) + PMFDDE_RATE(JLON, JK)*ZUD(JLON, IK))*ZMFA
          ZVD(JLON, JK) = (ZVD(JLON, IK)*PMFD(JLON, IK) - ZERATE*PVEN(JLON, IK) + PMFDDE_RATE(JLON, JK)*ZVD(JLON, IK))*ZMFA
        END IF
      END IF
    END DO
    
    !*    9.1          UPDATE TENDENCIES FOR U AND V IN SUBROUTINE CUDUDV
    !                  --------------------------------------------------
    
    ! for explicit/semi-implicit rescale massfluxes for stability in Momentum
    !------------------------------------------------------------------------
    
    ZMFS = 1.0_JPRB
    ! IF(RMFSOLUV<=0.5_JPRB) THEN
    IF (YDML_PHY_EC%YRECUMF%RMFSOLUV <= 1.0_JPRB) THEN
      DO JK=2,KLEV
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
          ZMFMAX = (PAPH(JLON, JK) - PAPH(JLON, JK - 1))*ZCONS
          IF (PMFU(JLON, JK) > ZMFMAX .and. JK >= KCTOP(JLON)) ZMFS = MIN(ZMFS, ZMFMAX / PMFU(JLON, JK))
        END IF
      END DO
    END IF
    DO JK=1,KLEV
      !DIR$ LOOP_INFO EST_TRIPS(16)
      ZMFUUS(JLON, JK) = PMFU(JLON, JK)
      ZMFDUS(JLON, JK) = PMFD(JLON, JK)
      IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
        ZMFUUS(JLON, JK) = PMFU(JLON, JK)*ZMFS
        ZMFDUS(JLON, JK) = PMFD(JLON, JK)*ZMFS
      END IF
    END DO
    
    ! recompute Draught properties below for Implicit
    ! based on linear flux profiles
    
    IF (YDML_PHY_EC%YRECUMF%RMFSOLUV > 0.0_JPRB) THEN
      
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LDCUM(JLON)) THEN
        JK = KCBOT(JLON)
        IK = JK - 1
        ZMFUUB = ZMFUUS(JLON, JK)*(ZUU(JLON, JK) - PUEN(JLON, IK))
        ZMFUVB = ZMFUUS(JLON, JK)*(ZVU(JLON, JK) - PVEN(JLON, IK))
      END IF
      
      DO JK=2,KLEV
        IK = JK - 1
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LDCUM(JLON) .and. JK > KCBOT(JLON)) THEN
          IKB = KCBOT(JLON)
          ZDZ = ((PAPH(JLON, KLEV + 1) - PAPH(JLON, JK)) / (PAPH(JLON, KLEV + 1) - PAPH(JLON, IKB)))
          IF (KTYPE(JLON) == 3) THEN
            ZDZ = ZDZ*ZDZ
          END IF
          ZMFA = 1.0_JPRB / MAX(YDML_PHY_EC%YRECUMF%RMFCMIN, ZMFUUS(JLON, JK))
          ZUU(JLON, JK) = PUEN(JLON, IK) + ZMFUUB*ZDZ*ZMFA
          ZVU(JLON, JK) = PVEN(JLON, IK) + ZMFUVB*ZDZ*ZMFA
          
          ZMFDUS(JLON, JK) = ZMFDUS(JLON, IKB)*ZDZ
          ZUD(JLON, JK) = PUEN(JLON, IK) + ZUD(JLON, IKB) - PUEN(JLON, IKB - 1)
          ZVD(JLON, JK) = PVEN(JLON, IK) + ZVD(JLON, IKB) - PVEN(JLON, IKB - 1)
        END IF
        ! add UV perturb to correct wind bias
        IF (LDCUM(JLON) .and. JK >= KCTOP(JLON)) THEN
          IF (LDTDKMF) THEN
            ZUU(JLON, JK) = ZUU(JLON, JK) - YDML_PHY_EC%YRECUMF%RUVPER*SIGN(1.0_JPRB, ZUU(JLON, JK))
            ZVU(JLON, JK) = ZVU(JLON, JK) - YDML_PHY_EC%YRECUMF%RUVPER*SIGN(1.0_JPRB, ZVU(JLON, JK))
          ELSE
            ZUU(JLON, JK) = ZUU(JLON, JK) - MIN(ABS(ZUU(JLON, JK)), YDML_PHY_EC%YRECUMF%RUVPER)*SIGN(1.0_JPRB, ZUU(JLON, JK))
            ZVU(JLON, JK) = ZVU(JLON, JK) - MIN(ABS(ZVU(JLON, JK)), YDML_PHY_EC%YRECUMF%RUVPER)*SIGN(1.0_JPRB, ZVU(JLON, JK))
          END IF
        END IF
      END DO
      
    END IF
    
    !-------------------------------------------------------------------
    ! End
    ! Intermediate Solution for stability in EPS:
    ! For original code replace line
    !  &, PUEN,     PVEN,     ZMFUUS,   ZMFDUS &
    !by
    !  &, PUEN,     PVEN,     PMFU,     PMFD
    
    CALL CUDUDV_OPENACC(YDCST, YDML_PHY_EC%YRECUMF, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, ITOPM2, KTYPE, KCBOT, KCTOP, LDCUM,  &
    & PTSPHY, PAPH, PAP, PUEN, PVEN, ZMFUUS, ZMFDUS, PMFU, PMFD, ZUU, ZUD, ZVU, ZVD, PGP2DSPP, PTENU, PTENV, YDSTACK=YLSTACK)
    
    IF (YDML_PHY_EC%YRECUMF%LMFUVDIS) THEN
      ! add KE dissipation
      !DIR$ LOOP_INFO EST_TRIPS(16)
      ZSUM12 = 0.0_JPRB
      ZSUM22 = 0.0_JPRB
      DO JK=1,KLEV
        !DIR$ LOOP_INFO EST_TRIPS(16)
        ZUV2(JLON, JK) = 0.0_JPRB
        IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
          ZDZ = (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
          ZDUTEN = PTENU(JLON, JK) - ZTENU(JLON, JK)
          ZDVTEN = PTENV(JLON, JK) - ZTENV(JLON, JK)
          ZUV2(JLON, JK) = SQRT(ZDUTEN**2 + ZDVTEN**2)
          ZSUM22 = ZSUM22 + ZUV2(JLON, JK)*ZDZ
          ZSUM12 = ZSUM12 - (PUEN(JLON, JK)*ZDUTEN + PVEN(JLON, JK)*ZDVTEN)*ZDZ
        END IF
      END DO
      ! Store vertically integrated dissipation rate (W m-2)
      !DIR$ LOOP_INFO EST_TRIPS(16)
      PVDISCU(JLON) = ZSUM12*ZRG
      ! Calculate local temperature tendency (K) due to dissipation
      DO JK=1,KLEV
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
          ZTDIS = ZORCPD*ZSUM12*ZUV2(JLON, JK) / MAX(1.E-15_JPRB, ZSUM22)
          PTENT(JLON, JK) = PTENT(JLON, JK) + ZTDIS
        END IF
      END DO
    END IF
    
  END IF
  
  !----------------------------------------------------------------------
  
  !*   10.           IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
  !                  NEED TO SET SOME VARIABLES A POSTERIORI TO ZERO
  !                  ---------------------------------------------------
  
  IF (.not.YDML_PHY_EC%YRECUMF%LMFSCV .or. .not.YDML_PHY_EC%YRECUMF%LMFPEN) THEN
    DO JK=2,KLEV
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LLO2 .and. JK >= KCTOP(JLON) - 1) THEN
        PMFUDE_RATE(JLON, JK) = 0.0_JPRB
        PMFDDE_RATE(JLON, JK) = 0.0_JPRB
      END IF
    END DO
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IF (LLO2) THEN
      KCTOP(JLON) = KLEV - 1
      KCBOT(JLON) = KLEV - 1
    END IF
  END IF
  
  !----------------------------------------------------------------------
  
  !*   11.0          CHEMICAL TRACER TRANSPORT
  !                  -------------------------
  
  IF (YDML_PHY_EC%YREPHY%LMFTRAC .and. KTRAC > 0) THEN
    
    ! transport switched off for mid-level convection
    !DIR$ LOOP_INFO EST_TRIPS(16)
    !IF( LDCUM(JL).AND.KTYPE(JL)/=3 ) THEN
    IF (LDCUM(JLON) .and. KTYPE(JLON) /= 3 .and. KCBOT(JLON) - KCTOP(JLON) >= 1) THEN
      LLDCUM(JLON) = .true.
      LLDDRAF3(JLON) = LLDDRAF(JLON)
    ELSE
      LLDCUM(JLON) = .false.
      LLDDRAF3(JLON) = .false.
    END IF
    
    ! check and correct mass fluxes for CFL criterium
    
    ZMFS = 1.0_JPRB
    IF (YDML_PHY_EC%YRECUMF%RMFSOLCT <= 3.0_JPRB) THEN
      DO JK=2,KLEV
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LLDCUM(JLON) .and. JK >= KCTOP(JLON)) THEN
          ZMFMAX = (PAPH(JLON, JK) - PAPH(JLON, JK - 1))*0.8_JPRB*ZCONS
          IF (PMFU(JLON, JK) > ZMFMAX) ZMFS = MIN(ZMFS, ZMFMAX / PMFU(JLON, JK))
        END IF
      END DO
    END IF
    DO JK=1,KLEV
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LLDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
        ZMFUUS(JLON, JK) = PMFU(JLON, JK)*ZMFS
        ZMFUDR(JLON, JK) = PMFUDE_RATE(JLON, JK)*ZMFS
        ZDMFUPC(JLON, JK) = (PMFLXR(JLON, JK) + PMFLXS(JLON, JK))*ZMFS
        IKB = KCBOT(JLON)
        IF (JK > IKB) THEN
          ZDMFUPC(JLON, JK) = ZDMFUPC(JLON, IKB)
        END IF
      ELSE
        ZMFUUS(JLON, JK) = 0._JPRB
        ZMFUDR(JLON, JK) = 0._JPRB
      END IF
      IF (LLDDRAF3(JLON) .and. JK >= IDTOP(JLON) - 1) THEN
        ZMFDUS(JLON, JK) = PMFD(JLON, JK)*ZMFS
        ZMFDDR(JLON, JK) = PMFDDE_RATE(JLON, JK)*ZMFS
        !     ZDMFDPC(JL,JK)=ZDMFDPC(JL,JK)*ZMFS(JL)
        ZDMFDPC(JLON, JK) = 0.0_JPRB
      ELSE
        ZMFDUS(JLON, JK) = 0._JPRB
        ZMFDDR(JLON, JK) = 0._JPRB
      END IF
    END DO
    
    IF (YDML_PHY_EC%YRECUMF%LMFSMOOTH) THEN
      ! smmoothing of mass fluxes (gradients) at top and bottom of draughts
      DO JK=2,KLEV - 1
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LLDDRAF3(JLON) .and. ZMFDUS(JLON, JK) < 0.0_JPRB .and. ZMFDUS(JLON, JK + 1) == 0.0_JPRB) THEN
          ZERATE = MIN(0._JPRB, ZMFDUS(JLON, JK) - 0.5*ZMFDUS(JLON, JK - 1))
          ZMFDUS(JLON, JK) = ZMFDUS(JLON, JK) - ZERATE
          ZMFDDR(JLON, JK) = ZMFDDR(JLON, JK) - ZERATE
          ZMFDDR(JLON, JK + 1) = -ZMFDUS(JLON, JK)
        END IF
        IF (LLDCUM(JLON) .and. JK == KCTOP(JLON)) THEN
          ZERATE = MAX(0.0_JPRB, ZMFUUS(JLON, JK) - 0.5_JPRB*ZMFUUS(JLON, JK + 1))
          ZMFUUS(JLON, JK) = ZMFUUS(JLON, JK) - ZERATE
          ZMFUDR(JLON, JK) = ZMFUDR(JLON, JK) + ZERATE
          ZMFUDR(JLON, JK - 1) = ZMFUUS(JLON, JK)
        END IF
      END DO
      DO JK=KLEV - 1,2,-1
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LLDCUM(JLON)) THEN
          IF (ZMFUDR(JLON, JK) == 0.0_JPRB .and. ZMFUDR(JLON, JK - 1) > 0.0_JPRB) THEN
            ZMFUDR(JLON, JK) = 0.5_JPRB*ZMFUDR(JLON, JK - 1)
          END IF
        END IF
      END DO
    END IF
    
    IF (YDML_PHY_EC%YREPHY%LMFSCAV) THEN
      CALL CUCTRACER_OPENACC(YDCST, YDML_PHY_SLIN%YRCUMFS, YDML_PHY_SLIN%YRECUMF2, YDML_PHY_EC%YRECUMF, KIDIA, KFDIA, KLON,  &
      & KLEV, KTRAC, KCTOP, IDTOP, KTYPE, LLDCUM, LLDDRAF3, PTSPHY, PAPH, PAP, ZMFUUS, ZMFDUS, PMFU, PMFD, ZMFUDR, ZMFDDR,  &
      & ZDMFUPC, ZDMFDPC, PCEN, PTENC, PSCAV, YDSTACK=YLSTACK)
    ELSE
      ZDMFUPC(JLON, :) = 0.0_JPRB
      CALL CUCTRACER_OPENACC(YDCST, YDML_PHY_SLIN%YRCUMFS, YDML_PHY_SLIN%YRECUMF2, YDML_PHY_EC%YRECUMF, KIDIA, KFDIA, KLON,  &
      & KLEV, KTRAC, KCTOP, IDTOP, KTYPE, LLDCUM, LLDDRAF3, PTSPHY, PAPH, PAP, ZMFUUS, ZMFDUS, PMFU, PMFD, ZMFUDR, ZMFDDR,  &
      & ZDMFUPC, ZDMFDPC, PCEN, PTENC, PSCAV0, YDSTACK=YLSTACK)
    END IF
    
  END IF
  
  !----------------------------------------------------------------------
  
  !*   12.           PUT DETRAINMENT RATES FROM MFLX UNITS IN UNITS MFLX/M
  !                  FOR ERA40, ESTIMATE VOLUME MEAN RAIN AND SNOW CONTENT
  !                  ---------------------------------------------------
  
  PDISS(JLON, 1) = 0.0_JPRB
  ZAR = 1.0_JPRB / 20.89_JPRB
  ZAS = 1.0_JPRB / 29.51_JPRB
  ZBR = 1.0_JPRB / 1.15_JPRB
  ZBS = 1.0_JPRB / 1.10_JPRB
  ZAK = 1.0_JPRB / 5.09E-3_JPRB
  DO JK=1,KLEV
    !DIR$ LOOP_INFO EST_TRIPS(16)
    IK = MIN(JK + 1, KLEV)
    IKD = MAX(JK - 1, 1)
    PDISS(JLON, JK) = 0.0_JPRB
    PRSUD(JLON, JK, 1) = 0.0_JPRB
    PRSUD(JLON, JK, 2) = 0.0_JPRB
    IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
      PLUDELI(JLON, JK, 3) =  &
      & PMFUDE_RATE(JLON, JK)*(PQU(JLON, IK) - ZQENH(JLON, IK)) + PMFDDE_RATE(JLON, JK)*(ZQD(JLON, IKD) - ZQENH(JLON, IKD))
      PLUDELI(JLON, JK, 4) =  &
      & PMFUDE_RATE(JLON, JK)*(PTU(JLON, IK) - ZTENH(JLON, IK)) + PMFDDE_RATE(JLON, JK)*(ZTD(JLON, IKD) - ZTENH(JLON, IKD))
      ! Change units for detrainment as for ERA
      ZRO = YDCST%RG / (PGEOH(JLON, JK) - PGEOH(JLON, JK + 1))        ! 1/dz
      PMFUDE_RATE(JLON, JK) = PMFUDE_RATE(JLON, JK)*ZRO
      PMFDDE_RATE(JLON, JK) = PMFDDE_RATE(JLON, JK)*ZRO
      ZRO = YDCST%RD*PTEN(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQEN(JLON, JK)) / PAP(JLON, JK)
      ! Dissipation definition used in stoch Backscatter
      ! PDISS(JL,JK)=MIN(PWMEAN(JL),17._JPRB)**2*PMFUDE_RATE(JL,JK)*ZRO*ZTAU(JL)
      ZDUTEN = PTENU(JLON, JK) - ZTENU(JLON, JK)
      ZDVTEN = PTENV(JLON, JK) - ZTENV(JLON, JK)
      ! Dissipation used in eddy dissipation(EDR) turbulence diagnostic
      PDISS(JLON, JK) = ABS(PUEN(JLON, JK)*ZDUTEN + PVEN(JLON, JK)*ZDVTEN)**0.3333
      
      ! grid-mean convective rain/snow parametrization (Geer et al. 2009) with factor 0.5 for snow
      PRSUD(JLON, JK, 1) = 1.E-3_JPRB*ZRO*(PMFLXR(JLON, JK)*3600._JPRB*ZAR)**ZBR
      PRSUD(JLON, JK, 2) = 0.5*1.E-3_JPRB*ZRO*(10.0_JPRB*PMFLXS(JLON, JK)*3600._JPRB*ZAS)**ZBS
      ! possible updraught fraction for incloud values
      ! ZFAC=PAPH(JL,NJKT5)/MIN(PAPH(JL,NJKT5),MAX(200.E2_JPRB,PAPH(JL,JK)))
      ! ZFAC=RCUCOV*ZFAC**2
      ! grid-mean convective rain/snow parametrization similar to that used in evaporation
      ! i.e. Kessler but without pressure factor and factor for snow
      ! PRSUD(JL,JK,1)=1.E-3_JPRB*ZRO*(PMFLXR(JL,JK)*ZAK)**0.888_JPRB
      ! PRSUD(JL,JK,2)=1.E-3_JPRB*ZRO*(PMFLXS(JL,JK)*2.5*ZAK)**0.888_JPRB
    ELSE
      PMFU(JLON, JK) = 0.0_JPRB
      PMFD(JLON, JK) = 0.0_JPRB
      PLUDE(JLON, JK) = 0.0_JPRB
      PLUDELI(JLON, JK, 1) = 0.0_JPRB
      PLUDELI(JLON, JK, 2) = 0.0_JPRB
      PLUDELI(JLON, JK, 3) = 0.0_JPRB
      PLUDELI(JLON, JK, 4) = 0.0_JPRB
      PSNDE(JLON, JK, 1) = 0.0_JPRB
      PSNDE(JLON, JK, 2) = 0.0_JPRB
      PTU(JLON, JK) = PTEN(JLON, JK)
      PQU(JLON, JK) = PQEN(JLON, JK)
      PLU(JLON, JK) = 0.0_JPRB
      PENTH(JLON, JK) = 0.0_JPRB
      PMFUDE_RATE(JLON, JK) = 0.0_JPRB
      PMFDDE_RATE(JLON, JK) = 0.0_JPRB
    END IF
  END DO
  
  !----------------------------------------------------------------------
  
END SUBROUTINE CUMASTRN_OPENACC
