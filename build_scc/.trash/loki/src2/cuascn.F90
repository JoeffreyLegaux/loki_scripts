SUBROUTINE CUASCN_OPENACC (YDTHF, YDCST, YDEPHLI, YDECLDP, YDECUMF, YDSPP_CONFIG, YGFL, KIDIA, KFDIA, KLON, KLEV, LDTDKMF,  &
& PTSPHY, PTENH, PQENH, PUEN, PVEN, PTEN, PQEN, PQSEN, PLITOT, PGEO, PGEOH, PAP, PAPH, PVERVEL, PWUBASE, PGP2DSPP, LDLAND,  &
& LDCUM, KTYPE, KLAB, LSCVFLAG, PTU, PQU, PLU, PLRAIN, PMFU, PMFUB, PLGLAC, PMFUS, PMFUQ, PMFUL, PLUDE, PLUDELI, PDMFUP,  &
& PLCRIT_AER, PDMFEN, KCBOT, KCTOP, KCTOP0, KDPL, PMFUDE_RATE, PKINEU, PWMEAN, YDSTACK)
  
  !          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
  !          FOR CUMULUS PARAMETERIZATION
  
  !          PURPOSE.
  !          --------
  !          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
  !          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
  !           FLUXES AS WELL AS PRECIPITATION RATES)
  
  !          INTERFACE
  !          ---------
  
  !          THIS ROUTINE IS CALLED FROM *CUMASTR*.
  
  !          METHOD.
  !          --------
  !          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
  !          AND THEN CALCULATE MOIST ASCENT FOR
  !          ENTRAINING/DETRAINING PLUME.
  !          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
  !          SHALLOW AND DEEP CUMULUS CONVECTION.
  !          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
  !          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
  !          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KTYPE*        TYPE OF CONVECTION
  !                       1 = PENETRATIVE CONVECTION
  !                       2 = SHALLOW CONVECTION
  !                       3 = MIDLEVEL CONVECTION
  !    *KCBOT*        CLOUD BASE LEVEL
  !    *KDPL*         DEPARTURE LEVEL FOR CONVECTION
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
  !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
  !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
  !    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
  !    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
  !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)      K
  !    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1) KG/KG
  !    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)  KG/KG
  !    *PGEO*         GEOPOTENTIAL                                 M2/S2
  !    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT           KG/KG
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
  !    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS           PA
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
  !    *PVERVEL*      VERTICAL VELOCITY                            PA/S
  !    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)
  !    *PLCRIT_AER*   CRITICAL LIQUID MMR FOR AUTOCONVERSION PROCESS KG/KG
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  !    *LDTDKMF*      Arpege tuning (if TRUE)
  
  !    UPDATED PARAMETERS (INTEGER):
  
  !    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
  !                        KLAB=2 FOR CLOUD LEVELS
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PTU*          TEMPERATURE IN UPDRAFTS                        K
  !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
  !    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
  !    *PLRAIN*       RAIN   WATER CONTENT IN UPDRAFTS             KG/KG
  
  !    OUTPUT PARAMETERS (INTEGER):
  
  !    *KCTOP*        CLOUD TOP LEVEL
  !    *KCTOP0*       FIRST GUESS OF CLOUD TOP LEVEL
  !    *LSCVFLAG*     MASK FOR LIQUID ONLY SHALLOW
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PMFU*         MASSFLUX IN UPDRAFTS                         KG/(M2*S)
  !    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
  !    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS         J/(M2*S)
  !    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS           KG/(M2*S)
  !    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS             KG/(M2*S)
  !    *PLUDE*        DETRAINED TOTAL CONDENSATE                   KG/(M2*S)
  !    *PLUDELI*      DETRAINED LIQUID, ICE                        KG/(M2*S)
  !    *PLGLAC*       FROZEN CLOUD WATER/RAIN CONTENT              KG/KG
  !    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS       KG/(M2*S)
  !    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                     KG/(M2*S)
  !    *PKINEU*       UPDRAFT KINETIC ENERGY                       M2/S2
  !    *PWMEAN*       MEAN UPDRAUGHT VELOCITY                      M/S
  
  !          EXTERNALS
  !          ---------
  !          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
  !          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
  !          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
  
  !          REFERENCE
  !          ---------
  !          (TIEDTKE,1989)
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      01-05-22 : Modified flux limiter M.CULLEN
  !      02-08-14 : Allow for departure level =/ KLEV  P.BECHTOLD
  !      03-08-28 : Clean-up detrainment rates         P.BECHTOLD
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      J.Hague       08-Dec-2005 Tuning: LLFLAG indexing
  !      07-06-01 : Organized entrainment based on RH  P.BECHTOLD
  !      08-02-11 : Simplify entrainment to org no turb P.BECHTOLD
  !      11-03-26 : Detrainment dependent on       RH  P.BECHTOLD
  !      11-09-19 : JJMorcrette effect of prognostic aerosols on autoconversion
  !      16-01-27 : Introduced SPP scheme (LSPP)   M. Leutbecher & S.-J. Lock
  !      30-Jan-20: Single precision fix    F. Vana
  !      20-10-12 : SPP abstraction                M. Leutbecher
  !      21-09-13 : Introduced LDTDKMF for Arpege
  !    2021-09-24 : Fix bug on PKINEU, this bug is active only in simple precision mode.
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUASCN_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOECUMF, ONLY: TECUMF
  USE YOEPHLI, ONLY: TEPHLI
  USE YOECLDP, ONLY: TECLDP
  USE YOM_YGFL, ONLY: TYPE_GFLD
  USE SPP_MOD, ONLY: TSPP_CONFIG
  USE SPP_GEN_MOD, ONLY: SPP_PERT
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TECLDP), INTENT(IN) :: YDECLDP
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  TYPE(TSPP_CONFIG), INTENT(IN) :: YDSPP_CONFIG
  TYPE(TYPE_GFLD), INTENT(IN) :: YGFL
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
  REAL(KIND=JPRB), INTENT(IN) :: PLCRIT_AER(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLITOT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PVERVEL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PWUBASE
  REAL(KIND=JPRB), INTENT(IN) :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
  LOGICAL, INTENT(IN) :: LDLAND
  LOGICAL, INTENT(INOUT) :: LDCUM
  LOGICAL, INTENT(IN) :: LDTDKMF
  LOGICAL, INTENT(OUT) :: LSCVFLAG
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KTYPE
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KLAB(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLRAIN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFUB
  REAL(KIND=JPRB), INTENT(OUT) :: PLGLAC(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLUDE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLUDELI(KLON, KLEV, 2)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFUP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFEN(KLON, KLEV)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KCBOT
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCTOP
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KCTOP0
  INTEGER(KIND=JPIM), INTENT(IN) :: KDPL
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUDE_RATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PKINEU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PWMEAN
  
  REAL(KIND=JPRB) :: ZDMFEN
  REAL(KIND=JPRB) :: ZDMFDE
  REAL(KIND=JPRB) :: ZQOLD
  REAL(KIND=JPRB) :: ZPRECIP
  temp (REAL (KIND=JPRB), ZBUO, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZLUOLD
  REAL(KIND=JPRB) :: ZDPMEAN
  REAL(KIND=JPRB) :: ZOENTR
  REAL(KIND=JPRB) :: ZPH
  LOGICAL :: LLFLAG
  LOGICAL :: LLFLAGUV
  LOGICAL :: LLO1
  LOGICAL :: LLO3
  
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IS
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  INTEGER(KIND=JPIM) :: IKB
  INTEGER(KIND=JPIM) :: JLL
  INTEGER(KIND=JPIM) :: JLM
  INTEGER(KIND=JPIM) :: JLX
  
  REAL(KIND=JPRB) :: Z_CLDMAX
  REAL(KIND=JPRB) :: Z_CPRC2
  REAL(KIND=JPRB) :: Z_CWDRAG
  REAL(KIND=JPRB) :: Z_CWIFRAC
  REAL(KIND=JPRB) :: ZALFAW
  REAL(KIND=JPRB) :: ZBC
  REAL(KIND=JPRB) :: ZBE
  REAL(KIND=JPRB) :: ZBUOC
  REAL(KIND=JPRB) :: ZC
  REAL(KIND=JPRB) :: ZCBF
  REAL(KIND=JPRB) :: ZCONS2
  REAL(KIND=JPRB) :: ZD
  REAL(KIND=JPRB) :: ZDFI
  REAL(KIND=JPRB) :: ZDKBUO
  REAL(KIND=JPRB) :: ZDKEN
  REAL(KIND=JPRB) :: ZDNOPRC
  REAL(KIND=JPRB) :: ZRG
  REAL(KIND=JPRB) :: ZORCPD
  REAL(KIND=JPRB) :: ZDT
  REAL(KIND=JPRB) :: ZFAC
  REAL(KIND=JPRB) :: ZFACBUO
  REAL(KIND=JPRB) :: ZINT
  REAL(KIND=JPRB) :: ZKEDKE
  REAL(KIND=JPRB) :: ZLCRIT
  REAL(KIND=JPRB) :: ZLEEN
  REAL(KIND=JPRB) :: ZLNEW
  REAL(KIND=JPRB) :: ZMFMAX
  REAL(KIND=JPRB) :: ZMFTEST
  REAL(KIND=JPRB) :: ZMFULK
  REAL(KIND=JPRB) :: ZMFUN
  REAL(KIND=JPRB) :: ZMFUQK
  REAL(KIND=JPRB) :: ZMFUSK
  REAL(KIND=JPRB) :: ZOEALFA
  REAL(KIND=JPRB) :: ZOEALFAP
  REAL(KIND=JPRB) :: ZPRCDGW
  REAL(KIND=JPRB) :: ZPRCON
  REAL(KIND=JPRB) :: ZQEEN
  REAL(KIND=JPRB) :: ZQUDE
  REAL(KIND=JPRB) :: ZRNEW
  REAL(KIND=JPRB) :: ZROLD
  REAL(KIND=JPRB) :: ZSCDE
  REAL(KIND=JPRB) :: ZSEEN
  REAL(KIND=JPRB) :: ZVI
  REAL(KIND=JPRB) :: ZVV
  REAL(KIND=JPRB) :: ZVW
  REAL(KIND=JPRB) :: ZWU
  REAL(KIND=JPRB) :: ZZCO
  REAL(KIND=JPRB) :: ZOCUDET
  REAL(KIND=JPRB) :: ZGLAC
  REAL(KIND=JPRB) :: ZXENTRORG
  REAL(KIND=JPRB) :: ZXENTSHALP
  REAL(KIND=JPRB) :: ZXPRCDGW
  
  ! A bunch of SPP variables
  LOGICAL :: LLPERT_ENTRORG  ! SPP perturbation on?
  LOGICAL :: LLPERT_ENTSHALP  ! SPP perturbation on?
  LOGICAL :: LLPERT_RPRCON  ! SPP perturbation on?
  INTEGER(KIND=JPIM) :: IPENTRORG  ! SPP random field pointer
  INTEGER(KIND=JPIM) :: IPENTSHALP  ! SPP random field pointer
  INTEGER(KIND=JPIM) :: IPRPRCON  ! SPP random field pointer
  INTEGER(KIND=JPIM) :: IPN  ! SPP perturbation pointer
  TYPE(SPP_PERT) :: PN1ENTRORG  ! SPP pertn. configs. for ENTRORG, ENTSHALP and RPRCON, respectively
  TYPE(SPP_PERT) :: PN1ENTSHALP  ! SPP pertn. configs. for ENTRORG, ENTSHALP and RPRCON, respectively
  TYPE(SPP_PERT) :: PN1RPRCON  ! SPP pertn. configs. for ENTRORG, ENTSHALP and RPRCON, respectively
  
  REAL(KIND=JPRB) :: ZCHANGE
  REAL(KIND=JPRB) :: ZXS
  REAL(KIND=JPRB) :: ZXE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  LOGICAL :: LLKLAB
  
#include "cuadjtq.intfb.h"
#include "cubasmcn.intfb.h"
#include "cuentr.intfb.h"
  
  !DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ZBUO)
  JL = KIDIA
  !----------------------------------------------------------------------
  
  !*    1.           SPECIFY PARAMETERS
  !                  ------------------
  
  
  ZCONS2 = YDECUMF%RMFCFL / (YDCST%RG*PTSPHY)
  ZRG = 1.0_JPRB / YDCST%RG
  ZORCPD = 1.0_JPRB / YDCST%RCPD
  ZFACBUO = 0.5_JPRB / (1.0_JPRB + 0.5_JPRB)
  ZPRCDGW = YDECUMF%RPRCON / YDCST%RG
  ZDNOPRC = 3.E-4_JPRB    !condensate threshold for precip
  Z_CLDMAX = 5.E-3_JPRB
  Z_CWIFRAC = 0.5_JPRB
  Z_CPRC2 = 0.5_JPRB
  Z_CWDRAG = (3._JPRB / 8._JPRB)*0.506_JPRB / 0.2_JPRB
  
  IF (YDECUMF%LMFGLAC) THEN
    ZGLAC = 1.0_JPRB
  ELSE
    ZGLAC = 0.0_JPRB
  END IF
  
  !----------------------------------------------------------------------
  
  !     2.           SET DEFAULT VALUES
  !                  ------------------
  
  LLO3 = .false.
  ZLUOLD = 0.0_JPRB
  IF (.not.LDCUM) THEN
    KCBOT = -1
    PMFUB = 0.0_JPRB
    PQU(JL, KLEV) = 0.0_JPRB
    KTYPE = 0
  END IF
  PWMEAN = 0.0_JPRB
  ZDPMEAN = 0.0_JPRB
  ZOENTR = 0.0_JPRB
  LSCVFLAG = .false.
  IF (KTYPE > 1 .and. YDECUMF%LSCVLIQ) LSCVFLAG = .true.
  
  !  Prepare SPP
  IF (YDSPP_CONFIG%LSPP) THEN
    
    IPN = YDSPP_CONFIG%PPTR%ENTRORG
    LLPERT_ENTRORG = IPN > 0
    IF (LLPERT_ENTRORG) THEN
      PN1ENTRORG = YDSPP_CONFIG%SM%PN(IPN)
      IPENTRORG = PN1ENTRORG%MP
    END IF
    
    IPN = YDSPP_CONFIG%PPTR%ENTSHALP
    LLPERT_ENTSHALP = IPN > 0
    IF (LLPERT_ENTSHALP) THEN
      PN1ENTSHALP = YDSPP_CONFIG%SM%PN(IPN)
      IPENTSHALP = PN1ENTSHALP%MP
    END IF
    
    IPN = YDSPP_CONFIG%PPTR%RPRCON
    LLPERT_RPRCON = IPN > 0
    IF (LLPERT_RPRCON) THEN
      PN1RPRCON = YDSPP_CONFIG%SM%PN(IPN)
      IPRPRCON = PN1RPRCON%MP
    END IF
    
  ELSE
    LLPERT_ENTRORG = .false.
    LLPERT_ENTSHALP = .false.
    LLPERT_RPRCON = .false.
  END IF
  
  
  ! initalize various quantities
  ! note that liquid water and kinetic energy at cloud base is
  ! preserved from cubase
  
  LLKLAB = .false.
  IF (.not.LDCUM .or. KTYPE == 3) LLKLAB = .true.
  
  DO JK=1,KLEV
    IF (JK /= KCBOT) THEN
      PLU(JL, JK) = 0.0_JPRB
    END IF
    PKINEU(JL, JK) = 0.0_JPRB
    PMFU(JL, JK) = 0.0_JPRB
    PMFUS(JL, JK) = 0.0_JPRB
    PMFUQ(JL, JK) = 0.0_JPRB
    PMFUL(JL, JK) = 0.0_JPRB
    PLUDE(JL, JK) = 0.0_JPRB
    PLUDELI(JL, JK, 1) = 0.0_JPRB
    PLUDELI(JL, JK, 2) = 0.0_JPRB
    PLGLAC(JL, JK) = 0.0_JPRB
    PDMFUP(JL, JK) = 0.0_JPRB
    PLRAIN(JL, JK) = 0.0_JPRB
    ZBUO(JL, JK) = 0.0_JPRB
    IF (LLKLAB) KLAB(JL, JK) = 0
    IF (.not.LDCUM .and. PAPH(JL, JK) < 4.E4_JPRB) KCTOP0 = JK
    PDMFEN(JL, JK) = 0.0_JPRB
    PMFUDE_RATE(JL, JK) = 0.0_JPRB
  END DO
  !DIR$ IVDEP
  !OCL NOVREC
  IF (KTYPE == 3) LDCUM = .false.
  
  !----------------------------------------------------------------------
  
  !     3.0          INITIALIZE VALUES AT cloud base LEVEL
  !                  -------------------------------------
  
  KCTOP = KCBOT
  IF (LDCUM) THEN
    IKB = KCBOT
    PKINEU(JL, IKB) = 0.5*PWUBASE**2
    PMFU(JL, IKB) = PMFUB
    PMFUS(JL, IKB) = PMFUB*(YDCST%RCPD*PTU(JL, IKB) + PGEOH(JL, IKB))
    PMFUQ(JL, IKB) = PMFUB*PQU(JL, IKB)
    PMFUL(JL, IKB) = PMFUB*PLU(JL, IKB)
  END IF
  
  !----------------------------------------------------------------------
  
  !     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
  !                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
  !                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
  !                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
  !                  -------------------------------------------------
  
  DO JK=KLEV - 1,3,-1
    
    !                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
    !                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
    !                  ----------------------------------------------------
    
    IK = JK
    CALL CUBASMCN_OPENACC(YDCST, YDECUMF, KIDIA, KFDIA, KLON, KLEV, IK, PTEN, PQEN, PQSEN, PVERVEL, PGEO, PGEOH, LDCUM, KTYPE,  &
    & KLAB, KCBOT, PMFU, PMFUB, PLRAIN, PTU, PQU, PLU, PMFUS, PMFUQ, PMFUL, PDMFUP, YDSTACK=YLSTACK)
    
    IS = 0
    JLM = 0
    ! also liquid only for ktype=3
    ! IF(KTYPE(JL)>1.AND.LSCVLIQ) LSCVFLAG(JL)=.TRUE.
    LLFLAG = .false.
    ZPRECIP = 0.0_JPRB
    LLO1 = .false.
    IS = IS + KLAB(JL, JK + 1)
    IF (KLAB(JL, JK + 1) == 0) KLAB(JL, JK) = 0
    IF (LDCUM .and. KLAB(JL, JK + 1) == 2 .or. KTYPE == 3 .and. KLAB(JL, JK + 1) == 1) THEN
      LLFLAG = .true.
      JLM = JLM + 1
      JLX = JL
    END IF
    IF (KLAB(JL, JK + 1) > 0) THEN
      LLFLAGUV = .true.
    ELSE
      LLFLAGUV = .false.
    END IF
    ZPH = PAPH(JL, JK)
    IF (KTYPE == 3 .and. JK == KCBOT) THEN
      ZMFMAX = (PAPH(JL, JK) - PAPH(JL, JK - 1))*ZCONS2
      IF (PMFUB > ZMFMAX) THEN
        ZFAC = ZMFMAX / PMFUB
        PMFU(JL, JK + 1) = PMFU(JL, JK + 1)*ZFAC
        PMFUS(JL, JK + 1) = PMFUS(JL, JK + 1)*ZFAC
        PMFUQ(JL, JK + 1) = PMFUQ(JL, JK + 1)*ZFAC
        PMFUB = ZMFMAX
      END IF
    END IF
    
    IF (IS > 0) LLO3 = .true.
    
    !*                  SPECIFY ENTRAINMENT RATES IN *CUENTR*
    !                   -------------------------------------
    
    IK = JK
    CALL CUENTR_OPENACC(YDCST, YDECUMF, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, IK, KCBOT, KTYPE, LDCUM, LLO3, PQSEN, PAPH,  &
    & PGEOH, PMFU, PGP2DSPP, ZDMFEN, ZDMFDE, YDSTACK=YLSTACK)
    
    !                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
    !                  ---------------------------------------------------
    
    IF (LLO3) THEN
      
      ZQOLD = 0.0_JPRB
      DO JLL=1,JLM
        JL = JLX
        ZDMFDE = MIN(ZDMFDE, 0.75_JPRB*PMFU(JL, JK + 1))
        IF (JK == KCBOT) THEN
          IF (YDSPP_CONFIG%LSPP .and. LLPERT_ENTRORG) THEN
            ZXENTRORG = YDECUMF%ENTRORG*EXP(PN1ENTRORG%MU(1) + PN1ENTRORG%XMAG(1)*PGP2DSPP(JL, IPENTRORG))
          ELSE
            ZXENTRORG = YDECUMF%ENTRORG
          END IF
          ZOENTR =  &
          & -ZXENTRORG*(MIN(1.0_JPRB, PQEN(JL, JK) / PQSEN(JL, JK)) - YDECUMF%ENTR_RH)*(PGEOH(JL, JK) - PGEOH(JL, JK + 1))*ZRG
          ZOENTR = MIN(0.4_JPRB, ZOENTR)*PMFU(JL, JK + 1)
        END IF
        IF (JK < KCBOT) THEN
          ZMFMAX = (PAPH(JL, JK) - PAPH(JL, JK - 1))*ZCONS2
          ZXS = MAX(PMFU(JL, JK + 1) - ZMFMAX, 0.0_JPRB)
          PWMEAN = PWMEAN + PKINEU(JL, JK + 1)*(PAP(JL, JK + 1) - PAP(JL, JK))
          ZDPMEAN = ZDPMEAN + PAP(JL, JK + 1) - PAP(JL, JK)
          ZDMFEN = ZOENTR
          IF (KTYPE >= 2) THEN
            IF (YDSPP_CONFIG%LSPP .and. LLPERT_ENTSHALP) THEN
              ZXENTSHALP = YDECUMF%ENTSHALP*EXP(PN1ENTSHALP%MU(1) + PN1ENTSHALP%XMAG(1)*PGP2DSPP(JL, IPENTSHALP))
            ELSE
              ZXENTSHALP = YDECUMF%ENTSHALP
            END IF
            ZDMFEN = ZXENTSHALP*ZDMFEN
            ZDMFDE = ZDMFEN
            ! ZDMFDE(JL)=MAX(ZDMFDE(JL),ZDMFEN(JL))
          END IF
          ZDMFDE = ZDMFDE*(1.6_JPRB - MIN(1.0_JPRB, PQEN(JL, JK) / PQSEN(JL, JK)))
          ZMFTEST = PMFU(JL, JK + 1) + ZDMFEN - ZDMFDE
          ZCHANGE = MAX(ZMFTEST - ZMFMAX, 0.0_JPRB)
          ZXE = MAX(ZCHANGE - ZXS, 0.0_JPRB)
          ZDMFEN = ZDMFEN - ZXE
          ZCHANGE = ZCHANGE - ZXE
          ZDMFDE = ZDMFDE + ZCHANGE
        END IF
        
        PDMFEN(JL, JK) = ZDMFEN - ZDMFDE
        
        PMFU(JL, JK) = PMFU(JL, JK + 1) + ZDMFEN - ZDMFDE
        ZQEEN = PQENH(JL, JK + 1)*ZDMFEN
        ZSEEN = (YDCST%RCPD*PTENH(JL, JK + 1) + PGEOH(JL, JK + 1))*ZDMFEN
        IF (PLITOT(JL, JK) > YDECLDP%RLMIN) THEN
          ZLEEN = PLITOT(JL, JK)*ZDMFEN
        ELSE
          ZLEEN = 0.0_JPRB
        END IF
        ZSCDE = (YDCST%RCPD*PTU(JL, JK + 1) + PGEOH(JL, JK + 1))*ZDMFDE
        ZQUDE = PQU(JL, JK + 1)*ZDMFDE
        PLUDE(JL, JK) = PLU(JL, JK + 1)*ZDMFDE
        ZMFUSK = PMFUS(JL, JK + 1) + ZSEEN - ZSCDE
        ZMFUQK = PMFUQ(JL, JK + 1) + ZQEEN - ZQUDE
        ZMFULK = PMFUL(JL, JK + 1) + ZLEEN - PLUDE(JL, JK)
        ZFAC = 1.0_JPRB / MAX(YDECUMF%RMFCMIN, PMFU(JL, JK))
        PLU(JL, JK) = ZMFULK*ZFAC
        PQU(JL, JK) = ZMFUQK*ZFAC
        PTU(JL, JK) = (ZMFUSK*ZFAC - PGEOH(JL, JK))*ZORCPD
        PTU(JL, JK) = MAX(100._JPRB, PTU(JL, JK))
        PTU(JL, JK) = MIN(400._JPRB, PTU(JL, JK))
        ZQOLD = PQU(JL, JK)
        PLRAIN(JL, JK) = PLRAIN(JL, JK + 1)*MAX(0.0_JPRB, PMFU(JL, JK + 1) - ZDMFDE)*ZFAC
        ZLUOLD = PLU(JL, JK)
      END DO
      ! reset to environmental values if below departure level
      IF (JK > KDPL) THEN
        PTU(JL, JK) = PTENH(JL, JK)
        PQU(JL, JK) = PQENH(JL, JK)
        PLU(JL, JK) = 0.0_JPRB
        ZLUOLD = PLU(JL, JK)
      END IF
      
      !                  DO CORRECTIONS FOR MOIST ASCENT
      !                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
      !                  -----------------------------------
      
      IK = JK
      IF (JLM > 0) THEN
        IF (YDECUMF%LSCVLIQ) THEN
          CALL CUADJTQ_OPENACC(YDTHF, YDCST, YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK, ZPH, PTU, PQU, LLFLAG, 6, LSCVFLAG,  &
          & YDSTACK=YLSTACK)
        ELSE
          CALL CUADJTQ_OPENACC(YDTHF, YDCST, YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK, ZPH, PTU, PQU, LLFLAG, 1, LSCVFLAG,  &
          & YDSTACK=YLSTACK)
        END IF
      END IF
      
      IF (YDEPHLI%LPHYLIN) THEN
        
        !DIR$ IVDEP
        !OCL NOVREC
        DO JLL=1,JLM
          JL = JLX
          IF (PQU(JL, JK) /= ZQOLD) THEN
            ZOEALFA = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(PTU(JL, JK) - YDEPHLI%RLPTRC)) + 1.0_JPRB))
            ZOEALFAP = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(PTU(JL, JK + 1) - YDEPHLI%RLPTRC)) + 1.0_JPRB))
            PLGLAC(JL, JK) = PLU(JL, JK)*((1.0_JPRB - ZOEALFA) - (1.0_JPRB - ZOEALFAP))
            ! add glaciation of rain
            ZFAC = 0.545_JPRB*(TANH(0.17_JPRB*(PTEN(JL, JK) - YDEPHLI%RLPTRC)) + 1.0_JPRB)
            PLGLAC(JL, JK) = PLGLAC(JL, JK) + ZFAC*PDMFUP(JL, JK + 1) / MAX(YDECUMF%RMFCMIN, PMFU(JL, JK + 1))*(0.5_JPRB +  &
            & SIGN(0.5_JPRB, YDCST%RTT - PTEN(JL, JK)))*ZGLAC
            PTU(JL, JK) = PTU(JL, JK) + YDTHF%RALFDCP*PLGLAC(JL, JK)
          END IF
        END DO
        
      ELSE
        
        !DIR$ IVDEP
        !OCL NOVREC
        DO JLL=1,JLM
          JL = JLX
          IF (PQU(JL, JK) /= ZQOLD) THEN
            PLGLAC(JL, JK) = PLU(JL, JK)*((1.0_JPRB - FOEALFCU(PTU(JL, JK))) - (1.0_JPRB - FOEALFCU(PTU(JL, JK + 1))))
            IF (LSCVFLAG) PLGLAC(JL, JK) = 0.0_JPRB
            ! add glaciation of rain, only fraction added to updraught heat
            ZFAC = FOEALFCU(PTEN(JL, JK))
            PLGLAC(JL, JK) = PLGLAC(JL, JK) + ZFAC*PDMFUP(JL, JK + 1) / MAX(YDECUMF%RMFCMIN, PMFU(JL, JK + 1))*(0.5_JPRB +  &
            & SIGN(0.5_JPRB, YDCST%RTT - PTEN(JL, JK)))*ZGLAC
            PTU(JL, JK) = PTU(JL, JK) + YDTHF%RALFDCP*PLGLAC(JL, JK)
          END IF
        END DO
        
      END IF
      
      DO JLL=1,JLM
        JL = JLX
        IF (PQU(JL, JK) /= ZQOLD) THEN
          KLAB(JL, JK) = 2
          PLU(JL, JK) = PLU(JL, JK) + ZQOLD - PQU(JL, JK)
          ZBC = PTU(JL, JK)*(1.0_JPRB + YDCST%RETV*PQU(JL, JK) - PLU(JL, JK + 1) - PLRAIN(JL, JK + 1))
          ZBE = PTENH(JL, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JL, JK))
          ZBUO(JL, JK) = ZBC - ZBE
          
          ! set flags in case of midlevel convection
          
          IF (KTYPE == 3 .and. KLAB(JL, JK + 1) == 1) THEN
            IF (ZBUO(JL, JK) > -0.5_JPRB) THEN
              LDCUM = .true.
              KCTOP = JK
              PKINEU(JL, JK) = 0.5_JPRB
            ELSE
              KLAB(JL, JK) = 0
              PMFU(JL, JK) = 0.0_JPRB
              PLUDE(JL, JK) = 0.0_JPRB
              PLU(JL, JK) = 0.0_JPRB
            END IF
          END IF
          
          IF (KLAB(JL, JK + 1) == 2) THEN
            
            !IF(ZBUO(JL,JK) < 0.0_JPRB.AND.KTYPE(JL)==1) THEN
            IF (ZBUO(JL, JK) < -0.1_JPRB) THEN
              PTENH(JL, JK) = 0.5_JPRB*(PTEN(JL, JK) + PTEN(JL, JK - 1))
              PQENH(JL, JK) = 0.5_JPRB*(PQEN(JL, JK) + PQEN(JL, JK - 1))
              ZBUO(JL, JK) = ZBC - PTENH(JL, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JL, JK))
            END IF
            ZBUOC = 0.5*(ZBUO(JL, JK) + ZBUO(JL, JK + 1)) / (PTENH(JL, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JL, JK)))
            ZDKBUO = (PGEOH(JL, JK) - PGEOH(JL, JK + 1))*ZFACBUO*ZBUOC
            
            ! either use entrainment rate or if zero
            ! use detrainmnet rate as a substitute for
            ! mixing and "pressure" gradient term in upper
            ! troposphere
            
            IF (ZDMFEN > 0.0_JPRB) THEN
              ZDKEN = MIN(1.0_JPRB, (1 + Z_CWDRAG)*ZDMFEN / MAX(YDECUMF%RMFCMIN, PMFU(JL, JK + 1)))
            ELSE
              ZDKEN = MIN(1.0_JPRB, (1 + Z_CWDRAG)*ZDMFDE / MAX(YDECUMF%RMFCMIN, PMFU(JL, JK + 1)))
            END IF
            
            IF (LDTDKMF) THEN
              PKINEU(JL, JK) = (PKINEU(JL, JK + 1)*(1 - ZDKEN) + ZDKBUO) / (1 + ZDKEN)
              ZBUOC = ZBUO(JL, JK)
            ELSE
              PKINEU(JL, JK) = MAX(-1.E3_JPRB, (PKINEU(JL, JK + 1)*(1 - ZDKEN) + ZDKBUO) / (1 + ZDKEN))
            END IF
            IF (ZBUOC < 0.0_JPRB) THEN
              ZKEDKE = PKINEU(JL, JK) / MAX(1.E-3_JPRB, PKINEU(JL, JK + 1))
              ZKEDKE = MAX(0.0_JPRB, MIN(1.0_JPRB, ZKEDKE))
              ZOCUDET = (1.6_JPRB - MIN(1.0_JPRB, PQEN(JL, JK) / PQSEN(JL, JK)))
              ZMFUN = ZOCUDET*SQRT(ZKEDKE)
              ZDMFDE = MAX(ZDMFDE, PMFU(JL, JK + 1)*(1.0_JPRB - ZMFUN))
              PLUDE(JL, JK) = PLU(JL, JK + 1)*ZDMFDE
              PMFU(JL, JK) = PMFU(JL, JK + 1) + ZDMFEN - ZDMFDE
            END IF
            
            IF (ZBUO(JL, JK) > -0.2_JPRB) THEN
              IKB = KCBOT
              IF (YDSPP_CONFIG%LSPP .and. LLPERT_ENTRORG) THEN
                ZXENTRORG = YDECUMF%ENTRORG*EXP(PN1ENTRORG%MU(1) + PN1ENTRORG%XMAG(1)*PGP2DSPP(JL, IPENTRORG))
              ELSE
                ZXENTRORG = YDECUMF%ENTRORG
              END IF
              ZOENTR = ZXENTRORG*(0.3_JPRB - (MIN(1.0_JPRB, PQEN(JL, JK - 1) / PQSEN(JL, JK - 1)) - YDECUMF%ENTR_RH))*(PGEOH(JL,  &
              & JK - 1) - PGEOH(JL, JK))*ZRG*MIN(1.0_JPRB, PQSEN(JL, JK) / PQSEN(JL, IKB))**3
              ZOENTR = MIN(0.4_JPRB, ZOENTR)*PMFU(JL, JK)
            ELSE
              ZOENTR = 0.0_JPRB
            END IF
            
            ! Erase values if below departure level
            IF (JK > KDPL) THEN
              PMFU(JL, JK) = PMFU(JL, JK + 1)
              PKINEU(JL, JK) = 0.5_JPRB
            END IF
            IF (PKINEU(JL, JK) > 0.0_JPRB .and. PMFU(JL, JK) > 0.0_JPRB .and. (ZBUO(JL, JK) > -2._JPRB .or. (PTEN(JL, JK - 1) -  &
            & PTEN(JL, JK)) / (ZRG*(PGEO(JL, JK - 1) - PGEO(JL, JK))) < -3.E-3_JPRB)) THEN
              ! add overshoot limiter for convective top evaluation
              KCTOP = JK
              LLO1 = .true.
            ELSE
              KLAB(JL, JK) = 0
              PMFU(JL, JK) = 0.0_JPRB
              PKINEU(JL, JK) = 0.0_JPRB
              ZDMFDE = PMFU(JL, JK + 1)
              PLUDE(JL, JK) = PLU(JL, JK + 1)*ZDMFDE
            END IF
            
            ! store detrainment rates for updraught
            
            IF (PMFU(JL, JK + 1) > 0.0_JPRB) THEN
              PMFUDE_RATE(JL, JK) = ZDMFDE
            END IF
            
          END IF
          
          !     ELSEIF(LLFLAG(JL).AND.KTYPE(JL)==2.AND.PQU(JL,JK) == ZQOLD(JL)) THEN
        ELSE IF (KTYPE == 2 .and. PQU(JL, JK) == ZQOLD) THEN
          KLAB(JL, JK) = 0
          PMFU(JL, JK) = 0.0_JPRB
          PKINEU(JL, JK) = 0.0_JPRB
          ZDMFDE = PMFU(JL, JK + 1)
          PLUDE(JL, JK) = PLU(JL, JK + 1)*ZDMFDE
          PMFUDE_RATE(JL, JK) = ZDMFDE
          
        END IF
      END DO
      
      !              CALCULATE PRECIPITATION RATE BY
      !              ANALYTIC INTEGRATION OF EQUATION FOR L
      
      IF (LLO1) THEN
        IF (PLU(JL, JK) > ZDNOPRC) THEN
          ZWU = MIN(15._JPRB, SQRT(2.0_JPRB*MAX(0.5_JPRB, PKINEU(JL, JK + 1))))
          !       increase conversion for liquid phase only
          IF (LDTDKMF) THEN
            ZZCO = FOEALFCU(PTU(JL, JK))*1.3_JPRB + (1.0_JPRB - FOEALFCU(PTU(JL, JK)))
          ELSE
            ZZCO = 1.0_JPRB + 0.3_JPRB*FOEALFCU(PTU(JL, JK))
          END IF
          ! IF(LSCVFLAG(JL)) ZZCO=1.3_JPRB
          
          IF (YDSPP_CONFIG%LSPP .and. LLPERT_RPRCON) THEN
            ZXPRCDGW = ZPRCDGW*EXP(PN1RPRCON%MU(1) + PN1RPRCON%XMAG(1)*PGP2DSPP(JL, IPRPRCON))
          ELSE
            ZXPRCDGW = ZPRCDGW
          END IF
          ZPRCON = ZXPRCDGW / (0.75_JPRB*ZWU)*ZZCO
          
          !           PARAMETERS FOR BERGERON-FINDEISEN PROCESS (T < -5C)
          
          ZDT = MIN(YDTHF%RTBERCU - YDTHF%RTICECU, MAX(YDTHF%RTBERCU - PTU(JL, JK), 0.0_JPRB))
          ZCBF = 1 + Z_CPRC2*SQRT(ZDT)
          ZZCO = ZPRCON*ZCBF
          
          IF (YGFL%NACTAERO > 0) THEN
            !         (Prognostic) aerosols in clouds, CCN effect only in pure warm phase
            !         then phases out in mixed phase - ***NOT*** linear physics
            ! CCN number taken at each height from environment
            IF (YDECLDP%LAERLIQAUTOCP) THEN
              ZALFAW = FOEALFCU(PTU(JL, JK))
              ZLCRIT = ZALFAW*PLCRIT_AER(JL, JK) + (1.0_JPRB - ZALFAW)*ZDNOPRC
              ! CCN number determined at cloud base (entrainment ignored)
            ELSE IF (YDECLDP%LAERLIQAUTOCPB) THEN
              ZALFAW = FOEALFCU(PTU(JL, JK))
              ZLCRIT = ZALFAW*PLCRIT_AER(JL, KCBOT) + (1.0_JPRB - ZALFAW)*ZDNOPRC
            ELSE
              ZLCRIT = ZDNOPRC / ZCBF
            END IF
          ELSE
            ZLCRIT = ZDNOPRC / ZCBF
          END IF
          
          ZDFI = PGEOH(JL, JK) - PGEOH(JL, JK + 1)
          ZC = (PLU(JL, JK) - ZLUOLD)
          ZD = ZZCO*(1.0_JPRB - EXP(-(PLU(JL, JK) / ZLCRIT)**2))*ZDFI
          ZINT = EXP(-ZD)
          ZLNEW = ZLUOLD*ZINT + ZC / ZD*(1.0_JPRB - ZINT)
          ZLNEW = MAX(0.0_JPRB, MIN(PLU(JL, JK), ZLNEW))
          ZLNEW = MIN(Z_CLDMAX, ZLNEW)
          ZPRECIP = MAX(0.0_JPRB, ZLUOLD + ZC - ZLNEW)
          PDMFUP(JL, JK) = ZPRECIP*PMFU(JL, JK)
          PLRAIN(JL, JK) = PLRAIN(JL, JK) + ZPRECIP
          PLU(JL, JK) = ZLNEW
        END IF
      END IF
      
      IF (YDEPHLI%LPHYLIN) THEN
        
        !DEC$ IVDEP
        IF (LLO1) THEN
          IF (PLRAIN(JL, JK) > 0.0_JPRB) THEN
            ZVW = 21.18_JPRB*PLRAIN(JL, JK)**0.2_JPRB
            ZVI = Z_CWIFRAC*ZVW
            ZALFAW = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(PTU(JL, JK) - YDEPHLI%RLPTRC)) + 1.0_JPRB))
            ZVV = ZALFAW*ZVW + (1.0_JPRB - ZALFAW)*ZVI
            ZROLD = PLRAIN(JL, JK) - ZPRECIP
            ZC = ZPRECIP
            ZWU = MIN(15._JPRB, SQRT(2.0_JPRB*MAX(0.1_JPRB, PKINEU(JL, JK))))
            ZD = ZVV / ZWU
            ZINT = EXP(-ZD)
            ZRNEW = ZROLD*ZINT + ZC / ZD*(1.0_JPRB - ZINT)
            ZRNEW = MAX(0.0_JPRB, MIN(PLRAIN(JL, JK), ZRNEW))
            PLRAIN(JL, JK) = ZRNEW
          END IF
        END IF
        
      ELSE
        
        IF (LLO1) THEN
          IF (PLRAIN(JL, JK) > 0.0_JPRB) THEN
            ZVW = 21.18_JPRB*PLRAIN(JL, JK)**0.2_JPRB
            ZVI = Z_CWIFRAC*ZVW
            ZALFAW = FOEALFCU(PTU(JL, JK))
            IF (LSCVFLAG) ZALFAW = 1.0_JPRB
            ZVV = ZALFAW*ZVW + (1.0_JPRB - ZALFAW)*ZVI
            ZROLD = PLRAIN(JL, JK) - ZPRECIP
            ZC = ZPRECIP
            ZWU = MIN(15._JPRB, SQRT(2.0_JPRB*MAX(0.1_JPRB, PKINEU(JL, JK))))
            ZD = ZVV / ZWU
            ZINT = EXP(-ZD)
            ZRNEW = ZROLD*ZINT + ZC / ZD*(1.0_JPRB - ZINT)
            ZRNEW = MAX(0.0_JPRB, MIN(PLRAIN(JL, JK), ZRNEW))
            PLRAIN(JL, JK) = ZRNEW
          END IF
        END IF
        
      END IF
      
      DO JLL=1,JLM
        JL = JLX
        PMFUL(JL, JK) = PLU(JL, JK)*PMFU(JL, JK)
        PMFUS(JL, JK) = (YDCST%RCPD*PTU(JL, JK) + PGEOH(JL, JK))*PMFU(JL, JK)
        PMFUQ(JL, JK) = PQU(JL, JK)*PMFU(JL, JK)
        ZALFAW = FOEALFCU(PTU(JL, JK))
        IF (LSCVFLAG) ZALFAW = 1.0_JPRB
        PLUDELI(JL, JK, 1) = ZALFAW*PLUDE(JL, JK)
        PLUDELI(JL, JK, 2) = (1.0_JPRB - ZALFAW)*PLUDE(JL, JK)
      END DO
      
    END IF
  END DO
  
  !----------------------------------------------------------------------
  
  !     5.           FINAL CALCULATIONS
  !                  ------------------
  
  IF (KCTOP == -1) LDCUM = .false.
  KCBOT = MAX(KCBOT, KCTOP)
  IF (LDCUM) THEN
    PWMEAN = MAX(1.E-2_JPRB, PWMEAN / MAX(1.0_JPRB, ZDPMEAN))
    PWMEAN = SQRT(2.0_JPRB*PWMEAN)
  END IF
  
END SUBROUTINE CUASCN_OPENACC
