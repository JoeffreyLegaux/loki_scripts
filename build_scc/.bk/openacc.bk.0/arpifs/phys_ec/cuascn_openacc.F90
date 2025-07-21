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
  REAL(KIND=JPRB), INTENT(IN) :: PWUBASE(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
  LOGICAL, INTENT(IN) :: LDLAND(KLON)
  LOGICAL, INTENT(INOUT) :: LDCUM(KLON)
  LOGICAL, INTENT(IN) :: LDTDKMF
  LOGICAL, INTENT(OUT) :: LSCVFLAG(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KTYPE(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KLAB(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLRAIN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFUB(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PLGLAC(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLUDE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLUDELI(KLON, KLEV, 2)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFUP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFEN(KLON, KLEV)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCTOP(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KCTOP0(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KDPL(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFUDE_RATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PKINEU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PWMEAN(KLON)
  
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
#include "abor1.intfb.h"
#include "cuadjtq.func.h"
  LOGICAL :: CUENTR_LLO1
  INTEGER(KIND=JPIM) :: CUENTR_JL
  LOGICAL :: LLPERT_DETRPEN  ! SPP perturbation on?
  INTEGER(KIND=JPIM) :: IPDETRPEN  ! SPP random field pointer
  INTEGER(KIND=JPIM) :: CUENTR_IPN  ! SPP perturbation pointer
  TYPE(SPP_PERT) :: PN1  ! SPP pertn. config. for RTAU
  REAL(KIND=JPRB) :: ZDZ
  REAL(KIND=JPRB) :: ZENTR
  REAL(KIND=JPRB) :: ZMF
  REAL(KIND=JPRB) :: CUENTR_ZRG
  REAL(KIND=JPRB) :: ZXDETRPEN
  REAL(KIND=JPHOOK) :: CUENTR_ZHOOK_HANDLE
  INTEGER(KIND=JPIM) :: CUBASMCN_JL
  REAL(KIND=JPRB) :: ZZZMB
  REAL(KIND=JPRB) :: CUBASMCN_ZRG
  REAL(KIND=JPRB) :: CUBASMCN_ZORCPD
  REAL(KIND=JPHOOK) :: CUBASMCN_ZHOOK_HANDLE
  INTEGER(KIND=JPIM) :: CUADJTQ_JL
  REAL(KIND=JPRB) :: Z1S
  REAL(KIND=JPRB) :: Z2S
  REAL(KIND=JPRB) :: ZCOND
  REAL(KIND=JPRB) :: ZCOND1
  REAL(KIND=JPRB) :: ZCOR
  REAL(KIND=JPRB) :: ZFOEEWI
  REAL(KIND=JPRB) :: ZFOEEWL
  REAL(KIND=JPRB) :: CUADJTQ_ZOEALFA
  REAL(KIND=JPRB) :: ZQMAX
  REAL(KIND=JPRB) :: ZQSAT
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZQP
  REAL(KIND=JPRB) :: ZL
  REAL(KIND=JPRB) :: ZI
  REAL(KIND=JPRB) :: ZF
  LOGICAL :: CUADJTQ_LLFLAG
  REAL(KIND=JPHOOK) :: CUADJTQ_ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ZBUO)
  JLON = KIDIA
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
  IF (.not.LDCUM(JLON)) THEN
    KCBOT(JLON) = -1
    PMFUB(JLON) = 0.0_JPRB
    PQU(JLON, KLEV) = 0.0_JPRB
    KTYPE(JLON) = 0
  END IF
  PWMEAN(JLON) = 0.0_JPRB
  ZDPMEAN = 0.0_JPRB
  ZOENTR = 0.0_JPRB
  LSCVFLAG(JLON) = .false.
  IF (KTYPE(JLON) > 1 .and. YDECUMF%LSCVLIQ) LSCVFLAG(JLON) = .true.
  
  !  Prepare SPP
  LLPERT_ENTRORG = .false.
  LLPERT_ENTSHALP = .false.
  LLPERT_RPRCON = .false.
  
  
  ! initalize various quantities
  ! note that liquid water and kinetic energy at cloud base is
  ! preserved from cubase
  
  LLKLAB = .false.
  IF (.not.LDCUM(JLON) .or. KTYPE(JLON) == 3) LLKLAB = .true.
  
  DO JK=1,KLEV
    IF (JK /= KCBOT(JLON)) THEN
      PLU(JLON, JK) = 0.0_JPRB
    END IF
    PKINEU(JLON, JK) = 0.0_JPRB
    PMFU(JLON, JK) = 0.0_JPRB
    PMFUS(JLON, JK) = 0.0_JPRB
    PMFUQ(JLON, JK) = 0.0_JPRB
    PMFUL(JLON, JK) = 0.0_JPRB
    PLUDE(JLON, JK) = 0.0_JPRB
    PLUDELI(JLON, JK, 1) = 0.0_JPRB
    PLUDELI(JLON, JK, 2) = 0.0_JPRB
    PLGLAC(JLON, JK) = 0.0_JPRB
    PDMFUP(JLON, JK) = 0.0_JPRB
    PLRAIN(JLON, JK) = 0.0_JPRB
    ZBUO(JLON, JK) = 0.0_JPRB
    IF (LLKLAB) KLAB(JLON, JK) = 0
    IF (.not.LDCUM(JLON) .and. PAPH(JLON, JK) < 4.E4_JPRB) KCTOP0(JLON) = JK
    PDMFEN(JLON, JK) = 0.0_JPRB
    PMFUDE_RATE(JLON, JK) = 0.0_JPRB
  END DO
  !DIR$ IVDEP
  !OCL NOVREC
  IF (KTYPE(JLON) == 3) LDCUM(JLON) = .false.
  
  !----------------------------------------------------------------------
  
  !     3.0          INITIALIZE VALUES AT cloud base LEVEL
  !                  -------------------------------------
  
  KCTOP(JLON) = KCBOT(JLON)
  IF (LDCUM(JLON)) THEN
    IKB = KCBOT(JLON)
    PKINEU(JLON, IKB) = 0.5*PWUBASE(JLON)**2
    PMFU(JLON, IKB) = PMFUB(JLON)
    PMFUS(JLON, IKB) = PMFUB(JLON)*(YDCST%RCPD*PTU(JLON, IKB) + PGEOH(JLON, IKB))
    PMFUQ(JLON, IKB) = PMFUB(JLON)*PQU(JLON, IKB)
    PMFUL(JLON, IKB) = PMFUB(JLON)*PLU(JLON, IKB)
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
    ! [Loki] inlined child subroutine: CUBASMCN
    ! =========================================
    
    !----------------------------------------------------------------------
    
    
    !*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    !                  -------------------------------------------
    
    !DIR$ IVDEP
    !OCL NOVREC
    
    CUBASMCN_ZRG = 1.0_JPRB / YDCST%RG
    CUBASMCN_ZORCPD = 1.0_JPRB / YDCST%RCPD
    IF (.not.LDCUM(JLON) .and. KLAB(JLON, IK + 1) == 0) THEN
      IF (YDECUMF%LMFMID .and. IK > YDECUMF%NJKT7 .and. PQEN(JLON, IK) > 0.80_JPRB*PQSEN(JLON, IK)) THEN
        PTU(JLON, IK + 1) = (YDCST%RCPD*PTEN(JLON, IK) + PGEO(JLON, IK) - PGEOH(JLON, IK + 1))*CUBASMCN_ZORCPD
        PQU(JLON, IK + 1) = PQEN(JLON, IK)
        PLU(JLON, IK + 1) = 0.0_JPRB
        ! ZZZMB=MAX(RMFCMIN,-PVERVEL(JL,KK)/RG)
        ZZZMB = MAX(1.E-6_JPRB, -PVERVEL(JLON, IK)*CUBASMCN_ZRG)
        ZZZMB = MIN(ZZZMB, YDECUMF%RMFLIA)
        PMFUB(JLON) = ZZZMB
        PMFU(JLON, IK + 1) = PMFUB(JLON)
        PMFUS(JLON, IK + 1) = PMFUB(JLON)*(YDCST%RCPD*PTU(JLON, IK + 1) + PGEOH(JLON, IK + 1))
        PMFUQ(JLON, IK + 1) = PMFUB(JLON)*PQU(JLON, IK + 1)
        PMFUL(JLON, IK + 1) = 0.0_JPRB
        PDMFUP(JLON, IK + 1) = 0.0_JPRB
        PLRAIN(JLON, IK + 1) = 0.0_JPRB
        KCBOT(JLON) = IK
        KLAB(JLON, IK + 1) = 1
        KTYPE(JLON) = 3
      END IF
    END IF
    
    ! =========================================
    
    IS = 0
    JLM = 0
    ! also liquid only for ktype=3
    ! IF(KTYPE(JL)>1.AND.LSCVLIQ) LSCVFLAG(JL)=.TRUE.
    LLFLAG = .false.
    ZPRECIP = 0.0_JPRB
    LLO1 = .false.
    IS = IS + KLAB(JLON, JK + 1)
    IF (KLAB(JLON, JK + 1) == 0) KLAB(JLON, JK) = 0
    IF (LDCUM(JLON) .and. KLAB(JLON, JK + 1) == 2 .or. KTYPE(JLON) == 3 .and. KLAB(JLON, JK + 1) == 1) THEN
      LLFLAG = .true.
      JLM = JLM + 1
      JLX = JLON
    END IF
    IF (KLAB(JLON, JK + 1) > 0) THEN
      LLFLAGUV = .true.
    ELSE
      LLFLAGUV = .false.
    END IF
    ZPH = PAPH(JLON, JK)
    IF (KTYPE(JLON) == 3 .and. JK == KCBOT(JLON)) THEN
      ZMFMAX = (PAPH(JLON, JK) - PAPH(JLON, JK - 1))*ZCONS2
      IF (PMFUB(JLON) > ZMFMAX) THEN
        ZFAC = ZMFMAX / PMFUB(JLON)
        PMFU(JLON, JK + 1) = PMFU(JLON, JK + 1)*ZFAC
        PMFUS(JLON, JK + 1) = PMFUS(JLON, JK + 1)*ZFAC
        PMFUQ(JLON, JK + 1) = PMFUQ(JLON, JK + 1)*ZFAC
        PMFUB(JLON) = ZMFMAX
      END IF
    END IF
    
    IF (IS > 0) LLO3 = .true.
    
    !*                  SPECIFY ENTRAINMENT RATES IN *CUENTR*
    !                   -------------------------------------
    
    IK = JK
    ! [Loki] inlined child subroutine: CUENTR
    ! =========================================
    
    !----------------------------------------------------------------------
    
    !*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    !                  -------------------------------------------
    
    IF (LLO3) THEN
      
      CUENTR_ZRG = 1.0_JPRB / YDCST%RG
      
      ! prepare SPP perturbations (just once)
      LLPERT_DETRPEN = .false.
      
      ZDMFEN = 0.0_JPRB
      ZDMFDE = 0.0_JPRB
      ZENTR = 0.0
      
      !*    1.1          SPECIFY ENTRAINMENT RATES
      !                  -------------------------
      
      IF (LDCUM(JLON)) THEN
        ZDZ = (PGEOH(JLON, IK) - PGEOH(JLON, IK + 1))*CUENTR_ZRG
        ZMF = PMFU(JLON, IK + 1)*ZDZ
        CUENTR_LLO1 = IK < KCBOT(JLON)
        IF (CUENTR_LLO1) THEN
          ZXDETRPEN = YDECUMF%DETRPEN
          ZDMFEN = ZENTR*ZMF
          ZDMFDE = ZXDETRPEN*ZMF
        END IF
      END IF
      
    END IF
    
    ! =========================================
    
    !                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
    !                  ---------------------------------------------------
    
    IF (LLO3) THEN
      
      ZQOLD = 0.0_JPRB
      DO JLL=1,JLM
        JLON = JLX
        ZDMFDE = MIN(ZDMFDE, 0.75_JPRB*PMFU(JLON, JK + 1))
        IF (JK == KCBOT(JLON)) THEN
          ZXENTRORG = YDECUMF%ENTRORG
          ZOENTR = -ZXENTRORG*(MIN(1.0_JPRB, PQEN(JLON, JK) / PQSEN(JLON, JK)) - YDECUMF%ENTR_RH)*(PGEOH(JLON, JK) - PGEOH(JLON,  &
          & JK + 1))*ZRG
          ZOENTR = MIN(0.4_JPRB, ZOENTR)*PMFU(JLON, JK + 1)
        END IF
        IF (JK < KCBOT(JLON)) THEN
          ZMFMAX = (PAPH(JLON, JK) - PAPH(JLON, JK - 1))*ZCONS2
          ZXS = MAX(PMFU(JLON, JK + 1) - ZMFMAX, 0.0_JPRB)
          PWMEAN(JLON) = PWMEAN(JLON) + PKINEU(JLON, JK + 1)*(PAP(JLON, JK + 1) - PAP(JLON, JK))
          ZDPMEAN = ZDPMEAN + PAP(JLON, JK + 1) - PAP(JLON, JK)
          ZDMFEN = ZOENTR
          IF (KTYPE(JLON) >= 2) THEN
            ZXENTSHALP = YDECUMF%ENTSHALP
            ZDMFEN = ZXENTSHALP*ZDMFEN
            ZDMFDE = ZDMFEN
            ! ZDMFDE(JL)=MAX(ZDMFDE(JL),ZDMFEN(JL))
          END IF
          ZDMFDE = ZDMFDE*(1.6_JPRB - MIN(1.0_JPRB, PQEN(JLON, JK) / PQSEN(JLON, JK)))
          ZMFTEST = PMFU(JLON, JK + 1) + ZDMFEN - ZDMFDE
          ZCHANGE = MAX(ZMFTEST - ZMFMAX, 0.0_JPRB)
          ZXE = MAX(ZCHANGE - ZXS, 0.0_JPRB)
          ZDMFEN = ZDMFEN - ZXE
          ZCHANGE = ZCHANGE - ZXE
          ZDMFDE = ZDMFDE + ZCHANGE
        END IF
        
        PDMFEN(JLON, JK) = ZDMFEN - ZDMFDE
        
        PMFU(JLON, JK) = PMFU(JLON, JK + 1) + ZDMFEN - ZDMFDE
        ZQEEN = PQENH(JLON, JK + 1)*ZDMFEN
        ZSEEN = (YDCST%RCPD*PTENH(JLON, JK + 1) + PGEOH(JLON, JK + 1))*ZDMFEN
        IF (PLITOT(JLON, JK) > YDECLDP%RLMIN) THEN
          ZLEEN = PLITOT(JLON, JK)*ZDMFEN
        ELSE
          ZLEEN = 0.0_JPRB
        END IF
        ZSCDE = (YDCST%RCPD*PTU(JLON, JK + 1) + PGEOH(JLON, JK + 1))*ZDMFDE
        ZQUDE = PQU(JLON, JK + 1)*ZDMFDE
        PLUDE(JLON, JK) = PLU(JLON, JK + 1)*ZDMFDE
        ZMFUSK = PMFUS(JLON, JK + 1) + ZSEEN - ZSCDE
        ZMFUQK = PMFUQ(JLON, JK + 1) + ZQEEN - ZQUDE
        ZMFULK = PMFUL(JLON, JK + 1) + ZLEEN - PLUDE(JLON, JK)
        ZFAC = 1.0_JPRB / MAX(YDECUMF%RMFCMIN, PMFU(JLON, JK))
        PLU(JLON, JK) = ZMFULK*ZFAC
        PQU(JLON, JK) = ZMFUQK*ZFAC
        PTU(JLON, JK) = (ZMFUSK*ZFAC - PGEOH(JLON, JK))*ZORCPD
        PTU(JLON, JK) = MAX(100._JPRB, PTU(JLON, JK))
        PTU(JLON, JK) = MIN(400._JPRB, PTU(JLON, JK))
        ZQOLD = PQU(JLON, JK)
        PLRAIN(JLON, JK) = PLRAIN(JLON, JK + 1)*MAX(0.0_JPRB, PMFU(JLON, JK + 1) - ZDMFDE)*ZFAC
        ZLUOLD = PLU(JLON, JK)
      END DO
      ! reset to environmental values if below departure level
      IF (JK > KDPL(JLON)) THEN
        PTU(JLON, JK) = PTENH(JLON, JK)
        PQU(JLON, JK) = PQENH(JLON, JK)
        PLU(JLON, JK) = 0.0_JPRB
        ZLUOLD = PLU(JLON, JK)
      END IF
      
      !                  DO CORRECTIONS FOR MOIST ASCENT
      !                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
      !                  -----------------------------------
      
      IK = JK
      IF (JLM > 0) THEN
        IF (YDECUMF%LSCVLIQ) THEN
          ! [Loki] inlined child subroutine: CUADJTQ
          ! =========================================
          
          !----------------------------------------------------------------------
          
          !     1.           DEFINE CONSTANTS
          !                  ----------------
          
          
          !IF (LHOOK) CALL DR_HOOK('CUADJTQ',0,ZHOOK_HANDLE)
          
          
          ZQMAX = 0.5_JPRB
          
          !*********************************************
          IF (.not.YDEPHLI%LPHYLIN) THEN
            !*********************************************
            
            !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
            !                  -----------------------------------------------------
            
            
            !   mixed phase saturation
            
            CUADJTQ_LLFLAG = LLFLAG
            IF (.true.) THEN
              CUADJTQ_LLFLAG = CUADJTQ_LLFLAG .and. .not.LSCVFLAG(JLON)
            ELSE
              CALL ABOR1_ACC('CUADJTQ: LDOFLAG has to be present when KCALL==6')
            END IF
            
            IF (CUADJTQ_LLFLAG) THEN
              ZQP = 1.0_JPRB / ZPH
              ZL = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4LES)
              ZI = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4IES)
              !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
              ZQSAT = YDTHF%R2ES*(FOEALFCU(PTU(JLON, IK))*EXP(YDTHF%R3LES*(PTU(JLON, IK) - YDCST%RTT)*ZL) + (1.0_JPRB -  &
              & FOEALFCU(PTU(JLON, IK)))*EXP(YDTHF%R3IES*(PTU(JLON, IK) - YDCST%RTT)*ZI))
              ZQSAT = ZQSAT*ZQP
              ZQSAT = MIN(0.5_JPRB, ZQSAT)
              ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
              ZF = FOEALFCU(PTU(JLON, IK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(PTU(JLON, IK)))*YDTHF%R5ALSCP*ZI**2
              ZCOND = (PQU(JLON, IK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
              ZCOND = MAX(ZCOND, 0.0_JPRB)
              PTU(JLON, IK) = PTU(JLON, IK) + FOELDCPMCU(PTU(JLON, IK))*ZCOND
              PQU(JLON, IK) = PQU(JLON, IK) - ZCOND
              ZL = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4LES)
              ZI = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4IES)
              !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
              ZQSAT = YDTHF%R2ES*(FOEALFCU(PTU(JLON, IK))*EXP(YDTHF%R3LES*(PTU(JLON, IK) - YDCST%RTT)*ZL) + (1.0_JPRB -  &
              & FOEALFCU(PTU(JLON, IK)))*EXP(YDTHF%R3IES*(PTU(JLON, IK) - YDCST%RTT)*ZI))
              ZQSAT = ZQSAT*ZQP
              ZQSAT = FMINJ(0.5_JPRB, ZQSAT)
              ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
              ZF = FOEALFCU(PTU(JLON, IK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(PTU(JLON, IK)))*YDTHF%R5ALSCP*ZI**2
              ZCOND1 = (PQU(JLON, IK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
              IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
              PTU(JLON, IK) = PTU(JLON, IK) + FOELDCPMCU(PTU(JLON, IK))*ZCOND1
              PQU(JLON, IK) = PQU(JLON, IK) - ZCOND1
            END IF
            
            
            IF (LLFLAG .and. LSCVFLAG(JLON)) THEN
              ZQP = 1.0_JPRB / ZPH
              ZL = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4LES)
              ZQSAT = YDTHF%R2ES*EXP(YDTHF%R3LES*(PTU(JLON, IK) - YDCST%RTT)*ZL)
              ZQSAT = ZQSAT*ZQP
              ZQSAT = MIN(0.5_JPRB, ZQSAT)
              ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
              ZF = YDTHF%R5ALVCP*ZL**2
              ZCOND = (PQU(JLON, IK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
              ZCOND = MAX(ZCOND, 0.0_JPRB)
              PTU(JLON, IK) = PTU(JLON, IK) + YDTHF%RALVDCP*ZCOND
              PQU(JLON, IK) = PQU(JLON, IK) - ZCOND
              ZL = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4LES)
              ZQSAT = YDTHF%R2ES*EXP(YDTHF%R3LES*(PTU(JLON, IK) - YDCST%RTT)*ZL)
              ZQSAT = ZQSAT*ZQP
              ZQSAT = FMINJ(0.5_JPRB, ZQSAT)
              ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
              ZF = YDTHF%R5ALVCP*ZL**2
              ZCOND1 = (PQU(JLON, IK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
              IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
              PTU(JLON, IK) = PTU(JLON, IK) + YDTHF%RALVDCP*ZCOND1
              PQU(JLON, IK) = PQU(JLON, IK) - ZCOND1
            END IF
            
            
            
            
            
            
            
            !*********************************************
          ELSE
            !*********************************************
            
            !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
            !                  -----------------------------------------------------
            
            
            
            
            
            !*********************************************
          END IF
          !*********************************************
          
          
          !IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
          
          ! =========================================
        ELSE
          ! [Loki] inlined child subroutine: CUADJTQ
          ! =========================================
          
          !----------------------------------------------------------------------
          
          !     1.           DEFINE CONSTANTS
          !                  ----------------
          
          
          !IF (LHOOK) CALL DR_HOOK('CUADJTQ',0,ZHOOK_HANDLE)
          
          
          ZQMAX = 0.5_JPRB
          
          !*********************************************
          IF (.not.YDEPHLI%LPHYLIN) THEN
            !*********************************************
            
            !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
            !                  -----------------------------------------------------
            
            
            !   mixed phase saturation
            
            CUADJTQ_LLFLAG = LLFLAG
            
            IF (CUADJTQ_LLFLAG) THEN
              ZQP = 1.0_JPRB / ZPH
              ZL = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4LES)
              ZI = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4IES)
              !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
              ZQSAT = YDTHF%R2ES*(FOEALFCU(PTU(JLON, IK))*EXP(YDTHF%R3LES*(PTU(JLON, IK) - YDCST%RTT)*ZL) + (1.0_JPRB -  &
              & FOEALFCU(PTU(JLON, IK)))*EXP(YDTHF%R3IES*(PTU(JLON, IK) - YDCST%RTT)*ZI))
              ZQSAT = ZQSAT*ZQP
              ZQSAT = MIN(0.5_JPRB, ZQSAT)
              ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
              ZF = FOEALFCU(PTU(JLON, IK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(PTU(JLON, IK)))*YDTHF%R5ALSCP*ZI**2
              ZCOND = (PQU(JLON, IK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
              ZCOND = MAX(ZCOND, 0.0_JPRB)
              PTU(JLON, IK) = PTU(JLON, IK) + FOELDCPMCU(PTU(JLON, IK))*ZCOND
              PQU(JLON, IK) = PQU(JLON, IK) - ZCOND
              ZL = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4LES)
              ZI = 1.0_JPRB / (PTU(JLON, IK) - YDTHF%R4IES)
              !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
              ZQSAT = YDTHF%R2ES*(FOEALFCU(PTU(JLON, IK))*EXP(YDTHF%R3LES*(PTU(JLON, IK) - YDCST%RTT)*ZL) + (1.0_JPRB -  &
              & FOEALFCU(PTU(JLON, IK)))*EXP(YDTHF%R3IES*(PTU(JLON, IK) - YDCST%RTT)*ZI))
              ZQSAT = ZQSAT*ZQP
              ZQSAT = FMINJ(0.5_JPRB, ZQSAT)
              ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
              ZF = FOEALFCU(PTU(JLON, IK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(PTU(JLON, IK)))*YDTHF%R5ALSCP*ZI**2
              ZCOND1 = (PQU(JLON, IK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
              IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
              PTU(JLON, IK) = PTU(JLON, IK) + FOELDCPMCU(PTU(JLON, IK))*ZCOND1
              PQU(JLON, IK) = PQU(JLON, IK) - ZCOND1
            END IF
            
            
            
            
            
            
            
            !*********************************************
          ELSE
            !*********************************************
            
            !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
            !                  -----------------------------------------------------
            
            
            !DIR$    IVDEP
            !OCL NOVREC
            !DIR$ LOOP_INFO EST_TRIPS(16)
            IF (LLFLAG) THEN
              ZQP = 1.0_JPRB / ZPH
              ZTARG = PTU(JLON, IK)
              CUADJTQ_ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
              ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
              ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
              ZQSAT = ZQP*(CUADJTQ_ZOEALFA*ZFOEEWL + (1.0_JPRB - CUADJTQ_ZOEALFA)*ZFOEEWI)
              Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
              ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
              
              ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
              ZQSAT = ZQSAT*ZCOR
              
              Z2S = CUADJTQ_ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - CUADJTQ_ZOEALFA) &
              & *YDTHF%R5ALSCP*(1.0_JPRB / (ZTARG - YDTHF%R4IES)**2)
              ZCOND = (PQU(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
              
              ZCOND = MAX(ZCOND, 0.0_JPRB)
              
              IF (ZCOND /= 0.0_JPRB) THEN
                
                PTU(JLON, IK) =  &
                & PTU(JLON, IK) + (CUADJTQ_ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - CUADJTQ_ZOEALFA)*YDTHF%RALSDCP)*ZCOND
                PQU(JLON, IK) = PQU(JLON, IK) - ZCOND
                ZTARG = PTU(JLON, IK)
                CUADJTQ_ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
                ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
                ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
                ZQSAT = ZQP*(CUADJTQ_ZOEALFA*ZFOEEWL + (1.0_JPRB - CUADJTQ_ZOEALFA)*ZFOEEWI)
                Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
                ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
                
                ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
                ZQSAT = ZQSAT*ZCOR
                
                Z2S = CUADJTQ_ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - CUADJTQ_ZOEALFA) &
                & *YDTHF%R5ALSCP*(1.0_JPRB / (ZTARG - YDTHF%R4IES)**2)
                ZCOND1 = (PQU(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
                
                PTU(JLON, IK) =  &
                & PTU(JLON, IK) + (CUADJTQ_ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - CUADJTQ_ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
                
                PQU(JLON, IK) = PQU(JLON, IK) - ZCOND1
              END IF
            END IF
            
            
            
            
            
            !*********************************************
          END IF
          !*********************************************
          
          
          !IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
          
          ! =========================================
        END IF
      END IF
      
      IF (YDEPHLI%LPHYLIN) THEN
        
        !DIR$ IVDEP
        !OCL NOVREC
        DO JLL=1,JLM
          JLON = JLX
          IF (PQU(JLON, JK) /= ZQOLD) THEN
            ZOEALFA = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(PTU(JLON, JK) - YDEPHLI%RLPTRC)) + 1.0_JPRB))
            ZOEALFAP = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(PTU(JLON, JK + 1) - YDEPHLI%RLPTRC)) + 1.0_JPRB))
            PLGLAC(JLON, JK) = PLU(JLON, JK)*((1.0_JPRB - ZOEALFA) - (1.0_JPRB - ZOEALFAP))
            ! add glaciation of rain
            ZFAC = 0.545_JPRB*(TANH(0.17_JPRB*(PTEN(JLON, JK) - YDEPHLI%RLPTRC)) + 1.0_JPRB)
            PLGLAC(JLON, JK) = PLGLAC(JLON, JK) + ZFAC*PDMFUP(JLON, JK + 1) / MAX(YDECUMF%RMFCMIN, PMFU(JLON, JK + 1))*(0.5_JPRB  &
            & + SIGN(0.5_JPRB, YDCST%RTT - PTEN(JLON, JK)))*ZGLAC
            PTU(JLON, JK) = PTU(JLON, JK) + YDTHF%RALFDCP*PLGLAC(JLON, JK)
          END IF
        END DO
        
      ELSE
        
        !DIR$ IVDEP
        !OCL NOVREC
        DO JLL=1,JLM
          JLON = JLX
          IF (PQU(JLON, JK) /= ZQOLD) THEN
            PLGLAC(JLON, JK) = PLU(JLON, JK)*((1.0_JPRB - FOEALFCU(PTU(JLON, JK))) - (1.0_JPRB - FOEALFCU(PTU(JLON, JK + 1))))
            IF (LSCVFLAG(JLON)) PLGLAC(JLON, JK) = 0.0_JPRB
            ! add glaciation of rain, only fraction added to updraught heat
            ZFAC = FOEALFCU(PTEN(JLON, JK))
            PLGLAC(JLON, JK) = PLGLAC(JLON, JK) + ZFAC*PDMFUP(JLON, JK + 1) / MAX(YDECUMF%RMFCMIN, PMFU(JLON, JK + 1))*(0.5_JPRB  &
            & + SIGN(0.5_JPRB, YDCST%RTT - PTEN(JLON, JK)))*ZGLAC
            PTU(JLON, JK) = PTU(JLON, JK) + YDTHF%RALFDCP*PLGLAC(JLON, JK)
          END IF
        END DO
        
      END IF
      
      DO JLL=1,JLM
        JLON = JLX
        IF (PQU(JLON, JK) /= ZQOLD) THEN
          KLAB(JLON, JK) = 2
          PLU(JLON, JK) = PLU(JLON, JK) + ZQOLD - PQU(JLON, JK)
          ZBC = PTU(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQU(JLON, JK) - PLU(JLON, JK + 1) - PLRAIN(JLON, JK + 1))
          ZBE = PTENH(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JLON, JK))
          ZBUO(JLON, JK) = ZBC - ZBE
          
          ! set flags in case of midlevel convection
          
          IF (KTYPE(JLON) == 3 .and. KLAB(JLON, JK + 1) == 1) THEN
            IF (ZBUO(JLON, JK) > -0.5_JPRB) THEN
              LDCUM(JLON) = .true.
              KCTOP(JLON) = JK
              PKINEU(JLON, JK) = 0.5_JPRB
            ELSE
              KLAB(JLON, JK) = 0
              PMFU(JLON, JK) = 0.0_JPRB
              PLUDE(JLON, JK) = 0.0_JPRB
              PLU(JLON, JK) = 0.0_JPRB
            END IF
          END IF
          
          IF (KLAB(JLON, JK + 1) == 2) THEN
            
            !IF(ZBUO(JL,JK) < 0.0_JPRB.AND.KTYPE(JL)==1) THEN
            IF (ZBUO(JLON, JK) < -0.1_JPRB) THEN
              PTENH(JLON, JK) = 0.5_JPRB*(PTEN(JLON, JK) + PTEN(JLON, JK - 1))
              PQENH(JLON, JK) = 0.5_JPRB*(PQEN(JLON, JK) + PQEN(JLON, JK - 1))
              ZBUO(JLON, JK) = ZBC - PTENH(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JLON, JK))
            END IF
            ZBUOC = 0.5*(ZBUO(JLON, JK) + ZBUO(JLON, JK + 1)) / (PTENH(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JLON, JK)))
            ZDKBUO = (PGEOH(JLON, JK) - PGEOH(JLON, JK + 1))*ZFACBUO*ZBUOC
            
            ! either use entrainment rate or if zero
            ! use detrainmnet rate as a substitute for
            ! mixing and "pressure" gradient term in upper
            ! troposphere
            
            IF (ZDMFEN > 0.0_JPRB) THEN
              ZDKEN = MIN(1.0_JPRB, (1 + Z_CWDRAG)*ZDMFEN / MAX(YDECUMF%RMFCMIN, PMFU(JLON, JK + 1)))
            ELSE
              ZDKEN = MIN(1.0_JPRB, (1 + Z_CWDRAG)*ZDMFDE / MAX(YDECUMF%RMFCMIN, PMFU(JLON, JK + 1)))
            END IF
            
            IF (LDTDKMF) THEN
              PKINEU(JLON, JK) = (PKINEU(JLON, JK + 1)*(1 - ZDKEN) + ZDKBUO) / (1 + ZDKEN)
              ZBUOC = ZBUO(JLON, JK)
            ELSE
              PKINEU(JLON, JK) = MAX(-1.E3_JPRB, (PKINEU(JLON, JK + 1)*(1 - ZDKEN) + ZDKBUO) / (1 + ZDKEN))
            END IF
            IF (ZBUOC < 0.0_JPRB) THEN
              ZKEDKE = PKINEU(JLON, JK) / MAX(1.E-3_JPRB, PKINEU(JLON, JK + 1))
              ZKEDKE = MAX(0.0_JPRB, MIN(1.0_JPRB, ZKEDKE))
              ZOCUDET = (1.6_JPRB - MIN(1.0_JPRB, PQEN(JLON, JK) / PQSEN(JLON, JK)))
              ZMFUN = ZOCUDET*SQRT(ZKEDKE)
              ZDMFDE = MAX(ZDMFDE, PMFU(JLON, JK + 1)*(1.0_JPRB - ZMFUN))
              PLUDE(JLON, JK) = PLU(JLON, JK + 1)*ZDMFDE
              PMFU(JLON, JK) = PMFU(JLON, JK + 1) + ZDMFEN - ZDMFDE
            END IF
            
            IF (ZBUO(JLON, JK) > -0.2_JPRB) THEN
              IKB = KCBOT(JLON)
              ZXENTRORG = YDECUMF%ENTRORG
              ZOENTR = ZXENTRORG*(0.3_JPRB - (MIN(1.0_JPRB, PQEN(JLON, JK - 1) / PQSEN(JLON, JK - 1)) - YDECUMF%ENTR_RH)) &
              & *(PGEOH(JLON, JK - 1) - PGEOH(JLON, JK))*ZRG*MIN(1.0_JPRB, PQSEN(JLON, JK) / PQSEN(JLON, IKB))**3
              ZOENTR = MIN(0.4_JPRB, ZOENTR)*PMFU(JLON, JK)
            ELSE
              ZOENTR = 0.0_JPRB
            END IF
            
            ! Erase values if below departure level
            IF (JK > KDPL(JLON)) THEN
              PMFU(JLON, JK) = PMFU(JLON, JK + 1)
              PKINEU(JLON, JK) = 0.5_JPRB
            END IF
            IF (PKINEU(JLON, JK) > 0.0_JPRB .and. PMFU(JLON, JK) > 0.0_JPRB .and. (ZBUO(JLON, JK) > -2._JPRB .or. (PTEN(JLON, JK  &
            & - 1) - PTEN(JLON, JK)) / (ZRG*(PGEO(JLON, JK - 1) - PGEO(JLON, JK))) < -3.E-3_JPRB)) THEN
              ! add overshoot limiter for convective top evaluation
              KCTOP(JLON) = JK
              LLO1 = .true.
            ELSE
              KLAB(JLON, JK) = 0
              PMFU(JLON, JK) = 0.0_JPRB
              PKINEU(JLON, JK) = 0.0_JPRB
              ZDMFDE = PMFU(JLON, JK + 1)
              PLUDE(JLON, JK) = PLU(JLON, JK + 1)*ZDMFDE
            END IF
            
            ! store detrainment rates for updraught
            
            IF (PMFU(JLON, JK + 1) > 0.0_JPRB) THEN
              PMFUDE_RATE(JLON, JK) = ZDMFDE
            END IF
            
          END IF
          
          !     ELSEIF(LLFLAG(JL).AND.KTYPE(JL)==2.AND.PQU(JL,JK) == ZQOLD(JL)) THEN
        ELSE IF (KTYPE(JLON) == 2 .and. PQU(JLON, JK) == ZQOLD) THEN
          KLAB(JLON, JK) = 0
          PMFU(JLON, JK) = 0.0_JPRB
          PKINEU(JLON, JK) = 0.0_JPRB
          ZDMFDE = PMFU(JLON, JK + 1)
          PLUDE(JLON, JK) = PLU(JLON, JK + 1)*ZDMFDE
          PMFUDE_RATE(JLON, JK) = ZDMFDE
          
        END IF
      END DO
      
      !              CALCULATE PRECIPITATION RATE BY
      !              ANALYTIC INTEGRATION OF EQUATION FOR L
      
      IF (LLO1) THEN
        IF (PLU(JLON, JK) > ZDNOPRC) THEN
          ZWU = MIN(15._JPRB, SQRT(2.0_JPRB*MAX(0.5_JPRB, PKINEU(JLON, JK + 1))))
          !       increase conversion for liquid phase only
          IF (LDTDKMF) THEN
            ZZCO = FOEALFCU(PTU(JLON, JK))*1.3_JPRB + (1.0_JPRB - FOEALFCU(PTU(JLON, JK)))
          ELSE
            ZZCO = 1.0_JPRB + 0.3_JPRB*FOEALFCU(PTU(JLON, JK))
          END IF
          ! IF(LSCVFLAG(JL)) ZZCO=1.3_JPRB
          
          ZXPRCDGW = ZPRCDGW
          ZPRCON = ZXPRCDGW / (0.75_JPRB*ZWU)*ZZCO
          
          !           PARAMETERS FOR BERGERON-FINDEISEN PROCESS (T < -5C)
          
          ZDT = MIN(YDTHF%RTBERCU - YDTHF%RTICECU, MAX(YDTHF%RTBERCU - PTU(JLON, JK), 0.0_JPRB))
          ZCBF = 1 + Z_CPRC2*SQRT(ZDT)
          ZZCO = ZPRCON*ZCBF
          
          IF (YGFL%NACTAERO > 0) THEN
            !         (Prognostic) aerosols in clouds, CCN effect only in pure warm phase
            !         then phases out in mixed phase - ***NOT*** linear physics
            ! CCN number taken at each height from environment
            IF (YDECLDP%LAERLIQAUTOCP) THEN
              ZALFAW = FOEALFCU(PTU(JLON, JK))
              ZLCRIT = ZALFAW*PLCRIT_AER(JLON, JK) + (1.0_JPRB - ZALFAW)*ZDNOPRC
              ! CCN number determined at cloud base (entrainment ignored)
            ELSE IF (YDECLDP%LAERLIQAUTOCPB) THEN
              ZALFAW = FOEALFCU(PTU(JLON, JK))
              ZLCRIT = ZALFAW*PLCRIT_AER(JLON, KCBOT(JLON)) + (1.0_JPRB - ZALFAW)*ZDNOPRC
            ELSE
              ZLCRIT = ZDNOPRC / ZCBF
            END IF
          ELSE
            ZLCRIT = ZDNOPRC / ZCBF
          END IF
          
          ZDFI = PGEOH(JLON, JK) - PGEOH(JLON, JK + 1)
          ZC = (PLU(JLON, JK) - ZLUOLD)
          ZD = ZZCO*(1.0_JPRB - EXP(-(PLU(JLON, JK) / ZLCRIT)**2))*ZDFI
          ZINT = EXP(-ZD)
          ZLNEW = ZLUOLD*ZINT + ZC / ZD*(1.0_JPRB - ZINT)
          ZLNEW = MAX(0.0_JPRB, MIN(PLU(JLON, JK), ZLNEW))
          ZLNEW = MIN(Z_CLDMAX, ZLNEW)
          ZPRECIP = MAX(0.0_JPRB, ZLUOLD + ZC - ZLNEW)
          PDMFUP(JLON, JK) = ZPRECIP*PMFU(JLON, JK)
          PLRAIN(JLON, JK) = PLRAIN(JLON, JK) + ZPRECIP
          PLU(JLON, JK) = ZLNEW
        END IF
      END IF
      
      IF (YDEPHLI%LPHYLIN) THEN
        
        !DEC$ IVDEP
        IF (LLO1) THEN
          IF (PLRAIN(JLON, JK) > 0.0_JPRB) THEN
            ZVW = 21.18_JPRB*PLRAIN(JLON, JK)**0.2_JPRB
            ZVI = Z_CWIFRAC*ZVW
            ZALFAW = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(PTU(JLON, JK) - YDEPHLI%RLPTRC)) + 1.0_JPRB))
            ZVV = ZALFAW*ZVW + (1.0_JPRB - ZALFAW)*ZVI
            ZROLD = PLRAIN(JLON, JK) - ZPRECIP
            ZC = ZPRECIP
            ZWU = MIN(15._JPRB, SQRT(2.0_JPRB*MAX(0.1_JPRB, PKINEU(JLON, JK))))
            ZD = ZVV / ZWU
            ZINT = EXP(-ZD)
            ZRNEW = ZROLD*ZINT + ZC / ZD*(1.0_JPRB - ZINT)
            ZRNEW = MAX(0.0_JPRB, MIN(PLRAIN(JLON, JK), ZRNEW))
            PLRAIN(JLON, JK) = ZRNEW
          END IF
        END IF
        
      ELSE
        
        IF (LLO1) THEN
          IF (PLRAIN(JLON, JK) > 0.0_JPRB) THEN
            ZVW = 21.18_JPRB*PLRAIN(JLON, JK)**0.2_JPRB
            ZVI = Z_CWIFRAC*ZVW
            ZALFAW = FOEALFCU(PTU(JLON, JK))
            IF (LSCVFLAG(JLON)) ZALFAW = 1.0_JPRB
            ZVV = ZALFAW*ZVW + (1.0_JPRB - ZALFAW)*ZVI
            ZROLD = PLRAIN(JLON, JK) - ZPRECIP
            ZC = ZPRECIP
            ZWU = MIN(15._JPRB, SQRT(2.0_JPRB*MAX(0.1_JPRB, PKINEU(JLON, JK))))
            ZD = ZVV / ZWU
            ZINT = EXP(-ZD)
            ZRNEW = ZROLD*ZINT + ZC / ZD*(1.0_JPRB - ZINT)
            ZRNEW = MAX(0.0_JPRB, MIN(PLRAIN(JLON, JK), ZRNEW))
            PLRAIN(JLON, JK) = ZRNEW
          END IF
        END IF
        
      END IF
      
      DO JLL=1,JLM
        JLON = JLX
        PMFUL(JLON, JK) = PLU(JLON, JK)*PMFU(JLON, JK)
        PMFUS(JLON, JK) = (YDCST%RCPD*PTU(JLON, JK) + PGEOH(JLON, JK))*PMFU(JLON, JK)
        PMFUQ(JLON, JK) = PQU(JLON, JK)*PMFU(JLON, JK)
        ZALFAW = FOEALFCU(PTU(JLON, JK))
        IF (LSCVFLAG(JLON)) ZALFAW = 1.0_JPRB
        PLUDELI(JLON, JK, 1) = ZALFAW*PLUDE(JLON, JK)
        PLUDELI(JLON, JK, 2) = (1.0_JPRB - ZALFAW)*PLUDE(JLON, JK)
      END DO
      
    END IF
  END DO
  
  !----------------------------------------------------------------------
  
  !     5.           FINAL CALCULATIONS
  !                  ------------------
  
  IF (KCTOP(JLON) == -1) LDCUM(JLON) = .false.
  KCBOT(JLON) = MAX(KCBOT(JLON), KCTOP(JLON))
  IF (LDCUM(JLON)) THEN
    PWMEAN(JLON) = MAX(1.E-2_JPRB, PWMEAN(JLON) / MAX(1.0_JPRB, ZDPMEAN))
    PWMEAN(JLON) = SQRT(2.0_JPRB*PWMEAN(JLON))
  END IF
  
  
  
  ! (C) Copyright 1988- ECMWF.
  !
  ! This software is licensed under the terms of the Apache Licence Version 2.0
  ! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
  !
  ! In applying this licence, ECMWF does not waive the privileges and immunities
  ! granted to it by virtue of its status as an intergovernmental organisation
  ! nor does it submit to any jurisdiction.
  
  
END SUBROUTINE CUASCN_OPENACC
