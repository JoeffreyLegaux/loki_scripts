SUBROUTINE CUBASEN_OPENACC (PPLDARE, PPLRG, YDTHF, YDCST, YDEPHLI, YDECLDP, YDECUMF, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV,  &
& KINDEX, LDMIXS, LDTDKMF, PTENH, PQENH, PGEOH, PAPH, PQHFL, PAHFS, PGP2DSPP, PKMFL, PTEN, PQEN, PQSEN, PGEO, PTU, PQU, PLU,  &
& PWU2H, PWUBASE, KLAB, LDCUM, LDSC, KCBOT, KBOTSC, KCTOP, KDPL, PCAPE, YDSTACK)
  
  !          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
  !          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT
  
  !          PURPOSE.
  !          --------
  !          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION
  
  !          INTERFACE
  !          ---------
  !          THIS ROUTINE IS CALLED FROM *CUMASTR*.
  !          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
  !          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
  !                 KLAB=0 FOR STABLE LAYERS
  !                 KLAB=1 FOR SUBCLOUD LEVELS
  !                 KLAB=2 FOR CLOUD LEVELS LEVEL
  
  !          METHOD.
  !          --------
  !          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
  !          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KINDEX*       TOP INDEX FOR OUTER VERTICAL LOOP
  
  !    INPUT PARAMETERS (LOGICAL):
  !    *LDMIXS*        WEAK (FALSE) OR STRONG (TRUE) CLOUD MIXING FOR SURFACE PARCEL ONLY
  !    *LDTDKMF*      Arpege tuning (if TRUE)
  
  !    INPUT PARAMETERS (REAL):
  
  ! not used at the moment because we want to use linear intepolation
  ! for fields on the half levels.
  
  !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
  !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
  
  !    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
  !    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
  !    *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2
  
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
  !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
  !    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
  !    *PQSEN*        PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  KG/KG
  !    *PGEO*         GEOPOTENTIAL                                  M2/S2
  !    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
  !    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
  !    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PTU*          TEMPERATURE IN UPDRAFTS                         K
  !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
  !    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
  !    *PWU2H*        KINETIC ENERGY IN UPDRAFTS  SURFACE PARCEL    M2/S2
  
  !    UPDATED PARAMETERS (INTEGER):
  
  !    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
  !                        KLAB=2 FOR CLOUD LEVELS
  
  !    OUTPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  !    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST
  
  !    OUTPUT PARAMETERS (INTEGER):
  
  !    *KCBOT*       CLOUD BASE LEVEL !
  !    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL
  !                  WITH A NON-ZERO CLOUD UPDRAFT.
  !    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
  !    *KDPL*        DEPARTURE LEVEL
  !    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)
  
  !          EXTERNALS
  !          ---------
  !          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
  
  !     AUTHOR.
  !     -------
  !      A. Pier Siebesma   KNMI ********
  
  !     MODIFICATIONS.
  !     --------------
  !      modified C Jakob (ECMWF) (01/2001)
  !      modified P Bechtold (ECMWF) (08/2002)
  !      02-11-02 : Use fixed last possible departure level and
  !                 last updraft computation level for bit-reproducibility
  !                 D.Salmond &  J. Hague
  !      03-07-03 : Tuning for p690     J. Hague
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      N.Semane+P.Bechtold    04-10-2012 Add RPLRG/RPLDARE factors for small planet
  !      M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
  !      S.-J. Lock (01 Nov 2016) SPP bug fix for ENTRORG perturbations
  !      M. Leutbecher (Oct 2020) SPP abstraction
  !      20210913 : Modifications for Arpege Y.Bouteloup (LDTDKMF)
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUBASEN_OPENACC ) seq
  
  USE YOEPHLI, ONLY: TEPHLI
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMCST, ONLY: TCST
  USE PARPHY, ONLY: RKAP
  USE YOECUMF, ONLY: TECUMF
  USE YOECLDP, ONLY: TECLDP
  USE YOETHF, ONLY: TTHF
  USE SPP_MOD, ONLY: TSPP_CONFIG
  USE SPP_GEN_MOD, ONLY: SPP_PERT
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  REAL(KIND=JPRB), INTENT(IN) :: PPLDARE
  REAL(KIND=JPRB), INTENT(IN) :: PPLRG
  TYPE(TECLDP), INTENT(IN) :: YDECLDP
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  TYPE(TSPP_CONFIG), INTENT(IN) :: YDSPP_CONFIG
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KINDEX
  LOGICAL, INTENT(IN) :: LDMIXS
  LOGICAL, INTENT(IN) :: LDTDKMF
  REAL(KIND=JPRB), INTENT(IN) :: PTENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PQHFL(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAHFS(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
  REAL(KIND=JPRB), INTENT(IN) :: PKMFL(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PWU2H(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PWUBASE(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KLAB(KLON, KLEV)
  LOGICAL, INTENT(INOUT) :: LDCUM(KLON)
  LOGICAL, INTENT(OUT) :: LDSC(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KBOTSC(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCTOP(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KDPL(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PCAPE(KLON)
  INTEGER(KIND=JPIM) :: ICTOP
  INTEGER(KIND=JPIM) :: ICBOT
  INTEGER(KIND=JPIM) :: IBOTSC
  temp (INTEGER (KIND=JPIM), ILAB, (KLON, KLEV))
  INTEGER(KIND=JPIM) :: IDPL
  
  LOGICAL :: LL_LDBASE
  LOGICAL :: LLGO_ON
  LOGICAL :: LLDEEP
  LOGICAL :: LLDCUM
  LOGICAL :: LLDSC
  LOGICAL :: LLFIRST
  LOGICAL :: LLRESET
  LOGICAL :: LLRESETJL
  LOGICAL :: LLPERT_ENTRORG  ! SPP perturbation on?
  LOGICAL :: LLPERT_ENTSTPC1  ! SPP perturbation on?
  
  
  
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IS
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  INTEGER(KIND=JPIM) :: JKK
  INTEGER(KIND=JPIM) :: JKT1
  INTEGER(KIND=JPIM) :: JKT2
  INTEGER(KIND=JPIM) :: JKT
  INTEGER(KIND=JPIM) :: JKB
  INTEGER(KIND=JPIM) :: IPENTRORG
  INTEGER(KIND=JPIM) :: IPENTSTPC1
  
  temp (REAL (KIND=JPRB), ZS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZSENH, (KLON, KLEV + 1))
  temp (REAL (KIND=JPRB), ZQENH, (KLON, KLEV + 1))
  temp (REAL (KIND=JPRB), ZSUH, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZBUOH, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZWU2H, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZQOLD
  REAL(KIND=JPRB) :: ZPH
  REAL(KIND=JPRB) :: ZMIX
  REAL(KIND=JPRB) :: ZDZ
  REAL(KIND=JPRB) :: ZCBASE
  
  temp (REAL (KIND=JPRB), ZLU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZQU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZTU, (KLON, KLEV))
  
  temp (REAL (KIND=JPRB), ZCAPE, (KLON, KLEV))
  
  REAL(KIND=JPRB) :: ZBUOF  ! BUOYANCY
  REAL(KIND=JPRB) :: ZRHO  ! DENSITY AT SURFACE (KG/M^3)
  REAL(KIND=JPRB) :: ZKHVFL  ! SURFACE BUOYANCY FLUX (K M/S)
  REAL(KIND=JPRB) :: ZWS  ! SIGMA_W AT LOWEST MODEL HALFLEVEL (M/S)
  REAL(KIND=JPRB) :: ZUST  ! U*
  REAL(KIND=JPRB) :: ZQEX  ! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
  REAL(KIND=JPRB) :: ZQEXC  ! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
  REAL(KIND=JPRB) :: ZTEX  ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
  REAL(KIND=JPRB) :: ZTEXC  ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
  REAL(KIND=JPRB) :: ZEPS  ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
  REAL(KIND=JPRB) :: ZTVENH  ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)
  REAL(KIND=JPRB) :: ZTVUH  ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
  REAL(KIND=JPRB) :: ZLGLAC  ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
  REAL(KIND=JPRB) :: ZQSU
  REAL(KIND=JPRB) :: ZCOR
  REAL(KIND=JPRB) :: ZDQ
  REAL(KIND=JPRB) :: ZALFAW
  REAL(KIND=JPRB) :: ZFACW
  REAL(KIND=JPRB) :: ZFACI
  REAL(KIND=JPRB) :: ZFAC
  REAL(KIND=JPRB) :: ZESDP
  REAL(KIND=JPRB) :: ZDQSDT
  REAL(KIND=JPRB) :: ZDTDP
  REAL(KIND=JPRB) :: ZDP
  REAL(KIND=JPRB) :: ZPDIFFTOP
  REAL(KIND=JPRB) :: ZPDIFFBOT
  REAL(KIND=JPRB) :: ZSF
  REAL(KIND=JPRB) :: ZQF
  REAL(KIND=JPRB) :: ZAW
  REAL(KIND=JPRB) :: ZBW
  REAL(KIND=JPRB) :: ZWORK1  ! work arrays for T and w perturbations
  REAL(KIND=JPRB) :: ZWORK2  ! work arrays for T and w perturbations
  REAL(KIND=JPRB) :: ZRCPD
  REAL(KIND=JPRB) :: ZRG
  REAL(KIND=JPRB) :: ZTMP
  
  REAL(KIND=JPRB) :: ZXENTRORG
  REAL(KIND=JPRB) :: ZXENTSTPC1
  INTEGER(KIND=JPIM) :: IPN  ! SPP perturbation pointer
  TYPE(SPP_PERT) :: PN1
  TYPE(SPP_PERT) :: PN2
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  
#include "cuadjtq.intfb.h"
#include "fcttre.func.h"
#include "abor1.intfb.h"
#include "cuadjtq.func.h"
  INTEGER(KIND=JPIM) :: CUADJTQ_JL
  REAL(KIND=JPRB) :: Z1S
  REAL(KIND=JPRB) :: Z2S
  REAL(KIND=JPRB) :: ZCOND
  REAL(KIND=JPRB) :: ZCOND1
  REAL(KIND=JPRB) :: CUADJTQ_ZCOR
  REAL(KIND=JPRB) :: ZFOEEWI
  REAL(KIND=JPRB) :: ZFOEEWL
  REAL(KIND=JPRB) :: ZOEALFA
  REAL(KIND=JPRB) :: ZQMAX
  REAL(KIND=JPRB) :: ZQSAT
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZQP
  REAL(KIND=JPRB) :: ZL
  REAL(KIND=JPRB) :: ZI
  REAL(KIND=JPRB) :: ZF
  LOGICAL :: LLFLAG
  REAL(KIND=JPHOOK) :: CUADJTQ_ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ILAB)
  alloc (ZS)
  alloc (ZSENH)
  alloc (ZQENH)
  alloc (ZSUH)
  alloc (ZBUOH)
  alloc (ZWU2H)
  alloc (ZLU)
  alloc (ZQU)
  alloc (ZTU)
  alloc (ZCAPE)
  JLON = KIDIA
  
  !----------------------------------------------------------------------
  !     0.           INITIALIZE CONSTANTS AND FIELDS
  !                  -------------------------------
  !----------------------------------------------------------------------
  
  ZAW = 1.0_JPRB
  ZBW = 1.0_JPRB
  
  PWUBASE(JLON) = 0.0_JPRB
  LLGO_ON = .true.
  LLFIRST = .true.
  KDPL(JLON) = KLEV
  
  JKT1 = KINDEX
  JKT2 = YDECUMF%NJKT2
  ZRG = 1.0_JPRB / YDCST%RG
  ZRCPD = 1.0_JPRB / YDCST%RCPD
  
  DO JK=1,KLEV
    ZTU(JLON, JK) = PTU(JLON, JK)
    ZQU(JLON, JK) = PQU(JLON, JK)
    ZLU(JLON, JK) = PLU(JLON, JK)
    ILAB(JLON, JK) = KLAB(JLON, JK)
    ZCAPE(JLON, JK) = 0.0_JPRB
  END DO
  
  ! prepare SPP perturbations
  LLPERT_ENTSTPC1 = .false.
  LLPERT_ENTRORG = .false.
  
  !----------------------------------------------------------------------
  !       -----------------------------------------------------------
  !       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
  !             OF SPECIFIC HUMIDITY AND STATIC ENERGY
  !       -----------------------------------------------------------
  
  DO JK=1,KLEV
    PWU2H(JLON, JK) = 0.0_JPRB
    ZWU2H(JLON, JK) = 0.0_JPRB
    ZS(JLON, JK) = YDCST%RCPD*PTEN(JLON, JK) + PGEO(JLON, JK)
    ZQENH(JLON, JK) = PQENH(JLON, JK)
    ZSENH(JLON, JK) = YDCST%RCPD*PTENH(JLON, JK) + PGEOH(JLON, JK)
  END DO
  
  DO JKK=KLEV,JKT1,-1
    ! Big external loop for level testing:
    ! find first departure level that produces deepest cloud top
    ! or take surface level for shallow convection and Sc
    !
    !        ---------------------------------------------------------
    !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
    !        ---------------------------------------------------------
    !
    IS = 0
    IF (LLGO_ON) THEN
      IS = IS + 1
      IDPL = JKK        ! departure level
      ICBOT = JKK        ! cloud base level for convection, (-1 if not found)
      IBOTSC = KLEV - 1        ! sc    base level for sc-clouds , (-1 if not found)
      ICTOP = KLEV - 1        ! cloud top for convection (-1 if not found)
      LLDCUM = .false.        ! on exit: true if cloudbase=found
      LLDSC = .false.        ! on exit: true if cloudbase=found
      LL_LDBASE = .false.        ! on exit: true if cloudbase=found
    END IF
    
    IF (IS /= 0) THEN
      
      IF (JKK == KLEV) THEN
        
        ZTEXC = 0.2_JPRB
        ZQEXC = 1.E-4_JPRB
        IF (LLGO_ON) THEN
          ZRHO = PAPH(JLON, JKK + 1) / (YDCST%RD*(PTEN(JLON, JKK)*(1.0_JPRB + YDCST%RETV*PQEN(JLON, JKK))))
          ZKHVFL = (PAHFS(JLON, JKK + 1)*ZRCPD + YDCST%RETV*PTEN(JLON, JKK)*PQHFL(JLON, JKK + 1)) / (ZRHO*PPLRG*PPLDARE)
          ZUST = MAX(SQRT(PKMFL(JLON)), 0.1_JPRB)
          ZWS = ZUST**3._JPRB - 1.5_JPRB*RKAP*ZKHVFL*(PGEOH(JLON, KLEV) - PGEOH(JLON, KLEV + 1)) / PTEN(JLON, KLEV)
          ZTEX = 0.0_JPRB
          ZQEX = 0.0_JPRB
          IF (ZKHVFL < 0.0_JPRB) THEN
            ZWS = 1.2_JPRB*ZWS**.3333_JPRB
            ILAB(JLON, JKK) = 1
            ZTEX = MAX(-1.5_JPRB*PAHFS(JLON, JKK + 1) / (ZRHO*ZWS*YDCST%RCPD*PPLRG*PPLDARE), ZTEXC)
            ZQEX = MAX(-1.5_JPRB*PQHFL(JLON, JKK + 1) / (ZRHO*ZWS*PPLRG*PPLDARE), ZQEXC)
            ZQU(JLON, JKK) = ZQENH(JLON, JKK) + ZQEX
            ZSUH(JLON, JKK) = ZSENH(JLON, JKK) + YDCST%RCPD*ZTEX
            ZTU(JLON, JKK) = (ZSENH(JLON, JKK) - PGEOH(JLON, JKK))*ZRCPD + ZTEX
            ZLU(JLON, JKK) = 0.0_JPRB
            ZWU2H(JLON, JKK) = ZWS**2 + 0.1_JPRB
            PWU2H(JLON, JKK) = ZWU2H(JLON, JKK)
            !
            !  determine buoyancy at lowest half level
            !
            ZTVENH = (1.0_JPRB + YDCST%RETV*ZQENH(JLON, JKK))*(ZSENH(JLON, JKK) - PGEOH(JLON, JKK))*ZRCPD
            ZTVUH = (1.0_JPRB + YDCST%RETV*ZQU(JLON, JKK))*ZTU(JLON, JKK)
            ZBUOH(JLON, JKK) = (ZTVUH - ZTVENH)*YDCST%RG / ZTVENH
          ELSE
            LLGO_ON = .false.              ! non-convective point
          END IF
        END IF
        
      ELSE
        
        IF (LLGO_ON) THEN
          ZRHO = PAPH(JLON, JKK + 1) / (YDCST%RD*(PTEN(JLON, JKK)*(1. + YDCST%RETV*PQEN(JLON, JKK))))
          ILAB(JLON, JKK) = 1
          ZTEXC = 0.2_JPRB
          ZQEXC = 1.E-4_JPRB
          IF (JKK == KLEV - 1) THEN
            ZTEXC = MAX(ZTEXC, ZTEX)
            ZQEXC = MAX(ZQEXC, ZQEX)
            IF (LDTDKMF) THEN
              ZTEXC = MIN(ZTEXC, 3.0_JPRB)
              ZQEXC = MIN(ZQEXC, 2.E-3_JPRB)
            ELSE
              ZTEXC = MIN(ZTEXC, 1.0_JPRB)
              ZQEXC = MIN(ZQEXC, 5.E-4_JPRB)
            END IF
          END IF
          ZQU(JLON, JKK) = ZQENH(JLON, JKK) + ZQEXC
          ZSUH(JLON, JKK) = ZSENH(JLON, JKK) + YDCST%RCPD*ZTEXC
          ZTU(JLON, JKK) = (ZSENH(JLON, JKK) - PGEOH(JLON, JKK))*ZRCPD + ZTEXC
          ZLU(JLON, JKK) = 0.0_JPRB
          ! construct mixed layer for parcels emanating in lowest 60 hPa
          IF (PAPH(JLON, KLEV + 1) - PAPH(JLON, JKK - 1) < 60.E2_JPRB) THEN
            ZQU(JLON, JKK) = 0.0_JPRB
            ZSUH(JLON, JKK) = 0.0_JPRB
            ZWORK1 = 0.0_JPRB
            DO JK=JKK + 1,JKK - 1,-1
              IF (ZWORK1 < 50.E2_JPRB) THEN
                ZWORK2 = PAPH(JLON, JK) - PAPH(JLON, JK - 1)
                ZWORK1 = ZWORK1 + ZWORK2
                ZQU(JLON, JKK) = ZQU(JLON, JKK) + ZQENH(JLON, JK)*ZWORK2
                ZSUH(JLON, JKK) = ZSUH(JLON, JKK) + ZSENH(JLON, JK)*ZWORK2
              END IF
            END DO
            ZQU(JLON, JKK) = ZQU(JLON, JKK) / ZWORK1 + ZQEXC
            ZSUH(JLON, JKK) = ZSUH(JLON, JKK) / ZWORK1 + YDCST%RCPD*ZTEXC
            ZTU(JLON, JKK) = (ZSUH(JLON, JKK) - PGEOH(JLON, JKK))*ZRCPD + ZTEXC
          END IF
          ZWU2H(JLON, JKK) = 1.0_JPRB
          ! PWU2H(JL,JKK) = ZWU2H(JL,JKK)
          
          !
          !  determine buoyancy at lowest half level
          !
          ZTVENH = (1.0_JPRB + YDCST%RETV*ZQENH(JLON, JKK))*(ZSENH(JLON, JKK) - PGEOH(JLON, JKK))*ZRCPD
          ZTVUH = (1.0_JPRB + YDCST%RETV*ZQU(JLON, JKK))*ZTU(JLON, JKK)
          ZBUOH(JLON, JKK) = (ZTVUH - ZTVENH)*YDCST%RG / ZTVENH
        END IF
        
      END IF
      
    END IF
    
    !----------------------------------------------------------------------
    !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
    !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
    !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
    !                  CHECK FOR BUOYANCY AND SET FLAGS
    !                  -------------------------------------
    !       ------------------------------------------------------------
    !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
    !       ------------------------------------------------------------
    DO JK=JKK - 1,JKT2,-1
      IS = 0
      
      IF (JKK == KLEV .or. KINDEX == KLEV - 1) THEN
        ! 1/z mixing for shallow
        
        
        ZXENTSTPC1 = YDECUMF%ENTSTPC1
        
        IF (LLGO_ON) THEN
          IS = IS + 1
          ZDZ = (PGEOH(JLON, JK) - PGEOH(JLON, JK + 1))*ZRG
          IF (LDTDKMF) THEN
            ZEPS = ZXENTSTPC1 / ((PGEOH(JLON, JK) - PGEOH(JLON, KLEV + 1))*ZRG*PPLRG) + YDECUMF%ENTSTPC2
          ELSE
            ZEPS = ZXENTSTPC1 / ((PGEO(JLON, JK) - PGEOH(JLON, KLEV + 1))*ZRG*PPLRG) + YDECUMF%ENTSTPC2
            IF (LDMIXS .and. KINDEX == KLEV .and. ZLU(JLON, JK + 1) > 0.0_JPRB) ZEPS = ZEPS*YDECUMF%ENTSTPC3
          END IF
          ZMIX = 0.5_JPRB*ZDZ*PPLRG*ZEPS
          IF (.not.LDTDKMF) ZMIX = MIN(1.0_JPRB, ZMIX)
          ZQF = (PQENH(JLON, JK + 1) + PQENH(JLON, JK))*0.5_JPRB
          ZSF = (ZSENH(JLON, JK + 1) + ZSENH(JLON, JK))*0.5_JPRB
          ZTMP = 1.0_JPRB / (1.0_JPRB + ZMIX)
          ZQU(JLON, JK) = (ZQU(JLON, JK + 1)*(1.0_JPRB - ZMIX) + 2.0_JPRB*ZMIX*ZQF)*ZTMP
          ZSUH(JLON, JK) = (ZSUH(JLON, JK + 1)*(1.0_JPRB - ZMIX) + 2.0_JPRB*ZMIX*ZSF)*ZTMP
          ZQOLD = ZQU(JLON, JK)
          ZTU(JLON, JK) = (ZSUH(JLON, JK) - PGEOH(JLON, JK))*ZRCPD
          ZPH = PAPH(JLON, JK)
        END IF
        
      ELSE
        
        ZXENTRORG = YDECUMF%ENTRORG
        IF (LLGO_ON) THEN
          !SPP: perturb ENTRORG parameter
          ZDZ = (PGEOH(JLON, JK) - PGEOH(JLON, JK + 1))*ZRG
          ZMIX = 0.4_JPRB*ZXENTRORG*ZDZ*MIN(1.0_JPRB, (PQSEN(JLON, JK) / PQSEN(JLON, KLEV))**3)
        END IF
        
        IF (LLGO_ON) THEN
          IS = IS + 1
          ZMIX = MIN(1.0_JPRB, ZMIX)
          ZQF = (PQENH(JLON, JK + 1) + PQENH(JLON, JK))*0.5_JPRB
          ZSF = (ZSENH(JLON, JK + 1) + ZSENH(JLON, JK))*0.5_JPRB
          ZQU(JLON, JK) = ZQU(JLON, JK + 1)*(1.0_JPRB - ZMIX) + ZQF*ZMIX
          ZSUH(JLON, JK) = ZSUH(JLON, JK + 1)*(1.0_JPRB - ZMIX) + ZSF*ZMIX
          ZQOLD = ZQU(JLON, JK)
          ZTU(JLON, JK) = (ZSUH(JLON, JK) - PGEOH(JLON, JK))*ZRCPD
          ZPH = PAPH(JLON, JK)
        END IF
        
      END IF
      
      IF (IS == 0) EXIT
      
      IK = JK
      
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
        
        LLFLAG = LLGO_ON
        
        IF (LLFLAG) THEN
          ZQP = 1.0_JPRB / ZPH
          ZL = 1.0_JPRB / (ZTU(JLON, IK) - YDTHF%R4LES)
          ZI = 1.0_JPRB / (ZTU(JLON, IK) - YDTHF%R4IES)
          !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
          ZQSAT = YDTHF%R2ES*(FOEALFCU(ZTU(JLON, IK))*EXP(YDTHF%R3LES*(ZTU(JLON, IK) - YDCST%RTT)*ZL) + (1.0_JPRB -  &
          & FOEALFCU(ZTU(JLON, IK)))*EXP(YDTHF%R3IES*(ZTU(JLON, IK) - YDCST%RTT)*ZI))
          ZQSAT = ZQSAT*ZQP
          ZQSAT = MIN(0.5_JPRB, ZQSAT)
          CUADJTQ_ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
          ZF = FOEALFCU(ZTU(JLON, IK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(ZTU(JLON, IK)))*YDTHF%R5ALSCP*ZI**2
          ZCOND = (ZQU(JLON, IK)*CUADJTQ_ZCOR**2 - ZQSAT*CUADJTQ_ZCOR) / (CUADJTQ_ZCOR**2 + ZQSAT*ZF)
          ZCOND = MAX(ZCOND, 0.0_JPRB)
          ZTU(JLON, IK) = ZTU(JLON, IK) + FOELDCPMCU(ZTU(JLON, IK))*ZCOND
          ZQU(JLON, IK) = ZQU(JLON, IK) - ZCOND
          ZL = 1.0_JPRB / (ZTU(JLON, IK) - YDTHF%R4LES)
          ZI = 1.0_JPRB / (ZTU(JLON, IK) - YDTHF%R4IES)
          !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
          ZQSAT = YDTHF%R2ES*(FOEALFCU(ZTU(JLON, IK))*EXP(YDTHF%R3LES*(ZTU(JLON, IK) - YDCST%RTT)*ZL) + (1.0_JPRB -  &
          & FOEALFCU(ZTU(JLON, IK)))*EXP(YDTHF%R3IES*(ZTU(JLON, IK) - YDCST%RTT)*ZI))
          ZQSAT = ZQSAT*ZQP
          ZQSAT = FMINJ(0.5_JPRB, ZQSAT)
          CUADJTQ_ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
          ZF = FOEALFCU(ZTU(JLON, IK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(ZTU(JLON, IK)))*YDTHF%R5ALSCP*ZI**2
          ZCOND1 = (ZQU(JLON, IK)*CUADJTQ_ZCOR**2 - ZQSAT*CUADJTQ_ZCOR) / (CUADJTQ_ZCOR**2 + ZQSAT*ZF)
          IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
          ZTU(JLON, IK) = ZTU(JLON, IK) + FOELDCPMCU(ZTU(JLON, IK))*ZCOND1
          ZQU(JLON, IK) = ZQU(JLON, IK) - ZCOND1
        END IF
        
        
        
        
        
        
        
        !*********************************************
      ELSE
        !*********************************************
        
        !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
        !                  -----------------------------------------------------
        
        
        !DIR$    IVDEP
        !OCL NOVREC
        !DIR$ LOOP_INFO EST_TRIPS(16)
        IF (LLGO_ON) THEN
          ZQP = 1.0_JPRB / ZPH
          ZTARG = ZTU(JLON, IK)
          ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
          ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
          ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
          ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
          Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
          ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
          
          CUADJTQ_ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
          ZQSAT = ZQSAT*CUADJTQ_ZCOR
          
          Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB /  &
          & (ZTARG - YDTHF%R4IES)**2)
          ZCOND = (ZQU(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*CUADJTQ_ZCOR*Z2S)
          
          ZCOND = MAX(ZCOND, 0.0_JPRB)
          
          IF (ZCOND /= 0.0_JPRB) THEN
            
            ZTU(JLON, IK) = ZTU(JLON, IK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND
            ZQU(JLON, IK) = ZQU(JLON, IK) - ZCOND
            ZTARG = ZTU(JLON, IK)
            ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
            ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
            ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
            ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
            Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
            ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
            
            CUADJTQ_ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
            ZQSAT = ZQSAT*CUADJTQ_ZCOR
            
            Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB /  &
            & (ZTARG - YDTHF%R4IES)**2)
            ZCOND1 = (ZQU(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*CUADJTQ_ZCOR*Z2S)
            
            ZTU(JLON, IK) = ZTU(JLON, IK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
            
            ZQU(JLON, IK) = ZQU(JLON, IK) - ZCOND1
          END IF
        END IF
        
        
        
        
        
        !*********************************************
      END IF
      !*********************************************
      
      
      !IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
      
      ! =========================================
      
      !DIR$ IVDEP
      !OCL NOVREC
      
      IF (LLGO_ON) THEN
        
        ! add condensation to water
        
        ZDQ = MAX(ZQOLD - ZQU(JLON, JK), 0.0_JPRB)
        ZLU(JLON, JK) = ZLU(JLON, JK + 1) + ZDQ
        
        ! freezing
        
        ZLGLAC = ZDQ*((1.0_JPRB - FOEALFCU(ZTU(JLON, JK))) - (1.0_JPRB - FOEALFCU(ZTU(JLON, JK + 1))))
        
        
        ! pseudo-microphysics
        
        ZLU(JLON, JK) = 0.5_JPRB*ZLU(JLON, JK)
        
        ! update dry static energy after condensation + freezing
        
        ZTU(JLON, JK) = ZTU(JLON, JK) + YDTHF%RALFDCP*ZLGLAC
        ZSUH(JLON, JK) = YDCST%RCPD*ZTU(JLON, JK) + PGEOH(JLON, JK)
        
        ! Buoyancy on half and full levels
        
        ZTVUH = (1.0_JPRB + YDCST%RETV*ZQU(JLON, JK) - ZLU(JLON, JK))*ZTU(JLON, JK) + YDTHF%RALFDCP*ZLGLAC
        ZTVENH = (1.0_JPRB + YDCST%RETV*ZQENH(JLON, JK))*(ZSENH(JLON, JK) - PGEOH(JLON, JK))*ZRCPD
        ZBUOH(JLON, JK) = (ZTVUH - ZTVENH)*YDCST%RG / ZTVENH
        
        ZBUOF = (ZBUOH(JLON, JK) + ZBUOH(JLON, JK + 1))*0.5_JPRB
        
        ! solve kinetic energy equation
        
        ZTMP = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB*ZBW*ZMIX)
        ZWU2H(JLON, JK) = (ZWU2H(JLON, JK + 1)*(1.0_JPRB - 2.0_JPRB*ZBW*ZMIX) + 2.0_JPRB*ZAW*ZBUOF*ZDZ)*ZTMP
        
        IF (JKK == KLEV) THEN
          PWU2H(JLON, JK) = ZWU2H(JLON, JK)            !save surface parcel KE
        END IF
        
        ! compute CAPE for diagnostics
        
        ZCAPE(JLON, JKK) = ZCAPE(JLON, JKK) + MAX(0.0_JPRB, ZBUOF*ZDZ)
        
        ! first layer with liquid water - find exact cloud base
        
        IF (ZLU(JLON, JK) > 0.0_JPRB .and. ILAB(JLON, JK + 1) == 1) THEN
          
          IK = JK + 1
          ZQSU = FOEEWM(ZTU(JLON, IK)) / PAPH(JLON, IK)
          ZESDP = ZQSU
          ZQSU = MIN(0.5_JPRB, ZQSU)
          ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSU)
          ZQSU = ZQSU*ZCOR
          ZDQ = MIN(0._JPRB, ZQU(JLON, IK) - ZQSU)
          ZALFAW = FOEALFA(ZTU(JLON, IK))
          ZFACW = YDTHF%R5LES / ((ZTU(JLON, IK) - YDTHF%R4LES)**2)
          ZFACI = YDTHF%R5IES / ((ZTU(JLON, IK) - YDTHF%R4IES)**2)
          ZFAC = ZALFAW*ZFACW + (1. - ZALFAW)*ZFACI
          ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZESDP)
          ZDQSDT = ZFAC*ZCOR*ZQSU
          ZDTDP = YDCST%RD*ZTU(JLON, IK) / (YDCST%RCPD*PAPH(JLON, IK))
          ZDP = ZDQ / (ZDQSDT*ZDTDP)
          ZCBASE = PAPH(JLON, IK) + ZDP
          
          ! chose nearest half level as cloud base
          
          ZPDIFFTOP = ZCBASE - PAPH(JLON, JK)
          ZPDIFFBOT = PAPH(JLON, JK + 1) - ZCBASE
          
          IF (ZPDIFFTOP > ZPDIFFBOT .and. ZWU2H(JLON, JK + 1) > 0.0_JPRB) THEN
            JKB = MIN(KLEV - 1, JK + 1)
            ILAB(JLON, JKB) = 2
            ILAB(JLON, JK) = 2
            LL_LDBASE = .true.
            LLDSC = .true.
            IBOTSC = JKB
            ICBOT = JKB
            ZLU(JLON, JK + 1) = YDECLDP%RLMIN
          ELSE IF (ZPDIFFTOP <= ZPDIFFBOT .and. ZWU2H(JLON, JK) > 0.0_JPRB) THEN
            ILAB(JLON, JK) = 2
            LL_LDBASE = .true.
            LLDSC = .true.
            IBOTSC = JK
            ICBOT = JK
          END IF
          
        END IF
        
        ! decide on presence of convection, cloud base and cloud top based on
        ! kinetic energy
        
        IF (ZWU2H(JLON, JK) < 0.0_JPRB) THEN
          LLGO_ON = .false.
          IF (ZLU(JLON, JK + 1) > 0.0_JPRB) THEN
            ICTOP = JK
            LLDCUM = .true.
          ELSE
            LLDCUM = .false.
          END IF
        ELSE
          IF (ZLU(JLON, JK) > 0.0_JPRB) THEN
            ILAB(JLON, JK) = 2
          ELSE
            ILAB(JLON, JK) = 1
          END IF
        END IF
      END IF
      
      !     IF (IS == 0) EXIT
    END DO
    
    IF (JKK == KLEV .or. JKK == KLEV - 1 .and. KINDEX == KLEV - 1) THEN
      
      ! set values for departure level for PBL clouds = first model level
      LDSC(JLON) = LLDSC
      IF (LDSC(JLON)) THEN
        KBOTSC(JLON) = IBOTSC
      ELSE
        KBOTSC(JLON) = -1
      END IF
      
      JKT = ICTOP
      JKB = ICBOT
      LLDEEP = PAPH(JLON, JKB) - PAPH(JLON, JKT) > YDECUMF%RDEPTHS
      IF (LLDEEP .and. KINDEX < KLEV - 1) LLDCUM = .false.
      ! no deep allowed for KLEV
      LLDEEP = .false.        ! for deep convection start only at level KLEV-1
      ! and form mixed layer, so go on
      ! test further for deep convective columns as not yet found
      ! IF ( LLDEEP(JL) ) LLFIRST(JL)=.FALSE.
      LLGO_ON = .not.LLDEEP
      IF (KINDEX == KLEV - 1) LLGO_ON = .not.LLDCUM
      IF (LLDCUM) THEN
        KCBOT(JLON) = ICBOT
        KCTOP(JLON) = ICTOP
        KDPL(JLON) = IDPL
        LDCUM(JLON) = LLDCUM
        PWUBASE(JLON) = SQRT(MAX(ZWU2H(JLON, JKB), 0.0_JPRB))
      ELSE
        KCTOP(JLON) = -1
        KCBOT(JLON) = -1
        KDPL(JLON) = KLEV - 1
        LDCUM(JLON) = .false.
        PWUBASE(JLON) = 0.0_JPRB
      END IF
      DO JK=KLEV,1,-1
        JKT = ICTOP
        IF (JK >= JKT .or. KINDEX >= KLEV - 1 .and. JKK == KLEV) THEN
          KLAB(JLON, JK) = ILAB(JLON, JK)
          PTU(JLON, JK) = ZTU(JLON, JK)
          PQU(JLON, JK) = ZQU(JLON, JK)
          PLU(JLON, JK) = ZLU(JLON, JK)
        END IF
      END DO
    END IF
    
    IF (JKK < KLEV) THEN
      LLRESET = .false.
      IF (.not.LLDEEP) THEN
        JKT = ICTOP
        JKB = ICBOT
        ! test on cloud thickness and buoyancy
        LLDEEP = PAPH(JLON, JKB) - PAPH(JLON, JKT) >= YDECUMF%RDEPTHS
      END IF
      LLRESETJL = LLDEEP .and. LLFIRST
      LLRESET = LLRESET .or. LLRESETJL
      
      
      IF (LLRESET) THEN
        DO JK=KLEV,1,-1
          IF (LLRESETJL) THEN
            JKT = ICTOP
            JKB = IDPL
            IF (JK <= JKB .and. JK >= JKT) THEN
              KLAB(JLON, JK) = ILAB(JLON, JK)
              PTU(JLON, JK) = ZTU(JLON, JK)
              PQU(JLON, JK) = ZQU(JLON, JK)
              PLU(JLON, JK) = ZLU(JLON, JK)
            ELSE
              KLAB(JLON, JK) = 1
              PTU(JLON, JK) = PTENH(JLON, JK)
              PQU(JLON, JK) = PQENH(JLON, JK)
              PLU(JLON, JK) = 0.0_JPRB
            END IF
            IF (JK < JKT) KLAB(JLON, JK) = 0
          END IF
        END DO
      END IF
      
      IF (LLDEEP .and. LLFIRST) THEN
        KDPL(JLON) = IDPL
        KCTOP(JLON) = ICTOP
        KCBOT(JLON) = ICBOT
        LDCUM(JLON) = LLDCUM
        LDSC(JLON) = .false.
        KBOTSC(JLON) = -1
        JKB = KCBOT(JLON)
        PWUBASE(JLON) = SQRT(MAX(ZWU2H(JLON, JKB), 0.0_JPRB))
        !  no initialization of wind for deep here, this is done in
        !  CUINI and CUASCN
        LLFIRST = .false.
      END IF
      LLGO_ON = .not.LLDEEP
    END IF
    
  END DO
  ! end of big loop for search of departure level
  
  ! chose maximum CAPE value
  PCAPE(JLON) = MAXVAL(ZCAPE(JLON, :))
  
  
  
  ! (C) Copyright 1988- ECMWF.
  !
  ! This software is licensed under the terms of the Apache Licence Version 2.0
  ! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
  !
  ! In applying this licence, ECMWF does not waive the privileges and immunities
  ! granted to it by virtue of its status as an intergovernmental organisation
  ! nor does it submit to any jurisdiction.
  
  
END SUBROUTINE CUBASEN_OPENACC
