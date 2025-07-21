SUBROUTINE CUENTR_OPENACC (YDCST, YDECUMF, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, KK, KCBOT, KTYPE, LDCUM, LDWORK, PQSEN, PAPH,  &
& PGEOH, PMFU, PGP2DSPP, PDMFEN, PDMFDE, YDSTACK)
  
  !          M.TIEDTKE         E.C.M.W.F.     12/89
  !          P.BECHTOLD        E.C.M.W.F.     06/07
  !          P.BECHTOLD        E.C.M.W.F.     03/12
  
  !          PURPOSE.
  !          --------
  !          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
  !          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
  
  !          INTERFACE
  !          ---------
  
  !          THIS ROUTINE IS CALLED FROM *CUASC*.
  !          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
  !          AND UPDRAFT VALUES T,Q ETC
  !          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
  
  !          METHOD.
  !          --------
  !          TURBULENT ENTRAINMENT IS SIMULATED BY A CONSTANT
  !          MULTIPLIED BY A VERTICAL SCALING FUNCTION
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KK*           CURRENT LEVEL
  !    *KCBOT*        CLOUD BASE LEVEL
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PQSEN*        SATURATION SPEC. HUMIDITY                   KG/KG
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS          PA
  !    *PGEOH*        PROVISIONAL GEOPOTENTIAL ON HALF LEVELS      PA
  !    *PMFU*         MASSFLUX IN UPDRAFTS                        KG/(M2*S)
  !    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PDMFEN*       ENTRAINMENT RATE                            KG/(M2*S)
  !    *PDMFDE*       DETRAINMENT RATE                            KG/(M2*S)
  
  !          EXTERNALS
  !          ---------
  !          NONE
  !
  !
  !    MODIFICATIONS.
  !    --------------
  !
  !    M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
  !    M. Leutbecher (Oct 2020) SPP abstraction
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUENTR_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOECUMF, ONLY: TECUMF
  USE SPP_MOD, ONLY: TSPP_CONFIG
  USE SPP_GEN_MOD, ONLY: SPP_PERT
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TSPP_CONFIG), INTENT(IN) :: YDSPP_CONFIG
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KK
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(KLON)
  LOGICAL, INTENT(IN) :: LDCUM(KLON)
  LOGICAL, INTENT(IN) :: LDWORK
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFEN(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFDE(KLON)
  
  LOGICAL :: LLO1
  
  INTEGER(KIND=JPIM) :: JL
  
  ! A bunch of SPP variables that are saved between calls for efficiency reasons
  LOGICAL :: LLPERT_DETRPEN  ! SPP perturbation on?
  INTEGER(KIND=JPIM) :: IPDETRPEN  ! SPP random field pointer
  INTEGER(KIND=JPIM) :: IPN  ! SPP perturbation pointer
  TYPE(SPP_PERT) :: PN1  ! SPP pertn. config. for RTAU
  
  REAL(KIND=JPRB) :: ZDZ
  REAL(KIND=JPRB) :: ZENTR
  REAL(KIND=JPRB) :: ZMF
  REAL(KIND=JPRB) :: ZRG
  REAL(KIND=JPRB) :: ZXDETRPEN
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !----------------------------------------------------------------------
  
  !*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
  !                  -------------------------------------------
  
  IF (LDWORK) THEN
    
    ZRG = 1.0_JPRB / YDCST%RG
    
    ! prepare SPP perturbations (just once)
    LLPERT_DETRPEN = .false.
    
    PDMFEN(JLON) = 0.0_JPRB
    PDMFDE(JLON) = 0.0_JPRB
    ZENTR = 0.0
    
    !*    1.1          SPECIFY ENTRAINMENT RATES
    !                  -------------------------
    
    IF (LDCUM(JLON)) THEN
      ZDZ = (PGEOH(JLON, KK) - PGEOH(JLON, KK + 1))*ZRG
      ZMF = PMFU(JLON, KK + 1)*ZDZ
      LLO1 = KK < KCBOT(JLON)
      IF (LLO1) THEN
        ZXDETRPEN = YDECUMF%DETRPEN
        PDMFEN(JLON) = ZENTR*ZMF
        PDMFDE(JLON) = ZXDETRPEN*ZMF
      END IF
    END IF
    
  END IF
  
END SUBROUTINE CUENTR_OPENACC
