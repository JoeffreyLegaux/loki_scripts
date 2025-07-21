SUBROUTINE CUADJTQS_OPENACC (YDTHF, YDCST, KIDIA, KFDIA, KLON, KLEV, KK, PSP, PT, PQ, LDFLAG, KCALL, YDSTACK)
  
  !**   *CUADJTQS* - SIMPLIFIED VERSION OF MOIST ADJUSTMENT
  
  !     PURPOSE.
  !     --------
  !     TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
  
  !     INTERFACE
  !     ---------
  !     THIS ROUTINE IS CALLED FROM SUBROUTINES:
  
  !       *COND*
  !       *CUBMADJ*
  !       *CUBMD*
  !       *CONDAD*
  !       *CUBMADJAD*
  !       *CUBMDAD*
  
  !     INPUT ARE UNADJUSTED T AND Q VALUES,
  !     IT RETURNS ADJUSTED VALUES OF T AND Q
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KK*           LEVEL
  !    *KCALL*        DEFINES CALCULATION AS
  !                      KCALL=0  ENV. T AND QS IN*CUINI*
  !                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
  !                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
  
  !     INPUT PARAMETERS (LOGICAL):
  
  !    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
  
  !     INPUT PARAMETERS (REAL):
  
  !    *PSP*          PRESSURE                                        PA
  
  !     UPDATED PARAMETERS (REAL):
  
  !    *PT*           TEMPERATURE                                     K
  !    *PQ*           SPECIFIC HUMIDITY                             KG/KG
  
  !     AUTHOR.
  !     -------
  !      J.F. MAHFOUF      ECMWF
  
  !     MODIFICATIONS.
  !     --------------
  !      M.Hamrud     01-Oct-2003 CY28 Cleaning
  !      20180303 : Gabor: Just a comment line to force recompilation due to
  !                        compiler wrapper optimation exception liat change
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  
  !----------------------------------------------------------------------
  
!$acc routine( CUADJTQS_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KK
  REAL(KIND=JPRB), INTENT(IN) :: PSP(KLON)
  REAL(KIND=JPRB), INTENT(INOUT) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQ(KLON, KLEV)
  LOGICAL, INTENT(IN) :: LDFLAG(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCALL
  REAL(KIND=JPRB) :: Z3ES
  REAL(KIND=JPRB) :: Z4ES
  REAL(KIND=JPRB) :: Z5ALCP
  REAL(KIND=JPRB) :: ZALDCP
  
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZQMAX
  REAL(KIND=JPRB) :: ZQP
  REAL(KIND=JPRB) :: ZCOND
  REAL(KIND=JPRB) :: ZCOND1
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZCOR
  REAL(KIND=JPRB) :: ZQSAT
  REAL(KIND=JPRB) :: ZFOEEW
  REAL(KIND=JPRB) :: Z2S
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  !DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  !----------------------------------------------------------------------
  
  !     1.           DEFINE CONSTANTS
  !                  ----------------
  
  
  ZQMAX = 0.5_JPRB
  
  !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
  !                  -----------------------------------------------------
  
  !*    ICE-WATER THERMODYNAMICAL FUNCTIONS
  
  IF (PT(JLON, KK) > YDCST%RTT) THEN
    Z3ES = YDTHF%R3LES
    Z4ES = YDTHF%R4LES
    Z5ALCP = YDTHF%R5ALVCP
    ZALDCP = YDTHF%RALVDCP
  ELSE
    Z3ES = YDTHF%R3IES
    Z4ES = YDTHF%R4IES
    Z5ALCP = YDTHF%R5ALSCP
    ZALDCP = YDTHF%RALSDCP
  END IF
  
  IF (KCALL == 1) THEN
    
    !DIR$    IVDEP
    !OCL NOVREC
    IF (LDFLAG(JLON)) THEN
      ZQP = 1.0_JPRB / PSP(JLON)
      ZTARG = PT(JLON, KK)
      ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
      ZQSAT = ZQP*ZFOEEW
      IF (ZQSAT > ZQMAX) THEN
        ZQSAT = ZQMAX
      END IF
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      Z2S = Z5ALCP / (ZTARG - Z4ES)**2
      ZCOND = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
      ZCOND = MAX(ZCOND, 0.0_JPRB)
      !     IF(ZCOND /= _ZERO_) THEN
      PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND
      PQ(JLON, KK) = PQ(JLON, KK) - ZCOND
      ZTARG = PT(JLON, KK)
      ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
      ZQSAT = ZQP*ZFOEEW
      IF (ZQSAT > ZQMAX) THEN
        ZQSAT = ZQMAX
      END IF
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      Z2S = Z5ALCP / (ZTARG - Z4ES)**2
      ZCOND1 = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
      IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
      PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND1
      PQ(JLON, KK) = PQ(JLON, KK) - ZCOND1
      !     ENDIF
    END IF
    
  END IF
  
  IF (KCALL == 2) THEN
    
    !DIR$    IVDEP
    !OCL NOVREC
    IF (LDFLAG(JLON)) THEN
      ZQP = 1.0_JPRB / PSP(JLON)
      ZTARG = PT(JLON, KK)
      ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
      ZQSAT = ZQP*ZFOEEW
      IF (ZQSAT > ZQMAX) THEN
        ZQSAT = ZQMAX
      END IF
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      Z2S = Z5ALCP / (ZTARG - Z4ES)**2
      ZCOND = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
      ZCOND = MIN(ZCOND, 0.0_JPRB)
      !     IF(ZCOND /= _ZERO_) THEN
      PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND
      PQ(JLON, KK) = PQ(JLON, KK) - ZCOND
      ZTARG = PT(JLON, KK)
      ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
      ZQSAT = ZQP*ZFOEEW
      IF (ZQSAT > ZQMAX) THEN
        ZQSAT = ZQMAX
      END IF
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      Z2S = Z5ALCP / (ZTARG - Z4ES)**2
      ZCOND1 = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
      IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
      PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND1
      PQ(JLON, KK) = PQ(JLON, KK) - ZCOND1
      !     ENDIF
    END IF
    
  END IF
  
  IF (KCALL == 0) THEN
    
    !DIR$    IVDEP
    !OCL NOVREC
    ZQP = 1.0_JPRB / PSP(JLON)
    ZTARG = PT(JLON, KK)
    ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
    ZQSAT = ZQP*ZFOEEW
    IF (ZQSAT > ZQMAX) THEN
      ZQSAT = ZQMAX
    END IF
    ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
    ZQSAT = ZQSAT*ZCOR
    Z2S = Z5ALCP / (ZTARG - Z4ES)**2
    ZCOND1 = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
    PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND1
    PQ(JLON, KK) = PQ(JLON, KK) - ZCOND1
    ZTARG = PT(JLON, KK)
    ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
    ZQSAT = ZQP*ZFOEEW
    IF (ZQSAT > ZQMAX) THEN
      ZQSAT = ZQMAX
    END IF
    ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
    ZQSAT = ZQSAT*ZCOR
    Z2S = Z5ALCP / (ZTARG - Z4ES)**2
    ZCOND1 = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
    PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND1
    PQ(JLON, KK) = PQ(JLON, KK) - ZCOND1
    
  END IF
  
  IF (KCALL == 4) THEN
    
    !DIR$    IVDEP
    !OCL NOVREC
    ZQP = 1.0_JPRB / PSP(JLON)
    ZTARG = PT(JLON, KK)
    ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
    ZQSAT = ZQP*ZFOEEW
    IF (ZQSAT > ZQMAX) THEN
      ZQSAT = ZQMAX
    END IF
    ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
    ZQSAT = ZQSAT*ZCOR
    Z2S = Z5ALCP / (ZTARG - Z4ES)**2
    ZCOND = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
    PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND
    PQ(JLON, KK) = PQ(JLON, KK) - ZCOND
    ZTARG = PT(JLON, KK)
    ZFOEEW = YDTHF%R2ES*EXP(Z3ES*(ZTARG - YDCST%RTT) / (ZTARG - Z4ES))
    ZQSAT = ZQP*ZFOEEW
    IF (ZQSAT > ZQMAX) THEN
      ZQSAT = ZQMAX
    END IF
    ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
    ZQSAT = ZQSAT*ZCOR
    Z2S = Z5ALCP / (ZTARG - Z4ES)**2
    ZCOND1 = (PQ(JLON, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
    PT(JLON, KK) = PT(JLON, KK) + ZALDCP*ZCOND1
    PQ(JLON, KK) = PQ(JLON, KK) - ZCOND1
    
  END IF
  
END SUBROUTINE CUADJTQS_OPENACC
