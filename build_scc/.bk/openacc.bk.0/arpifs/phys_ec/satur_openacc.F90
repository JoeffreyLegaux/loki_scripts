SUBROUTINE SATUR_OPENACC (YDTHF, YDCST, KIDIA, KFDIA, KLON, KTDIA, KLEV, LDPHYLIN, PAPRSF, PT, PQSAT, KFLAG, YDSTACK)
  
  !***
  
  ! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION
  
  !       J.F. MAHFOUF       E.C.M.W.F.     15/05/96
  
  !       Modified J. HAGUE          13/01/03 MASS Vector Functions
  
  !       PURPOSE.
  !       --------
  
  !       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
  !       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
  !       AND LIQUID WATER CONTENT
  
  !       INTERFACE
  !       ---------
  
  !       THIS ROUTINE IS CALLED FROM *CALLPAR*.
  
  !       PARAMETER     DESCRIPTION                                 UNITS
  !       ---------     -----------                                 -----
  !       INPUT PARAMETERS (INTEGER):
  
  !      *KIDIA*        START POINT
  !      *KFDIA*        END POINT
  !      *KLON*         NUMBER OF GRID POINTS PER PACKET
  !      *KTDIA*        START OF THE VERTICAL LOOP
  !      *KLEV*         NUMBER OF LEVELS
  
  !       INPUT PARAMETERS (REAL):
  
  !      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
  !      *PT*            TEMPERATURE AT T-DT                          K
  
  !       INPUT PARAMETERS (INTEGER):
  
  !      *KFLAG*         FLAG TO DETECT CALL FROM
  
  !                      CONVECTION  KFLAG=1
  !                      OTHER       KFLAG=2
  
  !       OUTPUT PARAMETER (REAL):
  
  !      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG
  
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !-------------------------------------------------------------------------
  
!$acc routine( SATUR_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  LOGICAL, INTENT(IN) :: LDPHYLIN
  REAL(KIND=JPRB), INTENT(IN) :: PAPRSF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQSAT(KLON, KLEV)
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLAG
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZCOR
  REAL(KIND=JPRB) :: ZEW
  REAL(KIND=JPRB) :: ZFOEEW
  REAL(KIND=JPRB) :: ZQMAX
  REAL(KIND=JPRB) :: ZQS
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZALFA
  REAL(KIND=JPRB) :: ZFOEEWL
  REAL(KIND=JPRB) :: ZFOEEWI
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  !DIR$ VFUNCTION EXPHF
  
#include "fcttre.func.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !----------------------------------------------------------------------
  
  !*    1.           DEFINE CONSTANTS
  !                  ----------------
  
  
  ZQMAX = 0.5_JPRB
  
  !     *
  !----------------------------------------------------------------------
  
  !     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
  !                       --------------------------------------
  
  IF (LDPHYLIN) THEN
    DO JK=KTDIA,KLEV
      ZTARG = PT(JLON, JK)
      ZALFA = FOEALFA(ZTARG)
      
      ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
      ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
      ZFOEEW = ZALFA*ZFOEEWL + (1.0_JPRB - ZALFA)*ZFOEEWI
      
      ZQS = ZFOEEW / PAPRSF(JLON, JK)
      IF (ZQS > ZQMAX) THEN
        ZQS = ZQMAX
      END IF
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQS)
      PQSAT(JLON, JK) = ZQS*ZCOR
    END DO
  ELSE
    
    DO JK=KTDIA,KLEV
      IF (KFLAG == 1) THEN
        ZEW = FOEEWMCU(PT(JLON, JK))
      ELSE
        ZEW = FOEEWM(PT(JLON, JK))
      END IF
      ZQS = ZEW / PAPRSF(JLON, JK)
      ZQS = MIN(ZQMAX, ZQS)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQS)
      PQSAT(JLON, JK) = ZQS*ZCOR
    END DO
    
  END IF
  
END SUBROUTINE SATUR_OPENACC
