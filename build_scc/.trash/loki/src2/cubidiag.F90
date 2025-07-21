SUBROUTINE CUBIDIAG_OPENACC (KIDIA, KFDIA, KLON, KLEV, KCTOP, LD_LCUMASK, PA, PB, PR, PU, YDSTACK)
  
  !          P. Bechtold         E.C.M.W.F.     07/03
  
  !          PURPOSE.
  !          --------
  !          SOLVES BIDIAGONAL SYSTEM
  !          FOR IMPLICIT SOLUTION OF ADVECTION EQUATION
  
  !          INTERFACE
  !          ---------
  
  !          THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
  !          IT RETURNS UPDATED VALUE OF QUANTITY
  
  !          METHOD.
  !          --------
  !          NUMERICAL RECIPES (Cambridge Press)
  !          DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
  !          (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
  !          M  x  U  = R
  !          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C
  
  !          (  B(kctop-1)    C(kctop-1)    0          0        )
  !          (  A(kctop)      B(kctop)      C(kctop)   0        )
  !          (  0             A(jk)         B(jk)      C(jk)    )
  !          (  0             0             A(klev)    B(klev)  )
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  
  !    INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KCTOP*        CLOUD TOP LEVELS
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS
  !    *PR*           RHS VECTOR CONTAINING "CONSTANTS"
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY
  
  !          EXTERNALS
  !          ---------
  !          NONE
  
  !----------------------------------------------------------------------
  
!$acc routine( CUBIDIAG_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  !     DUMMY INTEGER
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP
  LOGICAL, INTENT(IN) :: LD_LCUMASK(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PA(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PB(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PR(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PU(KLON, KLEV)
  !     DUMMY REALS
  !     DUMMY LOGICALS
  !     LOCALS
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  REAL(KIND=JPRB) :: ZBET
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JL = KIDIA
  
  !----------------------------------------------------------------------
  
  
  PU(JL, :) = 0.0_JPRB
  
  ! Forward Substitution
  
  DO JK=2,KLEV
    IF (LD_LCUMASK(JL, JK)) THEN
      IF (JK == KCTOP - 1) THEN
        ZBET = 1.0_JPRB / (PB(JL, JK) + 1.E-35_JPRB)
        PU(JL, JK) = PR(JL, JK)*ZBET
      ELSE IF (JK > KCTOP - 1) THEN
        ZBET = 1.0_JPRB / (PB(JL, JK) + 1.E-35_JPRB)
        PU(JL, JK) = (PR(JL, JK) - PA(JL, JK)*PU(JL, JK - 1))*ZBET
      END IF
    END IF
  END DO
  
END SUBROUTINE CUBIDIAG_OPENACC
