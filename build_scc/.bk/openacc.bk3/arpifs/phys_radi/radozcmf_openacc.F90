SUBROUTINE RADOZCMF_OPENACC (YDCST, YDEOZOC, KIDIA, KFDIA, KLON, KLEV, PAPRS, PGEMU, POZON, YDSTACK)
  
  !**** *RADOZCMF* - COMPUTES DISTRIBUTION OF OZONE FROM CLIMATOLOGY
  
  !     PURPOSE.
  !     --------
  
  !**   INTERFACE.
  !     ----------
  
  !        EXPLICIT ARGUMENTS :
  !        --------------------
  
  !        IMPLICIT ARGUMENTS :   NONE
  !        --------------------
  
  !     METHOD.
  !     -------
  
  !     EXTERNALS.
  !     ----------
  
  !          NONE
  
  !     REFERENCE.
  !     ----------
  
  !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"
  
  !     AUTHOR.
  !     -------
  !     F. BOUYSSEL  Meteo-France    2008/10/01
  !     Adaptation of radozc to work on rotated and stretched geometry
  
  !     MODIFICATIONS.
  !     --------------
  !     F. BOUYSSEL 2010-01-19   Optimisation
  
  !-----------------------------------------------------------------------
  
!$acc routine( RADOZCMF_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOEOZOC, ONLY: TEOZOC
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TEOZOC), INTENT(IN) :: YDEOZOC
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  REAL(KIND=JPRB), INTENT(IN) :: PAPRS(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PGEMU(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: POZON(KLON, KLEV)
  !     -----------------------------------------------------------------
  
  !*       0.1   ARGUMENTS.
  !              ----------
  
  !     -----------------------------------------------------------------
  
  !*       0.2   LOCAL ARRAYS.
  !              -------------
  
  temp (REAL (KIND=JPRB), ZOZLT, (KLON, 0:35))
  temp (REAL (KIND=JPRB), ZOZON, (KLON, KLEV + 1))
  temp (REAL (KIND=JPRB), ZRRR, (KLON, 0:34))
  REAL(KIND=JPRB) :: ZSILAT
  
  INTEGER(KIND=JPIM) :: INLA
  INTEGER(KIND=JPIM) :: JC
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: I1
  
  REAL(KIND=JPRB) :: Z1
  REAL(KIND=JPRB) :: Z2
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZOZLT) == 8) THEN
    alloc8 (ZOZLT)
  ELSE
    IF (KIND (ZOZLT) == 4) THEN
      alloc4 (ZOZLT)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZOZON) == 8) THEN
    alloc8 (ZOZON)
  ELSE
    IF (KIND (ZOZON) == 4) THEN
      alloc4 (ZOZON)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZRRR) == 8) THEN
    alloc8 (ZRRR)
  ELSE
    IF (KIND (ZRRR) == 4) THEN
      alloc4 (ZRRR)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KIDIA
  
  !     ------------------------------------------------------------------
  !     ------------------------------------------------------------------
  !     ------------------------------------------------------------------
  
  !*         1.     LATITUDE INDEX WITHIN OZONE CLIMATOLOGY
  !                 ---------------------------------------
  
  
  INLA = 1 + INT(18.0_JPRB / YDCST%RPI*ASIN(PGEMU(JLON)) + 9.0_JPRB)
  ZSILAT = (PGEMU(JLON) - YDEOZOC%RSINC(INLA)) / (YDEOZOC%RSINC(INLA + 1) - YDEOZOC%RSINC(INLA))
  
  !     ------------------------------------------------------------------
  
  !*         2.     LATITUDE INTERPOLATED FIELD
  !                 ---------------------------
  
  DO JC=0,35
    Z1 = REAL(INLA)
    Z2 = MIN(ABS((Z1 - 1._JPRB)*(Z1 - 18._JPRB)), 1.0_JPRB)
    I1 = INLA
    ZOZLT(JLON, JC) =  &
    & (1.0_JPRB - Z2)*YDEOZOC%ROZT(I1, JC) + Z2*(YDEOZOC%ROZT(I1, JC) + ZSILAT*(YDEOZOC%ROZT(I1 + 1, JC) - YDEOZOC%ROZT(I1, JC)))
  END DO
  
  !     ------------------------------------------------------------------
  
  !*         3.     VERTICAL INTERPOLATION
  !                 ----------------------
  
  DO JC=0,34
    ZRRR(JLON, JC) = (1.0_JPRB / (YDEOZOC%RPROC(JC) - YDEOZOC%RPROC(JC + 1)))*(ZOZLT(JLON, JC) - ZOZLT(JLON, JC + 1))
  END DO
  
  DO JLEV=1,KLEV + 1
    ZOZON(JLON, JLEV) = 0.0_JPRB
  END DO
  
  DO JC=0,34
    DO JLEV=1,KLEV + 1
      ZOZON(JLON, JLEV) = ZOZON(JLON, JLEV) + (ZOZLT(JLON, JC + 1) + (PAPRS(JLON, JLEV) - YDEOZOC%RPROC(JC + 1))*ZRRR(JLON, JC)) &
      & *MAX(0.0_JPRB, SIGN(1.0_JPRB, PAPRS(JLON, JLEV) - YDEOZOC%RPROC(JC)))*MAX(0.0_JPRB, SIGN(1.0_JPRB, YDEOZOC%RPROC(JC + 1)  &
      & - PAPRS(JLON, JLEV)))
    END DO
  END DO
  
  DO JLEV=1,KLEV + 1
    ZOZON(JLON, JLEV) = ZOZON(JLON, JLEV) + ZOZLT(JLON, 35)*MAX(0.0_JPRB, SIGN(1.0_JPRB, PAPRS(JLON, JLEV) - YDEOZOC%RPROC(35)))
  END DO
  
  ! INTEGRATION IN THE VERTICAL:
  
  DO JLEV=1,KLEV
    POZON(JLON, JLEV) = (PAPRS(JLON, JLEV + 1) - PAPRS(JLON, JLEV))*(ZOZON(JLON, JLEV) + ZOZON(JLON, JLEV + 1))*0.5_JPRB
  END DO
  
  !     -----------------------------------------------------------
  
END SUBROUTINE RADOZCMF_OPENACC
