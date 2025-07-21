SUBROUTINE CUCCDIA_OPENACC (YDERAD, YDEPHLI, YDEPHY, KIDIA, KFDIA, KLON, KLEV, KSTEP, KCBOT, KCTOP, LDCUM, PQU, PLU, PMFU,  &
& PRAIN, PARPRC, KTOPC, KBASEC, YDSTACK)
  
  !**** *CUCCDIA*- UPDATES PRECIPITAION, CLOUD BASE AND CLOUD TOP
  !                FOR DIAGNOSTIC SCHEME FOR CONVECTIVE CLOUDS
  
  !          M.TIEDTKE         E.C.M.W.F.    12/89
  
  !**   INTERFACE.
  !     ----------
  
  !          *CUCCDIA* IS CALLED FROM *CUCALL*
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KSTEP*        CURRENT TIME STEP INDEX
  !    *KCBOT*        CLOUD BASE LEVEL
  !    *KCTOP*        CLOUD TOP LEVEL
  
  !     INPUT PARAMETERS (LOGICAL)
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  
  !     INPUT PARAMETERS (REAL)
  
  !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
  !    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
  !    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PARPRC*       ACCUMULATED PRECIPITATION AMMOUNT             KG/(M2*S)
  !                   FOR RADIATION CALCULATION
  !    *KTOPC*        CONVECTIVE CLOUD TOP LEVEL FOR RADIATION
  !    *KBASEC*       CONVECTIVE CLOUD BASE LEVEL FOR RADIATION
  
  !-----------------------------------------------------------------------
  
!$acc routine( CUCCDIA_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOEPHY, ONLY: TEPHY
  USE YOERAD, ONLY: TERAD
  USE YOEPHLI, ONLY: TEPHLI
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  TYPE(TEPHY), INTENT(IN) :: YDEPHY
  TYPE(TERAD), INTENT(IN) :: YDERAD
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP(KLON)
  LOGICAL, INTENT(IN) :: LDCUM(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRAIN(KLON)
  REAL(KIND=JPRB), INTENT(INOUT) :: PARPRC(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KTOPC(KLON)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KBASEC(KLON)
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: IKB
  INTEGER(KIND=JPIM) :: IKT
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZDMFQ
  REAL(KIND=JPRB) :: ZNORMR
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  
  !---------------------------------------------------------------------
  
  !*    1.0          STORE CLOUD PARAMETERS FOR RADIATION CALCULATION
  !                  -----------------------------------------------
  
  !********************************
  IF (.not.YDEPHLI%LPHYLIN .and. YDERAD%NRADFR /= 1) THEN
    !********************************
    IF (YDEPHY%LECUMF .and. YDEPHY%LERADI .and. MOD(KSTEP + 1, YDERAD%NRADFR) /= 0) THEN
      ZNORMR = 1.0_JPRB / REAL(MAX(YDERAD%NRADFR - 1, 1), kind=JPRB)
      !DIR$ IVDEP
      !OCL NOVREC
      IF (LDCUM(JLON)) THEN
        KBASEC(JLON) = MAX(KBASEC(JLON), KCBOT(JLON))
        KTOPC(JLON) = MIN(KTOPC(JLON), KCTOP(JLON))
        IF (KTOPC(JLON) == 1) KTOPC(JLON) = KCTOP(JLON)
        IKB = KCBOT(JLON)
        IKT = KCTOP(JLON)
        ZDMFQ = PMFU(JLON, IKB)*MAX(PQU(JLON, IKB) + PLU(JLON, IKB) - PQU(JLON, IKT), PQU(JLON, IKB)*0.05_JPRB)
        PARPRC(JLON) = PARPRC(JLON) + MAX(ZDMFQ, PRAIN(JLON))*ZNORMR
      END IF
    END IF
    !********************************
  ELSE
    !********************************
    IF (YDEPHY%LECUMF) THEN
      !DIR$ IVDEP
      !OCL NOVREC
      IF (LDCUM(JLON)) THEN
        KBASEC(JLON) = KCBOT(JLON)
        KTOPC(JLON) = KCTOP(JLON)
        IF (KTOPC(JLON) == 1) KTOPC(JLON) = KCTOP(JLON)
        IKB = KCBOT(JLON)
        IKT = KCTOP(JLON)
        ZDMFQ = PMFU(JLON, IKB)*MAX(PQU(JLON, IKB) + PLU(JLON, IKB) - PQU(JLON, IKT), PQU(JLON, IKB)*0.05_JPRB)
        PARPRC(JLON) = MAX(ZDMFQ, PRAIN(JLON))
      END IF
    END IF
    !********************************
  END IF
  !********************************
  
END SUBROUTINE CUCCDIA_OPENACC
