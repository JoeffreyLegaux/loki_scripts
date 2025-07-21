SUBROUTINE CUININ_OPENACC (YDCST, YDTHF, YDEPHLI, YDECUMF, KIDIA, KFDIA, KLON, KLEV, PTEN, PQEN, PQSEN, PUEN, PVEN, PGEO, PAPH,  &
& KLAB, PTENH, PQENH, PQSENH, PGEOH, PTU, PQU, PTD, PQD, PUU, PVU, PUD, PVD, PLU, YDSTACK)
  
  !          PURPOSE
  !          -------
  
  !          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
  !          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
  !          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
  !          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
  
  !          INTERFACE
  !          ---------
  !          THIS ROUTINE IS CALLED FROM *CUMASTR*.
  
  !          METHOD.
  !          --------
  !          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
  !    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
  !    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
  !    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
  !    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
  !    *PGEO*         GEOPOTENTIAL                                  M2/S2
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
  
  !    OUTPUT PARAMETERS (INTEGER):
  
  !    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
  !                        KLAB=2 FOR CONDENSATION LEVEL
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS         K
  !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS    KG/KG
  !    *PQSENH*       ENV. SPEC. SATURATION HUMIDITY (T+1)
  !                   ON HALF LEVELS                              KG/KG
  !    *PTU*          TEMPERATURE IN UPDRAFTS                       K
  !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                  KG/KG
  !    *PTD*          TEMPERATURE IN DOWNDRAFTS                     K
  !    *PQU*          SPEC. HUMIDITY IN DOWNDRAFTS                KG/KG
  !    *PUU*          U-VELOCITY IN UPDRAFTS                       M/S
  !    *PVU*          V-VELOCITY IN UPDRAFTS                       M/S
  !    *PUD*          U-VELOCITY IN DOWNDRAFTS                     M/S
  !    *PVD*          V-VELOCITY IN DOWNDRAFTS                     M/S
  !    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
  
  !          EXTERNALS
  !          ---------
  !          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.     12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      05-02-11 : Optimisation (NJKT2) P. BECHTOLD
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUININ_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOECUMF, ONLY: TECUMF
  USE YOEPHLI, ONLY: TEPHLI
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(IN) :: PTEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KLAB(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQSENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(OUT) :: PTU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PTD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PUU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PVU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PUD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PVD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB) :: ZPH
  LOGICAL :: LLFLAG
  
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZORCPD  !, ZZS
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cuadjtq.intfb.h"
#include "cuadjtqs.intfb.h"
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JL = KIDIA
  
  !----------------------------------------------------------------------
  
  !*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
  !*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
  !*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
  !                  ----------------------------------------------
  
  ZORCPD = 1._JPRB / YDCST%RCPD
  DO JK=2,KLEV
    PTENH(JL, JK) =  &
    & (MAX(YDCST%RCPD*PTEN(JL, JK - 1) + PGEO(JL, JK - 1), YDCST%RCPD*PTEN(JL, JK) + PGEO(JL, JK)) - PGEOH(JL, JK))*ZORCPD
    PQENH(JL, JK) = PQEN(JL, JK - 1)
    PQSENH(JL, JK) = PQSEN(JL, JK - 1)
    ZPH = PAPH(JL, JK)
    LLFLAG = .true.
    
    !orig   IF(JK.GE.KLEV-1) GO TO 130
    IF (JK >= KLEV - 1 .or. JK < YDECUMF%NJKT2) CYCLE
    IK = JK
    IF (YDEPHLI%LPHYLIN) THEN
      CALL CUADJTQS_OPENACC(YDTHF, YDCST, KIDIA, KFDIA, KLON, KLEV, IK, ZPH, PTENH, PQSENH, LLFLAG, 0, YDSTACK=YLSTACK)
    ELSE
      CALL CUADJTQ_OPENACC(YDTHF, YDCST, YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK, ZPH, PTENH, PQSENH, LLFLAG, 3, YDSTACK=YLSTACK)
    END IF
    
    PQENH(JL, JK) = MIN(PQEN(JL, JK - 1), PQSEN(JL, JK - 1)) + (PQSENH(JL, JK) - PQSEN(JL, JK - 1))
    PQENH(JL, JK) = MAX(PQENH(JL, JK), 0.0_JPRB)
    !orig  130   continue
  END DO
  
  PTENH(JL, KLEV) = (YDCST%RCPD*PTEN(JL, KLEV) + PGEO(JL, KLEV) - PGEOH(JL, KLEV))*ZORCPD
  PQENH(JL, KLEV) = PQEN(JL, KLEV)
  PTENH(JL, 1) = PTEN(JL, 1)
  PQENH(JL, 1) = PQEN(JL, 1)
  
  !DO JK=KLEV-1,2,-1
  !  DO JL=KIDIA,KFDIA
  !    ZZS=MAX(RCPD*PTENH(JL,JK)+PGEOH(JL,JK),&
  !     & RCPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))
  !   PTENH(JL,JK)=(ZZS-PGEOH(JL,JK))*ZORCPD
  !  ENDDO
  !ENDDO
  
  !-----------------------------------------------------------------------
  
  !*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
  !*                 ---------------------------------------------
  
  DO JK=1,KLEV
    IK = JK - 1
    IF (JK == 1) IK = 1
    PTU(JL, JK) = PTENH(JL, JK)
    PTD(JL, JK) = PTENH(JL, JK)
    PQU(JL, JK) = PQENH(JL, JK)
    PQD(JL, JK) = PQENH(JL, JK)
    PLU(JL, JK) = 0.0_JPRB
    PUU(JL, JK) = PUEN(JL, IK)
    PUD(JL, JK) = PUEN(JL, IK)
    PVU(JL, JK) = PVEN(JL, IK)
    PVD(JL, JK) = PVEN(JL, IK)
    KLAB(JL, JK) = 0
  END DO
  
END SUBROUTINE CUININ_OPENACC
