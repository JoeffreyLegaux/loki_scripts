SUBROUTINE GPRT_OPENACC (LDSPRT, KPROMA, KST, KEND, KLEV, PRD, PRV, PR, PT, PTL, PTM, PQL, PQM, PRT, PRTL, PRTM, PRL, PRM,  &
& YDSTACK)
  
  !*    *GPRT*
  
  !     Purpose
  !     -------  To calculate RT and its derivates
  !     Interface
  !     ---------
  
  !     Explicit arguments
  !     ------------------
  !     Input:
  !    -------
  !              LDSPRT   : .TRUE. if PTL and PTM already contain
  !                         the derivatives of 'TV')
  !              KPROMA   : Horizontal dimension
  !              KSTART   : Start index
  !              KEND     : End index
  !              KLEV     : number of levels
  !              PRD      : Rd
  !              PRV      : Rv
  !              PR       : R
  !              PT       : T
  !              PTL ,PTM : Horizontal derivatives of T
  !              PQL, PQM : Horizontal derivatives of q
  !     Output:
  !    --------
  !              PRT       : RT
  !              PRTL,PRTM : Horizontal derivatives of RT
  
  !     Author
  !     ------
  !           J.Boutahar *MAROC-METEO*
  !     Modifications
  !     -------------
  !           Original: 97/06/06
  !     C. Fischer 02-06-27 : cdlock
  !     M.Hamrud      01-Oct-2003 CY28 Cleaning
  !     K. Yessad (Dec 2008): remove dummy CDLOCK
  !     H Petithomme (Dec 2020): tests simplification and hoisting
  !----------------------------------------------------------
  
!$acc routine( GPRT_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  !----------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN) :: LDSPRT
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  REAL(KIND=JPRB), INTENT(IN) :: PRD
  REAL(KIND=JPRB), INTENT(IN) :: PRV
  REAL(KIND=JPRB), INTENT(IN) :: PR(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTL(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTM(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQL(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQM(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PRT(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PRTL(KPROMA, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PRTM(KPROMA, KLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PRL(KPROMA, KLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PRM(KPROMA, KLEV)
  
  !----------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  REAL(KIND=JPRB) :: ZR
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !----------------------------------------------------------
  
  ZR = PRV - PRD
  
  !---------------------------------------------------------------
  !*       1. Compute RT and its derivatives
  
  DO JLEV=1,KLEV
    PRT(JLON, JLEV) = PR(JLON, JLEV)*PT(JLON, JLEV)
  END DO
  
  IF (LDSPRT) THEN
    DO JLEV=1,KLEV
      PRTL(JLON, JLEV) = PRD*PTL(JLON, JLEV)
      PRTM(JLON, JLEV) = PRD*PTM(JLON, JLEV)
    END DO
  ELSE IF (PRESENT(PRL)) THEN
    DO JLEV=1,KLEV
      PRTL(JLON, JLEV) = PRL(JLON, JLEV)*PT(JLON, JLEV) + PR(JLON, JLEV)*PTL(JLON, JLEV)
      PRTM(JLON, JLEV) = PRM(JLON, JLEV)*PT(JLON, JLEV) + PR(JLON, JLEV)*PTM(JLON, JLEV)
    END DO
  ELSE
    DO JLEV=1,KLEV
      PRTL(JLON, JLEV) = ZR*PT(JLON, JLEV)*PQL(JLON, JLEV) + PR(JLON, JLEV)*PTL(JLON, JLEV)
      PRTM(JLON, JLEV) = ZR*PT(JLON, JLEV)*PQM(JLON, JLEV) + PR(JLON, JLEV)*PTM(JLON, JLEV)
    END DO
  END IF
  
  !----------------------------------------------------------
END SUBROUTINE GPRT_OPENACC
