SUBROUTINE GPHLUV_EXPL_OPENACC (YDDIMV, KPROMA, KST, KEND, PU, PV, PWWI, PUH, PVH, YDSTACK)
  
  !**** *GPHLUV_EXPL* - wind components calculation in half-levels
  
  !     Purpose.
  !     --------
  !           Compute wind components in half-levels
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPHLUV_EXPL()
  
  !        Explicit arguments :
  !        --------------------
  !        * INPUT:
  !          KPROMA  - length of work
  !          KSTART  - start of work
  !          KPROF   - end of work
  !          PU      - U-wind at full levels
  !          PV      - V-wind at full levels
  
  !        * IN/OUT:
  !          PUVH    - horizontal wind and weights at half levels
  !                    (IN) for weights, (OUT) for half level wind
  
  !        Implicit arguments : none.
  !        --------------------
  
  !     Method.
  !     -------
  !        Interpolation from full-levels using logaritmic pressure profile
  !        Then modify values on the top and bottom using boundary condition
  !        (at the bottom free-slip condition)
  !        Store also the weights of vertical interpolation
  
  !     Externals.
  !     ----------
  !     Reference.
  !     ----------
  !        ARPEGE/ALADIN documentation
  
  !     Author.
  !     -------
  !      Radmila Bubnova & Martin Janousek,  CNRM/GMAP/EXT
  !      Original : February 1996
  
  ! Modifications
  ! -------------
  !   Modified 02-07-02 by C. Fischer : intents for dummy arrays
  !   Modified 08-2002 C. Smith and K. YESSAD:
  !    - make optional the vertical averaging.
  !    - some cleanings and optimisations.
  !   01-Oct-2003 M. Hamrud  CY28 Cleaning
  !   09-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
  !   21-Jan-2005 K. Yessad  Remove useless dummy arguments.
  !   07-Mar-2007 K. Yessad  Remove LVERAVE_HLUV.
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !   K. Yessad (Dec 2008): remove dummy CDLOCK
  !   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
  ! End Modifications
  !------------------------------------------------------------------
  
!$acc routine( GPHLUV_EXPL_OPENACC )
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, JPHOOK, DR_HOOK
  
  USE YOMDIMV, ONLY: TDIMV
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TDIMV), INTENT(IN) :: YDDIMV
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PU(KPROMA, YDDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PV(KPROMA, YDDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PWWI(KPROMA, 0:YDDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PUH(KPROMA, 0:YDDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PVH(KPROMA, 0:YDDIMV%NFLEVG)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JROF
  INTEGER(KIND=JPIM) :: JLEV
  REAL(KIND=JPRB) :: ZWW
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !     ------------------------------------------------------------------
  
  
  
  !     ------------------------------------------------------------------
  
  !*    1. General case:
  
  DO JLEV=1,YDDIMV%NFLEVG - 1
    ZWW = PWWI(JLON, JLEV)
    PUH(JLON, JLEV) = ZWW*PU(JLON, JLEV) + (1.0_JPRB - ZWW)*PU(JLON, JLEV + 1)
    PVH(JLON, JLEV) = ZWW*PV(JLON, JLEV) + (1.0_JPRB - ZWW)*PV(JLON, JLEV + 1)
  END DO
  
  !     ------------------------------------------------------------------
  
  !*    2. Top.
  
  PUH(JLON, 0) = PU(JLON, 1)
  PVH(JLON, 0) = PV(JLON, 1)
  
  !     ------------------------------------------------------------------
  
  !*    3. Surface.
  
  PUH(JLON, YDDIMV%NFLEVG) = PU(JLON, YDDIMV%NFLEVG)
  PVH(JLON, YDDIMV%NFLEVG) = PV(JLON, YDDIMV%NFLEVG)
  
  !     ------------------------------------------------------------------
  
  
  
END SUBROUTINE GPHLUV_EXPL_OPENACC
