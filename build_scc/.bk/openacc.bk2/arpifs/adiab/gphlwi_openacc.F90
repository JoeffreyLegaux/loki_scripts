SUBROUTINE GPHLWI_OPENACC (YDDIMV, KPROMA, KST, KEND, PLNPR, PALPH, PWW, LDVERINT, YDSTACK)
  
  !**** *GPHLWI* - to half-levels interpolation weights
  
  !     Purpose.
  !     --------
  !           Compute weights for interpolation of winds to
  !           half levels (in non-hydrostatic dynamics)
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPHLWI()
  
  !        Explicit arguments :
  !        --------------------
  !         INPUT:
  !          KPROMA  - dimensioning.
  !          KSTART  - start of work.
  !          KPROF   - depth of work.
  !          PLNPR   - "delta" (log(prehyd) depth) at full levels.
  !          PALPH   - "alpha" at full levels.
  
  !         OUTPUT:
  !          PWW     - vertical weight.
  
  !        Implicit arguments : none.
  !        --------------------
  
  !     Method.
  !     -------
  !        Interpolation from full-levels using logaritmic pressure profile
  !        Then modify values on the top and bottom using boundary condition
  !        (at the bottom free-slip condition)
  !        Store the weights of vertical interpolation
  
  !     Externals.
  !     ----------
  !     Reference.
  !     ----------
  !        ARPEGE/ALADIN documentation
  
  !     Author.
  !     -------
  !        Radmila Bubnova,  CNRM/GMAP/EXT
  
  !     Modifications.
  !     --------------
  !        Original : November 1997
  !        Modified 02-07-02 by C. Fischer : remove pww5 computation + intents
  !        M.Hamrud      01-Oct-2003 CY28 Cleaning
  !        K.Yessad      20-Mar-2006 Cleaning + optimisation
  !        K. Yessad (Dec 2008): remove dummy CDLOCK
  !     ----------------------------------------------------------------------
  
!$acc routine( GPHLWI_OPENACC )
  
  USE YOMDIMV, ONLY: TDIMV
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TDIMV), INTENT(IN) :: YDDIMV
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PLNPR(KPROMA, YDDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PALPH(KPROMA, YDDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PWW(KPROMA, YDDIMV%NFLEVG)
  LOGICAL, OPTIONAL, INTENT(IN) :: LDVERINT
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  LOGICAL :: LL_VERINT
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !     ------------------------------------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  IF (PRESENT(LDVERINT)) THEN
    LL_VERINT = LDVERINT
  ELSE
    LL_VERINT = .false.
  END IF
  
  IF (.not.LL_VERINT) THEN
    DO JLEV=1,YDDIMV%NFLEVG - 1
      PWW(JLON, JLEV) =  &
      & (PLNPR(JLON, JLEV + 1) - PALPH(JLON, JLEV + 1)) / (PLNPR(JLON, JLEV + 1) - PALPH(JLON, JLEV + 1) + PALPH(JLON, JLEV))
    END DO
  ELSE
    DO JLEV=1,YDDIMV%NFLEVG - 1
      PWW(JLON, JLEV) = 0.5_JPRB
    END DO
  END IF
  !     ------------------------------------------------------------------
  
  !* Bottom half level:
  PWW(JLON, YDDIMV%NFLEVG) = 0.0_JPRB
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE GPHLWI_OPENACC
