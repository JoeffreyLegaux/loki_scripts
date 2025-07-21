SUBROUTINE GPUVS_OPENACC (KFLEV, KPROMA, KST, KEND, LDER, PUF, PVF, PUS, PVS, PDIVF, PVORF, PUFL, PVFL, PUS_L, PVS_L, PUS_M,  &
& PVS_M, YDSTACK)
  ! --- INPUT -----------------------------------------------------------------
  ! --- OUTPUT ----------------------------------------------------------------
  ! --- OPTIONAL INPUT --------------------------------------------------------
  ! --- OPTIONAL OUTPUT -------------------------------------------------------
  
  ! GPUVS - Diagnoses "V_surf" and "grad(V_surf)".
  
  ! Purpose
  ! -------
  !   Diagnoses "V_surf" and "grad(V_surf)" (surface wind).
  !   The current assumptions which are done are:
  !   * V_surf=V(l=KFLEV)
  !   * grad(V_surf)=grad(V(l=KFLEV))
  !   The consistency must be kept with the content of routine GPHLUV.
  
  ! Interface
  ! ---------
  !   * INPUT:
  !   KFLEV   - number of levels.
  !   KPROMA  - horizontal dimension.
  !   KST  - start of work.
  !   KEND    - end of work.
  !   LDER    - T: treatment of the derivatives.
  !   PUF     - upper air U-wind at full levels.
  !   PVF     - upper air V-wind at full levels.
  
  !   * OUTPUT:
  !   PUS     - surface U wind.
  !   PVS     - surface V wind.
  
  !   * OPTIONAL INPUT:
  !   PDIVF   - upper air divergence at full levels.
  !   PVORF   - upper air vorticity at full levels.
  !   PUFL    - upper air zonal component of grad(U-wind) at full levels.
  !   PVFL    - upper air zonal component of grad(V-wind) at full levels.
  
  !   * OPTIONAL OUTPUT:
  !   PUS_L   - zonal component of grad(surface U wind).
  !   PVS_L   - zonal component of grad(surface V wind).
  !   PUS_M   - meridian component of grad(surface U wind).
  !   PVS_M   - meridian component of grad(surface V wind).
  
  ! Externals
  ! ---------
  
  ! Method
  ! ------
  
  ! Reference
  ! ---------
  
  ! Author
  ! ------
  !   K. Yessad, Dec 2004 (after GNHPDVD, GNHGRP, GPHLUV)
  
  ! Modifications
  ! -------------
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !------------------------------------------------------------------
  
!$acc routine( GPUVS_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  ! -----------------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  LOGICAL, INTENT(IN) :: LDER
  REAL(KIND=JPRB), INTENT(IN) :: PUF(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVF(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PUS(KPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PVS(KPROMA)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PDIVF(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PVORF(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PUFL(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PVFL(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PUS_L(KPROMA)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PVS_L(KPROMA)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PUS_M(KPROMA)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PVS_M(KPROMA)
  
  ! -----------------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JROF
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  ! -----------------------------------------------------------------------------
  
#include "abor1.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  
  ! -----------------------------------------------------------------------------
  
  
  ! -----------------------------------------------------------------------------
  
  !*      1. COMPUTES "V_surf" and "grad(V_surf)".
  !       ----------------------------------------
  
  PUS(JLON) = PUF(JLON, KFLEV)
  PVS(JLON) = PVF(JLON, KFLEV)
  
  IF (LDER) THEN
    IF (.not.(PRESENT(PDIVF) .and. PRESENT(PVORF) .and. PRESENT(PUFL) .and. PRESENT(PVFL) .and. PRESENT(PUS_L) .and.  &
    & PRESENT(PVS_L) .and. PRESENT(PUS_M) .and. PRESENT(PVS_M))) CALL ABOR1_ACC(' GPUVS: LDER=T => PDIVF to PVS_M should be  &
    & present!')
    PUS_L(JLON) = PUFL(JLON, KFLEV)
    PVS_L(JLON) = PVFL(JLON, KFLEV)
    PUS_M(JLON) = PVFL(JLON, KFLEV) - PVORF(JLON, KFLEV)
    PVS_M(JLON) = PDIVF(JLON, KFLEV) - PUFL(JLON, KFLEV)
  END IF
  
  ! -----------------------------------------------------------------------------
  
  
END SUBROUTINE GPUVS_OPENACC
