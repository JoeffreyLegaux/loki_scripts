SUBROUTINE GP_TNDLAGADIAB_UV_OPENACC (LDRPLANE, YDGEOMETRY, YDEPHY, YDDYN, KST, KEND, PRCORI, PGNORDL, PGNORDM, PSGRTL, PSGRTM,  &
& PU, PV, PTNDU, PTNDV, PTNDU_NOC, PTNDV_NOC, YDSTACK)
  
  !**** *GP_TNDLAGADIAB_UV*   Compute adiabatic Lagrangian tendency of horizontal wind.
  
  !     Purpose.
  !     --------
  !          Compute adiabatic Lagrangian tendency of horizontal wind, with the following assumptions:
  !          - explicit representation of Coriolis term.
  !          - no curvature term.
  !          - Rayleigh friction taken into account.
  !          - pressure gradient term taken into account.
  
  !**   Interface.
  !     ----------
  !        *CALL* *GP_TNDLAGADIAB_UV(..)
  
  !        Explicit arguments :
  !        --------------------
  
  !        INPUT:
  !          KST       - first element of work.
  !          KPROF     - depth of work.
  !          PRCORI    - Coriolis parameter "f = 2 Omega sin(theta)".
  !          PGNORDL   - zonal component ("Gnordl") of the unit vector
  !                      directed towards the true North pole.
  !          PGNORDM   - meridian component ("Gnordm") of the unit vector
  !                      directed towards the true North pole.
  !          PSGRTL    - zonal component of the pressure force grad.
  !          PSGRTM    - merid component of the pressure force grad.
  !          PGMV      - GMV variables at t-dt and t.
  
  !        OUTPUT:
  !          PTNDU     - Tendency for U-wind.
  !          PTNDV     - Tendency for V-wind.
  !          PTNDU_NOC - Tendency for U-wind without Coriolis term.
  !          PTNDV_NOC - Tendency for V-wind without Coriolis term.
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  !           none
  
  !     Reference.
  !     ----------
  !             Arpege documentation about model equations.
  
  !     Author.
  !     -------
  !        K. YESSAD (METEO FRANCE/CNRM/GMAP)
  !         after some code present in LAVENT, CPEULDYN, LATTEX.
  !        Original : SEPT 2010.
  
  ! Modifications
  ! -------------
  !  K. Yessad (Nov 2011): more flexible Rayleigh friction, and adaptations.
  !  K. Yessad (July 2014): Move some variables.
  !  K. Yessad (Feb 2018): remove deep-layer formulations.
  !  H. Petithomme (Dec 2020): test reordering and optimisation
  !------------------------------------------------------------------------------
  ! End Modifications
  
!$acc routine( GP_TNDLAGADIAB_UV_OPENACC ) seq
  
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK
  USE PARKIND1, ONLY: JPIM, JPRB
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMDYN, ONLY: TDYN
  USE YOEPHY, ONLY: TEPHY
  
  
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN) :: LDRPLANE
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(TEPHY), INTENT(IN) :: YDEPHY
  TYPE(TDYN), INTENT(IN) :: YDDYN
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PRCORI(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PGNORDL(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PGNORDM(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PU(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PV(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PTNDU(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PTNDV(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PTNDU_NOC(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PTNDV_NOC(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  INTEGER(KIND=JPIM) :: JROF
  INTEGER(KIND=JPIM) :: JLEV
  REAL(KIND=JPRB) :: ZGNORDLM
  REAL(KIND=JPRB) :: ZGNORDL2
  REAL(KIND=JPRB) :: ZGNORDM2
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  
  
  
  IF (YDEPHY%LEGWWMS .or. .not.YDDYN%LRFRIC) THEN
    ! no Rayleigh friction:
    PTNDU_NOC(JLON, 1:YDGEOMETRY%YRDIMV%NFLEVG) = -PSGRTL(JLON, 1:YDGEOMETRY%YRDIMV%NFLEVG)
    PTNDV_NOC(JLON, 1:YDGEOMETRY%YRDIMV%NFLEVG) = -PSGRTM(JLON, 1:YDGEOMETRY%YRDIMV%NFLEVG)
  ELSE IF (YDDYN%LRFRICISOTR) THEN
    ! isotropic Rayleigh friction:
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      PTNDU_NOC(JLON, JLEV) = -YDDYN%RKRF(JLEV)*PU(JLON, JLEV) - PSGRTL(JLON, JLEV)
      PTNDV_NOC(JLON, JLEV) = -YDDYN%RKRF(JLEV)*PV(JLON, JLEV) - PSGRTM(JLON, JLEV)
    END DO
  ELSE IF (.not.LDRPLANE .and. YDGEOMETRY%YRGEM%NSTTYP == 1) THEN
    ! non-isotropic Rayleigh friction in untilted spherical geometry:
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      PTNDU_NOC(JLON, JLEV) = -YDDYN%RKRF(JLEV)*PU(JLON, JLEV) - PSGRTL(JLON, JLEV)
      PTNDV_NOC(JLON, JLEV) = -PSGRTM(JLON, JLEV)
    END DO
  ELSE
    ! non-isotropic Rayleigh friction, other cases:
    ZGNORDLM = PGNORDM(JLON)*PGNORDL(JLON)
    ZGNORDL2 = PGNORDL(JLON)**2
    ZGNORDM2 = PGNORDM(JLON)**2
    
    ! simplified version from (with rkrfu=rkrf, rkrfv=0):
    ! ptndu=(rkrfu-rkrfv)*nordm*nordl*v-(rkrfu*nordm**2+rkrfv*nordl**2)*u
    ! ptndv=(rkrfu-rkrfv)*nordm*nordl*u-(rkrfu*nordl**2+rkrfv*nordm**2)*v
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      PTNDU_NOC(JLON, JLEV) = YDDYN%RKRF(JLEV)*(ZGNORDLM*PV(JLON, JLEV) - ZGNORDM2*PU(JLON, JLEV)) - PSGRTL(JLON, JLEV)
      PTNDV_NOC(JLON, JLEV) = YDDYN%RKRF(JLEV)*(ZGNORDLM*PU(JLON, JLEV) - ZGNORDL2*PV(JLON, JLEV)) - PSGRTM(JLON, JLEV)
    END DO
  END IF
  
  DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
    PTNDU(JLON, JLEV) = PTNDU_NOC(JLON, JLEV) + PRCORI(JLON)*PV(JLON, JLEV)
    PTNDV(JLON, JLEV) = PTNDV_NOC(JLON, JLEV) - PRCORI(JLON)*PU(JLON, JLEV)
  END DO
  
END SUBROUTINE GP_TNDLAGADIAB_UV_OPENACC
