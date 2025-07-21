SUBROUTINE GPXX_OPENACC (YDGEOMETRY, KFLEV, KPROMA, KST, KEND, PHIHL, PHIHM, PHIFL, PHIFM, PLNPR, PRT, PUF, PVF, PUH, PVH, PX,  &
& PNHPPI, LDVFE, YDSTACK)
  
  ! GPXX - Diagnose NHX-term
  
  ! Purpose
  ! -------
  !   Diagnose NHX-term
  !    NHX = (pre/(RT)) grad[gz] (d V / d prehyd)
  !   It is better to rewrite NHX as follows for the discretisation:
  !    NHX = (pre/prehyd) (1/(RT)) grad[gz] ( d V / d (log(prehyd)) )
  
  !   NHX is discretised as follows, at full levels, for FD discretisation:
  !   [X]_[l] = [pre/prehyd]_[l] * [ 1/(R_[l] T_[l] delta_[l]) ] *
  !             [ grad[gz]_[lbar-1] (V_[l] - V[lbar-1])
  !             + grad[gz]_[lbar] (V_[lbar] - V_[l]) ]
  !   VFE discretisation is different.
  
  !   This routine can be used in a NHEE model: PNHPPI should be present in this case.
  
  !   This routine can also be used in a hydrostatic model or in a NHQE model:
  !   in this case the ratio (pre/prehyd) is equal to 1 and PNHPPI should not be used.
  
  ! Interface
  ! ---------
  !   * INPUT:
  !   YDGEOMETRY   : structure containing all geometry
  !   KFLEV        : number of levels.
  !   KPROMA       : length of work
  !   KSTART       : start of work
  !   KEND         : end of work
  !   PHIHL        : zonal component of "grad[gz]" at half levels.
  !   PHIHM        : meridian component of "grad[gz]" at half levels.
  !   PHIFL        : zonal component of "grad[gz]" at full levels.
  !   PHIFM        : meridian component of "grad[gz]" at full levels.
  !   PLNPR        : "delta" at full levels.
  !   PRT          : (R*Temperature) at full levels.
  !   PUF          : U-wind at full levels.
  !   PVF          : V-wind at full levels.
  !   PUH          : U-wind at half levels.
  !   PVH          : V-wind at half levels.
  
  !   * OUTPUT:
  !   PX           : NHX-term at full levels.
  
  !   * OPTIONAL INPUT:
  !   PNHPPI       : [pre/prehyd] at full levels (NHEE model).
  !   LDVFE        : T if VFE discretisation is used in this routine.
  
  ! Externals
  ! ---------
  
  ! Method
  ! ------
  
  ! Reference
  ! ---------
  
  ! Author
  ! ------
  !   06 Dec 2004 K. Yessad (after GNHX).
  
  ! Modifications
  ! -------------
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !   K. Yessad (March 2009): correct false comments for LRWSDLG=T
  !   J. Vivoda (Oct 2013): new options for VFE-NH
  !   K. Yessad (June 2017): Introduce NHQE model.
  !   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
  !   K. Yessad (Feb 2018): remove deep-layer formulations.
  !   J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
  !   R. El Khatib 31-Jan-2022 Optional VFE
  !   C. Wastl and C. Wittmann (Feb 2022): Add updraft helicity
  ! End Modifications
  !---------------------------------------------------------------------
  
!$acc routine( GPXX_OPENACC )
  
  USE PARKIND1, ONLY: JPIM, JPRB, JPRD
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  
  
  
  ! -----------------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PHIHL(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PHIHM(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PHIFL(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PHIFM(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRT(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUF(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVF(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUH(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVH(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PX(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PNHPPI(KPROMA, KFLEV)
  LOGICAL, OPTIONAL, INTENT(IN) :: LDVFE
  
  ! -----------------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  temp (REAL (KIND=JPRB), ZDUF, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZDVF, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZF, (KPROMA, 0:KFLEV + 1))
  LOGICAL :: LLVFE
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  ! -----------------------------------------------------------------------------
  
#include "verdisint_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZDUF) == 8) THEN
    alloc8 (ZDUF)
  ELSE
    IF (KIND (ZDUF) == 8) THEN
      alloc4 (ZDUF)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDVF) == 8) THEN
    alloc8 (ZDVF)
  ELSE
    IF (KIND (ZDVF) == 8) THEN
      alloc4 (ZDVF)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZF) == 8) THEN
    alloc8 (ZF)
  ELSE
    IF (KIND (ZF) == 8) THEN
      alloc4 (ZF)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KST
  
  ! -----------------------------------------------------------------------------
  
  
  ! -----------------------------------------------------------------------------
  
  IF (PRESENT(LDVFE)) THEN
    LLVFE = LDVFE
  ELSE
    LLVFE = YDGEOMETRY%YRVERT_GEOM%YRCVER%LVERTFE
  END IF
  
  IF (LLVFE) THEN
    
    ! vertical derivatives of wind
    ZF(JLON, 0) = 0.0_JPRB
    ZF(JLON, KFLEV + 1) = 0.0_JPRB
    ZF(JLON, 1:KFLEV) = PUF(JLON, 1:KFLEV)
    CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, 'FDER', '11', KPROMA, KST, KEND, KFLEV,  &
    & ZF, ZDUF, YDSTACK=YLSTACK)
    ZF(JLON, 1:KFLEV) = PVF(JLON, 1:KFLEV)
    CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, 'FDER', '11', KPROMA, KST, KEND, KFLEV,  &
    & ZF, ZDVF, YDSTACK=YLSTACK)
    
    IF (PRESENT(PNHPPI)) THEN
      DO JLEV=1,KFLEV
        PX(JLON, JLEV) = REAL(PNHPPI(JLON, JLEV), kind=JPRD) / REAL(PRT(JLON, JLEV)*PLNPR(JLON, JLEV), kind=JPRD) &
        & *REAL(PHIFL(JLON, JLEV)*ZDUF(JLON, JLEV) + PHIFM(JLON, JLEV)*ZDVF(JLON, JLEV), kind=JPRD) /  &
        & YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
      END DO
    ELSE
      DO JLEV=1,KFLEV
        PX(JLON, JLEV) = 1._JPRD / REAL(PRT(JLON, JLEV)*PLNPR(JLON, JLEV), kind=JPRD)*REAL(PHIFL(JLON, JLEV)*ZDUF(JLON, JLEV) +  &
        & PHIFM(JLON, JLEV)*ZDVF(JLON, JLEV), kind=JPRD) / YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
      END DO
    END IF
    
  ELSE
    
    IF (PRESENT(PNHPPI)) THEN
      DO JLEV=1,KFLEV
        PX(JLON, JLEV) = (REAL(PUH(JLON, JLEV) - PUF(JLON, JLEV), kind=JPRD)*PHIHL(JLON, JLEV) + REAL(PVH(JLON, JLEV) - PVF(JLON, &
        &  JLEV), kind=JPRD)*PHIHM(JLON, JLEV) + REAL(PUF(JLON, JLEV) - PUH(JLON, JLEV - 1), kind=JPRD)*PHIHL(JLON, JLEV - 1) +  &
        & REAL(PVF(JLON, JLEV) - PVH(JLON, JLEV - 1), kind=JPRD)*PHIHM(JLON, JLEV - 1))*REAL(PNHPPI(JLON, JLEV), kind=JPRD) /  &
        & REAL(PRT(JLON, JLEV)*PLNPR(JLON, JLEV), kind=JPRD)
      END DO
    ELSE
      DO JLEV=1,KFLEV
        PX(JLON, JLEV) = (REAL(PUH(JLON, JLEV) - PUF(JLON, JLEV), kind=JPRD)*PHIHL(JLON, JLEV) + REAL(PVH(JLON, JLEV) - PVF(JLON, &
        &  JLEV), kind=JPRD)*PHIHM(JLON, JLEV) + REAL(PUF(JLON, JLEV) - PUH(JLON, JLEV - 1), kind=JPRD)*PHIHL(JLON, JLEV - 1) +  &
        & REAL(PVF(JLON, JLEV) - PVH(JLON, JLEV - 1), kind=JPRD)*PHIHM(JLON, JLEV - 1)) / REAL(PRT(JLON, JLEV)*PLNPR(JLON, JLEV), &
        &  kind=JPRD)
      END DO
    END IF
    
  END IF
  
  ! -----------------------------------------------------------------------------
  
  
END SUBROUTINE GPXX_OPENACC
