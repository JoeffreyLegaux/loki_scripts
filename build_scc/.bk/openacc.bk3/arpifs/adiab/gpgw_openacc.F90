SUBROUTINE GPGW_OPENACC (YDGEOMETRY, LDNHDYN, KFLEV, KPROMA, KST, KEND, LDGWF, LDGDWI, POROGL, POROGM, PLNPR, PALPH, PUS, PVS,  &
& PRT, PDVER, PGWH, PGWF, LDVFE, PRNHPPI, PGDW, YDSTACK)
  
  ! GPGW - Diagnoses "Gw" from the vertical divergence "dver" or from "-G dw".
  
  ! Purpose
  ! -------
  !   Diagnoses "Gw" from the vertical divergence "dver" if LDGDWI=F, from "-G dw" if LDGDWI=T.
  !   For the finite element vertical discretization (lvertfe=T + lvfe_gw=T), "Gw" is computed only at full levels.
  !   For the finite difference vertical discretization, "Gw" is computed at both half and full levels.
  !   Calculation is done by vertical integration of the formula:
  !    dver = - G/(RT) (pre/prehyd) (d w / d log(prehyd))
  !   with the following bottom condition:
  !    w_surf = V_surf grad[Phi_s].
  
  !   This routine can be used in a NHEE model: PRNHPPI should be present in this case.
  
  !   This routine can be used in a NHQE or a hydrostatic model: in this case
  !   the ratio (prehyd/pre) is equal to 1 and PRNHPPI should not be used.
  
  ! Interface
  ! ---------
  !   * INPUT:
  !   YDGEOMETRY   : structure containing all geometry.
  !   KFLEV        : number of levels.
  !   KPROMA       : horizontal dimension.
  !   KSTART       : start of work.
  !   KEND         : end of work.
  !   LDGWF        : calculation of "Gw" at full layers asked for
  !                  if finite difference vertical discretization.
  !   LDGDWI       : T vs F: input content of PDVER is "-G dw" vs "dver".
  !   POROGL       : zonal component of grad[Phi_s].
  !   POROGM       : meridian component of grad[Phi_s].
  !   PLNPR        : "delta" at full layers (computed in GPXYB).
  !   PALPH        : "alpha" at full layers (computed in GPXYB).
  !   PUS          : surface U wind.
  !   PVS          : surface V wind.
  !   PRT          : (RT) at full levels, with the version of R used to define vertical divergence "dver".
  !                  R may be Rdry or Rmoist according to definition of vertical divergence "dver".
  !   PDVER        : vertical divergence "dver" at full layers.
  
  !   * OUTPUT:
  !   PGWH         : G times vertical velocity w at half layers.
  !                  (computed only if lvertfe=F)
  !   PGWF         : G times vertical velocity w at full layers.
  
  !   * OPTIONAL INPUT:
  !   LDVFE        : T if VFE discretisation is used in this routine.
  !   PRNHPPI      : (prehyd/pre) at full layers, required for NHEE.
  
  !   * OPTIONAL OUTPUT:
  !   PGDW         : contains 'G (dw)'.
  
  ! Externals
  ! ---------
  
  ! Method
  ! ------
  
  ! Reference
  ! ---------
  
  ! Author
  ! ------
  !   K. Yessad, Dec 2004 (after GNHSVD2GW)
  
  ! Modifications
  ! -------------
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !   K. Yessad (Dec 2016): Prune obsolete options.
  !   K. Yessad (June 2017): Introduce NHQE model.
  !   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
  !   K. Yessad (Feb 2018): remove deep-layer formulations.
  !   J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
  !   H. Petithomme (Dec 2020): optimisation and test re-organization
  !------------------------------------------------------------------
  
!$acc routine( GPGW_OPENACC ) seq
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  ! -----------------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  LOGICAL, INTENT(IN) :: LDNHDYN
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  LOGICAL, INTENT(IN) :: LDGWF
  LOGICAL, INTENT(IN) :: LDGDWI
  REAL(KIND=JPRB), INTENT(IN) :: POROGL(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: POROGM(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PALPH(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUS(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PVS(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PRT(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDVER(KPROMA, KFLEV)
  REAL(KIND=JPRB), TARGET, INTENT(OUT) :: PGWH(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PGWF(KPROMA, KFLEV)
  LOGICAL, OPTIONAL, INTENT(IN) :: LDVFE
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PRNHPPI(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PGDW(KPROMA, KFLEV)
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  temp (REAL (KIND=JPRB), ZGDW0, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZIN, (KPROMA, 0:KFLEV + 1))
  temp (REAL (KIND=JPRB), ZGWH, (KPROMA, KFLEV + 1))
  temp (REAL (KIND=JPRB), ZGDW, (KPROMA, KFLEV))
  LOGICAL :: LLVFE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "abor1.intfb.h"
#include "verdisint_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZGDW0) == 8) THEN
    alloc8 (ZGDW0)
  ELSE
    IF (KIND (ZGDW0) == 4) THEN
      alloc4 (ZGDW0)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZIN) == 8) THEN
    alloc8 (ZIN)
  ELSE
    IF (KIND (ZIN) == 4) THEN
      alloc4 (ZIN)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZGWH) == 8) THEN
    alloc8 (ZGWH)
  ELSE
    IF (KIND (ZGWH) == 4) THEN
      alloc4 (ZGWH)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KST
  
  
  IF (PRESENT(LDVFE)) THEN
    LLVFE = LDVFE
  ELSE
    LLVFE = YDGEOMETRY%YRVERT_GEOM%YRCVER%LVERTFE .and. (YDGEOMETRY%YRVERT_GEOM%YRCVER%LVFE_GW .or. .not.LDNHDYN)
  END IF
  
  ! optim: use of pointer avoids copying, but dependencies may arise (use ivdep/nodep)
  IF (LLVFE .and. PRESENT(PGDW)) THEN
    assoc (ZGDW,PGDW)
  ELSE
    assoc (ZGDW,ZGDW0)
  END IF
  
  ! * Compute "Gw" at the surface (free slip boundary condition)
  ! optim: use pointer on last level
  PGWH(JLON, KFLEV) = PUS(JLON)*POROGL(JLON) + PVS(JLON)*POROGM(JLON)
  
  ! * Transform "dver" into "-G.dw"
  IF (LDGDWI) THEN
    DO JLEV=1,KFLEV
      ZGDW(JLON, JLEV) = PDVER(JLON, JLEV)
    END DO
  ELSE IF (PRESENT(PRNHPPI)) THEN
    DO JLEV=1,KFLEV
      ZGDW(JLON, JLEV) = PDVER(JLON, JLEV)*PRT(JLON, JLEV)*PLNPR(JLON, JLEV)*PRNHPPI(JLON, JLEV)
    END DO
  ELSE
    DO JLEV=1,KFLEV
      ZGDW(JLON, JLEV) = PDVER(JLON, JLEV)*PRT(JLON, JLEV)*PLNPR(JLON, JLEV)
    END DO
  END IF
  
  ! * Compute "Gw" at full (llvfe=T) or half levels
  IF (LLVFE) THEN
    ! * Store "G dw" at full levels.
    
    DO JLEV=1,KFLEV
      ZIN(JLON, JLEV) = -ZGDW(JLON, JLEV)*YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
    END DO
    
    IF (LDGDWI) THEN
      ZIN(JLON, 0) = ZGDW(JLON, 1)
      ZIN(JLON, KFLEV + 1) = ZGDW(JLON, KFLEV)
      ! Apply RINTBF00, constructed from bottom
      CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, 'ITOP', '00', KPROMA, KST, KEND,  &
      & KFLEV, ZIN, ZGWH, YDSTACK=YLSTACK)
    ELSE IF (LDNHDYN) THEN
      ZIN(JLON, 0) = ZGDW(JLON, 1)
      ZIN(JLON, KFLEV + 1) = ZGDW(JLON, KFLEV)        ! not applied with INGW
      CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, 'INGW', '00', KPROMA, KST, KEND,  &
      & KFLEV, ZIN, ZGWH, YDSTACK=YLSTACK)
    ELSE
      ZIN(JLON, 0) = 0.0_JPRB
      ZIN(JLON, KFLEV + 1) = 0.0_JPRB
      CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, 'IBOT', '11', KPROMA, KST, KEND,  &
      & KFLEV, ZIN, ZGWH, YDSTACK=YLSTACK)
    END IF
    
    DO JLEV=1,KFLEV
      PGWF(JLON, JLEV) = ZGWH(JLON, JLEV) + PGWH(JLON, KFLEV)
    END DO
  ELSE
    ! transform -G.dw into Gw
    DO JLEV=KFLEV,1,-1
      !CDIR NODEP
      !DIR$ IVDEP
      PGWH(JLON, JLEV - 1) = PGWH(JLON, JLEV) + ZGDW(JLON, JLEV)
    END DO
    
    IF (LDGWF) THEN
      IF (LDGDWI) CALL abor1_ACC(' GPGW: compute "Gw" at full levels: case not coded')
      
      ! * Also compute "Gw" at full levels
      ! k.y.: formula pgwf(jlev)=pgwh(jlev)(1-palph(jlev)/plnpr(jlev))
      ! +pgwh(jlev-1)(palph(jlev)/plnpr(jlev)) must be equivalent.
      IF (PRESENT(PRNHPPI)) THEN
        DO JLEV=1,KFLEV
          PGWF(JLON, JLEV) = PGWH(JLON, JLEV) + PDVER(JLON, JLEV)*PRT(JLON, JLEV)*PALPH(JLON, JLEV)*PRNHPPI(JLON, JLEV)
        END DO
      ELSE
        DO JLEV=1,KFLEV
          PGWF(JLON, JLEV) = PGWH(JLON, JLEV) + PDVER(JLON, JLEV)*PRT(JLON, JLEV)*PALPH(JLON, JLEV)
        END DO
      END IF
    END IF
  END IF
  
END SUBROUTINE GPGW_OPENACC
