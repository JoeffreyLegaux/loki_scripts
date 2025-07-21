SUBROUTINE GPGRGEO_EXPL_OPENACC (YDGEOMETRY, KPROMA, KST, KEND, KFLEV, PRT, PRTL, PRTM, PLNPR, PALPH, POROGL, POROGM, PHIFL,  &
& PHIFM, PHIHL, PHIHM, LDNHEE, LDNHHY, PRNHPPI, PQCHAL, PQCHAM, PEPS_NL, PNH1L, PNH1M, PLNPRL_DER, PLNPRM_DER, PALPHL_DER,  &
& PALPHM_DER, YDSTACK)
  
  !**** *GPGRGEO_EXPL* - Computes half and full level gradient of geopotential height "gz".
  
  !     Purpose.
  !     --------
  
  !      Expression of this term is:
  
  !       grad (gz) = grad Phi_s
  !       + grad {int[prehyd'=prehyds to prehyd] (-RT/pre) d prehyd'}
  
  !      where:
  !       - "Phi_s = g z[surf]" is the surface orography.
  !       - "prehyd" is the hydrostatic pressure.
  !       - "pre" is the total pressure including non-hydrostatic effects.
  !       - "prehyds" is the surface hydrostatic pressure.
  !       - "R" is the air constant (including water contribution).
  !       - "T" is the temperature.
  
  !      For the NHQE and the HYD model, pre is replaced by prehyd, some terms disappear.
  !      For the NHQE model, T is a modified temperature
  
  !      Discretisation of the gradient of geopotential height yields:
  
  !      * "grad (gz)" at half level "lbar":
  
  !        (grad (gz))[lbar] = (grad Phi_s)
  !        + sum[k=L to l+1] (prehyd/pre)[k] R[k] T[k] (grad (delta))[k]
  !        + sum[k=L to l+1] (prehyd/pre)[k] (grad (RT))[k] delta[k]
  !        + sum[k=L to l+1] (prehyd/pre)[k] R[k] T[k] delta[k]
  !          { (grad(prehyd)/prehyd)[k] - (grad(pre)/pre)[k] }
  
  !      * "grad (gz)" at full level "l":
  
  !        (grad (gz))[l] = (grad (gz))[lbar]
  !        + (prehyd/pre)[l] R[l] T[l] (grad (alpha))[l]
  !        + (prehyd/pre)[l] (grad (RT))[l] alpha[l]
  !        + (prehyd/pre)[l] R[l] T[l] alpha[l]
  !          { (grad(prehyd)/prehyd)[l] - (grad(pre)/pre)[l] }
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPGRGEO_EXPL(...)
  
  !        Explicit arguments :
  !        --------------------
  !         * INPUT:
  !           YDGEOMETRY   : structure containing all geometry
  !           KPROMA       : horizontal dimension
  !           KD           : start of work
  !           KF           : working length
  !           KFLEV        : number of levels
  !           PRT          : (RT) at full levels
  !           PRTL         : zonal component of "grad (RT)" at full levels
  !           PRTM         : meridian component of "grad (RT)" at full levels
  !           PLNPR        : "delta" at full levels
  !           PALPH        : "alpha" at full levels
  !           POROGL       : zonal component of "grad(surf orography)"
  !           POROGM       : meridian component of "grad(surf orography)"
  
  !         * OUTPUT:
  !           PHIFL        : zonal component of "grad (gz)" at full levels
  !           PHIFM        : merid component of "grad (gz)" at full levels
  !           PHIHL        : zonal component of "grad (gz)" at half levels
  !           PHIHM        : merid component of "grad (gz)" at half levels
  
  !         * INPUT OPTIONAL:
  !           LDNHEE       : .T.: fully elastic non hydrostatic (NHEE) model.
  !                          .F.: hydrostatic or NHQE model.
  !           PRNHPPI      : "prehyd/pre" at full levels.
  !           PQCHAL,PQCHAM: zonal and meridian components at full levels of
  !                          "grad(log(pre/prehyd))=(grad pre)/pre - (grad(prehyd))/prehyd"
  
  !         * OUTPUT OPTIONAL:
  !           PNH1L        : zonal comp of RT grad(log(prehyd/pre)) at full levels
  !           PNH1M        : merid comp of RT grad(log(prehyd/pre)) at full levels
  
  !        Implicit arguments :   None.
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.    None.
  !     ----------
  
  !     Reference.
  !     ----------
  !        See  documentation Arpege ALGORITHME CHAPITRE 6 paragraphe 6
  
  !     Author.
  !     -------
  !      K. YESSAD
  !      Original : 2000-08-11
  
  !     Modifications.
  !     --------------
  !      K. Yessad (Dec 2008): remove dummy CDLOCK
  !      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
  !      K. Yessad (June 2017): Introduce NHQE model.
  !      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
  !      H. Petithomme (Nov 2020): use of pointers for avoiding array copies
  !     ------------------------------------------------------------------
  
!$acc routine( GPGRGEO_EXPL_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, JPHOOK, DR_HOOK
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  
  
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  REAL(KIND=JPRB), INTENT(IN) :: PRT(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRTL(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRTM(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PALPH(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: POROGL(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: POROGM(KPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PHIFL(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PHIFM(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PHIHL(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PHIHM(KPROMA, 0:KFLEV)
  LOGICAL, OPTIONAL, INTENT(IN) :: LDNHEE
  LOGICAL, OPTIONAL, INTENT(IN) :: LDNHHY
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PRNHPPI(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PQCHAL(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PQCHAM(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PEPS_NL(KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PNH1L(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PNH1M(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPRL_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPRM_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PALPHL_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PALPHM_DER(KPROMA, KFLEV)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  LOGICAL :: LLNHEE
  LOGICAL :: LLNHHY
  CHARACTER(LEN=4) :: CLOPER
  temp (REAL (KIND=JPRB), ZNH1L, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZNH1M, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZINL, (KPROMA, 0:KFLEV + 1))
  temp (REAL (KIND=JPRB), ZINM, (KPROMA, 0:KFLEV + 1))
  temp (REAL (KIND=JPRB), ZRNHPPI, (KPROMA, KFLEV))
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  !     ------------------------------------------------------------------
  
#include "abor1.intfb.h"
#include "verdisint_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  IF (KIND (ZNH1L) == 8) THEN
    alloc8 (ZNH1L)
  ELSE
    IF (KIND (ZNH1L) == 4) THEN
      alloc4 (ZNH1L)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZNH1M) == 8) THEN
    alloc8 (ZNH1M)
  ELSE
    IF (KIND (ZNH1M) == 4) THEN
      alloc4 (ZNH1M)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZINL) == 8) THEN
    alloc8 (ZINL)
  ELSE
    IF (KIND (ZINL) == 4) THEN
      alloc4 (ZINL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZINM) == 8) THEN
    alloc8 (ZINM)
  ELSE
    IF (KIND (ZINM) == 4) THEN
      alloc4 (ZINM)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZRNHPPI) == 8) THEN
    alloc8 (ZRNHPPI)
  ELSE
    IF (KIND (ZRNHPPI) == 4) THEN
      alloc4 (ZRNHPPI)
    ELSE
      STOP 1
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  IF (PRESENT(LDNHEE)) THEN
    LLNHEE = LDNHEE
  ELSE
    LLNHEE = .false.
  END IF
  
  IF (PRESENT(LDNHHY)) THEN
    LLNHHY = LDNHHY
  ELSE
    LLNHHY = .false.
  END IF
  
  CLOPER = 'IBOT'
  IF (YDGEOMETRY%YRVERT_GEOM%YRCVER%LVFE_COMPATIBLE) CLOPER = 'INTG'
  
  !     ------------------------------------------------------------------
  
  !*    1. Computation of term:
  !        R[l] T[l] {(grad(prehyd)/prehyd)[l] - (grad(pre)/pre)[l]}
  !        according to "llnhee".
  !*    2.1 Calculation of "grad (gz)" at half levels.
  !         ("delta" and "grad delta" terms contributions.
  
  PHIHL(JLON, KFLEV) = POROGL(JLON)
  PHIHM(JLON, KFLEV) = POROGM(JLON)
  
  IF (LLNHEE) THEN
    IF (.not.(PRESENT(PRNHPPI) .and. PRESENT(PQCHAL) .and. PRESENT(PQCHAM))) CALL ABOR1_ACC(' GPGRGEO_EXPL: missing input  &
    & PRNHPPI, PQCHAL, PQCHAM !!!')
    
    IF (PRESENT(PEPS_NL) .and. LLNHHY) THEN
      DO JLEV=1,KFLEV
        ZRNHPPI(JLON, JLEV) = 1.0_JPRB + PEPS_NL(JLEV)*(PRNHPPI(JLON, JLEV) - 1.0_JPRB)
        ZNH1L(JLON, JLEV) = -PEPS_NL(JLEV)*PRT(JLON, JLEV)*PQCHAL(JLON, JLEV)
        ZNH1M(JLON, JLEV) = -PEPS_NL(JLEV)*PRT(JLON, JLEV)*PQCHAM(JLON, JLEV)
      END DO
    ELSE
      DO JLEV=1,KFLEV
        ZRNHPPI(JLON, JLEV) = PRNHPPI(JLON, JLEV)
        ZNH1L(JLON, JLEV) = -PRT(JLON, JLEV)*PQCHAL(JLON, JLEV)
        ZNH1M(JLON, JLEV) = -PRT(JLON, JLEV)*PQCHAM(JLON, JLEV)
      END DO
    END IF
    
    DO JLEV=KFLEV,1,-1
      PHIHL(JLON, JLEV - 1) = PHIHL(JLON, JLEV) + PLNPR(JLON, JLEV)*ZRNHPPI(JLON, JLEV)*(PRTL(JLON, JLEV) + ZNH1L(JLON, JLEV)) +  &
      & PLNPRL_DER(JLON, JLEV)*ZRNHPPI(JLON, JLEV)*PRT(JLON, JLEV)
      PHIHM(JLON, JLEV - 1) = PHIHM(JLON, JLEV) + PLNPR(JLON, JLEV)*ZRNHPPI(JLON, JLEV)*(PRTM(JLON, JLEV) + ZNH1M(JLON, JLEV)) +  &
      & PLNPRM_DER(JLON, JLEV)*ZRNHPPI(JLON, JLEV)*PRT(JLON, JLEV)
    END DO
    
    ! note: merge tests since both or none present (only present in cpg5_cp)
    IF (PRESENT(PNH1L) .and. PRESENT(PNH1M)) THEN
      PNH1L(JLON, 1:KFLEV) = ZNH1L(JLON, 1:KFLEV)
      PNH1M(JLON, 1:KFLEV) = ZNH1M(JLON, 1:KFLEV)
    END IF
  ELSE
    ! if not present, these arrays are unused
    IF (PRESENT(PNH1L) .and. PRESENT(PNH1M)) THEN
      PNH1L(JLON, 1:KFLEV) = 0.0_JPRB
      PNH1M(JLON, 1:KFLEV) = 0.0_JPRB
    END IF
    
    DO JLEV=KFLEV,1,-1
      PHIHL(JLON, JLEV - 1) = PHIHL(JLON, JLEV) + PLNPR(JLON, JLEV)*PRTL(JLON, JLEV) + PLNPRL_DER(JLON, JLEV)*PRT(JLON, JLEV)
      PHIHM(JLON, JLEV - 1) = PHIHM(JLON, JLEV) + PLNPR(JLON, JLEV)*PRTM(JLON, JLEV) + PLNPRM_DER(JLON, JLEV)*PRT(JLON, JLEV)
    END DO
  END IF
  
  !*    2.2 Calculation of "grad (gz)" at half levels.
  !         "alpha" and "grad alpha" terms contributions.
  
  IF (LLNHEE) THEN
    DO JLEV=1,KFLEV
      ZINL(JLON, JLEV) = (-PLNPR(JLON, JLEV)*ZRNHPPI(JLON, JLEV)*(PRTL(JLON, JLEV) + ZNH1L(JLON, JLEV)) - PLNPRL_DER(JLON, JLEV) &
      & *ZRNHPPI(JLON, JLEV)*PRT(JLON, JLEV))*YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
      ZINM(JLON, JLEV) = (-PLNPR(JLON, JLEV)*ZRNHPPI(JLON, JLEV)*(PRTM(JLON, JLEV) + ZNH1M(JLON, JLEV)) - PLNPRM_DER(JLON, JLEV) &
      & *ZRNHPPI(JLON, JLEV)*PRT(JLON, JLEV))*YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
    END DO
  ELSE
    DO JLEV=1,KFLEV
      ZINL(JLON, JLEV) = (-PLNPR(JLON, JLEV)*PRTL(JLON, JLEV) - PLNPRL_DER(JLON, JLEV)*PRT(JLON, JLEV)) &
      & *YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
      ZINM(JLON, JLEV) = (-PLNPR(JLON, JLEV)*PRTM(JLON, JLEV) - PLNPRM_DER(JLON, JLEV)*PRT(JLON, JLEV)) &
      & *YDGEOMETRY%YRVERT_GEOM%YRVETA%VFE_RDETAH(JLEV)
    END DO
  END IF
  
  ZINL(JLON, 0) = 0.0_JPRB
  ZINL(JLON, KFLEV + 1) = 0.0_JPRB
  ZINM(JLON, 0) = 0.0_JPRB
  ZINM(JLON, KFLEV + 1) = 0.0_JPRB
  CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, CLOPER, '11', KPROMA, KST, KEND, KFLEV,  &
  & ZINL, PHIFL, PINS=POROGL, YDSTACK=YLSTACK)
  CALL VERDISINT_OPENACC(YDGEOMETRY%YRVERT_GEOM%YRVFE, YDGEOMETRY%YRVERT_GEOM%YRCVER, CLOPER, '11', KPROMA, KST, KEND, KFLEV,  &
  & ZINM, PHIFM, PINS=POROGM, YDSTACK=YLSTACK)
  
  !     ------------------------------------------------------------------
  
  
END SUBROUTINE GPGRGEO_EXPL_OPENACC
