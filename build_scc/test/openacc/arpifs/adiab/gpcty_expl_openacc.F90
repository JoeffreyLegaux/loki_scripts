SUBROUTINE GPCTY_EXPL_OPENACC (YDVFE, YDCVER, KPROMA, KST, KEND, KFLEV, LDRUBC, YDVAB, YDVETA, PU, PV, PD, PEVT, PSPL, PSPM,  &
& PRPREF, PDPHYCTY, PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, PEVEL, PVVEL, PPSDIV, PPSDVBC, PDIVDP, YDSTACK)
  
  !**** *GPCTY_EXPL* - Computes vertical velocities.
  
  !     Purpose.
  !     --------
  
  !      Computes the vertical velocities
  !       "etadot (d prehyd / d eta)" at half levels
  !       "(omega / prehyd)" at full levels
  
  !      ----- The following discretisations are valid for lvertfe=.false. -----
  
  !      Omitting the "delta m=1" flux precipitation terms, which will be later
  !      added in routine "cpmvvps", the discretised expression of
  !      "etadot (d prehyd / d eta)" on the half level number "lbar" is:
  
  !       etadot (d prehyd / d eta) [lbar] =
  !        B[lbar] * { sum[k=1 to L]
  !        (Delta B)[k] * vec(V)[k] * (grad prehyds) }
  !        + B[lbar] * { sum[k=1 to L] (Delta prehyd)[k] * (grad vec(V)[k]) }
  !        - sum[k=1 to l] (Delta B)[k] * vec(V)[k] * (grad prehyds)
  !        - sum[k=1 to l] (Delta prehyd)[k] * (grad vec(V)[k])
  !        + { 1 - B[lbar] } * (etadot (d prehyd / d eta))[top]
  
  !      where:
  !       - "vec(V)" is the horizontal wind.
  !       - "prehyds" is the surface hydrostatic pressure.
  !       - "grad" is the "horizontal gradient" operator:
  !         grad X = vnabla X = M vnabla' X
  
  !      Omitting the "delta m=1" flux precipitation terms, which will be later
  !      added in routine "cpmvvps", the discretised expression of
  !      "(omega / prehyd)" on the full level number "l" is:
  
  !       (omega / prehyd)[l] =
  !        vec(V)[l] (grad prehyd / prehyd)[l]
  !        - delta[l]/(Delta prehyd)[l] * { sum[k=1 to l-1]
  !        (Delta B)[k] * vec(V)[k] * (grad prehyds) }
  !        - delta[l]/(Delta prehyd)[l] * { sum[k=1 to l-1]
  !        (Delta prehyd)[k] * (grad vec(V)[k]) }
  !        - alpha[l]/(Delta prehyd)[l] * (Delta B)[l] * vec(V)[l] * (grad prehyds)
  !        - alpha[l] * (grad vec(V)[l])
  !        + delta[l]/(Delta prehyd)[l] * (etadot (d prehyd / d eta))[top]
  
  !      This routine stores additional quantities:
  
  !      * vertical integral of divergence without the "lrubc" contribution:
  !        for the half level number "lbar" its discretised expression is:
  
  !        psdiv[lbar] = sum[k=1 to l]
  !        { (Delta B)[k] * vec(V)[k] * (grad prehyds)
  !        + (Delta prehyd)[k] * (grad vec(V)[k]) }
  
  !      * vertical integral of divergence with the "lrubc" contribution:
  !        for the half level number "lbar" its discretised expression is:
  
  !        psdvbc[lbar] = sum[k=1 to l]
  !        { (Delta B)[k] * vec(V)[k] * (grad prehyds)
  !        + (Delta prehyd)[k] * (grad vec(V)[k]) }
  !        - (etadot (d prehyd / d eta))[top]
  
  !      * divergence term
  !        "pdivdp=grad( vec(V) * (Delta prehyd) )"
  !        at full levels: for the full level number "l" its discretised expression is:
  
  !        pdivdp[l] =
  !        (Delta B)[l] * vec(V)[l] * (grad prehyds)
  !        + (Delta prehyd)[l] * (grad vec(V)[l])
  
  !      -----------------------------------------------------------------------
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPCTY_EXPL(...)
  
  !        Explicit arguments :
  !        --------------------
  !        * INPUT:
  !          KPROMA       : horizontal dimensioning
  !          KST           : start of work
  !          KEND           : depth of work
  !          KFLEV        : number of layers
  !          LDRUBC       : upper boundary condition switch
  !          YDVAB        : contains information about hybrid vertical coordinate
  !          YDVETA       : contains information about hybrid vertical coordinate: "eta"
  !          PU           : U component of the wind, at full levels
  !          PV           : V component of the wind, at full levels
  !          PD           : horizontal divergence of the wind, at full levels
  !          PEVT         : "etadot (d prehyd / d eta)" at the top
  !          PXYB         : contains pressure depth, "delta", "alpha".
  !          PSPL         : zonal component of "grad prehyds"
  !          PSPM         : meridian component of "grad prehyds"
  !          PRPREF       : inverse of full level pressures.
  
  !        * OUTPUT:
  !          PCTY         : contains vertical velocities, vertical integral of divergence.
  
  !        * OPTIONAL INPUT:
  !          PDPHYCTY     : mass source/sink from physics (previous dt in LAG physics)
  
  !        Implicit arguments :   None.
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.    None.
  !     ----------
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !      Mats Hamrud and Philippe Courtier  *ECMWF*
  !      Original : 88-02-04
  
  !     Modifications.
  !     --------------
  !   K. Yessad (Jan 2008): complete (LVERCOR,LVERTFE)=(T,T).
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !   K. Yessad (Dec 2008): remove dummy CDLOCK
  !   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
  !   K. Yessad (Dec 2011): use YDVAB.
  !   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
  !   K. Yessad (Feb 2018): remove deep-layer formulations.
  !   H Petithomme (Dec 2020): use of pointer, add directives for dependencies
  !     ------------------------------------------------------------------
  
!$acc routine( GPCTY_EXPL_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, JPHOOK, DR_HOOK
  
  USE YOMVERT, ONLY: TVFE, TVAB, TVETA
  USE YOMCVER, ONLY: TCVER
  
  
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TVFE), INTENT(IN) :: YDVFE
  TYPE(TCVER), INTENT(IN) :: YDCVER
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  LOGICAL, INTENT(IN) :: LDRUBC
  TYPE(TVAB), INTENT(IN) :: YDVAB
  TYPE(TVETA), INTENT(IN) :: YDVETA
  REAL(KIND=JPRB), INTENT(IN) :: PU(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PV(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PD(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PEVT(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPL(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPM(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PRPREF(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PDPHYCTY(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PALPH(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRTGR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRPRE(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRPP(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PEVEL(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PVVEL(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PPSDIV(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PPSDVBC(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PDIVDP(KPROMA, 0:KFLEV)
  
  !     ------------------------------------------------------------------
  
  temp (REAL (KIND=JPRB), ZSDIV, (KPROMA, 0:KFLEV + 1))
  temp (REAL (KIND=JPRB), ZPSDIVFE, (KPROMA))
  temp (REAL (KIND=JPRB), ZVP, (KPROMA, KFLEV))
  
  INTEGER(KIND=JPIM) :: JROF
  INTEGER(KIND=JPIM) :: JLEV
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  !     ------------------------------------------------------------------
  
#include "abor1.intfb.h"
#include "verdisint_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  IF (KIND (ZSDIV) == 8) THEN
    alloc8 (ZSDIV)
  ELSE
    IF (KIND (ZSDIV) == 4) THEN
      alloc4 (ZSDIV)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZPSDIVFE) == 8) THEN
    alloc8 (ZPSDIVFE)
  ELSE
    IF (KIND (ZPSDIVFE) == 4) THEN
      alloc4 (ZPSDIVFE)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZVP) == 8) THEN
    alloc8 (ZVP)
  ELSE
    IF (KIND (ZVP) == 4) THEN
      alloc4 (ZVP)
    ELSE
      STOP 1
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  !     check for non-supported configurations
  
  IF (LDRUBC .and. YDCVER%LVERTFE) CALL ABOR1_ACC('GPCTY_EXPL: BAD OPTIONS')
  
  !     ------------------------------------------------------------------
  
  !*    1. Computes "vec(V) * grad prehyds" at full levels.
  !     ---------------------------------------------------
  
  DO JLEV=1,KFLEV
    ZVP(JLON, JLEV) = PU(JLON, JLEV)*PSPL(JLON) + PV(JLON, JLEV)*PSPM(JLON)
  END DO
  
  !     ------------------------------------------------------------------
  
  !*    2. Sum divergence.
  !     ------------------
  
  IF (PRESENT(PDPHYCTY)) THEN
    DO JLEV=1,KFLEV
      PDIVDP(JLON, JLEV) =  &
      & PD(JLON, JLEV)*PDELP(JLON, JLEV) + YDVAB%VDELB(JLEV)*ZVP(JLON, JLEV) - PDPHYCTY(JLON, JLEV)*PDELP(JLON, JLEV)
    END DO
  ELSE
    DO JLEV=1,KFLEV
      PDIVDP(JLON, JLEV) = PD(JLON, JLEV)*PDELP(JLON, JLEV) + YDVAB%VDELB(JLEV)*ZVP(JLON, JLEV)
    END DO
  END IF
  
  
  
  PPSDIV(JLON, 0) = 0.0_JPRB
  
  IF (YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
      ZSDIV(JLON, JLEV) = PDIVDP(JLON, JLEV)*YDVETA%VFE_RDETAH(JLEV)
    END DO
    
    ZSDIV(JLON, 0) = 0.0_JPRB
    ZSDIV(JLON, KFLEV + 1) = 0.0_JPRB
    CALL VERDISINT_OPENACC(YDVFE, YDCVER, 'ITOP', '11', KPROMA, KST, KEND, KFLEV, ZSDIV, PPSDIV(:, 1:KFLEV), POUTS=ZPSDIVFE,  &
    & YDSTACK=YLSTACK)
  ELSE
    DO JLEV=1,KFLEV
      !dir$ ivdep
      PPSDIV(JLON, JLEV) = PPSDIV(JLON, JLEV - 1) + PDIVDP(JLON, JLEV)
    END DO
  END IF
  
  IF (LDRUBC) THEN
    IF (YDCVER%LVERTFE) THEN
      PPSDVBC(JLON, KFLEV) = ZPSDIVFE(JLON) - PEVT(JLON)
    ELSE
      DO JLEV=0,KFLEV
        !dir$ ivdep
        PPSDVBC(JLON, JLEV) = PPSDIV(JLON, JLEV) - PEVT(JLON)
      END DO
    END IF
  ELSE
    IF (YDCVER%LVERTFE) THEN
      PPSDVBC(JLON, KFLEV) = ZPSDIVFE(JLON)
    ELSE
      DO JLEV=0,KFLEV
        !dir$ ivdep
        PPSDVBC(JLON, JLEV) = PPSDIV(JLON, JLEV)
      END DO
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
  !*    3. Computes "etadot (d prehyd / d eta)".
  !     ----------------------------------------
  
  ! "etadot (d prehyd / d eta)" is computed at full levels if LVERTFE=T,
  !  at half levels otherwise.
  
  IF (YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
      !dir$ ivdep
      PEVEL(JLON, JLEV) = YDVAB%VBF(JLEV)*ZPSDIVFE(JLON) - PPSDIV(JLON, JLEV)
    END DO
  ELSE
    DO JLEV=1,KFLEV - 1
      !dir$ ivdep
      PEVEL(JLON, JLEV) = YDVAB%VBH(JLEV)*PPSDIV(JLON, KFLEV) - PPSDIV(JLON, JLEV)
    END DO
  END IF
  
  PEVEL(JLON, 0) = 0.0_JPRB
  IF (.not.YDCVER%LVERTFE) PEVEL(JLON, KFLEV) = 0.0_JPRB
  
  IF (LDRUBC) THEN
    IF (YDCVER%LVERTFE) THEN
      CALL ABOR1_ACC(' GPCTY_EXPL: VFE not coded for LRUBC.')
    ELSE
      DO JLEV=0,KFLEV
        PEVEL(JLON, JLEV) = PEVEL(JLON, JLEV) + (1.0_JPRB - YDVAB%VBH(JLEV))*PEVT(JLON)
      END DO
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
  !*    4. Computes "(omega/prehyd)" at full levels.
  !     --------------------------------------------
  
  IF (YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
      !dir$ ivdep
      PVVEL(JLON, JLEV) = (YDVAB%VBF(JLEV)*ZVP(JLON, JLEV) - PPSDIV(JLON, JLEV))*PRPREF(JLON, JLEV)
    END DO
  ELSE
    DO JLEV=1,KFLEV
      !dir$ ivdep
      PVVEL(JLON, JLEV) = PRTGR(JLON, JLEV)*ZVP(JLON, JLEV) - PRDELP(JLON, JLEV)*PALPH(JLON, JLEV)*YDVAB%VDELB(JLEV)*ZVP(JLON,  &
      & JLEV) - PALPH(JLON, JLEV)*PD(JLON, JLEV) - PRDELP(JLON, JLEV)*PPSDVBC(JLON, JLEV - 1)*PLNPR(JLON, JLEV)
    END DO
  END IF
  
  !     ------------------------------------------------------------------
  
  
END SUBROUTINE GPCTY_EXPL_OPENACC
