SUBROUTINE GPGRXYB_EXPL_OPENACC (YDCVER, KPROMA, KST, KEND, KFLEV, LDCOEF, YDVAB, PREL, PREM, PDELP, PLNPR, PRDELP, PALPH,  &
& PRTGR, PRPRE, PRPP, PCOEFD_DER, PLNPRL_DER, PLNPRM_DER, PCOEFA_DER, PCOEFAPL_DER, PALPHPLL_DER, PALPHPLM_DER, PALPHL_DER,  &
& PALPHM_DER, YDSTACK)
  
  !**** *GPGRXYB_EXPL* - Complement to routine "GPXYB".
  !                 Computation of the horizontal gradient of quantities
  !                 "alpha" and "delta" at model levels.
  
  !     Purpose.
  !     --------
  
  !     "alpha" and "delta" are computed at model levels in routine "GPXYB",
  !     but not their horizontal gradient. So this routine provides the
  !     horizontal gradients at full levels. Quantity
  !     "(grad(alpha)) + (grad(prehyd)/prehyd)"
  !     is also provided separately (for case LVERTFE=.F.
  !     its "NDLNPR=0" expression is simpler than the expressions
  !     of "grad(alpha)" and "grad(prehyd)/prehyd").
  !     Discretisation depends on variables "NDLNPR" and "LVERTFE".
  
  !     For LVERTFE=.F., NDLNPR=0, discretisations are:
  
  !      (grad(delta))[l] =
  !      - (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])/(prehyd[lbar]*prehyd[lbar-1])
  !      * (grad prehyds)
  
  !      (grad(alpha))[l] + (grad(prehyd)/prehyd)[l] =
  !      B[lbar]/prehyd[lbar] * (grad prehyds)
  
  !      Quantity "(grad(alpha))[l]" is computed by substracting
  !      "(grad(prehyd)/prehyd)[l]" from
  !      "(grad(alpha))[l] + (grad(prehyd)/prehyd)[l]"
  
  !     For LVERTFE=.F., NDLNPR=1 or 2, discretisations are:
  
  !      (grad(delta))[l] =
  !      - delta[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])
  !      * (1/sqrt(prehyd[lbar]*prehyd[lbar-1])) * (1/(delta prehyd[l]))
  !      * (grad prehyds)
  
  !      (grad(alpha))[l] =
  !      - alpha[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])
  !      * (1/sqrt(prehyd[lbar]*prehyd[lbar-1])) * (1/(delta prehyd[l]))
  !      * (grad prehyds)
  
  !      (grad(prehyd)/prehyd)[l] = prtgr[l] * (grad prehyds)
  !      where "prtgr[l]" is computed in routine "gpxyb" as:
  !      prtgr[l] = { (delta B)[l]
  !      + delta[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])/(delta prehyd[l]) }
  !      * { 1/(delta prehyd[l]) }
  
  !      In this case "(grad(alpha))[l]" is computed prior to
  !      "(grad(alpha))[l] + (grad(prehyd)/prehyd)[l]"
  
  !     For LVERTFE=.T., NDLNPR=0, discretisations are:
  
  !      (grad(delta))[l] =
  !      delta[l] * ((Delta B)[l]/(Delta prehyd)[l] - B[l]/prehyd[l])
  !      * (grad prehyds)
  
  !      grad(alpha) is useless in this case.
  
  !     Notations:
  !      - "grad" is the horizontal gradient operator.
  !        (grad X = vnabla X = M vnabla' X)
  !      - "prehyd" is the hydrostatic pressure.
  !      - "prehyds" is the surface hydrostatic pressure.
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPGRXYB(...)
  
  !        Explicit arguments :
  !        --------------------
  !         * INPUT:
  !           KPROMA       : horizontal dimension
  !           KD           : start of work
  !           KF           : working length
  !           KFLEV        : number of levels
  !           LDCOEF       : if T, stores ZCOEFD, ZCOEFA, ZCOEFAPL in PXYBDER.
  !           YDVAB        : contains information about hybrid vertical coordinate
  !           PREL         : zonal component of "grad prehyds"
  !           PREM         : meridian component of "grad prehyds"
  !           PXYB         : contains pressure depth, "delta", "alpha".
  
  !         * OUTPUT:
  !           PXYBDER      : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
  
  !        Implicit arguments :   None.
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.    None.
  !     ----------
  
  !     Reference.
  !     ----------
  
  !     Author.
  !     -------
  !        K. YESSAD
  !        Original : 00-08-11
  
  !     Modifications.
  !     --------------
  !        K. Yessad (Dec 2008): remove dummy CDLOCK
  !        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
  !        K. Yessad (Dec 2011): use YDVAB.
  !        K. Yessad (June 2017): introduce NDLNPR=2 (for NHQE model).
  !        H. Petithomme (Dec 2020): use of pointers for optimisation
  !     ------------------------------------------------------------------
  
!$acc routine( GPGRXYB_EXPL_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, JPHOOK, DR_HOOK
  USE YOMCVER, ONLY: TCVER
  USE YOMVERT, ONLY: TVAB
  
  
  
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCVER), INTENT(IN) :: YDCVER
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  LOGICAL, INTENT(IN) :: LDCOEF
  TYPE(TVAB), INTENT(IN) :: YDVAB
  REAL(KIND=JPRB), INTENT(IN) :: PREL(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PREM(KPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PALPH(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRTGR(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRPRE(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRPP(KPROMA, KFLEV)
  REAL(KIND=JPRB), TARGET, INTENT(OUT) :: PCOEFD_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLNPRL_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PLNPRM_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), TARGET, INTENT(OUT) :: PCOEFA_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), TARGET, INTENT(OUT) :: PCOEFAPL_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PALPHPLL_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PALPHPLM_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PALPHL_DER(KPROMA, KFLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PALPHM_DER(KPROMA, KFLEV)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  REAL(KIND=JPRB) :: ZCOEFDT
  REAL(KIND=JPRB) :: ZCOEFAT
  REAL(KIND=JPRB) :: ZCOEFAPLT
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !     ------------------------------------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  !*    1/ Calculation of "grad delta" at full levels.
  
  
  ! optim: compilers may report dependence between pxyb and ydvab, ignored with ivdep/nodep
  
  IF (YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      ZCOEFDT = (YDVAB%VDELB(JLEV)*PRDELP(JLON, JLEV) - PRTGR(JLON, JLEV))*PLNPR(JLON, JLEV)
      PLNPRL_DER(JLON, JLEV) = ZCOEFDT*PREL(JLON)
      PLNPRM_DER(JLON, JLEV) = ZCOEFDT*PREM(JLON)
      IF (LDCOEF) THEN
        PCOEFD_DER(JLON, JLEV) = ZCOEFDT
      END IF
    END DO
  ELSE
    DO JLEV=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      ZCOEFDT = -YDVAB%VC(JLEV)*PRPP(JLON, JLEV)
      PLNPRL_DER(JLON, JLEV) = ZCOEFDT*PREL(JLON)
      PLNPRM_DER(JLON, JLEV) = ZCOEFDT*PREM(JLON)
      IF (LDCOEF) THEN
        PCOEFD_DER(JLON, JLEV) = ZCOEFDT
      END IF
    END DO
    
    !*    2/ Calculation of "grad (alpha + log prehyd)" at full levels
    !        discretised as "grad alpha + (grad prehyd) / prehyd ",
    !        and calculation of "grad alpha" at full levels.
    
    IF (YDCVER%NDLNPR == 0) THEN
      DO JLEV=1,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        ZCOEFAPLT = YDVAB%VBH(JLEV)*PRPRE(JLON, JLEV)
        PALPHPLL_DER(JLON, JLEV) = ZCOEFAPLT*PREL(JLON)
        PALPHPLM_DER(JLON, JLEV) = ZCOEFAPLT*PREM(JLON)
        
        ZCOEFAT = ZCOEFAPLT - PRTGR(JLON, JLEV)
        PALPHL_DER(JLON, JLEV) = ZCOEFAT*PREL(JLON)
        PALPHM_DER(JLON, JLEV) = ZCOEFAT*PREM(JLON)
        
        IF (LDCOEF) THEN
          PCOEFA_DER(JLON, JLEV) = ZCOEFAT
          PCOEFAPL_DER(JLON, JLEV) = ZCOEFAPLT
        END IF
        
      END DO
    ELSE IF (YDCVER%NDLNPR == 1 .or. YDCVER%NDLNPR == 2) THEN
      DO JLEV=1,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        ZCOEFAT = -YDVAB%VC(JLEV)*PRPP(JLON, JLEV)*PALPH(JLON, JLEV) / PLNPR(JLON, JLEV)
        PALPHL_DER(JLON, JLEV) = ZCOEFAT*PREL(JLON)
        PALPHM_DER(JLON, JLEV) = ZCOEFAT*PREM(JLON)
        
        ZCOEFAPLT = ZCOEFAT + PRTGR(JLON, JLEV)
        PALPHPLL_DER(JLON, JLEV) = ZCOEFAPLT*PREL(JLON)
        PALPHPLM_DER(JLON, JLEV) = ZCOEFAPLT*PREM(JLON)
        
        IF (LDCOEF) THEN
          PCOEFA_DER(JLON, JLEV) = ZCOEFAT
          PCOEFAPL_DER(JLON, JLEV) = ZCOEFAPLT
        END IF
        
      END DO
    END IF
  END IF
  
  
END SUBROUTINE GPGRXYB_EXPL_OPENACC
