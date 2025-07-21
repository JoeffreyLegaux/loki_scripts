SUBROUTINE GPHPRE_EXPL_OPENACC (YDCVER, TOPPRES, YDCST, KPROMA, KFLEV, KST, KEND, YDVAB, PRESH, PRESF, LHSET, LDELP, LALPHA,  &
& LRTGR, LRPP, PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, YDSTACK)
  
  !**** *GPHPRE_EXPL* - Computes half and full level pressure
  !                Modern version of former GPPRE.
  !                Modern version of former GPPREH+GPXYB+GPPREF
  
  !     Purpose.
  !     --------
  !           Computes pressures at half and full model levels.
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPHPRE_EXPL(...)
  
  !        Explicit arguments :
  !        --------------------
  
  !          KPROMA    : horizontal dimensioning                                (in)
  !          KFLEV     : vertical dimensioning                                  (in)
  !          KST    : start of work                                          (in)
  !          KEND     : depth of work                                          (in)
  !          YDVAB     : contains information about hybrid vertical coordinate  (in)
  !          PRESH     : half level pressure                                    (inout)
  !          PRESF     : full level pressure                                    (opt out)
  !          LDELP,LALPHA,... : activation keys for partial computations        (opt in)
  
  !        Implicit arguments :  NONE.
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.  None.
  !     ----------
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !      K. YESSAD (Sep 2011) after GPPRE, GPPREH, GPXYB and GPPREF.
  
  !     Modifications.
  !     --------------
  !   K. Yessad (Dec 2016): Prune obsolete options.
  !   K. Yessad (Mar 2017): Introduce NDLNPR=2 for NHQE model.
  !   H Petithomme (Dec 2020): add options, use of pointers, group VFE tests
  !     ------------------------------------------------------------------
  
!$acc routine( GPHPRE_EXPL_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, JPHOOK, DR_HOOK
  
  USE YOMCST, ONLY: TCST
  USE YOMVERT, ONLY: TVAB
  USE YOMCVER, ONLY: TCVER
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCVER), INTENT(IN) :: YDCVER
  REAL(KIND=JPRB), INTENT(IN) :: TOPPRES
  TYPE(TCST), INTENT(IN) :: YDCST
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  TYPE(TVAB), INTENT(IN) :: YDVAB
  REAL(KIND=JPRB), INTENT(INOUT) :: PRESH(KPROMA, 0:KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PRESF(KPROMA, KFLEV)
  LOGICAL, OPTIONAL, INTENT(IN) :: LHSET
  LOGICAL, OPTIONAL, INTENT(IN) :: LDELP
  LOGICAL, OPTIONAL, INTENT(IN) :: LALPHA
  LOGICAL, OPTIONAL, INTENT(IN) :: LRTGR
  LOGICAL, OPTIONAL, INTENT(IN) :: LRPP
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PRDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PALPH(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PRTGR(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PRPRE(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PRPP(KPROMA, KFLEV)
  
#include "gphpre_expl_vertfe0_openacc.intfb.h"
#include "gphpre_expl_vertfe1_openacc.intfb.h"
  
  INTEGER(KIND=JPIM) :: JL
  LOGICAL :: LLHSET
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  
  
  LLHSET = .false.
  IF (PRESENT(LHSET)) LLHSET = LHSET
  
  IF (.not.LLHSET) THEN
    DO JL=0,KFLEV - 1
      PRESH(JLON, JL) = YDVAB%VAH(JL) + YDVAB%VBH(JL)*PRESH(JLON, KFLEV)
    END DO
  END IF
  
  IF (YDCVER%LVERTFE) THEN
    CALL GPHPRE_EXPL_VERTFE1_OPENACC(YDCVER, TOPPRES, YDCST, KPROMA, KFLEV, KST, KEND, YDVAB, PRESH, PRESF, LHSET, LDELP,  &
    & LALPHA, LRTGR, LRPP, PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, YDSTACK=YLSTACK)
  ELSE
    CALL GPHPRE_EXPL_VERTFE0_OPENACC(YDCVER, TOPPRES, YDCST, KPROMA, KFLEV, KST, KEND, YDVAB, PRESH, PRESF, LHSET, LDELP,  &
    & LALPHA, LRTGR, LRPP, PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, YDSTACK=YLSTACK)
  END IF
  
  
END SUBROUTINE GPHPRE_EXPL_OPENACC
