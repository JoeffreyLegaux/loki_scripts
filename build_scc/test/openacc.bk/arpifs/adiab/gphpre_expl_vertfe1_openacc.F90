SUBROUTINE GPHPRE_EXPL_VERTFE1_OPENACC (YDCVER, TOPPRES, YDCST, KPROMA, KFLEV, KST, KEND, YDVAB, PRESH, PRESF, LHSET, LDELP,  &
& LALPHA, LRTGR, LRPP, PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, YDSTACK)
  
  !**** *GPHPRE_EXPL_VERTFE1* - Computes half and full level pressure
  !                Modern version of former GPPRE.
  !                Modern version of former GPPREH+GPXYB+GPPREF
  
  !     Purpose.
  !     --------
  !           Computes pressures at half and full model levels.
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPHPRE_EXPL_VERTFE1(...)
  
  !        Explicit arguments :
  !        --------------------
  
  !          KPROMA    : horizontal dimensioning                                (in)
  !          KFLEV     : vertical dimensioning                                  (in)
  !          KSTART    : start of work                                          (in)
  !          KPROF     : depth of work                                          (in)
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
  
!$acc routine( GPHPRE_EXPL_VERTFE1_OPENACC ) seq
  
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
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PRESF(KPROMA, KFLEV)
  LOGICAL, OPTIONAL, INTENT(IN) :: LHSET
  LOGICAL, OPTIONAL, INTENT(IN) :: LDELP
  LOGICAL, OPTIONAL, INTENT(IN) :: LALPHA
  LOGICAL, OPTIONAL, INTENT(IN) :: LRTGR
  LOGICAL, OPTIONAL, INTENT(IN) :: LRPP
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PLNPR(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PRDELP(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PALPH(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PRTGR(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PRPRE(KPROMA, KFLEV)
  REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(OUT) :: PRPP(KPROMA, KFLEV)
  
  INTEGER(KIND=JPIM) :: IFIRST
  INTEGER(KIND=JPIM) :: JL
  INTEGER(KIND=JPIM) :: JROF
  LOGICAL :: LLDELP
  LOGICAL :: LLALPHA
  LOGICAL :: LLRTGR
  LOGICAL :: LLRPP
  LOGICAL :: LTEST
  LOGICAL :: LLXYB
  REAL(KIND=JPRB) :: ZPRE
  
  temp (REAL (KIND=JPRB), ZZPRESF, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZDELP, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZLNPR, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZRDELP, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZALPH, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZRTGR, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZRPRE, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZZRPP, (KPROMA, KFLEV))
  
  temp (REAL (KIND=JPRB), ZDELP, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZLNPR, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZRDELP, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZALPH, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZRTGR, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZRPRE, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZRPP, (KPROMA, KFLEV))
  temp (REAL (KIND=JPRB), ZPRESF, (KPROMA, KFLEV))
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  IF (KIND (ZZPRESF) == 8) THEN
    alloc8 (ZZPRESF)
  ELSE
    IF (KIND (ZZPRESF) == 4) THEN
      alloc4 (ZZPRESF)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZDELP) == 8) THEN
    alloc8 (ZZDELP)
  ELSE
    IF (KIND (ZZDELP) == 4) THEN
      alloc4 (ZZDELP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZLNPR) == 8) THEN
    alloc8 (ZZLNPR)
  ELSE
    IF (KIND (ZZLNPR) == 4) THEN
      alloc4 (ZZLNPR)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZRDELP) == 8) THEN
    alloc8 (ZZRDELP)
  ELSE
    IF (KIND (ZZRDELP) == 4) THEN
      alloc4 (ZZRDELP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZALPH) == 8) THEN
    alloc8 (ZZALPH)
  ELSE
    IF (KIND (ZZALPH) == 4) THEN
      alloc4 (ZZALPH)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZRTGR) == 8) THEN
    alloc8 (ZZRTGR)
  ELSE
    IF (KIND (ZZRTGR) == 4) THEN
      alloc4 (ZZRTGR)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZRPRE) == 8) THEN
    alloc8 (ZZRPRE)
  ELSE
    IF (KIND (ZZRPRE) == 4) THEN
      alloc4 (ZZRPRE)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZZRPP) == 8) THEN
    alloc8 (ZZRPP)
  ELSE
    IF (KIND (ZZRPP) == 4) THEN
      alloc4 (ZZRPP)
    ELSE
      STOP 1
    END IF
  END IF
  
  
  LLXYB = PRESENT(PDELP) .and. PRESENT(PLNPR) .and. PRESENT(PRDELP) .and. PRESENT(PALPH) .and. PRESENT(PRTGR) .and. PRESENT(PRPRE &
  & ) .and. PRESENT(PRPP)
  
  IF (LLXYB) THEN
    LLDELP = .true.
    IF (PRESENT(LDELP)) LLDELP = LDELP
    LLALPHA = .true.
    IF (PRESENT(LALPHA)) LLALPHA = LALPHA
    LLRTGR = .true.
    IF (PRESENT(LRTGR)) LLRTGR = LRTGR
    
    assoc (ZDELP,PDELP)
    assoc (ZLNPR,PLNPR)
    assoc (ZRDELP,PRDELP)
    assoc (ZALPH,PALPH)
    assoc (ZRTGR,PRTGR)
    assoc (ZRPRE,PRPRE)
    assoc (ZRPP,PRPP)
    
    IF (PRESENT(PRESF)) THEN
      assoc (ZPRESF,PRESF)
    ELSE
      assoc (ZPRESF,ZZPRESF)
    END IF
  ELSE
    nullptr(ZPRESF)
  END IF
  
  
  IF (LLXYB) THEN
    DO JL=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      ZDELP(JLON, JL) = YDVAB%VDELA(JL) + YDVAB%VDELB(JL)*PRESH(JLON, KFLEV)
      ZPRESF(JLON, JL) = YDVAB%VAF(JL) + YDVAB%VBF(JL)*PRESH(JLON, KFLEV)
      ZPRE = 1._JPRB / ZPRESF(JLON, JL)
      ZLNPR(JLON, JL) = ZDELP(JLON, JL)*ZPRE
      
      IF (LLDELP) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        ZRDELP(JLON, JL) = 1._JPRB / ZDELP(JLON, JL)
      END IF
      
      IF (LLALPHA) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        ZALPH(JLON, JL) = (PRESH(JLON, JL) - ZPRESF(JLON, JL))*ZPRE
      END IF
      
      IF (LLRTGR) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        ZRTGR(JLON, JL) = YDVAB%VBF(JL)*ZPRE
      END IF
    END DO
  ELSE IF (PRESENT(PRESF)) THEN
    DO JL=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      PRESF(JLON, JL) = YDVAB%VAF(JL) + YDVAB%VBF(JL)*PRESH(JLON, KFLEV)
    END DO
  END IF
  
  
END SUBROUTINE GPHPRE_EXPL_VERTFE1_OPENACC
