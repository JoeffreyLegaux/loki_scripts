SUBROUTINE GPHPRE_EXPL_VERTFE0_OPENACC (YDCVER, TOPPRES, YDCST, KPROMA, KFLEV, KST, KEND, YDVAB, PRESH, PRESF, LHSET, LDELP,  &
& LALPHA, LRTGR, LRPP, PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, YDSTACK)
  
  !**** *GPHPRE_EXPL_VERTFE0* - Computes half and full level pressure
  !                Modern version of former GPPRE.
  !                Modern version of former GPPREH+GPXYB+GPPREF
  
  !     Purpose.
  !     --------
  !           Computes pressures at half and full model levels.
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPHPRE_EXPL_VERTFE0(...)
  
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
  
!$acc routine( GPHPRE_EXPL_VERTFE0_OPENACC ) seq
  
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
  REAL(KIND=JPRB) :: ZCOUNT
  REAL(KIND=JPRB) :: ZSUM
  
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
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
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
    LLRPP = .true.
    IF (PRESENT(LRPP)) LLRPP = LRPP
    
    ! broader condition for computing alpha:
    IF (PRESENT(PRESF)) LLALPHA = LLALPHA .or. YDCVER%NDLNPR == 1 .or. YDCVER%NDLNPR == 2 .or. .not.YDCVER%LAPRXPK
    
    assoc (ZDELP,PDELP)
    assoc (ZLNPR,PLNPR)
    assoc (ZRDELP,PRDELP)
    assoc (ZALPH,PALPH)
    assoc (ZRTGR,PRTGR)
    assoc (ZRPRE,PRPRE)
    assoc (ZRPP,PRPP)
    
  ELSE
    LLDELP = .false.
    LLRTGR = .false.
    LLRPP = .false.
    
    ! reduced condition for computing alpha:
    LLALPHA = PRESENT(PRESF) .and. (YDCVER%NDLNPR == 1 .or. YDCVER%NDLNPR == 2 .or. .not.YDCVER%LAPRXPK)
    
    IF (LLALPHA) THEN
      assoc (ZDELP,ZZDELP)
      assoc (ZLNPR,ZZLNPR)
      assoc (ZRDELP,ZZRDELP)
      assoc (ZALPH,ZZALPH)
      assoc (ZRTGR,ZZRTGR)
      assoc (ZRPRE,ZZRPRE)
      assoc (ZRPP,ZZRPP)
    END IF
  END IF
  
  IF (LLXYB .or. LLALPHA) THEN
    ! pressure at top
    ZPRE = YDVAB%VAH(0) + YDVAB%VBH(0)*PRESH(JLON, KFLEV)
    
    ZCOUNT = 0._JPRB
    
    IF (ZPRE <= TOPPRES) THEN
      ZCOUNT = ZCOUNT + 1
    END IF
    
    ZSUM = ZCOUNT
    
    IF (ZSUM > 0._JPRB) THEN
      IFIRST = 2
    ELSE
      IFIRST = 1
    END IF
    
    IF (YDCVER%NDLNPR == 0) THEN
      IF (IFIRST == 2) THEN
        ZLNPR(JLON, 1) = LOG(PRESH(JLON, 1) / TOPPRES)
        
        IF (LLDELP .or. LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZDELP(JLON, 1) = PRESH(JLON, 1) - PRESH(JLON, 0)
          ZRDELP(JLON, 1) = 1._JPRB / ZDELP(JLON, 1)
        END IF
        
        IF (LLALPHA) ZALPH(JLON, 1) = YDCVER%RHYDR0
        
        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZRTGR(JLON, 1) = ZRDELP(JLON, 1)*YDVAB%VDELB(1)
        END IF
        
        IF (LLRPP) THEN
          ZRPRE(JLON, 1) = 1._JPRB / PRESH(JLON, 1)
          ZRPP(JLON, 1) = 1._JPRB / (PRESH(JLON, 1)*TOPPRES)
        END IF
      END IF
      
      ZPRE = 1._JPRB / PRESH(JLON, IFIRST - 1)
      
      DO JL=IFIRST,KFLEV
        ZLNPR(JLON, JL) = LOG(PRESH(JLON, JL)*ZPRE)
        
        IF (LLDELP .or. LLALPHA .or. LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZDELP(JLON, JL) = PRESH(JLON, JL) - PRESH(JLON, JL - 1)
          ZRDELP(JLON, JL) = 1._JPRB / ZDELP(JLON, JL)
        END IF
        
        IF (LLALPHA) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZALPH(JLON, JL) = 1._JPRB - PRESH(JLON, JL - 1)*ZRDELP(JLON, JL)*ZLNPR(JLON, JL)
        END IF
        
        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZRTGR(JLON, JL) = ZRDELP(JLON, JL)*(YDVAB%VDELB(JL) + YDVAB%VC(JL)*ZLNPR(JLON, JL)*ZRDELP(JLON, JL))
        END IF
        
        IF (LLRPP) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZRPRE(JLON, JL) = 1._JPRB / PRESH(JLON, JL)
          ZRPP(JLON, JL) = ZRPRE(JLON, JL)*ZPRE
          ZPRE = ZRPRE(JLON, JL)
        ELSE
          ZPRE = 1._JPRB / PRESH(JLON, JL)
        END IF
      END DO
    ELSE IF (YDCVER%NDLNPR == 1 .or. YDCVER%NDLNPR == 2) THEN
      LTEST = LLDELP .or. LLALPHA .or. LLRTGR
      
      DO JL=IFIRST,KFLEV
        ! optim: keep lnpr separate but consistent
        IF (LTEST) THEN
          ZDELP(JLON, JL) = PRESH(JLON, JL) - PRESH(JLON, JL - 1)
          ZRDELP(JLON, JL) = 1._JPRB / ZDELP(JLON, JL)
          ZRPP(JLON, JL) = 1._JPRB / (PRESH(JLON, JL)*PRESH(JLON, JL - 1))
          ZLNPR(JLON, JL) = ZDELP(JLON, JL)*SQRT(ZRPP(JLON, JL))
        ELSE
          ZLNPR(JLON, JL) = (PRESH(JLON, JL) - PRESH(JLON, JL - 1)) / SQRT(PRESH(JLON, JL)*PRESH(JLON, JL - 1))
        END IF
        
        IF (LLALPHA) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZALPH(JLON, JL) = 1._JPRB - PRESH(JLON, JL - 1)*ZRDELP(JLON, JL)*ZLNPR(JLON, JL)
        END IF
        
        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZRTGR(JLON, JL) = ZRDELP(JLON, JL)*(YDVAB%VDELB(JL) + YDVAB%VC(JL)*ZLNPR(JLON, JL)*ZRDELP(JLON, JL))
        END IF
        
        IF (LLRPP) ZRPRE(JLON, JL) = 1._JPRB / PRESH(JLON, JL)
      END DO
      
      IF (IFIRST == 2) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        ZDELP(JLON, 1) = PRESH(JLON, 1)
        ZRDELP(JLON, 1) = 1._JPRB / ZDELP(JLON, 1)
        
        IF (YDCVER%NDLNPR == 1) THEN
          ZLNPR(JLON, 1) = 2._JPRB + YDCST%RCVD / YDCST%RD
        ELSE
          !DIR$ IVDEP
          !CDIR NODEP
          ZLNPR(JLON, 1) = 1._JPRB + ZLNPR(JLON, 2)*(ZDELP(JLON, 1) / ZDELP(JLON, 2))*SQRT(PRESH(JLON, 2) / PRESH(JLON, 1))
        END IF
        
        IF (LLALPHA) ZALPH(JLON, 1) = 1._JPRB
        
        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZRTGR(JLON, 1) = ZRDELP(JLON, 1)*YDVAB%VDELB(1)
        END IF
        
        IF (LLRPP) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          ZRPRE(JLON, 1) = 1._JPRB / PRESH(JLON, 1)
          ZRPP(JLON, 1) = (ZLNPR(JLON, 1)*ZRDELP(JLON, 1))**2
        END IF
      END IF
    END IF
  END IF
  
  IF (PRESENT(PRESF)) THEN
    IF (YDCVER%NDLNPR == 0) THEN
      IF (YDCVER%LAPRXPK) THEN
        DO JL=1,KFLEV
          PRESF(JLON, JL) = (PRESH(JLON, JL - 1) + PRESH(JLON, JL))*0.5_JPRB
        END DO
      ELSE
        DO JL=1,KFLEV
          !DIR$ IVDEP
          !CDIR NODEP
          PRESF(JLON, JL) = EXP(-ZALPH(JLON, JL))*PRESH(JLON, JL)
        END DO
      END IF
    ELSE IF (YDCVER%NDLNPR == 1 .or. YDCVER%NDLNPR == 2) THEN
      DO JL=IFIRST,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        PRESF(JLON, JL) = (1._JPRB - ZALPH(JLON, JL))*PRESH(JLON, JL)
      END DO
      
      IF (IFIRST == 2) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        PRESF(JLON, 1) = PRESH(JLON, 1) / ZLNPR(JLON, 1)
      END IF
    END IF
  END IF
  
  
END SUBROUTINE GPHPRE_EXPL_VERTFE0_OPENACC
