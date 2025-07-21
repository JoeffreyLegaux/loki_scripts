SUBROUTINE CPMVVPS_OPENACC (YDCST, YDVAB, KLON, KIDIA, KFDIA, KFLEV, PDT, PFP, PRES, PFEVL, PFEVN, PEVEL, PSDVBC, PSPT1, YDSTACK)
  !  Variables 2D Input
  !  Variables 1D Input
  !  Variables 2D In/Out
  
  !     ------------------------------------------------------------------
  !     MODIFICATION DE LA VITESSE VERTICALE ET DE LA PRESSION DE SURFACE
  !     DANS LE CAS NDPSFI = 1.
  !     ------------------------------------------------------------------
  
  !      VOIR DOCUMENTATION, INTERFACE PHYSICO-DYNAMIQUE
  !                            -----------------
  
  !     ARGUMENTS D ENTREE.
  !     ------------------.
  !       KLON   : DIMENSION HORIZONTALE.
  !       KIDIA  : DEBUT DE LA BOUCLE HORIZONTALE.
  !       KFDIA : BORNE HORIZONTALE DES CALCULS.
  !       KFLEV : DIMENSION ET BORNE VERTICALE.
  !       PDT   : Delta t (SL2TL or first timestep) or 2 Delta t.
  
  ! --- INPUT 2D.
  !     --------.
  !       PFP : FLUX TOTAL DE PRECIPITATIONS LIQUIDES ET NEIGEUSES.
  
  ! --- INPUT 1D.
  !     --------.
  !       PRES : PRESSION DE SURFACE A L'INSTANT OU EST CALCULEE LA PHYSIQUE.
  !              Surface pressure at the same instant as non lagged physics.
  !       PFEVL, PFEVN ( IDEM) : FLUX D'EVAPORATION.
  
  !     ARGUMENTS IMPLICITES
  !     --------------------
  !       CONSTANTES UNIVERSELLES = COMMON /YOMCST/: RG
  !       DECOUPAGE VERTICAL      = YOMGEM: YRVAB%YRVAB%VBH
  
  !     SORTIES
  !     -------
  !       PEVEL (KLON,0:KFLEV) : VITESSE VERTICALE GENERALISEE.
  !       PSDVBC (IDEM) : Integral of divergence term, including
  !                       the "lrubc" and "delta m=1" contributions
  !                       of (etadot DP/DETA), but not the
  !                       "delta m=1" physics.
  !       PSPT1 (KLON) : surface pressure or log of pressure buffer.
  
  !     AUTEUR : E.BAZILE  JUIN 93.
  !     ------
  
  !     INSPIRE DE CPATY ECRIT EN SON TEMPS PAR A.JOLY
  
  !     Modifications:
  !     --------------
  !      K. Yessad (Dec 2008): remove dummy CDLOCK and useless dummy arg
  !      K. Yessad (Jan 2011): remove useless calculations.
  !      K. Yessad (Feb 2018): remove deep-layer formulations.
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !    ------------------------------------------------------------------
  
!$acc routine( CPMVVPS_OPENACC ) seq
  
  USE YOMVERT, ONLY: TVAB
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMCST, ONLY: TCST
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TVAB), INTENT(IN) :: YDVAB
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  REAL(KIND=JPRB), INTENT(IN) :: PDT
  REAL(KIND=JPRB), INTENT(IN) :: PFP(KLON, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRES(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PFEVL(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PFEVN(KLON)
  REAL(KIND=JPRB), INTENT(INOUT) :: PEVEL(KLON, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSDVBC(KLON, 0:KFLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSPT1(KLON)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JLON
  
  temp (REAL (KIND=JPRB), ZFE, (KLON, 0:KFLEV))
  REAL(KIND=JPRB) :: ZEVELS
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  IF (KIND (ZFE) == 8) THEN
    alloc8 (ZFE)
  ELSE
    IF (KIND (ZFE) == 4) THEN
      alloc4 (ZFE)
    ELSE
      STOP 1
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  !     1. - MODIFY :
  !          * The divergence vertical integral term.
  !          * etadot (d prehyd/d eta).
  !     ------------------------------------------------------------------
  
  ! * ZFE : FLUX D EVAPORATION TOTAL
  ZFE(JLON, KFLEV) = PFEVL(JLON) + PFEVN(JLON)
  ZEVELS = ZFE(JLON, KFLEV) + PFP(JLON, KFLEV)
  ZFE(JLON, 0:KFLEV - 1) = 0.0_JPRB
  
  ! * 1.1 Surface/bottom values:
  
  PEVEL(JLON, KFLEV) = PEVEL(JLON, KFLEV) + YDCST%RG*ZEVELS
  PSDVBC(JLON, KFLEV) = PSDVBC(JLON, KFLEV) + YDCST%RG*ZEVELS
  
  ! * 1.2 Other levels:
  
  DO JLEV=1,KFLEV - 1
    PEVEL(JLON, JLEV) = PEVEL(JLON, JLEV) + YDCST%RG*(YDVAB%VBH(JLEV)*ZEVELS)
    PSDVBC(JLON, JLEV) = PSDVBC(JLON, JLEV) + YDCST%RG*(YDVAB%VBH(JLEV)*ZEVELS)
  END DO
  
  !    ------------------------------------------------------------------
  !     2. - ADD PHYSICS TO PSPT1.
  !    ------------------------------------------------------------------
  
  PSPT1(JLON) = PSPT1(JLON) - PDT*YDCST%RG*(PFP(JLON, KFLEV) + ZFE(JLON, KFLEV)) / PRES(JLON)
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE CPMVVPS_OPENACC
