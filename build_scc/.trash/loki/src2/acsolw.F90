SUBROUTINE ACSOLW_OPENACC (YDPHY1, KIDIA, KFDIA, KLON, PARG, PD2, PLSM, PIVEG, PSAB, LDHMT, PWFC, PWPMX, PWSAT, PWSMX, PWWILT,  &
& YDSTACK)
  !-----------------------------------------------------------------------
  ! - INPUT  1D
  ! - INPUT  LOGIQUE
  ! - OUTPUT 1D .
  
  !**** *ACSOLW  * - DETERMINATION DES CONTENUS EN EAU CARACTERISTIQUES DU SOL
  
  !     Sujet.
  !     ------
  
  !     - ROUTINE DE CALCUL INTERMEDIAIRE.
  !       DETERMINATION DES CONTENUS EN EAU CARACTERISTIQUES DU SOL.
  
  !    SUR MER      (PLSM=0. , PIVEG=NTVMER    )
  !       wwilt=0.        wfc=1.        wsat=1.
  !    SUR BANQUISE (PLSM=1. , PIVEG=NTVGLA*1.1) (MODE CLIMAT SEULEMENT)
  !       wwilt=0.        wfc=1.        wsat=1.
  !    SUR GLACIER  (PLSM=1. , PIVEG=NTVGLA    )
  !       wwilt=f1(arg)   wfc=f2(arg)   wsat=wfc
  !    SUR TERRE    (PLSM=1. , AUTRES PIVEG    )
  !       wwilt=f1(arg)   wfc=f2(arg)   wsat=f3(sab)
  
  !**   Interface.
  !     ----------
  !        *CALL* *ACSOLW*
  
  !-----------------------------------------------------------------------
  ! WARNING: THE ENGLISH VERSION OF VARIABLES' NAMES IS TO BE READ IN THE
  !          "APLPAR" CODE.
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS D'ENTREE.
  !     -------------------
  
  ! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.
  
  ! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
  ! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
  ! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
  
  ! - NOM DES CLES LOGIQUES
  
  ! LDHMT      : CALCULS REDUITS (POUR L'ANALYSE OU FULL-POS)
  !              SEULS LES TABLEAUX (*) SONT UTILISES/CALCULES EN E/S
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (GEOGRAPHIQUE) .
  
  ! PARG (*)   : POURCENTAGE D'ARGILE DANS LA MAILLE.
  ! PD2        : EPAISSEUR DU RESERVOIR PROFOND.
  ! PLSM (*)   : INDICE TERRE/MER.
  ! PIVEG      : TYPE DE SURFACE (MER, BANQUISE OU GLACIER, TERRE).
  ! PSAB       : POURCENTAGE DE SABLE DANS LA MAILLE.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS DE SORTIE.
  !     --------------------
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (DIAGNOSTIQUE) .
  
  ! PWFC   (*) : TENEUR EN EAU A LA CAPACITE AUX CHAMPS.
  ! PWPMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR PROFOND.
  ! PWSAT      : TENEUR EN EAU A LA SATURATION
  ! PWSMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR DE SURFACE.
  ! PWWILT (*) : TENEUR EN EAU CORRESPONDANT AU POINT DE FLETRISSEMENT.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS IMPLICITES.
  !     ---------------------
  
  ! COMMON/YOMPHY1/
  
  !-----------------------------------------------------------------------
  
  !     Externes.
  !     ---------
  
  !     Methode.
  !     --------
  
  !     Auteur.
  !     -------
  !        99-02, D. Giard, d'apres ACSOL
  
  !     Modifications.
  !     --------------
  !        M.Hamrud      01-Oct-2003 CY28 Cleaning
  !        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
  !        2011-06: M. Jerczynski - some cleaning to meet norms
  !-----------------------------------------------------------------------
  
!$acc routine( ACSOLW_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMPHY1, ONLY: TPHY1
  
  !-----------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TPHY1), INTENT(IN) :: YDPHY1
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(IN) :: PARG
  REAL(KIND=JPRB), INTENT(IN) :: PD2
  REAL(KIND=JPRB), INTENT(IN) :: PLSM
  REAL(KIND=JPRB), INTENT(IN) :: PIVEG
  REAL(KIND=JPRB), INTENT(IN) :: PSAB
  LOGICAL, INTENT(IN) :: LDHMT
  REAL(KIND=JPRB), INTENT(INOUT) :: PWFC
  REAL(KIND=JPRB), INTENT(OUT) :: PWPMX
  REAL(KIND=JPRB), INTENT(INOUT) :: PWSAT
  REAL(KIND=JPRB), INTENT(OUT) :: PWSMX
  REAL(KIND=JPRB), INTENT(OUT) :: PWWILT
  
  !-----------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLON
  
  REAL(KIND=JPRB) :: ZARG
  REAL(KIND=JPRB) :: ZEPS
  REAL(KIND=JPRB) :: ZSAB
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  ZEPS = 1.E-1_JPRB
  
  !     ------------------------------------------------------------------
  !     I - CALCULS REDUITS DANS LE CAS LDHMT=.TRUE.
  !         ----------------------------------------
  !         REDUCED COMPUTATIONS IN CASE LDHMT=.T.
  !         --------------------------------------
  
  ! ** WFC , WWILT **
  
  !DEC$ IVDEP
  IF (PLSM <= 0.5_JPRB) THEN
    PWWILT = 0.0_JPRB
    PWFC = 1.0_JPRB
  ELSE
    ZARG = MAX(ZEPS, PARG)
    PWWILT = YDPHY1%GWWILT*(ZARG**YDPHY1%EWWILT)
    PWFC = YDPHY1%GWFC*(ZARG**YDPHY1%EWFC)
  END IF
  
  IF (.not.LDHMT) THEN
    
    !     ------------------------------------------------------------------
    !     II - CALCUL DES AUTRES CHAMPS
    !          ------------------------
    !          COMPUTING OTHER FIELDS
    !          ----------------------
    
    !DEC$ IVDEP
    
    ! ** CAS DE LA BANQUISE EN MODE CLIMAT - SEA-ICE - (CLIMAT) **
    
    IF (NINT(10._JPRB*PIVEG) == (10*YDPHY1%NTVGLA + 1)) THEN
      PWWILT = 0.0_JPRB
      PWFC = 1.0_JPRB
    END IF
    
    ! ** WSAT **
    
    IF (PLSM <= 0.5_JPRB .or. NINT(PIVEG) == YDPHY1%NTVGLA) THEN
      PWSAT = PWFC
    ELSE
      ZSAB = MAX(ZEPS, PSAB)
      PWSAT = YDPHY1%G1WSAT*ZSAB + YDPHY1%G2WSAT
    END IF
    
    ! **  WSMX , WPMX **
    
    PWPMX = PWSAT*PD2*YDPHY1%GCONV
    PWSMX = PWSAT*YDPHY1%RD1*YDPHY1%GCONV
    
    
  END IF
  
  
END SUBROUTINE ACSOLW_OPENACC
