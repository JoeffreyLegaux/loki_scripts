SUBROUTINE ACVEG_OPENACC (YDPHY, YDPHY1, KIDIA, KFDIA, KLON, KLEV, PFRSO, PQ, PQSAT, PT, PD2, PLSM, PIVEG, PLAI, PNEIJ, PVEG,  &
& PRSMIN, PCHROV, PGWDCS, PWFC, PWL, PWWILT, PWP, PQSATS, PHQ, PHTR, PHU, PHV, PWLMX, YDSTACK)
  !-----------------------------------------------------------------------
  ! - INPUT  2D .
  ! - INPUT  1D .
  ! - OUTPUT 1D .
  
  !**** *ACVEG   * - DETERMINATION DES CARACTERISTIQUES DE LA VEGETATION.
  
  !     Sujet.
  !     ------
  
  !     - ROUTINE DE CALCUL ACTIF .
  !       DETERMINATION DES CARACTERISTIQUES DE SURFACE CONCERNANT LA VEGE
  !       CALCUL DE WLMX , HV , HTR , HU , HQ.
  !     -  VERSION ANGLAISE.
  
  !**   Interface.
  !     ----------
  !        *CALL* *ACVEG*
  
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
  ! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 2D (0:KLEV) .
  
  ! PFRSO      : FLUX DE RAYONNEMENT SOLAIRE.
  
  ! - 2D (1:KLEV) .
  
  ! PQ         : HUMIDITE SPECIFIQUE DE LA VAPEUR D'EAU.
  ! PQSAT      : HUMIDITE SPECIFIQUE SATURANTE.
  ! PT         : TEMPERATURE.
  
  ! - 1D (PRONOSTIQUE) .
  
  ! PWL        : CONTENU EN EAU DU RESERVOIR D'INTERCEPTION.
  ! PWP        : CONTENU EN EAU DU RESERVOIR TOTAL (SOL).
  
  ! - 1D (GEOGRAPHIQUE) .
  
  ! PD2        : PROFONDEUR DU SOL.
  ! PLSM       : INDICE TERRE/MER.
  ! PIVEG      : TABLEAU PERMETTANT D'IDENTIFIER LE TYPE DE VEGETATION DAN
  !              CHAQUE MAILLE DU DOMAINE SIMULE.
  ! PLAI       : INDICE FOLIAIRE.
  ! PRSMIN     : RESISTANCE STOMATIQUE MINIMALE.
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (DIAGNOSTIQUE) .
  
  ! PNEIJ      : PROPORTION DE SOL ENNEIGE.
  ! PVEG       : FRACTION DE VEGETATION APPARENTE.
  ! PCHROV     : COEFFICIENT D'ECHANGE EN SURFACE POUR T ET Q
  !              RENORME EN DENSITE FOIS VITESSE.
  ! PGWDCS     : DENSITE EN SURFACE
  ! PWFC       : TENEUR EN EAU A LA CAPACITE AUX CHAMPS.
  ! PWWILT     : TENEUR EN EAU AU POINT DE FLETRISSEMENT.
  ! PQSATS     : HUMIDITE SPECIFIQUE DE SATURATION EN SURFACE.
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS DE SORTIE.
  !     --------------------
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (DIAGNOSTIQUE) .
  
  ! PHQ        : POIDS DE L'HUMIDITE DE L'AIR DANS L'HUMIDITE DE SURFACE.
  ! PHU        : POIDS DE L'HUMIDITE SATURANTE DANS L'HUMIDITE DE SURFACE.
  ! PHV        : RESISTANCE A L'EVAPOTRANSPIRATION DU COUVERT VEGETAL.
  ! PHTR       : RESISTANCE A LA TRANSPIRATION DU COUVERT VEGETAL.
  ! PWLMX      : CONTENU EN EAU MAXIMAL DU RESERVOIR D'INTERCEPTION.
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS IMPLICITES.
  !     ---------------------
  
  !-----------------------------------------------------------------------
  
  !     Externes.
  !     ---------
  
  !     Methode.
  !     --------
  
  !     Auteur.
  !     -------
  !      91-12, J. Noilhan.
  
  !     Modifications.
  !     --------------
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      K. Yessad (Jul 2009): remove CDLOCK + some cleanings
  !-----------------------------------------------------------------------
  
!$acc routine( ACVEG_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMPHY, ONLY: TPHY
  USE YOMPHY1, ONLY: TPHY1
  
  !-----------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TPHY), INTENT(IN) :: YDPHY
  TYPE(TPHY1), INTENT(IN) :: YDPHY1
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(IN) :: PFRSO(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSAT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PD2(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PLSM(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PIVEG(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PLAI(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PNEIJ(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PVEG(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PRSMIN(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PCHROV(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PGWDCS(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PWFC(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PWL(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PWWILT(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PWP(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PQSATS(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PHQ(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PHTR(KLON)
  REAL(KIND=JPRB), INTENT(INOUT) :: PHU(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PHV(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PWLMX(KLON)
  
  !-----------------------------------------------------------------------
  
  REAL(KIND=JPRB) :: ZDELTA
  REAL(KIND=JPRB) :: ZRSTO
  
  INTEGER(KIND=JPIM) :: IVEG
  INTEGER(KIND=JPIM) :: JLON
  
  REAL(KIND=JPRB) :: ZEPS
  REAL(KIND=JPRB) :: ZF
  REAL(KIND=JPRB) :: ZF1
  REAL(KIND=JPRB) :: ZF2
  REAL(KIND=JPRB) :: ZF3
  REAL(KIND=JPRB) :: ZF4
  REAL(KIND=JPRB) :: ZLAI
  REAL(KIND=JPRB) :: ZWL
  REAL(KIND=JPRB) :: ZWP
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  ZEPS = 1.E-5_JPRB
  
  !     ------------------------------------------------------------------
  !     I - CALCUL DE WLMX.
  !         ---------------
  
  !DEC$ IVDEP
  
  ZLAI = MAX(ZEPS, PLAI(JLON))
  PWLMX(JLON) = YDPHY1%GWLMX*ZLAI*PVEG(JLON)
  ZWL = MAX(0.0_JPRB, MIN(PWLMX(JLON), PWL(JLON)))
  ZDELTA = 0.0_JPRB
  IF (ZWL > 0.0_JPRB) THEN
    ZDELTA = (ZWL / PWLMX(JLON))**YDPHY1%GWLEX
  END IF
  
  
  !     ------------------------------------------------------------------
  !     II- CALCUL DE LA RESISTANCE STOMATIQUE.
  !         -----------------------------------
  
  
  ZLAI = MAX(ZEPS, PLAI(JLON))
  IVEG = NINT(PIVEG(JLON))
  ZF = YDPHY1%GF1*PFRSO(JLON, KLEV) / YDPHY1%RGL(IVEG)*2.0_JPRB / ZLAI
  ZF1 = (1.0_JPRB + ZF) / (ZF + PRSMIN(JLON) / YDPHY1%RSMAX)
  ZWP = PWP(JLON) / (PD2(JLON)*YDPHY1%GCONV)
  ZF2 = MIN(1.0_JPRB, MAX(ZEPS, (ZWP - PWWILT(JLON)) / (PWFC(JLON) - PWWILT(JLON))))
  ZF3 = MAX(ZEPS, (1.0_JPRB - YDPHY1%GF3(IVEG)*(PQSAT(JLON, KLEV) - PQ(JLON, KLEV))))
  ZF4 = MAX(ZEPS, (1.0_JPRB - YDPHY1%GF4(IVEG)*(YDPHY1%TREF4(IVEG) - PT(JLON, KLEV))**2))
  ZRSTO = PRSMIN(JLON) / ZLAI*ZF1 / (ZF2*ZF3*ZF4)
  ZRSTO = MIN(ZRSTO, YDPHY1%RSMAX)
  
  
  !     ------------------------------------------------------------------
  !     III - CALCUL DES HUMIDITES HV ET HTR.
  !           -------------------------------
  !           COMPUTING OF HUMIDITY HV AND HTR.
  !           ---------------------------------
  
  
  PHV(JLON) = 1.0_JPRB - MAX(0.0_JPRB, SIGN(1.0_JPRB, PQSATS(JLON) - PQ(JLON, KLEV)))*PLSM(JLON)*ZRSTO*PCHROV(JLON)*(1.0_JPRB -  &
  & ZDELTA) / (PGWDCS(JLON) + ZRSTO*PCHROV(JLON))
  
  PHTR(JLON) = (PHV(JLON) - ZDELTA)*PLSM(JLON)
  IF (YDPHY%LNEIGE) THEN
    PHV(JLON) = (1.0_JPRB - PNEIJ(JLON))*PHV(JLON) + PNEIJ(JLON)
  END IF
  
  PHU(JLON) = PVEG(JLON)*PHV(JLON) + (1.0_JPRB - PVEG(JLON))*PHU(JLON)
  PHQ(JLON) = PVEG(JLON)*(1.0_JPRB - PHV(JLON))
  
  
  !-----------------------------------------------------------------------
END SUBROUTINE ACVEG_OPENACC
