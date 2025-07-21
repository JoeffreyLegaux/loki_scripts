SUBROUTINE ACCLPH_OPENACC (YDCST, YDPHY0, YDPHY2, KIDIA, KFDIA, KLON, KTDIA, KLEV, PTHETAV, PAPHI, PAPHIF, PU, PV, PTHETAVS,  &
& KCLPH, PCLPH, PVEIN, PUGST, PVGST, YDSTACK)
  !-----------------------------------------------------------------------
  ! - INPUT  2D .
  ! - INPUT  1D .
  ! - OUTPUT 1D .
  
  !**** *ACCLPH * - PHI-PHIS,U,V AU SOMMET DIAGNOSTIQUE DE LA CLP.
  
  !     Sujet.
  !     ------
  
  !     - ROUTINE DE CALCUL ACTIF .
  !          SOMMET DE LA CLA DIAGNOSTIQUE D'APRES AYOTTE BLM 1996.
  !          VENTS AU SOMMET DE LA COUCHE LIMITE PLANETAIRE,
  
  !          TOP OF PBL DIAGNOSED VIA FOLLOWING AYOTTE BLM 1996.
  !          WIND AT THE TOP OF THE PLANETARY BOUNDARY LAYER
  
  !**   Interface.
  !     ----------
  !        *CALL* *ACCLPH*
  
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
  ! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
  ! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 2D (0:KLEV) .
  
  ! PAPHI      : GEOPOTENTIEL AUX DEMI-NIVEAUX.
  
  ! - 2D (1:KLEV) .
  
  ! PAPHIF     : GEOPOTENTIEL AUX NIVEAUX.
  ! PTHETAV    : THETAV.
  ! PU         : COMPOSANTE EN X DU VENT.
  ! PV         : COMPOSANTE EN Y DU VENT.
  
  ! - 1D (DIAGNOSTIQUE) .
  
  ! PTHETAVS   : THETAV EN SURFACE.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS DE SORTIE.
  !     --------------------
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (DIAGNOSTIQUE) .
  
  ! KCLPH      : INDICE DU NIVEAU MODELE CONTENANT LE SOMMET DE LA CLP.
  ! PCLPH      : ALTITUDE EN METRES (PAR RAPPORT A LA SURFACE) AU SOMMET DE LA CLP.
  ! PVEIN      : INDEX DE VENTILATION DANS LA CLP.
  ! PUGST      : U COMPONENT OF GUSTS (DIAGNOSTIC)
  ! PVGST      : V COMPONENT OF GUSTS (DIAGNOSTIC)
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS IMPLICITES.
  !     ---------------------
  
  ! COMMON/YOMCST /
  ! COMMON/YOMPHY0 /
  
  !-----------------------------------------------------------------------
  
  !     Externes.
  !     ---------
  
  !     Methode.
  !     --------
  
  !     Auteurs.
  !     -------
  !        2002-08, J.M. Piriou et Jean-Francois Geleyn.
  
  !     Modifications.
  !     --------------
  !     2003-04. M. Bellus: wind gusts in case of LRAFTUR=.F.
  !        M.Hamrud      01-Oct-2003 CY28 Cleaning
  !     03-03-2006  R.Brozkova Shear linked convection after M. Tudor
  !                          and harmonisation with Troen-Mahrt method
  !     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !-----------------------------------------------------------------------
  
!$acc routine( ACCLPH_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOMPHY0, ONLY: TPHY0
  USE YOMPHY2, ONLY: TPHY2
  
  !-----------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TPHY0), INTENT(IN) :: YDPHY0
  TYPE(TPHY2), INTENT(IN) :: YDPHY2
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  REAL(KIND=JPRB), INTENT(IN) :: PTHETAV(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHI(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PV(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTHETAVS(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCLPH(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PCLPH(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PVEIN(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PUGST(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PVGST(KLON)
  
  !-----------------------------------------------------------------------
  
  REAL(KIND=JPRB) :: ZINT
  REAL(KIND=JPRB) :: ZX_PREV
  temp (REAL (KIND=JPRB), ZTHETAV, (KLON, KLEV + 1))
  temp (REAL (KIND=JPRB), ZPHI, (KLON, KLEV + 1))
  REAL(KIND=JPRB) :: ZLCH_PREV
  temp (REAL (KIND=JPRB), ZU, (KLON, KLEV + 1))
  temp (REAL (KIND=JPRB), ZV, (KLON, KLEV + 1))
  temp (REAL (KIND=JPRB), ZTHETAVS, (KLON, KLEV + 1))
  
  ! ZCLPU :: VENT "U" AU SOMMET DE LA CLP
  ! ZCLPV :: VENT "V" AU SOMMET DE LA CLP
  REAL(KIND=JPRB) :: ZCLPU
  REAL(KIND=JPRB) :: ZCLPV
  
  ! ZVITM   :: VITESSE DU VENT MOYEN VERTICAL DANS LA CLP
  REAL(KIND=JPRB) :: ZVITM
  
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: IBIN
  INTEGER(KIND=JPIM) :: IC
  INTEGER(KIND=JPIM) :: INDI
  
  REAL(KIND=JPRB) :: ZPBLHK1
  REAL(KIND=JPRB) :: ZX
  REAL(KIND=JPRB) :: ZPSI
  REAL(KIND=JPRB) :: ZEPS
  REAL(KIND=JPRB) :: ZBIG
  REAL(KIND=JPRB) :: ZDELTA
  REAL(KIND=JPRB) :: ZARG1
  REAL(KIND=JPRB) :: ZARG2
  REAL(KIND=JPRB) :: ZBIN
  REAL(KIND=JPRB) :: ZLCH
  REAL(KIND=JPRB) :: ZVAL_LCH
  REAL(KIND=JPRB) :: ZLOG2
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZTHETAV) == 8) THEN
    alloc8 (ZTHETAV)
  ELSE
    IF (KIND (ZTHETAV) == 4) THEN
      alloc4 (ZTHETAV)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZPHI) == 8) THEN
    alloc8 (ZPHI)
  ELSE
    IF (KIND (ZPHI) == 4) THEN
      alloc4 (ZPHI)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZU) == 8) THEN
    alloc8 (ZU)
  ELSE
    IF (KIND (ZU) == 4) THEN
      alloc4 (ZU)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZV) == 8) THEN
    alloc8 (ZV)
  ELSE
    IF (KIND (ZV) == 4) THEN
      alloc4 (ZV)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZTHETAVS) == 8) THEN
    alloc8 (ZTHETAVS)
  ELSE
    IF (KIND (ZTHETAVS) == 4) THEN
      alloc4 (ZTHETAVS)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KIDIA
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  !     PREPARATION.
  
  !         INITIAL SETTING.
  
  ZPBLHK1 = YDPHY0%GPBLHK0 / YDPHY0%GPBLHRA
  ZEPS = 0.001_JPRB
  ZBIG = 30._JPRB
  ZLOG2 = LOG(2.0_JPRB)
  
  ZARG1 = -YDPHY0%GPBLHK0 / ZPBLHK1
  ZVAL_LCH = MAX(LOG(COSH(MIN(ABS(ZARG1), ZBIG))), ABS(ZARG1) - ZLOG2)
  PCLPH(JLON) = 0.0_JPRB
  PVEIN(JLON) = 0.0_JPRB
  ZCLPU = 0.0_JPRB
  ZCLPV = 0.0_JPRB
  ZVITM = SQRT(PU(JLON, KLEV)**2 + PV(JLON, KLEV)**2)
  KCLPH(JLON) = KLEV
  ZINT = 0.0_JPRB
  ZX_PREV = 0.0_JPRB
  ZLCH_PREV = ZVAL_LCH
  ZTHETAV(JLON, KLEV + 1) = MAX(PTHETAVS(JLON), PTHETAV(JLON, KLEV))
  ZTHETAVS(JLON, KLEV + 1) = ZTHETAV(JLON, KLEV + 1)
  ZPHI(JLON, KLEV + 1) = PAPHI(JLON, KLEV)
  ZU(JLON, KLEV + 1) = 0.0_JPRB
  ZV(JLON, KLEV + 1) = 0.0_JPRB
  DO JLEV=KTDIA,KLEV
    ZTHETAV(JLON, JLEV) = PTHETAV(JLON, JLEV)
    ZPHI(JLON, JLEV) = PAPHIF(JLON, JLEV)
    ZU(JLON, JLEV) = PU(JLON, JLEV)
    ZV(JLON, JLEV) = PV(JLON, JLEV)
  END DO
  !MT shear linked convection modification
  !cdir unroll=8
  DO JLEV=KLEV,KTDIA,-1
    ZTHETAVS(JLON, JLEV) = ZTHETAVS(JLON, JLEV + 1) + ZTHETAV(JLON, JLEV) - ZTHETAV(JLON, JLEV + 1) - 0.5_JPRB*(ZTHETAV(JLON,  &
    & JLEV) + ZTHETAV(JLON, JLEV + 1))*((ZU(JLON, JLEV) - ZU(JLON, JLEV + 1))**2 + (ZV(JLON, JLEV) - ZV(JLON, JLEV + 1))**2) /  &
    & (ZPHI(JLON, JLEV) - ZPHI(JLON, JLEV + 1))
  END DO
  
  !-------------------------------------------------
  ! INTEGRALE ASCENDANTE.
  ! UPWARD INTEGRATIONS.
  !-------------------------------------------------
  
  DO JLEV=KLEV,KTDIA,-1
    !DEC$ IVDEP
    !
    !-------------------------------------------------
    ! INTEGRALE DE THETAV.
    ! THETAV INTEGRAL.
    !-------------------------------------------------
    !
    ZINT = ZINT + (ZPHI(JLON, JLEV) - ZPHI(JLON, JLEV + 1))*0.5_JPRB*(ZTHETAVS(JLON, JLEV) + ZTHETAVS(JLON, JLEV + 1))
    !
    !-------------------------------------------------
    ! ECART THETAV' ENTRE THETAV DU NIVEAU COURANT
    ! ET LA MOYENNE DE THETAV ENTRE LA SURFACE ET LE NIVEAU COURANT.
    ! COMPUTE THETAV', DIFFERENCE BETWEEN CURRENT THETAV VALUE
    ! AND ITS INTEGRAL BETWEEN SURFACE AND CURRENT LEVEL.
    !-------------------------------------------------
    !
    ZX = ZTHETAVS(JLON, JLEV) - ZINT / (ZPHI(JLON, JLEV) - ZPHI(JLON, KLEV + 1))
    !
    !-------------------------------------------------
    ! DIFFERENCE ENTRE THETAV' DU NIVEAU COURANT ET CELUI DU NIVEAU DU DESSOUS.
    ! DIFFERENCE BETWEEN CURRENT LEVEL THETAV' AND LEVEL BELOW.
    !-------------------------------------------------
    !
    ZDELTA = ZX - ZX_PREV
    !
    !-------------------------------------------------
    ! SECURITE EN DIVISION: ZDELTA VA SERVIR AU DENOMINATEUR.
    ! ON LE BORNE POUR EVITER LA DIVISION PAR ZERO, EN CONSERVANT SON SIGNE.
    ! SECURE DIVISION: ZDELTA IS AT THE DENOMINATOR.
    ! ONE PREVENTS A ZERO VALUE, WHILE KEEPING THE SAME SIGN.
    !-------------------------------------------------
    !
    ZBIN = MAX(0.0_JPRB, SIGN(1.0_JPRB, ABS(ZDELTA) - ZEPS))
    ZDELTA = ZBIN*ZDELTA + (1.0_JPRB - ZBIN)*ZEPS*SIGN(1.0_JPRB, ZDELTA)
    ZX = ZX_PREV + ZDELTA
    !
    !-------------------------------------------------
    ! SECURITE EN DEBORDEMENT DE L'EXPONENTIELLE: ON UTILISE
    ! LA FONCTION LOG(COSH) POUR LES ARGUMENTS INFERIEURS A ZBIG, ET SON ASYMPTOTE AU-DELA.
    ! SECURE EXPONENTIAL OVERFLOWS: ONE USES
    ! THE LOG(COSH) FUNCTION FOR ARGUMENTS BELOW ZBIG, AND ITS ASYMPTOT ELSE CASE.
    !-------------------------------------------------
    !
    ZARG2 = (ZX - YDPHY0%GPBLHK0) / ZPBLHK1
    ZLCH = MAX(LOG(COSH(MIN(ABS(ZARG2), ZBIG))), ABS(ZARG2) - ZLOG2)
    ZPSI = 0.5_JPRB*(1.0_JPRB - ZPBLHK1 / ZDELTA*(ZLCH - ZLCH_PREV))
    !
    !-------------------------------------------------
    ! INTEGRALE DE L'EPAISSEUR GEOPOTENTIELLE, AVEC POUR POIDS PSI.
    ! INTEGRATE GEOPOTENTIAL DEPTH, WITH PSI WEIGHTING FUNCTION.
    !-------------------------------------------------
    !
    PCLPH(JLON) = PCLPH(JLON) + (ZPHI(JLON, JLEV) - ZPHI(JLON, JLEV + 1))*ZPSI
    !
    !-------------------------------------------------
    ! ON REMPLACE LES VALEURS DU NOUVEAU INFERIEUR PAR CELLES DU NIVEAU COURANT.
    ! ONE REPLACES VALUES FROM LEVEL BELOW BY THOSE OF CURRENT LEVEL.
    !-------------------------------------------------
    !
    ZX_PREV = ZX
    ZLCH_PREV = ZLCH
  END DO
  
  PCLPH(JLON) = MAX(PCLPH(JLON), PAPHIF(JLON, KLEV) - PAPHI(JLON, KLEV))
  
  !-------------------------------------------------
  ! SECONDE BOUCLE VERTICALE: VENT AU SOMMET DE LA CLA.
  ! SECOND VERTICAL LOOP: WIND AT THE TOP OF PBL.
  !-------------------------------------------------
  
  !cdir unroll=8
  DO JLEV=KLEV,KTDIA,-1
    ZPSI = MAX(0.0_JPRB, MIN(1.0_JPRB, (PCLPH(JLON) + PAPHI(JLON, KLEV) - ZPHI(JLON, JLEV + 1)) / (ZPHI(JLON, JLEV) - ZPHI(JLON,  &
    & JLEV + 1))))
    ZCLPU = ZPSI*PU(JLON, JLEV) + (1.0_JPRB - ZPSI)*ZCLPU
    ZCLPV = ZPSI*PV(JLON, JLEV) + (1.0_JPRB - ZPSI)*ZCLPV
    IBIN = NINT(ZPSI)
    KCLPH(JLON) = IBIN*JLEV + (1 - IBIN)*KCLPH(JLON)
  END DO
  
  !cdir unroll=8
  IC = 1
  DO JLEV=KLEV - 1,KTDIA,-1
    INDI = INT(MAX(0.0_JPRB, SIGN(1.0_JPRB, PCLPH(JLON) + PAPHI(JLON, KLEV) - ZPHI(JLON, JLEV + 1))))
    ZVITM = ZVITM + INDI*SQRT(PU(JLON, JLEV)**2 + PV(JLON, JLEV)**2)
    IC = IC + INDI
  END DO
  
  PVEIN(JLON) = PCLPH(JLON)*(ZVITM / FLOAT(IC)) / YDCST%RG
  !-------------------------------------------------
  ! WIND GUSTS IN CASE OF LRAFTUR=.F.
  !-------------------------------------------------
  
  IF (.not.YDPHY2%LRAFTUR) THEN
    PUGST(JLON) = ZCLPU
    PVGST(JLON) = ZCLPV
  END IF
  
  !-------------------------------------------------
  ! GEOPOTENTIELS > ALTITUDES.
  ! GEOPOTENTIALS > ALTITUDES.
  !-------------------------------------------------
  
  PCLPH(JLON) = PCLPH(JLON) / YDCST%RG
  
  !-----------------------------------------------------------------------
END SUBROUTINE ACCLPH_OPENACC
