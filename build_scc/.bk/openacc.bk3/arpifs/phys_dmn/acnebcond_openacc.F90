SUBROUTINE ACNEBCOND_OPENACC (YDCST, YDRIP, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, LDREDPR, YDSTA, PAPHI, PAPHIF, PAPRSF,  &
& PCP, PR, PDELP, PRH, PBLH, PQ, PQI, PQL, PQW, PT, PNCV, PGM, PTS, PQCS, PNEBCOND, PHCRICS, PRHOUT, PQSATS, PRMF, PQCS0,  &
& PNEBS0, YDSTACK)
  !-----------------------------------------------------------------------
  ! - INPUT  1D VERTICAL.
  !-----------------------------------------------------------------------
  ! - INPUT  2D .
  ! - INPUT-OUTPUT  2D .
  ! - INPUT  1D .
  ! - OUTPUT 2D .
  
  !**** *ACNEBCOND * - CALCUL DE NEBULOSITE ET HUMIDITE CRITIQUE POUR LA
  !                    CONDENSATION RESOLUE
  
  !     Sujet.
  !     ------
  !     - ROUTINE DE CALCUL ACTIF .
  !       CALCUL DE NEBULOSITE ET HUMIDITE CRITIQUE POUR LA CONDENSATION RESOLUE.
  !     - COMPUTATION OF CLOUDINESS AND CRITICAL HUMIDITY FOR THE RESOLVED
  !              CONDENSATION .
  
  !**   Interface.
  !     ----------
  !        *CALL* *ACNEBCOND*
  
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
  
  ! LOGICAL KEY
  ! LDREDPR    : T FOR REDUCED PROTECTION, F for FULL PROTECTION OF PNCV
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (1:KLEV) .
  
  ! PHUC       : PROFIL DE BASE POUR L'HUMIDITE RELATIVE CRITIQUE.
  ! PVETAF     : LA COORDONNEE VERTICALE AUX COUCHES.
  
  ! - 2D (0:KLEV) .
  
  ! PAPHI      : HALF LEVEL GEOPOTENTIAL.
  
  ! - 2D (1:KLEV) .
  
  ! PAPHIF     : FULL LEVEL GEOPOTENTIAL.
  ! PAPRSF     : PRESSION AUX NIVEAUX PLEIN.
  ! PCP        : CHALEUR MASSIQUE A PRESSION CONSTANTE.
  ! PR         : CONSTANTE DE L'AIR HUMIDE.
  ! PDELP      : EPAISSEUR EN PRESSION DE LA COUCHE.
  ! PQ         : HUMIDITE SPECIFIQUE DE LA VAPEUR D'EAU.
  ! PQW        : HUMIDITE SPECIFIQUE DU THERMOMETRE MOUILLE.
  ! PT         : TEMPERATURE.
  ! PNCV       : CONVECTIVE CLOUDINESS (HISTORIC, INPUT)
  !               -> EFFECTIVELY PROTECTED FRACTION (OUTPUT)
  
  
  ! ADDITIONAL PROGNOSTIC VARIABLES
  
  ! PQI        : RATIO OF SUSPENDED ICE.
  ! PQL        : RATIO OF SUSPENDED LIQUID WATER.
  
  ! - 1D .
  
  ! PGM        : FACTEUR D'ECHELLE.
  ! PTS        : TEMPERATURE DE SURFACE.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS DE SORTIE.
  !     --------------------
  
  ! - 2D (KLEV) .
  
  ! PQCS       : CONTENU STRATIFORME EN CONDENSAT NUAGEUX (option SMITH).
  ! PNEBCOND   : NEBULOSITE POUR LA CONDENSATION RESOLUE.
  ! PHCRICS    : HUMIDITE CRITIQUE RELATIVE POUR LA CONDENSATION RESOLUE.
  ! PQSATS     : HUMIDITE SPECIFIQUE DE SATURATION.
  ! PRMF       : PROPORTION DE LA GLACE.
  ! PQCS0      : CONTENU STRATIFORME EN CONDENSAT NUAGEUX POUR RAYONNEMENT.
  ! PNEBS0     : NEBULOSITE PARTIELLE STRATIFORME POUR RAYONNEMENT.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS IMPLICITES.
  !     ---------------------
  
  ! COMMON/YOMPHY /
  ! COMMON/YOMCST /
  ! COMMON/YOMPHY0/
  ! COMMON/YOMPHY2/
  
  !-----------------------------------------------------------------------
  
  !     Externes.
  !     ---------
  
  !     Methode.
  !     --------
  
  !     Auteur.
  !     -------
  !        07-02, R. Brozkova.
  
  !     Modifications.
  !     --------------
  !        09-10, L. Bengtsson: Rasch Kristjansson (RK) scheme
  !        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
  !        2011-06: M. Jerczynski - some cleaning to meet norms
  !        K-I Ivarsson 2011-05: Some modifications to RK-scheme
  !        2012-02: R. Brozkova: Smith scheme and output for radiation.
  !        2012-06: R. Brozkova: XR scheme - better RHcrit dependency on dx.
  !        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
  !        R. Brozkova (Oct 2014) parameters QXRAL_ADJ and ADJTAU for XR option.
  !        L. Gerard (Apr 2016) XR: reduced protection LDREDPR still allowing
  !                                 condensation; simplified PQCS calculation
  !                                 (remove unnecessary phase separation).
  !                             SMG: Fix uninitialized ZLV, ZLS.
  !        P. Marguinaud (Oct 2016) : Port to single precision
  !        R. Brozkova (Sep 2018): Fixes in thermodynamic adjustment - deep
  !                                convective condensates protection.
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !-------------------------------------------------------------------------
  
!$acc routine( ACNEBCOND_OPENACC ) seq
  
  USE MODEL_PHYSICS_MF_MOD, ONLY: MODEL_PHYSICS_MF_TYPE
  USE PARKIND1, ONLY: JPIM, JPRB, JPRD
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOMRIP, ONLY: TRIP
  USE YOMSTA, ONLY: TSTA
  !-----------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(MODEL_PHYSICS_MF_TYPE), INTENT(IN) :: YDML_PHY_MF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TRIP), INTENT(IN) :: YDRIP
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  LOGICAL, INTENT(IN) :: LDREDPR
  TYPE(TSTA), INTENT(IN) :: YDSTA
  REAL(KIND=JPRB), INTENT(IN) :: PAPHI(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPRSF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PCP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PR(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDELP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQI(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQW(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PNCV(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGM(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTS(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PRH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQCS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PNEBCOND(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PHCRICS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQSATS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PRMF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PBLH(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PRHOUT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQCS0(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PNEBS0(KLON, KLEV)
  
  temp (REAL (KIND=JPRB), ZQTOT, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZRHL, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZSHALTEMP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZHCUTMP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZLSTMP, (KLON, KLEV))
  
  REAL(KIND=JPRB) :: ZDC
  REAL(KIND=JPRB) :: ZAW
  REAL(KIND=JPRB) :: ZX0
  REAL(KIND=JPRB) :: ZSITER
  REAL(KIND=JPRB) :: ZNEB1
  REAL(KIND=JPRB) :: ZA
  REAL(KIND=JPRB) :: ZB
  REAL(KIND=JPRB) :: ZC
  REAL(KIND=JPRB) :: ZRMF
  REAL(KIND=JPRB) :: ZS1
  REAL(KIND=JPRB) :: ZS2
  REAL(KIND=JPRB) :: ZMESHEXP
  
  INTEGER(KIND=JPIM) :: ILIQ
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: IITER
  INTEGER(KIND=JPIM) :: JITER
  INTEGER(KIND=JPIM) :: ILEV
  INTEGER(KIND=JPIM) :: JLEV2
  LOGICAL :: LLWHOLE
  
  REAL(KIND=JPRB) :: ZEPS1
  REAL(KIND=JPRB) :: ZEPS2
  REAL(KIND=JPRB) :: ZEPS3
  REAL(KIND=JPRB) :: ZEPS4
  REAL(KIND=JPRB) :: ZEPS5
  REAL(KIND=JPRB) :: ZESAT
  REAL(KIND=JPRB) :: ZESATL
  REAL(KIND=JPRB) :: ZESP
  REAL(KIND=JPRB) :: ZNI
  REAL(KIND=JPRB) :: ZNPI
  REAL(KIND=JPRB) :: ZDNII
  REAL(KIND=JPRB) :: ZFAC
  REAL(KIND=JPRB) :: ZXN
  REAL(KIND=JPRB) :: ZDXN
  REAL(KIND=JPRB) :: ZQXRAL
  REAL(KIND=JPRB) :: ZALFI
  REAL(KIND=JPRB) :: ZTKLEV
  REAL(KIND=JPRB) :: ZZ
  REAL(KIND=JPRB) :: ZLEN0
  REAL(KIND=JPRB) :: ZLESEFR
  REAL(KIND=JPRB) :: ZLESEFS
  REAL(KIND=JPRB) :: ZMESH
  REAL(KIND=JPRB) :: ZQCST
  REAL(KIND=JPRB) :: ZRHC
  REAL(KIND=JPRB) :: ZRATQ
  REAL(KIND=JPRB) :: ZNEB0
  REAL(KIND=JPRB) :: ZRAT2
  REAL(KIND=JPRB) :: ZDNEB
  REAL(KIND=JPRB) :: ZDCRI
  REAL(KIND=JPRB) :: ZFACT
  REAL(KIND=JPRB) :: ZFACT0
  REAL(KIND=JPRB) :: ZFACTA
  REAL(KIND=JPRB) :: ZFACTB
  REAL(KIND=JPRB) :: ZFACTC
  REAL(KIND=JPRB) :: ZDS
  REAL(KIND=JPRB) :: ZDELTAP
  REAL(KIND=JPRB) :: ZNCVNEB
  REAL(KIND=JPRB) :: ZNIM
  REAL(KIND=JPRB) :: ZZRHC
  REAL(KIND=JPRB) :: ZTEST
  REAL(KIND=JPRB) :: ZSIGMA
  REAL(KIND=JPRB) :: ZSIGMASTAB
  REAL(KIND=JPRB) :: ZEPSILO
  REAL(KIND=JPRB) :: ZCLD
  REAL(KIND=JPRB) :: ZRHD
  REAL(KIND=JPRB) :: ZRHIN
  REAL(KIND=JPRB) :: ZDRHDZ
  REAL(KIND=JPRB) :: ZSIGMAZ
  REAL(KIND=JPRB) :: ZONEMRHCRIT
  REAL(KIND=JPRB) :: ZSTIMEINC
  REAL(KIND=JPRB) :: ZRHDIF
  REAL(KIND=JPRB) :: ZDESDTM
  REAL(KIND=JPRB) :: ZESN
  REAL(KIND=JPRB) :: ZDENOM
  REAL(KIND=JPRB) :: ZNEBT
  REAL(KIND=JPRB) :: ZRNEB
  REAL(KIND=JPRB) :: ZIRNEB
  REAL(KIND=JPRB) :: ZQVCS
  REAL(KIND=JPRB) :: ZQ
  REAL(KIND=JPRB) :: ZQL
  REAL(KIND=JPRB) :: ZQI
  REAL(KIND=JPRB) :: ZALF
  REAL(KIND=JPRB) :: ZNEBM1
  REAL(KIND=JPRB) :: ZQVN
  REAL(KIND=JPRB) :: ZDQVN
  REAL(KIND=JPRB) :: ZRMLOC
  REAL(KIND=JPRB) :: ZTESTB
  REAL(KIND=JPRB) :: ZCONL
  REAL(KIND=JPRB) :: ZCONI
  REAL(KIND=JPRB) :: ZADJL
  REAL(KIND=JPRB) :: ZADJI
  REAL(KIND=JPRB) :: ZWEIGHT
  REAL(KIND=JPRB) :: ZFRACON
  REAL(KIND=JPRB) :: ZSQRT6
  REAL(KIND=JPRB) :: ZRDSRV
  REAL(KIND=JPRB) :: ZTLIQ
  REAL(KIND=JPRB) :: ZSIGS
  REAL(KIND=JPRB) :: ZLEFF
  REAL(KIND=JPRB) :: ZLV
  REAL(KIND=JPRB) :: ZLS
  REAL(KIND=JPRB) :: ZZDIST2
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  !-----------------------------------------------------------------------
  
#include "fctdoi.func.h"
#include "fcttrm.func.h"
  
  !-----------------------------------------------------------------------
  
#include "acnebsm_openacc.intfb.h"
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZQTOT) == 8) THEN
    alloc8 (ZQTOT)
  ELSE
    IF (KIND (ZQTOT) == 4) THEN
      alloc4 (ZQTOT)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZRHL) == 8) THEN
    alloc8 (ZRHL)
  ELSE
    IF (KIND (ZRHL) == 4) THEN
      alloc4 (ZRHL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZS) == 8) THEN
    alloc8 (ZS)
  ELSE
    IF (KIND (ZS) == 4) THEN
      alloc4 (ZS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZSHALTEMP) == 8) THEN
    alloc8 (ZSHALTEMP)
  ELSE
    IF (KIND (ZSHALTEMP) == 4) THEN
      alloc4 (ZSHALTEMP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZHCUTMP) == 8) THEN
    alloc8 (ZHCUTMP)
  ELSE
    IF (KIND (ZHCUTMP) == 4) THEN
      alloc4 (ZHCUTMP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZLSTMP) == 8) THEN
    alloc8 (ZLSTMP)
  ELSE
    IF (KIND (ZLSTMP) == 4) THEN
      alloc4 (ZLSTMP)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KIDIA
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  !*
  !     ------------------------------------------------------------------
  !     I - CALCUL DE PARAMETRES DERIVES ET CONSTANTES DE SECURITE.
  
  !         COMPUTATION OF DERIVED PARAMETERS AND SECURITY CONSTANTS.
  
  IF (JPRB == JPRD) THEN
    ZEPS1 = 1.E-14_JPRB
    ZEPS5 = 1.E-12_JPRB
  ELSE
    ZEPS1 = 1.E-06_JPRB
    ZEPS5 = 1.E-06_JPRB
  END IF
  ZEPS2 = 1.E-02_JPRB
  ZEPS3 = 1.E-20_JPRB
  ZEPS4 = 1.E-10_JPRB
  ZSQRT6 = SQRT(6._JPRB)
  ZRDSRV = YDCST%RD / YDCST%RV
  ZWEIGHT = 1.0_JPRB - EXP(-YDML_PHY_MF%YRPHY2%TSPHY / YDML_PHY_MF%YRPHY0%ADJTAU)
  
  ZQXRAL = YDML_PHY_MF%YRPHY0%QXRAL_ADJ
  ZLESEFR = 1._JPRB / YDML_PHY_MF%YRPHY0%SCLESPR**YDML_PHY_MF%YRPHY0%RHCEXPDX
  ZLESEFS = 1._JPRB / YDML_PHY_MF%YRPHY0%SCLESPS**YDML_PHY_MF%YRPHY0%RHCEXPDX
  IITER = 3
  
  !*
  !     ------------------------------------------------------------------
  !     II - CALCULS PRELIMINAIRES DE LA PROPORTION DE LA GLACE.
  
  !          PRELIMINARY COMPUTATIONS OF ICE PROPORTION.
  
  
  DO JLEV=KTDIA,KLEV
    !DEC$ IVDEP
    PRMF(JLON, JLEV) = FONICE(PT(JLON, JLEV), YDML_PHY_MF%YRPHY0%RDTFAC)
  END DO
  
  !*
  !     ------------------------------------------------------------------
  !     III - CALCULS DE CONDENSATION-EVAPORATION.
  
  !           CONDENSATION-EVAPORATION COMPUTATIONS.
  
  IF (YDML_PHY_MF%YRPHY%LXRCDEV) THEN
    !*
    !     ------------------------------------------------------------------
    !     IIIa - CALCUL DE L'HUMIDITE CRITIQUE (FORMULE ACPLUIE_PROG).
    
    !            CRITICAL HUMIDITY (ACPLUIE_PROG).
    
    ZMESHEXP = (YDML_PHY_MF%YRPHY0%REFLRHC / (YDML_PHY_MF%YRPHY0%TEQH*PGM(JLON)))**YDML_PHY_MF%YRPHY0%RHCEXPDX
    
    DO JLEV=KTDIA,KLEV
      
      !DEC$ IVDEP
      
      !     CALCUL DE L'INFLUENCE DE LA TAILLE DE MAILLE.
      !     COMPUTATION OF THE MESH-SIZE'S INFLUENCE.
      
      ZLEN0 = 1.0_JPRB / (PRMF(JLON, JLEV)*ZLESEFS + (1.0_JPRB - PRMF(JLON, JLEV))*ZLESEFR)
      PHCRICS(JLON, JLEV) = ((YDML_PHY_MF%YRPHY0%HUCRED*YDML_PHY_MF%YRPHY0%RHUC(JLEV) + 1.0_JPRB - YDML_PHY_MF%YRPHY0%HUCRED) &
      & *ZMESHEXP + ZLEN0) / (ZMESHEXP + ZLEN0)
      
      
    END DO
    ! JLEV
    
  END IF
  ! LXRCDEV
  
  IF (YDML_PHY_MF%YRPHY%LSMGCDEV) THEN
    !*
    !     ------------------------------------------------------------------
    !     IIIb - CALCUL DE L'HUMIDITE CRITIQUE (APRES LOPEZ).
    
    !            CRITICAL HUMIDITY (AFTER LOPEZ).
    
    ! CRITICAL RELATIVE HUMIDITY PROFILE
    ! DEPENDENT ON THE LOCAL RESOLUTION.
    ! ----------------------------------
    
    ZFACT0 = (YDML_PHY_MF%YRPHY0%RETAMIN*(YDML_PHY_MF%YRPHY0%RETAMIN - 1.0_JPRB))**2
    ZFACTA = (2.0_JPRB*YDML_PHY_MF%YRPHY0%RETAMIN - 1.0_JPRB)
    ZFACTB = (1.0_JPRB - 3._JPRB*YDML_PHY_MF%YRPHY0%RETAMIN**2)
    ZFACTC = (3._JPRB*YDML_PHY_MF%YRPHY0%RETAMIN - 2.0_JPRB)*YDML_PHY_MF%YRPHY0%RETAMIN
    ZDCRI = (YDML_PHY_MF%YRPHY0%RHCRIT2 - YDML_PHY_MF%YRPHY0%RHCRIT1) / ZFACT0
    
    ZFACT = ZDCRI*(1.0_JPRB - YDML_PHY_MF%YRPHY0%GRHCMOD*EXP(-1.0_JPRB / (YDML_PHY_MF%YRPHY0%TEQH*PGM(JLON))))
    ZA = ZFACTA*ZFACT
    ZB = ZFACTB*ZFACT
    ZC = ZFACTC*ZFACT
    
    DO JLEV=1,KLEV
      !DEC$ IVDEP
      PHCRICS(JLON, JLEV) =  &
      & YDML_PHY_MF%YRPHY0%RHCRIT2 + YDSTA%SVETAF(JLEV)*(ZC + YDSTA%SVETAF(JLEV)*(ZB + YDSTA%SVETAF(JLEV)*ZA))
    END DO
    
  END IF
  ! LSMGCDEV
  
  IF (YDML_PHY_MF%YRPHY%LRKCDEV) THEN
    !     ------------------------------------------------------------------
    !     IIIc - CALCUL DE L'HUMIDITE CRITIQUE (FORMULE RASCH-KRISTJANSSON).
    
    !            CRITICAL HUMIDITY (RASCH-KRISTJANSSON).
    
    
    ZSIGMA = 4.0E-04_JPRB
    ZSIGMASTAB = 1.6E-02_JPRB
    ZEPSILO = YDCST%RMV / YDCST%RMD
    ZSIGMAZ = 1.6E-2_JPRB
    
    DO JLEV=1,KLEV
      ZLSTMP(JLON, JLEV) = 0._JPRB
    END DO
    
    DO JLEV=1,KLEV
      !DEC$ IVDEP
      
      ZMESH = YDML_PHY_MF%YRPHY0%REFLRHC / (YDML_PHY_MF%YRPHY0%TEQH*PGM(JLON))
      
      ZESN = FOEW(PT(JLON, JLEV), 1._JPRB)*PRMF(JLON, JLEV) + FOEW(PT(JLON, JLEV), 0._JPRB)*(1._JPRB - PRMF(JLON, JLEV))
      
      ZDENOM = PAPRSF(JLON, JLEV) - (1._JPRB - ZEPSILO)*ZESN
      
      PQSATS(JLON, JLEV) = ZEPSILO*ZESN / ZDENOM
      
      PRHOUT(JLON, JLEV) = PQ(JLON, JLEV) / PQSATS(JLON, JLEV)
      
      ZRHIN = MAX(PRHOUT(JLON, JLEV), 0.05_JPRB)
      
      ZDRHDZ = ZRHIN*YDCST%RG / (PT(JLON, JLEV)*YDCST%RD)*(ZEPSILO*(YDCST%RLVTT + PRMF(JLON, JLEV)*(YDCST%RLSTT - YDCST%RLVTT)) / &
      &  (YDCST%RCPD*PT(JLON, JLEV)) - 1._JPRB)
      
      ZTEST = MAX(0.0_JPRB, SIGN(1.0_JPRB, PAPHIF(JLON, JLEV) - PAPHI(JLON, KLEV) - YDCST%RG*PBLH(JLON)))
      
      ZZDIST2 = (PAPHI(JLON, JLEV - 1) - PAPHI(JLON, JLEV)) / YDCST%RG
      
      PHCRICS(JLON, JLEV) =  &
      & 1.0_JPRB - 0.5_JPRB*SQRT(2.0_JPRB*ZMESH*ZSIGMA**2 + (1.0_JPRB - ZTEST)*(ZZDIST2*ZDRHDZ)**2 + ZTEST*ZZDIST2*ZSIGMAZ**2)
      
      PHCRICS(JLON, JLEV) = MAX(PHCRICS(JLON, JLEV), 0.5_JPRB)
      
      !   Decrease spinup of phcrics in the beginning of a forecast. No
      !   physical relevance:
      
      ZONEMRHCRIT = MIN(0.075_JPRB, 1._JPRB - PHCRICS(JLON, JLEV))
      ZSTIMEINC = MAX(0._JPRB, MIN(1._JPRB, 1._JPRB - YDRIP%RSTATI / (12._JPRB*3600._JPRB)))
      ZSTIMEINC = ZSTIMEINC**4
      
      PHCRICS(JLON, JLEV) = PHCRICS(JLON, JLEV) - ZSTIMEINC*ZONEMRHCRIT
      
      
    END DO
    
  END IF
  ! LRKCDEV
  
  IF (YDML_PHY_MF%YRPHY%LXRCDEV) THEN
    !*
    !     ------------------------------------------------------------------
    !     IVa - CALCUL DE LA NEBULOSITE (FORMULE XU-RANDAL MODIFIE).
    
    !           CLOUDINESS DIAGNOSED BY MODIFIED XU-RANDAL SCHEME.
    
    ZALFI = 1._JPRB / ZQXRAL
    
    DO JLEV=KTDIA,KLEV
      
      ! - TEMPORAIRE(S) 1D .
      
      ! ZAW       : CONSTANTE CENTRALE DU CALCUL NEBULOSITE/CONDENSAT.
      !           : CENTRAL CONSTANT OF THE CLOUD-COVER/CONDENSATE COMP..
      ! ZDC       : UN MOINS L'HUMIDITE RELATIVE CRITIQUE LOCALE.
      !           : ONE MINUS THE LOCAL CRITICAL RELATIVE HUMIDITY.
      ! ZSITER    : VALEUR ITERATIVE DE LA VARIABLE DE CONTROLE (1/(1-N)).
      !           : ITERATIVE VALUE OF THE CONTROL VARIABLE (1/(1-N)).
      ! ZX0       : VALEUR CIBLE POUR Q/QW-HUC = FONCTION DE (1/(1-N)).
      !           : TARGET VALUE FOR Q/QW-HUC = FUNCTION OF (1/(1-N)).
      !     PREPARATION DES PARAMETRES DE LA BOUCLE DE NEWTON ET PREMIER PAS
      !     ANALYTIQUE.
      !     PREPARATION OF THE NEWTON LOOP'S PARAMETERS AND ANALYTICAL FIRST
      !     STEP.
      
      ZDC = 1._JPRB - PHCRICS(JLON, JLEV)
      ZAW = ZALFI*SQRT(ZDC / PQW(JLON, JLEV))
      ZX0 = MAX(0._JPRB, (PQ(JLON, JLEV) + PQL(JLON, JLEV) + PQI(JLON, JLEV)) / PQW(JLON, JLEV) - PHCRICS(JLON, JLEV))
      ZDNII = 1._JPRB / SQRT(SQRT(PHCRICS(JLON, JLEV)))
      ZXN = PNCV(JLON, JLEV)*(ZDC + ZAW*ZDNII)
      ZDXN = (1._JPRB - PNCV(JLON, JLEV))*ZDC + ZAW*ZDNII*(1._JPRB + PNCV(JLON, JLEV)*(ZDNII - 0.25_JPRB*ZDC / PHCRICS(JLON, JLEV &
      & ) - 1.5_JPRB))
      ZSITER = MAX(1._JPRB, MAX(0._JPRB, (1._JPRB + 0.5_JPRB*(ZX0 - ZXN) / MAX(ZEPS1, ZDXN)))**2)
      
      !     BOUCLE DE NEWTON; LA VARIABLE DE CONTROLE EST L'INVERSE DE UN
      !     MOINS LA NEBULOSITE.
      !     NEWTON LOOP; THE CONTROL VARIABLE IS THE INVERSE OF ONE MINUS THE
      !     CLOUD COVER.
      
      DO JITER=1,IITER
        
        ZNCVNEB = 1._JPRB - PNCV(JLON, JLEV)
        ZNI = 1._JPRB - (1._JPRB / ZSITER)
        ZNIM = 1._JPRB - (ZNCVNEB / ZSITER)
        ZNPI = 1._JPRB - ZDC*(1._JPRB / ZSITER)
        ZDNII = 1._JPRB / (SQRT(SQRT(ZNPI)) - ZNI)
        ZFAC = ZAW / (ZDC*SQRT(ZSITER))
        ZXN = ZDC*ZNIM*(1._JPRB + (ZFAC*ZDNII))
        ZDXN = (ZDC / ZSITER**2)*(ZNCVNEB + (ZFAC*ZDNII)*(ZNCVNEB + ZNIM*(ZDNII*(1._JPRB - 0.25_JPRB*ZDC*SQRT(SQRT(ZNPI)) / ZNPI) &
        &  - 0.5_JPRB*ZSITER)))
        ZSITER = MAX(1._JPRB, MAX(0._JPRB, (ZSITER + 0.5_JPRB*(ZX0 - ZXN) / MAX(ZEPS1, ZDXN)))**2 / ZSITER)
        
      END DO
      
      !     CALCULS FINAUX.
      !     FINAL COMPUTATIONS.
      
      PNEBCOND(JLON, JLEV) = 1._JPRB - (1._JPRB / ZSITER)
      PQSATS(JLON, JLEV) = PQW(JLON, JLEV)
      PRHOUT(JLON, JLEV) = PQ(JLON, JLEV) / PQSATS(JLON, JLEV)
      PNEBS0(JLON, JLEV) = MIN(1.0_JPRB - ZEPS1, MAX(PNEBCOND(JLON, JLEV), ZEPS1))
      
    END DO
    ! JLEV
    
    DO JLEV=KTDIA,KLEV
      ZNEBT = PNEBCOND(JLON, JLEV) + PNCV(JLON, JLEV) - PNEBCOND(JLON, JLEV)*PNCV(JLON, JLEV)
      ZZ = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZNEBT - ZEPS4))
      ZRNEB = ZZ*PNEBCOND(JLON, JLEV) / (ZNEBT + (1.0_JPRB - ZZ))
      
      ZZ = MAX(0.0_JPRB, SIGN(1.0_JPRB, PNEBCOND(JLON, JLEV) - ZEPS4))
      ZIRNEB = ZZ*ZNEBT / (PNEBCOND(JLON, JLEV) + (1.0_JPRB - ZZ))
      
      ZQVCS = MAX(0.0_JPRB, PQ(JLON, JLEV) - PQW(JLON, JLEV)*(PNEBCOND(JLON, JLEV) + PNCV(JLON, JLEV) - PNEBCOND(JLON, JLEV) &
      & *PNCV(JLON, JLEV)))
      ZQ = (PQ(JLON, JLEV) - ZQVCS)*ZRNEB + ZQVCS / MAX(ZEPS4, 1.0_JPRB - PNCV(JLON, JLEV))
      ZADJL = ZWEIGHT*((1.0_JPRB - PRMF(JLON, JLEV))*(PQL(JLON, JLEV) + PQI(JLON, JLEV)) - PQL(JLON, JLEV))
      ZQL = (PQL(JLON, JLEV) + ZADJL)*ZRNEB
      ZADJI = ZWEIGHT*(PRMF(JLON, JLEV)*(PQL(JLON, JLEV) + PQI(JLON, JLEV)) - PQI(JLON, JLEV))
      ZQI = (PQI(JLON, JLEV) + ZADJI)*ZRNEB
      ZALF = ZQI / MAX(ZEPS4, ZQL + ZQI)
      ZNEBM1 = 1.0_JPRB - PNEBCOND(JLON, JLEV)
      ZQVN = PQW(JLON, JLEV)*(1.0_JPRB + (PHCRICS(JLON, JLEV) - 1.0_JPRB)*ZNEBM1)
      ZDQVN = ZQ - ZQVN
      ZTEST = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZQVN - (ZQ + ZQL + ZQI)))
      ZTESTB = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZDQVN))
      ZRMLOC = (1.0_JPRB - ZTESTB)*ZALF + ZTESTB*PRMF(JLON, JLEV)
      ZCONL = (1.0_JPRB - ZRMLOC)*(1.0_JPRB - ZTEST)*ZDQVN - ZTEST*ZQL
      ZCONI = ZRMLOC*(1.0_JPRB - ZTEST)*ZDQVN - ZTEST*ZQI
      
      ZCONL = (1.0_JPRB - PNCV(JLON, JLEV))*ZCONL
      ZCONI = (1.0_JPRB - PNCV(JLON, JLEV))*ZCONI
      
      PQCS(JLON, JLEV) = (MAX(0.0_JPRB, ZQL + ZCONL) + MAX(0.0_JPRB, ZQI + ZCONI))*ZIRNEB
      PQCS0(JLON, JLEV) = PQCS(JLON, JLEV)
      
    END DO
    
  END IF
  ! LXRCDEV
  
  IF (YDML_PHY_MF%YRPHY%LSMGCDEV) THEN
    !*
    !     ------------------------------------------------------------------
    !     IVb - CALCUL DE LA NEBULOSITE (FORMULE SMITH-GERARD).
    
    !           CLOUDINESS DIAGNOSED BY SMITH-GERARD SCHEME.
    
    ! Initialisations for the Smith scheme.
    
    ! Compute the saturation with respect to ice and liquid.
    
    DO JLEV=KTDIA,KLEV - 1
      ZQTOT(JLON, JLEV) = PQI(JLON, JLEV) + PQL(JLON, JLEV) + PQ(JLON, JLEV)
      ZS(JLON, JLEV) = PCP(JLON, JLEV)*PT(JLON, JLEV) + PAPHIF(JLON, JLEV)
      
      ZESATL = FOEW(PT(JLON, JLEV), 0.0_JPRB)
      ZESAT = ZESATL*(1.0_JPRB - PRMF(JLON, JLEV)) + FOEW(PT(JLON, JLEV), 1.0_JPRB)*PRMF(JLON, JLEV)
      
      ZESP = ZESATL / PAPRSF(JLON, JLEV)
      ZRHL(JLON, JLEV) = ZQTOT(JLON, JLEV) / FOQS(ZESP)
      
      ZESP = ZESAT / PAPRSF(JLON, JLEV)
      PQSATS(JLON, JLEV) = FOQS(ZESP)
      
    END DO
    
    ! Lowest level: use optionally the surface T in the estimation.
    
    !DEC$ IVDEP
    
    IF (YDML_PHY_MF%YRPHY%NSMTBOT == 0) THEN
      ZTKLEV = (PTS(JLON) + PT(JLON, KLEV))*0.5_JPRB
    ELSE
      ZTKLEV = PT(JLON, KLEV)
    END IF
    
    ZRMF = PRMF(JLON, KLEV)
    PRMF(JLON, KLEV) = FONICE(ZTKLEV, YDML_PHY_MF%YRPHY0%RDTFAC)
    
    ZQTOT(JLON, KLEV) = PQI(JLON, KLEV) + PQL(JLON, KLEV) + PQ(JLON, KLEV)
    ZS(JLON, KLEV) = PCP(JLON, KLEV)*ZTKLEV + PAPHIF(JLON, KLEV)
    
    ZESATL = FOEW(ZTKLEV, 0.0_JPRB)
    ZESAT = ZESATL*(1.0_JPRB - PRMF(JLON, KLEV)) + FOEW(ZTKLEV, 1.0_JPRB)*PRMF(JLON, KLEV)
    
    ZESP = ZESATL / PAPRSF(JLON, KLEV)
    ZRHL(JLON, KLEV) = ZQTOT(JLON, KLEV) / FOQS(ZESP)
    
    ZESP = ZESAT / PAPRSF(JLON, KLEV)
    PQSATS(JLON, KLEV) = FOQS(ZESP)
    
    ! Treble point smoothing: on PQSATS only.
    
    IF (YDML_PHY_MF%YRPHY%LSMTPS) THEN
      
      ! Computation done only for liquid containing levels.
      
      DO JLEV=KTDIA,KLEV
        ZS1 = (ZQTOT(JLON, JLEV) / PQSATS(JLON, JLEV))*PDELP(JLON, JLEV)
        ZS2 = PDELP(JLON, JLEV)
        ILIQ = MAX(0, SIGN(1, 1 - NINT(1000._JPRB*PRMF(JLON, JLEV))))
        
        ! Use NSMTPB layers below when dry adiabatic unstable.
        
        ILEV = MIN(JLEV + YDML_PHY_MF%YRPHY0%NSMTPB, KLEV)
        
        DO JLEV2=JLEV + 1,ILEV
          ZDS = ZS(JLON, JLEV2) - ZS(JLON, JLEV)
          ZDELTAP = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZDS))*PDELP(JLON, JLEV2)
          ZS1 = ZS1 + (ZDELTAP*ILIQ)*ZRHL(JLON, JLEV2)
          ZS2 = ZS2 + (ZDELTAP*ILIQ)
        END DO
        
        ! Use NSMTPA layers above when dry adiabatic unstable (reversal of ZDS sign).
        
        ILEV = MAX(JLEV - YDML_PHY_MF%YRPHY0%NSMTPA, KTDIA)
        
        DO JLEV2=ILEV,JLEV - 1
          ZDS = ZS(JLON, JLEV) - ZS(JLON, JLEV2)
          ZDELTAP = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZDS))*PDELP(JLON, JLEV2)
          ZS1 = ZS1 + (ZDELTAP*ILIQ)*ZRHL(JLON, JLEV2)
          ZS2 = ZS2 + (ZDELTAP*ILIQ)
        END DO
        
        ! Keep the original if ZS1 small but ALSO if ZQTOT is too small (desert)!
        
        ZZ = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZS1 - ZEPS1))*MAX(0.0_JPRB, SIGN(1.0_JPRB, ZQTOT(JLON, JLEV) - ZEPS1))
        PQSATS(JLON, JLEV) = ZZ*ZS2 / (ZS1 + (1.0_JPRB - ZZ))*ZQTOT(JLON, JLEV) + (1.0_JPRB - ZZ)*PQSATS(JLON, JLEV)
        
      END DO
      
    END IF
    ! LSMTPS
    
    ! Reinstalling the true ice proportion at the bottom (for modular use).
    
    PRMF(JLON, KLEV) = ZRMF
    
    ! Stratiform cloud cover, based on PQSATS.
    
    DO JLEV=KTDIA,KLEV
      !DEC$ IVDEP
      
      ZQCST = PQI(JLON, JLEV) + PQL(JLON, JLEV)
      ZRHC = PHCRICS(JLON, JLEV)
      ZRATQ = ZQTOT(JLON, JLEV) / PQSATS(JLON, JLEV)
      ZZRHC = MIN(ZRHC, 1.0_JPRB - ZEPS1)
      ZRATQ = MAX(ZZRHC, MIN(2.0_JPRB - ZZRHC, ZRATQ))
      ZRAT2 = (ZRATQ - 1.0_JPRB) / (1.0_JPRB - ZZRHC)
      ZLV = FOLH(PT(JLON, JLEV), 0._JPRB)
      ZLS = FOLH(PT(JLON, JLEV), 1._JPRB)
      ZLEFF = ZLV*(1.0_JPRB - PRMF(JLON, JLEV)) + ZLS*PRMF(JLON, JLEV)
      ZTLIQ = PT(JLON, JLEV) - (ZLV*PQL(JLON, JLEV) + ZLS*PQI(JLON, JLEV)) / PCP(JLON, JLEV)
      ZSIGS = (1.0_JPRB - ZRHC) / ZSQRT6*PQSATS(JLON, JLEV) / (1.0_JPRB + ZLEFF*ZLEFF*PQSATS(JLON, JLEV)*ZRDSRV / PR(JLON, JLEV)  &
      & / ZTLIQ / ZTLIQ / PCP(JLON, JLEV))
      
      ZTEST = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZRATQ - 1.0_JPRB))
      PNEBCOND(JLON, JLEV) =  &
      & ZTEST*(1.0_JPRB - 0.5_JPRB*(1.0_JPRB - ZRAT2)**2) + (1.0_JPRB - ZTEST)*(0.5_JPRB*(1.0_JPRB + ZRAT2)**2)
      
      PQCS(JLON, JLEV) = (ZTEST*(6.0_JPRB*ZRAT2 + (1.0_JPRB - ZRAT2)**3) + (1.0_JPRB - ZTEST)*(1.0_JPRB + ZRAT2)**3) / ZSQRT6
      
      PQCS(JLON, JLEV) = PQCS(JLON, JLEV)*ZSIGS
      
      ZNEBT = PNEBCOND(JLON, JLEV) + PNCV(JLON, JLEV) - PNEBCOND(JLON, JLEV)*PNCV(JLON, JLEV)
      ZZ = MAX(0.0_JPRB, SIGN(1.0_JPRB, ZNEBT - ZEPS1))
      ZFRACON = ZZ*PNCV(JLON, JLEV) / (ZNEBT + (1.0_JPRB - ZZ))
      
      PQCS(JLON, JLEV) = (PQCS(JLON, JLEV)**2 + (ZQCST*ZFRACON)**2) / MAX(ZEPS1, (PQCS(JLON, JLEV) + ZQCST*ZFRACON))
      PQCS0(JLON, JLEV) = PQCS(JLON, JLEV)
      
      ! Where there is condensate impose a minimum of cloudiness.
      
      PNEBCOND(JLON, JLEV) = MAX(PNEBCOND(JLON, JLEV), ZEPS2*MAX(0.0_JPRB, SIGN(1.0_JPRB, ZQCST - ZEPS3)))
      
      ! Cloud cover reduction to avoid excessive values
      ! when overlapping is performed in radiation.
      
      PNEBCOND(JLON, JLEV) =  &
      & PNEBCOND(JLON, JLEV) / (1.0_JPRB + (PAPHI(JLON, JLEV - 1) - PAPHI(JLON, JLEV)) / YDML_PHY_MF%YRPHY0%RDPHIC)
      
      PNEBS0(JLON, JLEV) = MIN(1.0_JPRB - ZEPS1, MAX(PNEBCOND(JLON, JLEV), ZEPS1))
      PRHOUT(JLON, JLEV) = PQ(JLON, JLEV) / PQSATS(JLON, JLEV)
    END DO
    
    ! Smoothing of the cloudiness where there are too big jumps
    ! diluting an excess condensate and preventing its autoconversion.
    
    IF (YDML_PHY_MF%YRPHY%NSMDNEB == 1) THEN
      
      ! Smoothing the profile.
      
      ZNEB1 = 0._JPRB
      
      DO JLEV=KTDIA,KLEV - 1
        ZNEB0 = 0.5_JPRB*(PNEBCOND(JLON, JLEV) + PNEBCOND(JLON, JLEV + 1))
        PNEBCOND(JLON, JLEV) = (ZNEB1 + ZNEB0 + PNEBCOND(JLON, JLEV)*2.0_JPRB)*0.25_JPRB
        ZNEB1 = ZNEB0
        PNEBS0(JLON, JLEV) = MIN(1.0_JPRB - ZEPS1, MAX(PNEBCOND(JLON, JLEV), ZEPS1))
      END DO
      
      ZNEB0 = PNEBCOND(JLON, KLEV)
      PNEBCOND(JLON, KLEV) = (ZNEB1 + ZNEB0 + PNEBCOND(JLON, KLEV)*2.0_JPRB)*0.25_JPRB
      PNEBS0(JLON, KLEV) = MIN(1.0_JPRB - ZEPS1, MAX(PNEBCOND(JLON, KLEV), ZEPS1))
      
    ELSE IF (YDML_PHY_MF%YRPHY%NSMDNEB == 2) THEN
      
      ! Limitation of the gradient.
      
      ZNEB1 = 0._JPRB
      
      DO JLEV=KTDIA,KLEV
        !DEC$ IVDEP
        ZNEB0 = PNEBCOND(JLON, JLEV)
        ZDNEB = ZNEB0 - ZNEB1
        ZNEB1 = ZNEB1 + SIGN(MIN(ABS(ZDNEB), YDML_PHY_MF%YRPHY0%RSMDNEBX), ZDNEB)
        PNEBCOND(JLON, JLEV) = ZNEB1
        PNEBS0(JLON, JLEV) = MIN(1.0_JPRB - ZEPS1, MAX(PNEBCOND(JLON, JLEV), ZEPS1))
      END DO
      
    END IF
    ! NSMDNEB
    
  END IF
  ! LSMGCDEV
  IF (YDML_PHY_MF%YRPHY%LRKCDEV) THEN
    !  ------------------------------------------------------------------
    !           IVc - CALCUL DE LA NEBULOSITE (FORMULE RASCH-KRISTJANSSON).
    
    !           CLOUDINESS DIAGNOSED BY MODIFIED RASCH-KRISTJANSSON CAM3.
    ! ---------------------------------------------------------------------
    
    !    1. Define some constants for the scheme depending on height
    
    !   ---------------------------------------------------------------------------
    !    2. In HIRLAM this is where shallow and convective cloudiness is computed,
    !       Need to know what to do with 3MT, and Meteo-France physics for Cu and
    !       shallow convection. Now I use historic PNCV for convective cloudiness.
    !   ---------------------------------------------------------------------------
    
    DO JLEV=KTDIA,KLEV
      
      ZSHALTEMP(JLON, JLEV) = 0._JPRB
      ZHCUTMP(JLON, JLEV) = 0._JPRB
      
    END DO
    !   ---------------------------------------------------------------------------
    !    3. Compute stratiform cloudiness and merge it which shallow convective
    !       clouds from (2.) into "total" clouds, excluding deep convective clouds.
    !   ---------------------------------------------------------------------------
    
    DO JLEV=2,KLEV
      
      !      Following lines currently not needed, but may be later on if shallow
      !      convection enters from TOUCANS:
      !      ZCLD = MIN(0.8_JPRB,ZHCUTMP(JLON,JLEV)+ZSHALTEMP(JLON,JLEV))
      !      ZRHD = (PRH(JLON,JLEV)-ZCLD)/(1._JPRB-ZCLD)
      
      ZRHD = MIN(1._JPRB, PRHOUT(JLON, JLEV))
      ZRHDIF = (1._JPRB - ZRHD) / (1._JPRB - PHCRICS(JLON, JLEV))
      ZRHDIF = 1._JPRB - SQRT(MAX(0._JPRB, ZRHDIF))
      ZLSTMP(JLON, JLEV) = MIN(0.999_JPRB, MAX(ZRHDIF, 0.0_JPRB))
      
      !      Following lines currently not needed, but may be later on if shallow
      !      convection enters from TOUCANS:
      !      ZLSTMP(JLON,JLEV)=ZLSTMP(JLON,JLEV)*(1._JPRB-ZSHALTEMP(JLON,JLEV))&
      !     & + ZSHALTEMP(JLON,JLEV)
      
    END DO
    
    
    !     ----------------------------------------------------------------------------
    !     4. Merge deep convective and layered cloud fraction for total cloud.
    !     ----------------------------------------------------------------------------
    
    DO JLEV=KTDIA,KLEV
      
      !      Following lines currently not needed, but may be later on:
      !      IF((ZHCUTMP(JLON,JLEV)+ZLSTMP(JLON,JLEV)) > 1.0_JPRB)THEN
      !        ZHCUTMP(JLON,JLEV)=ZHCUTMP(JLON,JLEV)/(ZHCUTMP(JLON,JLEV)&
      !       &+ZLSTMP(JLON,JLEV))
      !        ZLSTMP(JLON,JLEV)=ZLSTMP(JLON,JLEV)/(ZHCUTMP(JLON,JLEV)&
      !       &+ZLSTMP(JLON,JLEV))
      !      ENDIF
      
      ! Other output variables:
      
      PNEBCOND(JLON, JLEV) = MAX(0.0_JPRB, MIN(0.99_JPRB, ZLSTMP(JLON, JLEV)))
      
    END DO
    
  END IF
  ! LRKCDEV
  
  IF (YDML_PHY_MF%YRPHY%LSMITH_CDEV) THEN
    !  ------------------------------------------------------------------
    !           IVd - CALCUL DE LA NEBULOSITE (FORMULE SMITH).
    
    !           CLOUDINESS DIAGNOSED BY SMITH SCHEME.
    ! ---------------------------------------------------------------------
    
    !    1. Computation following Smith QJRMS 1990 paper: routine ACNEBSM
    
    CALL ACNEBSM_OPENACC(YDCST, YDML_PHY_MF%YRPHY0, KIDIA, KFDIA, KLON, KTDIA, KLEV, PT, PQ, PQL, PQI, PAPHI, PAPRSF, PCP, PR,  &
    & PGM, YDSTA, PQCS, PNEBCOND, PHCRICS, PRMF, PQSATS, YDSTACK=YLSTACK)
    
    IF (YDML_PHY_MF%YRPHY%L3MT) THEN
      
      DO JLEV=KTDIA,KLEV
        !DEC$ IVDEP
        
        ! Where there is condensate impose a minimum of cloudiness.
        
        PNEBCOND(JLON, JLEV) = MAX(PNEBCOND(JLON, JLEV), ZEPS2*MAX(0.0_JPRB, SIGN(1.0_JPRB, PQCS(JLON, JLEV) - ZEPS3)))
        
        ! Cloud cover reduction to avoid excessive values
        ! when overlapping is performed in radiation.
        
        PNEBCOND(JLON, JLEV) =  &
        & PNEBCOND(JLON, JLEV) / (1.0_JPRB + (PAPHI(JLON, JLEV - 1) - PAPHI(JLON, JLEV)) / YDML_PHY_MF%YRPHY0%RDPHIC)
        
        PNEBS0(JLON, JLEV) = MIN(1.0_JPRB - ZEPS1, MAX(PNEBCOND(JLON, JLEV), ZEPS1))
        
        PNEBCOND(JLON, JLEV) = MIN(1.0_JPRB, MAX(PNEBCOND(JLON, JLEV), 0.0_JPRB))
        
        ZNEBT = PNEBCOND(JLON, JLEV) + PNCV(JLON, JLEV) - PNEBCOND(JLON, JLEV)*PNCV(JLON, JLEV)
        ZQCST = PQI(JLON, JLEV) + PQL(JLON, JLEV)
        
        PQCS(JLON, JLEV) =  &
        & (ZQCST*PNCV(JLON, JLEV) + PQCS(JLON, JLEV)*(1._JPRB - PNCV(JLON, JLEV))*PNEBCOND(JLON, JLEV)) / MAX(ZEPS4, ZNEBT)
        
        PQCS0(JLON, JLEV) = PQCS(JLON, JLEV)
        
        PRHOUT(JLON, JLEV) = PQ(JLON, JLEV) / PQSATS(JLON, JLEV)
        
      END DO
      
    ELSE
      
      DO JLEV=KTDIA,KLEV
        PNEBCOND(JLON, JLEV) = MIN(1.0_JPRB - ZEPS5, MAX(PNEBCOND(JLON, JLEV), ZEPS5))
        PNEBS0(JLON, JLEV) = PNEBCOND(JLON, JLEV)
        PQCS0(JLON, JLEV) = PQCS(JLON, JLEV)
        PRHOUT(JLON, JLEV) = PQ(JLON, JLEV) / PQSATS(JLON, JLEV)
      END DO
      
    END IF
    ! L3MT
    
  END IF
  ! LSMITH_CDEV
  !-----------------------------------------------------------------------
END SUBROUTINE ACNEBCOND_OPENACC
