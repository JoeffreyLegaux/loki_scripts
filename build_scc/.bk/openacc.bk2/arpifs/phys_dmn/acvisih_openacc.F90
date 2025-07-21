SUBROUTINE ACVISIH_OPENACC (YDCST, YDVISI, KIDIA, KFDIA, KLON, KTDIA, KLEV, PAPHI, PAPHIF, PAPRSF, PT, PR, PQL, PQI, PQR, PQS,  &
& PQG, PVISICLD, PVISIHYD, PMXCLWC, YDSTACK)
  !-----------------------------------------------------------------------
  ! - INPUT  2D .
  ! - OUTPUT 1D .
  
  !**** *ACVISIH * - HORIZONTAL VISIBILITY DUE TO FOG AND PRECIPITATIONS.
  
  !     Sujet.
  !     ------
  
  !     - ROUTINE DE CALCUL ACTIF.
  !          VISIBILITES LIEES AUX NUAGES ET AUX PRECIPITATIONS.
  !          CONTENU DE L'EAU NUAGEUSE POUR LE CALCUL DE MAXIMUM SUR UNE PERIODE.
  
  !          VISIBILITIES DUE TO CLOUD WATER (FOG) AND PRECIPITATIONS.
  !          LIQUID CLOUD WATER CONTENT TO COMPUTE MAXIMUM OVER A PERIOD.
  
  !**   Interface.
  !     ----------
  !        *CALL* *ACVISIH*
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS D'ENTREE.
  ! -   INPUT ARGUMENTS.
  !     -------------------
  
  ! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.
  ! - DIMENSIONS.
  
  ! KIDIA, KFDIA : BORNES BOUCLES HORIZONTALES   (IST,IEND DANS CPG).
  ! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
  ! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
  ! KLON : HORIZONTAL DIMENSION                  (NPROMA IN *CPG*).
  ! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
  ! KTDIA : START OF THE VERTICAL LOOP IN THE PHYSICS
  ! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
  ! KLEV : END OF VERTICAL LOOP AND VERTICAL DIMENSION(NFLEVG IN *CPG*).
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 2D (0:KLEV) .
  
  ! PAPHI      : GEOPOTENTIEL AUX DEMI-NIVEAUX.
  ! PAPHI      : GEOPOTENTIAL ON HALF-LEVELS.
  
  ! - 2D (1:KLEV) .
  
  ! PAPHIF     : GEOPOTENTIEL AUX NIVEAUX.
  ! PAPHIF     : GEOPOTENTIAL ON FULL LEVELS.
  ! PAPRSF     : PRESSION AUX NIVEAUX.
  ! PAPRSF     : PRESSURE ON FULL LEVELS.
  ! PT         : TEMPERATURE.
  ! PR         : CONSTANTE DE GAZ - L'AIR.
  ! PR         : GAS CONSTANT FOR AIR.
  ! PQL        : HUMIDITE SPECIFIQUE DE L'EAU LIQUIDE NUAGEUSE.
  ! PQL        : SPECIFIC HUMIDITY OF CLOUD WATER.
  ! PQI        : HUMIDITE SPECIFIQUE DE LA GLACE NUAGEUSE.
  ! PQI        : SPECIFIC HUMIDITY OF ICE.
  ! PQR        : HUMIDITE SPECIFIQUE DE LA PLUIE TOTALE.
  ! PQR        : SPECIFIC HUMIDITY OF RAIN.
  ! PQS        : HUMIDITE SPECIFIQUE DE LA NEIGE TOTALE.
  ! PQS        : SPECIFIC HUMIDITY OF SNOW.
  ! PQG        : HUMIDITE SPECIFIQUE DU GRAUPEL TOTAL.
  ! PQG        : SPECIFIC HUMIDITY OF GRAUPEL.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS DE SORTIE.
  !     --------------------
  
  ! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
  !   CATEGORIE).
  
  ! - 1D (DIAGNOSTIQUE) .
  
  ! PVISICLD   : VISIBILITE LIEE AUX NUAGES.
  ! PVISICLD   : VISIBILITE DUE TO CLOUDS (FOG)
  ! PVISIHYD   : VISIBILITE LIEE AUX PRECIPITATIONS.
  ! PVISIHYD   : VISIBILITE DUE TO PRECIPITATIONS
  ! PMXCLWC    : MAX DE L'HUMIDITE SPECIFIQUE DE L'EAU LIQUIDE NUAGEUSE.
  ! PMXCLWC    : MAX OF SPECIFIC HUMIDITY OF CLOUD WATER.
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS IMPLICITES.
  !     ---------------------
  
  ! COMMON/YOMCST /
  ! COMMON/YOMPHY2 /
  
  !-----------------------------------------------------------------------
  
  !     Externes.
  !     ---------
  
  !     Methode.
  !     --------
  
  !     Auteurs.
  !     -------
  !        2019-01, R. Brozkova according to the code of I. Etchevers.
  !
  !     Modifications.
  !     --------------
  !        2019-10, I. Etchevers : optimization, cleaning and new formula
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !
  !-----------------------------------------------------------------------
  
!$acc routine( ACVISIH_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOMDVISI, ONLY: TDVISI
  
  
  !-----------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TDVISI), INTENT(IN) :: YDVISI
  
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  REAL(KIND=JPRB), INTENT(IN) :: PAPHI(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPRSF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PR(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQI(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQR(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQG(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PVISICLD(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PVISIHYD(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PMXCLWC(KLON)
  
  !-----------------------------------------------------------------------
  
  ! ZCONTENT : (cloud liquid water, cloud ice, rain, snow or graupel) Content in g/m3
  REAL(KIND=JPRB) :: ZLCONTENT
  REAL(KIND=JPRB) :: ZICONTENT
  REAL(KIND=JPRB) :: ZRCONTENT
  REAL(KIND=JPRB) :: ZSCONTENT
  REAL(KIND=JPRB) :: ZGCONTENT
  temp (REAL (KIND=JPRB), ZRHOAIR, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZC1
  REAL(KIND=JPRB) :: ZC2
  
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: ILEV
  INTEGER(KIND=JPIM) :: ILEVH
  INTEGER(KIND=JPIM) :: ILEVB
  
  REAL(KIND=JPRB) :: ZBETARAYL
  REAL(KIND=JPRB) :: ZCONTRAST
  REAL(KIND=JPRB) :: ZZVISI
  REAL(KIND=JPRB) :: ZEPS
  REAL(KIND=JPRB) :: ZEPSQ
  REAL(KIND=JPRB) :: ZSCALE
  REAL(KIND=JPRB) :: ZDELTA
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZRHOAIR) == 8) THEN
    alloc8 (ZRHOAIR)
  ELSE
    IF (KIND (ZRHOAIR) == 4) THEN
      alloc4 (ZRHOAIR)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KIDIA
  
  !-----------------------------------------------------------------------
  !ASSOCIATE(HVISI=>YRPHY2%HVISI,COEFFEXTQ=>YRPHY2%COEFFEXTQ,COEFFPWRQ=>YRPHY2%COEFFPWRQ)
  
  !-----------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  !     PREPARATION.
  
  !         INITIAL SETTINGS.
  
  ZCONTRAST = -LOG(0.05_JPRB)    ! log du seuil du contraste de visibilite
  ZBETARAYL = 0.013_JPRB    ! Rayleigh diffusion moleculaire ATTENTION en km^(-1)
  ZEPS = EPSILON(1._JPRB)*100._JPRB    ! protection contre la division par zero
  ZEPSQ = 1.E-6_JPRB    ! contenu minimum EN g/m3 !
  ZSCALE = 1.E+3_JPRB    ! facteur pour passer en g/kg !
  
  
  !------------------------------------------------------------------------------
  !         CALCUL DU NIVEAU AU DESSUS DE LA HAUTEUR HVISI ET DE LA DENSITE D'AIR.
  !         COMPUTATION OF LEVEL ABOVE HVISI HEIGHT AND AIR DENSITY.
  !------------------------------------------------------------------------------
  
  
  ILEV = 1
  
  ! Search level ILEV above HVISI
  
  DO JLEV=KLEV,1,-1
    ZZVISI = PAPHIF(JLON, JLEV) - PAPHI(JLON, KLEV) - YDVISI%HVISI*YDCST%RG
    ILEVH = MAX(0._JPRB, SIGN(1._JPRB, ZZVISI))*JLEV
    ILEV = MAX(ILEVH, ILEV)
    
    ZRHOAIR(JLON, JLEV) = PAPRSF(JLON, JLEV) / (PR(JLON, JLEV)*PT(JLON, JLEV))
  END DO
  
  !-----------------------------------------------------------
  !         CALCUL DU CONTENU DES HYDROMETEORES.
  !         COMPUTATION OF HYDROMETEORES CONTENT.
  !-----------------------------------------------------------
  
  ILEVH = ILEV
  ILEVB = MIN(ILEV + 1, KLEV)
  ZDELTA =  &
  & (PAPHIF(JLON, ILEVH) - YDVISI%HVISI*YDCST%RG - PAPHI(JLON, KLEV)) / MAX(ZEPS, (PAPHIF(JLON, ILEVH) - PAPHIF(JLON, ILEVB)))
  
  IF (ILEVH == ILEVB) ZDELTA = 1._JPRB
  
  ZLCONTENT =  &
  & MAX(ZEPSQ, ZSCALE*(PQL(JLON, ILEVB)*ZRHOAIR(JLON, ILEVB)*ZDELTA + PQL(JLON, ILEVH)*ZRHOAIR(JLON, ILEVH)*(1 - ZDELTA)))
  ZICONTENT =  &
  & MAX(ZEPSQ, ZSCALE*(PQI(JLON, ILEVB)*ZRHOAIR(JLON, ILEVB)*ZDELTA + PQI(JLON, ILEVH)*ZRHOAIR(JLON, ILEVH)*(1 - ZDELTA)))
  ZRCONTENT =  &
  & MAX(ZEPSQ, ZSCALE*(PQR(JLON, ILEVB)*ZRHOAIR(JLON, ILEVB)*ZDELTA + PQR(JLON, ILEVH)*ZRHOAIR(JLON, ILEVH)*(1 - ZDELTA)))
  ZSCONTENT =  &
  & MAX(ZEPSQ, ZSCALE*(PQS(JLON, ILEVB)*ZRHOAIR(JLON, ILEVB)*ZDELTA + PQS(JLON, ILEVH)*ZRHOAIR(JLON, ILEVH)*(1 - ZDELTA)))
  ZGCONTENT =  &
  & MAX(ZEPSQ, ZSCALE*(PQG(JLON, ILEVB)*ZRHOAIR(JLON, ILEVB)*ZDELTA + PQG(JLON, ILEVH)*ZRHOAIR(JLON, ILEVH)*(1 - ZDELTA)))
  PMXCLWC(JLON) = PQL(JLON, ILEVB)*ZRHOAIR(JLON, ILEVB)*ZDELTA + PQL(JLON, ILEVH)*ZRHOAIR(JLON, ILEVH)*(1 - ZDELTA)
  
  !-----------------------------------------------------------
  !         CALCUL DES VISIBILITES.
  !-----------------------------------------------------------
  
  
  ZC1 = ZSCALE*EXP(YDVISI%COEF_CM2*(LOG(ZLCONTENT)**2))    !mise en km-1
  ZC2 = EXP(YDVISI%COEF_CM3*(LOG(ZLCONTENT)**3))
  
  PVISICLD(JLON) =  &
  & ZCONTRAST / (YDVISI%COEF_CM1*ZC1*ZC2*(ZLCONTENT**YDVISI%COEF_CM4) + YDVISI%COEF_IM1*ZICONTENT**YDVISI%COEF_IM2)
  ! conversion from km to m !
  PVISICLD(JLON) = ZSCALE*MIN(PVISICLD(JLON), 20.0_JPRB)
  
  PVISIHYD(JLON) = ZCONTRAST / (ZBETARAYL + YDVISI%COEF_RM1*ZRCONTENT**YDVISI%COEF_RM2 +  &
  & YDVISI%COEF_SM1*ZSCONTENT**YDVISI%COEF_SM2 + YDVISI%COEF_GM1*ZGCONTENT**YDVISI%COEF_GM2)
  ! conversion from km to m !
  PVISIHYD(JLON) = ZSCALE*MIN(PVISIHYD(JLON), 20.0_JPRB)
  
  !-----------------------------------------------------------------------
END SUBROUTINE ACVISIH_OPENACC
