SUBROUTINE ACMICRO_OPENACC (YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, PNEBST, PT, PQL, PQI, PTS, PNEIJ, PLSM, PAUTOL,  &
& PAUTOI, YDSTACK)
  !-----------------------------------------------------------------------
  ! - INPUT -
  ! - OUTPUT -
  
  !**** *ACMICRO * - CALCULS D'AUTOCONVERSION.
  !                  AUTOCONVERSION COMPUTATIONS.
  
  !     Sujet.
  !     ------
  
  !**   Interface.
  !     ----------
  !        *CALL* *ACMICRO*
  
  !-----------------------------------------------------------------------
  
  ! -   ARGUMENTS D'ENTREE./INPUT ARGUMENTS.
  !     ------------------------------------
  
  ! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.
  
  ! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
  !            : START OF HORIZONTAL LOOP
  ! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
  !            : END OF HORIZONTAL LOOP
  ! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
  !            : HORIZONTAL DIMENSION
  ! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES.
  !            : START OF THE VERTICAL LOOP IN THE PHYSICS.
  ! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
  !            : END OF VERTICAL LOOP AND VERTICAL DIMENSION
  
  ! - 2D (1:KLEV)
  
  ! PNEBST     : NEBULOSITE PARTIELLE STRATIFORME.
  !            : STRATIFORM FRACTIONAL CLOUDINESS.
  ! PT         : TEMPERATURE.
  !            : TEMPERATURE.
  ! PQL        : CONTENU EN CONDENSAT NUAGEUX LIQUIDE.
  !            : CLOUD WATER LIQUID.
  ! PQI        : CONTENU EN CONDENSAT NUAGEUX SOLIDE.
  !            : CLOUD WATER SOLID.
  
  ! - 1D (1:KLON) .
  
  ! PTS        : TEMPERATURE DE SURFACE.
  !            : SURFACE TEMPERATURE.
  ! PNEIJ      : PROPORTION DE LA MAILLE RECOUVERTE DE NEIGE.
  !            : SNOW FRACTION.
  ! PLSM       : INDICE TERRE/MER.
  !            : LAND/SEA MASK.
  
  ! -   ARGUMENTS EN SORTIE.
  !     ---------------------------
  
  ! PAUTOL     : GENERATION DE PRECIPITATIONS A PARTIR DE L'EAU NUAGEUSE LIQ.
  !            : GENERATION OF PRECIPITATION FROM LIQUID CLOUD WATER.
  ! PAUTOI     : GENERATION DE PRECIPITATIONS A PARTIR DE L'EAU NUAGEUSE SOL.
  !            : GENERATION OF PRECIPITATION FROM SOLID CLOUD WATER.
  
  !-----------------------------------------------------------------------
  
  !     Auteur.
  !     -------
  !        99-10, Philippe Lopez
  !               CNRM/GMME/RECYF
  
  !     Modifications.
  !     --------------
  !      2002-01, P. Marquet Comments :
  !               here PAUTO is >0 when conversion occur, and the
  !               contribution is "-PAUTO" in the "qc" equation.
  !      2004-06, P. Marquet : introduce RAUTEFR for liquid->rain, with
  !                            the old RAUTEFS valid for ice->snow
  !        M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      2004-10, F. Bouyssel : Cleaning
  !      2005-05, F. Bouyssel : More tunings in autoconvertion thresholds
  !      2005-07, F. Bouyssel : Modification of ice autoconvertion thresholds
  !                             over snow and sea ice (RQICRSN)
  !      2006-01, F. Bouyssel : Remove ZSURF
  !      2006-04, F. Bouyssel : Replace PQC by PQL,PQI
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  
  !-----------------------------------------------------------------------
  
!$acc routine( ACMICRO_OPENACC ) seq
  
  USE MODEL_PHYSICS_MF_MOD, ONLY: MODEL_PHYSICS_MF_TYPE
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(MODEL_PHYSICS_MF_TYPE), INTENT(IN) :: YDML_PHY_MF
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  REAL(KIND=JPRB), INTENT(IN) :: PNEBST(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQI(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTS(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PNEIJ(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PLSM(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PAUTOL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PAUTOI(KLON, KLEV)
  
  REAL(KIND=JPRB) :: ZQL
  REAL(KIND=JPRB) :: ZQI
  REAL(KIND=JPRB) :: ZCAUT
  REAL(KIND=JPRB) :: ZALPH
  REAL(KIND=JPRB) :: ZDUM
  REAL(KIND=JPRB) :: ZQCR
  REAL(KIND=JPRB) :: ZARG1
  REAL(KIND=JPRB) :: ZARG2
  REAL(KIND=JPRB) :: ZBETA
  REAL(KIND=JPRB) :: ZFACICE
  
  temp (REAL (KIND=JPRB), ZAUTOL, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZAUTOI, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDELT, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZEFFA, (KLON, KLEV))
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: JLEV
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  IF (KIND (ZAUTOL) == 8) THEN
    alloc8 (ZAUTOL)
  ELSE
    IF (KIND (ZAUTOL) == 4) THEN
      alloc4 (ZAUTOL)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZAUTOI) == 8) THEN
    alloc8 (ZAUTOI)
  ELSE
    IF (KIND (ZAUTOI) == 4) THEN
      alloc4 (ZAUTOI)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDELT) == 8) THEN
    alloc8 (ZDELT)
  ELSE
    IF (KIND (ZDELT) == 4) THEN
      alloc4 (ZDELT)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZEFFA) == 8) THEN
    alloc8 (ZEFFA)
  ELSE
    IF (KIND (ZEFFA) == 4) THEN
      alloc4 (ZEFFA)
    ELSE
      STOP 1
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
  !     CHECK RELIABILITY OF INPUT ARGUMENTS.
  
  !     ------------------------------------------------------------------
  
  ! Define threshold for ice autoconversion as a function of temperature
  ! --------------------------------------------------------------------
  
  ZARG1 = 2.0_JPRB*YDML_PHY_MF%YRPHY0%RQICRMAX*(1.0_JPRB - 0.999_JPRB) / (YDML_PHY_MF%YRPHY0%RQICRMAX -  &
  & YDML_PHY_MF%YRPHY0%RQICRMIN) - 1.0_JPRB
  ZARG2 = 2.0_JPRB*(YDML_PHY_MF%YRPHY0%RQICRMAX - 1.5_JPRB*YDML_PHY_MF%YRPHY0%RQICRMIN) / (YDML_PHY_MF%YRPHY0%RQICRMAX -  &
  & YDML_PHY_MF%YRPHY0%RQICRMIN) - 1.0_JPRB
  ZARG1 = 0.5_JPRB*LOG(ABS((1.0_JPRB + ZARG1) / (1.0_JPRB - ZARG1)))
  ZARG2 = 0.5_JPRB*LOG(ABS((1.0_JPRB + ZARG2) / (1.0_JPRB - ZARG2)))
  ZALPH = (ZARG1 - ZARG2) / (YDML_PHY_MF%YRPHY0%RQICRT2 - YDML_PHY_MF%YRPHY0%RQICRT1)
  ZBETA = ZARG1 - YDML_PHY_MF%YRPHY0%RQICRT2*ZALPH
  
  DO JLEV=KTDIA,KLEV
    ZDELT(JLON, JLEV) = PT(JLON, JLEV) - YDCST%RTT
    
    ! Efficiency for ice conversion as a function of temperature.
    ! -----------------------------------------------------------
    ZEFFA(JLON, JLEV) = EXP(YDML_PHY_MF%YRPHY0%RAUTSBET*ZDELT(JLON, JLEV))
    
  END DO
  
  ! ---------------------------------------------------
  ! MICROPHYSICAL AUTOCONVERSION IN THE STRATIFORM CASE
  ! ---------------------------------------------------
  
  DO JLEV=KTDIA,KLEV
    !DEC$ IVDEP
    
    ! Compute in-cloud values
    ! -----------------------
    
    ZQL = MAX(0.0_JPRB, PQL(JLON, JLEV) / PNEBST(JLON, JLEV))
    ZQI = MAX(0.0_JPRB, PQI(JLON, JLEV) / PNEBST(JLON, JLEV))
    
    ! AUTOCONVERSION OF CLOUD LIQUID WATER INTO RAIN
    ! ----------------------------------------------
    
    ZCAUT = YDML_PHY_MF%YRPHY0%RAUTEFR
    ZDUM = (1.0_JPRB - EXP(-ZCAUT*YDML_PHY_MF%YRPHY2%TSPHY))*(ZQL - YDML_PHY_MF%YRPHY0%RQLCR) / YDML_PHY_MF%YRPHY2%TSPHY
    ZAUTOL(JLON, JLEV) = MAX(0.0_JPRB, ZDUM)
    
    ! AUTOCONVERSION OF CLOUD ICE INTO PRECIPITATING ICE
    ! --------------------------------------------------
    
    ZCAUT = YDML_PHY_MF%YRPHY0%RAUTEFS*ZEFFA(JLON, JLEV)
    ZQCR = YDML_PHY_MF%YRPHY0%RQICRMAX - (YDML_PHY_MF%YRPHY0%RQICRMAX - YDML_PHY_MF%YRPHY0%RQICRMIN)*0.5_JPRB*(1.0_JPRB +  &
    & TANH(ZALPH*ZDELT(JLON, JLEV) + ZBETA))
    ZFACICE =  &
    & PLSM(JLON)*PNEIJ(JLON) + (1.0_JPRB - PLSM(JLON))*MAX(0.0_JPRB, SIGN(1.0_JPRB, YDML_PHY_MF%YRPHY1%TMERGL - PTS(JLON)))
    ZQCR = ZQCR*(1.0_JPRB - ZFACICE*(1.0_JPRB - YDML_PHY_MF%YRPHY0%RQICRSN))
    ZDUM = (1.0_JPRB - EXP(-ZCAUT*YDML_PHY_MF%YRPHY2%TSPHY))*(ZQI - ZQCR) / YDML_PHY_MF%YRPHY2%TSPHY
    ZAUTOI(JLON, JLEV) = MAX(0.0_JPRB, ZDUM)
    
    ! TOTAL AUTOCONVERSION TERM
    ! -------------------------
    
    PAUTOL(JLON, JLEV) = ZAUTOL(JLON, JLEV)*PNEBST(JLON, JLEV)
    PAUTOI(JLON, JLEV) = ZAUTOI(JLON, JLEV)*PNEBST(JLON, JLEV)
    
  END DO
  
END SUBROUTINE ACMICRO_OPENACC
