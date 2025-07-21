SUBROUTINE ACCLDIA_OPENACC (YDCST, LDXCLP, LDXTGST, LDXXGST, YDPHY, YDPHY2, YDTOPH, KIDIA, KFDIA, KLON, KLEV, PUCLS, PVCLS, PU,  &
& PV, PCAPE, PDCAPE, PTKE, PAPHIFM, POROG, PUGST, PVGST, PBLH, KCLPH, YDSTACK)
  
!$acc routine( ACCLDIA_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMPHY, ONLY: TPHY
  USE YOMPHY2, ONLY: TPHY2
  USE YOMCST, ONLY: TCST
  USE YOMTOPH, ONLY: TTOPH
  
  !**** *ACCLDIA*  - Compute some PBL Diags  -
  
  !     PURPOSE.
  !     --------
  !        To Compute some PBl diags
  !                   * Wind gusts from PBL wind and Tke
  !                   * Height of PBL
  !**   INTERFACE.
  !     ----------
  !       *CALL* *ACCLDIA*
  
  !        EXPLICIT ARGUMENTS
  !        --------------------
  !            INPUT :
  !        KIDIA   : start of work
  !        KFDIA   : end of work
  !        KLON    : depth of work
  !        KLEV    : number of levels
  !        PUCLS   : x-CLS wind              (KLON)
  !        PVCLS   : y-CLS wind              (KLON)
  !        PU      : x-wind                  (KLON,KLEV)
  !        PV      : y-wind                  (KLON,KLEV)
  !        PTKE    : TKE                     (KLON,KLEV)
  !        PCAPE   : CAPE                    (KLON)
  !        PDCAPE  : downward CAPE           (KLON)
  !        PAPHIFM : full level geopotential (KLON,KLEV)
  !        POROG   : orography times g       (KLON)
  !            OUTPUT:
  !        PUGST   : x-wind gust             (KLON)
  !        PVGST   : y-wind gust             (KLON)
  !        PBLH    : PBL height              (KLON)
  !        KCLPH   : level of PBL            (KLON)
  
  !        IMPLICIT ARGUMENTS
  !        --------------------
  !           NONE
  
  !     METHOD.
  !     -------
  !        Consider That winds gust are // to cls wind. Consider Tke to
  !        increment the cls wind to obtain wind gust. Formulation first
  !        implemented and tested in the meso-Nh diags for the 1999
  !        storms.
  
  !     EXTERNALS.
  !     ----------
  
  !     REFERENCE.
  !     ----------
  !        None
  
  !     AUTHOR.
  !     -------
  !        Gwenaelle Hello *METEO-FRANCE*
  
  !     MODIFICATIONS.
  !     --------------
  !        ORIGINAL : 07-04-18
  !        S. Riette    : 2009-03-25 HTKERAF and FACRAF are introduced
  !        E. Bazile et Y. Seity : 2009-07-16 vectorisation
  !        E. Bazile et Y. Seity : 2010-05-16 add PBL Height
  !        2011-06: M. Jerczynski - some cleaning to meet norms
  !        2018-06: J.M. Piriou : convective wind gusts (FACRAFCV, GCAPERAF, FACRAFDCAPE).
  !        2020-10: R. Brozkova : usage with TOUCANS not touching PBL height
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  
  ! End modifications
  !-------------------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TPHY), INTENT(IN) :: YDPHY
  TYPE(TPHY2), INTENT(IN) :: YDPHY2
  TYPE(TTOPH), INTENT(IN) :: YDTOPH
  TYPE(TCST), INTENT(IN) :: YDCST
  LOGICAL, INTENT(IN) :: LDXCLP
  LOGICAL, INTENT(IN) :: LDXTGST
  LOGICAL, INTENT(IN) :: LDXXGST
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  REAL(KIND=JPRB), INTENT(IN) :: PUCLS(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PVCLS(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: POROG(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PUGST(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PVGST(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PBLH(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTKE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIFM(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PV(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PCAPE(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PDCAPE(KLON)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KCLPH(KLON)
  
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: ILEVM1
  INTEGER(KIND=JPIM) :: JJLEVM1
  INTEGER(KIND=JPIM) :: JLEVM1
  temp (INTEGER (KIND=JPIM), ICM, (KLON, KLEV))
  INTEGER(KIND=JPIM) :: ILEVBI
  REAL(KIND=JPRB) :: ZCVGUSTS
  REAL(KIND=JPRB) :: ZALPHA
  REAL(KIND=JPRB) :: ZVCLS
  REAL(KIND=JPRB) :: ZZ
  REAL(KIND=JPRB) :: ZEPS
  REAL(KIND=JPRB) :: ZECTBLH
  REAL(KIND=JPRB) :: ZHBOT
  REAL(KIND=JPRB) :: ZHTOP
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  IF (KIND (ICM) == 8) THEN
    alloc8 (ICM)
  ELSE
    IF (KIND (ICM) == 4) THEN
      alloc4 (ICM)
    ELSE
      STOP 1
    END IF
  END IF
  
  !      ----------------------------------------------------------------
  
  ! 1. Computation of gusts
  
  IF ((LDXTGST .or. LDXXGST) .and. YDPHY2%LRAFTKE) THEN
    ZEPS = 1.E-6_JPRB
    JLEVM1 = 0
    
    ! Recherche du niveau JLEV juste au dessus de HTKERAF
    DO JLEV=KLEV,1,-1
      ZZ = PAPHIFM(JLON, JLEV) - POROG(JLON) - YDPHY2%HTKERAF*YDCST%RG
      JJLEVM1 = MAX(0._JPRB, SIGN(1._JPRB, ZZ))*JLEV
      JLEVM1 = MAX(JJLEVM1, JLEVM1)
    END DO
    
    
    !  Cas ou HTKERAF est en dessous du niveau KLEV donc extrapolation
    JLEVM1 = MIN(JLEVM1, KLEV - 1)
    ! Turbulent gusts.
    ZVCLS = MAX(ZEPS, PUCLS(JLON)**2 + PVCLS(JLON)**2)
    ZALPHA = PTKE(JLON, JLEVM1)*(YDPHY2%HTKERAF*YDCST%RG + POROG(JLON) - PAPHIFM(JLON, JLEVM1 + 1)) - PTKE(JLON, JLEVM1 + 1) &
    & *(YDPHY2%HTKERAF*YDCST%RG + POROG(JLON) - PAPHIFM(JLON, JLEVM1))
    ZALPHA = MAX(ZEPS, ZALPHA / (PAPHIFM(JLON, JLEVM1) - PAPHIFM(JLON, JLEVM1 + 1)))
    ZALPHA = ZALPHA / ZVCLS
    ! Convective gusts, if CAPE > given threshold.
    ZCVGUSTS = MAX(0._JPRB, SQRT(PU(JLON, YDTOPH%NT850)**2 + PV(JLON, YDTOPH%NT850)**2) - SQRT(PU(JLON, YDTOPH%NT950)**2 +  &
    & PV(JLON, YDTOPH%NT950)**2))*MAX(0._JPRB, SIGN(1._JPRB, PCAPE(JLON) - YDPHY2%GCAPERAF))*MIN(1._JPRB, SQRT(MAX(0._JPRB,  &
    & PCAPE(JLON) / (5._JPRB*YDPHY2%GCAPERAF))))
    ! Gust from 3 processes: turbulence (FACRAF), convective transport (FACRAFCV),
    ! precipitation evaporation in the PBL (FACRAFDCAPE).
    ZALPHA = 1.0_JPRB + YDPHY2%FACRAF*SQRT(ZALPHA) + YDPHY2%FACRAFCV*ZCVGUSTS / SQRT(ZVCLS) +  &
    & YDPHY2%FACRAFDCAPE*SQRT(2._JPRB*MAX(0._JPRB, -PDCAPE(JLON))) / SQRT(ZVCLS)
    PUGST(JLON) = ZALPHA*PUCLS(JLON)
    PVGST(JLON) = ZALPHA*PVCLS(JLON)
  END IF
  
  ! 2. Computation of PBL Height
  
  ! do not overwrite diagnosed PBL height for TOUCANS (LCOEFKTKE=.T.)
  IF (LDXCLP .and. .not.YDPHY%LCOEFKTKE) THEN
    
    ! CALCUL DE LA HAUTEUR DE COUCHE LIMITE:
    ! PREMIER NIVEAU EN PARTANT DE LA SURFACE OU LA TKE <0.01
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Compute the INDEX array ILEVBI from the lowest
    ! half level (KLEV-1) to the "Top-PBL" half level :
    !- - - - - - - - - - - - - - - - - - - - - - - - - -
    ! ILEVBI(JLON) = 1 from KLEV-1 to the "Top-PBL"
    ! ILEVBI(JLON) = 0 above the "Top-PBL" half level
    !- - - - - - - - - - - - - - - - - - - - - - - - - -
    ! the case "MAX(JLEV, ILEVBI(JLON))" avoid the
    ! detection of the other "PBL" located above
    ! the first one close to the ground.
    !- - - - - - - - - - - - - - - - - - - - - - - - - -
    ILEVBI = 0
    ICM(JLON, :) = 0
    ZECTBLH = 0.01_JPRB
    !NTCOET = 1 by default in the setup
    DO JLEV=YDTOPH%NTCOET,KLEV
      ICM(JLON, JLEV) = INT(MAX(0.0_JPRB, SIGN(1.0_JPRB, PTKE(JLON, JLEV) - ZECTBLH)))
    END DO
    DO JLEV=KLEV,YDTOPH%NTCOET,-1
      ILEVM1 = MAX(YDTOPH%NTCOET, JLEV - 1)
      IF (ICM(JLON, JLEV) == 1 .and. ICM(JLON, ILEVM1) == 0) THEN
        ILEVBI = MAX(JLEV, ILEVBI)
      END IF
    END DO
    ILEVBI = ILEVBI*MAX(ICM(JLON, KLEV), ICM(JLON, KLEV - 1))
    IF (ICM(JLON, KLEV) == 0 .and. ILEVBI == 0) ILEVBI = KLEV
    KCLPH(JLON) = MAX(1, ILEVBI)
    IF (ILEVBI > 1 .and. ILEVBI < KLEV) THEN
      ZHBOT = (PAPHIFM(JLON, ILEVBI) - POROG(JLON)) / YDCST%RG
      ZHTOP = (PAPHIFM(JLON, ILEVBI - 1) - POROG(JLON)) / YDCST%RG
      PBLH(JLON) = ZHBOT + (ZHTOP - ZHBOT) / (PTKE(JLON, ILEVBI - 1) - PTKE(JLON, ILEVBI))*(ZECTBLH - PTKE(JLON, ILEVBI))
    ELSE IF (ILEVBI == KLEV) THEN
      PBLH(JLON) = (PAPHIFM(JLON, KLEV) - POROG(JLON)) / YDCST%RG
    ELSE IF (ILEVBI == 0) THEN
      PBLH(JLON) = (PAPHIFM(JLON, YDTOPH%NTCOET) - POROG(JLON)) / YDCST%RG
    ELSE
      PBLH(JLON) = (PAPHIFM(JLON, ILEVBI) - POROG(JLON)) / YDCST%RG
    END IF
    
    !WRITE(NULOUT,*)'sous accldia PBLH',MINVAL(PBLH),MAXVAL(PBLH)
    
  END IF
  
  
END SUBROUTINE ACCLDIA_OPENACC
