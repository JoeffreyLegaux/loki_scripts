INTERFACE
  SUBROUTINE ACCLDIA_OPENACC (YDCST, LDXCLP, LDXTGST, LDXXGST, YDPHY, YDPHY2, YDTOPH, KIDIA, KFDIA, KLON, KLEV, PUCLS, PVCLS,  &
  & PU, PV, PCAPE, PDCAPE, PTKE, PAPHIFM, POROG, PUGST, PVGST, PBLH, KCLPH, YDSTACK)
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
    
    IMPLICIT NONE
    
    TYPE(TPHY), INTENT(IN) :: YDPHY
    TYPE(TPHY2), INTENT(IN) :: YDPHY2
    TYPE(TTOPH), INTENT(IN) :: YDTOPH
    TYPE(TCST), INTENT(IN) :: YDCST
    LOGICAL, INTENT(IN) :: LDXCLP
    LOGICAL, INTENT(IN) :: LDXTGST
    LOGICAL, INTENT(IN) :: LDXXGST
    
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV
    REAL(KIND=JPRB), INTENT(IN) :: PUCLS(KLON), PVCLS(KLON), POROG(KLON)
    REAL(KIND=JPRB), INTENT(OUT) :: PUGST(KLON), PVGST(KLON), PBLH(KLON)
    REAL(KIND=JPRB), INTENT(IN) :: PTKE(KLON, KLEV), PAPHIFM(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PU(KLON, KLEV), PV(KLON, KLEV), PCAPE(KLON), PDCAPE(KLON)
    INTEGER(KIND=JPIM), INTENT(OUT) :: KCLPH(KLON)
    
    
    TYPE(STACK), INTENT(IN) :: YDSTACK
  END SUBROUTINE ACCLDIA_OPENACC
END INTERFACE
