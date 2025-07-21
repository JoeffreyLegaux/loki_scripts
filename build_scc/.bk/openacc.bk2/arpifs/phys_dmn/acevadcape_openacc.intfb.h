INTERFACE
  SUBROUTINE ACEVADCAPE_OPENACC (YDPHY2, YDCST, KIDIA, KFDIA, KLON, KLEV, PFPLSL, PFPLSN, PFPLCL, PFPLCN, PAPRS, PAPRSF, PT,  &
  & PCP, PAPHIF, PAPHI, PDCAPE, YDSTACK)
!$acc routine( ACEVADCAPE_OPENACC )
    
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE YOMCST, ONLY: TCST
    USE YOMPHY2, ONLY: TPHY2
    
    !**** *ACEVADCAPE*  - Compute DCAPE due to precipitation evaporation.
    
    !     PURPOSE.
    !     --------
    !**   INTERFACE.
    !     ----------
    !       *CALL* *ACEVADCAPE*
    
    !        EXPLICIT ARGUMENTS
    !        --------------------
    !            INPUT :
    !        KIDIA   : start of work
    !        KFDIA   : end of work
    !        KLON    : depth of work
    !        KLEV    : number of levels
    !            OUTPUT:
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !     REFERENCE.
    !     ----------
    !        None
    
    !     AUTHOR.
    !     -------
    !        J.M. Piriou.
    
    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 2018-10-13
    !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
    !     ------------------------------------------------------------------
    
    USE STACK_MOD
    
    IMPLICIT NONE
    
    TYPE(TPHY2), INTENT(IN) :: YDPHY2
    TYPE(TCST), INTENT(IN) :: YDCST
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
    REAL(KIND=JPRB), INTENT(IN) :: PFPLSL(KLON, 0:KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PFPLSN(KLON, 0:KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PFPLCL(KLON, 0:KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PFPLCN(KLON, 0:KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAPRS(KLON, 0:KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAPRSF(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PCP(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAPHI(KLON, 0:KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: PDCAPE(KLON)
    
    
    
    !-----------------------------------------------------------------------
    
    TYPE(STACK), INTENT(IN) :: YDSTACK
  END SUBROUTINE ACEVADCAPE_OPENACC
END INTERFACE
