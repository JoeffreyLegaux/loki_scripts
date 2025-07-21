INTERFACE
  SUBROUTINE RADOZCMF_OPENACC (YDCST, YDEOZOC, KIDIA, KFDIA, KLON, KLEV, PAPRS, PGEMU, POZON, YDSTACK)
!$acc routine( RADOZCMF_OPENACC )
    
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    
    USE YOMCST, ONLY: TCST
    USE YOEOZOC, ONLY: TEOZOC
    
    USE STACK_MOD
    
    IMPLICIT NONE
    
    TYPE(TCST), INTENT(IN) :: YDCST
    TYPE(TEOZOC), INTENT(IN) :: YDEOZOC
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
    REAL(KIND=JPRB), INTENT(IN) :: PAPRS(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(IN) :: PGEMU(KLON)
    REAL(KIND=JPRB), INTENT(OUT) :: POZON(KLON, KLEV)
    !     -----------------------------------------------------------------
    
    !*       0.1   ARGUMENTS.
    !              ----------
    
    !     -----------------------------------------------------------------
    
    !*       0.2   LOCAL ARRAYS.
    !              -------------
    
    
    
    TYPE(STACK), INTENT(IN) :: YDSTACK
  END SUBROUTINE RADOZCMF_OPENACC
END INTERFACE
