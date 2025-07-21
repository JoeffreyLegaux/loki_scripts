INTERFACE
  SUBROUTINE VERDISINT_OPENACC (YDVFE, YDCVER, CDOPER, CDBC, KPROMA, KST, KEND, KFLEV, PIN, POUT, KCHUNK, YDSTACK)
!$acc routine( VERDISINT_OPENACC )
    
    USE PARKIND1, ONLY: JPIM, JPRB, JPRD
    USE YOMCVER, ONLY: TCVER
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE YOMLUN, ONLY: NULERR
    USE YOMVERT, ONLY: TVFE
    USE STACK_MOD
    
    
    !     ---------------------------------------------------------------------------
    
    IMPLICIT NONE
    TYPE(TVFE), TARGET, INTENT(IN) :: YDVFE
    TYPE(TCVER), INTENT(IN) :: YDCVER
    CHARACTER(LEN=*), INTENT(IN) :: CDOPER, CDBC
    INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA, KST, KEND, KFLEV
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KCHUNK
    REAL(KIND=JPRB), INTENT(IN) :: PIN(KPROMA, 0:KFLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: POUT(KPROMA, KFLEV + 1)
    
    !     ------------------------------------------------------------------
    
    
    !     ------------------------------------------------------------------
    
    TYPE(STACK), INTENT(IN) :: YDSTACK
  END SUBROUTINE VERDISINT_OPENACC
END INTERFACE
