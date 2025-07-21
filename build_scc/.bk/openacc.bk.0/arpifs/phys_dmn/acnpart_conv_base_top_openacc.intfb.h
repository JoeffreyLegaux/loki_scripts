INTERFACE
  SUBROUTINE ACNPART_CONV_BASE_TOP_OPENACC (YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, PNEBC, PAPRSF, PAPHIF, PTOPC,  &
  & YDSTACK)
!$acc routine( ACNPART_CONV_BASE_TOP_OPENACC ) seq
    
    USE MODEL_PHYSICS_MF_MOD, ONLY: MODEL_PHYSICS_MF_TYPE
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE STACK_MOD
    
    USE YOMCST, ONLY: TCST
    
    ! TOP  of convective cloud
    ! If BASE is needed, add PBASEC in the arguments of this subroutine.
    
    ! Interface:
    ! ----------
    ! INPUT:
    !   PNEBC     - convective cloud cover on levels (protected from 0 and 1,
    !               missing in AROME)
    !   PAPRSF    - full level pressure
    ! OUTPUT :
    !!   PBASEC    - BASE of convective cloud [Pa]
    !   PTOPC     - TOP of convective cloud  [Pa]
    
    TYPE(TCST), INTENT(IN) :: YDCST
    TYPE(MODEL_PHYSICS_MF_TYPE), INTENT(IN) :: YDML_PHY_MF
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON
    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
    
    REAL(KIND=JPRB), INTENT(IN) :: PNEBC(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAPRSF(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
    !REAL(KIND=JPRB)  ,INTENT(OUT)   :: PBASEC(KLON) ! Base desactivated
    REAL(KIND=JPRB), INTENT(OUT) :: PTOPC(KLON)
    
    
    TYPE(STACK), INTENT(IN) :: YDSTACK
  END SUBROUTINE ACNPART_CONV_BASE_TOP_OPENACC
END INTERFACE
