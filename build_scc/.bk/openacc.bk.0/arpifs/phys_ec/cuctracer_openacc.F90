SUBROUTINE CUCTRACER_OPENACC (YDCST, YDCUMFS, YDECUMF2, YDECUMF, KIDIA, KFDIA, KLON, KLEV, KTRAC, KCTOP, KDTOP, KTYPE, LDCUM,  &
& LDDRAF, PTSPHY, PAPH, PAP, PMFU, PMFD, PMFUO, PMFDO, PUDRATE, PDDRATE, PDMFUP, PDMFDP, PCEN, PTENC, PSCAV, YDSTACK)
  
  !**** *CUCTRACER* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
  !                   IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES
  !                              SCAVENGING IS DONE IF SCAVENGING COEFFICIENT>0
  
  !          P.BECHTOLD        E.C.M.W.F.              11/02/2004
  
  !**   INTERFACE.
  !     ----------
  
  !          *CUTRACER* IS CALLED FROM *CUMASTR*
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KTRAC*        NUMBER OF CHEMICAL TRACERS
  
  !    *KCTOP*        CLOUD TOP  LEVEL
  !    *KDTOP*        DOWNDRAFT TOP LEVEL
  !    *KTYPE*        TYPE OF CONVECTION
  !                       1 = PENETRATIVE CONVECTION
  !                       2 = SHALLOW CONVECTION
  !                       3 = MIDLEVEL CONVECTION
  
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  !    *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTSPHY*       PHYSICS TIME-STEP                              S
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
  !    *PAPH          PROVISIONAL PRESSURE ON FULL LEVELS            PA
  !    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
  !    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
  !    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
  !    *PUDRATE*      UPDRAFT DETRAINMENT                           KG/(M2*S)
  !    *PDDRATE*      DOWNDRAFT DETRAINMENT                         KG/(M2*S)
  !    *PDMFUP*       UPDRAUGT PRECIP (PRODUCTION) IN LAYER         KG/(M2*S)
  !    *PDMFDP*       DDRAUGHT PRECIP (PRODUCTION) IN LAYER         KG/(M2*S)
  !    *PSCAV*        WET SCAVENGING COEFFICIENTS                   UNITLESS
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S
  
  !          METHOD
  !          -------
  !     EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
  !     DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT
  
  !     FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
  !     FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
  !                            TO CORRECT TENDENCIES BELOW CLOUD BASE
  
  
  !------------------------------------------------------------------------------------
  !     COMMENTS FOR OFFLINE USERS IN CHEMICAL TRANSPORT MODELS
  !     (i.e. reading mass fluxes and detrainment rates from ECMWF archive:
  !      ------------------------------------------------------------------
  !     KCTOP IS FIRST LEVEL FROM TOP WHERE PMFU>0
  !     KDTOP IS FIRST LEVEL FROM TOP WHERE PMFD<0
  !     ATTENTION: ON ARCHIVE DETRAINMENT RATES HAVE UNITS KG/(M3*S), SO FOR USE
  !                IN CURRENT ROUTINE YOU HAVE TO MULTIPLY ARCHIVED VALUES BY DZ
  !     LDCUM  IS TRUE IF CONVECTION EXISTS, i.e. IF PMFU>0 IN COLUMN OR IF
  !                       KCTOP>0 AND KCTOP<KLEV
  !     LDDRAF IS TRUE IF DOWNDRAUGHTS EXIST IF PMFD<0 IN COLUMN OR IF
  !                       KDTOP>0 AND KDTOP<KLEV
  !     IF MASSFLUX SATISFIES CFL CRITERIUM M<=DP/Dt IT IS SUFFICIENT TO
  !     ONLY CONSIDER EXPLICIT SOLUTION (RMFSOLCT=0), IN THIS CASE
  !     YOU CAN IGNORE IMPLICIT PART 7.0 OF CURRENT ROUTINE
  !------------------------------------------------------------------------------------
  
  !          EXTERNALS
  !          ---------
  !          CUBIDIAG
  
  !          MODIFICATIONS
  !          -------------
  !        M.Hamrud      01-Oct-2003 CY28 Cleaning
  !        P.Bechtold    07-Oct-2009 add scavenging
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  
  !----------------------------------------------------------------------
  
!$acc routine( CUCTRACER_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOECUMF, ONLY: TECUMF
  USE YOECUMF2, ONLY: TECUMF2
  USE YOMCUMFS, ONLY: TCUMFS
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TCUMFS), INTENT(IN) :: YDCUMFS
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TECUMF2), INTENT(IN) :: YDECUMF2
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KTRAC
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KDTOP(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(KLON)
  LOGICAL, INTENT(IN) :: LDCUM(KLON)
  LOGICAL, INTENT(IN) :: LDDRAF(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFUO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFDO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUDRATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDDRATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDMFUP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDMFDP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PCEN(KLON, KLEV, KTRAC)
  REAL(KIND=JPRB), INTENT(IN) :: PSCAV(KTRAC)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENC(KLON, KLEV, KTRAC)
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IM
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  INTEGER(KIND=JPIM) :: JN
  
  REAL(KIND=JPRB) :: ZZP
  REAL(KIND=JPRB) :: ZMFA
  REAL(KIND=JPRB) :: ZIMP
  REAL(KIND=JPRB) :: ZERATE
  REAL(KIND=JPRB) :: ZPOSI
  REAL(KIND=JPRB) :: ZTSPHY
  REAL(KIND=JPRB) :: ZC
  REAL(KIND=JPRB) :: ZADVW
  REAL(KIND=JPRB) :: ZRMFCMIN
  REAL(KIND=JPRB) :: ZRMFSOLCT
  
  !     ALLOCATABLE ARAYS
  temp (REAL (KIND=JPRB), ZCEN, (KLON, KLEV, KTRAC))
  temp (REAL (KIND=JPRB), ZCU, (KLON, KLEV, KTRAC))
  temp (REAL (KIND=JPRB), ZCD, (KLON, KLEV, KTRAC))
  temp (REAL (KIND=JPRB), ZTENC, (KLON, KLEV, KTRAC))
  temp (REAL (KIND=JPRB), ZMFC, (KLON, KLEV, KTRAC))
  temp (REAL (KIND=JPRB), ZDP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZB, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZR1, (KLON, KLEV))
  temp (LOGICAL, LLCUMASK, (KLON, KLEV))
  temp (LOGICAL, LLCUMBAS, (KLON, KLEV))
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cubidiag_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ZCEN)
  alloc (ZCU)
  alloc (ZCD)
  alloc (ZTENC)
  alloc (ZMFC)
  alloc (ZDP)
  alloc (ZB)
  alloc (ZR1)
  alloc (LLCUMASK)
  alloc (LLCUMBAS)
  JLON = KIDIA
  !----------------------------------------------------------------------
  
  IF (YDCUMFS%LECUMFS) THEN
    ZRMFSOLCT = YDECUMF2%RMFSOLCT2
    ZRMFCMIN = YDECUMF2%RMFCMIN2
  ELSE
    ZRMFSOLCT = YDECUMF%RMFSOLCT
    ZRMFCMIN = YDECUMF%RMFCMIN
  END IF
  
  ZIMP = 1.0_JPRB - ZRMFSOLCT
  ZTSPHY = 1.0_JPRB / PTSPHY
  
  ZADVW = 0.0_JPRB
  IF (KTYPE(JLON) == 1) ZADVW = YDECUMF%RMFADVW
  
  ! Initialize Cumulus mask + some setups
  
  DO JK=2,KLEV
    LLCUMASK(JLON, JK) = .false.
    IF (LDCUM(JLON)) THEN
      ZDP(JLON, JK) = YDCST%RG / (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
      IF (JK >= KCTOP(JLON) - 1) THEN
        LLCUMASK(JLON, JK) = .true.
      END IF
    END IF
  END DO
  !----------------------------------------------------------------------
  
  DO JN=1,KTRAC
    
    !*    1.0          DEFINE TRACERS AT HALF LEVELS
    !                  -----------------------------
    
    DO JK=2,KLEV
      IK = JK - 1
      ZCEN(JLON, JK, JN) = PCEN(JLON, JK, JN)
      ZCD(JLON, JK, JN) = PCEN(JLON, IK, JN)
      ZCU(JLON, JK, JN) = PCEN(JLON, IK, JN)
      ZMFC(JLON, JK, JN) = 0.0_JPRB
      ZTENC(JLON, JK, JN) = 0.0_JPRB
    END DO
    ZCU(JLON, KLEV, JN) = PCEN(JLON, KLEV, JN)
    
    !*    2.0          COMPUTE UPDRAFT VALUES
    !                  ----------------------
    
    DO JK=KLEV - 1,3,-1
      IK = JK + 1
      IF (LLCUMASK(JLON, JK)) THEN
        ZERATE = PMFU(JLON, JK) - PMFU(JLON, IK) + PUDRATE(JLON, JK)
        ZMFA = 1.0_JPRB / MAX(ZRMFCMIN, PMFU(JLON, JK))
        IF (JK >= KCTOP(JLON)) THEN
          ZCU(JLON, JK, JN) = (PMFU(JLON, IK)*ZCU(JLON, IK, JN) + ZERATE*PCEN(JLON, JK, JN) - (PUDRATE(JLON, JK) + PDMFUP(JLON,  &
          & JK)*PSCAV(JN))*ZCU(JLON, IK, JN))*ZMFA
          ! if you have a source term dc/dt=dcdt write
          !             ZCU(JL,JK,JN)=( PMFU(JL,IK)*ZCU(JL,IK,JN)+ZERATE*PCEN(JL,JK,JN) &
          !                           -PUDRATE(JL,JK)*ZCU(JL,IK,JN) )*ZMFA
          !                           +dcdt(jl,ik,jn)*ptsphy
        END IF
      END IF
    END DO
    
    
    !*    3.0          COMPUTE DOWNDRAFT VALUES
    !                  ------------------------
    
    DO JK=3,KLEV
      IK = JK - 1
      IF (LDDRAF(JLON) .and. JK == KDTOP(JLON)) THEN
        !Nota: in order to avoid final negative Tracer values at LFS the allowed value of ZCD
        !      depends on the jump in mass flux at the LFS
        !ZCD(JL,JK,JN)=0.5_JPRB*ZCU(JL,JK,JN)+0.5_JPRB*PCEN(JL,IK,JN)
        ZCD(JLON, JK, JN) = 0.1_JPRB*ZCU(JLON, JK, JN) + 0.9_JPRB*PCEN(JLON, IK, JN)
      ELSE IF (LDDRAF(JLON) .and. JK > KDTOP(JLON)) THEN
        ZERATE = -PMFD(JLON, JK) + PMFD(JLON, IK) + PDDRATE(JLON, JK)
        ZMFA = 1._JPRB / MIN(-ZRMFCMIN, PMFD(JLON, JK))
        ZCD(JLON, JK, JN) = (PMFD(JLON, IK)*ZCD(JLON, IK, JN) - ZERATE*PCEN(JLON, IK, JN) + (PDDRATE(JLON, JK) + PDMFDP(JLON, JK) &
        & *PSCAV(JN))*ZCD(JLON, IK, JN))*ZMFA
        ! if you have a source term dc/dt=dcdt write
        !             ZCD(JL,JK,JN)=( PMFD(JL,IK)*ZCD(JL,IK,JN)-ZERATE*PCEN(JL,IK,JN) &
        !                           &+PDDRATE(JL,JK)*ZCD(JL,IK,JN) &
        !                           &+dcdt(jl,ik,jn)*ptsphy
      END IF
    END DO
    
    ! In order to avoid negative Tracer at KLEV adjust ZCD
    JK = KLEV
    IK = JK - 1
    IF (LDDRAF(JLON)) THEN
      ZPOSI = -ZDP(JLON, JK)*(PMFU(JLON, JK)*ZCU(JLON, JK, JN) + PMFD(JLON, JK)*ZCD(JLON, JK, JN) - (PMFU(JLON, JK) + PMFD(JLON,  &
      & JK))*PCEN(JLON, IK, JN))
      IF (PCEN(JLON, JK, JN) + ZPOSI*PTSPHY < 0.0_JPRB) THEN
        ZMFA = 1._JPRB / MIN(-ZRMFCMIN, PMFD(JLON, JK))
        ZCD(JLON, JK, JN) = ((PMFU(JLON, JK) + PMFD(JLON, JK))*PCEN(JLON, IK, JN) - PMFU(JLON, JK)*ZCU(JLON, JK, JN) + PCEN(JLON, &
        &  JK, JN) / (PTSPHY*ZDP(JLON, JK)))*ZMFA
      END IF
    END IF
    
  END DO
  
  !----------------------------------------------------------------------
  
  DO JN=1,KTRAC
    
    !*    4.0          COMPUTE FLUXES
    !                  --------------
    
    DO JK=2,KLEV
      IK = JK - 1
      IF (LLCUMASK(JLON, JK)) THEN
        ZMFA = PMFU(JLON, JK) + PMFD(JLON, JK)
        ZMFC(JLON, JK, JN) = PMFU(JLON, JK)*ZCU(JLON, JK, JN) + PMFD(JLON, JK)*ZCD(JLON, JK, JN) - ZIMP*ZMFA*ZCEN(JLON, IK, JN)
      END IF
    END DO
    
    !*    5.0          COMPUTE TENDENCIES = RHS
    !                  ------------------------
    
    DO JK=2,KLEV - 1
      IK = JK + 1
      IF (LLCUMASK(JLON, JK)) THEN
        ZTENC(JLON, JK, JN) = ZDP(JLON, JK)*(ZMFC(JLON, IK, JN) - ZMFC(JLON, JK, JN))
        !ZTENC(JL,JK,JN)=ZTENC(JL,JK,JN)-ZDP(JL,JK)*PDMFUP(JL,JK)*PSCAV(JN)*ZCU(JL,JK,JN)
      END IF
    END DO
    
    JK = KLEV
    IF (LDCUM(JLON)) THEN
      ZTENC(JLON, JK, JN) = -ZDP(JLON, JK)*ZMFC(JLON, JK, JN)
      !ZTENC(JL,JK,JN)=ZTENC(JL,JK,JN)-ZDP(JL,JK)*PDMFUP(JL,JK)*PSCAV(JN)*ZCU(JL,JK,JN)
    END IF
    
  END DO
  
  IF (ZRMFSOLCT == 0.0_JPRB) THEN
    
    
    !*    6.0          UPDATE TENDENCIES
    !                  -----------------
    
    DO JN=1,KTRAC
      DO JK=2,KLEV
        IF (LLCUMASK(JLON, JK)) THEN
          PTENC(JLON, JK, JN) = PTENC(JLON, JK, JN) + ZTENC(JLON, JK, JN)
        END IF
      END DO
    END DO
    
  ELSE
    
    !---------------------------------------------------------------------------
    
    !*    7.0          IMPLICIT SOLUTION
    !                  -----------------
    
    ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
    ! reuse ZMFC=A and ZB=B;
    ! ZTENC corresponds to the RHS ("constants") of the equation
    ! The solution is in ZR1
    
    LLCUMBAS(JLON, :) = .false.
    ZB(JLON, :) = 1._JPRB
    
    DO JN=1,KTRAC
      
      ! Fill vectors A, B and RHS
      
      DO JK=2,KLEV
        IK = JK + 1
        IM = JK - 1
        ! LLCUMBAS(JL,JK)=LLCUMASK(JL,JK).AND.JK<=KCBOT(JL)
        LLCUMBAS(JLON, JK) = LLCUMASK(JLON, JK)
        IF (LLCUMBAS(JLON, JK)) THEN
          ZZP = ZRMFSOLCT*ZDP(JLON, JK)*PTSPHY
          ZMFC(JLON, JK, JN) = -ZZP*(PMFU(JLON, JK) + PMFD(JLON, JK))
          IF (JK < KLEV) THEN
            ZB(JLON, JK) = 1.0_JPRB + ZZP*(PMFU(JLON, IK) + PMFD(JLON, IK))
          ELSE
            ZB(JLON, JK) = 1.0_JPRB
          END IF
          ZZP = YDCST%RG*(PMFUO(JLON, JK) + YDECUMF%RMFADVWDD*PMFDO(JLON, JK)) / (PAP(JLON, JK) - PAP(JLON, IM))*PTSPHY*ZADVW
          ZC = ZZP*(PCEN(JLON, IM, JN) - PCEN(JLON, JK, JN))
          ZTENC(JLON, JK, JN) = ZTENC(JLON, JK, JN)*PTSPHY + PCEN(JLON, JK, JN) - ZC
        END IF
      END DO
      
      CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFC(:, :, JN), ZB, ZTENC(:, :, JN), ZR1, YDSTACK=YLSTACK)
      
      ! Compute tendencies
      
      DO JK=2,KLEV
        IF (LLCUMBAS(JLON, JK)) THEN
          PTENC(JLON, JK, JN) = PTENC(JLON, JK, JN) + (ZR1(JLON, JK) - PCEN(JLON, JK, JN))*ZTSPHY
          !  for implicit solution including tendency source term
          !  PTENC(JL,JK,JN)=(ZR1(JL,JK)-PCEN(JL,JK,JN))*ZTSPHY
        END IF
      END DO
      
    END DO
    
  END IF
  !---------------------------------------------------------------------------
  
END SUBROUTINE CUCTRACER_OPENACC
