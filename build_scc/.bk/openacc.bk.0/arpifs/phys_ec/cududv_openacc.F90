SUBROUTINE CUDUDV_OPENACC (YDCST, YDECUMF, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, KTOPM2, KTYPE, KCBOT, KCTOP, LDCUM, PTSPHY,  &
& PAPH, PAP, PUEN, PVEN, PMFU, PMFD, PMFUO, PMFDO, PUU, PUD, PVU, PVD, PGP2DSPP, PTENU, PTENV, YDSTACK)
  
  !**** *CUDUDV* - UPDATES U AND V TENDENCIES,
  !                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
  
  !**   INTERFACE.
  !     ----------
  
  !          *CUDUDV* IS CALLED FROM *CUMASTR*
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KTYPE*        TYPE OF CONVECTION
  !                       1 = PENETRATIVE CONVECTION
  !                       2 = SHALLOW CONVECTION
  !                       3 = MIDLEVEL CONVECTION
  !    *KCBOT*        CLOUD BASE LEVEL
  !    *KCTOP*        CLOUD TOP LEVEL
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
  !    *PAP *         PROVISIONAL PRESSURE ON FULL LEVELS            PA
  !    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
  !    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
  !    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
  !    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
  !    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
  !    *PUD*          U-VELOCITY IN DOWNDRAFTS                       M/S
  !    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
  !    *PVD*          V-VELOCITY IN DOWNDRAFTS                       M/S
  !    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
  !    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
  
  !            METHOD
  !            -------
  !       EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
  !       DEPENDING ON VALUE OF RMFSOLUV:
  !       0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT
  
  !       FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
  !       FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
  !                              TO CORRECT TENDENCIES BELOW CLOUD BASE
  
  !            EXTERNALS
  !            ---------
  !            CUBIDIAG
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      P.BECHTOLD        E.C.M.W.F.    11/02/05 IMPLICIT SOLVER
  !      M. Leutbecher&S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
  !      S.-J. Lock    Jul 2017    Extended SPP perturbations to include shallow convection
  !      L. Descamps 2020-02-25 Introduced perturbed parameter option for PEARP (LPERTPAR)
  !      M. Leutbecher & S. Lang Oct 2020    SPP abstraction and revision
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUDUDV_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOECUMF, ONLY: TECUMF
  USE SPP_MOD, ONLY: TSPP_CONFIG
  USE YOMPERTPAR, ONLY: LPERTPAR, LPERT_CUDUV, CUDU_MOD, CUDV_MOD
  USE SPP_GEN_MOD, ONLY: SPP_PERT
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TSPP_CONFIG), INTENT(IN) :: YDSPP_CONFIG
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KTOPM2
  INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP(KLON)
  LOGICAL, INTENT(IN) :: LDCUM(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFUO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFDO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PUD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PVD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENV(KLON, KLEV)
  temp (REAL (KIND=JPRB), ZMFUU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDU, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUV, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDV, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZADVW
  
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IM
  INTEGER(KIND=JPIM) :: IKB
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  INTEGER(KIND=JPIM) :: IPCUDU
  INTEGER(KIND=JPIM) :: IPCUDV
  INTEGER(KIND=JPIM) :: IPCUDUS
  INTEGER(KIND=JPIM) :: IPCUDVS
  INTEGER(KIND=JPIM) :: IPN1
  INTEGER(KIND=JPIM) :: IPN2
  TYPE(SPP_PERT) :: PN1  ! SPP pertn. config.
  TYPE(SPP_PERT) :: PN2  ! SPP pertn. config.
  
  REAL(KIND=JPRB) :: ZZP
  REAL(KIND=JPRB) :: ZIMP
  REAL(KIND=JPRB) :: ZTSPHY
  REAL(KIND=JPRB) :: ZU
  REAL(KIND=JPRB) :: ZV
  REAL(KIND=JPRB) :: ZMDU
  REAL(KIND=JPRB) :: ZMDV
  REAL(KIND=JPRB) :: ZLIMN2
  REAL(KIND=JPRB) :: ZXU
  REAL(KIND=JPRB) :: ZXV
  REAL(KIND=JPRB) :: ZXU_MAG
  REAL(KIND=JPRB) :: ZXV_MAG
  REAL(KIND=JPRB) :: ZN2
  REAL(KIND=JPRB) :: ZFAC
  REAL(KIND=JPRB) :: ZRDU
  REAL(KIND=JPRB) :: ZRDV
  temp (REAL (KIND=JPRB), ZDUDT, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDVDT, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZB, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZR1, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZR2, (KLON, KLEV))
  temp (LOGICAL, LLCUMBAS, (KLON, KLEV))
  LOGICAL :: LLPPAR_CUDUV
  LOGICAL :: LLPERT_CUDUDV  ! SPP perturbation on?
  LOGICAL :: LLPERT_CUDUDVS  ! SPP perturbation on?
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cubidiag_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ZMFUU)
  alloc (ZMFDU)
  alloc (ZMFUV)
  alloc (ZMFDV)
  alloc (ZDUDT)
  alloc (ZDVDT)
  alloc (ZDP)
  alloc (ZB)
  alloc (ZR1)
  alloc (ZR2)
  alloc (LLCUMBAS)
  JLON = KIDIA
  !----------------------------------------------------------------------
  
  ZIMP = 1.0_JPRB - YDECUMF%RMFSOLUV
  ZTSPHY = 1.0_JPRB / PTSPHY
  
  LLPPAR_CUDUV = .false.
  LLPERT_CUDUDV = .false.
  LLPERT_CUDUDVS = .false.
  
  ZADVW = 0.0_JPRB
  IF (KTYPE(JLON) == 1) ZADVW = YDECUMF%RMFADVW
  
  
  DO JK=1,KLEV
    IF (LDCUM(JLON)) THEN
      ZDP(JLON, JK) = YDCST%RG / (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
    END IF
  END DO
  
  !----------------------------------------------------------------------
  
  
  !*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
  !                  ----------------------------------------------
  
  DO JK=KTOPM2,KLEV
    IK = JK - 1
    IF (LDCUM(JLON)) THEN
      ZMFUU(JLON, JK) = PMFU(JLON, JK)*(PUU(JLON, JK) - ZIMP*PUEN(JLON, IK))
      ZMFUV(JLON, JK) = PMFU(JLON, JK)*(PVU(JLON, JK) - ZIMP*PVEN(JLON, IK))
      ZMFDU(JLON, JK) = PMFD(JLON, JK)*(PUD(JLON, JK) - ZIMP*PUEN(JLON, IK))
      ZMFDV(JLON, JK) = PMFD(JLON, JK)*(PVD(JLON, JK) - ZIMP*PVEN(JLON, IK))
    END IF
  END DO
  
  ! linear fluxes below cloud
  IF (YDECUMF%RMFSOLUV == 0.0_JPRB) THEN
    DO JK=KTOPM2,KLEV
      !DIR$ IVDEP
      !OCL NOVREC
      IF (LDCUM(JLON) .and. JK > KCBOT(JLON)) THEN
        IKB = KCBOT(JLON)
        ZZP = ((PAPH(JLON, KLEV + 1) - PAPH(JLON, JK)) / (PAPH(JLON, KLEV + 1) - PAPH(JLON, IKB)))
        IF (KTYPE(JLON) == 3) THEN
          ZZP = ZZP*ZZP
        END IF
        ZMFUU(JLON, JK) = ZMFUU(JLON, IKB)*ZZP
        ZMFUV(JLON, JK) = ZMFUV(JLON, IKB)*ZZP
        ZMFDU(JLON, JK) = ZMFDU(JLON, IKB)*ZZP
        ZMFDV(JLON, JK) = ZMFDV(JLON, IKB)*ZZP
      END IF
    END DO
  END IF
  
  !*    1.2          COMPUTE TENDENCIES
  !                  ------------------
  
  DO JK=KTOPM2,KLEV
    
    IF (JK < KLEV) THEN
      IK = JK + 1
      IF (LDCUM(JLON)) THEN
        ZDUDT(JLON, JK) = ZDP(JLON, JK)*(ZMFUU(JLON, IK) - ZMFUU(JLON, JK) + ZMFDU(JLON, IK) - ZMFDU(JLON, JK))
        ZDVDT(JLON, JK) = ZDP(JLON, JK)*(ZMFUV(JLON, IK) - ZMFUV(JLON, JK) + ZMFDV(JLON, IK) - ZMFDV(JLON, JK))
      END IF
      
    ELSE
      IF (LDCUM(JLON)) THEN
        ZDUDT(JLON, JK) = -ZDP(JLON, JK)*(ZMFUU(JLON, JK) + ZMFDU(JLON, JK))
        ZDVDT(JLON, JK) = -ZDP(JLON, JK)*(ZMFUV(JLON, JK) + ZMFDV(JLON, JK))
      END IF
    END IF
    
  END DO
  
  IF (YDECUMF%RMFSOLUV == 0.0_JPRB) THEN
    
    !*    1.3          UPDATE TENDENCIES
    !                  -----------------
    
    DO JK=KTOPM2,KLEV
      IF (LDCUM(JLON)) THEN
        PTENU(JLON, JK) = PTENU(JLON, JK) + ZDUDT(JLON, JK)
        PTENV(JLON, JK) = PTENV(JLON, JK) + ZDVDT(JLON, JK)
      END IF
    END DO
    
  ELSE
    !----------------------------------------------------------------------
    
    !*      1.6          IMPLICIT SOLUTION
    !                    -----------------
    
    ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
    ! reuse ZMFUU=A and ZB=B;
    ! ZDUDT and ZDVDT correspond to the RHS ("constants") of the equation
    ! The solution is in ZR1 and ZR2
    
    LLCUMBAS(JLON, :) = .false.
    ZB(JLON, :) = 1.0_JPRB
    ZMFUU(JLON, :) = 0.0_JPRB
    
    ! Fill vectors A, B and RHS
    
    DO JK=KTOPM2,KLEV
      IK = JK + 1
      IM = JK - 1
      LLCUMBAS(JLON, JK) = LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1
      IF (LLCUMBAS(JLON, JK)) THEN
        ZZP = YDECUMF%RMFSOLUV*ZDP(JLON, JK)*PTSPHY
        ZMFUU(JLON, JK) = -ZZP*(PMFU(JLON, JK) + PMFD(JLON, JK))
        IF (JK < KLEV) THEN
          ZB(JLON, JK) = 1.0_JPRB + ZZP*(PMFU(JLON, IK) + PMFD(JLON, IK))
        ELSE
          ZB(JLON, JK) = 1.0_JPRB
        END IF
        ZZP = YDCST%RG*(PMFUO(JLON, JK) + YDECUMF%RMFADVWDD*PMFDO(JLON, JK)) / (PAP(JLON, JK) - PAP(JLON, IM))*PTSPHY*ZADVW
        ZU = ZZP*(PUEN(JLON, IM) - PUEN(JLON, JK))
        ZV = ZZP*(PVEN(JLON, IM) - PVEN(JLON, JK))
        ZDUDT(JLON, JK) = (ZDUDT(JLON, JK) + PTENU(JLON, JK)*YDECUMF%RMFSOLRHS)*PTSPHY + PUEN(JLON, JK) - ZU
        ZDVDT(JLON, JK) = (ZDVDT(JLON, JK) + PTENV(JLON, JK)*YDECUMF%RMFSOLRHS)*PTSPHY + PVEN(JLON, JK) - ZV
      END IF
    END DO
    
    CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFUU, ZB, ZDUDT, ZR1, YDSTACK=YLSTACK)
    
    CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFUU, ZB, ZDVDT, ZR2, YDSTACK=YLSTACK)
    
    ! prepare SPP perturbations
    LLPERT_CUDUDV = .false.
    LLPERT_CUDUDVS = .false.
    
    
    
    ! Compute tendencies
    
    DO JK=KTOPM2,KLEV
      IF (LLCUMBAS(JLON, JK)) THEN
        IF (LLPERT_CUDUDV .or. LLPPAR_CUDUV) THEN
          !Apply SPP perturbations or RP perturbation
          PTENU(JLON, JK) = PTENU(JLON, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + ZRDU*(ZR1(JLON, JK) - PUEN(JLON, JK))*ZTSPHY
          PTENV(JLON, JK) = PTENV(JLON, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + ZRDV*(ZR2(JLON, JK) - PVEN(JLON, JK))*ZTSPHY
        ELSE
          PTENU(JLON, JK) = PTENU(JLON, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + (ZR1(JLON, JK) - PUEN(JLON, JK))*ZTSPHY
          PTENV(JLON, JK) = PTENV(JLON, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + (ZR2(JLON, JK) - PVEN(JLON, JK))*ZTSPHY
        END IF
      END IF
    END DO
    
    
  END IF
  !----------------------------------------------------------------------
  
END SUBROUTINE CUDUDV_OPENACC
