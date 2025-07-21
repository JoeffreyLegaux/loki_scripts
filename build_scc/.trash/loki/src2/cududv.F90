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
  INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP
  LOGICAL, INTENT(IN) :: LDCUM
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
  
#include "cubidiag.intfb.h"
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
  JL = KIDIA
  !----------------------------------------------------------------------
  
  ZIMP = 1.0_JPRB - YDECUMF%RMFSOLUV
  ZTSPHY = 1.0_JPRB / PTSPHY
  
  LLPPAR_CUDUV = .false.
  LLPERT_CUDUDV = .false.
  LLPERT_CUDUDVS = .false.
  
  ZADVW = 0.0_JPRB
  IF (KTYPE == 1) ZADVW = YDECUMF%RMFADVW
  
  
  DO JK=1,KLEV
    IF (LDCUM) THEN
      ZDP(JL, JK) = YDCST%RG / (PAPH(JL, JK + 1) - PAPH(JL, JK))
    END IF
  END DO
  
  !----------------------------------------------------------------------
  
  
  !*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
  !                  ----------------------------------------------
  
  DO JK=KTOPM2,KLEV
    IK = JK - 1
    IF (LDCUM) THEN
      ZMFUU(JL, JK) = PMFU(JL, JK)*(PUU(JL, JK) - ZIMP*PUEN(JL, IK))
      ZMFUV(JL, JK) = PMFU(JL, JK)*(PVU(JL, JK) - ZIMP*PVEN(JL, IK))
      ZMFDU(JL, JK) = PMFD(JL, JK)*(PUD(JL, JK) - ZIMP*PUEN(JL, IK))
      ZMFDV(JL, JK) = PMFD(JL, JK)*(PVD(JL, JK) - ZIMP*PVEN(JL, IK))
    END IF
  END DO
  
  ! linear fluxes below cloud
  IF (YDECUMF%RMFSOLUV == 0.0_JPRB) THEN
    DO JK=KTOPM2,KLEV
      !DIR$ IVDEP
      !OCL NOVREC
      IF (LDCUM .and. JK > KCBOT) THEN
        IKB = KCBOT
        ZZP = ((PAPH(JL, KLEV + 1) - PAPH(JL, JK)) / (PAPH(JL, KLEV + 1) - PAPH(JL, IKB)))
        IF (KTYPE == 3) THEN
          ZZP = ZZP*ZZP
        END IF
        ZMFUU(JL, JK) = ZMFUU(JL, IKB)*ZZP
        ZMFUV(JL, JK) = ZMFUV(JL, IKB)*ZZP
        ZMFDU(JL, JK) = ZMFDU(JL, IKB)*ZZP
        ZMFDV(JL, JK) = ZMFDV(JL, IKB)*ZZP
      END IF
    END DO
  END IF
  
  !*    1.2          COMPUTE TENDENCIES
  !                  ------------------
  
  DO JK=KTOPM2,KLEV
    
    IF (JK < KLEV) THEN
      IK = JK + 1
      IF (LDCUM) THEN
        ZDUDT(JL, JK) = ZDP(JL, JK)*(ZMFUU(JL, IK) - ZMFUU(JL, JK) + ZMFDU(JL, IK) - ZMFDU(JL, JK))
        ZDVDT(JL, JK) = ZDP(JL, JK)*(ZMFUV(JL, IK) - ZMFUV(JL, JK) + ZMFDV(JL, IK) - ZMFDV(JL, JK))
      END IF
      
    ELSE
      IF (LDCUM) THEN
        ZDUDT(JL, JK) = -ZDP(JL, JK)*(ZMFUU(JL, JK) + ZMFDU(JL, JK))
        ZDVDT(JL, JK) = -ZDP(JL, JK)*(ZMFUV(JL, JK) + ZMFDV(JL, JK))
      END IF
    END IF
    
  END DO
  
  IF (YDECUMF%RMFSOLUV == 0.0_JPRB) THEN
    
    !*    1.3          UPDATE TENDENCIES
    !                  -----------------
    
    DO JK=KTOPM2,KLEV
      IF (LDCUM) THEN
        PTENU(JL, JK) = PTENU(JL, JK) + ZDUDT(JL, JK)
        PTENV(JL, JK) = PTENV(JL, JK) + ZDVDT(JL, JK)
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
    
    LLCUMBAS(JL, :) = .false.
    ZB(JL, :) = 1.0_JPRB
    ZMFUU(JL, :) = 0.0_JPRB
    
    ! Fill vectors A, B and RHS
    
    DO JK=KTOPM2,KLEV
      IK = JK + 1
      IM = JK - 1
      LLCUMBAS(JL, JK) = LDCUM .and. JK >= KCTOP - 1
      IF (LLCUMBAS(JL, JK)) THEN
        ZZP = YDECUMF%RMFSOLUV*ZDP(JL, JK)*PTSPHY
        ZMFUU(JL, JK) = -ZZP*(PMFU(JL, JK) + PMFD(JL, JK))
        IF (JK < KLEV) THEN
          ZB(JL, JK) = 1.0_JPRB + ZZP*(PMFU(JL, IK) + PMFD(JL, IK))
        ELSE
          ZB(JL, JK) = 1.0_JPRB
        END IF
        ZZP = YDCST%RG*(PMFUO(JL, JK) + YDECUMF%RMFADVWDD*PMFDO(JL, JK)) / (PAP(JL, JK) - PAP(JL, IM))*PTSPHY*ZADVW
        ZU = ZZP*(PUEN(JL, IM) - PUEN(JL, JK))
        ZV = ZZP*(PVEN(JL, IM) - PVEN(JL, JK))
        ZDUDT(JL, JK) = (ZDUDT(JL, JK) + PTENU(JL, JK)*YDECUMF%RMFSOLRHS)*PTSPHY + PUEN(JL, JK) - ZU
        ZDVDT(JL, JK) = (ZDVDT(JL, JK) + PTENV(JL, JK)*YDECUMF%RMFSOLRHS)*PTSPHY + PVEN(JL, JK) - ZV
      END IF
    END DO
    
    CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFUU, ZB, ZDUDT, ZR1, YDSTACK=YLSTACK)
    
    CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFUU, ZB, ZDVDT, ZR2, YDSTACK=YLSTACK)
    
    ! prepare SPP perturbations
    IF (YDSPP_CONFIG%LSPP) THEN
      IPN1 = YDSPP_CONFIG%PPTR%CUDUDV
      LLPERT_CUDUDV = IPN1 > 0
      IF (LLPERT_CUDUDV .or. LLPPAR_CUDUV) THEN
        PN1 = YDSPP_CONFIG%SM%PN(IPN1)
        IPCUDU = PN1%MP
        IF (PN1%NRF == 1) THEN
          IPCUDV = IPCUDU
        ELSE
          IPCUDV = IPCUDU + 1
        END IF
      END IF
      IPN2 = YDSPP_CONFIG%PPTR%CUDUDVS
      LLPERT_CUDUDVS = IPN2 > 0
      IF (LLPERT_CUDUDVS) THEN
        PN2 = YDSPP_CONFIG%SM%PN(IPN2)
        IPCUDUS = PN2%MP
        IF (PN2%NRF == 1) THEN
          IPCUDVS = IPCUDUS
        ELSE
          IPCUDVS = IPCUDUS + 1
        END IF
      END IF
    ELSE
      LLPERT_CUDUDV = .false.
      LLPERT_CUDUDVS = .false.
    END IF
    
    IF (YDSPP_CONFIG%LSPP) THEN
      IF (LLPERT_CUDUDV .or. LLPERT_CUDUDVS) THEN
        ZMDU = 1.0_JPRB          ! mean      zonal momentum transport pdf
        ZMDV = 1.0_JPRB          ! mean meridional momentum transport pdf
        ZLIMN2 = 9.0_JPRB*PN1%SDEV**2          ! limit for norm squared (3 stdev)
        
        IF (KTYPE == 1 .and. LLPERT_CUDUDV .or. KTYPE == 2 .and. LLPERT_CUDUDVS) THEN
          !perturb deep/shallow convection
          
          IF (KTYPE == 1) THEN
            IF (.not.LLPPAR_CUDUV) THEN
              ZXU = PGP2DSPP(JL, IPCUDU)
              ZXV = PGP2DSPP(JL, IPCUDV)
              ZXU_MAG = PN1%XMAG(1)
              ZXV_MAG = PN1%XMAG(2)
            ELSE
              ZXU = CUDU_MOD
              ZXV = CUDV_MOD
            END IF
          END IF
          
          IF (KTYPE == 2) THEN
            ZXU = PGP2DSPP(JL, IPCUDUS)
            ZXV = PGP2DSPP(JL, IPCUDVS)
            ZXU_MAG = PN2%XMAG(1)
            ZXV_MAG = PN2%XMAG(2)
          END IF
          
          ZN2 = ZXU**2 + ZXV**2
          IF (ZN2 > ZLIMN2) THEN
            ZFAC = SQRT(ZLIMN2 / ZN2)
            ZXU = ZFAC*ZXU
            ZXV = ZFAC*ZXV
          END IF
          IF (.not.LLPPAR_CUDUV) THEN
            ZRDU = ZMDU + ZXU_MAG*ZXU
            ZRDV = ZMDV + ZXV_MAG*ZXV
          ELSE
            ZRDU = ZMDU + ZXU
            ZRDV = ZMDV + ZXV
          END IF
        ELSE
          !other convection - unperturbed
          ZRDU = ZMDU
          ZRDV = ZMDV
        END IF
      END IF
    END IF
    
    
    ! Compute tendencies
    
    DO JK=KTOPM2,KLEV
      IF (LLCUMBAS(JL, JK)) THEN
        IF (LLPERT_CUDUDV .or. LLPPAR_CUDUV) THEN
          !Apply SPP perturbations or RP perturbation
          PTENU(JL, JK) = PTENU(JL, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + ZRDU*(ZR1(JL, JK) - PUEN(JL, JK))*ZTSPHY
          PTENV(JL, JK) = PTENV(JL, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + ZRDV*(ZR2(JL, JK) - PVEN(JL, JK))*ZTSPHY
        ELSE
          PTENU(JL, JK) = PTENU(JL, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + (ZR1(JL, JK) - PUEN(JL, JK))*ZTSPHY
          PTENV(JL, JK) = PTENV(JL, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + (ZR2(JL, JK) - PVEN(JL, JK))*ZTSPHY
        END IF
      END IF
    END DO
    
    
  END IF
  !----------------------------------------------------------------------
  
END SUBROUTINE CUDUDV_OPENACC
