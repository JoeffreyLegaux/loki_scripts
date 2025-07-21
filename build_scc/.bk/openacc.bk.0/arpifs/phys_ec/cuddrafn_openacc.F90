SUBROUTINE CUDDRAFN_OPENACC (YDCST, YDTHF, YDEPHLI, YDECUMF, KIDIA, KFDIA, KLON, KLEV, LDDRAF, PTENH, PQENH, PQSEN, PGEO, PGEOH,  &
& PAPH, PRFL, PTD, PQD, PMFU, PMFD, PMFDS, PMFDQ, PDMFDP, PDMFDE, PMFDDE_RATE, PKINED, YDSTACK)
  
  !          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
  
  !          PURPOSE.
  !          --------
  !          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
  !          (I.E. T,Q,U AND V AND FLUXES)
  
  !          INTERFACE
  !          ---------
  
  !          THIS ROUTINE IS CALLED FROM *CUMASTR*.
  !          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
  !          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
  !          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
  
  !          METHOD.
  !          --------
  !          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
  !          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
  !          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
  !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
  !    *PQSEN*        ENV. SATUR. SPEC. HUMIDITY (T+1)             KG/KG
  !    *PGEO*         GEOPOTENTIAL                                  M2/S2
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
  !    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
  !    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
  !    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
  !    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
  !    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
  !    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)
  !    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                   KG/(M2*S)
  !    *PKINED*       DOWNDRAFT KINETIC ENERGY                     M2/S2
  
  !          EXTERNALS
  !          ---------
  !          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
  !          SATURATED DESCENT
  
  !          REFERENCE
  !          ---------
  !          (TIEDTKE,1989)
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      03-08-28 : Clean-up detrainment rates   P. BECHTOLD
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUDDRAFN_OPENACC ) seq
  
  USE YOEPHLI, ONLY: TEPHLI
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOECUMF, ONLY: TECUMF
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  LOGICAL, INTENT(IN) :: LDDRAF(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(INOUT) :: PRFL(KLON)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFDS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFDQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFDP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFDE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFDDE_RATE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PKINED(KLON, KLEV)
  REAL(KIND=JPRB) :: ZDMFEN
  REAL(KIND=JPRB) :: ZDMFDE
  REAL(KIND=JPRB) :: ZCOND
  REAL(KIND=JPRB) :: ZOENTR
  REAL(KIND=JPRB) :: ZBUOY
  REAL(KIND=JPRB) :: ZPH
  LOGICAL :: LLO2
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IS
  INTEGER(KIND=JPIM) :: ITOPDE
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZBUO
  REAL(KIND=JPRB) :: ZBUOYZ
  REAL(KIND=JPRB) :: ZBUOYV
  REAL(KIND=JPRB) :: ZDMFDP
  REAL(KIND=JPRB) :: ZDZ
  REAL(KIND=JPRB) :: ZENTR
  REAL(KIND=JPRB) :: ZMFDQK
  REAL(KIND=JPRB) :: ZMFDSK
  REAL(KIND=JPRB) :: ZQDDE
  REAL(KIND=JPRB) :: ZQEEN
  REAL(KIND=JPRB) :: ZRAIN
  REAL(KIND=JPRB) :: ZSDDE
  REAL(KIND=JPRB) :: ZSEEN
  REAL(KIND=JPRB) :: ZZENTR
  REAL(KIND=JPRB) :: ZRG
  REAL(KIND=JPRB) :: ZFACBUO
  REAL(KIND=JPRB) :: Z_CWDRAG
  REAL(KIND=JPRB) :: ZDKBUO
  REAL(KIND=JPRB) :: ZDKEN
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cuadjtq.intfb.h"
#include "abor1.intfb.h"
#include "fcttre.func.h"
#include "cuadjtq.func.h"
  INTEGER(KIND=JPIM) :: CUADJTQ_JL
  REAL(KIND=JPRB) :: Z1S
  REAL(KIND=JPRB) :: Z2S
  REAL(KIND=JPRB) :: CUADJTQ_ZCOND
  REAL(KIND=JPRB) :: ZCOND1
  REAL(KIND=JPRB) :: ZCOR
  REAL(KIND=JPRB) :: ZFOEEWI
  REAL(KIND=JPRB) :: ZFOEEWL
  REAL(KIND=JPRB) :: ZOEALFA
  REAL(KIND=JPRB) :: ZQMAX
  REAL(KIND=JPRB) :: ZQSAT
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZQP
  REAL(KIND=JPRB) :: ZL
  REAL(KIND=JPRB) :: ZI
  REAL(KIND=JPRB) :: ZF
  LOGICAL :: LLFLAG
  REAL(KIND=JPHOOK) :: CUADJTQ_ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  ITOPDE = YDECUMF%NJKT3
  ZRG = 1.0_JPRB / YDCST%RG
  ZFACBUO = 0.5_JPRB / (1.0_JPRB + 0.5_JPRB)
  Z_CWDRAG = (3._JPRB / 8._JPRB)*0.506_JPRB / 0.2_JPRB
  !----------------------------------------------------------------------
  
  !     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
  !                     (A) CALCULATING ENTRAINMENT/DETRAINMENT RATES,
  !                         INCLUDING ORGANIZED ENTRAINMENT DEPENDENT ON
  !                         NEGATIVE BUOYANCY AND ASSUMING
  !                         LINEAR DECREASE OF MASSFLUX IN PBL
  !                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
  !                         AND MOISTENING IS CALCULATED IN *CUADJTQ*
  !                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
  !                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
  !                    -------------------------------------------------
  
  ZOENTR = 0.0_JPRB
  ZBUOY = 0.0_JPRB
  ZDMFEN = 0.0_JPRB
  ZDMFDE = 0.0_JPRB
  PDMFDE(JLON, :) = 0.0_JPRB
  PMFDDE_RATE(JLON, :) = 0.0_JPRB
  PKINED(JLON, :) = 0.0_JPRB
  
  DO JK=3,KLEV
    IS = 0
    ZPH = PAPH(JLON, JK)
    LLO2 = LDDRAF(JLON) .and. PMFD(JLON, JK - 1) < 0.0_JPRB
    IF (LLO2) THEN
      IS = IS + 1
    END IF
    IF (IS == 0) CYCLE
    
    IF (LLO2) THEN
      ZENTR = YDECUMF%ENTRDD*PMFD(JLON, JK - 1)*(PGEOH(JLON, JK - 1) - PGEOH(JLON, JK))*ZRG
      ! &*(1.6_JPRB-MIN(1.0_JPRB,PQENH(JL,JK-1)/PQSEN(JL,JK)))
      ZDMFEN = ZENTR
      ZDMFDE = ZENTR
    END IF
    
    IF (JK > ITOPDE) THEN
      IF (LLO2) THEN
        ZDMFEN = 0.0_JPRB
        ZDMFDE = PMFD(JLON, ITOPDE)*(PAPH(JLON, JK) - PAPH(JLON, JK - 1)) / (PAPH(JLON, KLEV + 1) - PAPH(JLON, ITOPDE))
      END IF
    END IF
    
    IF (JK <= ITOPDE) THEN
      IF (LLO2) THEN
        ZDZ = -(PGEOH(JLON, JK - 1) - PGEOH(JLON, JK))*ZRG
        ZZENTR = ZOENTR*ZDZ*PMFD(JLON, JK - 1)
        ZDMFEN = ZDMFEN + ZZENTR
        ZDMFEN = MAX(ZDMFEN, 0.3_JPRB*PMFD(JLON, JK - 1))
        ZDMFEN = MAX(ZDMFEN, -0.75_JPRB*PMFU(JLON, JK) - (PMFD(JLON, JK - 1) - ZDMFDE))
        ZDMFEN = MIN(ZDMFEN, 0.0_JPRB)
      END IF
      
      PDMFDE(JLON, JK) = ZDMFEN - ZDMFDE
      
    END IF
    IF (LLO2) THEN
      PMFD(JLON, JK) = PMFD(JLON, JK - 1) + ZDMFEN - ZDMFDE
      ZSEEN = (YDCST%RCPD*PTENH(JLON, JK - 1) + PGEOH(JLON, JK - 1))*ZDMFEN
      ZQEEN = PQENH(JLON, JK - 1)*ZDMFEN
      ZSDDE = (YDCST%RCPD*PTD(JLON, JK - 1) + PGEOH(JLON, JK - 1))*ZDMFDE
      ZQDDE = PQD(JLON, JK - 1)*ZDMFDE
      ZMFDSK = PMFDS(JLON, JK - 1) + ZSEEN - ZSDDE
      ZMFDQK = PMFDQ(JLON, JK - 1) + ZQEEN - ZQDDE
      PQD(JLON, JK) = ZMFDQK*(1.0_JPRB / MIN(-YDECUMF%RMFCMIN, PMFD(JLON, JK)))
      PTD(JLON, JK) = (ZMFDSK*(1.0_JPRB / MIN(-YDECUMF%RMFCMIN, PMFD(JLON, JK))) - PGEOH(JLON, JK)) / YDCST%RCPD
      PTD(JLON, JK) = MIN(400._JPRB, PTD(JLON, JK))
      PTD(JLON, JK) = MAX(100._JPRB, PTD(JLON, JK))
      ZCOND = PQD(JLON, JK)
    END IF
    
    IK = JK
    ! [Loki] inlined child subroutine: CUADJTQ
    ! =========================================
    
    !----------------------------------------------------------------------
    
    !     1.           DEFINE CONSTANTS
    !                  ----------------
    
    
    !IF (LHOOK) CALL DR_HOOK('CUADJTQ',0,ZHOOK_HANDLE)
    
    
    ZQMAX = 0.5_JPRB
    
    !*********************************************
    IF (.not.YDEPHLI%LPHYLIN) THEN
      !*********************************************
      
      !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
      !                  -----------------------------------------------------
      
      
      
      !DIR$    IVDEP
      !OCL NOVREC
      IF (LLO2) THEN
        ZQP = 1.0_JPRB / ZPH
        ZQSAT = FOEEWMCU(PTD(JLON, IK))*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        CUADJTQ_ZCOND = (PQD(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEMCU(PTD(JLON, IK)))
        CUADJTQ_ZCOND = MIN(CUADJTQ_ZCOND, 0.0_JPRB)
        PTD(JLON, IK) = PTD(JLON, IK) + FOELDCPMCU(PTD(JLON, IK))*CUADJTQ_ZCOND
        PQD(JLON, IK) = PQD(JLON, IK) - CUADJTQ_ZCOND
        ZQSAT = FOEEWMCU(PTD(JLON, IK))*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        ZCOND1 = (PQD(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEMCU(PTD(JLON, IK)))
        IF (CUADJTQ_ZCOND == 0.0_JPRB) ZCOND1 = MIN(ZCOND1, 0.0_JPRB)
        PTD(JLON, IK) = PTD(JLON, IK) + FOELDCPMCU(PTD(JLON, IK))*ZCOND1
        PQD(JLON, IK) = PQD(JLON, IK) - ZCOND1
      END IF
      
      
      
      
      
      !*********************************************
    ELSE
      !*********************************************
      
      !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
      !                  -----------------------------------------------------
      
      
      
      !DIR$    IVDEP
      !OCL NOVREC
      IF (LLO2) THEN
        ZQP = 1.0_JPRB / ZPH
        
        ZTARG = PTD(JLON, IK)
        ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
        ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
        ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
        ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
        Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
        ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
        
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        
        Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB / (ZTARG &
        &  - YDTHF%R4IES)**2)
        CUADJTQ_ZCOND = (PQD(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
        
        CUADJTQ_ZCOND = MIN(CUADJTQ_ZCOND, 0.0_JPRB)
        
        IF (CUADJTQ_ZCOND /= 0.0_JPRB) THEN
          
          PTD(JLON, IK) = PTD(JLON, IK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*CUADJTQ_ZCOND
          PQD(JLON, IK) = PQD(JLON, IK) - CUADJTQ_ZCOND
          ZTARG = PTD(JLON, IK)
          ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
          ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
          ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
          ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
          Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
          ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
          
          ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
          ZQSAT = ZQSAT*ZCOR
          
          Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB /  &
          & (ZTARG - YDTHF%R4IES)**2)
          ZCOND1 = (PQD(JLON, IK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
          
          PTD(JLON, IK) = PTD(JLON, IK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
          
          PQD(JLON, IK) = PQD(JLON, IK) - ZCOND1
        END IF
      END IF
      
      
      
      
      !*********************************************
    END IF
    !*********************************************
    
    
    !IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
    
    ! =========================================
    
    IF (LLO2) THEN
      ZCOND = ZCOND - PQD(JLON, JK)
      ZBUO = PTD(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQD(JLON, JK)) - PTENH(JLON, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JLON, JK))
      IF (PRFL(JLON) > 0.0_JPRB .and. PMFU(JLON, JK) > 0.0_JPRB) THEN
        ZRAIN = PRFL(JLON) / PMFU(JLON, JK)
        ZBUO = ZBUO - PTD(JLON, JK)*ZRAIN
      END IF
      IF (ZBUO >= 0.0_JPRB .or. PRFL(JLON) <= (PMFD(JLON, JK)*ZCOND)) THEN
        PMFD(JLON, JK) = 0.0_JPRB
        ZBUO = 0.0_JPRB
      END IF
      PMFDS(JLON, JK) = (YDCST%RCPD*PTD(JLON, JK) + PGEOH(JLON, JK))*PMFD(JLON, JK)
      PMFDQ(JLON, JK) = PQD(JLON, JK)*PMFD(JLON, JK)
      ZDMFDP = -PMFD(JLON, JK)*ZCOND
      PDMFDP(JLON, JK - 1) = ZDMFDP
      PRFL(JLON) = PRFL(JLON) + ZDMFDP
      
      ! COMPUTE ORGANIZED ENTRAINMENT FOR USE AT NEXT LEVEL
      
      ZBUOYZ = ZBUO / PTENH(JLON, JK)
      ZBUOYV = ZBUOYZ
      ZBUOYZ = MIN(ZBUOYZ, 0.0_JPRB)
      ZDZ = -(PGEO(JLON, JK - 1) - PGEO(JLON, JK))
      ZBUOY = ZBUOY + ZBUOYZ*ZDZ
      ZOENTR = YDCST%RG*ZBUOYZ*0.5_JPRB / (1.0_JPRB + ZBUOY)
      
      ! STORE DOWNDRAUGHT DETRAINMENT RATES
      
      PMFDDE_RATE(JLON, JK) = -ZDMFDE
      
      ! COMPUTE KINETIC ENERGY
      
      ZDKBUO = ZDZ*ZBUOYV*ZFACBUO
      IF (ZDMFEN < 0.0_JPRB) THEN
        ZDKEN = MIN(1.0_JPRB, (1 + Z_CWDRAG)*ZDMFEN / MIN(-YDECUMF%RMFCMIN, PMFD(JLON, JK - 1)))
      ELSE
        ZDKEN = MIN(1.0_JPRB, (1 + Z_CWDRAG)*ZDMFDE / MIN(-YDECUMF%RMFCMIN, PMFD(JLON, JK - 1)))
      END IF
      PKINED(JLON, JK) = MAX(0.0_JPRB, (PKINED(JLON, JK - 1)*(1 - ZDKEN) + ZDKBUO) / (1 + ZDKEN))
      
    END IF
    
  END DO
  
  
  
  ! (C) Copyright 1988- ECMWF.
  !
  ! This software is licensed under the terms of the Apache Licence Version 2.0
  ! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
  !
  ! In applying this licence, ECMWF does not waive the privileges and immunities
  ! granted to it by virtue of its status as an intergovernmental organisation
  ! nor does it submit to any jurisdiction.
  
  
END SUBROUTINE CUDDRAFN_OPENACC
