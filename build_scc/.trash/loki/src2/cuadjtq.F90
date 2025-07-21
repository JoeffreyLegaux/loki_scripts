SUBROUTINE CUADJTQ_OPENACC (YDTHF, YDCST, YDEPHLI, KIDIA, KFDIA, KLON, KLEV, KK, PSP, PT, PQ, LDFLAG, KCALL, LDOFLAG, YDSTACK)
  
  !          PURPOSE.
  !          --------
  !          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
  
  !          INTERFACE
  !          ---------
  !          THIS ROUTINE IS CALLED FROM SUBROUTINES:
  !              *COND*     (T AND Q AT CONDENSATION LEVEL)
  !              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
  !              *CUASC*    (T AND Q AT CLOUD LEVELS)
  !              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
  !              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
  !          INPUT ARE UNADJUSTED T AND Q VALUES,
  !          IT RETURNS ADJUSTED VALUES OF T AND Q
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KK*           LEVEL
  !    *KCALL*        DEFINES CALCULATION AS
  !                      KCALL=0  ENV. T AND QS IN*CUINI*
  !                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
  !                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
  
  !     INPUT PARAMETERS (LOGICAL):
  
  !    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
  
  !     INPUT PARAMETERS (REAL):
  
  !    *PSP*          PRESSURE                                        PA
  
  !     UPDATED PARAMETERS (REAL):
  
  !    *PT*           TEMPERATURE                                     K
  !    *PQ*           SPECIFIC HUMIDITY                             KG/KG
  
  !          EXTERNALS
  !          ---------
  !          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
  !          FOR CONDENSATION CALCULATIONS.
  !          THE TABLES ARE INITIALISED IN *SUPHEC*.
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.     12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      J.HAGUE               03-01-13   MASS Vector Functions
  !      J.HAGUE               03-07-07   More MASS V.F.
  !      M.Hamrud              01-Oct-2003 CY28 Cleaning
  !      J.Hague & D.Salmond   22-Nov-2005 Optimisations
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUADJTQ_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOEPHLI, ONLY: TEPHLI
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KK
  REAL(KIND=JPRB), INTENT(IN) :: PSP
  REAL(KIND=JPRB), INTENT(INOUT) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PQ(KLON, KLEV)
  LOGICAL, INTENT(IN) :: LDFLAG
  INTEGER(KIND=JPIM), INTENT(IN) :: KCALL
  LOGICAL, OPTIONAL, INTENT(IN) :: LDOFLAG
  
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: Z1S
  REAL(KIND=JPRB) :: Z2S
  REAL(KIND=JPRB) :: ZCOND
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
  
#include "abor1.intfb.h"
  
  !DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "cuadjtq.func.h"
  
  !     STATEMENT FUNCTIONS
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JL = KIDIA
  
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
    
    IF (KCALL == 1 .or. KCALL == 6) THEN
      
      !   mixed phase saturation
      
      LLFLAG = LDFLAG
      IF (KCALL == 6) THEN
        IF (PRESENT(LDOFLAG)) THEN
          LLFLAG = LLFLAG .and. .not.LDOFLAG
        ELSE
          CALL ABOR1_OPENACC('CUADJTQ: LDOFLAG has to be present when KCALL==6', YDSTACK=YLSTACK)
        END IF
      END IF
      
      IF (LLFLAG) THEN
        ZQP = 1.0_JPRB / PSP
        ZL = 1.0_JPRB / (PT(JL, KK) - YDTHF%R4LES)
        ZI = 1.0_JPRB / (PT(JL, KK) - YDTHF%R4IES)
        !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
        ZQSAT = YDTHF%R2ES*(FOEALFCU(PT(JL, KK))*EXP(YDTHF%R3LES*(PT(JL, KK) - YDCST%RTT)*ZL) + (1.0_JPRB - FOEALFCU(PT(JL, KK))) &
        & *EXP(YDTHF%R3IES*(PT(JL, KK) - YDCST%RTT)*ZI))
        ZQSAT = ZQSAT*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
        ZF = FOEALFCU(PT(JL, KK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(PT(JL, KK)))*YDTHF%R5ALSCP*ZI**2
        ZCOND = (PQ(JL, KK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
        ZCOND = MAX(ZCOND, 0.0_JPRB)
        PT(JL, KK) = PT(JL, KK) + FOELDCPMCU(PT(JL, KK))*ZCOND
        PQ(JL, KK) = PQ(JL, KK) - ZCOND
        ZL = 1.0_JPRB / (PT(JL, KK) - YDTHF%R4LES)
        ZI = 1.0_JPRB / (PT(JL, KK) - YDTHF%R4IES)
        !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
        ZQSAT = YDTHF%R2ES*(FOEALFCU(PT(JL, KK))*EXP(YDTHF%R3LES*(PT(JL, KK) - YDCST%RTT)*ZL) + (1.0_JPRB - FOEALFCU(PT(JL, KK))) &
        & *EXP(YDTHF%R3IES*(PT(JL, KK) - YDCST%RTT)*ZI))
        ZQSAT = ZQSAT*ZQP
        ZQSAT = FMINJ(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
        ZF = FOEALFCU(PT(JL, KK))*YDTHF%R5ALVCP*ZL**2 + (1.0_JPRB - FOEALFCU(PT(JL, KK)))*YDTHF%R5ALSCP*ZI**2
        ZCOND1 = (PQ(JL, KK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
        IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
        PT(JL, KK) = PT(JL, KK) + FOELDCPMCU(PT(JL, KK))*ZCOND1
        PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      END IF
      
      IF (KCALL == 6) THEN
        
        IF (LDFLAG .and. LDOFLAG) THEN
          ZQP = 1.0_JPRB / PSP
          ZL = 1.0_JPRB / (PT(JL, KK) - YDTHF%R4LES)
          ZQSAT = YDTHF%R2ES*EXP(YDTHF%R3LES*(PT(JL, KK) - YDCST%RTT)*ZL)
          ZQSAT = ZQSAT*ZQP
          ZQSAT = MIN(0.5_JPRB, ZQSAT)
          ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
          ZF = YDTHF%R5ALVCP*ZL**2
          ZCOND = (PQ(JL, KK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
          ZCOND = MAX(ZCOND, 0.0_JPRB)
          PT(JL, KK) = PT(JL, KK) + YDTHF%RALVDCP*ZCOND
          PQ(JL, KK) = PQ(JL, KK) - ZCOND
          ZL = 1.0_JPRB / (PT(JL, KK) - YDTHF%R4LES)
          ZQSAT = YDTHF%R2ES*EXP(YDTHF%R3LES*(PT(JL, KK) - YDCST%RTT)*ZL)
          ZQSAT = ZQSAT*ZQP
          ZQSAT = FMINJ(0.5_JPRB, ZQSAT)
          ZCOR = 1.0_JPRB - YDCST%RETV*ZQSAT
          ZF = YDTHF%R5ALVCP*ZL**2
          ZCOND1 = (PQ(JL, KK)*ZCOR**2 - ZQSAT*ZCOR) / (ZCOR**2 + ZQSAT*ZF)
          IF (ZCOND == 0.0_JPRB) ZCOND1 = 0.0_JPRB
          PT(JL, KK) = PT(JL, KK) + YDTHF%RALVDCP*ZCOND1
          PQ(JL, KK) = PQ(JL, KK) - ZCOND1
        END IF
        
      END IF
      
    END IF
    
    IF (KCALL == 2) THEN
      
      !DIR$    IVDEP
      !OCL NOVREC
      IF (LDFLAG) THEN
        ZQP = 1.0_JPRB / PSP
        ZQSAT = FOEEWMCU(PT(JL, KK))*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        ZCOND = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEMCU(PT(JL, KK)))
        ZCOND = MIN(ZCOND, 0.0_JPRB)
        PT(JL, KK) = PT(JL, KK) + FOELDCPMCU(PT(JL, KK))*ZCOND
        PQ(JL, KK) = PQ(JL, KK) - ZCOND
        ZQSAT = FOEEWMCU(PT(JL, KK))*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEMCU(PT(JL, KK)))
        IF (ZCOND == 0.0_JPRB) ZCOND1 = MIN(ZCOND1, 0.0_JPRB)
        PT(JL, KK) = PT(JL, KK) + FOELDCPMCU(PT(JL, KK))*ZCOND1
        PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      END IF
      
    END IF
    
    IF (KCALL == 0) THEN
      
      !DIR$    IVDEP
      !OCL NOVREC
      
      !DIR$ LOOP_INFO EST_TRIPS(16)
      ZQP = 1.0_JPRB / PSP
      ZQSAT = FOEEWM(PT(JL, KK))*ZQP
      ZQSAT = MIN(0.5_JPRB, ZQSAT)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(PT(JL, KK)))
      PT(JL, KK) = PT(JL, KK) + FOELDCPM(PT(JL, KK))*ZCOND1
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      ZQSAT = FOEEWM(PT(JL, KK))*ZQP
      ZQSAT = MIN(0.5_JPRB, ZQSAT)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(PT(JL, KK)))
      PT(JL, KK) = PT(JL, KK) + FOELDCPM(PT(JL, KK))*ZCOND1
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      
    END IF
    
    IF (KCALL == 4) THEN
      
      IF (LDFLAG) THEN
        ZQP = 1.0_JPRB / PSP
        ZQSAT = FOEEWM(PT(JL, KK))*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        ZCOND = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(PT(JL, KK)))
        PT(JL, KK) = PT(JL, KK) + FOELDCPM(PT(JL, KK))*ZCOND
        PQ(JL, KK) = PQ(JL, KK) - ZCOND
        ZQSAT = FOEEWM(PT(JL, KK))*ZQP
        ZQSAT = MIN(0.5_JPRB, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(PT(JL, KK)))
        PT(JL, KK) = PT(JL, KK) + FOELDCPM(PT(JL, KK))*ZCOND1
        PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      END IF
      
    END IF
    
    IF (KCALL == 5) THEN
      ! Same as 4 but with LDFLAG all true
      
      !OCL NOVREC
      !DIR$    IVDEP
      !DIR$ LOOP_INFO EST_TRIPS(16)
      ZQP = 1.0_JPRB / PSP
      ZQSAT = FOEEWM(PT(JL, KK))*ZQP
      ZQSAT = MIN(0.5_JPRB, ZQSAT)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      ZCOND = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(PT(JL, KK)))
      PT(JL, KK) = PT(JL, KK) + FOELDCPM(PT(JL, KK))*ZCOND
      PQ(JL, KK) = PQ(JL, KK) - ZCOND
      ZQSAT = FOEEWM(PT(JL, KK))*ZQP
      ZQSAT = MIN(0.5_JPRB, ZQSAT)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(PT(JL, KK)))
      PT(JL, KK) = PT(JL, KK) + FOELDCPM(PT(JL, KK))*ZCOND1
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
    END IF
    
    IF (KCALL == 3) THEN
      !DIR$ LOOP_INFO EST_TRIPS(16)
      ZQP = 1.0_JPRB / PSP
      ZQSAT = FOEEWMCU(PT(JL, KK))*ZQP
      ZQSAT = MIN(0.5_JPRB, ZQSAT)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEMCU(PT(JL, KK)))
      PT(JL, KK) = PT(JL, KK) + FOELDCPMCU(PT(JL, KK))*ZCOND1
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      ZQSAT = FOEEWMCU(PT(JL, KK))*ZQP
      ZQSAT = MIN(0.5_JPRB, ZQSAT)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEMCU(PT(JL, KK)))
      PT(JL, KK) = PT(JL, KK) + FOELDCPMCU(PT(JL, KK))*ZCOND1
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      
    END IF
    !*********************************************
  ELSE
    !*********************************************
    
    !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
    !                  -----------------------------------------------------
    
    IF (KCALL == 1) THEN
      
      !DIR$    IVDEP
      !OCL NOVREC
      !DIR$ LOOP_INFO EST_TRIPS(16)
      IF (LDFLAG) THEN
        ZQP = 1.0_JPRB / PSP
        ZTARG = PT(JL, KK)
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
        ZCOND = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
        
        ZCOND = MAX(ZCOND, 0.0_JPRB)
        
        IF (ZCOND /= 0.0_JPRB) THEN
          
          PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND
          PQ(JL, KK) = PQ(JL, KK) - ZCOND
          ZTARG = PT(JL, KK)
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
          ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
          
          PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
          
          PQ(JL, KK) = PQ(JL, KK) - ZCOND1
        END IF
      END IF
      
    END IF
    
    IF (KCALL == 2) THEN
      
      !DIR$    IVDEP
      !OCL NOVREC
      IF (LDFLAG) THEN
        ZQP = 1.0_JPRB / PSP
        
        ZTARG = PT(JL, KK)
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
        ZCOND = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
        
        ZCOND = MIN(ZCOND, 0.0_JPRB)
        
        IF (ZCOND /= 0.0_JPRB) THEN
          
          PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND
          PQ(JL, KK) = PQ(JL, KK) - ZCOND
          ZTARG = PT(JL, KK)
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
          ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
          
          PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
          
          PQ(JL, KK) = PQ(JL, KK) - ZCOND1
        END IF
      END IF
      
    END IF
    
    IF (KCALL == 0) THEN
      
      !DIR$    IVDEP
      !OCL NOVREC
      ZQP = 1.0_JPRB / PSP
      
      ZTARG = PT(JL, KK)
      ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
      ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
      ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
      ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
      Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
      ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
      
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      
      Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB / (ZTARG - &
      &  YDTHF%R4IES)**2)
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
      
      PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
      
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      
      ZTARG = PT(JL, KK)
      ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
      ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
      ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
      ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
      Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
      ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
      
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
      ZQSAT = ZQSAT*ZCOR
      
      Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB / (ZTARG - &
      &  YDTHF%R4IES)**2)
      ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
      
      PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
      
      PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      
    END IF
    
    IF (KCALL == 4) THEN
      
      !DIR$    IVDEP
      !OCL NOVREC
      IF (LDFLAG) THEN
        ZQP = 1.0_JPRB / PSP
        
        ZTARG = PT(JL, KK)
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
        ZCOND = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
        
        PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND
        
        PQ(JL, KK) = PQ(JL, KK) - ZCOND
        
        ZTARG = PT(JL, KK)
        ZOEALFA = 0.5_JPRB*(TANH(YDEPHLI%RLPAL1*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB)
        ZFOEEWL = YDTHF%R2ES*EXP(YDTHF%R3LES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4LES))
        ZFOEEWI = YDTHF%R2ES*EXP(YDTHF%R3IES*(ZTARG - YDCST%RTT) / (ZTARG - YDTHF%R4IES))
        ZQSAT = ZQP*(ZOEALFA*ZFOEEWL + (1.0_JPRB - ZOEALFA)*ZFOEEWI)
        Z1S = TANH(YDEPHLI%RLPAL2*(ZQSAT - ZQMAX))
        ZQSAT = 0.5_JPRB*((1.0_JPRB - Z1S)*ZQSAT + (1.0_JPRB + Z1S)*ZQMAX)
        
        ZQSAT = MIN(ZQMAX, ZQSAT)
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
        ZQSAT = ZQSAT*ZCOR
        
        Z2S = ZOEALFA*YDTHF%R5ALVCP*(1.0_JPRB / (ZTARG - YDTHF%R4LES)**2) + (1.0_JPRB - ZOEALFA)*YDTHF%R5ALSCP*(1.0_JPRB / (ZTARG &
        &  - YDTHF%R4IES)**2)
        ZCOND1 = (PQ(JL, KK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*Z2S)
        
        PT(JL, KK) = PT(JL, KK) + (ZOEALFA*YDTHF%RALVDCP + (1.0_JPRB - ZOEALFA)*YDTHF%RALSDCP)*ZCOND1
        
        PQ(JL, KK) = PQ(JL, KK) - ZCOND1
      END IF
      
    END IF
    
    !*********************************************
  END IF
  !*********************************************
  
  
  !IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
  
END SUBROUTINE CUADJTQ_OPENACC
