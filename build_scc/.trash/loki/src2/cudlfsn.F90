SUBROUTINE CUDLFSN_OPENACC (YDTHF, YDCST, YDEPHLI, YDECUMF, KIDIA, KFDIA, KLON, KLEV, KCBOT, KCTOP, LDCUM, PTENH, PQENH, PTEN,  &
& PQSEN, PGEO, PGEOH, PAPH, PTU, PQU, PMFUB, PRFL, PTD, PQD, PMFD, PMFDS, PMFDQ, PDMFDP, KDTOP, LDDRAF, YDSTACK)
  
  !          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
  !          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
  
  !          PURPOSE.
  !          --------
  !          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
  !          FOR MASSFLUX CUMULUS PARAMETERIZATION
  
  !          INTERFACE
  !          ---------
  !          THIS ROUTINE IS CALLED FROM *CUMASTR*.
  !          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
  !          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
  !          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
  !          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
  
  !          METHOD.
  
  !          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
  !          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  !    *KCBOT*        CLOUD BASE LEVEL
  !    *KCTOP*        CLOUD TOP LEVEL
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
  !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
  !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
  !    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
  !    *PGEO*         GEOPOTENTIAL                                  M2/S2
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
  !    *PTU*          TEMPERATURE IN UPDRAFTS                        K
  !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
  !    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
  !    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
  !    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
  !    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
  !    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
  !    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)
  
  !    OUTPUT PARAMETERS (INTEGER):
  
  !    *KDTOP*        TOP LEVEL OF DOWNDRAFTS
  
  !    OUTPUT PARAMETERS (LOGICAL):
  
  !    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST
  
  !          EXTERNALS
  !          ---------
  !          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      P. Lopez      20-Jun-2007 CY32R2 Bug correction in latent heat when LPHYLIN=T.
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !----------------------------------------------------------------------
  
!$acc routine( CUDLFSN_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOECUMF, ONLY: TECUMF
  USE YOEPHLI, ONLY: TEPHLI
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP
  LOGICAL, INTENT(IN) :: LDCUM
  REAL(KIND=JPRB), INTENT(IN) :: PTENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PTU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFUB
  REAL(KIND=JPRB), INTENT(INOUT) :: PRFL
  REAL(KIND=JPRB), INTENT(OUT) :: PTD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PQD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PMFD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFDS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PMFDQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PDMFDP(KLON, KLEV)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KDTOP
  LOGICAL, INTENT(OUT) :: LDDRAF
  INTEGER(KIND=JPIM) :: IKHSMIN
  temp (REAL (KIND=JPRB), ZTENWB, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZQENWB, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZCOND
  REAL(KIND=JPRB) :: ZPH
  REAL(KIND=JPRB) :: ZHSMIN
  LOGICAL :: LLO2
  
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IKE
  INTEGER(KIND=JPIM) :: IS
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  
  REAL(KIND=JPRB) :: ZBUO
  REAL(KIND=JPRB) :: ZHSK
  REAL(KIND=JPRB) :: ZMFTOP
  REAL(KIND=JPRB) :: ZOEALFA
  REAL(KIND=JPRB) :: ZOELHM
  REAL(KIND=JPRB) :: ZQTEST
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZTTEST
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cuadjtq.intfb.h"
  
#include "fcttre.func.h"
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ZTENWB)
  alloc (ZQENWB)
  JL = KIDIA
  !----------------------------------------------------------------------
  
  !     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
  !                  ---------------------------------
  
  LDDRAF = .false.
  KDTOP = KLEV + 1
  IKHSMIN = KLEV + 1
  ZHSMIN = 1.E8_JPRB
  
  !orig IF(.NOT.LMFDD) GO TO 300
  IF (YDECUMF%LMFDD) THEN
    
    !----------------------------------------------------------------------
    
    !     2.           DETERMINE LEVEL OF FREE SINKING:
    !                  DOWNDRAFTS SHALL START AT MODEL LEVEL OF MINIMUM
    !                  OF SATURATION MOIST STATIC ENERGY OR BELOW
    !                  RESPECTIVELY
    
    !                  FOR EVERY POINT AND PROCEED AS FOLLOWS:
    
    !                    (1) DETERMINE LEVEL OF MINIMUM OF HS
    !                    (2) DETERMINE WET BULB ENVIRONMENTAL T AND Q
    !                    (3) DO MIXING WITH CUMULUS CLOUD AIR
    !                    (4) CHECK FOR NEGATIVE BUOYANCY
    !                    (5) IF BUOYANCY>0 REPEAT (2) TO (4) FOR NEXT
    !                        LEVEL BELOW
    
    !                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
    !                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
    !                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
    !                  EVAPORATION OF RAIN AND CLOUD WATER)
    !                  ----------------------------------------------------
    
    DO JK=3,KLEV - 2
      
      IF (YDEPHLI%LPHYLIN) THEN
        
        ZTARG = PTEN(JL, JK)
        ZOEALFA = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB))
        ZOELHM = ZOEALFA*YDCST%RLVTT + (1.0_JPRB - ZOEALFA)*YDCST%RLSTT
        ZHSK = YDCST%RCPD*PTEN(JL, JK) + PGEO(JL, JK) + ZOELHM*PQSEN(JL, JK)
        IF (ZHSK < ZHSMIN) THEN
          ZHSMIN = ZHSK
          IKHSMIN = JK
        END IF
        
      ELSE
        
        ZHSK = YDCST%RCPD*PTEN(JL, JK) + PGEO(JL, JK) + FOELHMCU(PTEN(JL, JK))*PQSEN(JL, JK)
        IF (ZHSK < ZHSMIN) THEN
          ZHSMIN = ZHSK
          IKHSMIN = JK
        END IF
        
      END IF
      
    END DO
    IKE = KLEV - 3
    DO JK=3,IKE
      
      !     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
      !                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
      !                  -------------------------------------------
      
      IS = 0
      ZTENWB(JL, JK) = PTENH(JL, JK)
      ZQENWB(JL, JK) = PQENH(JL, JK)
      ZPH = PAPH(JL, JK)
      LLO2 = LDCUM .and. PRFL > 0.0_JPRB .and. .not.LDDRAF .and. JK < KCBOT .and. JK > KCTOP .and. JK >= IKHSMIN
      IF (LLO2) THEN
        IS = IS + 1
      END IF
      !orig   IF(IS.EQ.0) GO TO 290
      IF (IS == 0) CYCLE
      
      IK = JK
      CALL CUADJTQ_OPENACC(YDTHF, YDCST, YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK, ZPH, ZTENWB, ZQENWB, LLO2, KCALL=2,  &
      & YDSTACK=YLSTACK)
      
      !     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
      !                  AND CHECK FOR NEGATIVE BUOYANCY.
      !                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
      !                  ----------------------------------------
      
      !DIR$ IVDEP
      !OCL NOVREC
      IF (LLO2) THEN
        ZTTEST = 0.5_JPRB*(PTU(JL, JK) + ZTENWB(JL, JK))
        ZQTEST = 0.5_JPRB*(PQU(JL, JK) + ZQENWB(JL, JK))
        ZBUO = ZTTEST*(1.0_JPRB + YDCST%RETV*ZQTEST) - PTENH(JL, JK)*(1.0_JPRB + YDCST%RETV*PQENH(JL, JK))
        ZCOND = PQENH(JL, JK) - ZQENWB(JL, JK)
        ZMFTOP = -YDECUMF%RMFDEPS*PMFUB
        IF (ZBUO < 0.0_JPRB .and. PRFL > 10._JPRB*ZMFTOP*ZCOND) THEN
          KDTOP = JK
          LDDRAF = .true.
          PTD(JL, JK) = ZTTEST
          PQD(JL, JK) = ZQTEST
          PMFD(JL, JK) = ZMFTOP
          PMFDS(JL, JK) = PMFD(JL, JK)*(YDCST%RCPD*PTD(JL, JK) + PGEOH(JL, JK))
          PMFDQ(JL, JK) = PMFD(JL, JK)*PQD(JL, JK)
          PDMFDP(JL, JK - 1) = -0.5_JPRB*PMFD(JL, JK)*ZCOND
          PRFL = PRFL + PDMFDP(JL, JK - 1)
        END IF
      END IF
      
      ! 290   continue
    END DO
    
    !300  CONTINUE
  END IF
  
END SUBROUTINE CUDLFSN_OPENACC
