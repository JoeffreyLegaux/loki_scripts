SUBROUTINE CUDTDQN_OPENACC (YDTHF, YDCST, YDEPHLI, YDPHNC, YDECUMF, YDEPHY, KIDIA, KFDIA, KLON, KLEV, KTOPM2, KTYPE, KCTOP,  &
& KCBOT, KDTOP, LDTDKMF, LDCUM, LDDRAF, LSCVFLAG, PTSPHY, PAPH, PAP, PGEOH, PGEO, PTEN, PTENH, PQEN, PQENH, PQSEN, PLGLAC,  &
& PLUDE, PLUDELI, PSNDE, PMFU, PMFD, PMFUS, PMFDS, PMFUQ, PMFDQ, PMFUL, PDMFUP, PDPMEL, PMFLXR, PMFLXS, PTENT, PTENQ, PENTH,  &
& YDSTACK)
  
  !**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
  !                DOES GLOBAL DIAGNOSTICS
  
  !**   INTERFACE.
  !     ----------
  
  !          *CUDTDQ* IS CALLED FROM *CUMASTR*
  
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
  !    *KCTOP*        CLOUD TOP LEVEL
  !    *KCBOT*        CLOUD BASE LEVEL
  !    *KDTOP*        TOP LEVEL OF DOWNDRAFTS
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  !    *LDDRAF*       FLAG: .TRUE. FOR DOWNDRAFT LEVEL
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
  !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
  !    *PAP *         PROVISIONAL PRESSURE ON FULL LEVELS            PA
  !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
  !    *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2
  !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
  !    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
  !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
  !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
  !    *PQSEN*        SATURATION ENV. SPEC. HUMIDITY (T+1)          KG/KG
  !    *PLGLAC*       FLUX OF FROZEN CLOUDWATER IN UPDRAFTS         KG/(M2*S)
  !    *PLUDE*        DETRAINED CONDENSATE                          KG/(M2*S)
  !    *PLUDELI*      DETRAINED LIQUID, ICE                         KG/(M2*S)
  !    *PSNDE*        DETRAINED SNOW/RAIN                           KG/(M2*S)
  !    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
  !    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
  !    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
  !    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
  !    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
  !    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
  !    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
  !    *PDMFUP*       FLUX DIFFERENCE OF PRECIP.                    KG/(M2*S)
  !    *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)
  
  !    UPDATED PARAMETERS (REAL):
  
  !    *PTENT*        TEMPERATURE TENDENCY                           K/S
  !    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
  
  !----------------------------------------------------------------------
  
  !     AUTHOR.
  !     -------
  !      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
  
  !     MODIFICATIONS.
  !     --------------
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      03-08-28       : Clean up LINPHYS  P.BECHTOLD
  !      05-10-13       : implicit solution P.BECHTOLD
  !      10-11-17       : allow to substract subbsidence (dyn) P. Bechtold
  !      28-12-21       : correction for enthalpy conservation ice phase P. Bechtold
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !---------------------------------------------------------------------------------
  
!$acc routine( CUDTDQN_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOECUMF, ONLY: TECUMF
  USE YOEPHY, ONLY: TEPHY
  USE YOEPHLI, ONLY: TEPHLI
  USE YOPHNC, ONLY: TPHNC
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  TYPE(TEPHLI), INTENT(IN) :: YDEPHLI
  TYPE(TEPHY), INTENT(IN) :: YDEPHY
  TYPE(TPHNC), INTENT(IN) :: YDPHNC
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KTOPM2
  INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KDTOP(KLON)
  LOGICAL, INTENT(INOUT) :: LDCUM(KLON)
  LOGICAL, INTENT(IN) :: LDDRAF(KLON)
  LOGICAL, INTENT(IN) :: LDTDKMF
  LOGICAL, INTENT(IN) :: LSCVFLAG(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEO(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PGEOH(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PTEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PTENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQENH(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQSEN(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLGLAC(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PLUDELI(KLON, KLEV, 2)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSNDE(KLON, KLEV, 2)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFD(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFUS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFDS(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFUQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFDQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFUL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDMFUP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDPMEL(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFLXR(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN) :: PMFLXS(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PTENQ(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PENTH(KLON, KLEV)
  
  LOGICAL :: LLTEST
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: IK
  INTEGER(KIND=JPIM) :: IM
  INTEGER(KIND=JPIM) :: JL
  REAL(KIND=JPRB) :: ZTSPHY
  REAL(KIND=JPRB) :: ZIMP
  REAL(KIND=JPRB) :: ZORCPD
  REAL(KIND=JPRB) :: ZOEALFA
  REAL(KIND=JPRB) :: ZTARG
  REAL(KIND=JPRB) :: ZZP
  REAL(KIND=JPRB) :: ZGQ
  REAL(KIND=JPRB) :: ZGS
  REAL(KIND=JPRB) :: ZGH
  REAL(KIND=JPRB) :: ZS
  REAL(KIND=JPRB) :: ZQ
  temp (REAL (KIND=JPRB), ZALV, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZALV2, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZINT
  REAL(KIND=JPRB) :: ZINT2
  
  temp (REAL (KIND=JPRB), ZMFUS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFUQ, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDS, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZMFDQ, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZADVW
  
  temp (REAL (KIND=JPRB), ZDTDT, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDQDT, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDP, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZB, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZR1, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZR2, (KLON, KLEV))
  temp (LOGICAL, LLCUMBAS, (KLON, KLEV))
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "cubidiag_openacc.intfb.h"
#include "fcttre.func.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  IF (KIND (ZALV) == 8) THEN
    alloc8 (ZALV)
  ELSE
    IF (KIND (ZALV) == 4) THEN
      alloc4 (ZALV)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZALV2) == 8) THEN
    alloc8 (ZALV2)
  ELSE
    IF (KIND (ZALV2) == 4) THEN
      alloc4 (ZALV2)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUS) == 8) THEN
    alloc8 (ZMFUS)
  ELSE
    IF (KIND (ZMFUS) == 4) THEN
      alloc4 (ZMFUS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFUQ) == 8) THEN
    alloc8 (ZMFUQ)
  ELSE
    IF (KIND (ZMFUQ) == 4) THEN
      alloc4 (ZMFUQ)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFDS) == 8) THEN
    alloc8 (ZMFDS)
  ELSE
    IF (KIND (ZMFDS) == 4) THEN
      alloc4 (ZMFDS)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZMFDQ) == 8) THEN
    alloc8 (ZMFDQ)
  ELSE
    IF (KIND (ZMFDQ) == 4) THEN
      alloc4 (ZMFDQ)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDTDT) == 8) THEN
    alloc8 (ZDTDT)
  ELSE
    IF (KIND (ZDTDT) == 4) THEN
      alloc4 (ZDTDT)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDQDT) == 8) THEN
    alloc8 (ZDQDT)
  ELSE
    IF (KIND (ZDQDT) == 4) THEN
      alloc4 (ZDQDT)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZDP) == 8) THEN
    alloc8 (ZDP)
  ELSE
    IF (KIND (ZDP) == 4) THEN
      alloc4 (ZDP)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZB) == 8) THEN
    alloc8 (ZB)
  ELSE
    IF (KIND (ZB) == 4) THEN
      alloc4 (ZB)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZR1) == 8) THEN
    alloc8 (ZR1)
  ELSE
    IF (KIND (ZR1) == 4) THEN
      alloc4 (ZR1)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZR2) == 8) THEN
    alloc8 (ZR2)
  ELSE
    IF (KIND (ZR2) == 4) THEN
      alloc4 (ZR2)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (LLCUMBAS) == 8) THEN
    alloc8 (LLCUMBAS)
  ELSE
    IF (KIND (LLCUMBAS) == 4) THEN
      alloc4 (LLCUMBAS)
    ELSE
      STOP 1
    END IF
  END IF
  
  !----------------------------------------------------------------------
  
  !*    1.0          SETUP AND INITIALIZATIONS
  !                  -------------------------
  
  ZIMP = 1.0_JPRB - YDECUMF%RMFSOLTQ
  ZTSPHY = 1.0_JPRB / PTSPHY
  ZORCPD = 1.0_JPRB / YDCST%RCPD
  ZADVW = 0.0_JPRB
  IF (KTYPE(JLON) == 1) ZADVW = YDECUMF%RMFADVW
  
  DO JK=1,KLEV
    PENTH(JLON, JK) = 0.0_JPRB
  END DO
  
  !         MASS-FLUX APPROACH SWITCHED ON FOR DEEP CONVECTION ONLY
  !         IN THE TANGENT-LINEAR AND ADJOINT VERSIONS
  
  IF (KTYPE(JLON) /= 1 .and. YDEPHLI%LPHYLIN) LDCUM(JLON) = .false.
  
  ! zero detrained liquid water if diagnostic cloud scheme to be used
  
  ! this means that detrained liquid water will be evaporated in the
  ! cloud environment and not fed directly into a cloud liquid water
  ! variable
  
  LLTEST = .not.YDEPHY%LEPCLD .and. .not.YDPHNC%LENCLD2 .and. .not.YDPHNC%LEPCLD2 .or. YDEPHLI%LPHYLIN .and. .not.YDPHNC%LENCLD2  &
  & .and. .not.YDPHNC%LEPCLD2
  
  IF (LLTEST) THEN
    DO JK=1,KLEV
      PLUDE(JLON, JK) = 0.0_JPRB
      PLUDELI(JLON, JK, 1) = 0.0_JPRB
      PLUDELI(JLON, JK, 2) = 0.0_JPRB
    END DO
  END IF
  
  DO JK=1,KLEV
    IF (LDCUM(JLON)) THEN
      ZDP(JLON, JK) = YDCST%RG / (PAPH(JLON, JK + 1) - PAPH(JLON, JK))
      ZMFUS(JLON, JK) = PMFUS(JLON, JK)
      ZMFDS(JLON, JK) = PMFDS(JLON, JK)
      ZMFUQ(JLON, JK) = PMFUQ(JLON, JK)
      ZMFDQ(JLON, JK) = PMFDQ(JLON, JK)
    END IF
  END DO
  
  !------------------------------------------------------------------------------
  
  IF (YDECUMF%RMFSOLTQ > 0.0_JPRB) THEN
    
    !*    2.0          RECOMPUTE CONVECTIVE FLUXES IF IMPLICIT
    
    DO JK=KTOPM2,KLEV
      IK = JK - 1
      !DIR$ IVDEP
      !OCL NOVREC
      IF (LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1) THEN
        ! compute interpolating coefficients ZGS and ZGQ for half-level values
        ZGQ = (PQENH(JLON, JK) - PQEN(JLON, IK)) / PQSEN(JLON, JK)
        ZGH = YDCST%RCPD*PTEN(JLON, JK) + PGEO(JLON, JK)
        ZGS = (YDCST%RCPD*(PTENH(JLON, JK) - PTEN(JLON, IK)) + PGEOH(JLON, JK) - PGEO(JLON, IK)) / ZGH
        
        !half-level environmental values for S and Q
        ZS = YDCST%RCPD*(ZIMP*PTEN(JLON, IK) + ZGS*PTEN(JLON, JK)) + PGEO(JLON, IK) + ZGS*PGEO(JLON, JK)
        ZQ = ZIMP*PQEN(JLON, IK) + ZGQ*PQSEN(JLON, JK)
        ZMFUS(JLON, JK) = PMFUS(JLON, JK) - PMFU(JLON, JK)*ZS
        ZMFUQ(JLON, JK) = PMFUQ(JLON, JK) - PMFU(JLON, JK)*ZQ
        IF (LDDRAF(JLON) .and. JK >= KDTOP(JLON)) THEN
          ZMFDS(JLON, JK) = PMFDS(JLON, JK) - PMFD(JLON, JK)*ZS
          ZMFDQ(JLON, JK) = PMFDQ(JLON, JK) - PMFD(JLON, JK)*ZQ
        END IF
      END IF
    END DO
    
  END IF
  
  !*    3.0          COMPUTE TENDENCIES
  !                  ------------------
  
  ZINT = 0._JPRB
  DO JK=KTOPM2,KLEV
    ZALV(JLON, JK) = FOELHMCU(PTEN(JLON, JK))
    ZALV2(JLON, JK) = ZALV(JLON, JK)
  END DO
  IF (YDECUMF%LMFENTHCONS) THEN
    ! Integrals for conservation correction ice-phase
    DO JK=KTOPM2,KLEV - 1
      IF (LSCVFLAG(JLON)) ZALV(JLON, JK) = YDCST%RLVTT
      ZINT = ZINT + YDCST%RLMLT*(PLGLAC(JLON, JK) - PDPMEL(JLON, JK)) - ZALV(JLON, JK)*(PMFUL(JLON, JK + 1) - PMFUL(JLON, JK)) +  &
      & ZALV2(JLON, JK)*PDMFUP(JLON, JK)
    END DO
    JK = KLEV
    IF (LSCVFLAG(JLON)) ZALV(JLON, JK) = YDCST%RLVTT
    ZINT = ZINT + YDCST%RLMLT*(PLGLAC(JLON, JK) - PDPMEL(JLON, JK)) + ZALV(JLON, JK)*PMFUL(JLON, JK) + ZALV2(JLON, JK) &
    & *PDMFUP(JLON, JK)
    ZINT2 = ZINT - YDCST%RLVTT*PMFLXR(JLON, KLEV + 1) - YDCST%RLSTT*PMFLXS(JLON, KLEV + 1)
    IF (ABS(ZINT) > 0.1_JPRB) THEN
      ZINT = ZINT2 / (ZINT + 1.E-3_JPRB)
    ELSE
      ZINT = 0.0_JPRB
    END IF
    IF (ABS(ZINT) > 2.0_JPRB) ZINT = 0.0_JPRB
  END IF
  
  DO JK=KTOPM2,KLEV
    
    IF (JK < KLEV) THEN
      IF (LDCUM(JLON)) THEN
        IF (YDEPHLI%LPHYLIN) THEN
          ZTARG = PTEN(JLON, JK)
          ZOEALFA = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB))
          ZALV(JLON, JK) = ZOEALFA*YDCST%RLVTT + (1.0_JPRB - ZOEALFA)*YDCST%RLSTT
          ZALV2(JLON, JK) = ZALV(JLON, JK)
        END IF
        IF (LDTDKMF) THEN
          ZDTDT(JLON, JK) = ZDP(JLON, JK)*ZORCPD*(ZMFUS(JLON, JK + 1) - ZMFUS(JLON, JK) + ZMFDS(JLON, JK + 1) - ZMFDS(JLON, JK) + &
          &  YDCST%RLMLT*PLGLAC(JLON, JK) - YDCST%RLMLT*PDPMEL(JLON, JK) - ZALV(JLON, JK)*(PMFUL(JLON, JK + 1) - PMFUL(JLON, JK)  &
          & - PLUDE(JLON, JK) - PDMFUP(JLON, JK)))
          
          ZDQDT(JLON, JK) = ZDP(JLON, JK)*(ZMFUQ(JLON, JK + 1) - ZMFUQ(JLON, JK) + ZMFDQ(JLON, JK + 1) - ZMFDQ(JLON, JK) +  &
          & PMFUL(JLON, JK + 1) - PMFUL(JLON, JK) - PLUDE(JLON, JK) - PDMFUP(JLON, JK))
        ELSE
          ZDTDT(JLON, JK) = ZDP(JLON, JK)*ZORCPD*(ZMFUS(JLON, JK + 1) - ZMFUS(JLON, JK) + ZMFDS(JLON, JK + 1) - ZMFDS(JLON, JK) + &
          &  (YDCST%RLMLT*(PLGLAC(JLON, JK) - PDPMEL(JLON, JK)) - ZALV(JLON, JK)*(PMFUL(JLON, JK + 1) - PMFUL(JLON, JK)) +  &
          & ZALV2(JLON, JK)*PDMFUP(JLON, JK))*(1.0_JPRB - ZINT) + YDCST%RLVTT*(PLUDELI(JLON, JK, 1) + PSNDE(JLON, JK, 1)) +  &
          & YDCST%RLSTT*(PLUDELI(JLON, JK, 2) + PSNDE(JLON, JK, 2)))
          
          ZDQDT(JLON, JK) = ZDP(JLON, JK)*(ZMFUQ(JLON, JK + 1) - ZMFUQ(JLON, JK) + ZMFDQ(JLON, JK + 1) - ZMFDQ(JLON, JK) +  &
          & PMFUL(JLON, JK + 1) - PMFUL(JLON, JK) - PLUDE(JLON, JK) - PDMFUP(JLON, JK) - PSNDE(JLON, JK, 1) - PSNDE(JLON, JK, 2))
        END IF
      END IF
      
    ELSE
      IF (LDCUM(JLON)) THEN
        IF (YDEPHLI%LPHYLIN) THEN
          ZTARG = PTEN(JLON, JK)
          ZOEALFA = MIN(1.0_JPRB, 0.545_JPRB*(TANH(0.17_JPRB*(ZTARG - YDEPHLI%RLPTRC)) + 1.0_JPRB))
          ZALV(JLON, JK) = ZOEALFA*YDCST%RLVTT + (1.0_JPRB - ZOEALFA)*YDCST%RLSTT
          ZALV2(JLON, JK) = ZALV(JLON, JK)
        END IF
        ZDTDT(JLON, JK) = -ZDP(JLON, JK)*ZORCPD*(ZMFUS(JLON, JK) + ZMFDS(JLON, JK) + (YDCST%RLMLT*(PDPMEL(JLON, JK) -  &
        & PLGLAC(JLON, JK)) - ZALV(JLON, JK)*PMFUL(JLON, JK) - ZALV2(JLON, JK)*PDMFUP(JLON, JK))*(1.0_JPRB - ZINT))
        
        ZDQDT(JLON, JK) = -ZDP(JLON, JK)*(ZMFUQ(JLON, JK) + ZMFDQ(JLON, JK) + (PMFUL(JLON, JK) + PDMFUP(JLON, JK)))
      END IF
    END IF
    
  END DO
  
  IF (YDECUMF%RMFSOLTQ == 0.0_JPRB) THEN
    
    !*    3.1          UPDATE TENDENCIES
    !                  -----------------
    
    DO JK=KTOPM2,KLEV
      IF (LDCUM(JLON)) THEN
        PTENT(JLON, JK) = PTENT(JLON, JK) + ZDTDT(JLON, JK)
        PTENQ(JLON, JK) = PTENQ(JLON, JK) + ZDQDT(JLON, JK)
        PENTH(JLON, JK) = ZDTDT(JLON, JK)*YDCST%RCPD
      END IF
    END DO
    
  ELSE
    !----------------------------------------------------------------------
    
    !*    3.2          IMPLICIT SOLUTION
    !                  -----------------
    
    ! Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1;
    ! reuse ZMFUS=A
    ! ZDTDT and ZDQDT correspond to the RHS ("constants") of the equation
    ! The solution is in ZR1 and ZR2
    
    LLCUMBAS(JLON, :) = .false.
    ZB(JLON, :) = 1.0_JPRB
    ZMFUS(JLON, :) = 0.0_JPRB
    
    ! Fill vectors A, B and RHS
    
    DO JK=KTOPM2,KLEV
      IK = JK + 1
      IM = JK - 1
      LLCUMBAS(JLON, JK) = LDCUM(JLON) .and. JK >= KCTOP(JLON) - 1
      IF (LLCUMBAS(JLON, JK)) THEN
        ZZP = YDECUMF%RMFSOLTQ*ZDP(JLON, JK)*PTSPHY
        ZMFUS(JLON, JK) = -ZZP*(PMFU(JLON, JK) + PMFD(JLON, JK))
        IF (JK < KLEV) THEN
          ZB(JLON, JK) = 1.0_JPRB + ZZP*(PMFU(JLON, IK) + PMFD(JLON, IK))
        ELSE
          ZB(JLON, JK) = 1.0_JPRB
        END IF
        IF (LDTDKMF) THEN
          ZDTDT(JLON, JK) = ZDTDT(JLON, JK)*PTSPHY + PTEN(JLON, JK)
          ZDQDT(JLON, JK) = ZDQDT(JLON, JK)*PTSPHY + PQEN(JLON, JK)
        ELSE
          ZZP = YDCST%RG*(PMFU(JLON, JK) + YDECUMF%RMFADVWDD*PMFD(JLON, JK)) / (PAP(JLON, JK) - PAP(JLON, IM))*PTSPHY*ZADVW
          ZS = ZZP*(PTEN(JLON, IM) - PTEN(JLON, JK) + ZORCPD*(PGEO(JLON, IM) - PGEO(JLON, JK)))
          ZQ = ZZP*(PQEN(JLON, IM) - PQEN(JLON, JK))
          ZDTDT(JLON, JK) = (ZDTDT(JLON, JK) + PTENT(JLON, JK)*YDECUMF%RMFSOLRHS)*PTSPHY + PTEN(JLON, JK) - ZS
          ZDQDT(JLON, JK) = (ZDQDT(JLON, JK) + PTENQ(JLON, JK)*YDECUMF%RMFSOLRHS)*PTSPHY + PQEN(JLON, JK) - ZQ
        END IF
      END IF
    END DO
    
    CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFUS, ZB, ZDTDT, ZR1, YDSTACK=YLSTACK)
    
    CALL CUBIDIAG_OPENACC(KIDIA, KFDIA, KLON, KLEV, KCTOP, LLCUMBAS, ZMFUS, ZB, ZDQDT, ZR2, YDSTACK=YLSTACK)
    
    ! Compute tendencies
    
    DO JK=KTOPM2,KLEV
      IF (LLCUMBAS(JLON, JK)) THEN
        IF (LDTDKMF) THEN
          PTENT(JLON, JK) = PTENT(JLON, JK) + (ZR1(JLON, JK) - PTEN(JLON, JK))*ZTSPHY
          PTENQ(JLON, JK) = PTENQ(JLON, JK) + (ZR2(JLON, JK) - PQEN(JLON, JK))*ZTSPHY
        ELSE
          PTENT(JLON, JK) = PTENT(JLON, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + (ZR1(JLON, JK) - PTEN(JLON, JK))*ZTSPHY
          PTENQ(JLON, JK) = PTENQ(JLON, JK)*(1.0_JPRB - YDECUMF%RMFSOLRHS) + (ZR2(JLON, JK) - PQEN(JLON, JK))*ZTSPHY
        END IF
        PENTH(JLON, JK) = (ZR1(JLON, JK) - PTEN(JLON, JK))*ZTSPHY
      END IF
    END DO
    
    !----------------------------------------------------------------------
  END IF
  
END SUBROUTINE CUDTDQN_OPENACC
