SUBROUTINE CULIGHT_OPENACC (PPLDARE, PPLRG, YDTHF, YDCST, YDEPHY, YGFL, YDECUMF, KIDIA, KFDIA, KLON, KLEV, PGAW, PGELAT, PAP,  &
& PAPH, PAPHI, PAPHIF, LDLAND, PT, PLU, PMFU, PCAPE, PFPLCL, PFPLCN, PQPFROZ, LDCUM, KCBOT, KCTOP, LDLINOX, PLIGH_TOT,  &
& PLIGH_CTG, PCTOPH, PPRECMX, PICE, PCDEPTH, PWMFU, PCHARGE, YDSTACK)
  ! Outputs
  
  !    THIS ROUTINE CALCULATES LIGHTNING FLASH RATES.
  
  !    P. Lopez     ECMWF   (03/2005)
  !
  !
  !    PURPOSE.
  !    --------
  !    TO CALCULATE LIGHTNING FLASH RATES (FLASHES/KM2/DAY)
  !    FROM INFORMATION COMING FROM THE CONVECTION SCHEME.
  !    TOTAL AND CLOUD-TO-GROUND FLASH RATES ARE COMPUTED.
  !    Partly inspired from Kurz and Grewe (2002).
  !
  !    The following lightning parameterizations are available:
  !
  !      NLIMODE = 1  -> Lopez, 2005, version 1.
  !                2  -> Lopez, 2005, version 2.
  !                3  -> Lopez, 2005, version 3.
  !                4  -> TM5 Meijer, 2001.
  !                5  -> Price and Rind, 1993.
  !                6  -> Lopez, 2015, latest version.
  !                7  -> Grewe et al., 2001.
  !                8  -> Allen and Prickering, 2002.
  
  !    INTERFACE
  !    ---------
  !    THIS ROUTINE IS CALLED FROM *LIGHTNING_LAYER*.
  
  !    METHOD.
  !    --------
  
  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):
  !
  !    *KIDIA*        START POINT
  !    *KFDIA*        END POINT
  !    *KLON*         NUMBER OF GRID POINTS PER PACKET
  !    *KLEV*         NUMBER OF LEVELS
  
  !    INPUT PARAMETERS (REAL):
  
  !    *PGAW*         GAUSSIAN WEIGHTS Reduced ~ grid box area
  !    *PGELAT*       Latitude (radians)
  !    *PAP*          PRESSURE ON FULL LEVELS                PA
  !    *PAPH*         PRESSURE ON HALF LEVELS                PA
  !    *PAPHI*        GEOPOTENTIAL ON HALF LEVELS            M2/S2
  !    *PAPHIF*       GEOPOTENTIAL ON FULL LEVELS            M2/S2
  !    *PT*           TEMPERATURE                            K
  !    *PLU*          CONVECTIVE CLOUD CONDENSATE            KG/KG
  !    *PMFU*         CONVECTIVE UPDRAUGHT MASS FLUX         KG/(SM2)
  !    *PFPLCL*       CONVECTIVE RAIN FLUX                   KG/(SM2)
  !    *PFPLCN*       CONVECTIVE SNOW FLUX                   KG/(SM2)
  !    *PQPFROZ*      CONVECTIVE FROZEN PRECIP CONTENT       KG/KG
  
  !    INPUT PARAMETERS (LOGICAL):
  
  !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
  !    *LDLAND*       LAND-SEA MASK
  
  !    INPUT PARAMETERS (INTEGER):
  
  !    *KCBOT*       CONVECTIVE CLOUD BASE LEVEL
  !    *KCTOP*       CONVECTIVE CLOUD TOP LEVEL
  
  !    OUTPUT PARAMETERS (LOGICAL):
  
  !    *LDLINOX*     GRID-POINT FLAG: .TRUE. FOR LIGHTNING NOx COMPUTATIONS IN CULINOX
  
  !    OUTPUT PARAMETERS (REAL):
  
  !    *PLIGH_TOT*   TOTAL LIGHTNING FLASH RATES               FL/KM2/DAY
  !    *PLIGH_CTG*   CLOUD-TO-GROUND LIGHTNING FLASH RATES     FL/KM2/DAY
  !    *PCTOPH*      CONVECTIVE CLOUD TOP HEIGHT               KM
  !    *PPRECMX*     MAXIMUM CONVECTIVE PRECIPITATION IN THE   KG/(SM2)
  !                  VERTICAL.
  !    *PICE*        TOTAL CONVECTIVE CLOUD ICE CONTENT        KG/M2
  !    *PCDEPTH*     DEPTH OF CLOUD ABOVE FREEZING LEVEL       KM
  !    *PWMFU*       MEAN UPDRAUGHT MASS FLUX                  M/S
  !    *PCHARGE*     ELECTRIC CHARGE PRODUCTION RATE
  
  !    EXTERNALS
  !    ---------
  
  !    MODIFICATIONS
  !    -------------
  !    J Flemming   ECMWF   (04/2010) ! Prince & Rind 1993 and Meijer 2001 added
  !    P Bechtold 18/05/2012   Use RDAYI for 86400
  !    J  Flemming 03/07/2013  Remove difference in CTG and IC lighting , introdruce new profile according to Ott et al., 2010
  !    N. Semane+P.Bechtold 04/10/2012   add RPLRG for small planet
  !    P. Lopez 24/07/2015     New ECMWF lightning parameterization plus added Grewe (2001) and Allen-Pickering (2002).
  !    P. Lopez 08/10/2015     Moved lightning NOx computations to CULINOX.
  !    P. Lopez 25/03/2021     Extended usage to the case of explicit convection (i.e. at very high resolutions).
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !
  !----------------------------------------------------------------------
  
!$acc routine( CULIGHT_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE YOMCST, ONLY: TCST
  USE YOETHF, ONLY: TTHF
  USE YOEPHY, ONLY: TEPHY
  USE YOM_YGFL, ONLY: TYPE_GFLD
  USE YOECUMF, ONLY: TECUMF
  ! USE YOMLUN    , ONLY : NULOUT
  
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  REAL(KIND=JPRB), INTENT(IN) :: PPLDARE
  REAL(KIND=JPRB), INTENT(IN) :: PPLRG
  TYPE(TTHF), INTENT(IN) :: YDTHF
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TEPHY), INTENT(IN) :: YDEPHY
  TYPE(TYPE_GFLD), INTENT(IN) :: YGFL
  TYPE(TECUMF), INTENT(IN) :: YDECUMF
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(IN) :: PGAW(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PGELAT(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHI(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PLU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PCAPE(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PFPLCL(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PFPLCN(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PQPFROZ(KLON, KLEV)
  LOGICAL, INTENT(IN) :: LDCUM(KLON)
  LOGICAL, INTENT(IN) :: LDLAND(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT(KLON)
  INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP(KLON)
  LOGICAL, INTENT(OUT) :: LDLINOX(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PLIGH_TOT(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PLIGH_CTG(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PCTOPH(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PPRECMX(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PICE(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PCDEPTH(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PWMFU(KLON)
  REAL(KIND=JPRB), INTENT(OUT) :: PCHARGE(KLON)
  
#include "abor1.intfb.h"
  
  !             LOCAL STORAGE
  !             ----- -------
  
  INTEGER(KIND=JPIM) :: JK
  INTEGER(KIND=JPIM) :: JL
  INTEGER(KIND=JPIM) :: IFREEZ
  INTEGER(KIND=JPIM) :: IFREEZ25
  INTEGER(KIND=JPIM) :: IFREEZ15
  INTEGER(KIND=JPIM) :: ILEVTOPG
  INTEGER(KIND=JPIM) :: ILEVBOTG
  
  REAL(KIND=JPRB) :: ZGEO2KM
  REAL(KIND=JPRB) :: ZBASHG
  REAL(KIND=JPRB) :: ZPCTG
  REAL(KIND=JPRB) :: Z1G
  REAL(KIND=JPRB) :: ZDEPTH
  REAL(KIND=JPRB) :: ZRHO
  REAL(KIND=JPRB) :: ZAREA_REF
  REAL(KIND=JPRB) :: ZAREA
  REAL(KIND=JPRB) :: ZCDEPTH
  REAL(KIND=JPRB) :: ZPREC_CV
  REAL(KIND=JPRB) :: ZCBASEH
  REAL(KIND=JPRB) :: ZQGRAUP
  REAL(KIND=JPRB) :: ZQSNOW
  REAL(KIND=JPRB) :: ZBETA
  REAL(KIND=JPRB) :: ZPFROZEN
  REAL(KIND=JPRB) :: ZCOLDDPT
  REAL(KIND=JPRB) :: ZHBMAX
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  
#include "fcttre.func.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !----------------------------------------------------------------------
  
  
  !----------------------------------------------------------------------
  !     0.           INITIALIZE CONSTANTS AND FIELDS
  !----------------------------------------------------------------------
  
  ZHBMAX = 1.8_JPRB
  Z1G = 1.0_JPRB / YDCST%RG
  ZGEO2KM = 0.001_JPRB*Z1G
  
  ! Initializations
  PLIGH_TOT(JLON) = 0.0_JPRB
  PLIGH_CTG(JLON) = 0.0_JPRB
  PCTOPH(JLON) = -99.0_JPRB
  PPRECMX(JLON) = 0.0_JPRB
  PICE(JLON) = 0.0_JPRB
  PCDEPTH(JLON) = 0.0_JPRB
  PWMFU(JLON) = -99.0_JPRB
  ZCDEPTH = 0.0_JPRB
  PCHARGE(JLON) = 0.0_JPRB
  ZPREC_CV = 0.0_JPRB
  ZCBASEH = -99.0_JPRB
  
  LDLINOX(JLON) = .false.
  
  ! grid box size in km**2
  ZAREA = 1.0E-6_JPRB*4.0_JPRB*YDCST%RPI*YDCST%RA*YDCST%RA*PGAW(JLON)
  ! see Allen and Pickering JGR, 2002, 2.0 x 2.5 reference grid box centred at 30 deg latitude.
  !LOP:Incorrect ZAREA_REF=(20000.0/FLOAT(90))*(20000.0/FLOAT(72))
  ZAREA_REF = (20015.1_JPRB / 90.0_JPRB)*(17333.6_JPRB / 72.0_JPRB)
  
  !----------------------------------------------------------------------
  !     1.           CALCULATE LIGHTNING FLASH RATES
  !----------------------------------------------------------------------
  
  IF (YGFL%NCHEM /= 0 .and. YDEPHY%NLIMODE /= 4) THEN
    CALL ABOR1_ACC('CULIGHT: NCHEM /= 0 .and. NLIMODE /= 4 NOT SUPPORTED')
  END IF
  
  
  IF (LDCUM(JLON) .and. KCTOP(JLON) >= 1 .and. KCTOP(JLON) <= KLEV) THEN
    
    ! Compute lightning occurrences only if cloud base height is below
    ! 4 km and if convective precipitation is present at least at one level
    ! between cloud base and cloud top.
    
    ZBASHG = (PAPHIF(JLON, KCBOT(JLON)) - PAPHIF(JLON, KLEV))*ZGEO2KM
    PPRECMX(JLON) = MAXVAL(PFPLCL(JLON, KCTOP(JLON):KCBOT(JLON)) + PFPLCN(JLON, KCTOP(JLON):KCBOT(JLON)))
    
    ! naj is this cloud base <  4 km OK for TM5's and MOZART's parameterisations ?
    IF (ZBASHG < (4.0_JPRB / PPLRG) .and. PPRECMX(JLON) > 1.E-10_JPRB) THEN
      
      ! Set flag for grid-points where NOx computations should be performed in CULINOX.
      LDLINOX(JLON) = .true.
      
      IFREEZ = -1
      IFREEZ25 = -1
      IFREEZ15 = -1
      DO JK=KLEV,1,-1
        IF (PT(JLON, JK) <= YDCST%RTT .and. IFREEZ == -1) IFREEZ = JK
        ! below 0C
        IF (PT(JLON, JK) <= YDCST%RTT - 25.0_JPRB .and. IFREEZ25 == -1) IFREEZ25 = JK
        ! below -25C
        IF (PT(JLON, JK) <= YDCST%RTT - 15.0_JPRB .and. IFREEZ15 == -1) IFREEZ15 = JK
        ! below -15C
      END DO
      
      ! Total cloud depth
      PCTOPH(JLON) = MAX((PAPHI(JLON, KCTOP(JLON)) - PAPHI(JLON, KCBOT(JLON)))*ZGEO2KM, 0.0_JPRB)
      
      ! Cold cloud depth
      PCDEPTH(JLON) = MAX((PAPHI(JLON, KCTOP(JLON)) - PAPHI(JLON, IFREEZ))*ZGEO2KM, 0.0_JPRB)
      PCDEPTH(JLON) = MIN(PCDEPTH(JLON), PCTOPH(JLON))
      
      ! Diagnose updraft velocity for lightning formulae.
      
      IF (YDEPHY%NLIMODE <= 3 .or. YDEPHY%NLIMODE == 7) THEN
        PWMFU(JLON) = 0.0_JPRB
        ZDEPTH = 0.0_JPRB
        DO JK=KCTOP(JLON),KCBOT(JLON)
          IF (PMFU(JLON, JK) > 0.0_JPRB) THEN
            ZRHO = PAP(JLON, JK) / (YDCST%RD*PT(JLON, JK))
            ZDEPTH = ZDEPTH + (PAPHI(JLON, JK - 1) - PAPHI(JLON, JK))
            PWMFU(JLON) = PWMFU(JLON) + PMFU(JLON, JK)*(PAPHI(JLON, JK - 1) - PAPHI(JLON, JK)) / ZRHO
          END IF
        END DO
        IF (ZDEPTH > 0.0_JPRB) THEN
          PWMFU(JLON) = PWMFU(JLON) / ZDEPTH
        ELSE
          PWMFU(JLON) = -99.0_JPRB
        END IF
      END IF
      
      IF (YDEPHY%NLIMODE <= 3) THEN
        ! Total convective cloud ice.
        DO JK=1,KLEV
          PICE(JLON) = PICE(JLON) + MAX(0.0_JPRB, PLU(JLON, JK)*(1.0_JPRB - FOEALFCU(PT(JLON, JK))))*(PAPH(JLON, JK) - PAPH(JLON, &
          &  JK - 1))*Z1G
        END DO
      END IF
      
      IF (YDEPHY%NLIMODE == 6) THEN
        ILEVTOPG = MAX(IFREEZ25, KCTOP(JLON))
        ILEVBOTG = MIN(IFREEZ, KCBOT(JLON))
        
        ! Convective cloud base height.
        ZCBASEH = MAX(0.0_JPRB, (PAPHI(JLON, KCBOT(JLON)) - PAPHI(JLON, KLEV))*ZGEO2KM)
        
        ! Convective cloud depth.
        ZCDEPTH = MAX(0.0_JPRB, (PAPHI(JLON, KCTOP(JLON)) - PAPHI(JLON, KCBOT(JLON)))*ZGEO2KM)
        
        ! Land/sea-dependent factor for the partitioning of frozen precipitation into graupel and snow.
        IF (LDLAND(JLON)) THEN
          ZBETA = 0.70_JPRB
        ELSE
          ZBETA = 0.45_JPRB
        END IF
        
        ! Electric charge production rate.
        DO JK=ILEVTOPG,ILEVBOTG
          ! Fraction of frozen precipitation flux in the form of graupel inside this layer.
          IF (YDECUMF%LMFPEN) THEN
            ! With parameterized convection, first convert input frozen precipitation flux into content
            ! and then partition into graupel and snow.
            ZRHO = PAP(JLON, JK) / (YDCST%RD*PT(JLON, JK))
            ZPFROZEN = MAX(0.0_JPRB, PFPLCN(JLON, JK))
            ZQGRAUP = ZBETA*ZPFROZEN / (ZRHO*3.0_JPRB)
            ZQSNOW = (1.0_JPRB - ZBETA)*ZPFROZEN / (ZRHO*0.5_JPRB)
          ELSE
            ! With explicit convection, partition input frozen precipitation content into
            ! graupel and snow.
            ZQGRAUP = ZBETA*PQPFROZ(JLON, JK)
            ZQSNOW = (1.0_JPRB - ZBETA)*PQPFROZ(JLON, JK)
          END IF
          PCHARGE(JLON) =  &
          & PCHARGE(JLON) + ZQGRAUP*(MAX(0.0_JPRB, PLU(JLON, JK)) + ZQSNOW)*(PAPH(JLON, JK) - PAPH(JLON, JK - 1))*Z1G
        END DO
        
      END IF
      
      IF (YDEPHY%NLIMODE == 8) THEN
        ! For Allen and Prickering (2002) param.
        ! Total convective surface precipitation (in mm/day; capped at 90 mm/day).
        ZPREC_CV = MIN((PFPLCL(JLON, KLEV) + PFPLCN(JLON, KLEV))*86400.0_JPRB, 90.0_JPRB)
      END IF
      
      ! Lightning flash rates.
      
      IF (YDEPHY%NLIMODE == 1) THEN
        ! PL seems TechMemo code.
        IF (PCDEPTH(JLON) > 1.0_JPRB .and. (PCTOPH(JLON) - PCDEPTH(JLON)) > 0.5_JPRB .and. PWMFU(JLON) > 0.0_JPRB) THEN
          IF (LDLAND(JLON)) THEN
            PLIGH_TOT(JLON) = 3.1E-01_JPRB*PICE(JLON)*PWMFU(JLON)*PCDEPTH(JLON)
          ELSE
            PLIGH_TOT(JLON) = 1.3E-02_JPRB*PICE(JLON)*PWMFU(JLON)*PCDEPTH(JLON)
          END IF
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 2) THEN
        ! PL was commented out.
        IF (PCDEPTH(JLON) > 1.0_JPRB .and. (PCTOPH(JLON) - PCDEPTH(JLON)) > 0.7_JPRB .and. PCTOPH(JLON) > 5.0_JPRB .and.  &
        & PWMFU(JLON) > 0.05_JPRB) THEN
          IF (LDLAND(JLON)) THEN
            PLIGH_TOT(JLON) = 6.5E-01_JPRB*PICE(JLON)*PWMFU(JLON)*PCDEPTH(JLON)
          ELSE
            PLIGH_TOT(JLON) = 8.0E-02_JPRB*PICE(JLON)*PWMFU(JLON)*PCDEPTH(JLON)
          END IF
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 3) THEN
        ! PL version in code as is.
        IF (PCDEPTH(JLON) > 1.0_JPRB .and. (PCTOPH(JLON) - PCDEPTH(JLON)) > 0.7_JPRB .and. PCTOPH(JLON) > 5.0_JPRB .and.  &
        & PWMFU(JLON) > 0.05_JPRB) THEN
          IF (LDLAND(JLON)) THEN
            PLIGH_TOT(JLON) = 4.5E-03_JPRB*(PICE(JLON)*PCDEPTH(JLON))**2.0_JPRB*PWMFU(JLON)**0.5_JPRB
          ELSE
            PLIGH_TOT(JLON) = 6.6E-04_JPRB*(PICE(JLON)*PCDEPTH(JLON))**2.0_JPRB*PWMFU(JLON)**0.5_JPRB
          END IF
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 4) THEN
        ! TM5 Meijer 2001.
        IF (PCTOPH(JLON) > 5.0_JPRB .and. IFREEZ25 > 0) THEN
          IF (LDLAND(JLON)) THEN
            PLIGH_TOT(JLON) = 5.0_JPRB*4.0E6_JPRB*PFPLCL(JLON, KLEV)*1.0E-9_JPRB*YDCST%RDAYI
          ELSE
            PLIGH_TOT(JLON) = 0.1_JPRB*(5.0_JPRB*4.0E6_JPRB*PFPLCL(JLON, KLEV)*1.0E-9_JPRB*YDCST%RDAYI)
          END IF
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 5) THEN
        ! Price and Rind 1993 as coded in mozart/ Base/mo_hook.F90
        ! Flash rate in flashes / min -> convert to /(day*km2)
        IF (IFREEZ > 0) THEN
          IF (LDLAND(JLON)) THEN
            PLIGH_TOT(JLON) = 3.44E-5_JPRB*PCTOPH(JLON)**4.9_JPRB*(24.0_JPRB*60.0_JPRB / ZAREA)
          ELSE
            PLIGH_TOT(JLON) = 6.4E-4_JPRB*PCTOPH(JLON)**1.73_JPRB*(24.0_JPRB*60.0_JPRB / ZAREA)
          END IF
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 6) THEN
        ! New lightning parameterization (Ph. Lopez, 2015; flash density in flashes/km2/day).
        IF (ZCDEPTH > 4.0_JPRB .and. PCAPE(JLON) > 175._JPRB .and. ZCBASEH > 0.1_JPRB .and. KCTOP(JLON) < IFREEZ15 .and.  &
        & PCDEPTH(JLON) > 1.0_JPRB) THEN
          PLIGH_TOT(JLON) =  &
          & 36.0_JPRB*MIN(1.0_JPRB, (SQRT(ZAREA) / 30.0_JPRB)**0.15_JPRB)*PCHARGE(JLON)*SQRT(PCAPE(JLON))*MIN(ZCBASEH, ZHBMAX)**2
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 7) THEN
        ! Grewe et al. (2001).
        ! Flash rate in flashes/min  -> convert to flashes/(day*km2).
        IF (PWMFU(JLON) > 0.0_JPRB) THEN
          PLIGH_TOT(JLON) = 1.54E-5_JPRB*(PWMFU(JLON)*SQRT(PCTOPH(JLON)*1000.0_JPRB))**4.9_JPRB*(1440.0_JPRB / ZAREA)
        END IF
        
      ELSE IF (YDEPHY%NLIMODE == 8) THEN
        ! Allen and Pickering (2002).
        ! Note: This parameterization actually computes cloud-to-ground lightning.
        IF (ZPREC_CV > 7.0_JPRB) THEN
          IF (LDLAND(JLON)) THEN
            PLIGH_CTG(JLON) = 3.75E-2_JPRB - 4.76E-2_JPRB*ZPREC_CV + 5.41E-3_JPRB*ZPREC_CV**2 + 3.21E-4_JPRB*ZPREC_CV**3 -  &
            & 2.93E-6_JPRB*ZPREC_CV**4
          ELSE
            PLIGH_CTG(JLON) = 5.23E-2_JPRB - 4.80E-2_JPRB*ZPREC_CV + 5.45E-3_JPRB*ZPREC_CV**2 + 3.68E-5_JPRB*ZPREC_CV**3 -  &
            & 2.42E-7_JPRB*ZPREC_CV**4
          END IF
          PLIGH_CTG(JLON) = MAX(PLIGH_CTG(JLON), 0.0_JPRB)
          ! Cloud-to-ground flash rate in flashes/min  -> convert to flashes/(day*km2).
          PLIGH_CTG(JLON) = PLIGH_CTG(JLON)*(1440.0_JPRB / ZAREA_REF)
          ! Convert to total lightning using Price and Rind (1994).
          ZCOLDDPT = MAX(PCDEPTH(JLON), 5.5_JPRB)
          IF (ZCOLDDPT < 14.0_JPRB) THEN
            ZPCTG = 1.0_JPRB / (0.021_JPRB*(ZCOLDDPT**4) - 0.648_JPRB*(ZCOLDDPT**3) + 7.49_JPRB*(ZCOLDDPT**2) -  &
            & 36.54_JPRB*ZCOLDDPT + 64.09_JPRB)
            ZPCTG = MAX(ZPCTG, 0.0_JPRB)
          ELSE
            ZPCTG = 0.02_JPRB
          END IF
          ! Total flash rate in flashes/(day*km2).
          PLIGH_TOT(JLON) = PLIGH_CTG(JLON) / ZPCTG
        END IF
      END IF
      
      ! Cloud-to-ground lightning frequency as a function of depth of cloud
      ! cold sector (if higher than 5.5 km and less than 14 km).
      ! Note: if using Allen and Pickering (2002), CTG should already be available.
      
      IF (YDEPHY%NLIMODE /= 8) THEN
        IF (PCDEPTH(JLON) > 5.5_JPRB .and. PCDEPTH(JLON) < 14.0_JPRB) THEN
          ZPCTG = 1.0_JPRB / (0.021_JPRB*(PCDEPTH(JLON)**4.0_JPRB) - 0.648_JPRB*(PCDEPTH(JLON)**3.0_JPRB) +  &
          & 7.49_JPRB*(PCDEPTH(JLON)**2.0_JPRB) - 36.54_JPRB*PCDEPTH(JLON) + 64.09_JPRB)
          PLIGH_CTG(JLON) = ZPCTG*PLIGH_TOT(JLON)
        ELSE
          PLIGH_CTG(JLON) = PLIGH_TOT(JLON)
        END IF
      END IF
      
    END IF
  END IF
  
  
END SUBROUTINE CULIGHT_OPENACC
