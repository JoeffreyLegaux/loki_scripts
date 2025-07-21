SUBROUTINE CPPHINP_OPENACC (YDGEOMETRY, YDMODEL, KIDIA, KFDIA, PGEMU, PGELAM, PUT0, PVT0, PQT0, PQT0L, PQT0M, PQSLT0L, PQSLT0M,  &
& PRDELP0, PEVEL0, PCVGQSL, PMU0, PMU0LU, PMU0M, PMU0N, PCVGQ, YDSTACK)
  
  !**** *CPPHINP*  - ComPute PHysical INPut.
  
  !     Purpose.
  !     --------
  
  !**   Interface.
  !     ----------
  !        *CALL* *CPPHINP*
  
  !        Explicit arguments :
  !        --------------------
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  
  !     Reference.
  !     ----------
  
  !     Author.
  !     -------
  
  !     Modifications.
  !     --------------
  !      2002-05, K. YESSAD: consistency with finite element scheme.
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      2004-11, Y. Seity: Do not compute Convergence of Humidity for Arome
  !      2005-10  Y. Bouteloup : Modification for computation of CVGQ
  !        2011-03  Y. Seity: add PMU0N from ARPEGE-climat
  !        Y.Bouteloup 21-Jan-2011 : Compute pmu0m as average of pmu0 instead of
  !                                  pmu0(mean time)
  !      2011-05-10 E. Bazile: PMU0M=PMU0 for NRADFR=1
  !     ------------------------------------------------------------------
  
!$acc routine( CPPHINP_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE TYPE_MODEL, ONLY: MODEL
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(MODEL), INTENT(IN) :: YDMODEL
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(IN) :: PGEMU(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PGELAM(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PUT0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PVT0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PQT0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PQT0L(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PQT0M(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PQSLT0L(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PQSLT0M(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PRDELP0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PEVEL0(YDGEOMETRY%YRDIM%NPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PCVGQSL(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PMU0(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PMU0LU(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PMU0M(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PMU0N(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PCVGQ(YDGEOMETRY%YRDIM%NPROMM, YDGEOMETRY%YRDIMV%NFLEVG)
  
  
  !     ------------------------------------------------------------------
  
#include "abor1.intfb.h"
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JLON
  
  REAL(KIND=JPRB) :: ZED
  REAL(KIND=JPRB) :: ZEPS
  REAL(KIND=JPRB) :: ZABSTSPHY
  REAL(KIND=JPRB) :: ZUSTSPHY
  REAL(KIND=JPRB) :: Z1MU0M
  REAL(KIND=JPRB) :: Z2MU0M
  REAL(KIND=JPRB) :: ZCOS
  REAL(KIND=JPRB) :: ZSIN
  
  LOGICAL :: LLCOMPUTE_CVGQ
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !     ------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  
  !*       1.1   Astronomy.
  
  ZCOS = COS(PGELAM(JLON))
  ZSIN = SIN(PGELAM(JLON))
  
  !DEC$ IVDEP
  PMU0(JLON) = MIN(1.0_JPRB, MAX(YDMODEL%YRML_GCONF%YRRIP%RSIDEC*PGEMU(JLON) -  &
  & YDMODEL%YRML_GCONF%YRRIP%RCODEC*YDMODEL%YRML_GCONF%YRRIP%RCOVSR*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZCOS +  &
  & YDMODEL%YRML_GCONF%YRRIP%RCODEC*YDMODEL%YRML_GCONF%YRRIP%RSIVSR*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZSIN, 0.0_JPRB))
  
  IF (YDMODEL%YRML_PHY_MF%YRARPHY%LMSE) THEN
    !DEC$ IVDEP
    PMU0N(JLON) = MIN(1.0_JPRB, MAX(YDMODEL%YRML_GCONF%YRRIP%RSIDECN*PGEMU(JLON) -  &
    & YDMODEL%YRML_GCONF%YRRIP%RCODECN*YDMODEL%YRML_GCONF%YRRIP%RCOVSRN*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZCOS +  &
    & YDMODEL%YRML_GCONF%YRRIP%RCODECN*YDMODEL%YRML_GCONF%YRRIP%RSIVSRN*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZSIN, 0.0_JPRB))
  ELSE
    PMU0N(JLON) = 0.0
  END IF
  
  
  !*       Lunar astronomy.
  
  IF (YDMODEL%YRML_PHY_MF%YRPHY%LRAYLU) THEN
    !DEC$ IVDEP
    PMU0LU(JLON) = MAX(YDMODEL%YRML_GCONF%YRRIP%RSIDECLU*PGEMU(JLON) -  &
    & YDMODEL%YRML_GCONF%YRRIP%RCODECLU*YDMODEL%YRML_GCONF%YRRIP%RCOVSRLU*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZCOS +  &
    & YDMODEL%YRML_GCONF%YRRIP%RCODECLU*YDMODEL%YRML_GCONF%YRRIP%RSIVSRLU*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZSIN, 0.0_JPRB)
  END IF
  
  IF (YDMODEL%YRML_PHY_EC%YREPHY%LEPHYS .or. YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS .and. YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM .and.  &
  & .not.YDMODEL%YRML_PHY_MF%YRPHY%LRMU0M) THEN
    !DEC$ IVDEP
    PMU0M(JLON) = MAX(YDMODEL%YRML_PHY_RAD%YRERIP%RSIDECM*PGEMU(JLON) -  &
    & YDMODEL%YRML_PHY_RAD%YRERIP%RCODECM*YDMODEL%YRML_PHY_RAD%YRERIP%RCOVSRM*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZCOS +  &
    & YDMODEL%YRML_PHY_RAD%YRERIP%RCODECM*YDMODEL%YRML_PHY_RAD%YRERIP%RSIVSRM*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZSIN, 0.0_JPRB)
    IF (YDMODEL%YRML_PHY_RAD%YRERAD%NRADFR == 1) THEN
      PMU0M(JLON) = PMU0(JLON)
    END IF
  ELSE IF (YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS .and. YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM .and. YDMODEL%YRML_PHY_MF%YRPHY%LRMU0M) THEN
    !DEC$ IVDEP
    Z1MU0M = MAX(YDMODEL%YRML_PHY_RAD%YRERIP%RSIDECM*PGEMU(JLON) -  &
    & YDMODEL%YRML_PHY_RAD%YRERIP%RCODECM*YDMODEL%YRML_PHY_RAD%YRERIP%RCOVSRM*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZCOS +  &
    & YDMODEL%YRML_PHY_RAD%YRERIP%RCODECM*YDMODEL%YRML_PHY_RAD%YRERIP%RSIVSRM*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZSIN, 0.0_JPRB)
    Z2MU0M = 0.5_JPRB*(MAX(YDMODEL%YRML_GCONF%YRRIP%RSIDECF*PGEMU(JLON) -  &
    & YDMODEL%YRML_GCONF%YRRIP%RCODECF*YDMODEL%YRML_GCONF%YRRIP%RCOVSRF*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZCOS +  &
    & YDMODEL%YRML_GCONF%YRRIP%RCODECF*YDMODEL%YRML_GCONF%YRRIP%RSIVSRF*SQRT(1.0_JPRB - PGEMU(JLON)**2)*ZSIN, 0.0_JPRB) +  &
    & PMU0(JLON))
    PMU0M(JLON) = MAX(Z1MU0M, Z2MU0M)
    
  ELSE IF (YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS .and. YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM15) THEN
    CALL ABOR1_ACC('RADIATION SCHEME FROM CYCLE 15 NO LONGER AVAILABLE')
  ELSE
    PMU0M(JLON) = 0.0_JPRB
  END IF
  
  !*    1.2   Convergence of humidity and physical wind components.
  
  ! ky: variable to be put later in a module, saying if the CVGQ
  ! calculation is required or not (its definition must be the same
  ! everywhere in the code).
  ! For the time being this calculation is needed if MF physics
  ! (other than AROME physics) is activated.
  LLCOMPUTE_CVGQ =  &
  & (YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS .or. YDMODEL%YRML_PHY_MF%YRSIMPHL%LSIMPH) .and. .not.YDMODEL%YRML_PHY_MF%YRARPHY%LMPA
  
  IF (LLCOMPUTE_CVGQ) THEN
    
    IF (YDMODEL%YRML_DYN%YRDYN%NCOMP_CVGQ == 0 .or. YDMODEL%YRML_DYN%YRDYN%NCOMP_CVGQ == 1) THEN
      
      ! * Eulerian convergence:
      
      ! a/ Horizontal part of the moisture convergence:
      IF (YDMODEL%YRML_DYN%YRDYN%NCOMP_CVGQ == 0) THEN
        DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
          PCVGQ(JLON, JLEV) = -PQT0L(JLON, JLEV)*PUT0(JLON, JLEV) - PQT0M(JLON, JLEV)*PVT0(JLON, JLEV)
        END DO
      ELSE
        DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
          PCVGQ(JLON, JLEV) = -PQSLT0L(JLON, JLEV)*PUT0(JLON, JLEV) - PQSLT0M(JLON, JLEV)*PVT0(JLON, JLEV)
        END DO
      END IF
      
      ! b/ Vertical part of the moisture convergence:
      
      ! * "pevel0" contains full-layer "etadot (d prehyd / d eta)"
      !   The current discretisation of the vertical advection of a variable X
      !   which is proposed is:
      !   [etadot d X/d eta](l) = 0.5 * [1 / Delta eta](l)
      !    * [etadot d prehyd/d eta](l) * (X(l+1)-X(l-1))
      !   Remark K.Y.: this discretisation is not fully consistent with the
      !    discretisation of a vertical advection provided by ECMWF between
      !    CY29 and CY29R2 (first coded in CPDYN, then transferred to CPEULDYN),
      !    and it is desirable (but not urgent) to update this discretisation
      !    in the future to make it consistent with CPEULDYN (in this case
      !    one must use routine VERDER to compute the vertical derivative
      !    [d X/d eta](l)).
      
      ! * Layer 1:
      
      PCVGQ(JLON, 1) = PCVGQ(JLON, 1) + (PQT0(JLON, 1) - PQT0(JLON, 2))*PEVEL0(JLON, 1)*PRDELP0(JLON, 1)
      
      ! * Layers 2 to L-1:
      
      DO JLEV=2,YDGEOMETRY%YRDIMV%NFLEVG - 1
        ZED = 0.5_JPRB*PEVEL0(JLON, JLEV)
        
        PCVGQ(JLON, JLEV) = PCVGQ(JLON, JLEV) + (PQT0(JLON, JLEV - 1) - PQT0(JLON, JLEV + 1))*ZED*PRDELP0(JLON, JLEV)
      END DO
      
      ! * Layer L:
      
      PCVGQ(JLON, YDGEOMETRY%YRDIMV%NFLEVG) = PCVGQ(JLON, YDGEOMETRY%YRDIMV%NFLEVG) + (PQT0(JLON, YDGEOMETRY%YRDIMV%NFLEVG - 1) - &
      &  PQT0(JLON, YDGEOMETRY%YRDIMV%NFLEVG))*PEVEL0(JLON, YDGEOMETRY%YRDIMV%NFLEVG)*PRDELP0(JLON, YDGEOMETRY%YRDIMV%NFLEVG)
      
      
    ELSE
      
      ! * case NCOMP_CVGQ=2: Semi-Lagrangian convergence:
      
      ZEPS = 1.E-12
      ZABSTSPHY = MAX(ZEPS, ABS(YDMODEL%YRML_PHY_MF%YRPHY2%TSPHY))
      IF (YDMODEL%YRML_PHY_MF%YRPHY2%TSPHY >= 0) THEN
        ZUSTSPHY = 1.0_JPRB / ZABSTSPHY
      ELSE
        ZUSTSPHY = -1.0_JPRB / ZABSTSPHY
      END IF
      
      DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
        PCVGQ(JLON, JLEV) = PCVGQSL(JLON, JLEV)*ZUSTSPHY
      END DO
    END IF
    
  END IF
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE CPPHINP_OPENACC
