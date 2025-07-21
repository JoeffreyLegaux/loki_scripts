SUBROUTINE LATTEX_TNT_OPENACC (YDGEOMETRY, YDLDDH, YDDYN, KST, KEND, KXLAG, PESGP, PESGM, PXT9, PMOY1X, PXSI9, PXSI0, PXT1,  &
& PXL0, PXL9, PXLF9, PSIDDHXT1, PSIDDHXL0, YDSTACK)
  
  !**** *LATTEX_TNT*   Semi-Lagrangian scheme.
  !                Computation of the t and t-dt useful quantities
  !                 at grid-points. Equations for tri-dimensional
  !                 variables for 3TL scheme
  
  !     Purpose.
  !     --------
  
  !**   Interface.
  !     ----------
  !        *CALL* *LATTEX_TNT(..)
  
  !        Explicit arguments :
  !        --------------------
  
  !        INPUT:
  !          KST     - first element of work.
  !          KPROF   - depth of work.
  !          KXLAG   - type of SL discretisation.
  !          PESGP   - (1 + uncentering factor).
  !          PESGM   - (1 - uncentering factor).
  !          PXT9    - prognostic variable time t-dt
  !          PMOY1X  - full nonlinear model at time t [(Delta t/2) "cursive" A]
  !          PXSI9   - semi-implicit linear model at time t-dt
  !                    [- (Delta t/2) "cursive" B]
  
  !        INPUT/OUTPUT:
  !          PXSI0      - semi-implicit linear model at time t
  !                       [- (Delta t/2) "cursive" B]
  !          PXT1       - t+dt term and other final point terms
  !          PXL0       - second SL quantity to be interpolated (linearly) at O
  !                       if KXLAG>3
  !          PXL9       - SL quantity to be interpolated at O
  !                       (if NSPLTHOI=1 with diffusive interpolation)
  !          PXLF9      - SL quantity to be interpolated at O with high order
  !                       interpolation (only meaningfull with NSPLTHOI /= 0)
  !          PSIDDHXT1  - semi-implicit  term at final point terms
  !          PSIDDHXL0  - second Semi-implicit linear quantity to be interpolated
  !                       (linearly) at O if KXLAG>3
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  !           none
  !           Called by LATTEX.
  
  !     Reference.
  !     ----------
  !             Arpege documentation about semi-Lagrangian scheme.
  
  !     Author.
  !     -------
  !        J. VIVODA (SHMI/CHMI/LACE) - rationalization of LATTEX
  
  !     Modifications.
  !     --------------
  !     Original : March 2002
  !     M.Hamrud      01-Oct-2003 CY28 Cleaning
  !     K. Yessad Aug 2008: rationalisation of dummy arg interfaces + cleanings
  !     F. Vana   15-Oct-2009 : optin NSPLTHOI
  !     F. Vana   22-Feb-2011 : optin N[X]LAG=4
  !     ------------------------------------------------------------------
  
!$acc routine( LATTEX_TNT_OPENACC )
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMDYN, ONLY: TDYN
  USE YOMLDDH, ONLY: TLDDH
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(TDYN), INTENT(IN) :: YDDYN
  TYPE(TLDDH), INTENT(IN) :: YDLDDH
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  INTEGER(KIND=JPIM), INTENT(IN) :: KXLAG
  REAL(KIND=JPRB), INTENT(IN) :: PESGP
  REAL(KIND=JPRB), INTENT(IN) :: PESGM
  REAL(KIND=JPRB), INTENT(IN) :: PXT9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PMOY1X(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PXSI9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXSI0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXT1(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXL0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXL9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXLF9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSIDDHXT1(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSIDDHXL0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  
  !     ------------------------------------------------------------------
  
  temp (REAL (KIND=JPRB), ZMOYXNL, (YDGEOMETRY%YRDIM%NPROMA))
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZMOYXNL) == 8) THEN
    alloc8 (ZMOYXNL)
  ELSE
    IF (KIND (ZMOYXNL) == 8) THEN
      alloc4 (ZMOYXNL)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KST
  
  !     ------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  
  DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
    
    ZMOYXNL(JLON) = PMOY1X(JLON, JLEV) + PXSI0(JLON, JLEV)
    IF (KXLAG == 2) THEN
      IF (YDDYN%NSPLTHOI /= 0) THEN
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZMOYXNL(JLON) - PESGP*PXSI0(JLON, JLEV)
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT9(JLON, JLEV)
        PXLF9(JLON, JLEV) = PXLF9(JLON, JLEV) + PESGM*ZMOYXNL(JLON) - PESGM*PXSI9(JLON, JLEV)
      ELSE
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZMOYXNL(JLON) - PESGP*PXSI0(JLON, JLEV)
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT9(JLON, JLEV) + PESGM*ZMOYXNL(JLON) - PESGM*PXSI9(JLON, JLEV)
      END IF
    ELSE IF (KXLAG >= 3) THEN
      PXL0(JLON, JLEV) = PXL0(JLON, JLEV) + PESGM*ZMOYXNL(JLON) - PESGM*PXSI9(JLON, JLEV)
      PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZMOYXNL(JLON) - PESGP*PXSI0(JLON, JLEV)
      PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT9(JLON, JLEV)
    END IF
    
    PXSI0(JLON, JLEV) = PESGP*PXSI0(JLON, JLEV)
    
  END DO
  
  !     ------------------------------------------------------------------
  
  IF (YDLDDH%LRSIDDH) THEN
    IF (KXLAG >= 3) THEN
      DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
        PSIDDHXL0(JLON, JLEV) = PSIDDHXL0(JLON, JLEV) + PESGM*PXSI0(JLON, JLEV) - PESGM*PXSI9(JLON, JLEV)
      END DO
    END IF
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      PSIDDHXT1(JLON, JLEV) = PSIDDHXT1(JLON, JLEV) + PESGP*PXSI0(JLON, JLEV)
    END DO
  END IF
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE LATTEX_TNT_OPENACC
