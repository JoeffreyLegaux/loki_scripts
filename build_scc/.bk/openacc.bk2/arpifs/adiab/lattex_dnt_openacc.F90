SUBROUTINE LATTEX_DNT_OPENACC (KSTEP, YDGEOMETRY, YDLDDH, YDRIP, YDDYN, YDDYNA, KST, KEND, LDSETTLS, KXLAG, PESGP, PESGM, PXT0,  &
& PXT9, PMOY1X, PMIXNL, PXSI, PXNLT9, PXT1, PXL0, PXL9, PXLF9, PCXNLT9, PSIDDHXT1, PSIDDHXT9, PSIDDHXL0, PXLF0, LDNESC, LDVD5,  &
& PDYT0, PDYT9, YDSTACK)
  
  !------------------------------------------------------------------------------
  ! LATTEX_DNT - Semi-Lagrangian scheme.
  !              Computation of the t and t-dt useful quantities
  !              at grid-points. Equations for tri-dimensional
  !              variables for 2TL scheme.
  
  ! Purpose
  ! -------
  
  ! Interface
  ! ---------
  !   CALL LATTEX_DNT(..)
  
  ! Explicit arguments :
  ! --------------------
  
  ! * INPUT:
  !        KST       - first element of work.
  !        KPROF     - depth of work.
  !        LDSETTLS  - .T./.F.: Stable/Conventional extrapolations for SL2TL.
  !        KXLAG     - type of SL discretisation
  !        PESGP     - (1 + uncentering factor).
  !        PESGM     - (1 - uncentering factor).
  !        PXT0      - prognostic variable time t (predictor or SI),
  !                    preliminary t+dt (corrector)
  !        PXT9      - prognostic variable time t (corrector),
  !                    not used for predictor or SI
  !        PMOY1X    - full nonlinear model at time t [(Delta t/2) "cursive" A]
  !        PMIXNL    - extrapolation control variable for mixed NESC/SETTLS scheme
  
  ! * INPUT/OUTPUT:
  !        PXSI      - semi-implicit linear model at time t
  !                    [- (Delta t/2) "cursive" B]
  !        PXNLT9    - buffer used during predictor resp. SI time step
  !        PXT1      - t+dt term and other final point terms
  !        PXL0      - second SL quantity to be interpolated (linearly) at O if KXLAG>3
  !        PXLF0     - third SL quantity to be interpolated (linearly) at O if KXLAG==3
  !        PXL9      - SL quantity to be interpolated at O
  !                    (if NSPLTHOI=1 with diffusive interpolation)
  !        PXLF9     - SL quantity to be interpolated at O with high order
  !                    interpolation (only meaningfull with NSPLTHOI /= 0)
  !        PCXNLT9   - buffer used during corrector
  !        PSIDDHX.. - buffers for DDH (linear terms).
  
  ! * OPTIONAL INPUT:
  !        LDNESC    - cf. LNESC in YOMDYNA.
  !        LDVDWY    - activate extra term in gw equation for g[W]
  !        PDYT0     - in case of  g[W] used as vertical motion variable this
  !                    represents Yterm at t (predictor) resp. t+dt (corrector)
  !        PDYT9     - buffer used during corrector step for Yterm
  ! Externals
  ! ---------
  !   none
  !   Called by LATTEX.
  
  ! Method
  ! ------
  
  ! Reference
  ! ---------
  !   Arpege documentation about semi-Lagrangian scheme.
  
  ! Author
  ! ------
  !      Mar-2002 J. VIVODA (SHMI/CHMI/LACE) - rationalization of LATTEX
  
  ! Modifications
  ! -------------
  !   12-Oct-2002 J. Masek   PC bugfix
  !   01-Oct-2003 M. Hamrud  CY28 Cleaning
  !   09-Jun-2004 J. Masek   NH cleaning (LPC_NOTR, LDFULLIMP)
  !   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !   K. Yessad (Aug 2008): simplify XIDT treatment with PC + cleanings
  !   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
  !   F. Vana  15-Oct-2009: option NSPLTHOI
  !   F. Vana  22-Feb-2011: option K[X]LAG=4
  !   K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
  !   K. Yessad (July 2014): Move some variables, rename some variables.
  !   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
  ! End Modifications
  !------------------------------------------------------------------------------
  
!$acc routine( LATTEX_DNT_OPENACC )
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: DR_HOOK, JPHOOK, LHOOK
  USE YOMDYNA, ONLY: TDYNA
  
  USE YOMDYN, ONLY: TDYN
  USE YOMRIP, ONLY: TRIP
  USE YOMLDDH, ONLY: TLDDH
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(TLDDH), INTENT(IN) :: YDLDDH
  TYPE(TRIP), INTENT(IN) :: YDRIP
  TYPE(TDYN), INTENT(IN) :: YDDYN
  TYPE(TDYNA), INTENT(IN) :: YDDYNA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  LOGICAL, INTENT(IN) :: LDSETTLS
  INTEGER(KIND=JPIM), INTENT(IN) :: KXLAG
  REAL(KIND=JPRB), INTENT(IN) :: PESGP
  REAL(KIND=JPRB), INTENT(IN) :: PESGM
  REAL(KIND=JPRB), INTENT(IN) :: PXT0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PXT9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PMOY1X(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXSI(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXNLT9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXT1(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXL0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXL9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXLF9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(INOUT) :: PCXNLT9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSIDDHXT1(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSIDDHXT9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(INOUT) :: PSIDDHXL0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  REAL(KIND=JPRB), INTENT(OUT) :: PXLF0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
  LOGICAL, OPTIONAL, INTENT(IN) :: LDNESC
  LOGICAL, OPTIONAL, INTENT(IN) :: LDVD5
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PDYT0(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PDYT9(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  
  !     * ZXNLT0 (resp. ZXNLT1) the non linear term at t (resp. t+dt).
  temp (REAL (KIND=JPRB), ZXNLT1, (YDGEOMETRY%YRDIM%NPROMA))
  temp (REAL (KIND=JPRB), ZXNLT0, (YDGEOMETRY%YRDIM%NPROMA))
  REAL(KIND=JPRB) :: ZXIDT0
  REAL(KIND=JPRB) :: ZXIDT9
  REAL(KIND=JPRB) :: ZXIGP
  REAL(KIND=JPRB) :: ZNESC
  REAL(KIND=JPRB) :: ZSETTLS
  REAL(KIND=JPRB) :: ZMIXNL
  
  LOGICAL :: LLCT
  LOGICAL :: LLNESC
  LOGICAL :: LLVD5
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  IF (KIND (ZXNLT1) == 8) THEN
    alloc8 (ZXNLT1)
  ELSE
    IF (KIND (ZXNLT1) == 8) THEN
      alloc4 (ZXNLT1)
    ELSE
      STOP 1
    END IF
  END IF
  IF (KIND (ZXNLT0) == 8) THEN
    alloc8 (ZXNLT0)
  ELSE
    IF (KIND (ZXNLT0) == 8) THEN
      alloc4 (ZXNLT0)
    ELSE
      STOP 1
    END IF
  END IF
  JLON = KST
  
  !     ------------------------------------------------------------------
  !     ------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  
  !*      1. AUXILIARY VARIABLES.
  !       -----------------------
  
  IF (PRESENT(LDNESC)) THEN
    LLNESC = LDNESC
  ELSE
    LLNESC = YDDYNA%LNESC
  END IF
  
  LLCT = YDDYNA%LPC_FULL .and. YDDYN%NCURRENT_ITER > 0
  
  IF (PRESENT(LDVD5)) THEN
    LLVD5 = LDVD5
  ELSE
    LLVD5 = .false.
  END IF
  
  
  ZXIDT0 = 1.0_JPRB + YDDYN%XIDT
  ZXIDT9 = 1.0_JPRB + YDDYN%XIDT
  ZXIGP = 1.0_JPRB + YDDYN%XIDT
  
  !     ------------------------------------------------------------------
  
  !*      2. MAIN CALCULATIONS.
  !       ---------------------
  
  
  !############################################
  ! 2.1 Predictor for LPC_FULL,
  !     or case NSITER=0.
  !############################################
  
  DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
    
    ! nonlinear residual time t
    ZXNLT0(JLON) = PMOY1X(JLON, JLEV) + ZXIDT0*PXSI(JLON, JLEV)
    
    ! Fill PXL9,PXLF9,PXL0,PXT1.
    IF (KXLAG == 2 .and. YDDYN%NSPLTHOI /= 0) THEN
      IF (KSTEP <= YDRIP%NFOST .or. LLNESC) THEN
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT0(JLON, JLEV)
        PXLF9(JLON, JLEV) = PXLF9(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV)
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
      ELSE IF (LDSETTLS) THEN
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT0(JLON, JLEV)
        ZNESC = PESGM*PMOY1X(JLON, JLEV)
        ZSETTLS = PESGM*PMOY1X(JLON, JLEV) + (ZXNLT0(JLON) - PXNLT9(JLON, JLEV))
        ZMIXNL = PMIXNL(JLON, JLEV)
        PXLF9(JLON, JLEV) = PXLF9(JLON, JLEV) + ZMIXNL*ZSETTLS + (1.0_JPRB - ZMIXNL)*ZNESC
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT0(JLON, JLEV)
        PXLF9(JLON, JLEV) = PXLF9(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV) + 0.5_JPRB*PESGM*(ZXNLT0(JLON) - PXNLT9(JLON, JLEV))
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*(1.5_JPRB*ZXNLT0(JLON) - 0.5_JPRB*PXNLT9(JLON, JLEV))
      END IF
    ELSE IF (KXLAG == 2 .and. YDDYN%NSPLTHOI == 0) THEN
      IF (KSTEP <= YDRIP%NFOST .or. LLNESC) THEN
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT0(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV)
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
      ELSE IF (LDSETTLS) THEN
        ZNESC = PXT0(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV)
        ZSETTLS = PXT0(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV) + (ZXNLT0(JLON) - PXNLT9(JLON, JLEV))
        ZMIXNL = PMIXNL(JLON, JLEV)
        PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + ZMIXNL*ZSETTLS + (1.0_JPRB - ZMIXNL)*ZNESC
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        PXL9(JLON, JLEV) =  &
        & PXL9(JLON, JLEV) + PXT0(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV) + 0.5_JPRB*PESGM*(ZXNLT0(JLON) - PXNLT9(JLON, JLEV))
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*(1.5_JPRB*ZXNLT0(JLON) - 0.5_JPRB*PXNLT9(JLON, JLEV))
      END IF
    ELSE IF (KXLAG >= 3) THEN
      IF (KSTEP <= YDRIP%NFOST .or. LLNESC) THEN
        PXL0(JLON, JLEV) = PXL0(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV)
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
      ELSE IF (LDSETTLS) THEN
        IF (YDDYNA%LPC_CHEAP) THEN
          ZNESC = PESGM*PMOY1X(JLON, JLEV)
          ZSETTLS = ZXNLT0(JLON) - PXNLT9(JLON, JLEV)
          ZMIXNL = PMIXNL(JLON, JLEV)
          PXL0(JLON, JLEV) = PXL0(JLON, JLEV) + ZNESC
          PXLF0(JLON, JLEV) = ZMIXNL*ZSETTLS
          PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
        ELSE
          ZNESC = PESGM*PMOY1X(JLON, JLEV)
          ZSETTLS = ZXNLT0(JLON) - PXNLT9(JLON, JLEV)
          ZMIXNL = PMIXNL(JLON, JLEV)
          PXL0(JLON, JLEV) = PXL0(JLON, JLEV) + ZNESC + ZMIXNL*ZSETTLS
          PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*ZXNLT0(JLON)
        END IF
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        PXL0(JLON, JLEV) = PXL0(JLON, JLEV) + PESGM*PMOY1X(JLON, JLEV) + 0.5_JPRB*PESGM*(ZXNLT0(JLON) - PXNLT9(JLON, JLEV))
        PXT1(JLON, JLEV) = PXT1(JLON, JLEV) + PESGP*(1.5_JPRB*ZXNLT0(JLON) - 0.5_JPRB*PXNLT9(JLON, JLEV))
      END IF
      PXL9(JLON, JLEV) = PXL9(JLON, JLEV) + PXT0(JLON, JLEV)
    END IF
    
    IF (LLVD5) THEN
      PXL9(JLON, JLEV) = PXL9(JLON, JLEV) - PDYT0(JLON, JLEV)
    END IF
    
    ! save quantities for corrector step
    IF (YDDYNA%LPC_FULL) THEN
      ! save nonlinear model at time t
      PCXNLT9(JLON, JLEV) = PMOY1X(JLON, JLEV)
    END IF
    
    IF (.not.LLNESC) THEN
      ! save of nonlinear residual at time t
      ! to be used as nonlinear residual at time t-dt next time step
      PXNLT9(JLON, JLEV) = PMOY1X(JLON, JLEV) + ZXIDT9*PXSI(JLON, JLEV)
    END IF
    
  END DO
  
  
  !########################################################
  ! 2.3 Addition of preliminary quantity for LAGPHY physics
  !########################################################
  
  DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
    IF (YDDYN%XIDT > 0.0_JPRB) THEN
      PXT1(JLON, JLEV) = PXT1(JLON, JLEV) - ZXIGP*PXSI(JLON, JLEV)
      PXSI(JLON, JLEV) = ZXIGP*PXSI(JLON, JLEV)
    ELSE
      PXT1(JLON, JLEV) = PXT1(JLON, JLEV) - PESGP*PXSI(JLON, JLEV)
      PXSI(JLON, JLEV) = PESGP*PXSI(JLON, JLEV)
    END IF
  END DO
  
  !########################################################
  ! 2.4  DDH computations for SI correction
  !########################################################
  
  IF (YDLDDH%LRSIDDH) THEN
    IF (KXLAG >= 3) THEN
      IF (KSTEP <= YDRIP%NFOST .or. LLNESC) THEN
        DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
          PSIDDHXT1(JLON, JLEV) = PSIDDHXT1(JLON, JLEV) + PESGP*ZXIDT0*PXSI(JLON, JLEV)
        END DO
      ELSE IF (LDSETTLS) THEN
        DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
          PSIDDHXL0(JLON, JLEV) = PSIDDHXL0(JLON, JLEV) + (ZXIDT0*PXSI(JLON, JLEV) - PSIDDHXT9(JLON, JLEV))
          PSIDDHXT1(JLON, JLEV) = PSIDDHXT1(JLON, JLEV) + PESGP*ZXIDT0*PXSI(JLON, JLEV)
        END DO
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
          PSIDDHXL0(JLON, JLEV) = PSIDDHXL0(JLON, JLEV) + 0.5_JPRB*PESGM*(ZXIDT0*PXSI(JLON, JLEV) - PSIDDHXT9(JLON, JLEV))
          PSIDDHXT1(JLON, JLEV) =  &
          & PSIDDHXT1(JLON, JLEV) + PESGP*(1.5_JPRB*ZXIDT0*PXSI(JLON, JLEV) - 0.5_JPRB*PSIDDHXT9(JLON, JLEV))
        END DO
      END IF
    END IF
    IF (.not.LLNESC) THEN
      ! save of semi-implicit linear term  at time t
      ! to be used at time t-dt next time step
      DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
        PSIDDHXT9(JLON, JLEV) = ZXIDT9*PXSI(JLON, JLEV)
      END DO
    END IF
  END IF
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE LATTEX_DNT_OPENACC
