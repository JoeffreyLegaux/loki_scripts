SUBROUTINE LASURE_OPENACC (YDGEOMETRY, YDEPHY, YDDYN, YDEDYN, YDPHY, KST, KEND, PMAPPA, PGM, PBT, PBDT, PREDIV, YDSTACK)
  ! ----- INPUT ---------------------------------------------------------------
  ! ----- OUTPUT --------------------------------------------------------------
  
  !**** *LASURE*   Semi-Lagrangian scheme.
  !                Set-up for other subroutines called by LACDYN.
  
  !     Purpose.
  !     --------
  !       Computes some intermediate quantities necessary in the other
  !       subroutines called by LACDYN.
  
  !**   Interface.
  !     ----------
  !        *CALL* *LASURE(..)
  
  !        Explicit arguments :
  !        --------------------
  
  !        INPUT:
  !          KSTART  : first element of work.
  !          KPROF   : depth of work.
  !          PBETADT : BETADT or 0 according to configuration.
  !          PDT     : time step for the first time-integration step of
  !                    a leap-frog scheme or all time-integration steps of
  !                    a two-time level scheme; 2*time step for the following
  !                    time-integration steps of a leap-frog scheme.
  !          KIBL    : index into YDGSGEOM in YDGEOMETRY
  
  !        OUTPUT:
  !          PDTS2   : 0.5*PDT.
  !          PBT     : PDTS2*PBETADT.
  !          LD2TLFF1: .T./.F.: Refined treatement of (2*Omega Vec r) at
  !                    the origin point when there is t-dt (or t in SL2TL)
  !                    physics / Other cases.
  !          PBDT    : PBT or PBT*(c**2/GM**2) according to LSIDG.
  !          PREDIV  : 1. or c**2/GM**2 according to LSIDG.
  !          PESGP   : (1 + uncentering factor).
  !          PESGM   : (1 - uncentering factor).
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation about semi-lagrangian scheme.
  
  !     Externals.
  !     ----------
  !           none.
  !           Called by LACDYN, LACDYNTL and LACDYNAD.
  
  !     Reference.
  !     ----------
  !             Arpege documentation about semi-lagrangian scheme.
  
  !     Author.
  !     -------
  !      K. YESSAD (METEO FRANCE/CNRM/GMAP) after old part one of LACDYN.
  !      Loops are rewritten according to F90 norms.
  !      Original : JULY 1995.
  
  !     Modifications.
  !     --------------
  !      J.Vivoda 03-2002 PC schemes for NH dynamics (LPC_XXXX keys)
  !      M.Hamrud      01-Oct-2003 CY28 Cleaning
  !      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
  !      K. Yessad (Dec 2011): use YDGSGEOM.
  !      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
  !     ------------------------------------------------------------------
  
!$acc routine( LASURE_OPENACC ) seq
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOEPHY, ONLY: TEPHY
  USE YOMPHY, ONLY: TPHY
  USE YEMDYN, ONLY: TEDYN
  USE YOMDYN, ONLY: TDYN
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(TEPHY), INTENT(IN) :: YDEPHY
  TYPE(TDYN), INTENT(IN) :: YDDYN
  TYPE(TEDYN), INTENT(IN) :: YDEDYN
  TYPE(TPHY), INTENT(IN) :: YDPHY
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PMAPPA(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PGM(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PBT
  REAL(KIND=JPRB), INTENT(OUT) :: PBDT(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PREDIV(YDGEOMETRY%YRDIM%NPROMA)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JROF
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !     ------------------------------------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  !*       1.    PRELIMINARY INITIALISATIONS:
  !              ----------------------------
  
  !     * Time step.
  
  IF (YDDYN%LSIDG) THEN
    PBDT(JLON) = PBT
    PREDIV(JLON) = 1.0_JPRB
  ELSE IF (YDEDYN%LESIDG) THEN
    PBDT(JLON) = PBT*PMAPPA(JLON) / (PGM(JLON)*PGM(JLON))
    PREDIV(JLON) = PMAPPA(JLON) / (PGM(JLON)*PGM(JLON))
  ELSE
    PBDT(JLON) = PBT*YDGEOMETRY%YRGEM%RSTRET*YDGEOMETRY%YRGEM%RSTRET / (PGM(JLON)*PGM(JLON))
    PREDIV(JLON) = YDGEOMETRY%YRGEM%RSTRET*YDGEOMETRY%YRGEM%RSTRET / (PGM(JLON)*PGM(JLON))
  END IF
  
END SUBROUTINE LASURE_OPENACC
