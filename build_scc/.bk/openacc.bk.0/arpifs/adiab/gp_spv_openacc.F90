SUBROUTINE GP_SPV_OPENACC (YDGEOMETRY, YDDYN, YDDYNA, YDSIMPHL, LDTL, KST, KEND, PSPT0, PSPT0L, PSPT0M, PSPT9, PSPT9L, PSPT9M,  &
& PRE0, PRE0L, PRE0M, PRE9, PRE9L, PRE9M, YDSTACK)
  
  !     ------------------------------------------------------------------
  !**** *GP_SPV* -
  
  !     Purpose.
  !     --------
  
  !**   Interface.
  !     ----------
  !        *CALL* *GP_SPV(...)*
  
  !        Explicit arguments :
  !        --------------------
  
  !        abbreviation "prehyds" means "hydrostatic surface pressure".
  
  !        INPUT:
  !        ------
  !        LDTL               : true if gp_spv is called from TL/AD models
  !        KST                : start of work
  !        KEND               : depth of work
  !        PSPT0,PSPT0L,PSPT0M: log(prehyds) and derivatives at t.
  !        PSPT9,PSPT9L,PSPT9M: log(prehyds) and derivatives at t-dt.
  
  !        OUTPUT:
  !        -------
  !        PRE0,PRE0L,PRE0M   : prehyds and derivatives at t.
  !        PRE9,PRE9L,PRE9M   : prehyds and derivatives at t-dt.
  
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
  
  ! Modifications
  ! -------------
  !   07-Aug-2001 R. El Khatib  Pruning options
  !      Mar-2002 J. Vivoda     PC schemes for NH dynamics (LPC_XXXX keys)
  !   27-Jun-2002 C. Fischer    cdlock & ldtl
  !   01-Oct-2003 M. Hamrud     CY28 Cleaning
  !      Dec-2003 K. Yessad     multiplication by GM has moved in GPMPFC_GMVS.
  !   09-Jun-2004 J. Masek      NH cleaning (LFULLIMP)
  !   01-Jul-2004 K. Yessad     Make clearer the tests for PC scheme.
  !   K. Yessad (Dec 2008): remove dummy CDLOCK
  !   K. Yessad (Nov 2012): simplify testings.
  !   K. Yessad (July 2014): Move some variables.
  ! End Modifications
  !------------------------------------------------------------------------------
  
!$acc routine( GP_SPV_OPENACC ) seq
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMDYN, ONLY: TDYN
  USE YOMDYNA, ONLY: TDYNA
  USE YOMSIMPHL, ONLY: TSIMPHL
  
  !------------------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(TDYN), INTENT(IN) :: YDDYN
  TYPE(TDYNA), INTENT(IN) :: YDDYNA
  TYPE(TSIMPHL), INTENT(IN) :: YDSIMPHL
  LOGICAL, INTENT(IN) :: LDTL
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PSPT0(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPT0L(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPT0M(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPT9(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPT9L(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PSPT9M(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PRE0(YDGEOMETRY%YRDIM%NPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PRE0L(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PRE0M(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PRE9(YDGEOMETRY%YRDIM%NPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PRE9L(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(OUT) :: PRE9M(YDGEOMETRY%YRDIM%NPROMA)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JROF
  LOGICAL :: LLPRE9
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !     ------------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  
  IF (.not.((YDSIMPHL%LSIMPH .or. YDDYNA%LNHDYN) .and. LDTL .or. .not.LDTL)) THEN
    LLPRE9 = .false.
  ELSE
    LLPRE9 = .not.YDDYNA%LTWOTL
  END IF
  
  PRE0(JLON, YDGEOMETRY%YRDIMV%NFLEVG) = EXP(PSPT0(JLON))
  PRE0L(JLON) = PSPT0L(JLON)*PRE0(JLON, YDGEOMETRY%YRDIMV%NFLEVG)
  PRE0M(JLON) = PSPT0M(JLON)*PRE0(JLON, YDGEOMETRY%YRDIMV%NFLEVG)
  IF (LLPRE9) THEN
    PRE9(JLON, YDGEOMETRY%YRDIMV%NFLEVG) = EXP(PSPT9(JLON))
    PRE9L(JLON) = PSPT9L(JLON)*PRE9(JLON, YDGEOMETRY%YRDIMV%NFLEVG)
    PRE9M(JLON) = PSPT9M(JLON)*PRE9(JLON, YDGEOMETRY%YRDIMV%NFLEVG)
  END IF
  
  !     ------------------------------------------------------------------
  
  
END SUBROUTINE GP_SPV_OPENACC
