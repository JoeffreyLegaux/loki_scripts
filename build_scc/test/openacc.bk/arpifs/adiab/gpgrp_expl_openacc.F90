SUBROUTINE GPGRP_EXPL_OPENACC (YDGEOMETRY, KST, KEND, PRT, PRTL, PRTM, PREL, PREM, PRTGR, PALPH, PHIHL, PHIHM, PHIFL, PHIFM,  &
& PSGRTL, PSGRTM, PALPHPLL_DER, PALPHPLM_DER, YDSTACK)
  
  !**** *GPGRP_EXPL* - Computation of the pressure gradient force term used in the
  !               RHS of the horizontal wind equation in hydrostatic model.
  !               This term is computed at full levels.
  
  !     Purpose.
  !     --------
  
  !        The pressure gradient force term writes:
  !          - grad[gz] - RT grad[log(prehyds)]
  
  !**   Interface.
  !     ----------
  !        *CALL* *GPGRP_EXPL(...)
  
  !        Explicit arguments :
  !        --------------------
  !         * INPUT:
  !           KST         : start of work
  !           KEND        : working length
  !           PRT         : "R(air)*temperature" at full levels
  !           PRTL        : zonal component of "grad (RT)" at full levels
  !           PRTM        : meridian component of "grad (RT)" at full levels
  !           PREL        : zonal component of "grad (prehyds)"
  !           PREM        : meridian component of "grad (prehyds)"
  !           PXYB        : contains pressure depth, "delta", "alpha".
  !           PXYBDER     : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
  !           PHIHL       : zonal component of "grad[gz]" at half levels
  !           PHIHM       : merid component of "grad[gz]" at half levels
  !           PHIFL       : zonal component of "grad[gz]" at full levels
  !           PHIFM       : merid component of "grad[gz]" at full levels
  
  !         * OUTPUT:
  !           PSGRTL      : zonal comp. of pressure gradient force
  !           PSGRTM      : merid comp. of pressure gradient force
  
  !        Implicit arguments :   None.
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.    None.
  !     ----------
  
  !     Reference.
  !     ----------
  !        See  documentation Arpege ALGORITHME CHAPITRE 6 paragraphe 6
  
  !     Author.
  !     -------
  !      Florence Rabier
  !      Original : 91-11-01
  
  !     Modifications.
  !     --------------
  !   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
  !   K. Yessad (Dec 2008): remove dummy CDLOCK
  !   K. Yessad (March 2009): correct false comments for LRWSDLG=T
  !   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
  !   N. Wedi (Nov 2011) : cleaning
  !   K. Yessad (Mar 2017): move NHEE code in new GNHEE_GRP.
  !   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
  !   K. Yessad (Feb 2018): remove deep-layer formulations.
  !   H Petithomme (Dec 2020): merge VFD loops
  !     ------------------------------------------------------------------
  
!$acc routine( GPGRP_EXPL_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, JPHOOK, DR_HOOK
  
  USE GEOMETRY_MOD, ONLY: GEOMETRY
  
  
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  REAL(KIND=JPRB), INTENT(IN) :: PRT(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PRTL(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PRTM(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PREL(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PREM(YDGEOMETRY%YRDIM%NPROMA)
  REAL(KIND=JPRB), INTENT(IN) :: PRTGR(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PALPH(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PHIHL(YDGEOMETRY%YRDIM%NPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PHIHM(YDGEOMETRY%YRDIM%NPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PHIFL(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PHIFM(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(OUT) :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PALPHPLL_DER(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  REAL(KIND=JPRB), INTENT(IN) :: PALPHPLM_DER(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
  
  !     ------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JROF
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KST
  YLSTACK = YDSTACK
  
  !     ------------------------------------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  !*       1.    COMPUTATION OF THE PRESSURE GRADIENT TERM.
  !              ------------------------------------------
  
  
  ! * The pressure gradient term is directly computed at full levels in
  !   this case and no quantity at half layers is computed.
  
  ! * First transfer " grad[gz] + RT [grad (log(prehyd))] " in (PSGRTL,PSGRTM).
  DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
    PSGRTL(JLON, JLEV) = PHIFL(JLON, JLEV) + PRT(JLON, JLEV)*PRTGR(JLON, JLEV)*PREL(JLON)
    PSGRTM(JLON, JLEV) = PHIFM(JLON, JLEV) + PRT(JLON, JLEV)*PRTGR(JLON, JLEV)*PREM(JLON)
  END DO
  
  
  !     ------------------------------------------------------------------
  
  
  
END SUBROUTINE GPGRP_EXPL_OPENACC
