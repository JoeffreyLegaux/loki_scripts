SUBROUTINE ACNPART_OPENACC (YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, PAPHI, PAPHIF, PAPRSF, PDECRD, PNEB, PCLCH,  &
& PCLCM, PCLCL, PCLCT, PCLCT_RAD, PCLCC, PNEBC, PTOPC, YDSTACK)
  ! optional arguments (convective cloud cover)
  
  ! Purpose:
  ! --------
  !   ACNPART - computes high/medium/low, convective and total cloud cover.
  !   Several overlap options are implemented.
  
  ! Interface:
  ! ----------
  ! INPUT:
  !   KIDIA     - initial index for horizontal loops
  !   KFDIA     - final index for horizontal loops
  !   KLON      - horizontal dimension of arrays
  !   KTDIA     - initial index for vertical loops (usually 1)
  !   KLEV      - vertical dimension of full level arrays
  !   PAPHI     - half level geopotential
  !   PAPHIF    - full level geopotential
  !   PAPRSF    - full level pressure
  !   PDECRD    - decorrelation depth for cloud overlaps [Pa]
  !   PNEB      - total cloud cover on levels (protected from 0 and 1)
  
  ! OUTPUT:
  !   PCLCH     - high cloud cover
  !   PCLCM     - medium cloud cover
  !   PCLCL     - low cloud cover
  !   PCLCT     - total cloud cover
  !   PCLCT_RAD - total cloud cover for radiation
  
  ! INPUT, OPTIONAL:
  !   PNEBC     - convective cloud cover on levels (protected from 0 and 1,
  !               missing in AROME)
  
  ! OUTPUT, OPTIONAL:
  !   PCLCC     - convective cloud cover (missing in AROME)
  !   PTOPC     - TOP of convective cloud  [Pa] (missing in AROME)
  
  
  ! Externals:
  ! ----------
  
  ! Method:
  ! -------
  
  ! Reference:
  ! ----------
  
  ! Author:
  ! -------
  !   2007-02, R. Brozkova
  
  ! Modifications:
  ! --------------
  !   2009-03, C. Wittmann
  !   Introduction of LACPANMX and WMXOV.
  !
  !   2009-07, K. Yessad
  !   Remove CDLOCK + some cleaning.
  !
  !   2009-10, L. Bengtsson
  !   Introduction of LWMOCLOUD.
  !
  !   2016-04, J. Masek
  !   Introduction of LRNUEXP (exponential-random overlap), fix for LWMOCLOUD,
  !   modularization, reordering of arguments, PCLCC optional (missing in AROME).
  !
  !   2016-09, J. Masek
  !   Introduction of radiative cloud cover PCLCT_RAD, needed for consistent
  !   calculation of sunshine duration.
  !
  !   2018-09, J. Masek
  !   Fix of convective cloud cover when WMXOV or RDECRDRED differ from 1.
  !
  !   2018-07, O. Jaron
  !   Introduction of CONV_BASE_TOP to diagnose base and top of convective
  !   clouds. (Base desactivated)
  !
  ! End Modifications
  !-------------------------------------------------------------------------------
  
!$acc routine( ACNPART_OPENACC ) seq
  
  USE MODEL_PHYSICS_MF_MOD, ONLY: MODEL_PHYSICS_MF_TYPE
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMCST, ONLY: TCST
  
  !-----------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(MODEL_PHYSICS_MF_TYPE), INTENT(IN) :: YDML_PHY_MF
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  REAL(KIND=JPRB), INTENT(IN) :: PAPHI(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PAPRSF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDECRD
  REAL(KIND=JPRB), INTENT(IN) :: PNEB(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(OUT) :: PCLCH
  REAL(KIND=JPRB), INTENT(OUT) :: PCLCM
  REAL(KIND=JPRB), INTENT(OUT) :: PCLCL
  REAL(KIND=JPRB), INTENT(OUT) :: PCLCT
  REAL(KIND=JPRB), INTENT(OUT) :: PCLCT_RAD
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PCLCC
  REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PNEBC(KLON, KLEV)
  REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PTOPC
  
#include "acnpart_cloud_cover_wmo.intfb.h"
#include "acnpart_cloud_cover.intfb.h"
#include "acnpart_conv_base_top.intfb.h"
  
  !-----------------------------------------------------------------------
  
  REAL(KIND=JPRB) :: ZDECRDRED
  REAL(KIND=JPRB) :: ZWMXOV
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=jpim) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  ! settings for diagnostic cloud cover
  ZDECRDRED = YDML_PHY_MF%YRPHY0%RDECRDRED
  ZWMXOV = YDML_PHY_MF%YRPHY0%WMXOV
  
  IF (YDML_PHY_MF%YRPHY2%LWMOCLOUD) THEN
    
    ! high/medium/low cloud cover according to WMO heights
    CALL ACNPART_CLOUD_COVER_WMO_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV, PDECRD, PNEB,  &
    & PAPHI, PAPRSF, PCLCH, PCLCM, PCLCL, YDSTACK=YLSTACK)
    
  ELSE
    
    ! high/medium/low cloud cover according to fixed model levels
    CALL ACNPART_CLOUD_COVER_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV, KTDIA,  &
    & YDML_PHY_MF%YRPHY2%NTSHM, PDECRD, PNEB, PAPRSF, PCLCH, YDSTACK=YLSTACK)
    ! high
    CALL ACNPART_CLOUD_COVER_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV,  &
    & YDML_PHY_MF%YRPHY2%NTSHM + 1, YDML_PHY_MF%YRPHY2%NTSML, PDECRD, PNEB, PAPRSF, PCLCM, YDSTACK=YLSTACK)
    ! medium
    CALL ACNPART_CLOUD_COVER_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV,  &
    & YDML_PHY_MF%YRPHY2%NTSML + 1, KLEV, PDECRD, PNEB, PAPRSF, PCLCL, YDSTACK=YLSTACK)
    ! low
    
  END IF
  
  ! total cloud cover
  CALL ACNPART_CLOUD_COVER_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV, KTDIA, KLEV, PDECRD,  &
  & PNEB, PAPRSF, PCLCT, YDSTACK=YLSTACK)
  
  ! convective cloud cover
  IF (PRESENT(PCLCC)) THEN
    CALL ACNPART_CLOUD_COVER_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV, KTDIA, KLEV,  &
    & PDECRD, PNEBC, PAPRSF, PCLCC, YDSTACK=YLSTACK)
  END IF
  
  ! total cloud cover for radiation
  IF (YDML_PHY_MF%YRPHY%LRNUMX .and. (YDML_PHY_MF%YRPHY%LACPANMX .or. YDML_PHY_MF%YRPHY%LRNUEXP)) THEN
    ZDECRDRED = 1._JPRB      ! do not reduce decorrelation depth
    ZWMXOV = 1._JPRB      ! ignore LACPANMX
    CALL ACNPART_CLOUD_COVER_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, ZDECRDRED, ZWMXOV, KTDIA, KLEV,  &
    & PDECRD, PNEB, PAPRSF, PCLCT_RAD, YDSTACK=YLSTACK)
  ELSE
    PCLCT_RAD = PCLCT
  END IF
  
  ! convective top
  IF (PRESENT(PTOPC) .and. YDML_PHY_MF%YRPHY%LPTOPC) THEN
    CALL ACNPART_CONV_BASE_TOP_OPENACC(YDCST, YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV, PNEBC, PAPRSF, PAPHIF, PTOPC,  &
    & YDSTACK=YLSTACK)
  END IF
  
  
END SUBROUTINE ACNPART_OPENACC
