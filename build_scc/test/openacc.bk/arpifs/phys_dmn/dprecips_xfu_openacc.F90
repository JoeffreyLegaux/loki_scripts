SUBROUTINE DPRECIPS_XFU_OPENACC (KIDIA, KFDIA, KLON, KDTPREC, KSTATS, PDPRECIPS, PXPTYPE, LDRESET, YDSTACK)
  
  !**** *DPRECIPS_XFU*   -  Compute precipitation type diagnostic
  
  !     Purpose.
  !     --------
  !           Compute precipitation type diagnostics for fullpos
  
  !**   Interface.
  !     ----------
  !        *CALL* *DPRECIPS_XFU(...)
  
  !        Explicit arguments :
  !        --------------------
  !----
  ! 0D :
  !----
  ! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
  ! KLON : HORIZONTAL DIMENSION                  (ILONMNH IN *APL_AROME*)
  ! KSTATS : 0 if PTYPE FREQUENT, 2 if PTYPE SEVERE
  ! ------
  ! INOUT :
  ! ------
  ! PDPRECIPS   : precipitation type diagnostic :
  !    0: no precipitation
  !    1: rain   / pluie
  !    3: freezing rain / pluie verglacante
  !    5: dry snow / neige seche
  !    6: wet snow / neige humide
  !    7: rain now mixture / pluie et neige melees
  !    8: ice pellets/ granules de glace
  !    9: graupel   / gresil ou petite grele
  !   10: hail      / grele
  !   11: drizzle/ bruine
  !   12: freezing drizzle / bruine verglacante
  !  193: moist snow / neige mouillee
  !  201: Pluie intermittente
  !  205: Neige sèche intermittente
  !  206: Neige humide intermittente
  !  207: Pluie et neige mêlées intermittentes
  !  213: Neige mouillée intermittente
  
  !
  !        Implicit arguments :
  !        --------------------
  !        COMMON YOMDPRECIPS
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  
  !     Reference.
  !     ----------
  !        Documentation ARPEGE/AROME
  
  !     Authors.
  !     -------
  !      I.Etchevers Y. Seity.
  !      Original : 2018-09-14
  
  !     Modifications.
  !     --------------
  !      R. El Khatib 06-Mar-2023 Enhance vectorization
  !     ------------------------------------------------------------------
  
!$acc routine( DPRECIPS_XFU_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KDTPREC
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTATS
  REAL(KIND=JPRB), INTENT(IN) :: PDPRECIPS(KLON, KDTPREC)
  REAL(KIND=JPRB), INTENT(INOUT) :: PXPTYPE(KLON)
  LOGICAL, INTENT(IN) :: LDRESET
  
  REAL(KIND=JPRB) :: ZTYPE
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: JPREC
  INTEGER(KIND=JPIM) :: JTYPE
  INTEGER(KIND=JPIM) :: INB
  INTEGER(KIND=JPIM) :: IDPRECIPS(11)
  INTEGER(KIND=JPIM) :: INDTOT
  INTEGER(KIND=JPIM) :: ITYPE
  INTEGER(KIND=JPIM), PARAMETER :: JPTYPES = 11
  INTEGER(KIND=JPIM), PARAMETER :: INDPTYPE(JPTYPES) = (/ 11, 1, 7, 8, 9, 5, 6, 193, 10, 12, 3 /)
  !ordering important : from less dangerous to most dangerous
  !INDPTYPE(1)=11
  !INDPTYPE(2)=1
  !INDPTYPE(3)=7
  !INDPTYPE(4)=8
  !INDPTYPE(5)=9
  !INDPTYPE(6)=5
  !INDPTYPE(7)=6
  !INDPTYPE(8)=193
  !INDPTYPE(9)=10
  !INDPTYPE(10)=12
  !INDPTYPE(11)=3
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  JLON = KIDIA
  YLSTACK = YDSTACK
  
  !     ------------------------------------------------------------------
  
  IF (LDRESET) THEN
    PXPTYPE(JLON) = 0._JPRB
  END IF
  
  
  ! Number of ptype occurences
  DO JTYPE=1,JPTYPES
    ZTYPE = REAL(INDPTYPE(JTYPE), kind=JPRB)
    IDPRECIPS(JTYPE) = 0_JPIM
    DO JPREC=1,KDTPREC
      IF (PDPRECIPS(JLON, JPREC) == ZTYPE) THEN
        IDPRECIPS(JTYPE) = IDPRECIPS(JTYPE) + 1_JPIM
      END IF
    END DO
  END DO
  
  INDTOT = SUM(IDPRECIPS(:))
  
  INB = 0
  
  ITYPE = 0_JPIM
  IF (KSTATS == 0) THEN
    ! Most Frequent ptype
    DO JTYPE=1,JPTYPES
      IF (IDPRECIPS(JTYPE) >= INB) THEN
        ITYPE = INDPTYPE(JTYPE)
        INB = IDPRECIPS(JTYPE)
      END IF
    END DO
  ELSE IF (KSTATS == 2) THEN
    ! Most severe ptype
    DO JTYPE=1,JPTYPES
      IF (IDPRECIPS(JTYPE) > 0_JPIM) THEN
        ITYPE = INDPTYPE(JTYPE)
      END IF
    END DO
  END IF
  
  
  !set to 0 when only 1/12 NDTPERIOD is concerned by precipitations (except for hail)
  IF (INDTOT < INT(REAL(KDTPREC, kind=JPRB) / 12._JPRB) .and. ITYPE /= 10_JPIM) THEN
    ITYPE = 0_JPIM
  END IF
  
  ! Intermittent character
  IF (INDTOT < INT(REAL(KDTPREC, kind=JPRB)*5._JPRB / 6._JPRB)) THEN
    IF (ITYPE == 1_JPIM) ITYPE = 201_JPRB
    IF (ITYPE == 7_JPIM) ITYPE = 207_JPRB
    IF (ITYPE == 6_JPIM) ITYPE = 206_JPRB
    IF (ITYPE == 5_JPIM) ITYPE = 205_JPRB
    IF (ITYPE == 193_JPIM) ITYPE = 213_JPRB
    IF (ITYPE == 11_JPIM) ITYPE = 201_JPRB
    ! to avoid drizzle (=very light precipitation) around thunderstorm cells
  END IF
  
  PXPTYPE(JLON) = REAL(ITYPE, kind=JPRB)
  
  
END SUBROUTINE DPRECIPS_XFU_OPENACC
