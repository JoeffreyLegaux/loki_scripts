SUBROUTINE ACDAYD_OPENACC (YDCST, YDRIP, KIDIA, KFDIA, KLON, KLEV, KTDIA, KSGST, PGEMU, PMU0, PAPHIF, PDELP, PFRSO, YDSTACK)
  
  !**** *ACDAYD * - Compute DAY Duration, as a function of altitude, and correct solar fluxes accordingly.
  
  !-----------------------------------------------------------------------
  ! -   INPUT ARGUMENTS.
  !     ----------------
  
  ! - PHYSICS DIMENSIONNING PARAMETERS
  
  ! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
  ! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
  ! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
  ! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
  ! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
  ! KSGST      : DIMENSION POUR LE RAYONNEMENT.
  
  ! - PHYSICS VARIABLES.
  
  ! - 2D (0:KLEV) .
  
  ! - 2D (1:KLEV) .
  
  ! PAPHIF     : FULL LEVEL GEOPOTENTIAL
  ! PDELP      : PRESSURE DIFFERENCE OVER THE LAYER
  
  ! - 1D .
  
  ! PGEMU      : sinus de la latitude.
  ! PMU0       : sinus de la hauteur du Soleil au-dessus de l'horizon.
  
  !-----------------------------------------------------------------------
  
  ! -   OUTPUT ARGUMENTS
  !     ----------------
  
  ! - 2D (0:KLEV) .
  
  ! - 2D (1:KLEV) .
  
  ! - 1D (DIAGNOSTIQUE) .
  
  !-----------------------------------------------------------------------
  
  ! -   INPUT/OUTPUT ARGUMENTS
  !     ----------------------
  
  ! - 2D (0:KLEV) .
  
  ! PFRSO      : SOLAR RADIATION FLUX.
  
  ! - 2D (1:KLEV) .
  
  !-----------------------------------------------------------------------
  
  ! METHOD.
  !     --------
  
  ! AUTHOR.
  !     -------
  ! 2014-05-15, J.M. Piriou.
  
  ! MODIFICATIONS.
  !     --------------
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  
  !-----------------------------------------------------------------------
  
  
!$acc routine( ACDAYD_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMCST, ONLY: TCST
  USE YOMRIP, ONLY: TRIP
  USE YOMLSFORC, ONLY: LMUSCLFA, NMUSCLFA
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(TRIP), INTENT(IN) :: YDRIP
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KSGST
  
  REAL(KIND=JPRB), INTENT(IN) :: PGEMU(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PMU0(KLON)
  REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PDELP(KLON, KLEV)
  REAL(KIND=JPRB), INTENT(INOUT) :: PFRSO(KLON, 0:KLEV, KSGST + 1)
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM) :: JLON
  INTEGER(KIND=JPIM) :: JLEV
  temp (REAL (KIND=JPRB), ZTEND, (KLON, KLEV))
  temp (REAL (KIND=JPRB), ZDAYDUR, (KLON, KLEV))
  REAL(KIND=JPRB) :: ZZ
  REAL(KIND=JPRB) :: ZGEOM
  REAL(KIND=JPRB) :: ZSINLAT
  REAL(KIND=JPRB) :: ZCOSLAT
  REAL(KIND=JPRB) :: ZCOSPHI
  REAL(KIND=JPRB) :: ZEPSILON
  REAL(KIND=JPRB) :: ZPI
  
#include "wrscmr.intfb.h"
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  alloc (ZTEND)
  alloc (ZDAYDUR)
  JLON = KIDIA
  
  !-------------------------------------------------
  ! Day duration depending on altitude.
  !-------------------------------------------------
  
  ZPI = 4._JPRB*ATAN(1._JPRB)
  ZDAYDUR(JLON, :) = 0._JPRB
  DO JLEV=KTDIA,KLEV
    ZZ = PAPHIF(JLON, JLEV) / YDCST%RG
    ZGEOM = SQRT(MAX(0._JPRB, 1._JPRB - 1._JPRB / (1._JPRB + ZZ / YDCST%RA)**2))
    ZSINLAT = PGEMU(JLON)
    ZCOSLAT = SQRT(MAX(0._JPRB, 1._JPRB - ZSINLAT*ZSINLAT))
    ZCOSPHI = MAX(-0.999_JPRB, MIN(0.999_JPRB, -(ZGEOM + YDRIP%RSIDEC*ZSINLAT) / YDRIP%RCODEC / MAX(0.0001_JPRB, ZCOSLAT)))
    ZEPSILON = ZPI - ACOS(ZCOSPHI)
    ZDAYDUR(JLON, JLEV) = 86400._JPRB*(1._JPRB - ZEPSILON / ZPI)
  END DO
  
  ! The tendency is multiplied by a factor depending on day duration.
  DO JLEV=KTDIA,KLEV
    ZTEND(JLON, JLEV) = -YDCST%RG*(PFRSO(JLON, JLEV, 1) - PFRSO(JLON, JLEV - 1, 1)) / PDELP(JLON, JLEV)
    ZTEND(JLON, JLEV) = ZTEND(JLON, JLEV)*ZDAYDUR(JLON, JLEV) / ZDAYDUR(JLON, KLEV)
  END DO
  
  ! Recompute flux from modified tendency. This is done upwards in order to let
  ! surface flux unchanged.
  DO JLEV=KLEV,KTDIA,-1
    PFRSO(JLON, JLEV - 1, 1) = PFRSO(JLON, JLEV, 1) + PDELP(JLON, JLEV) / YDCST%RG*ZTEND(JLON, JLEV)
  END DO
  
END SUBROUTINE ACDAYD_OPENACC
