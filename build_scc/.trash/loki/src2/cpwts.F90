SUBROUTINE CPWTS_OPENACC (YDCST, YDCPG_OPTS, YDMCC, YDPHY, YDPHY1, KLON, KIDIA, KFDIA, KCSS, PDT, PTDTS, PTDTP, PTDWS, PTDWSI,  &
& PTDWP, PTDWPI, PTDWL, PTDSNS, PTDALBNS, PTDRHONS, PTPCLI, PWPCLI, PLSM, PIVEG, PTS1, PTP1, PWS1, PWSI1, PWP1, PWPI1, PWL1,  &
& PSNS1, PALBNS1, PRHONS1, YDSTACK)
  !----------------------------------------------------------------------
  ! - INPUT 1D
  ! - IN/OUT 1D
  
  !DEC$ NOOPTIMIZE
  
  !****
  !     ------------------------------------------------------------------
  
  !     EVOLUTION DES VARIABLES SOLS
  !     ------------------------------------------------------------------
  
  !     ARGUMENTS D ENTREE
  !     ------------------
  !       KLON   : DIMENSION HORIZONTALE
  !       KIDIA  : BORNE INITIALE HORIZONTALE DES CALCULS
  !       KFDIA  : BORNE FINALE HORIZONTALE DES CALCULS
  !       KCSS   : NBRE DE COUCHES DANS LE SOL PROFOND
  !       PDT    : PAS DE TEMPS EFFECTIF
  
  ! --- INPUT 1D
  !     --------
  !       PTDTS   (KLON): TENDANCE PHYSIQUE - TEMPERATURE DE SURFACE
  !       PTDTP(KLON,KCSS):  "     "        - TEMPERATURE PROFONDE
  !       PTDWS   (KLON):    "     "        - RESERVOIR D'EAU SUPERFICIEL
  !       PTDWSI  (KLON):    "     "        - RESERVOIR D'EAU SUPERFICIEL GELEE
  !                                       SI LSOLV ET LFGELS
  !       PTDWP   (KLON):    "     "        - RESERVOIR D'EAU PROFOND
  !       PTDWPI  (KLON):    "     "        - RESERVOIR D'EAU PROFOND GELEE
  !                                       SI LSOLV ET LFGEL
  !       PTDWL   (KLON):    "     "        - RESERVOIR D'EAU D'INTERCEPTION
  !                                       SI LSOLV
  !       PTDSNS  (KLON):    "     "        - RESERVOIR DE NEIGE
  !                                       SI LNEIGE
  !       PTDALBNS(KLON):    "     "        - ALBEDO DE LA NEIGE
  !       PTDRHONS(KLON):    "     "        - DENSITE DE LA NEIGE
  !                                       SI LSNV (ET LNEIGE)
  !       PTPCLI  (KLON): TEMPERATURE CLIMATOLOGIQUE DU SOL PROFOND
  !       PWPCLI  (KLON): RESERVOIR EAU CLIMATOLOGIQUE DU SOL PROFOND
  !       PLSM    (KLON): INDICATEUR TERRE (1.)/ MER (0.)
  !       PIVEG   (KLON): INDICE DE VEGETATION DOMINANTE (LSOLV)
  
  !     ARGUMENTS IMPLICITES
  !     --------------------
  !       INDICATEURS LOGIQUES    = COMMON /YOMPHY /: LSOLV,LNEIGE,LSNV,LFGEL,
  !                                                   LFGELS
  !       CONSTANTES UNIVERSELLES = COMMON /YOMCST /: RTT
  !       PROPRIETES DU SOL       = COMMON /YOMPHY1/: OMTPRO,OMWPRO,TMERGL
  !       CLIMAT                  = COMMON /YOMMCC/ : LMCCO2
  !       ACTION ON FIELDS        = YSD_VP%Y(field)%LSET
  
  !     SORTIES
  !     -------
  !       PTS1, PTP1, PWS1, PWSI1, PWP1, PWPI1, PWL1, PSNS1, PALBNS1, PRHONS1:
  !       CHAMPS DE SURFACE A FAIRE EVOLUER
  
  !     IMPROVISE PAR
  !     ------------- ERIC BAZILE 92/11 INSPIRE DE CPSP
  
  ! --- MODIFICATIONS.
  !     --------------
  !      R. El Khatib : 01-08-07 Pruning options
  !      Modified 02-09-03 by E. Bazile    Schema de neige LVGSN.
  !      Modified 06-09-10 by A. Alias     relaxation reservoir profond (Michel Deque)
  !      R. El Khatib : 09-Mar-2012 Cleaning for the sake of bounds checking controls
  !     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
  !     ------------------------------------------------------------------
  
  
!$acc routine( CPWTS_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: DR_HOOK, JPHOOK, LHOOK
  USE CPG_OPTS_TYPE_MOD, ONLY: CPG_OPTS_TYPE
  
  USE YOMCST, ONLY: TCST
  USE YOMPHY, ONLY: TPHY
  USE YOMPHY1, ONLY: TPHY1
  USE YOMMCC, ONLY: TMCC
  
  !     ------------------------------------------------------------------
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  TYPE(TCST), INTENT(IN) :: YDCST
  TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
  TYPE(TMCC), INTENT(IN) :: YDMCC
  TYPE(TPHY), INTENT(IN) :: YDPHY
  TYPE(TPHY1), INTENT(IN) :: YDPHY1
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KCSS
  REAL(KIND=JPRB), INTENT(IN) :: PDT
  REAL(KIND=JPRB), INTENT(IN) :: PTDTS
  REAL(KIND=JPRB), INTENT(IN) :: PTDTP(KLON, KCSS)
  REAL(KIND=JPRB), INTENT(IN) :: PTDWS
  REAL(KIND=JPRB), INTENT(IN) :: PTDWSI
  REAL(KIND=JPRB), INTENT(IN) :: PTDWP
  REAL(KIND=JPRB), INTENT(IN) :: PTDWPI
  REAL(KIND=JPRB), INTENT(IN) :: PTDWL
  REAL(KIND=JPRB), INTENT(IN) :: PTDSNS
  REAL(KIND=JPRB), INTENT(IN) :: PTDALBNS
  REAL(KIND=JPRB), INTENT(IN) :: PTDRHONS
  REAL(KIND=JPRB), INTENT(IN) :: PTPCLI
  REAL(KIND=JPRB), INTENT(IN) :: PWPCLI
  REAL(KIND=JPRB), INTENT(IN) :: PLSM
  REAL(KIND=JPRB), INTENT(IN) :: PIVEG
  REAL(KIND=JPRB), INTENT(INOUT) :: PTS1
  REAL(KIND=JPRB), INTENT(INOUT) :: PTP1(KLON, KCSS)
  REAL(KIND=JPRB), INTENT(INOUT) :: PWS1
  REAL(KIND=JPRB), INTENT(INOUT) :: PWSI1
  REAL(KIND=JPRB), INTENT(INOUT) :: PWP1
  REAL(KIND=JPRB), INTENT(INOUT) :: PWPI1
  REAL(KIND=JPRB), INTENT(INOUT) :: PWL1
  REAL(KIND=JPRB), INTENT(INOUT) :: PSNS1
  REAL(KIND=JPRB), INTENT(INOUT) :: PALBNS1
  REAL(KIND=JPRB), INTENT(INOUT) :: PRHONS1
  !     ------------------------------------------------------------------
  INTEGER(KIND=JPIM) :: ITVGLAM
  INTEGER(KIND=JPIM) :: JCSS
  INTEGER(KIND=JPIM) :: JROF
  
  REAL(KIND=JPRB) :: ZEPSN
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=jpim) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !     ------------------------------------------------------------------
  
  !*    ------------------------------------------------------------------
  
  !     1.- CALCUL DES EVOLUTIONS
  
  ZEPSN = 1.E-3_JPRB
  DO JCSS=1,KCSS
    DO JROF=KIDIA,KFDIA
      
      !        TEMPERATURE PROFONDE
      !        --------------------
      PTP1(JROF, JCSS) = PTP1(JROF, JCSS) + PDT*PTDTP(JROF, JCSS)
    END DO
  END DO
  
  !DEC$ IVDEP
  DO JROF=KIDIA,KFDIA
    
    !        TEMPERATURE DE SURFACE
    !        ----------------------
    PTS1 = PTS1 + PDT*PTDTS
    
    !        TEMPERATURE PROFONDE
    !        --------------------
    IF (YDCPG_OPTS%YRSURF_OPTS%YSD_VP%YTPC%LSET) THEN
      PTP1(JROF, KCSS) = PTP1(JROF, KCSS) + PDT*PLSM*YDPHY1%OMTPRO*(PTPCLI - PTP1(JROF, KCSS))
    END IF
    
    !        RESERVOIR SUPERFICIEL
    !        ---------------------
    PWS1 = PWS1 + PDT*PTDWS
    PWS1 = MAX(0.0_JPRB, PWS1)
    
  END DO
  
  IF (YDPHY%LSOLV) THEN
    !DEC$ IVDEP
    DO JROF=KIDIA,KFDIA
      
      !        RESERVOIR D'INTERCEPTION
      !        ------------------------
      PWL1 = PWL1 + PDT*PTDWL
      PWL1 = MAX(0.0_JPRB, PWL1)
      
      !        RESERVOIR PROFOND
      !        -----------------
      PWP1 = PWP1 + PDT*PTDWP
      IF (YDCPG_OPTS%YRSURF_OPTS%YSD_VP%YWPC%LSET) THEN
        PWP1 = PWP1 + PDT*PLSM*YDPHY1%OMWPRO*(PWPCLI - PWP1)
      END IF
      PWP1 = MAX(0.0_JPRB, PWP1)
      
      !        RESERVOIR PROFOND GELE
      !        ----------------------
      IF (YDPHY%LFGEL) THEN
        PWPI1 = PWPI1 + PDT*PTDWPI
        PWPI1 = MAX(0.0_JPRB, PWPI1)
        PWSI1 = PWSI1 + PDT*PTDWSI
        PWSI1 = MAX(0.0_JPRB, PWSI1)
      END IF
    END DO
  ELSE
    !DEC$ IVDEP
    DO JROF=KIDIA,KFDIA
      
      !        RESERVOIR PROFOND
      !        -----------------
      PWP1 = PWP1 + PDT*PTDWP
      IF (YDCPG_OPTS%YRSURF_OPTS%YSD_VP%YWPC%LSET) THEN
        PWP1 = PWP1 + PDT*PLSM*YDPHY1%OMWPRO*(PWPCLI - PWP1)
      END IF
    END DO
  END IF
  
  !        VARIABLES NEIGE EVENTUELLES
  !        ---------------------------
  IF (YDPHY%LNEIGE) THEN
    !DEC$ IVDEP
    DO JROF=KIDIA,KFDIA
      PSNS1 = PSNS1 + PDT*PTDSNS
      IF (YDPHY%LSOLV) THEN
        IF (NINT(PIVEG) == YDPHY1%NTVGLA) THEN
          PSNS1 = MAX(ZEPSN, PSNS1)
        END IF
      END IF
      PSNS1 = MAX(0.0_JPRB, PSNS1)
      IF (YDPHY%LSNV .or. YDPHY%LVGSN) THEN
        IF (PSNS1 > 0.0_JPRB) THEN
          PALBNS1 = PALBNS1 + PTDALBNS*PDT
          PALBNS1 = MIN(MAX(PALBNS1, YDPHY1%ALBMIN), YDPHY1%ALBMAX)
        ELSE
          PALBNS1 = YDPHY1%ALBMAX
        END IF
        IF (PSNS1 > 0) THEN
          PRHONS1 = PRHONS1 + PTDRHONS*PDT
          PRHONS1 = MIN(MAX(PRHONS1, YDPHY1%RHOMIN), YDPHY1%RHOMAX)
        ELSE
          PRHONS1 = YDPHY1%RHOMIN
        END IF
      END IF
    END DO
  END IF
  
  !      PSEUDO FONTE BANQUISE
  
  IF (YDMCC%LMCC02 .and. YDPHY%LSOLV) THEN
    ITVGLAM = 10*YDPHY1%NTVGLA + 1
    DO JCSS=1,KCSS
      DO JROF=KIDIA,KFDIA
        IF (NINT(10._JPRB*PIVEG) == ITVGLAM) THEN
          !               ON SEA-ICE THE TEMPERATURE CANNOT EXCEED 273.16 K
          !               OTHERWISE THE ICE MELTS TO MAINTAIN IT
          PTP1(JROF, JCSS) = MIN(PTP1(JROF, JCSS), YDCST%RTT)
        END IF
      END DO
    END DO
    DO JROF=KIDIA,KFDIA
      IF (NINT(10._JPRB*PIVEG) == ITVGLAM) THEN
        !               ON SEA-ICE THE TEMPERATURE CANNOT EXCEED 273.16 K
        !               OTHERWISE THE ICE MELTS TO MAINTAIN IT
        PTS1 = MIN(PTS1, YDCST%RTT)
      END IF
    END DO
  END IF
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE CPWTS_OPENACC
