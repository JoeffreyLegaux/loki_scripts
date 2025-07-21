SUBROUTINE SUOZON_OPENACC (KIDIA, KFDIA, KLON, KLEV, PROFO3, LDQINT, PRESI, PRDELP, LD_LO3ABC, PVO3A, PVO3B, PVO3C, YDSTACK)
  
  !     ROUTINE DE CALCUL ACTIF .
  !**** *SUOZON* INITIALISATION DE PROFILS VERTICAUX D'OZONE
  
  !     But.   INITIALISATION DE UN OU PLUSIEURS PROFILS VERTICAUX D'OZONE
  !     ----   A PARTIR D'UN PROFIL CLIMATOLOGIQUE,
  !           SOIT EN QUANTITE INTEGREE DOUBLE AU DESSUS D'UNE COUCHE (LDQINT=.TRU
  !           SOIT EN RAPPORT DE MELANGE (LDQINT=.FALSE.)
  
  !**   Interface.
  !     ----------
  !        *CALL* *SUOZON ( KIDIA,KFDIA,KLON,KLEV,
  !    S                    PROFO3,LDQINT,PRESI,PRDELP)
  
  !        Arguments explicites :
  !        ----------------------
  
  !       KIDIA    : DEBUT DES BOUCLES HORIZONTALES
  !                  BEGINNING OF HORIZONTAL LOOPS
  !       KFDIA    : FIN DES BOUCLES HORIZONTALES
  !                  END OF HORIZONTAL LOOPS
  !       KLON      : DIMENSION HORIZONTALE                       (input)
  !       KLEV      : DIMENSION VERTICALE                         (input)
  !       PROFO3   : PROFIL VERTICAL D'OZONE
  !              (0): VALEUR AU-DESSUS DU MODELE                  (output)
  !       LDQINT     : INDICATEUR DU TYPE DE PROFIL
  !                   .TRUE.  QUANTITE TOTALE AU DESSUS D'UNE COUCHE x2
  !                   .FALSE. RAPPORT DE MELANGE MASSIQUE         (input)
  !       PRESI    : PRESSION DE L'INTER-COUCHE                  (input)
  !       PRDELP    : INVERSE DE L'EPAISSEUR EN PRESSION
  !                   DE LA COUCHE                                (input)
  !       LO3ABC    : Switch to use climatological profile        (input)
  !       PVO3ABC   : Climatological coef                         (input)
  
  !        Arguments implicites :
  !        ----------------------
  !       Les points de grille du modele.
  
  !     Methode.
  !     --------
  !         L'ozone est une fonction de la pression P
  
  !        |0                a
  !        | qO3 dp =  --------------
  !        |P           1 + (b/P)**3/2
  
  !           a = ZQO31
  !           b = ZQO32
  
  !     On peut donc approximer le rapport de melange d'ozone
  
  !                        1    |Pj
  !              QO3 =  ------- |   qO3 dp
  !                      Pi-Pj  |Pi
  
  !     Reference.
  !    -----------
  !         AUCUNE
  
  !     Auteur.
  !    --------
  !         A. Lasserre-Bigorry
  
  !     Modifications :
  !    ----------------
  !         Original : 91-07-22
  !         M. Deque 91-09-25
  !         Y. Bouteloup 02-03-30  : Use of climatological profiles
  !        M.Hamrud      01-Oct-2003 CY28 Cleaning
  
  !     ------------------------------------------------------------------
  
!$acc routine( SUOZON_OPENACC ) seq
  
  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  
  USE STACK_MOD
#include "stack.h"
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  REAL(KIND=JPRB), INTENT(OUT) :: PROFO3(KLON, 0:KLEV)
  LOGICAL, INTENT(IN) :: LDQINT
  REAL(KIND=JPRB), INTENT(IN) :: PRESI(KLON, 0:KLEV)
  REAL(KIND=JPRB), INTENT(IN) :: PRDELP(KLON, KLEV)
  LOGICAL, INTENT(IN) :: LD_LO3ABC
  REAL(KIND=JPRB), INTENT(IN) :: PVO3A
  REAL(KIND=JPRB), INTENT(IN) :: PVO3B
  REAL(KIND=JPRB), INTENT(IN) :: PVO3C
  INTEGER(KIND=JPIM) :: JLEV
  INTEGER(KIND=JPIM) :: JLON
  
  REAL(KIND=JPRB) :: ZEPSO
  REAL(KIND=JPRB) :: ZP
  REAL(KIND=JPRB) :: ZQO3A
  REAL(KIND=JPRB) :: ZQO3B
  REAL(KIND=JPRB) :: ZQO3C
  REAL(KIND=JPRB) :: ZQO31
  REAL(KIND=JPRB) :: ZQO32
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KIDIA
  
  !*
  !     ------------------------------------------------------------------
  !     1. CALCULS PRELIMINAIRES ET DIVERSES INITIALISATIONS
  
  ZEPSO = 1.E-03_JPRB
  ZQO31 = 0.6012E-01_JPRB
  ZQO32 = 0.3166E+04_JPRB
  ZQO32 = SQRT(ZQO32*ZQO32*ZQO32)
  
  !*
  !     ------------------------------------------------------------------
  !     2.  INITIALISATION DE L'OZONE
  
  DO JLEV=0,KLEV
    IF (LD_LO3ABC) THEN
      ZQO3A = PVO3A
      ZQO3B = PVO3B
      ZQO3C = PVO3C / 2._JPRB
      ZP = MAX(ZEPSO, PRESI(JLON, JLEV))
      PROFO3(JLON, JLEV) = (2.0_JPRB*ZQO3A) / (1.0_JPRB + EXP(ZQO3C*LOG(ZQO3B / ZP)))
    ELSE
      ZP = MAX(ZEPSO, PRESI(JLON, JLEV))
      PROFO3(JLON, JLEV) = (2.0_JPRB*ZQO31) / (1.0_JPRB + ZQO32 / SQRT(ZP*ZP*ZP))
    END IF
  END DO
  IF (.not.LDQINT) THEN
    DO JLEV=KLEV,1,-1
      PROFO3(JLON, JLEV) = (PROFO3(JLON, JLEV) - PROFO3(JLON, JLEV - 1))*PRDELP(JLON, JLEV)*0.5_JPRB
    END DO
    PROFO3(JLON, 0) = PROFO3(JLON, 0) / MAX(ZEPSO, PRESI(JLON, 0))*0.5_JPRB
  END IF
  
END SUBROUTINE SUOZON_OPENACC
