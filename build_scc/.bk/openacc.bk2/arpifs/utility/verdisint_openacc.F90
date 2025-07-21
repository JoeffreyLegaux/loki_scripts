SUBROUTINE VERDISINT_OPENACC (YDVFE, YDCVER, CDOPER, CDBC, KPROMA, KST, KEND, KFLEV, PIN, POUT, KCHUNK, YDSTACK)
  
  !**** *VERDISINT*   VERtical DIScretization -
  !                INTerface for finite element type vertical operations:
  !                derivative or integral
  
  !     Purpose.
  !     --------
  !          This subroutine prepares an interface to VERINT
  !          computing vertical integral with respect to eta
  !          and VERDER computing vertical derivative with
  !          respect to eta of a function given at full or
  !          half model levels using a general scheme.
  
  !**   Interface.
  !     ----------
  !        *CALL* *VERDISINT(..)
  
  !        Explicit arguments :
  !        --------------------
  
  !        INPUT:
  !          CDOPER    - type of integral or derivative applied:
  !                      'ITOP' - integral from top
  !                      'IBOT' - integral from bottom
  !                      'INGW' - invertible integral operator
  !                      'HDER' - first derivative at half levels
  !                      'FDER' - first derivative on full levels
  !                      'DDER' - second derivative on full levels
  !                      'DEGW' - invertible derivative operator
  !          CDBC      - boundary conditions used ('00','01','10','11',
  !                      first digit for top, second digit for bottom,
  !                      0 for value prescribed, 1 for derivative prescribed)
  !          KPROMA    - horizontal dimension.
  !          KSTART    - first element of work.
  !          KPROF     - depth of work.
  !          KFLEV     - vertical dimension for array PIN.
  !          PIN       - input field
  !        OPTIONAL INPUT:
  !          KCHUNK    - Use NPROMA as blocking factor when called outside
  !                      OpenMP threaded region
  
  !        OUTPUT:
  !          POUT      - integral or derivative of PIN according to CDOPER
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  !        none
  
  !     Reference.
  !     ----------
  
  !     Author.
  !     -------
  !        P. Smolikova (CHMI/LACE/ALADIN)
  
  !     Modifications.
  !     --------------
  !        Original : Sep 2017
  !        P.Smolikova (Sep 2020): VFE pruning.
  !     ------------------------------------------------------------------
  
!$acc routine( VERDISINT_OPENACC )
  
  USE PARKIND1, ONLY: JPIM, JPRB, JPRD
  USE YOMCVER, ONLY: TCVER
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN, ONLY: NULERR
  USE YOMVERT, ONLY: TVFE
  USE STACK_MOD
#include "stack.h"
  
  
  !     ---------------------------------------------------------------------------
  
  IMPLICIT NONE
  TYPE(TVFE), TARGET, INTENT(IN) :: YDVFE
  TYPE(TCVER), INTENT(IN) :: YDCVER
  CHARACTER(LEN=*), INTENT(IN) :: CDOPER
  CHARACTER(LEN=*), INTENT(IN) :: CDBC
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
  INTEGER(KIND=JPIM), INTENT(IN) :: KST
  INTEGER(KIND=JPIM), INTENT(IN) :: KEND
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV
  INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KCHUNK
  REAL(KIND=JPRB), INTENT(IN) :: PIN(KPROMA, 0:KFLEV + 1)
  REAL(KIND=JPRB), INTENT(OUT) :: POUT(KPROMA, KFLEV + 1)
  
  !     ------------------------------------------------------------------
  
  CHARACTER(LEN=2) :: CLBC
  INTEGER(KIND=JPIM) :: ILEVIN
  INTEGER(KIND=JPIM) :: ILEVOUT
  INTEGER(KIND=JPIM) :: IND
  INTEGER(KIND=JPIM) :: IEND
  INTEGER(KIND=JPIM) :: ITYPE
  INTEGER(KIND=JPIM) :: JCHUNK
  REAL(KIND=JPRD), POINTER, CONTIGUOUS :: ZOPER(:, :)
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  !     ------------------------------------------------------------------
  
#include "abor1.intfb.h"
#include "verder_openacc.intfb.h"
#include "verint_openacc.intfb.h"
  INTEGER(KIND=    JPIM) :: JLON
  TYPE(STACK), INTENT(IN) :: YDSTACK
  TYPE(STACK) :: YLSTACK
  YLSTACK = YDSTACK
  JLON = KST
  
  !     ------------------------------------------------------------------
  !     ------------------------------------------------------------------
  
  !*   1. Initialization and control
  !    -----------------------------
  
  IF (YDCVER%LVFE_ECMWF) THEN
    IF (CDOPER /= 'FDER' .or. YDCVER%LVFE_NOBC) THEN
      CLBC = 'XX'        ! no BC applied for derivatives
    ELSE
      CLBC = 'BC'        ! BC according to Hortal applied on derivatives
    END IF
  ELSE
    CLBC = CDBC
  END IF
  
  
  !*   2. Set operator size according to boundary conditions applied
  !    --------------------------------------------------------------
  
  IF (CDOPER == 'FDER' .or. CDOPER == 'DDER') THEN
    ILEVOUT = KFLEV
  ELSE
    ILEVOUT = KFLEV + 1
  END IF
  
  IF (CDOPER == 'IBOT') THEN
    ITYPE = 1
  ELSE
    ITYPE = 0
  END IF
  
  IF (CDOPER == 'INGW' .or. CDOPER == 'DEGW') THEN
    IND = 0
    ILEVIN = KFLEV + 1
  ELSE IF (CLBC == 'XX') THEN
    IND = 1
    ILEVIN = KFLEV
  ELSE
    IND = 0
    ILEVIN = KFLEV + 2
  END IF
  IEND = ILEVIN - 1 + IND
  
  !*   3. Set operator according to boundary conditions applied
  !    --------------------------------------------------------
  
  !*   3.1 Vertical integral
  !    ---------------------
  IF (CDOPER == 'INGW') THEN
    ZOPER => YDVFE%RINTGW(:, :)
  ELSE IF (ANY(CDOPER == (/ 'ITOP', 'IBOT' /))) THEN
    IF (CLBC == 'XX') THEN
      ZOPER => YDVFE%RINTE(:, :)
    ELSE IF (CLBC == '00') THEN
      ZOPER => YDVFE%RINTBF00(:, :)
    ELSE IF (CLBC == '11') THEN
      ZOPER => YDVFE%RINTBF11(:, :)
    ELSE
      PRINT *,  'VERDISINT: ITOP/IBOT NOT IMPLEMENTED FOR CDBC=', CLBC
      CALL ABOR1_ACC(' VERDISINT: ABOR1 CALLED')
    END IF
  END IF
  
  !*   3.2 Vertical derivative
  !    -----------------------
  
  IF (CDOPER == 'DEGW') THEN
    ZOPER => YDVFE%RDERGW(:, :)
  ELSE IF (CDOPER == 'HDER') THEN
    IF (CLBC == '00') THEN
      ZOPER => YDVFE%RDERBH00(:, :)
    ELSE IF (CLBC == '01') THEN
      ZOPER => YDVFE%RDERBH01(:, :)
    ELSE
      PRINT *,  'VERDISINT: HDER NOT IMPLEMENTED FOR CDBC=', CLBC
      CALL ABOR1_ACC(' VERDISINT: ABOR1 CALLED')
    END IF
  ELSE IF (CDOPER == 'FDER') THEN
    IF (CLBC == 'XX') THEN
      ZOPER => YDVFE%RDERI(:, :)
    ELSE IF (CLBC == 'BC') THEN
      ZOPER => YDVFE%RDERB(:, :)
    ELSE IF (CLBC == '00') THEN
      ZOPER => YDVFE%RDERBF00(:, :)
    ELSE IF (CLBC == '01') THEN
      ZOPER => YDVFE%RDERBF01(:, :)
    ELSE IF (CLBC == '10') THEN
      ZOPER => YDVFE%RDERBF10(:, :)
    ELSE IF (CLBC == '11') THEN
      ZOPER => YDVFE%RDERBF11(:, :)
    ELSE
      PRINT *,  'VERDISINT: FDER NOT IMPLEMENTED FOR CDBC=', CLBC
      CALL ABOR1_ACC(' VERDISINT: ABOR1 CALLED')
    END IF
  ELSE IF (CDOPER == 'DDER') THEN
    IF (CLBC == 'XX') THEN
      ZOPER => YDVFE%RDDERI(:, :)
    ELSE IF (CLBC == '01') THEN
      ZOPER => YDVFE%RDDERBF01(:, :)
    ELSE
      PRINT *,  'VERDISINT: DDER NOT IMPLEMENTED FOR CDBC=', CLBC
      CALL ABOR1_ACC(' VERDISINT: ABOR1 CALLED')
    END IF
  END IF
  
  !*   4. Apply the required operation
  !    --------------------------------
  
  IF (ANY(CDOPER == (/ 'INGW', 'ITOP', 'IBOT' /))) THEN
    IF (PRESENT(KCHUNK)) THEN
      JCHUNK = KCHUNK
    ELSE
      JCHUNK = 1
    END IF
    CALL VERINT_OPENACC(KPROMA, KST, KEND, ILEVIN, ILEVOUT, ZOPER, PIN(:, IND:IEND), POUT, ITYPE, KCHUNK=JCHUNK, YDSTACK=YLSTACK)
  ELSE IF (ANY(CDOPER == (/ 'FDER', 'HDER', 'DEGW', 'DDER' /))) THEN
    CALL VERDER_OPENACC(KPROMA, KST, KEND, ILEVIN, ILEVOUT, ZOPER, PIN(:, IND:IEND), POUT, YDSTACK=YLSTACK)
  ELSE
    PRINT *,  'VERDISINT: UNKNOWN CDOPER=', CDOPER
    CALL ABOR1_ACC(' VERDISINT: ABOR1 CALLED')
  END IF
  
  !     ------------------------------------------------------------------
  
END SUBROUTINE VERDISINT_OPENACC
