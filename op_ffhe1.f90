  SUBROUTINE OP_FFHE1 (IL,T,FRE,ZHE1, AN311)

!     free-free He I cross-section

      USE OP_PRECISIONS
      USE OP_GLOBAL, ONLY : AT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IL
      REAL, INTENT(IN) :: T
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL, INTENT(IN) :: ZHE1
      REAL :: AN311

      REAL(KIND=QQ) :: F

      F = FRE/1.E14
      AN311= F**(-3.04) * AT(16,IL) / ZHE1

  END SUBROUTINE OP_FFHE1
