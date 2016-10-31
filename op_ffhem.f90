SUBROUTINE OP_FFHEM (T,FRE,PE,ZHE1, AN5)

!     free-free absorption coefficient of helium minus
!     T.L. John, MNRAS 138, 137 (1968), 

      USE OP_PRECISIONS

      IMPLICIT NONE

      REAL, INTENT(IN) :: T, PE
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL, INTENT(IN) :: ZHE1
      REAL, INTENT(OUT) :: AN5

      REAL :: A, B, C

      A =  3.397E-36 + (-5.216E-21 + 7.039E-05/FRE)/FRE
      B = -4.116E-32 + ( 1.067E-16 + 8.135E-01/FRE)/FRE
      C =  5.081E-27 + (-8.724E-13 - 5.659E+02/FRE)/FRE
      AN5 = (A*T +B +C/T)*PE/1.38E-16/T/ZHE1/1.0E+10

END SUBROUTINE OP_FFHEM
