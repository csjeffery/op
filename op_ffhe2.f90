SUBROUTINE OP_FFHE2 (IL,T,FRE,ZHE2, AN47)

   !     Total free-free continuous absorption coefficient per atom of helium II
   !     The gaunt factors are taken from Gingerich (1964)

   USE OP_PRECISIONS
   USE OP_GLOBAL, ONLY : AT

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: IL
   REAL, INTENT(IN) :: T
   REAL(KIND=QQ), INTENT(IN) :: FRE
   REAL, INTENT(IN) :: ZHE2
   REAL :: AN47

   REAL, PARAMETER :: HH = 4.79926E-11 
   REAL, PARAMETER :: RC = 3.28940E+15
   REAL(KIND=QQ) :: F
   REAL :: GNFF

   !.

   F = FRE/1.0000E+14
   GNFF = 1.+(0.1728*(FRE/(4.*RC))**0.33333333*(1.+(2.*T)/(HH*FRE)))
   AN47 = AT(22,IL) * GNFF/F**3 * T*3.5677E-18
   AN47 = AN47*2./ZHE2

END SUBROUTINE OP_FFHE2
