  SUBROUTINE OP_FFH1 (IL,T,FRE,ZH1, AN17)

!     Total free-free absorption coefficient per atom of hydrogen 
!     The gaunt factors are taken from Gingerich(1964)

      USE OP_PRECISIONS
      USE OP_GLOBAL, ONLY : AT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IL
      REAL, INTENT(IN) :: T
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL, INTENT(IN) :: ZH1
      REAL, INTENT(OUT) :: AN17

      REAL, PARAMETER :: HH = 4.7993E-11 
      REAl, PARAMETER :: RC = 3.2894E+15 

      REAL :: GNFF
      REAL(KIND=QQ) :: F, AF

      F = FRE/1.0000E+14
      GNFF = 1.+(0.1728*(FRE/RC)**0.33333333)*(1.+(2.*T)/(HH*FRE))

      AF = 2.81592E-13/F**3
      AN17 = GNFF*AF*AT(6,IL) * 2./ZH1

      ! write (*,*), 'op_ffh1',fre,f,gnff,af,at(6,il),an17,zh1

  END SUBROUTINE OP_FFH1
