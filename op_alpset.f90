SUBROUTINE OP_ALPSET (ND, ALPHA)

  !+
  !
  !  Name:
  !     OP_ALPSET
  !
  !  Purpose:
  !     Load composition (homogeneous) into OP global
  !  
  !  Language:
  !     Fortran 90
  !
  !  Authors:
  !     CSJ: C. S. Jeffery (Armagh Observatory)
  !
  !  History:
  !     24-FEB-1999 (CSJ): New program (called by sterne2.4/main/star.f)
  !     04-SEP-2012 (CSJ): Upgraded for SterneOS and Spectrum
  !
  !-

  !  Global OP library variables 
  USE OP_GLOBAL, ONLY : OP_ALPNN

  !  Implicit typing.
  IMPLICIT NONE

  !  Subroutine argument (Given).
  INTEGER, INTENT(IN) :: ND     ! Number of depth points
  DOUBLE PRECISION, INTENT(IN) :: ALPHA(*)     ! Composition vector

  !  Local array indices
  INTEGER :: JD   !  loop over depth

  !  Sterne default abundance vector      
  INTEGER, DIMENSION(11) :: STN__K = (/ 1,2,6,7,8,12,13,14,16,20,26 /)

  !.

  !  Loop through depth points
  DO JD=1,ND
     OP_ALPNN(STN__K,JD) = ALPHA(STN__K)
  ENDDO

END SUBROUTINE OP_ALPSET
