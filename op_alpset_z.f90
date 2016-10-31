SUBROUTINE OP_ALPSET_Z (ND, NZ, ALPHA_Z)

  !+
  !
  !  Load abundance fractions for each atomic species 
  !  into the OP common blocks.
  !  This version currently assumes homogeneous composition
  !  but envisages a version for inhomogeneous structures. 
  !
  !  
  !  Language:
  !     Fortran 77
  !
  !  Input arguments:
  !     ND                 !  Number of depth points in model
  !     NZ                 !  Atomic number
  !     ALPHA_Z(NJ)        !  Abundances (number fraction) of NZ
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


  !  Subroutine Arguments (Given).
  INTEGER, INTENT(IN)        ::  ND        !  Number of depth points	       
  INTEGER, INTENT(IN)        ::  NZ        !  Atomic number
  DOUBLE PRECISION, INTENT(IN) :: ALPHA_Z(*)  !  Abundances (number fraction) of NZ 

  !.

  !  Set composition data
  OP_ALPNN(NZ,1:ND) = ALPHA_Z(1:ND)

END SUBROUTINE OP_ALPSET_Z
