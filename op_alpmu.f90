SUBROUTINE OP_ALPMU (ND)

  !+
  !
  !  Name:
  !     OP_ALPMU
  !
  !  Purpose:
  !     Calculate mean atomic weight in each zone of model 
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
  USE OP_GLOBAL, ONLY : OP_MAXNA 
  USE OP_GLOBAL, ONLY : OP__ALPMU, OP_ALPNN, OP_ALPAM 

  !  Implicit typing.
  IMPLICIT NONE

  !  Subroutine argument (Given).
  INTEGER :: ND     ! Number of depth points

  !  Local array indices
  INTEGER :: JD   !  loop over depth
  INTEGER :: IZ   !  loop over atomic number

  !.

  !  Loop through depth points
  DO JD = 1,ND
     OP__ALPMU(JD) = 0.0
     DO IZ=1,OP_MAXNA
        OP__ALPMU(JD) = OP__ALPMU(JD) + OP_ALPNN(IZ,JD)*OP_ALPAM(IZ)
     ENDDO
  ENDDO

  !  Trap abundance errors
  IF (OP__ALPMU(1) > 10.) WRITE (*,*) 'op_alpmu: ** warning ** mu=',OP__ALPMU(1:ND)

  ! write(*,*) 'OP_ALPMU -- OP__ALPMU :',op__alpmu(1:nd)

END SUBROUTINE OP_ALPMU
