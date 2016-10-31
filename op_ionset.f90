 
  subroutine op_ionset(JD,IZ,NI,X,UU)

  !+
  !
  !  Load ionisation fractions partition functions 
  !  for individual ions into the OP common blocks
  !
  !  
  !  Language:
  !     Fortran 77
  !
  !  Input arguments:
  !     JD                 !  Depth pointer
  !     IZ                 !  Atomic species pointer
  !     NI                 !  Number of ions for atom IZ
  !     X(*)               !  Ion fractions for IZ
  !     UU(*)    	         !  Partition functions for IZ 
  !
  !  Authors:
  !     CSJ: C. S. Jeffery (Armagh Observatory)
  !
  !  History:
  !     24-FEB-1999 (CSJ): New program (called by sterne2.4/stn2/press.f)
  !
  !-

  !  Global OP library variables
    USE OP_GLOBAL

    implicit none

  !  Subroutine Arguments (Given).
    INTEGER :: JD              !  Depth pointer
    INTEGER :: IZ              !  Atomic species pointer
    INTEGER :: NI              !  Number of ions for atom IZ

    real :: X(*)  !  Ion fractions for IZ
    real :: UU(*) !  Partition functions for IZ 

    integer :: i

  !  Loop through ions
    DO i = 1,NI
      OP_XION(i,IZ,JD) = X(I)
      OP_XUU(i,IZ,JD) = UU(I)
    ENDDO

  end subroutine op_ionset
