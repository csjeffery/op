
  real FUNCTION SEATON(FREQ,FREQ0,XSECT,POWER,A)

      use op_precisions
      implicit none

      real(kind=qq), intent(in) :: FREQ
      real(kind=qq), intent(in) :: FREQ0, XSECT, POWER, A
    
      real :: IPWR

      IPWR   = 2.*POWER + 0.01
      SEATON = XSECT * (A+(1.-A) * (FREQ0/FREQ))* &
               & SQRT((FREQ0/FREQ)**IPWR) 

  END function SEATON
