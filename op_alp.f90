SUBROUTINE OP_ALP

  !+
  !
  !  Set up atomic data and ion limits, save data to OP_GLOBAL module
  !  OP_ALPZ(:)     : atomic number
  !  OP_ALPAM(:)    : atomic weights
  !  OP_ALPIMIN(:) : ion range, lowest ion to be treated
  !  OP_ALPIMAX(:) : ion range, highest ion to be treated
  !
  !  NB -- in Sterne the old coposition array indices were ...
  !  :1=Hydrogen, 2=Helium, 3=Carbon, 4=Nitrogen,
  !  5=Oxygen, 6=Mg, 7=Al, 8=Si, 9=S, 10=Ca, 11=Fe
  !  This has been changed so that index coincides with the atomic number
  !
  !-

  USE OP_GLOBAL

  IMPLICIT NONE

  !.

  ! atomic numbers
  OP_ALPZ(1:40) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, &
       & 23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40/)

  ! atomic weights
  OP_ALPAM(1:40) = (/1.00794,4.002602,6.94,9.01,10.81,12.01,14.01,16.00,19.00,20.18,22.99,24.31, &
       & 26.98,28.09,30.97,32.06,35.45,39.95,39.10,40.08,44.96,47.87,50.94,52.00, &
       & 54.94,55.85,58.933,58.693,63.546,65.39,69.723,72.61,74.9216,78.96, &
       & 79.904,83.80,85.4678,87.62,88.90585,91.224/)

  ! lowest ion
  OP_ALPIMIN(1:40) = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, & 
       & 1,1,1,1,1,1,1,1,1,1/)

  ! highest ion
  OP_ALPIMAX(1:40) = (/1,2,3,4,4,6,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,6,6,6,6,6, &
       & 6,6,6,6,6,6,6,6,6,6/)

END SUBROUTINE OP_ALP
