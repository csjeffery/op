real function U (IZ,IJ,T,PE)

  !+
  !
  !  Purpose.
  !=    Calculate partition functions
  !
  !  Input.
  !     IZ    Atomic number  (1=H, 2=He, etc...)
  !     IJ    Ion identifier (1=neutral, 2=1st ion, etc...)
  !     T     Temperature  (K)
  !     PE    Electron pressure (cgs)
  !
  !  Output.
  !     U     Partition function
  !
  !  History.
  !     02-MAR-1993 (CSJ)
  !        Now uses TAP version of Traving et al. algorithm
  !
  !-

  implicit none

  !  External function declaration
  double precision, external :: TAP_PARFN

  !  Subroutine Arguments:
  integer, intent(in) :: IZ
  integer, intent(in) :: IJ
  real, intent(in) :: T, PE

  !  Local variables.
  double precision :: P_EL, TEMP
  integer :: ION_Z
  integer :: ION_N

  !.


  !  Convert arguments to SI units.
  ION_Z = IZ
  ION_N = IJ - 1
  TEMP  = T
  P_EL  = PE / 10.0

  !  Call the TAP routine.
  U = TAP_PARFN ( ION_Z, ION_N, TEMP, P_EL )
  ! write (*,*) ion_z, ion_n, temp, p_el, u

end function U
