MODULE OP_GLOBAL

  !+
  !
  !  Module:
  !     OP_GLOBAL
  !
  !  Purpose:
  !     Global variable holders and constants 
  !     for Opacity Project continuous opacities
  !
  !     ns	number of states
  !     alpimin,alpimax	min,max ion numbers being treated
  !     xwt	weights for states
  !     xchi	ionisation energy
  !     xion	fractional abundance for each ion
  !     xpop      fractional population for each state
  !     xs        cross sections
  !
  !  Authors:
  !     PMH: P.M.Harrison (Armagh)
  !     CSJ: C.S.Jeffery  (Armagh)
  !
  !  History:
  !     11-Feb-1998 (PMH): Created as an include file with common blocks
  !     25-Feb-2003 (CSJ): Converted to module  
  !
  !-

  ! User defined types
  USE OP_PRECISIONS

  ! No implicit typing
  IMPLICIT NONE

  PRIVATE :: QQ, RR

  !  Array index parameters.
  INTEGER, PARAMETER :: OP_MAXNA=40        !  atoms with x-sectioins
  INTEGER, PARAMETER :: OP_MAXNZ=90        !  atoms overall
  INTEGER, PARAMETER :: OP_MAXNI=6         !  ion species
  INTEGER, PARAMETER :: OP_MAXNS=50        !  ion states
  INTEGER, PARAMETER :: OP_MAXND=201       !  depth points
  INTEGER, PARAMETER :: OP_MAXNW=10000     !  wavelengths
  INTEGER, PARAMETER :: OP_MAXNX=10000     !  wavelengths in cross-section data

  !  Array of selected atoms . 
  INTEGER :: OP_NATOMS                       ! number of atoms selected for inclusion
  INTEGER, DIMENSION(OP_MAXNZ) :: OP_SATOMS  ! array of atomic numbers selected
  INTEGER :: OP_IAS                          ! array index for indirect referencing 

  !  Global variables.
  INTEGER :: OP_NS(OP_MAXNI,OP_MAXNA)            ! number of levels
  INTEGER :: OP_NX(OP_MAXNS,OP_MAXNI,OP_MAXNA)   ! number of x-sections

  INTEGER :: OP_ALPIMIN(OP_MAXNA)
  INTEGER :: OP_ALPIMAX(OP_MAXNA)                ! ion range
  INTEGER :: OP_ALPZ(OP_MAXNA)                   ! atomic numbers
  REAL :: OP_ALPAM(OP_MAXNA)                     ! atomic weights
  REAL :: OP_ALPNN(OP_MAXNA,OP_MAXND)            ! atomic abundances
  REAL :: OP__ALPMU(OP_MAXND)                     ! mean at. weights
  REAL :: OP_ALPMF(OP_MAXNA,OP_MAXND)            ! mass fractions


  REAL :: OP_XCHI(OP_MAXNS,OP_MAXNI,OP_MAXNZ)    ! ionisation energies
  REAL :: OP_XION(OP_MAXNI,OP_MAXNZ,OP_MAXND)    ! fractional ion abundances
  REAL :: OP_XUU(OP_MAXNI,OP_MAXNZ,OP_MAXND)        ! partition functions
  REAL :: OP_XWT(OP_MAXNS,OP_MAXNI,OP_MAXNA)     ! weights for levels
  REAL :: OP_XS(OP_MAXNW+1,OP_MAXNS,OP_MAXNI,OP_MAXNA) ! cross-sections
  REAL :: OP_XPOP(OP_MAXNS,OP_MAXNI,OP_MAXNA,OP_MAXND) ! level populations
  REAL :: OP_XSSUM(OP_MAXNW+1,OP_MAXNI,OP_MAXNA) ! summed cross-sections


  INTEGER :: OP_NW = 0                           ! number of wavelength points
  REAL, ALLOCATABLE :: OP_WL(:)         ! wavelength grid assoc. with model
  REAL, ALLOCATABLE :: OP_NU(:)         ! frequency grid assoc. with model

  !  Holders for OP energy level information

  REAL :: OP_ECONT ( OP_MAXNA, OP_MAXNI, OP_MAXNS )     !  Ionisation energies
  INTEGER :: OP_ISLP  ( OP_MAXNA, OP_MAXNI, OP_MAXNS )  !  Series
  INTEGER :: OP_ILV   ( OP_MAXNA, OP_MAXNI, OP_MAXNS )  !  Level
  CHARACTER(LEN=15):: OP_ICONF ( OP_MAXNA, OP_MAXNI, OP_MAXNS )  !  Configurations

  !! Holders for free-free continuous opacities

  REAL :: OPFFH1, OPFFHM, OPFFHE1, OPFFHE2, OPFFHEM
  REAL :: OPFFC1, OPFFC2, OPFFC2X, OPFFC3, OPFFCM
  REAL :: OPFFN1, OPFFN2, OPFFN3
  REAL :: OPFFMG1, OPFFAL1, OPFFSI1, OPFFO1
  REAL :: OPFFMG2, OPFFSI2, OPFFCA2

  !! Holders for individual bound-free opacities

  REAL :: KAPPABF(OP_MAXNA,OP_MAXNI)   

  REAL :: OPAC_REF
  REAL :: OPAC_NU

  !! depth-dependent atomic data for ff-opacities
  DOUBLE PRECISION, DIMENSION(22,OP_MAXND) :: AT     ! T dependent 


END MODULE OP_GLOBAL
