SUBROUTINE OP_INXS

  !+
  !
  !  Name:
  !     OP_INXS
  !
  !  Purpose:
  !     Interpolate cross-section for current wavelength
  !
  !  Module Type:
  !     Subroutine
  !
  !  Language:
  !     Fortran 77
  !
  !  Output arguments:
  !     IWL    wavelength index (set to dummy reference value)
  !
  !  Authors:
  !     PMH: P.M.Harrison (St.Andrews)
  !     CSJ: C.S.Jeffery (Armagh Observatory)
  !
  !  History:
  !     19-Mar-1998 (PMH): Created.
  !     24-Feb-1999 (CSJ): Revised 
  !     13-Nov-2003 (NTB): modified for sterne90
  !-


  !  Global variables
  USE LTE_WLGRID
  USE OP_GLOBAL

  IMPLICIT NONE

  !  Subroutine Arguments.
  INTEGER :: IWL       !  wavelength pointer  (returned)

  !  Local constants.
  INTEGER, parameter :: NORDER = 2            ! order of polynomial interpolation

  !  Local variables.
  INTEGER :: I, II, IA, IS
  INTEGER :: IWR
  REAL :: DUMXS2(NORDER+1)
  REAL :: DUMXS, DUMDY

  !.

  !  Set dummy wavelength index
  IWL = OP_MAXNW+1

  !  Find index for interpolation
  IWR=0
  DO I=1,OP_NW-1
     IF (OP_WL(I) < REFWL) then
        IWR = I
     ENDIF
  ENDDO
  ! write(*,*) 'op_inxs: ',iwl,refwl,iwr,iwr+1,op_nw
  IF ((IWR == 0).OR.(IWR == OP_NW-1)) THEN
     WRITE(*,*) 'OP_INXS ERROR: REFWL OUT OF WL RANGE',REFWL
     STOP
  ENDIF

  !  Interpolate cross sections (linear)
  DO IA = 1,OP_MAXNA                     ! Loop over atoms
     IF (OP_ALPIMIN(IA).NE.0) THEN
        DO II = OP_ALPIMIN(IA),OP_ALPIMAX(IA)     ! Loop over ions
           DO IS = 1,OP_NS(II,IA)                 ! loop over states

              DUMXS = ( (REFWL-OP_WL(IWR))   * OP_XS(IWR+1,IS,II,IA) +  & 
                   (OP_WL(IWR+1)-REFWL) * OP_XS(IWR,IS,II,IA)     ) &
                   / ( OP_WL(IWR+1)-OP_WL(IWR) )

              OP_XS(IWL,IS,II,IA)=DUMXS

           ENDDO
        ENDDO
     ENDIF
  ENDDO


END SUBROUTINE OP_INXS
