SUBROUTINE OP_EDGEFIND ( WS, WE, EDGE_M, EDGE_WL, EDGE_ID, EDGE_N )

  !+
  !
  !  Subroutine:
  !     OP_EDGEFIND
  !
  !  Purpose:
  !     Return a list of opacity edges within a given wavelength interval
  !
  !  Method:
  !     The OP library contains a list of ionization energies. This
  !     is polled to locate edges within th egiven range and these are
  !     returned, with an identifier.
  !     The number of edges detected is limited by the dimension of EDGE_WL,
  !     no warning is given if this limit is reached, but the subrouine will
  !     terminate safely. 
  !
  !  Arguments (input):
  !     WS           start wavelength of interval to search for edges
  !     WE           end wavelength
  !     EDGE_M       Dimension of edge arrays
  !
  !  Arguments (returned):
  !     EDGE_WL(*)   list of wavelngths of opacity edges (in Angstrom)
  !     EDGE_ID(*)   Identification for each opacity edge (forma; zzee.sssll)
  !                     zz = number of protons
  !                     ee = number of electron
  !                    sss = state identifier
  !                     ll = level identifier
  !     EDGE_N       Number of edges
  !
  !  Authors:
  !     CSJ : C.S.Jeffery (Armagh Observatory)
  !
  !  History:
  !     26 Feb 2003 (CSJ) : Original version
  !
  !-

  ! OP library global variables
  USE OP_GLOBAL
  USE TAPP_PAR

  ! No implicit typing
  IMPLICIT NONE


  ! Subroutine Arguments (given)
  DOUBLE PRECISION, INTENT(IN) :: WS !  start wavelength of interval 
  DOUBLE PRECISION, INTENT(IN) :: WE !  end wavelength  
  INTEGER, INTENT(IN) :: EDGE_M      !  dimension of explicit shape arrays

  ! Subroutine Arguments (returned)

  !  list of wavelngths of opacity edges (in Angstrom)
  DOUBLE PRECISION, DIMENSION(EDGE_M), INTENT(OUT) :: EDGE_WL 
  !  identification for each opacity edge
  DOUBLE PRECISION, DIMENSION(EDGE_M), INTENT(OUT) :: EDGE_ID 
  INTEGER, INTENT(OUT) :: EDGE_N    ! Number of edges


  !  Local constants.
  DOUBLE PRECISION, PARAMETER :: RINV = 1.E10 / TAPP__RYD ! 1/R (A)Shell No. 3

  ! Local variables
  INTEGER :: I, J, K  ! Loop indices
  INTEGER :: ZZ  ! Number of protons
  INTEGER :: EE  ! Number of electrons
  INTEGER :: LL  ! Level identifier
  DOUBLE PRECISION :: WL  ! Wavelngth of opacity edge
  INTEGER :: EDGE_MAX     ! Max number of edges allowed

  !.

  ! Set counter
  EDGE_N = 0

  ! Loop over atoms
  ATOMS : DO I = 1, OP_MAXNA
     ZZ = I

     ! Loop over ions
     DO J = 1, OP_MAXNI
        EE = ZZ - J + 1

        ! Loop over states
        DO K = 1, OP_MAXNS
           WL = - OP_ECONT ( I,J,K ) * RINV

           ! Identify opacity edges and store all in interval
           IF ( WS < WL .AND. WL < WE ) THEN
              EDGE_N          = EDGE_N + 1
              EDGE_WL(EDGE_N) = WL
              EDGE_ID(EDGE_N) = 100*ZZ + EE + &
                   OP_ISLP(I,J,K)*0.001 + OP_ILV(I,J,K)*0.00001
              IF ( EDGE_N == EDGE_M ) EXIT ATOMS
           ENDIF

        ENDDO
     ENDDO
  ENDDO ATOMS

END SUBROUTINE OP_EDGEFIND
