SUBROUTINE OP_POP(IL,TL)

  !+

  !  Name:
  !     OP_POP

  !  Purpose:
  !     To calculate fractional populations for individual states
  !     for all species

  !  Module Type:
  !     Subroutine

  !  Language:
  !     Fortran 90

  !  Description:
  !     Note: individual states are not necessarily identified with a unique
  !     energy level, but also with a given slp state.

  !  Input arguments:
  !     IL        current depth point
  !     TL        current temperature
  !
  !  Authors:
  !     UNK: unknown
  !     PMH: P.M.Harrison (St.Andrews)
  !
  !  History:
  !     30-Nov-1994 (PMH): Modified for Sterne 3 from Sterne 2.1 routine PRESS
  !     09-Feb-1998 (PMH): Hacked from stn_presel.
  !     26-Feb-2003 (CSJ): Converted to F90.
  !
  !.

  !  Global OP library variables
  USE OP_GLOBAL
  USE TAPP_PAR

  IMPLICIT NONE

  !  Subroutine Arguments (Given).
  INTEGER :: IL
  REAL ::  TL

  !  Local variables.
  INTEGER :: I,J,K
  REAL :: BKDT

  !-

  BKDT =  (TAPP__EV/TAPP__K) / TL

  !  Calculate fractional populations for individual states:
  DO K=1,OP_MAXNA                         ! loop over species
     IF (OP_ALPIMIN(K).NE.0) THEN
        DO I=OP_ALPIMIN(K),OP_ALPIMAX(K)         ! loop over ions
           DO J=1,OP_NS(I,K)                      ! loop over states
              IF (OP_ALPNN(K,IL) == 0.0) THEN
                 OP_XPOP(J,I,K,IL) = 0.0
              ELSE
                 OP_XPOP(J,I,K,IL)= OP_ALPNN(K,IL)*OP_XION(I,K,IL) &
                      *OP_XWT(J,I,K)* &
                      EXP(-OP_XCHI(J,I,K)*BKDT)/OP_XUU(I,K,IL)
              END IF
           ENDDO
        ENDDO
     ENDIF
  ENDDO

END SUBROUTINE OP_POP

