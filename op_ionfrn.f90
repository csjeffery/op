SUBROUTINE OP_IONFRN (ND)

  !+
  !
  !  Name:
  !  OP_IONFRN
  !
  !  Purpose: 
  !  Compute ionisation fractions throughout atmosphere using
  !  model atmosphere in LTE_MODEL
  !  And stores results in OP_GLOBAL (OP_XUU, OP_XION)
  !
  !-

  !  Global variables.
  USE LTE_MODEL, ONLY : TEMP, P_TOTAL, P_E
  USE OP_GLOBAL, ONLY : OP_MAXNA, OP_MAXNZ, OP_MAXNI
  USE OP_GLOBAL, ONLY : OP_XUU, OP_XION, OP_ALPNN

  IMPLICIT NONE

  !  Subroutine Arguments (Given).
  INTEGER, INTENT(IN) :: ND

  !  Local constants.
  REAL, PARAMETER ::  S0=0.666698474
  REAL, PARAMETER ::  BKX = 1./8.6167E-05

  INTEGER :: J, I, K
  REAL :: TJ, PGJ, PEJ
  REAL :: SAHA, BKDT, PE2, PE3

  REAL :: FH0, FHM, FHEM, FHE1, FHE2, GHE1, GHE2, GHEM

  ! External function
  REAL, EXTERNAL :: U

  !.

  !     Loop through atmosphere

  DO J = 1,ND
     TJ  = TEMP(J)
     PGJ = P_TOTAL(J)
     PEJ = P_E(J)

     SAHA = S0  * TJ**2.5
     BKDT = BKX / TJ
     PE2  = PEJ**2
     PE3  = PEJ**3

     !     Hydrogen
     OP_XUU(1,1,J) = U  (1,1,TJ,PEJ)
     OP_XUU(2,1,J) = U  (1,2,TJ,PEJ)
     FH0 = SAHA * EXP(-13.598*BKDT) / OP_XUU(1,1,J) / PEJ
     OP_XION(1,1,J) = 1. / (1.+FH0)
     OP_XION(2,1,J) = 1. - OP_XION(1,1,J) 

     FHM = SAHA * EXP(- 0.75*BKDT) * OP_XUU(1,1,J) / PEJ
     OP_XION(3,1,J) = OP_XION(1,1,J) / FHM 
     OP_XION(2,1,J) = 1. - OP_XION(1,1,J) - OP_XION(3,1,J)

     !        IF (J == 1) write (*,*) 'ionfrn: t,p,pe',TJ,PGJ,PEJ
     !        IF (J == 1) write (*,*) 'ionfrn: u',OP_XUU(1,1,J),  OP_XUU(2,1,J),  OP_XUU(3,1,J) 
     !        IF (J == 1) write (*,*) 'ionfrn: x',OP_XION(1,1,J), OP_XION(2,1,J), OP_XION(3,1,J) 
     !        IF (J == 1) write (*,*) 'ionfrn: f',FH0, FHM         

     !     Helium
     OP_XUU(1,2,J) = U  (2,1,TJ,PEJ)
     OP_XUU(2,2,J) = U  (2,2,TJ,PEJ)

     FHEM = SAHA * EXP(+ 0.075*BKDT)
     FHE1 = SAHA * EXP(-24.587*BKDT)
     FHE2 = SAHA * EXP(-54.416*BKDT)

     GHE1 = FHE1 * OP_XUU(2,2,J) / OP_XUU(1,2,J)
     GHE2 = FHE2 / OP_XUU(2,2,J)
     GHEM = FHEM * OP_XUU(1,2,J) / 6.

     OP_XION(2,2,J) = PEJ*GHE1*GHEM / &
          &(PE3 + PE2*GHEM + PEJ*GHE1*GHEM + GHE1*GHE2*GHEM)
     OP_XION(3,2,J) = OP_XION(2,2,J) * GHE2/PEJ
     OP_XION(4,2,J) = OP_XION(2,2,J)*PE2 / (GHE1*GHEM)
     !         OP_XION(1,2,J) = 1.0 - OP_XION(2,2,J) - OP_XION(3,2,J) - OP_XION(4,2,J)!OP_XION(2,2,J)*PEJ/GHE1
     OP_XION(1,2,J) = OP_XION(2,2,J)*PEJ/GHE1



     !  Z > 2
     DO I = 3, OP_MAXNZ  !atoms
     !   IF (OP_ALPNN(I,J) > 0.0) THEN

           !  partition functions
           DO K = 1, OP_MAXNI !loop over ions
              OP_XUU(K,I,J) = U (I,K,TJ,PEJ)
           ENDDO

           !  ionisation fractions
           CALL OP_ION (I,TJ,PEJ,OP_XUU(1:OP_MAXNI,I,J),OP_XION(1:OP_MAXNI,I,J)) 

     !   ENDIF
     ! IF ( J == 1 ) write (*,*) j, i, op_maxni, op_xion(1:op_maxni,i,j)
     ENDDO

     !  level populations
     CALL OP_POP(J,TJ)

  ENDDO

END SUBROUTINE OP_IONFRN
