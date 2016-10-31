SUBROUTINE OP_ATFILL (ND)

  !! This subroutine needs to be reviwed. Explanations and comments for the 
  !! calculations made are needed. Explicitly state where these values will be 
  !! used, and the significance of the values used.

  !+
  !
  !       Fill array AT with temperature dependent constants
  !       for H,He,He(+) continuous opacity calculations
  !
  !-

  USE LTE_MODEL, ONLY : TEMP
  USE OP_GLOBAL, ONLY : AT
  USE TAPP_PAR

  IMPLICIT NONE

! Local constants
  REAL,PARAMETER :: KEV = TAPP__K / TAPP__EV 

! Subroutine argument (given)
  INTEGER, INTENT(IN) :: ND    ! number of depth zones

! Local variables
  INTEGER :: J
  REAL :: AF

  DO J = 1,ND

     AT(1,J) = EXP(-11.8403E+04/TEMP(J)) /8.0
     AT(2,J) = EXP(-14.0330E+04/TEMP(J)) /27.0
     AT(3,J) = EXP(-14.8004E+04/TEMP(J)) /64.0
     AT(4,J) = EXP(-15.1556E+04/TEMP(J)) /125.0
     AT(5,J) = EXP(-15.3486E+04/TEMP(J)) /216.0
     AT(6,J) = EXP(-15.4649E+04/TEMP(J))*TEMP(J)/31.5742E+04

     AF = KEV * TEMP(J)
     AT(7,J)  =  3.*EXP(-19.816/AF)*10.**(-17.387)
     AT(8,J)  =     EXP(-20.613/AF)*10.**(-15.09)
     AT(9,J)  =  9.*EXP(-20.961/AF)
     AT(10,J) =  3.*EXP(-21.215/AF)
     AT(11,J) =  3.*EXP(-22.716/AF)*10.**(-16.05)
     AT(12,J) =     EXP(-22.918/AF)*10.**(-15.68)
     AT(13,J) =  9.*EXP(-23.005/AF)*10.**(-14.99)
     AT(14,J) = 20.*EXP(-23.072/AF)*10.**(-14.66)
     AT(15,J) =  3.*EXP(-23.085/AF)*10.**(-14.92)
     AT(16,J) = 1.5508E-13*EXP(-23.723/AF)*AF/13.60*10.**0.55

     AT(17,J) = EXP(-47.3580E+04/TEMP(J))/8.0
     AT(18,J) = EXP(-56.1284E+04/TEMP(J))/27.0
     AT(19,J) = EXP(-59.1977E+04/TEMP(J))/64.0
     AT(20,J) = EXP(-60.6164E+04/TEMP(J))/125.0
     AT(21,J) = EXP(-61.3902E+04/TEMP(J))/216.0
     AT(22,J) = EXP(-61.8556E+04/TEMP(J))

  ENDDO

END SUBROUTINE OP_ATFILL
