      SUBROUTINE OP_KIEL_H1 (J,T,FRE,ZH1, AN1)

!+
!
!  opk_oph1
!
!
!  Total continuous absorption coefficient per atom of neutral hydrogen
!
!-

      USE OP_PRECISIONS
      USE OP_GLOBAL, ONLY : AT

      IMPLICIT NONE

      REAL, INTENT(IN) :: T
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL, INTENT(OUT) :: AN1
      REAL :: ZH1  

      INTEGER, INTENT(IN) :: J

!  Local constants.
      REAL, PARAMETER :: HH = 4.7993E-11 
      REAL, PARAMETER :: RC = 3.2894E+15 

      REAL(KIND=QQ) :: F, FRE2, GN1, GN2, GN3, GNFF, AF

!.
      F = FRE/1.0000E+14

!  Gaunt factors are taken from Gingerich (1964)
      FRE2 = FRE*FRE
      GN1 = 0.9916+2.719E+13/FRE-2.2685E+30/FRE2
      GN2 = 1.105 -2.375E+14/FRE+4.0770E+28/FRE2
      GN3 = 1.101 -9.863E+13/FRE+1.0350E+28/FRE2
      GNFF = 1.+(0.1728*(FRE/RC)**0.33333333)*(1.+(2.*T)/(HH*FRE))

      AN1 = 0.
      AF = 2.81592E-13/F**3
      IF (FRE > 3.28806E+15) THEN
        AN1  = AN1+GN1*AF
        AN1  = AN1+GN2*AF*AT(1,J)
        AN1  = AN1+GN3*AF*AT(2,J)
        AN1  = AN1+AF*AT(3,J)
        AN1  = AN1+AF*AT(4,J)
        AN1  = AN1+AF*AT(5,J)
      ELSEIF (FRE > 0.82225E+15) THEN
        AN1  = AN1+GN2*AF*AT(1,J)
        AN1  = AN1+GN3*AF*AT(2,J)
        AN1  = AN1+AF*AT(3,J)
        AN1  = AN1+AF*AT(4,J)
        AN1  = AN1+AF*AT(5,J)
      ELSEIF (FRE > 0.36554E+15)THEN
        AN1  = AN1+GN3*AF*AT(2,J)
        AN1  = AN1+AF*AT(3,J)
        AN1  = AN1+AF*AT(4,J)
        AN1  = AN1+AF*AT(5,J)
      ELSEIF (FRE > 0.20562E+15) THEN
        AN1  = AN1+AF*AT(3,J)
        AN1  = AN1+AF*AT(4,J)
        AN1  = AN1+AF*AT(5,J)
      ELSEIF (FRE > 0.13159E+15) THEN
        AN1  = AN1+AF*AT(4,J)
        AN1  = AN1+AF*AT(5,J)
      ELSEIF (FRE > 0.09134E+15) THEN
        AN1  = AN1+AF*AT(5,J)
      ENDIF

! free-free
      AN1  = AN1 + GNFF*AF*AT(6,J)

      AN1  = AN1*2./ZH1

!     write (*,*) 'op_kiel_h1: ',f,gn1,gn2,gn3,gnff,af,at(6,j),an1

  END SUBROUTINE OP_KIEL_H1
