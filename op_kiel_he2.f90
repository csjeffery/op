      SUBROUTINE OP_KIEL_HE2 (J,T,FRE,ZHE2, AN4)

!     TOTAL CONTINUOUS ABSORPTION COEFFICIENT PER ATOM OF HELIUM II
!     THE GAUNT FACTORS ARE TAKEN FROM GINGERICH (1964)

      USE OP_PRECISIONS
      USE OP_GLOBAL

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: J
      REAL, INTENT(IN) :: T
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL,INTENT(IN) :: ZHE2

      REAL, INTENT(OUT) :: AN4

!  Local constants.
      REAL, PARAMETER :: HH = 4.79926E-11 
      REAL, PARAMETER :: RC = 3.28940E+15 

!  Local variables.
      REAL ::  F, AF
      REAL(KIND=QQ) :: FRE2
      !real :: AN41, AN42, AN43, AN44, AN45, AN46, AN47
      REAL, DIMENSION(7) :: AN4P(7)
      REAL :: GN1, GN2, GN3, GNFF 

!.

      F = FRE/1.0000E+14
      AF = 4.5056E-12/F**3
      FRE2 = FRE*FRE
      AN4P = 0.
      AN4 = 0.

      IF (1.3158E+16-FRE) 6,6,1
    1 IF (3.2894E+15-FRE) 7,7,2
    2 IF (1.4620E+15-FRE) 8,8,3
    3 IF (8.2235E+14-FRE) 9,9,4
    4 IF (5.2630E+14-FRE) 10,10,5
    5 IF (3.6549E+14-FRE) 11,11,12

    6 GN1     = 0.9916+1.087408E+14/FRE-3.62955E+31/FRE2
      AN4P(1) = AF*GN1
      AN4     = AN4+AN4P(1)
    7 GN2     = 1.105-9.49984E+14/FRE+6.52288E+29/FRE2
      AN4P(2) = AF*GN2*AT(17,J)
      AN4     = AN4+AN4P(2)
    8 GN3     = 1.101- 3.94528E+14/FRE+1.65659E+29/FRE2
      AN4P(3) = AF*GN3*AT(18,J)
      AN4     = AN4+AN4P(3)
    9 AN4P(4) = AF*AT(19,J)
      AN4     = AN4+AN4P(4)
   10 AN4P(5) = AF*AT(20,J)
      AN4     = AN4+AN4P(5)
   11 AN4P(6) = AF*AT(21,J)
      AN4     = AN4+AN4P(6)
   12 GNFF    = 1.+(0.1728*(FRE/(4.*RC))**0.33333333*(1.+(2.*T)/(HH*FRE)))
      AN4P(7) = AT(22,J) * GNFF/F**3 * T*3.5677E-18
      AN4     = AN4+AN4P(7)

      AN4     = AN4*2./ZHE2

!     write (*,*) 'op_kiel_he2: ',fre,f,af,an4p,an4 

  END SUBROUTINE OP_KIEL_HE2
