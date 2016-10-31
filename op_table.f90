      SUBROUTINE OP_TABLE (TEM,WL,TAB,NL,WLS,T,ANEL,MM)

!**** AUFSUCHEN DER RICHTIGEN TABELLENWERTE FUER DIE INTERPOLATION
!**** MM = -1, ABSTEIGENDE WERTE IN X
!**** MM = 0, KEINE INTERPOLATION
!**** MM = 1, AUFSTEIGENDE WERTE IN X

      use op_precisions
      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: WLS
      real, intent(inout) :: ANEL
      real :: TEM, WL, TAB
      integer :: NL, MM

      real :: X, Y, Z
      integer :: K1, K2, J, I

      DIMENSION TEM(6),WL(50),TAB(6,50),X(2),Y(4),Z(4,2)
!**** POSITION IN X
      IF (MM) 100,300,200
  100 IF (WLS .LE. WL(1)) GOTO 1
         X(1) = WL(1)
         X(2) = WL(2)
         K1 = 1
         K2 = 2
         GOTO 3
    1 DO 2 J = 2,NL
      IF (WLS .LE. WL(NL)) GOTO 22
         IF (WLS .LT. WL(J)) GOTO 2
         X(1) = WL(J-1)
         X(2) = WL(J)
         K1 = J-1
         K2 = J
         GOTO 3
    2         CONTINUE
   22         X(1) = WL(NL-1)
         X(2) = WL(NL)
         K1 = NL-1
         K2 = NL
         GOTO 3
  200 IF (WLS .GE. WL(1)) GOTO 10
         X(1) = WL(2)
         X(2) = WL(1)
         K1 = 2
         K2 = 1
         GOTO 3
   10 DO 12 J = 2,NL
      IF (WLS .GE. WL(NL)) GOTO 32
         IF (WLS .GT. WL(J)) GOTO 12
         X(1) = WL(J)
         X(2) = WL(J-1)
         K1 = J
         K2 = J-1
         GOTO 3
   12         CONTINUE
   32         X(1) = WL(NL)
         X(2) = WL(NL-1)
         K1 = NL
         K2 = NL-1
!**** POSITION IN Y
    3 IF (T .GT. TEM(2)) GOTO 5
         DO 4 I = 1,4
         Y(I) = TEM(I)
         Z(I,1) = TAB(I,K1)
    4         Z(I,2) = TAB(I,K2)
         GOTO 9
    5 IF (T .GT. TEM(4)) GOTO 7
         DO 6 I = 2,5
         Y(I-1) = TEM(I)
         Z(I-1,1) = TAB(I,K1)
    6         Z(I-1,2) = TAB(I,K2)
         GOTO 9
    7 DO 8 I = 3,6
         Y(I-2) = TEM(I)
         Z(I-2,1) = TAB(I,K1)
    8         Z(I-2,2) = TAB(I,K2)
    9 CALL OP_DINPOL(X,Y,Z,WLS,T,ANEL)
      RETURN
  300 ANEL = 0.
      RETURN

      END SUBROUTINE OP_TABLE
