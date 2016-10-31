      SUBROUTINE OP_KIEL_N3 (T,FRE,WLS,ZN3, ANEL)

!**** TOTAL CONTINOUS ABSORPTION COEEFICIENT OF NITROGEN III ,
!**** FROM G.PEACH  MEM.R.A.S. 73,1

!
!      Peach opacities for C and N
!
      use op_precisions
      use op_peach

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FRE
      real, intent(in) :: WLS
      real :: ZN3
      real, intent(out) :: ANEL

      real :: AN

      real :: TEM(6), TEM1(6), TEM2(6), TEM3(6)
      real :: WL(27)
      real :: TAB (6,27)
      integer :: I, J

      DATA TEM1 / 22000., 23000., 24000., 25000., 26000., 27000. /
      DATA TEM2 / 28000., 29000., 30000., 32000., 34000., 36000. /
      DATA TEM3 / 38000., 40000., 42000., 44000., 46000., 48000. /
      DATA  WL / 2000., 1800., 1605., 1604., 1542., 1541., 1411., 1410.,&
               & 1227., 1226., 1046.,  867.,  866.,  798.,  730.,  729.,&
               &  675.,  620.,  619.,  530.,  441.,  351.,  261.,  260.,&
               &  225.,  200.,  175. /

      ANEL = 0.0
      AN = 0.0
      IF (WLS .GT. 2000.)  GOTO 17
      IF (WLS .LT. 175.)   GOTO 17
      IF (T .LT. TEM(1) .OR. T .GT. TEM3(6)) GOTO 17

!**** FIND THE CORRECT TABLE AND THEN INTERPOLATE
      IF (T .GT. TEM1(6)) GOTO 1
      CALL OP_LINPOL(WL,TEM1,WLS,T,TAB1N3,AN,27,6)
      GOTO 17
    1 IF (T .GE. TEM2(1)) GOTO 4
         DO 3 J = 1,27
         DO 2 I = 1,3
              TEM(I) = TEM1(I+3)
    2              TAB(I,J) = TAB1N3(I+3,J)
         DO 3 I = 4,6
              TEM(I) = TEM2(I-3)
    3              TAB(I,J) = TAB2N3(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,27,WLS,T,AN,-1)
      GOTO 17
    4 IF (T .GT. TEM2(6)) GOTO 5
      CALL OP_TABLE(TEM2,WL,TAB2N3,27,WLS,T,AN,-1)
      GOTO 17
    5 IF ( T .GE. TEM3(1)) GOTO 8
         DO 7 J = 1,27
         DO 6 I = 1,3
              TEM(I) = TEM2(I+3)
    6              TAB(I,J) = TAB2N3(I+3,J)
         DO 7 I = 4,6
              TEM(I) = TEM3(I-3)
    7              TAB(I,J) = TAB3N3(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,27,WLS,T,AN,-1)
      GOTO 17
    8 CALL OP_TABLE(TEM3,WL,TAB3N3,27,WLS,T,AN,-1)
      GOTO 17
 17   ANEL = AN*6./ZN3
      RETURN

      END subroutine OP_KIEL_N3
