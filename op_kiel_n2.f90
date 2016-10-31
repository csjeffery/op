      SUBROUTINE OP_KIEL_N2 (T,FRE,WLS,ZN2, ANEL)

!**** TOTAL CONTINOUS ABSORPTION COEEFICIENT OF NITROGEN II,
!**** FROM G.PEACH  MEM.R.A.S. 73,1

!
!        Peach opacities for C and N
!
      use op_precisions
      use op_peach

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FRE
      real, intent(in) :: WLS
      real :: ZN2
      real, intent(out) :: ANEL

      real :: AN

      real :: TEM(6), TEM1(6), TEM2(6), TEM3(6), TEM4(6)
      real :: WL(40)
      real :: TAB (6,40)
      integer :: I, J

      DATA TEM1 / 16000., 17000., 18000., 19000., 20000., 21000./
      DATA TEM2 / 22000., 23000., 24000., 25000., 26000., 27000./
      DATA TEM3 / 28000., 29000., 30000., 32000., 34000., 36000./
      DATA TEM4 / 38000., 40000., 42000., 44000., 46000., 48000./
      DATA  WL / 2400., 2052., 2051., 2019., 2018., 2001., 2000., 1946.,&
               & 1945., 1932., 1931., 1913., 1912., 1651., 1650., 1547.,&
               & 1546., 1466., 1465., 1430., 1429., 1385., 1384., 1347.,&
               & 1346., 1116., 1115., 1113., 1112.,  903.,  694.,  485.,&
               &  484.,  447.,  446.,  419.,  418.,  350.,  300.,  250./

      ANEL = 0.0
      AN = 0.0
      IF (WLS .GT. 2400.)  GOTO 17
      IF (WLS .LT.  250.)  GOTO 17
      IF (T .LT. TEM1(1) .OR. T .GT. TEM4(6))  GOTO 17

!**** FIND THE CORRECT TABLE AND THEN INTERPOLATE
      IF (T .GT. TEM1(6)) GOTO 1
      CALL OP_LINPOL(WL,TEM1,WLS,T,TAB1N2,AN,40,6)
      GOTO 17
    1 IF (T .GE. TEM2(1)) GOTO 4
         DO 3 J = 1,40
         DO 2 I = 1,3
              TEM(I) = TEM1(I+3)
    2              TAB(I,J) = TAB1N2(I+3,J)
         DO 3 I = 4,6
              TEM(I) = TEM2(I-3)
    3              TAB(I,J) = TAB2N2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,40,WLS,T,AN,-1)
      GOTO 17
    4 IF (T .GT. TEM2(6)) GOTO 5
      CALL OP_TABLE(TEM2,WL,TAB2N2,40,WLS,T,AN,-1)
      GOTO 17
    5 IF ( T .GE. TEM3(1)) GOTO 8
         DO 7 J = 1,40
         DO 6 I = 1,3
              TEM(I) = TEM2(I+3)
    6              TAB(I,J) = TAB2N2(I+3,J)
         DO 7 I = 4,6
              TEM(I) = TEM3(I-3)
    7              TAB(I,J) = TAB3N2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,40,WLS,T,AN,-1)
      GOTO 17
    8 IF (T .GT. TEM3(6)) GOTO 9
      CALL OP_TABLE(TEM3,WL,TAB3N2,40,WLS,T,AN,-1)
      GOTO 17
    9 IF (T .GE. TEM4(1)) GOTO 12
         DO 11 J = 1,40
         DO 10 I = 1,3
              TEM(I) = TEM3(I+3)
   10              TAB(I,J) = TAB3N2(I+3,J)
         DO 11 I = 4,6
              TEM(I) = TEM4(I-3)
   11              TAB(I,J) = TAB4N2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,40,WLS,T,AN,-1)
      GOTO 17
   12 CALL OP_TABLE (TEM4,WL,TAB4N2,40,WLS,T,AN,-1)
      GOTO 17
   17 ANEL = AN*9./ZN2
      RETURN

      END subroutine OP_KIEL_N2
