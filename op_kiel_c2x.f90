      SUBROUTINE OP_KIEL_C2X (T,FRE,WLS,ZC2, AN9)

!**** TOTAL CONTINUOUS ABSORPTION COEFFICIENT PER ATOM OF CARBON II,
!**** CII(2S2P 3P)NL COEFFICIENTS

      use op_precisions
      use op_peach

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FRE, WLS
      real :: ZC2
      real, intent(out) :: AN9

      real :: ANEL

      real :: TEM, TEM1, TEM2, TEM3, TEM4, TEM5, TEM6
      real :: WL
      real :: TAB

      real :: HH, CC, BKT
      integer :: I, J

      DIMENSION TAB(6,48)

      DIMENSION TEM1(6),TEM2(6),TEM3(6),TEM4(6),TEM5(6),TEM(6),WL(48)
      DATA TEM1 / 10000.,11000.,12000.,13000.,14000.,15000./
      DATA TEM2 / 16000.,17000.,18000.,19000.,20000.,21000./
      DATA TEM3 / 22000.,23000.,24000.,25000.,26000.,27000./
      DATA TEM4 / 28000.,29000.,30000.,32000.,34000.,36000./
      DATA TEM5 / 38000.,40000.,42000.,44000.,46000.,48000./
      DATA WL / 4000.,3700.,3400.,3100.,2800.,2500.,2332.,2133.,2132.,&
              & 2037.,2036.,1992.,1991.,1975.,1974.,1905.,1904.,1877.,&
              & 1876.,1835.,1834.,1654.,1653.,1597.,1596.,1553.,1552.,&
              & 1492.,1491.,1485.,1484.,1412.,1411.,1219.,1218.,1053.,&
              &  888., 723., 722., 656., 655., 575., 574., 486., 485.,&
              &  442., 400., 350./
      DATA  HH,CC/ 4.7993E-11, 2.99793E+18 /

      ANEL = 0.
      BKT = 8.6167E-05*T
      IF (WLS .GT. 4000.) GOTO 17
      IF (WLS .LT. 350.)  GOTO 17
      IF (T .LT. TEM1(1) .OR. T .GT. TEM5(6))  GOTO 17
!**** FIND THE CORRECT TABLE AND THEN INTERPOLATE
      IF (T .GT. TEM1(6)) GOTO 1
      CALL OP_TABLE(TEM1,WL,TAB1C4,48,WLS,T,ANEL,-1)
      GOTO 17
    1 IF (T .GE. TEM2(1)) GOTO 4
         DO 3 J = 1,48
         DO 2 I = 1,3
              TEM(I) = TEM1(I+3)
    2              TAB(I,J) = TAB1C4(I+3,J)
         DO 3 I = 4,6
              TEM(I) = TEM2(I-3)
    3              TAB(I,J) = TAB2C4(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,48,WLS,T,ANEL,-1)
      GOTO 17
    4 IF (T .GT. TEM2(6)) GOTO 5
      CALL OP_TABLE(TEM2,WL,TAB2C4,48,WLS,T,ANEL,-1)
      GOTO 17
    5 IF (T .GE. TEM3(1)) GOTO 8
         DO 7 J = 1,48
         DO 6 I = 1,3
              TEM(I) = TEM2(I+3)
    6              TAB(I,J) = TAB2C4(I+3,J)
         DO 7 I = 4,6
              TEM(I) = TEM3(I-3)
    7              TAB(I,J) = TAB3C4(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,48,WLS,T,ANEL,-1)
      GOTO 17
    8 IF (T .GT. TEM3(6)) GOTO 9
      CALL OP_TABLE(TEM3,WL,TAB3C4,48,WLS,T,ANEL,-1)
      GOTO 17
    9 IF (T .GE. TEM4(1)) GOTO 12
         DO 11 J = 1,48
         DO 10 I = 1,3
              TEM(I) = TEM3(I+3)
   10              TAB(I,J) = TAB3C4(I+3,J)
         DO 11 I = 4,6
              TEM(I) = TEM4(I-3)
   11              TAB(I,J) = TAB4C4(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,48,WLS,T,ANEL,-1)
      GOTO 17
   12 IF (T .GT. TEM4(6)) GOTO 13
      CALL OP_TABLE(TEM4,WL,TAB4C4,48,WLS,T,ANEL,-1)
      GOTO 17
   13 IF (T .GE. TEM5(1)) GOTO 16
         DO 15 J = 1,48
         DO 14 I = 1,3
              TEM(I) = TEM4(I+3)
   14              TAB(I,J) = TAB4C4(I+3,J)
         DO 15 I = 4,6
              TEM(I) = TEM5(I-3)
   15              TAB(I,J) = TAB5C4(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,48,WLS,T,ANEL,-1)
      GOTO 17
   16 CALL OP_TABLE(TEM5,WL,TAB5C4,48,WLS,T,ANEL,-1)
 17   AN9 = ANEL*(12.*EXP(-5.329/BKT)+10.*EXP(-9.284/BKT)+ &
                     &  2.*EXP(-11.956/BKT)+6.)    /ZC2
      AN9 = AN9*2.*EXP(-5.329/BKT)
      RETURN

      END subroutine OP_KIEL_C2X
