      SUBROUTINE OP_KIEL_C2 (T,FRE,WLS,ZC2, AN7)

!**** TOTAL CONTINUOUS ABSORPTION COEFFICIENT PER ATOM OF CARBON II

      use op_precisions
      use op_peach

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FRE, WLS
      real :: ZC2
      real, intent(out) :: AN7
      real :: ANEL

      real :: TEM, TEM1, TEM2, TEM3, TEM4, TEM5
      real :: WL
      real :: TAB

      real :: HH, CC
      integer :: I, J

      DIMENSION TAB(6,33)

      DIMENSION TEM1(6),TEM2(6),TEM3(6),TEM4(6),TEM5(6),TEM(6),WL(33)
      DATA TEM1 /10000.,11000.,12000.,13000.,14000.,15000./
      DATA TEM2 /16000.,17000.,18000.,19000.,20000.,21000./
      DATA TEM3 /22000.,23000.,24000.,25000.,26000.,27000./
      DATA TEM4 /28000.,29000.,30000.,32000.,34000.,36000./
      DATA TEM5 /38000.,40000.,42000.,44000.,46000.,48000./
      DATA WL /4200.,4000.,3800.,3613.,3612.,3550.,3505.,3504.,&
             & 3300.,3100.,2930.,2929.,2730.,2537.,2536.,2340.,&
             & 2150.,1957.,1956.,1750.,1541.,1540.,1410.,1380.,&
             & 1249.,1248.,1100.,900.,700.,509.,508.,450.,400./
      DATA  HH,CC/ 4.7993E-11, 2.99793E+18 /

      ANEL = 0.
      IF (WLS .GT. 4200. .OR. WLS .LT. 400.) GOTO 17
      IF (T .LT. TEM1(1) .OR. T .GT. TEM5(6))  GOTO 17
!**** FIND THE CORRECT TABLE AND THEN INTERPOLATE
      IF (T .GT. TEM1(6)) GOTO 1
      CALL OP_TABLE(TEM1,WL,TAB1C2,33,WLS,T,ANEL,-1)
      GOTO 17
    1 IF (T .GE. TEM2(1)) GOTO 4
         DO 3 J = 1,33
         DO 2 I = 1,3
              TEM(I) = TEM1(I+3)
    2              TAB(I,J) = TAB1C2(I+3,J)
         DO 3 I = 4,6
              TEM(I) = TEM2(I-3)
    3              TAB(I,J) = TAB2C2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,33,WLS,T,ANEL,-1)
      GOTO 17
    4 IF (T .GT. TEM2(6)) GOTO 5
      CALL OP_TABLE(TEM2,WL,TAB2C2,33,WLS,T,ANEL,-1)
      GOTO 17
    5 IF (T .GE. TEM3(1)) GOTO 8
         DO 7 J = 1,33
         DO 6 I = 1,3
              TEM(I) = TEM2(I+3)
    6              TAB(I,J) = TAB2C2(I+3,J)
         DO 7 I = 4,6
              TEM(I) = TEM3(I-3)
    7              TAB(I,J) = TAB3C2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,33,WLS,T,ANEL,-1)
      GOTO 17
    8 IF (T .GT. TEM3(6)) GOTO 9
      CALL OP_TABLE(TEM3,WL,TAB3C2,33,WLS,T,ANEL,-1)
      GOTO 17
    9 IF (T .GE. TEM4(1)) GOTO 12
         DO 11 J = 1,33
         DO 10 I = 1,3
              TEM(I) = TEM3(I+3)
   10              TAB(I,J) = TAB3C2(I+3,J)
         DO 11 I = 4,6
              TEM(I) = TEM4(I-3)
   11              TAB(I,J) = TAB4C2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,33,WLS,T,ANEL,-1)
      GOTO 17
   12 IF (T .GT. TEM4(6)) GOTO 13
      CALL OP_TABLE(TEM4,WL,TAB4C2,33,WLS,T,ANEL,-1)
      GOTO 17
   13 IF (T .GE. TEM5(1)) GOTO 16
         DO 15 J = 1,33
         DO 14 I = 1,3
              TEM(I) = TEM4(I+3)
   14              TAB(I,J) = TAB4C2(I+3,J)
         DO 15 I = 4,6
              TEM(I) = TEM5(I-3)
   15              TAB(I,J) = TAB5C2(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,33,WLS,T,ANEL,-1)
      GOTO 17
   16 CALL OP_TABLE(TEM5,WL,TAB5C2,33,WLS,T,ANEL,-1)
 17   AN7 = ANEL*6./ZC2
      RETURN

      END SUBROUTINE OP_KIEL_C2

