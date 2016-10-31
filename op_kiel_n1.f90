      SUBROUTINE OP_KIEL_N1 (T,FRE,WLS,ZN1, ANEL)

!**** TOTAL CONTINOUS ABSORPTION COEEFICIENT OF NITROGEN I ,
!**** FROM G.PEACH  MEM.R.A.S. 73,1

!
!      Peach opacities for N
!
      use op_precisions
      use op_peach

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FRE
      real, intent(in) :: WLS
      real :: ZN1
      real, intent(out) :: ANEL

      real :: AN

      real :: TEM(6), TEM1(6), TEM2(6), TEM3(6), TEM4(6), TEM5(6)
      real :: WL(45)
      real :: TAB (6,45)
      integer :: I, J
    
      DATA TEM1 / 10000., 11000., 12000., 13000., 14000., 15000./
      DATA TEM2 / 16000., 17000., 18000., 19000., 20000., 21000./
      DATA TEM3 / 22000., 23000., 24000., 25000., 26000., 27000./
      DATA TEM4 / 28000., 29000., 30000., 32000., 34000., 36000./
      DATA TEM5 / 38000., 40000., 42000., 44000., 46000., 48000./
      DATA  WL / 6300., 5900., 5500., 5092., 5091., 4960., 4855., 4854.,&
               & 4844., 4835., 4834., 4650., 4562., 4561., 4490., 4426.,&
               & 4425., 4800., 4193., 4192., 3800., 3500., 3201., 3200.,&
               & 3050., 2933., 2932., 2480., 2030., 1580., 1129., 1128.,&
               & 1070., 1019., 1018.,  963.,  962.,  882.,  881.,  852.,&
               &  851.,  825.,  824.,  800.,  750. /

      ANEL = 0.0
      AN = 0.0
      IF (WLS .GT. 6300. .OR. WLS .LT. 750.)    GOTO 17
      IF (T .LT. TEM1(1) .OR. T .GT. TEM5(6)) GOTO 17

!**** FIND THE CORRECT TABLE AND THEN INTERPOLATE
      IF (T .GT. TEM1(6)) GOTO 1
      CALL OP_LINPOL(WL,TEM1,WLS,T,TAB1N1,AN,45,6)
      GOTO 17
    1 IF (T .GE. TEM2(1)) GOTO 4
         DO 3 J = 1,45
         DO 2 I = 1,3
              TEM(I) = TEM1(I+3)
    2              TAB(I,J) = TAB1N1(I+3,J)
         DO 3 I = 4,6
              TEM(I) = TEM2(I-3)
    3              TAB(I,J) = TAB2N1(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,45,WLS,T,AN,-1)
      GOTO 17
    4 IF (T .GT. TEM2(6)) GOTO 5
      CALL OP_TABLE(TEM2,WL,TAB2N1,45,WLS,T,AN,-1)
      GOTO 17
    5 IF ( T .GE. TEM3(1)) GOTO 8
         DO 7 J = 1,45
         DO 6 I = 1,3
              TEM(I) = TEM2(I+3)
    6              TAB(I,J) = TAB2N1(I+3,J)
         DO 7 I = 4,6
              TEM(I) = TEM3(I-3)
    7              TAB(I,J) = TAB3N1(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,45,WLS,T,AN,-1)
      GOTO 17
    8 IF (T .GT. TEM3(6)) GOTO 9
      CALL OP_TABLE(TEM3,WL,TAB3N1,45,WLS,T,AN,-1)
      GOTO 17
    9 IF (T .GE. TEM4(1)) GOTO 12
         DO 11 J = 1,45
         DO 10 I = 1,3
              TEM(I) = TEM3(I+3)
   10              TAB(I,J) = TAB3N1(I+3,J)
         DO 11 I = 4,6
              TEM(I) = TEM4(I-3)
   11              TAB(I,J) = TAB4N1(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,45,WLS,T,AN,-1)
      GOTO 17
   12 IF (T .GT.TEM4(6)) GOTO 13
      CALL OP_TABLE (TEM4,WL,TAB4N1,45,WLS,T,AN,-1)
      GOTO 17
   13 IF (T .GE. TEM5(1)) GOTO 16
         DO 15 J= 1,45
         DO 14 I = 1,3
              TEM(I) = TEM4(I+3)
   14              TAB(I,J) = TAB4N1(I+3,J)
         DO 15 I = 4,6
              TEM(I) = TEM5(I-3)
   15              TAB(I,J) = TAB5N1(I-3,J)
      CALL OP_TABLE(TEM,WL,TAB,45,WLS,T,AN,-1)
      GOTO 17
   16 CALL OP_TABLE(TEM5,WL,TAB5N1,45,WLS,T,AN,-1)
   17 ANEL = AN*4./ZN1
      RETURN

      END  subroutine OP_KIEL_N1
