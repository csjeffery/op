SUBROUTINE OP_KIEL_C1(T,FRE,WLS,ZC1,AN6)

  !**** TOTAL CONTINUOUS ABSORPTION COEFFICIENT PER ATOM OF CARBON I

  use op_precisions
  use op_peach

  implicit none

  real, intent(in) :: T
  real(kind=qq), intent(in) :: FRE
  real, intent(in) :: WLS
  real :: ZC1
  real, intent(out) :: AN6
  real :: ANEL

  real :: TEM, TEM1, TEM2, TEM3, TEM4, TEM5
  real :: WL
  real :: TAB

  real :: HH, CC, U
  integer :: I, J

  DIMENSION TAB(6,37)

  DIMENSION TEM1(6),TEM2(6),TEM3(6),TEM4(6),TEM5(6),TEM(6),WL(37)
  DATA TEM1 /  4000., 5000., 6000., 7000., 8000., 9000./
  DATA TEM2 / 10000.,11000.,12000.,13000.,14000.,15000./
  DATA TEM3 / 16000.,17000.,18000.,19000.,20000.,21000./
  DATA TEM4 / 22000.,23000.,24000.,25000.,26000.,27000./
  DATA TEM5 / 28000.,29000.,30000.,32000.,34000.,36000./
  DATA WL / 7100.,6700.,6300.,5922.,5921.,5700.,5479.,5478.,5305.,&
       & 5132.,5131.,4971.,4970.,4730.,4729.,4544.,4543.,4183.,3823.,&
       & 3463.,3462.,3280.,3279.,2912.,2545.,2178.,1811.,1445.,1444.,&
       & 1240.,1239.,1101.,1100., 900., 700., 600., 500./
  DATA  HH,CC/ 4.7993E-11, 2.99793E+18 /

  ANEL = 0.
  IF (WLS > 10000.) THEN
     U = 4.79927E-11*FRE/T
     !****  ABSORPTIONSKOEFFIZIENT WASSERSTOFFAEHNLICH, VERGL. UNSOELD
     AN6 = 9.8558E-9*EXP(-1.3075E+5/T)/T/T*EXP(U)/U/U/U*0.83
     AN6 = AN6*9./ZC1
  ELSE

     IF (WLS >= 500.) THEN
        !**** FIND THE CORRECT TABLE AND THEN INTERPOLATE
        IF (T < TEM1(6)) THEN
           CALL OP_LINPOL(WL,TEM1,WLS,T,TAB1C1,ANEL,37,6)

        ELSEIF (T <= TEM2(1)) THEN
           DO  J = 1,37
              DO  I = 1,3
                 TEM(I) = TEM1(I+3)
                 TAB(I,J) = TAB1C1(I+3,J)
              END DO
              DO  I = 4,6
                 TEM(I) = TEM2(I-3)
                 TAB(I,J) = TAB2C1(I-3,J)
              ENDDO
           ENDDO
           CALL OP_TABLE(TEM,WL,TAB,37,WLS,T,ANEL,-1)

        ELSEIF (T < TEM2(6)) THEN
           CALL OP_TABLE(TEM2,WL,TAB2C1,37,WLS,T,ANEL,-1)

        ELSEIF (T < TEM3(1)) THEN
           DO J = 1,37
              DO I = 1,3
                 TEM(I) = TEM2(I+3)
                 TAB(I,J) = TAB2C1(I+3,J)
              ENDDO
              DO I = 4,6
                 TEM(I) = TEM3(I-3)
                 TAB(I,J) = TAB3C1(I-3,J)
              ENDDO
           ENDDO
           CALL OP_TABLE(TEM,WL,TAB,37,WLS,T,ANEL,-1)

        ELSEIF (T < TEM3(6)) THEN
           CALL OP_TABLE(TEM3,WL,TAB3C1,37,WLS,T,ANEL,-1)

        ELSEIF (T < TEM4(1)) THEN
           DO J = 1,37
              DO I = 1,3
                 TEM(I) = TEM3(I+3)
                 TAB(I,J) = TAB3C1(I+3,J)
              ENDDO
              DO I = 4,6
                 TEM(I) = TEM4(I-3)
                 TAB(I,J) = TAB4C1(I-3,J)
              ENDDO
           ENDDO
           CALL OP_TABLE(TEM,WL,TAB,37,WLS,T,ANEL,-1)

        ELSEIF (T < TEM4(6)) THEN
           CALL OP_TABLE(TEM4,WL,TAB4C1,37,WLS,T,ANEL,-1)

        ELSEIF (T < TEM5(1)) THEN
           DO  J= 1,37
              DO  I = 1,3
                 TEM(I) = TEM4(I+3)
                 TAB(I,J) = TAB4C1(I+3,J)
              ENDDO
              DO  I = 4,6
                 TEM(I) = TEM5(I-3)
                 TAB(I,J) = TAB5C1(I-3,J)
              ENDDO
           ENDDO
           CALL OP_TABLE(TEM,WL,TAB,37,WLS,T,ANEL,-1)
        ELSE
           CALL OP_TABLE(TEM5,WL,TAB5C1,37,WLS,T,ANEL,-1)
        ENDIF

     ENDIF
     AN6 = ANEL*9./ZC1

  ENDIF

END SUBROUTINE OP_KIEL_C1

