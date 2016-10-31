      SUBROUTINE OP_FFCM (T,FRE,PE,WLS,ZC1,ANEL)

! Free-free and bound-free absorption coefficients for C minus, 
! including stimutaled emission.
!*****        TOTALER KONTINUIERLICHER ABSORPTIONSKOEFFIZIENT VON C MINUS
!*****        MYERSCOUGH, MC DOWELL, MON.NOT.ROY.ASTR.SOC., 132, 451 (1966)

      use op_precisions

      implicit none

      real, intent(in) :: T, PE
      real, intent(in) :: wls
      real(kind=qq), intent(in) :: FRE
      real, intent(in) :: ZC1
      real, intent(out) :: ANEL

      real(kind=qq) :: AN6
      real, dimension (7,13) :: TAB
      real :: TEM, WL

      DIMENSION  TEM(7),WL(13)
      DATA  TEM/ 4000., 5000., 6000., 7000., 8000., 9000., 10000. /
      DATA  WL/ 22782., 18225., 15188., 13018., 11391., 10125., 9113.,&
               & 4556.,  3037.,  1822.,  1302.,  1012.,   911. /
      DATA  TAB/ 9.39, 5.30, 3.25, 2.12, 1.45, 1.03, .751,&
              &  6.68, 3.84, 2.39, 1.58, 1.09, .781, .576,&
              &  5.12, 1.22, .853, .621, .462, .462, .467,&
              &  4.19, 2.46, 1.57, 1.06, .742, .539, .403,&
              &  3.30, 2.07, 1.32, .894, .643, .462, .351,&
              &  3.03, 1.79, 1.15, .782, .565, .413, .318,&
              & 10.70, 3.81, 1.90, 1.12, .732, .509, .370,&
              & 34.38, 9.87, 4.02, 2.01, 1.16, .729, .491,&
              & 34.34, 9.69, 3.86, 1.89, 1.07, .661, .439,&
              & 32.68, 9.16, 3.62, 1.76, .984, .606, .400,&
              & 28.75, 8.03, 3.16, 1.53, .853, .523, .344,&
              & 24.62, 6.87, 2.70, 1.31, .727, .445, .292,&
              & 22.68, 6.33, 2.49, 1.20, .668, .408, .268 /

      ANEL =0.0
      AN6 = 0.0
      IF (T .GT. TEM(7))  RETURN
      CALL OP_LINPOL (WL,TEM,WLS,T,TAB,AN6,13,7)
      AN6 = AN6*PE*9./ZC1*1.E-26
!***** STIMULIERTE EMISSION MUSS RUECKGAENGIG GEMACHT WERDEN, SIEHE OPA
      ANEL = AN6/(1.- EXP(-4.7993E-11*FRE/T))

      END subroutine OP_FFCM
