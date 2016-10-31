      SUBROUTINE OP_KIEL_HM (T,FRE,PE,ZH1, AN2)

!     TOTAL CONTINUOUS ABSORPTION COEFFICIENT PER ATOM OF NEGATIVE
!     HYDROGEN ION

      USE OP_PRECISIONS

      IMPLICIT NONE

      REAL, INTENT(IN) :: T, PE
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL, INTENT(OUT) :: AN2
      REAL :: ZH1

      REAL :: HH, THET, AN21, AN22
      REAL(KIND=QQ) :: F, A, B

      DATA  HH / 4.7993E-11 /

      F = FRE/1.0000E+14
      THET = 5040.2/T
      AN2 = 0.
      IF (2.11122E+14-FRE) 2,2,1
    1 IF (1.82589E+14-FRE) 3,3,5

!     BOUND FREE ABSORPTION
    2 A = 29.9793/F
      B = 0.00680133 + 0.178708*A + 0.164790*A*A - 0.0204842*A*A*A + &
              & 5.95244E-04*A**4

      GOTO 4
    3 A = 16.419-29.9793/F
      B = 0.269818*A+0.220190*A**2-0.0411288*A**3+0.00273236*A**4
    4 AN21 = 0.4158E-26*PE*THET**2.5*EXP(1.726*THET)*(1.- &
                &  EXP(-(HH*FRE/T)))*B
      AN2 = AN2+MAX(AN21,0.)

!     FREE FREE ABSORPTION
    5 AN22 = 1.0000E-26*PE*(0.0053666-0.011493*THET+0.027029*THET**2- &
                &  (3.2062-11.924*THET+5.9390*THET**2)*2.99793E-02/F- &
                &  (0.40192-7.0355*THET+0.34592*THET**2)*0.8987584/F**2)
      AN2 = AN2+MAX(AN22,0.)
      AN2 = AN2*2./ZH1

      !  write (*,*) 'op_kiel_hm: ',fre,f,thet,a,b,zh1,an21,an22,an2

  END SUBROUTINE OP_KIEL_HM
