      SUBROUTINE OP_KIEL_HE1 (J,T,FRE,ZHE1, AN3)

!     TOTAL CONTINUOUS ABSORPTION COEFFICIENT PER ATOM OF HELIUM I

      USE OP_PRECISIONS
      USE OP_GLOBAL

      IMPLICIT NONE

      REAL, INTENT(IN) :: T
      REAL(KIND=QQ), INTENT(IN) :: FRE
      REAL, INTENT(OUT) :: AN3
      REAL :: ZHE1
      INTEGER, INTENT(IN) :: J

      REAL :: HH = 4.79926E-11
      REAL :: RC = 3.28940E+15
      REAL :: F, WLOG, WLOG2
      INTEGER :: IZG = 0
      INTEGER :: K
      INTEGER, SAVE :: IFG, KANTE

      REAL, DIMENSION(10) :: FREK
      REAL, SAVE, DIMENSION(11) :: AN3P = 0.

      DATA FREK /5.9447E15, 1.1526E15, 9.6026E14, 8.7607E14, 8.1465E14,&
                & 4.5191E14, 4.0310E14, 3.8212E14, 3.6600E14, 3.6278E14/

!**** HE1-PHOTOIONISATIONSQUERSCHNITTE NACH V.L.JACOBS FUER DIE
!**** ABSORPTIONSKANTEN:  LAMBDA=  2601, 3122, 3422, 3680 A
!**** PHYS.REV.A3,289(1971)/PHYS.REV.A4,939(1972)/PHYS.REV.A9,1938(1938)

      IFG = 0
      IF (IZG /= 0) GOTO 110
      IF (IFG /= 0) WRITE(6,1010)
      IZG = 1
110   CONTINUE

      KANTE = 1
      DO K=1,10
         IF (FREK(K) > FRE) KANTE = K + 1
      ENDDO

      F = FRE/1.E14
      WLOG = LOG(2.99793E18/FRE)
      WLOG2 = WLOG*WLOG

      GOTO (40,41,42,43,44,45,46,47,48,49,50) ,KANTE
40    AN3P(1) = 2.9512E-14/F/F
41    IF (IFG.EQ.0) GOTO 411
      AN3P(2) = F**0.679*10.**(-0.727*LOG10(F)**2)
      GOTO 42
411   AN3P(2) = -33.72341 + 8.21263*WLOG - 0.47310*WLOG2
      AN3P(2) = EXP(AN3P(2)) * 0.243781
42    IF (IFG.EQ.0) GOTO 421
      AN3P(3) = F**(-1.91)
      GOTO 43
421   AN3P(3) = -25.16025 + 5.19756*WLOG - 0.22242*WLOG2
      AN3P(3) = EXP(AN3P(3)) * 1.23027E-3
43    IF (IFG.EQ.0) GOTO 431
      AN3P(4) = 10.**(-14.03)/F**2.9 + 10.**(-15.14)/F**3.3
      GOTO 44
431   AN3P(4) = -28.26795 + 4.85866*WLOG - 0.12895*WLOG2
      AN3P(4) = EXP(AN3P(4)) * 1.0E-18
44    IF (IFG.EQ.0) GOTO 441
      AN3P(5) = 10.**(-13.69)/F**3.5 + 10.**(-14.92)/F**3.6
      GOTO 45
441   AN3P(5) = -53.49472 +10.82112*WLOG - 0.48589*WLOG2
      AN3P(5) = EXP(AN3P(5)) * 1.0E-18
45    AN3P(6) = F**(-1.54)
46    AN3P(7) = F**(-1.86)
47    AN3P(8) = F**(-2.60)
48    AN3P(9) = F**(-3.69)
49    AN3P(10)= F**(-2.89)
50    AN3P(11)= F**(-3.04)

      AN3 = 0.0
      GOTO (10,11,12,13,14,15,16,17,18,19,20) ,KANTE
10    AN3 = AN3 + AN3P(1)
11    AN3 = AN3 + AN3P(2) * AT(7,J)
12    AN3 = AN3 + AN3P(3) * AT(8,J)
13    AN3 = AN3 + AN3P(4) * AT(9,J)
14    AN3 = AN3 + AN3P(5) * AT(10,J)
15    AN3 = AN3 + AN3P(6) * AT(11,J)
16    AN3 = AN3 + AN3P(7) * AT(12,J)
17    AN3 = AN3 + AN3P(8) * AT(13,J)
18    AN3 = AN3 + AN3P(9) * AT(14,J)
19    AN3 = AN3 + AN3P(10)* AT(15,J)
20    AN3 = AN3 + AN3P(11)* AT(16,J)
      AN3 = AN3/ZHE1

!      write (*,*) 'op_kiel_he1: ',fre,f, &
!       an31,an32,an33,an34,an35,an36,an37,an38,an39,an310,an311
!      write (*,*) 'op_kiel_he1: at ',   at(7:16,j)

1000  FORMAT(/,3X,'HE1-PHOTOIONISATIONSQUERSCHNITTE NACH', &
                 & 'JACOBS,PHYS.REV.,A9,1938(1974)'       )
1010  FORMAT(/,3X,'HE1-PHOTOIONISATIONSQUERSCHNITTE NACH GINGERICH')

  END SUBROUTINE OP_KIEL_HE1
