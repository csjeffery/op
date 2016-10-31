      SUBROUTINE OP_KIEL_Z1 (T,FREQ,ASI,AMG,AAL)

!****  << AVMET >> NAME CHANGED TO AVOID CLASH WITH QUB
!***** SILIZIUM, MAGNESIUM U. ALUMINIUM OPAZITAET PRO NEUTRALES ATOM
!****  PROGRAMM VON CARBON U. GINGERICH, PROC. 3. HARVARD-SMITH. CONF.
!*****        ON STELL. ATMOSPHERES, S. 397

      use op_precisions

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FREQ

      real, intent(out) :: ASI, AMG, AAL

      real :: FLIM, TK, A, B, C, D, E, G, V, W, X, Y, Z
      real(kind=qq) :: HV, WFR, CFR, FAK
      integer :: I

      DIMENSION FLIM(11)
      DATA  FLIM/ 19.63E14, 18.47E14, 17.81E14, 15.08E14, 14.43E14, &
       & 11.91E14,  7.99E14,  7.73E14,  6.00E14,  4.29E14,  3.5E14 /

      TK = (.43429/5040.)*T
      HV = 4.1352E-15*FREQ
      A=0.
      B=0.
      C=0.
      E=0.
      G=0.
      V=0.
      W=0.
      X=0.
      Y=0.
      Z=0.
! **** ALLOW FOR NU**-14 DEPENDENCE OF MG TRIPLET P, 2515 A
      WFR = (11.91E14/FREQ)**11
! ***** ALOW FOR NU**-1.5 DEPENDENC OF SI SINGLET S, 1986 A
      CFR = SQRT((FREQ/15.08E14)**3)
      DO 40  I = 1,10
        IF (FLIM(I) - FREQ)  50,50,40
 40   CONTINUE
      I = 11
      GOTO 11
 50   GOTO (1,2,3,4,44,5,6,7,8,10,11), I
! ***** GROUND STATE 3P2 TRIPLET P FROM SI, 1526 A
 1    IF (FREQ - 22.20E14)  101,101,102
! ***** RICH EXPERIMENTAL VALUE, 40 MGBARNS AT LIMIT, CONSTANT WITH FREQ
 101  A = (FREQ/19.63E14)**3
        GOTO 2
 102  A = (22.20E14/FREQ)**2*1.446
! ***** RICH GUESS OF NU**-5 FALL OFF BELOW 1350 A
! ***** GROUND STATE 3S2 SINGLET S FROM MG, 1620 A
! ***** DITCHBURN A. MARR EXPERIMENTAL VALUE, 1.18 MGBARNS AT LIMIT
 2    V = 0.0282
      I = 3
! ***** 3P2 SINGLET D FROM SI, 1682 A
! ***** RICH EXPERIM. VALUE , 35 MGBARNS AT LIMIT, CALCULATED NU**3
 3    B = 0.444/EXP(.78/TK)
! ***** 3P2 SINGLET S FROM SI, 1986 A
! ***** RICH QUANTUM DEFECT CALCUL., 46.5 MGBARNS AT LIMIT, NU**-1.5
 4    C = CFR*0.0628/EXP(1.91/TK)
! ***** GROUND STATE OF AL, 2076 A
! ***** PARKINSON AND REEVES EXPERIMENT. VALUE, 22 MGBARNS
 44   G = 0.230
! ***** 3P TRIPLET & ODD FROM MG, 2515 A
! ***** BOTTICHER EXPERIMENT. VALUE, 45 MGBARNS, NU**-14 DEP.
 5    W = WFR*2.32/EXP(2.71/TK)
! ***** 3P SINGLET P ODD FROM MG, 3750 A
 6    X = 0.065/EXP(4.33/TK)
! ***** CONTINUUM FROM 4S FOR SI, 3880 A
 7    D = .97*TK/EXP(5.0/TK)
! ***** CONTINUUM FROM 4S FOR MG, 5000 A
 8    Y = 1.90*TK/EXP(5.24/TK)
! **** 3D TERMS FROM SI
      E = .069/EXP(5.7/TK)
! ***** 3D TERMS FROM MG, 7000 A
 10   Z = .090/EXP(5.9/TK)
      GOTO (13,13,13,13,13,13,13,13,12,11,11), I
! ***** HIGHER TERM FOR MG
 11   Y = .962*TK/EXP((7.64-HV)/TK)
! ***** HIGHER TERMS FOR SI
 12   D = .493*TK/EXP((8.14-HV)/TK)
 13   FAK =          2.815E29/FREQ/FREQ/FREQ
      ASI = FAK*(A+B+C+D+E)
      AMG = FAK*(V+W+X+Y+Z)
      AAL = FAK*G

      RETURN

      END subroutine OP_KIEL_Z1
