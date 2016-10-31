SUBROUTINE OP_HLIN (JD,IW, TT,PG,PE,WL,DOP, CAPLIN,CAPCON,CAPHOP, IF2)

  !+
  !
  !  Name:
  !     OP_HLIN
  !
  !  Purpose:
  !     To calculate the total monochromatic opacity or
  !     line absorption coefficient per gram for Hydrogen 
  !     for line-blanketed model atmosphere calculations.
  !
  !  Language:
  !     Fortran 90
  !
  !  Description:
  !
  !     The original version of this code used formulae given by Griem to give
  !     fast approximations to the H line opacity. CSJ failed to implement
  !     these sucessfully, so turned to another set of public domain routines
  !     due to Paul Barklem and Nik Piskunov (Uppsala), with input from 
  !     Deane Peterson, Bob Kurucz and Kjell Eriksson.
  ! 
  !     The routine now returns the b-b and b-f absorption coefficients due to 
  !     Lyman and Balmer series H bound levels at WAVE employing the occupation probability 
  !     formalism (OPF) of Daeppen, Anderson & Mihalas ApJ 1987 (DAM).
  !
  !     The extension of the occupation probability formalism to non-LTE by 
  !     Hubeny, Hummer and Lanz A&A 1994 (HHL) is available, but not used in this
  !     implementation. 
  !
  !     The routine can be used to compute only the H line opacity alone (cases 2,3)
  !     or the combined line and H contnuous opacity. The latter is correct if 
  !     using the OPF, as the depression of the continuum correctly results in a 
  !     changeover from line to continuum absorption at increasing density, along 
  !     with the wavelength blurring of the ionization edges. 
  !
  !     At present, the use of the continuous H opacity is slightly cludged. In "FLXOS", 
  !     SterneOS computes the H opacity from Opacity Project x-sections, with no 
  !     OPF. It also uses strict LTE occupation numbers for the H levels. In HBOP, 
  !     The occupation numbers are changed by the EPF, while I haven't examined 
  !     where the x-sections come from. Therefore, if the continuous H opacity 
  !     is computed here -- including EPF (CAPCON), the OP equivalent must also be 
  !     isolated (CAPHOP). 
  !
  !     Thus, FLXOS can then correct the total continuous opacity 
  !     CAC = CAC + CAPCON - CAPHOP. 
  !     This remains corect if both CAPCON and CAPHOP are zero. 
  !     
  !
  !  Input arguments:
  !     JD        Depth zone indicator
  !     IW        wavelength index (set to maxnw+1 or 0 for refwl)
  !     TT        Temperature (K)
  !     PG        Gas pressure (dyn / cm^2)
  !     PE        Electron pressure (dyn / cm^2)
  !     WL        Wavelength (A)
  !     DOP       Reduced Doppler width
  !
  !     IF2       switch choice of calculation
  !         1:    all done by HBOP -- Barklem code, including dissoultion of continuum
  !         2:    line opacity from HLINOP, level populations from OP 
  !         3:    line opacity from HLINOP, level populations as in Spectrum. 
  !         0:    computes all three (used for debugging) but returns (1). 
  !
  !  Output arguments:
  !     CAPLIN    Monochromatic b-b opacity  (cm^2 / g)
  !     CAPCON    Continuous (bf) opacity due to H (Barklem x-sections) (cm^2 / g)
  !     CAPHOP    Continuous (bf) opacity due to H (OP x-sections) (cm^2 / g)
  !
  !  Authors:
  !     UH:   Uli Heber (Kiel / Bamberg)
  !     NTB:  Natalie Behara (Armagh / Meudon)
  !     CSJ:  Simon Jeffery (Armagh Observatory)
  !
  !  History:
  !       -   -1982 (UH) : Original version based on Griem theory
  !       -   -2006 (NTB): Adapted for STERNE
  !       -May-2007 (CSJ): Barklem algorithm (HBOP) incorporated
  !     20-Jul-2007 (CSJ): Testing completed
  !
  !  Bugs:
  !     probably lots -- tell me about any you don't see here
  ! 
  !     1) see above for mismatch between OP and HBOP continuous opacity
  !     2) do better at avoiding repeat calcs of occ. numbers in HBOP and HLINOP
  !
  !-


  USE LTE_MODEL
  USE OP_PRECISIONS
  USE OP_GLOBAL
  USE TAPM_PAR
  USE TAPP_PAR

  IMPLICIT NONE

  !! input
  INTEGER :: JD            ! depth in the atmosphere
  INTEGER :: IW            ! wavelength pointer
  REAL :: TT               ! temperature
  REAL :: PG               ! total (gas) pressure
  REAL :: PE               ! electron pressure
  REAL(KIND=QQ) :: WL      ! wavelength
  REAL :: DOP              ! reduced dopler width delta_lambda / lambda_0

  !! output
  REAL(KIND=QQ) :: CAPLIN  ! Line opacity due to H
  REAL(KIND=QQ) :: CAPCON  ! Continuous opacity due to H, incl. lines
  REAL(KIND=QQ) :: CAPHOP  ! Continuous opacity from OP, excl. lines

  !! switch
  INTEGER :: IF2



  !! constants
  REAL, PARAMETER :: PI  = TAPM__PI
  REAL, PARAMETER :: AN0 = TAPP__E0 / 10.    ! .8852554E-12  ! = E**2*PI/(EM*C*C)
  REAL, PARAMETER :: SF  = 0.0265 ! pi*e*e/m/c [cgs]  sigma = sf * f * norm-prof
  REAL, PARAMETER :: AMP = TAPP__MP * 1.E03  ! 1.673E-24
  REAL, PARAMETER :: E   = TAPP__E * TAPP__C * 1.E01 ! e in e.s.u 
  REAL, PARAMETER :: EM  = TAPP__ME * 1.E03
  REAL, PARAMETER :: C   = TAPP__C * 1.E02
  REAL, PARAMETER :: RYD = TAPP__RYD / ( 1. + 1./TAPP__MPBME ) / 100. ! 109677.576
  REAL, PARAMETER :: BC  = TAPP__K * 1.E07 ! 1.38046E-16 
  REAL, PARAMETER :: H   = TAPP__H * 1.E07 
  REAL, PARAMETER :: HCK = TAPP__H * TAPP__C * 1.E+10 / TAPP__K

  INTEGER :: J
  INTEGER, PARAMETER :: NLINES = 100  !! effective series limit


  ! Atomic data for hydrogen lines
  ! Hydrogen Lyman series wavelengths in vacuum
  REAL*8, DIMENSION(NLINES), SAVE ::  WL0L = & 
       (/ 0.0, 1215.67, 1025.72, 972.537, 949.734, 937.803, &
       &  930.748, 926.226, 923.150, 920.963, 919.352, 918.129, &
       &  917.181, 916.429, 915.824, 915.329, 914.919, 914.576, &
       &  914.286, 914.039, 913.826, 913.641, 913.480, 913.339, &
       &  913.215, 913.104, 913.006, 912.918, 912.839, 912.768, &
       &  912.703, 912.645, 912.592, 912.543, 912.499, 912.458, &
       &  912.420, 912.385, 912.353, 912.324, (0.0,J=41,NLINES) /)
  ! Hydrogen Balmer series wavelengths in air
  REAL*8, DIMENSION(NLINES), SAVE ::  WL0B = & 
       (/ 0.0, 0.0, 6562.80, 4861.32, 4340.46, 4101.73, &
       &  3970.07, 3889.05, 3835.38, 3797.90, 3770.63, 3750.15, 3734.37, &
       &  3721.94, 3711.97, 3703.85, 3697.15, 3691.55, 3686.83, 3682.81, &
       &  3679.35, 3676.36, 3673.76, 3671.48, 3669.46, 3667.68, 3666.10, &
       &  3664.68, 3663.40, 3662.26, 3661.22, 3660.28, 3659.42, 3658.64, &
       &  3657.92, 3657.27, 3656.66, 3656.11, 3655.59, 3655.12,(0.0,J=41,NLINES)  /)
  ! Hydrogen Lyman series f - values
  REAL, DIMENSION(NLINES), SAVE ::  FL = & 
       (/ 0.0, 4.162E-1, 7.910E-2, 2.899E-2, 1.394E-2, &
       &  7.799E-3, 4.814E-3, 3.183E-3, 2.216E-3, 1.605E-3, 1.201E-3, &
       &  9.214E-4, 7.227E-4, 5.774E-4, 4.686E-4, 3.856E-4, 3.211E-4, &
       &  2.702E-4, 2.296E-4, 1.967E-4, 1.698E-4, 1.476E-4, 1.291E-4, &
       &  1.136E-4, 1.005E-4, 8.928E-5, 7.970E-5, 7.144E-5, 6.429E-5, &
       &  5.806E-5, 5.261E-5, 4.782E-5, 4.360E-5, 3.986E-5, 3.653E-5, &
       &  3.357E-5, 3.092E-5, 2.854E-5, 2.640E-5, 2.446E-5, (0.0,J=41,NLINES) /)
  ! Hydrogen Balmer series f - values
  REAL, DIMENSION(NLINES), SAVE ::  FB = & 
       (/ 0.0, 0.0, 6.407E-1, 1.193E-1, 4.467E-2, 2.209E-2, &
       & 1.270E-2, 8.036E-3, 5.429E-3, 3.851E-3, 2.835E-3, 2.151E-3, &
       & 1.672E-3, 1.326E-3, 1.070E-3, 8.764E-4, 7.270E-4, 6.099E-4, &
       & 5.167E-4, 4.416E-4, 3.805E-4, 3.302E-4, 2.884E-4, 2.534E-4, &
       & 2.238E-4, 1.987E-4, 1.772E-4, 1.587E-4, 1.427E-4, 1.288E-4, &
       & 1.167E-4, 1.060E-4, 9.658E-5, 8.825E-5, 8.086E-5, 7.427E-5, &
       & 6.837E-5, 6.309E-5, 5.834E-5, 5.405E-5, (0.0,J=41,NLINES) /)

  !   ! Hydrogen statistical weights
  REAL, DIMENSION(3) :: GH = (/ 2., 8., 18. /)
  !   ! Hydrogen energy levels (eV)
  REAL, DIMENSION(3) :: EH = (/ 0., 10.1985, 12.0885 /)


  ! External functions
  REAL, EXTERNAL :: HLINOP
  REAL, EXTERNAL :: HFNM


  !! local
  REAL :: FF  ! temporary holder for oscillator strength
  INTEGER :: A, B
  LOGICAL, SAVE :: LHDAT_INIT = .FALSE.

  REAL :: BKT        ! Boltzmann constant (erg/K) * temperature (K)
  REAL :: STIM       ! Correction for stimulated emission
  REAL :: AM         ! mean mass per particle (g)


  ! Variables for use with Barklem codes
  INTEGER, PARAMETER :: NHBOP=100
  REAL, DIMENSION(nhbop) :: NPOP   !  array of level populations to pass to Barklem code

  REAL*8 :: WAVE  !  interrogation wavelength in A
  REAL :: NI      !  ion number density in cm-3
  REAL :: NH1     !  number density of H I in cm-3
  REAL :: NHE1    !  number density of He I in cm-3
  REAL :: NE      !  electron number density in cm-3

  REAL :: hbop_tot  !  the total (line + continuous) absorption coefficient
  REAL :: hbop_con  !  the continuous absorption coefficient
  REAL :: hbop_lin, hbop_lin2, hbop_lin3 !  the line absorption coefficient

  ! Variables for computing occupation numbers and local line opacity
  REAL :: U_K, U_J, EQUIL, FEX, FION, HPHI

  ! H continuum opacities from OP
  REAL :: KAPPAX, KBF
  INTEGER :: II, IA, IS


  !.


  ! Initialise wavelength and f-values for high-order Ly and Ba lines

  ! but see fn HFNM in Barklem's codes

  IF ( .NOT. LHDAT_INIT ) THEN
     DO J = 41,NLINES
        ! vacuum
        WL0L(J) = 1.E08 / ( RYD * ( 1. - 1./(J*J) ) ) 
        FL(J)   = HFNM(1,J)
        ! air
        WL0B(J) = 1.E08 / ( RYD * ( 1./4. - 1./(J*J) ) ) / 1.000292
        FB(J)   = HFNM(2,J)
     END DO
     LHDAT_INIT = .TRUE.
  END IF


  !.

  ! Ensure returned values are set
  CAPLIN = 0.0
  CAPCON = 0.0
  CAPHOP = 0.0

  ! Lyman series ?
  IF  (WL > 800 .AND. WL < 1450) THEN
     A = 1
     ! Balmer series ?
  ELSEIF (WL > 3400 .AND. WL < 8000) THEN
     A = 2
     ! Forget it !
  ELSE 
     RETURN
  ENDIF

  ! H-abundance = 0 ?
  IF ( op_alpnn(1,JD) == 0.)  RETURN


  ! Basic quantities

  BKT = BC*TT
  NE  = PE/BKT
  NI  = (PG-PE)/BKT
  !  division by OP__ALPMU(JD) not understtod but consistnet with 
  !  correction in option (2). 
  NH1  = NI * OP_XION(1,1,JD) * OP_ALPNN(1,JD) ! / OP__ALPMU(JD)
  NHE1 = NI * OP_XION(1,2,JD) * OP_ALPNN(2,JD) ! / OP__ALPMU(JD)
  AM  = 1.0 / (OP__ALPMU(JD) * AMP)
  ! Stimulated emission correction
  STIM = ( 1.D0 - EXP( - HCK / (TT * WL) ) )   
  WAVE = WL


  ! Use Barklem's code including occupation probability formalism (HBOP)
  ! what we do here is get the dissolution of the continuum bit 
  ! from HBOP and stick it onto CAPHOP and send CAPHOP back -- rather than
  ! use the continuous opacity from HBOP. 
  IF ( IF2 == 1 .OR. IF2 == 0 ) THEN

     SELECT CASE (A)
     CASE (1)
        CALL HBOP(WAVE, NHBOP-1, (/(1,J=2,NHBOP)/), (/(J,J=2,NHBOP)/), WL0L(2:NHBOP), &
             NH1, NHE1, NE, TT, DOP, &
             NPOP, 0, HBOP_TOT, HBOP_CON) 
     CASE (2)
        CALL HBOP(WAVE, NHBOP-2, (/(2,J=3,NHBOP)/), (/(J,J=3,NHBOP)/), WL0B(3:NHBOP), &
             NH1, NHE1, NE, TT, DOP, &
             NPOP, 0, HBOP_TOT, HBOP_CON) 
     END SELECT

     HBOP_LIN = ( HBOP_TOT - HBOP_CON ) * AM / NI
     HBOP_CON = HBOP_CON * AM  /  NI

     CAPLIN = HBOP_LIN 
     CAPCON = HBOP_CON

     ! compute the OP contribution from H to the continuous opacity
     ! CAPHOP can then be subtracted from CAC and replaced by CAPCON.
     KBF = 0.D0
     IA=1
     DO II=OP_ALPIMIN(IA),OP_ALPIMAX(IA)
        KAPPAX = 0.D0
        DO IS=1,OP_NS(II,IA)
           KAPPAX = KAPPAX + OP_XPOP(IS,II,IA,JD)*OP_XS(IW,IS,II,IA)
        ENDDO
        KAPPABF(IA,II) = KAPPAX*STIM*AM
        KBF = KBF + KAPPABF(IA,II)
     ENDDO
     CAPHOP = KBF

  ENDIF



  ! Use HLINOP code alone
  IF ( IF2 == 2 .OR. IF2 == 0  ) THEN

     ! Notes
     ! 1) In principle, occupation numbers are already available in OP_XPOP
     ! but the total number for n=1 , n=2, ... must be computed by summming over the
     ! correct slp states. At present, OP cannot tell which these are automatically
     ! although it is possible to hardwire them.
     ! 2) These are strict LTE occupation numbers and therefore different from the 
     ! occupation numbers used in option 1), which use the OPF. 

     FION = NI * OP_ALPNN(1,JD) * OP_XION(1,1,JD) / OP_XUU(1,1,JD) ! / OP__ALPMU(JD)
     FION = SF * FION * GH(A) * EXP(-EH(A) * (TAPP__EV/TAPP__K) / TT )

     HBOP_LIN2 = 0.
     SELECT CASE (A)
     CASE (1)
        DO J = 2,NLINES
           HPHI = HLINOP(WAVE,A,J,WL0L(J),TT,NE,NH1,NHE1,DOP)
           HBOP_LIN2 = HBOP_LIN2 + FION * FL(J) * HPHI
        ENDDO
     CASE (2)
        DO J = 3,NLINES
           HPHI = HLINOP(WAVE,A,J,WL0B(J),TT,NE,NH1,NHE1,DOP)
           HBOP_LIN2 = HBOP_LIN2 + FION * FB(J) * HPHI
        ENDDO
     END SELECT
     HBOP_LIN2 = HBOP_LIN2 * STIM  * AM  /  NI
     CAPLIN = HBOP_LIN2 

  ENDIF



  ! Use HLINOP code alone with SPECTRUM derivation of FION.
  ! this works and shows that the populations in op_xpop are correct. 
  IF ( IF2 == 3 .OR. IF2 == 0  ) THEN

     !    Boltzmann equation
     U_K = OP_XUU(1,1,JD)
     FEX = GH(A) / U_K * EXP(-EH(A) * (TAPP__EV/TAPP__K) / TT )

     !    Ionization equilibria
     U_J = OP_XUU(2,1,JD)
     EQUIL = 6.6671737E-01 * (TT**2.5) / PE &
          *  U_J/U_K * EXP(-13.598 * (TAPP__EV/TAPP__K) / TT)

     !    FION is the CONST*(fractional population of the lower level)
     !                        *(stimulated emission correction)
     FION = 5.3313200E+03 * (TAPP__C*1.E+10) * OP_ALPNN(1,JD) / OP__ALPMU(JD) & 
            * FEX / (1.0+EQUIL)  

     ! The factor NI / AM is wrapped up above as (x * T / PE). I've taken
     ! it out here to show that the following treatment is the same as (2). 
     !     FION = FION * NI / AM
     HBOP_LIN3 = 0.
     SELECT CASE (A)
     CASE (1)
        DO J = A+1,NLINES
           HPHI = HLINOP(WAVE,A,J,WL0L(J),TT,NE,NH1,NHE1,DOP)
           HBOP_LIN3 = HBOP_LIN3 + FION * FL(J) * HPHI
        ENDDO
     CASE (2)
        DO J = A+1,NLINES
           HPHI = HLINOP(WAVE,A,J,WL0B(J),TT,NE,NH1,NHE1,DOP)
           HBOP_LIN3 = HBOP_LIN3 + FION * FB(J) * HPHI 
        ENDDO
     END SELECT
     HBOP_LIN3 = HBOP_LIN3 * STIM !  * AM / NI
     CAPLIN = HBOP_LIN3 

  ENDIF


  !  IF ( jd==25 .AND. IF2 == 0) THEN
  !     write(*,'(i6,i3,1p,12e11.3)') &
  !          int(wl),jd,hbop_lin,hbop_lin2,hbop_lin3,caplin,capcon,caphop
  !  ENDIF


  IF (IF2 == 0) CAPLIN = HBOP_LIN


END SUBROUTINE OP_HLIN
