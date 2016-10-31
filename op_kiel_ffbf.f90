SUBROUTINE OP_KIEL_FFBF (TT,PGJ,PEJ,WL,OPAC,CAP,SIG,J)

  !+
  !
  !       Calculate total monochromatic continuous opacity
  !       absorption coefficient per gram of stellar material
  !       for the elements H and He in various ionization stages
  !       alpha gives the relative abundance of atoms by number
  !       Uses Peach and Kiel opacities.
  !
  !
  !
  !   Input.
  !        TT     - double precision -  temperature 
  !        PGJ    - double precision - gas pressure 
  !        PEJ    - double precision - electron pressure 
  !        WLS    - double precision - wavelength
  !        J      - integer          - current depth zone pointer
  !
  !   Output.
  !        OPAC   - total opacity
  !        CAP    - absorption coefficient (KAP in opac, RBF in hydro)
  !        SIG    - scattering coefficient (SIG in opac, SIGM in hydro) 
  !
  !-

  ! LTE Precisions
  USE OP_PRECISIONS, ONLY : QQ

  ! Global constants.
  USE TAPM_PAR
  USE TAPP_PAR

  ! Global variables.
  USE OP_GLOBAL
  USE OP_PEACH

  ! No implicit typing.
  IMPLICIT NONE

  !  Local constants.
  DOUBLE PRECISION, PARAMETER :: HH = 4.7992155E-11
  DOUBLE PRECISION, PARAMETER :: AMP = TAPP__MP * 1.E3
  DOUBLE PRECISION, PARAMETER :: C = TAPP__C*100.          ! cm /s 
  DOUBLE PRECISION, PARAMETER :: SIGMA = 5.6705E-5         ! erg cm^-3 K^-4
  DOUBLE PRECISION, PARAMETER :: BK = 1.38054E-16
  DOUBLE PRECISION, PARAMETER :: TENLOG = 2.302585093      ! ln(10)

  !  Local constants.
  DOUBLE PRECISION, PARAMETER :: THOM1 = 8.0 *TAPM__PI * (TAPP__E*TAPP__C*10.)**4
  DOUBLE PRECISION, PARAMETER :: THOM2 = 3.0 * (TAPP__ME*1.E3)**2 
  DOUBLE PRECISION, PARAMETER :: THOM3 = (TAPP__C*100)**2 
  DOUBLE PRECISION, PARAMETER :: THOM = THOM1/(THOM2*THOM3**2) 

  ! Subroutine arguments (Given).
  REAL, INTENT(IN) :: TT, PGJ, PEJ 
  REAL, INTENT(IN) :: WL
  INTEGER, INTENT(IN) :: J

  ! Subroutine arguments (Returned).
  REAL, INTENT(OUT) :: OPAC
  REAL, INTENT(OUT) :: SIG
  REAL, INTENT(OUT) :: CAP

  ! Local variables
  REAL(KIND=QQ) :: FRE
  REAL(KIND=QQ) :: WLS
  REAL :: AMUMPI, STIM, SIGE
  REAL :: STUMPI

  !! outputs of opk2 subroutine:
  REAL :: ASI1, AMG1, AAL1, AH1, AHM
  REAL :: AHE1, AHE2, AHEM, AC1, AC2
  REAL :: AC3, AC2X, ACM, AN1, AN2, AN3 

  !! External function declarations
  REAL, EXTERNAL :: OP_KIEL_SI2
  REAL, EXTERNAL :: OP_KIEL_MG2
  REAL, EXTERNAL :: OP_KIEL_CA2
  REAL, EXTERNAL :: OP_KIEL_O1 

  REAL :: RAYH,RAYHE,RAYC,RAYN,RAYO

  REAL :: OPK_H1, OPK_HM, OPK_HE1, OPK_HE2, OPK_HEM
  REAL :: OPK_C1, OPK_C2, OPK_C3, OPK_C2X, OPK_CM
  REAL :: OPK_N1, OPK_N2, OPK_N3, OPK_MG1, OPK_MG2
  REAL :: OPK_SI1, OPK_SI2, OPK_AL1, OPK_CA2, OPK_O1

  !.

  ! write (*,*) tt,pgj,pej,wls,j

  !  Collect frequency
  WLS = WL
  FRE = 2.997927E18 / WLS

  !  Set summation to zero
  CAP = 0. 
  SIG = 0.

  ! mean atomic mass factor = 1. / (mu * mH)
  AMUMPI = 1.0 / (OP__ALPMU(J) * AMP)

  ! Stimulated emission
  STIM = (1. - EXP(-HH*FRE/TT))

  !  Thomson scattering
  !      SIGE = 0.6655E-24 * PEJ / (PGJ-PEJ)  !!typo? --> change:
  SIGE = THOM * PEJ / (PGJ-PEJ)

  !  Rayleigh scattering
  CALL OP_KIEL_RAY (WLS,J,RAYH,RAYHE,RAYC,RAYN,RAYO)

  ! write (*,*) wls, rayh, rayhe, rayc, rayn, rayo

  SIG =       RAYH *OP_XION(1,1,J)*OP_ALPNN(1,J)
  SIG = SIG + RAYHE*OP_XION(1,2,J)*OP_ALPNN(2,J)
  IF (OP_ALPNN(6,J) > 0.0)   SIG = SIG + RAYC *OP_XION(1,6,J)*OP_ALPNN(6,J)
  IF (OP_ALPNN(7,J) > 0.0)   SIG = SIG + RAYN *OP_XION(1,7,J)*OP_ALPNN(7,J)
  IF (OP_ALPNN(8,J) > 0.0)   SIG = SIG + RAYO *OP_XION(1,8,J)*OP_ALPNN(8,J)

  !  Nett contribution to opacity by scattering
  SIG = (SIG + SIGE) * AMUMPI

  ! write(*,*) 'sig', SIG, SIGE, AMUMPI

  STUMPI = STIM * AMUMPI

  !  Opacity due to neutral metals (Si, Mg and Al)
  CALL OP_KIEL_Z1   (TT,FRE,ASI1,AMG1,AAL1)
  OPK_MG1 = AMG1 * OP_XION(1,12,J) * OP_ALPNN(12,J) * STUMPI
  OPK_AL1 = AAL1 * OP_XION(1,13,J) * OP_ALPNN(13,J) * STUMPI
  OPK_SI1 = ASI1 * OP_XION(1,14,J) * OP_ALPNN(14,J) * STUMPI
  CAP = CAP + OPK_MG1 + OPK_AL1 + OPK_SI1

  !  write(*,*)  'neutral metals', CAP

  !  Opacity due to singly ionized metals (SiII, MgII, CaII and OI)
  OPK_SI2 = OP_KIEL_SI2(FRE,TT) * OP_XION(2,14,J) * OP_ALPNN(14,J) * STUMPI
  OPK_MG2 = OP_KIEL_MG2(FRE,TT) * OP_XION(2,12,J) * OP_ALPNN(12,J) * STUMPI
  OPK_CA2 = OP_KIEL_CA2(FRE,TT) * OP_XION(2,20,J) * OP_ALPNN(20,J) * STUMPI
  OPK_O1  = OP_KIEL_O1(FRE,TT) * OP_XION(1,8,J) * OP_ALPNN(8,J)  * STUMPI
  CAP  = CAP + OPK_SI2 + OPK_MG2 + OPK_CA2 + OPK_O1

  !  write(*,*) 'op_kiel_ffbf: metals', J, OPK_SI2, OPK_MG2, OPK_CA2, OPK_O1 
  !  write(*,*) 'ionized metals', CAP

  !  Opacity due to hydrogen (H-,H)
  IF (OP_ALPNN(1,J) .NE. 0.) THEN
     CALL OP_KIEL_H1 (J,TT,FRE,    OP_XUU(1,1,J), AH1)
     CALL OP_KIEL_HM (TT,FRE,PEJ,OP_XUU(1,1,J), AHM)
     OPK_H1 = AH1 * OP_XION(1,1,J) * OP_ALPNN(1,J) * STUMPI
     OPK_HM = AHM * OP_XION(1,1,J) * OP_ALPNN(1,J) * AMUMPI ! does AHM include STIM ?
     CAP = CAP + OPK_H1 + OPK_HM
  ENDIF

  !  write(*,*) 'op_kiel_ffbf: H', J, op_xuu(1,1,j)
  !  write(*,*) 'op_kiel_ffbf: H', J, AH1, AHM, PEJ
  !  write(*,*) 'op_kiel_ffbf: H', J, OPK_H1, OPK_HM

  !  Opacity due to helium (He-,He,HeII)
  IF (OP_ALPNN(2,J) .NE. 0.) THEN
     CALL OP_KIEL_HE1 (J,TT,FRE,    OP_XUU(1,2,J), AHE1)
     CALL OP_KIEL_HE2 (J,TT,FRE,    OP_XUU(2,2,J), AHE2)
     CALL OP_KIEL_HEM (TT,FRE,PEJ,OP_XUU(1,2,J), AHEM)
     OPK_HE1 = AHE1 * OP_XION(1,2,J) * OP_ALPNN(2,J) * STUMPI
     OPK_HE2 = AHE2 * OP_XION(2,2,J) * OP_ALPNN(2,J) * STUMPI
     OPK_HEM = AHEM * OP_XION(1,2,J) * OP_ALPNN(2,J) * STUMPI
     CAP = CAP + OPK_HE1 + OPK_HE2 + OPK_HEM
  ENDIF

  !  write(*,*) 'op_kiel_ffbf: He', J, op_xuu(1,2,j), op_xuu(2,2,j)
  !  write(*,*) 'op_kiel_ffbf: He', J, AHE1, AHE2, AHEM
  !  write(*,*) 'op_kiel_ffbf: He', J, OPK_HE1, OPK_HE2, OPK_HEM

  !  Opacity due to carbon (C-,C,CII,CIII)
  IF (OP_ALPNN(6,J) .NE. 0.) THEN
     CALL OP_KIEL_C1  (TT,FRE,    WLS,OP_XUU(1,6,J), AC1)
     CALL OP_KIEL_C2  (TT,FRE,    WLS,OP_XUU(2,6,J), AC2)
     CALL OP_KIEL_C3  (TT,FRE,    WLS,OP_XUU(3,6,J), AC3)
     CALL OP_KIEL_C2X (TT,FRE,    WLS,OP_XUU(2,6,J), AC2X)
     CALL OP_KIEL_CM  (TT,FRE,PEJ,WLS,OP_XUU(1,6,J), ACM)
     OPK_C1 = AC1 * OP_XION(1,6,J) * OP_ALPNN(6,J) * STUMPI
     OPK_C2 = AC2 * OP_XION(2,6,J) * OP_ALPNN(6,J) * STUMPI
     OPK_C3 = AC3 * OP_XION(3,6,J) * OP_ALPNN(6,J) * STUMPI
     OPK_C2X= AC2X* OP_XION(2,6,J) * OP_ALPNN(6,J) * STUMPI
     OPK_CM = ACM * OP_XION(1,6,J) * OP_ALPNN(6,J) * STUMPI
     CAP = CAP + OPK_C1 + OPK_C2 + OPK_C3 + OPK_C2X + OPK_CM
  ENDIF

  !  write (*,*) 'op_kiel_ffbf: C', J, AC1, AC2, AC3, AC2X, ACM
  !  write (*,*) 'op_kiel_ffbf: C', J, OPK_C1, OPK_C2, OPK_C3, OPK_C2X, OPK_CM

  !  Opacity due to nitrogen (N,NII,NIII)
  IF (OP_ALPNN(7,J) .NE. 0.) THEN
     CALL OP_KIEL_N1 (TT,FRE,WLS,OP_XUU(1,7,J), AN1)
     CALL OP_KIEL_N2 (TT,FRE,WLS,OP_XUU(2,7,J), AN2)
     CALL OP_KIEL_N3 (TT,FRE,WLS,OP_XUU(3,7,J), AN3)
     OPK_N1 = AN1 * OP_XION(1,7,J) * OP_ALPNN(7,J) * STUMPI
     OPK_N2 = AN2 * OP_XION(2,7,J) * OP_ALPNN(7,J) * STUMPI
     OPK_N3 = AN3 * OP_XION(3,7,J) * OP_ALPNN(7,J) * STUMPI
     CAP = CAP + OPK_N1 + OPK_N2 + OPK_N3
  ENDIF

  !  write (*,*) 'op_kiel_ffbf: N', J, AN1, AN2, AN3
  !  write (*,*) 'op_kiel_ffbf: N', J, OPK_N1, OPK_N2, OPK_N3

  ! write (*,*) 'op_kiel_ffbf: cap,sig ', cap, sig

  !  All opacities combined
  OPAC = CAP + SIG

END SUBROUTINE OP_KIEL_FFBF
