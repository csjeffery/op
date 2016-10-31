SUBROUTINE OP_OP_FFBF (JD,I_OPWL,TL,PGL,PEL,WL,OPAC,KAP,SIG)

  !+
  !
  !  Name:
  !     OP_OP_FFBF
  !
  !  Purpose:
  !     To calculate the total monochromatic opacity due to 
  !     continuum (b-f and f-f transitions) and electron 
  !     scattering sources.
  !
  !  Language:
  !     Fortran 90
  !
  !  Input arguments:
  !     JD        Depth zone indicator
  !     I_OPWL    wavelength index (set to maxnw+1 or 0 for refwl)
  !     TL        Temperature (K)
  !     PGL       Gas pressure (dyn / cm^2)
  !     PEL       Electron pressure (dyn / cm^2)
  !     WL        Wavelength (A)
  !
  !  Output arguments:
  !     OPAC      Total monochromatic opacity...
  !     KAP       ...continuous (bf & ff).
  !     SIG       ...electron scattering.
  !
  !  Authors:
  !     PMH: P.M.Harrison (St.Andrews)
  !
  !  History:
  !     08-Dec-1994 (PMH): Modified for Sterne 3 from Sterne 2.1 routine OPKC.
  !     25-Aug-1997 (PMH): Bound-Free opacities for all species replaced 
  !                        using Opacity Project data
  !     12-Feb-1997 (PMH): Modified from Sterne 3 to Sterne 2.3 
  !        Feb-2004 (NTB): Modified and updated for SterneOS
  !     19-Jul-2007 (CSJ): Cleaned up 
  !
  !-

  !  LTE Precisions
  USE OP_PRECISIONS, ONLY : QQ

  !  TAP Global constants.
  USE TAPM_PAR
  USE TAPP_PAR

  !  Global OP library variables
  USE OP_GLOBAL

  IMPLICIT NONE

  !  Local constants.
  DOUBLE PRECISION, PARAMETER :: THOM1 = 8.0 *TAPM__PI * (TAPP__E*TAPP__C*10.)**4
  DOUBLE PRECISION, PARAMETER :: THOM2 = 3.0 * (TAPP__ME*1.E3)**2 
  DOUBLE PRECISION, PARAMETER :: THOM3 = (TAPP__C*100)**2 
  DOUBLE PRECISION, PARAMETER :: THOM = THOM1/(THOM2*THOM3**2) 
  DOUBLE PRECISION, PARAMETER :: HDK = TAPP__H / TAPP__K 
  DOUBLE PRECISION, PARAMETER :: AMP = TAPP__MP * 1.E3

  !! input:
  INTEGER, INTENT(IN) :: JD       ! Depth zone indicator
  INTEGER, INTENT(IN) :: I_OPWL   ! OP wavelength index
  REAL, INTENT(IN) :: TL          ! Temperature
  REAL, INTENT(IN) :: PGL         ! Gas pressure
  REAL, INTENT(IN) :: PEL         ! Electron pressure
  REAL, INTENT(IN) :: WL          ! Wavelength (A)

  !! output:
  REAL, INTENT(OUT) :: OPAC        ! Total cont. + scat. opacity     
  REAL, INTENT(OUT) :: KAP         ! Total Continuous opacity (bf & ff)
  REAL, INTENT(OUT) :: SIG         ! Total scattering opacity



  !  Local variables.
  INTEGER :: I

  !  Macroscopic.
  INTEGER :: IWL  = 0             ! wavelength index (local)
  INTEGER :: IA,II,IS             ! Atom, Ion and State indices

  REAL(KIND=QQ) :: FRE
  REAL :: STIM

  !  Atomic.
  REAL :: AM

  !  Scattering.
  REAL :: RAYH,RAYHE,RAYC,RAYN,RAYO ! Rayleigh
  REAL :: SIGE                      ! Thomson 

  !  Opacities.
  REAL :: KFF         ! Continuous opacity (free-free)
  REAL :: KBF         ! Continuous opacity (bound-free)

  !  Continuous Opacities.
  DOUBLE PRECISION :: KAPPAX                ! b-f op. (individual)
  REAL :: AH1,AHM,AHE1, AHE2, AHEM
  REAL :: AC1, AC2, AC2X, AC3
  REAL :: ACM, AN2, AN3
  REAL :: AN1
  REAL :: AMG1, AAL1, ASI1
  REAL :: CH, CHE

  REAL :: OP_MG1, OP_AL1, OP_SI1, OP_SI2, OP_MG2, OP_CA2, OP_O1


  ! Holders for intperolation in bf opacity x-sections
  REAL :: WTA, WTB
  REAL :: OP_XSA, OP_XSB

  !.

  ! write(*,*) 'op_op_ffbf: jd,iw, wl,t,pg,pe: ', jd, i_opwl, wl,tl,pgl,pel
  ! write(*,*) 'op_op_ffbf: xion : ', op_xion(1,1:8,jd)
  ! write(*,*) 'op_op_ffbf: alpha: ', op_alpnn(1:8,jd)

  !  Set reference wavelength pointer
  IWL = I_OPWL
  ! -- if I_OPWL == 0,  use default x-sections for ref wavelength (see op_inxs)
  IF (I_OPWL == 0) IWL = OP_MAXNW + 1 
  ! -- if I_OPWL < 0 , interpolate xsections to WL. (see below)

  !  Local values for wavlength and frequency
  FRE = TAPP__C * 1.E10 / WL

  !! mean atomic mass factor:
  AM  = 1.0 / AMP / OP__ALPMU(JD)

  ! Set summation to zero
  OPAC= 0.D0
  KAP = 0.D0 

  ! Stimulated emission correction
  STIM = (1.D0 - EXP(-HDK*FRE/TL))


  ! Electron scattering
  SIG = 0.D0

  !   Thomson
  SIGE = THOM * PEL / (PGL-PEL)
  ! write(*,*) 'op_op_ffbf: am,thom,sige: ', am,thom,sige

  !   Rayleigh
  CALL OP_KIEL_RAY (WL,JD,RAYH,RAYHE,RAYC,RAYN,RAYO)
  SIG =       RAYH *OP_XION(1,1,JD)*OP_ALPNN(1,JD)
  SIG = SIG + RAYHE*OP_XION(1,2,JD)*OP_ALPNN(2,JD)
  SIG = SIG + RAYC *OP_XION(1,6,JD)*OP_ALPNN(6,JD)
  SIG = SIG + RAYN *OP_XION(1,7,JD)*OP_ALPNN(7,JD)
  SIG = SIG + RAYO *OP_XION(1,8,JD)*OP_ALPNN(8,JD)
  ! write(*,*) 'op_op_ffbf: wl,sig,sige,rayh,rayhe: ', jd, wl,sig,rayh,rayhe,rayc,rayn,rayo

  !   Total
  SIG = (SIG + SIGE) * AM
  ! write(*,*) 'op_op_ffbf: wl,sig: ', jd, wl,sig

  ! Free-Free opacity 
  KFF = 0.D0

  ! H   
  IF (OP_ALPNN(1,JD) > 0.0)   THEN 

     !   H I
     CALL OP_FFH1 (JD,TL,FRE,OP_XUU(1,1,JD),AH1)
     OPFFH1 = AH1*STIM*OP_XION(1,1,JD)*OP_ALPNN(1,JD)*AM
     !   H minus
     CALL OP_FFHM (TL,FRE,PEL,OP_XUU(1,1,JD),AHM)
     OPFFHM = AHM*OP_XION(1,1,JD)*OP_ALPNN(1,JD)*AM

     KFF=KFF+OPFFH1+OPFFHM
  ENDIF
  ! write(*,*) 'op_op_ffbf: wl,ah1,ahm,pe: ', jd, wl,ah1,ahm,pel
  ! write(*,*) 'op_op_ffbf: wl,kff,kff_h: ', jd, wl,kff,opffh1,opffhm

  !  He
  IF (OP_ALPNN(2,JD) > 0.0) THEN 

     !   He I
     CALL OP_FFHE1 (JD,TL,FRE,OP_XUU(1,2,JD),AHE1)
     OPFFHE1 = AHE1*OP_XION(1,2,JD)*OP_ALPNN(2,JD)*AM*STIM
     !   He II
     CALL OP_FFHE2 (JD,TL,FRE,OP_XUU(2,2,JD),AHE2)
     OPFFHE2 = AHE2*OP_XION(2,2,JD)*OP_ALPNN(2,JD)*AM*STIM
     !   He minus
     CALL OP_FFHEM (TL,FRE,PEL,OP_XUU(1,2,JD),AHEM)
     OPFFHEM = AHEM*OP_XION(1,2,JD)*OP_ALPNN(2,JD)*AM*STIM

     KFF = KFF+OPFFHE1+OPFFHE2+OPFFHEM 
  ENDIF
  ! write(*,*) 'op_op_ffbf: wl,ahe1,ahe2,ahem: ', jd, wl,ahe1,ahe2,ahem
  ! write(*,*) 'op_op_ffbf: wl,kff,kff_he: ', jd, wl,kff,opffhe1,opffhe2,opffhem

  ! Carbon minus (free-free and bound-free from: 132, 457 MNRAS 1966)
  IF (OP_ALPNN(6,JD) > 0.0)   THEN 

     !  C minus
     CALL OP_FFCM (TL,FRE,PEL,WL,OP_XUU(1,6,JD),ACM)
     OPFFCM = ACM*OP_XION(1,6,JD)*OP_ALPNN(6,JD)*AM*STIM

     KFF = KFF+OPFFCM
  ENDIF
  ! write(*,*) 'op_op_ffbf: wl,kff,kff_cm: ', jd, wl,kff,opffcm


  ! Bound-Free opacity for all species being treated
  KBF = 0.D0

  IF (I_OPWL < 0) THEN 
     CALL GRID_LOC ( OP_WL, OP_NW, WL, IWL )
     IWL = MIN ( MAX (IWL,1), OP_NW-1 )
     WTA = (OP_WL(IWL+1) - WL) / (OP_WL(IWL+1)-OP_WL(IWL))
     WTB = (WL - OP_WL(IWL))   / (OP_WL(IWL+1)-OP_WL(IWL))
  ENDIF
  ! write (*,*) 'op_op_ffbf - grid_loc: ',wl, iwl, op_wl(iwl:iwl+1), wta, wtb
  DO OP_IAS = 1,OP_NATOMS
     IA = OP_SATOMS(OP_IAS)
     IF (OP_ALPNN(IA,JD) > 0.0) THEN
        IF (OP_ALPIMIN(IA) /= 0) THEN
           DO II=OP_ALPIMIN(IA),OP_ALPIMAX(IA)
              KAPPAX = 0.D0
              DO IS=1,OP_NS(II,IA)
                 IF (I_OPWL >= 0) THEN 
                 KAPPAX = KAPPAX + OP_XPOP(IS,II,IA,JD) * OP_XS(IWL,IS,II,IA)
                 ! write (*,*) is,  OP_XPOP(IS,II,IA,JD), OP_XS(IWL,IS,II,IA), kappax
                 ELSE
                 ! write (*,*) op_nw, i_opwl, iwl, op_maxnw, wl
                 OP_XSA = OP_XS(IWL,IS,II,IA); OP_XSB = OP_XS(IWL+1,IS,II,IA)
                 KAPPAX = KAPPAX + OP_XPOP(IS,II,IA,JD) * (WTA*OP_XSA + WTB*OP_XSB) 
                 ! write(*,*) is,  OP_XPOP(IS,II,IA,JD), op_xsa, op_xsb, kappax
                 ENDIF
                 ! write(*,*) op_xs(:,is,ii,ia)
              ENDDO
              KAPPABF(IA,II) = KAPPAX*STIM*AM
              ! write(*,*) 'op_op_ffbf: jd,iwl,wl,op_wl,iz,ii,kbf_ij ', jd,iwl, wl,op_wl(iwl),op_ias,ii,kappabf(ia,ii)
              KBF = KBF + KAPPABF(IA,II)
           ENDDO
        ENDIF
     ENDIF
  ENDDO


  ! All opacities combined
  KAP  =       KFF + KBF
  OPAC = SIG + KFF + KBF

  ! write(*,*) 'op_op_ffbf: sig,kff,kbf: ', jd, sig,kff,kbf

END SUBROUTINE OP_OP_FFBF
