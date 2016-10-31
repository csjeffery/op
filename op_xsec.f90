SUBROUTINE OP_XSEC (NOP,OPWL,OPXS,LIMIT,LL,WL,XSL)

  !+
  !
  !  Name:
  !     OP_XSEC
  !
  !  Purpose:
  !     Resamples OP cross-sections from the given OP wavelength grid
  !     to a new set of GRID wavelengths 
  !
  !  Input:
  !     NOP  - number of OP crossection points 
  !     OPWL - wavelength grid of OP crossection
  !     OPXS - OP crosssection values
  !
  !  Output:
  !     XSL  - cross-sections on new GRID
  !
  !  Authors:
  !     PMH: P.M.Harrison (St.Andrews)
  !     CSJ: C.S.Jeffery (Armagh)
  !
  !  History:
  !     16-Jul-97: Created (PMH)
  !     02-Mar-99: Tidied up (CSJ)
  !     13-Nov-03: modified for sterne90 (NTB)
  !-

  USE OP_PRECISIONS, ONLY : QQ
  USE OP_GLOBAL, ONLY : OP_MAXNX

  IMPLICIT NONE

  !  Subroutine Arguments (Given).
  INTEGER, INTENT(IN) :: NOP         ! number of OP crossection points
  REAL, INTENT(IN) :: OPWL(OP_MAXNX)     ! wavelength grid of OP crossection
  REAL, INTENT(IN) :: OPXS(OP_MAXNX)     ! OP crosssection values
  REAL, INTENT(IN) :: LIMIT
  INTEGER, INTENT(IN) :: LL          ! number of valid wavelengths in WL
  REAL, INTENT(IN) :: WL(*)          ! Wavelength grid for required x-sections

  !  Subroutine Arguments (Returned).
  REAL, INTENT(OUT) :: XSL(*)         ! average cross section on new grid

  !  Local variables.
  REAL ::  WLS, WLF, XWLS(OP_MAXNX), XWLF(OP_MAXNX)
  REAL :: YY 
  REAL :: SUM
  INTEGER :: IX, IW, I, IXS, IXF
  REAL :: DXDW

  REAL :: SM3
  REAL :: ENER_I, ENER3, ENEREXP, ENER0, ENER0_3
  REAL :: FREQ_I, FREQ3

  REAL(KIND=QQ), PARAMETER :: RINV = 911.2670747D0 ! 1/R (A)
  REAL(KIND=QQ), PARAMETER :: REV = 13.60559024D0  ! 1 Ryd (eV)

  !.

  !  Set up histogram wavelength limits of OP cross-sections.
  XWLS(1)       = OPWL(1)
  XWLF(1)       = 0.5D0*(OPWL(1)+OPWL(2))
  XWLS(2:NOP-1) = XWLF(1:NOP-2)
  XWLF(2:NOP-1) = 0.5D0*(OPWL(2:NOP-1)+OPWL(3:NOP))
  XWLS(NOP)     = XWLF(NOP-1)
  XWLF(NOP)     = OPWL(NOP)


  !  Cycle through new wavelength grid establishing average cross-section
  !  in each interval

  IX=1              ! OP cross-section grid counter
  DO IW=1,LL       ! cycle through new wavelength grid

     !  Set up histogram wavelength limits for Sterne.
     IF (IW==1) THEN 
        WLS=WL(IW)
        WLF=0.5D0*(WL(IW)+WL(IW+1))
     ELSEIF (IW < LL) THEN 
        WLS=0.5D0*(WL(IW)+WL(IW-1))
        WLF=0.5D0*(WL(IW)+WL(IW+1))
     ELSE
        WLS=0.5D0*(WL(IW)+WL(IW-1))
        WLF=WL(IW)
     ENDIF

     !  Assign values for Sterne's cross-sections.

     !  - For Sterne wavelengths less than OP wavelengths
     !    extrapolate towards zero cross-section (not working ...)
     IF (WLF < XWLS(1)) THEN
        XSL(IW) = 0.0

        !  - For Sterne wavelengths greater than OP wavelengths.
     ELSEIF (WLS > XWLF(NOP)) THEN
        XSL(IW)=0.D0

        !  - For GRID wavelengths within OP wavelength range:
     ELSE

        !write(*,*) 'XWLF', XWLF(IX), WLS, IX, NOP, WLF, XWLS(IX), WL(IW)
        DO WHILE((XWLF(IX) < WLS).AND.(IX < NOP))
           IX=IX+1
        ENDDO
        IXS=MAX(MIN(IX,NOP-1),2)
        DO WHILE((XWLF(IX) <= WLF).AND.(IX < NOP))
           IX=IX+1
        ENDDO
        IXF=MAX(MIN(IX,NOP),IXS+1)


        !    If GRID bin < OP bin => linear interpolation in wavelength 
        !
        !    Schematic of overlap of OP and GRID wavelengths
        !
        !  	   			   WL(IW)
        !  	   		   WLS  	       WLF
        !  	   	    x	    |	      x 	|	 x
        !
        !
        !      o            |             o                 |              o
        !                   
        !  XWLS          (IXS=IXF)                 
        !  XWLF                                         (IXS=IXF) 
        !  OPWL                        IXS=IXF

        IF (IXS >= IXF-1) THEN
           IF (WL(IW) < OPWL(IXS)) THEN
              XSL(IW) = ( (WL(IW)-OPWL(IXS-1)) * OPXS(IXS) + &
                   (OPWL(IXS)-WL(IW)) * OPXS(IXS-1) ) &
                   / (OPWL(IXS)-OPWL(IXS-1))
           ELSE
              IXS = MAX(2,IXS)
              XSL(IW) = ( (WL(IW)-OPWL(IXS)) * OPXS(IXS+1) + &
                   (OPWL(IXS+1)-WL(IW)) * OPXS(IXS) ) &
                   / (OPWL(IXS+1)-OPWL(IXS))
           ENDIF


           !    If GRID bin ~ OP bin => polynomial interpolation in wavelength 
           !    not used right now....
           !
           !	  ELSEIF (IXS.GE.IXF-3) THEN
           !            CALL NR_POLINT (OPWL(IXS),OPXS(IXS),3,WL(IW),XSL(IW),DXDW)


           !    If GRID bin > OP bin => simple average in wavelength
           !
           !    Schematic of overlap of OP and GRID wavelengths
           !
           !                                WL(IW)
           !                        WLS                 WLF
           !                 x       |         x         |        x
           !
           !
           !              o   |   o   |   o   |   o   |   o   |   o   |   o
           !                   
           !  XWLS          (IXS)                   (IXF)
           !  XWLF                  (IXS)                   (IXF) 
           !  OPWL               IXS                     IXF

        ELSE
           SUM=0.D0
           SUM=SUM+(XWLF(IXS)-WLS)*OPXS(IXS)
           SUM=SUM+(WLF-XWLS(IXF))*OPXS(IXF)
           IF (IXS < IXF-1) THEN
              DO I=IXS+1,IXF-1
                 SUM=SUM+(XWLF(I)-XWLS(I))*OPXS(I)
              ENDDO
           ENDIF
           XSL(IW)=SUM/(WLF-WLS)
        ENDIF

        !     Prevent extrapolation to negative opacity
        XSL(IW) = MAX(XSL(IW),0.D0)

        !! check that wavelengths do not fall under ionisation limit:

        IF (WL(IW) > LIMIT) XSL(IW) = 0.0 

     ENDIF

!     write(*,*)'op_xsec: ',iw,WL(iw),WLs,WLf
!     write(*,*)ixs,ixf,xWLs(ixs:ixf),xWLf(ixs:ixf)
!     write(*,*)opxs(ixs:ixf),xsl(iw)
!     write(*,*)

  ENDDO

END SUBROUTINE OP_XSEC
