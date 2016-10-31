SUBROUTINE OP_BFIN

  !! Reads in Opacity Project data for continuous opacities.
  !! Subroutine is a modified version of op_bfopin.f90.
  !!  History:
  !!     21-Mar-97: Created (PMH)
  !!     12-Feb-98: modified for Sterne 2 (PMH)
  !!     02-Mar-99: cleaned up (CSJ)
  !!     05-Nov-03: modified (NTB)

  USE LTE_WLGRID
  USE OP_PRECISIONS
  USE OP_GLOBAL

  IMPLICIT NONE 

  INTEGER :: I, J, K, IX

  !  Local variables.
  INTEGER :: FLAG, count

  INTEGER :: NSC, ISC, ISCX              !  no of levels in ion (+counters)
  INTEGER :: NZ                          !  no of protons in atom
  INTEGER :: NE                          !  no of electron in target ion
  INTEGER :: ISLP                        !  state identifier (slp)
  INTEGER :: ILV                         !  state identifier (n)
  REAL :: XCHIC              !  energy of level

  INTEGER :: NXC                         !  no of wavelengths in level (+counters)
  REAL :: OPXWL(OP_MAXNX)      !  OP x-s wavelengths (A)
  REAL :: OPXXS(OP_MAXNX)      !  OP x-s (cm^2)

  INTEGER :: DUM1, DUM2
  REAL :: DUMWT, DUMIE
  REAL :: DUMXWL, DUMXS, DUMXS2(OP_MAXNW)
  REAL :: STORE, LIMIT(0:OP_MAXNX), A_limit(0:OP_MAXNX), A_limit1(0:OP_MAXNX)
  REAL :: XCHIG             !  energy of ground state wrt continuum


  INTEGER :: ALPK                      !  atomic number
  INTEGER :: ALPI                      !  electron number on ion

  !  Holders for opacity data file names.
  CHARACTER(LEN=2) :: ION
  CHARACTER(LEN=2) :: ATOM
  CHARACTER(LEN=33) :: PFILE                  !  file holding energy level data
  CHARACTER(LEN=33) :: XFILE                  !  file holding x-section data

  ! unit numbers for files
  INTEGER :: pf=12, xf=13, out=10, xchi=14

!!$  LOGICAL  LXSOUT                     !  flag whether to printout x-s
!!$  DATA LXSOUT / .FALSE. /
  CHARACTER(LEN=4) :: INSLP                  !  ILV + ISLP
  CHARACTER(LEN=37) :: XSOUT                  !  file for transformed x-sections

  !  Local constants.
  REAL(KIND=QQ), PARAMETER :: RINV = 911.2670747D0 ! 1/R (A)
  REAL(KIND=QQ), PARAMETER :: REV = 13.60559024D0  ! 1 Ryd (eV)

  !! locator for OP data
  CHARACTER(LEN=7), PARAMETER :: PATH='OPDATA/'       

  INTEGER :: TOP = 14, CARB3=15, CARB1=16
  INTEGER :: IA, II, IS

  !           IONLIMIT controls lower energy limit used for cross-sections
  !           1 = use ionisation limit, 0 = use cross-section data limit
  REAL :: IONLIMIT=1.0


  ! Store the wavelength array for the OP xsections
  OP_NW = LL
  ALLOCATE(OP_WL(OP_NW))
  ALLOCATE(OP_NU(OP_NW))
  OP_WL(1:LL) = WL(1:LL)

  !write (*,*) 'OP_BFIN'
  !write (*,*) 'op_bfin: ',op_nw
  !write (*,*) 'op_bfin: ',op_wl(1:10)
  !write (*,*) 'op_bfin: ',op_wl(ll-9:ll)

  OP_XS = 0.0

  ! Set up array of selected atoms -- to be read in from TopBase
  ! Array used for indirect referencing of OP arrays. 
  ! currenty only uses elements in TopBase.
  OP_NATOMS = 18
  OP_SATOMS(1:OP_NATOMS) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14, 16, 18, 20, 26 /)

  !  Loop over atoms and ions
  DO OP_IAS = 1,OP_NATOMS
     K = OP_SATOMS(OP_IAS)
     ALPK=OP_ALPZ(k)                       ! atomic number of current atom


     IF (OP_ALPIMIN(k) /= 0) THEN
        DO j = OP_ALPIMIN(k),OP_ALPIMAX(k)    ! Loop over ions
           ALPI = ALPK-j+1              ! ion no. (j) -> elec. no.
           !       Open files to be read for current atom/ion:
           ATOM = '00'
           IF (ALPK .LT. 10) THEN
              WRITE (ATOM,'(A1,I1)') '0',ALPK
           ELSE
              WRITE (ATOM,'(I2)') ALPK
           END IF
           ION = '00'
           IF (ALPI .LT. 10) THEN
              WRITE (ION,'(A1,I1)') '0',ALPI
           ELSE
              WRITE (ION,'(I2)') ALPI
           END IF

           PFILE = PATH // 'p' // ATOM // ION // '.op'
           XFILE = PATH // 'x' // ATOM // ION // '.op'

           OPEN (unit=pf,file=PFILE,status='UNKNOWN')
           OPEN (unit=xf,file=XFILE,status='UNKNOWN')

           !            write (out,*) PFILE, XFILE

           !       Read in weights and ionisation energies:
           FLAG=1
           NSC=0
           LIMIT(0)=0.D0
           A_limit = 0.0 

           ! loop over states
           DO count=1,OP_MAXNS+1
              IF (FLAG == 1) THEN  
                 READ (pf,*) DUM1,DUM2,DUMWT,DUMIE
                 !! reads in slp, ilv, statistical weight and ionization potential(ryd)
                 IF ((NSC == OP_MAXNS) .OR. (DUM1 == 0)) THEN  ! end of file
                    FLAG=0
                    OP_NS(j,k)=NSC
                    !                  write(*,*) 'number of states', k, j, NSC
                 ELSE
                    NSC=NSC+1
                    OP_XWT(NSC,j,k)=DUMWT                 ! stat. weight
                    LIMIT(NSC)=-DUMIE*IONLIMIT            ! store ionisation limit
                    A_limit(NSC) = RINV/LIMIT(NSC)        ! Rydberg -> angstrom
                    A_limit1(NSC) = A_limit(NSC)  

                    ! correct OP limits...should be done in input file...
                    IF (k == 2 .AND. j == 1 .AND. NSC == 1) THEN
                       A_limit(NSC) = 504.26
                    END IF
                    IF (k == 2 .AND. j == 1 .AND. NSC == 3) THEN
                       A_limit(NSC) = 3124.15
                    END IF
                    IF (k == 6 .AND. j == 1 .AND. NSC == 1) THEN
                       A_limit(NSC) = 1101.09
                    END IF
                    IF (k == 6 .AND. j == 1 .AND. NSC == 2) THEN
                       A_limit(NSC) = 1239.0
                    END IF
                    IF (k == 6 .AND. j == 1 .AND. NSC == 3) THEN
                       A_limit(NSC) = 1444.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 1) THEN
                       A_limit(NSC) = 508.482
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 2) THEN
                       A_limit(NSC) = 655.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 4) THEN
                       A_limit(NSC) = 1002.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 5) THEN
                       A_limit(NSC) = 1218.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 6) THEN
                       A_limit(NSC) = 1248.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 7) THEN
                       A_limit(NSC) = 1540.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 8) THEN
                       A_limit(NSC) = 1877.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 9) THEN
                       A_limit(NSC) = 1956.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 11) THEN
                       A_limit(NSC) = 2536.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 12) THEN
                       A_limit(NSC) = 2929.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 14) THEN
                       A_limit(NSC) = 3504.0
                    END IF
                    IF (k == 6 .AND. j == 2 .AND. NSC == 15) THEN
                       A_limit(NSC) = 3612.0
                    END IF

                    IF (k == 6 .AND. j == 3 .AND. NSC == 1) THEN
                       A_limit(NSC) = 258.91
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 2) THEN
                       A_limit(NSC) = 299.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 3) THEN
                       A_limit(NSC) = 353.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 7) THEN
                       A_limit(NSC) = 675.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 8) THEN
                       A_limit(NSC) = 719.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 9) THEN
                       A_limit(NSC) = 785.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 10) THEN
                       A_limit(NSC) = 790.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 11) THEN
                       A_limit(NSC) = 860.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 12) THEN
                       A_limit(NSC) = 911.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 14) THEN
                       A_limit(NSC) = 1302.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 16) THEN
                       A_limit(NSC) = 1342.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 17) THEN
                       A_limit(NSC) = 1461.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 19) THEN
                       A_limit(NSC) = 1543.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 20) THEN
                       A_limit(NSC) = 1557.0
                    END IF
                    IF (k == 6 .AND. j == 3 .AND. NSC == 24) THEN
                       A_limit(NSC) = 1612.0
                    END IF


                    IF (k == 7 .AND. j == 1 .AND. NSC == 1) THEN
                       A_limit(NSC) = 853.05
                    END IF
                    IF (k == 7 .AND. j == 1 .AND. NSC == 2) THEN
                       A_limit(NSC) = 1020.438
                    END IF
                    IF (k == 7 .AND. j == 1 .AND. NSC == 3) THEN
                       A_limit(NSC) = 1131.234
                    END IF

                    IF (k == 7 .AND. j == 2 .AND. NSC == 1) THEN
                       A_limit(NSC) = 418.84
                    END IF
                    IF (k == 7 .AND. j == 2 .AND. NSC == 2) THEN
                       A_limit(NSC) = 447.0
                    END IF
                    IF (k == 7 .AND. j == 2 .AND. NSC == 3) THEN
                       A_limit(NSC) = 485.0
                    END IF
                    IF (k == 7 .AND. j == 2 .AND. NSC == 8) THEN
                       A_limit(NSC) = 1112.0
                    END IF
                    IF (k == 7 .AND. j == 2 .AND. NSC == 9) THEN
                       A_limit(NSC) = 1115.0
                    END IF
                    IF (k == 7 .AND. j == 2 .AND. NSC == 11) THEN
                       A_limit(NSC) = 1346.0
                    END IF
                    IF (k == 7 .AND. j == 2 .AND. NSC == 12) THEN
                       A_limit(NSC) = 1384.0
                    END IF


                    OP_XCHI(NSC,j,k)= - DUMIE * REV       ! Rydberg -> eV
                    IF (NSC == 1) XCHIG=OP_XCHI(NSC,j,k)  ! energy of ground state
                    ! relative to continuum
                    OP_XCHI(NSC,j,k)=XCHIG-OP_XCHI(NSC,j,k)  ! energy relative 
                    ! to ground state
                 END IF
              END IF
           END DO

           !         Read in wavelengths and cross sections: 
           IF (NSC .GT. 0) THEN                      ! only if x-sections exist 
              DO ISC = 1, OP_NS(J,K)                  ! loop over states

                 READ (XF,*) ISCX, NZ, NE, ISLP, ILV, XCHIC, NXC 

                 IF (NXC .NE. 0) THEN

                    WRITE (INSLP,'(I1,I3)') ILV,ISLP
                    DO IX = NXC, 1, -1                     ! loop over wavelengths
                       READ (XF,*) DUMXWL, DUMXS
                       ! write (*,*) dumxwl, dumxs
                       OPXWL(IX)=RINV/DUMXWL            ! Rydberg -> Angstrom
                       OPXXS(IX)=1.e-18*DUMXS           ! Mbarns -> cm^2
                    END DO

                    !! taking into account the duplicate wavelengths 
                    !! (shift wavelengths slightly so no overlapping points with different xsec)
                    DO IX = 1,NXC
                       IF (ABS(OPXWL(IX+1)-OPXWL(IX)) .LT. 1.D-35) THEN
                          IF (OPXXS(IX+1) .LT. 1.D-35) THEN
                             OPXWL(IX+1)=OPXWL(IX+1)*1.00001
                          ELSE IF (OPXXS(IX) .LT. 1.D-35) THEN
                             OPXWL(IX)  =OPXWL(IX)  *0.999995
                          ELSE
                             OPXWL(IX)  =OPXWL(IX)  *0.999995
                             OPXWL(IX+1)=OPXWL(IX+1)*1.000005
                          END IF
                       END IF
                    END DO

                    !           Find average OP cross-section value across 
                    !           Sterne wavelength interval
                    IF (NXC .GT. 0) CALL OP_XSEC(NXC,OPXWL,OPXXS,A_limit(ISC),LL,WL,DUMXS2)
                    OP_XS(1:LL,ISC,J,K) = DUMXS2(1:LL)
                    !write (*,*) wl
                    !write (*,*) opxwl(1:nxc),opxxs(1:nxc)
                    !write (*,*) isc,j,k,op_xs(:,isc,j,k)
                    !stop
                 ENDIF
              ENDDO                               ! end state loop
           ENDIF

           CLOSE(PF)
           CLOSE(XF)

        ENDDO                                   ! end ion loop
     ENDIF


  END DO                                     ! end atom loop

  !  Interpolate x-s at reference wavelength for opt.depth scale
  CALL OP_INXS 


END SUBROUTINE OP_BFIN
