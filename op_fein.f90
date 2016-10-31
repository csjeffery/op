  SUBROUTINE OP_FEIN

 !+
 !! Reads in Iron Project data for FE I, II and III continuous opacities.
 !!  History:
 !!     05-Jul-05: Created (NTB)
 !-

    USE LTE_WLGRID
    USE OP_PRECISIONS
    USE OP_GLOBAL

 !.

    IMPLICIT NONE 

 !. 
    INTEGER :: I, J, K, IX

  !  Local variables.
    INTEGER :: FLAG, count, count2

    INTEGER :: NSC, ISC             !  no of levels in ion (+counter)
    INTEGER :: NXC                  !  no of wavelengths in level
    REAL :: OPXWL(OP_MAXNX)      !  OP x-s wavelengths (A)
    REAL :: OPXXS(OP_MAXNX)      !  OP x-s (cm^2)

    REAL :: IS, IL, IP, NS

    INTEGER :: DUM1, DUM2, DUM4
    REAL :: DUMIE, DUM5
    REAL :: DUMXWL, DUMXS, DUMXS2(OP_MAXNW)
    REAL :: A_LIMIT(0:OP_MAXNX)

  ! unit numbers for files
    INTEGER :: XF=12, XF2=13

  !  Local constants.
    REAL(KIND=QQ), PARAMETER :: RINV = 911.2670747D0 ! 1/R (A)
    REAL(KIND=QQ), PARAMETER :: REV = 13.60559024D0  ! 1 Ryd (eV)

 !! locator for OP data
    CHARACTER(LEN=7), PARAMETER :: PATH='FEDATA/'       
 !! files holding xesc data
    CHARACTER(LEN=33) :: FE1aFILE, FE1bFILE, FE2FILE, FE3FILE

 !.

    FE1aFILE = PATH // 'px.fe1.spt'
    FE1bFILE = PATH // 'px.fe1.qnt'
    FE2FILE = PATH // 'px.fe2'
    FE3FILE = PATH // 'px.fe3.cropped'

    K=26                       ! atomic number of current atom

 ! Loop over ions
    DO J=1,3


      IF (J==1) THEN
        OPEN(UNIT=XF2,FILE=FE1aFILE,STATUS='UNKNOWN')
        OPEN(UNIT=XF,FILE=FE1bFILE,STATUS='UNKNOWN')
      ELSEIF (J==2) THEN
        OPEN(UNIT=XF,FILE=FE2FILE,STATUS='UNKNOWN')
      ELSE
        OPEN(UNIT=XF,FILE=FE3FILE,STATUS='UNKNOWN')
      ENDIF

 

 ! Read in weights and ionisation energies:
      FLAG=1
      NSC=0
      A_limit = 0.0 

 ! count the states as we go
      READ (XF,*) DUM1,DUM2


      DO COUNT=1,OP_MAXNS+1
        IF (FLAG==1) THEN
          READ(XF,*) IS, IL, IP, NS
          IF (((IS == 0.0) .AND. (IL == 0.0) .AND. (IP == 0.0)) .OR. (COUNT == OP_MAXNS+1)) THEN  ! end of file
            FLAG=0
            OP_NS(j,k)=NSC     
        !    write(*,*) 'number of states', k, j, NSC 
          ELSE
            READ(XF,*) DUM4, NXC
            NSC=NSC+1
            OP_XWT(NSC,j,k) = is*((2*il) +1)     ! stat. weight
        ! store ionisation limit, angstrom
            IF (J == 1) A_limit(NSC) = 1575.39
            IF (J == 2) A_limit(NSC) = 767.22
            IF (J == 3) A_limit(NSC) = 404.50
            READ(XF,*) DUMIE, DUM5
            OP_XCHI(NSC,j,k)= DUMIE * REV       ! Rydberg -> eV, energy rel. to ground
    ! read in wavelengths and xsecs:
            DO IX = NXC, 1, -1
              READ (XF,*) DUMXWL, DUMXS
              OPXWL(IX)=RINV/DUMXWL            ! Rydberg -> Angstrom
              OPXXS(IX)=1.e-18*DUMXS           ! Mbarns -> cm^2
            ENDDO

 !! taking into account the duplicate wavelengths 
 !! (shift wavelengths slightly so no overlapping points with different xsec)
              DO IX = 1,NXC
                IF (ABS(OPXWL(IX+1)-OPXWL(IX)) < 1.D-35) THEN
                  IF (OPXXS(IX+1) < 1.D-35) THEN
                    OPXWL(IX+1)= OPXWL(IX+1)*1.00001
                  ELSE IF (OPXXS(IX) < 1.D-35) THEN
                    OPXWL(IX)  = OPXWL(IX)  *0.999995
                  ELSE
                    OPXWL(IX)  = OPXWL(IX)  *0.999995
                    OPXWL(IX+1)= OPXWL(IX+1)*1.000005
                  ENDIF
                ENDIF
              ENDDO
              CALL OP_XSEC(NXC,OPXWL,OPXXS,A_limit(NSC),LL,WL,DUMXS2)
              DO I=1,LL
                OP_XS(I,NSC,J,K)=DUMXS2(I)
              ENDDO

          ENDIF
        ENDIF
      ENDDO


     CLOSE(XF)
     IF (J == 1) CLOSE(XF2)

   ENDDO                                   ! end ion loop


  !  Interpolate x-s at reference wavelength for opt.depth scale
    CALL OP_INXS 


  END SUBROUTINE OP_FEIN
