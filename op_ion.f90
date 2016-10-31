SUBROUTINE OP_ION (K,T,PE,U,X)

  !+
  !
  !   Name:
  !   OP_ION
  !
  !   Purpose:
  !   Calculate fractional abundances of ions of a given atomic species
  !   from the Saha ionization equilibrium formula. This formulation is
  !   not currently suitable for use with -ve ions, eg H and He.
  !
  !
  !       Input:
  !       K        Atomic number (Z) of species required
  !       T        Temperature (degrees K)
  !       PE       Electron pressure
  !       U(6)     Partition functions for first 6 ionisation states (I - VI)
  !
  !       Output:
  !       X(6)     Fractional abundance of ions in each state (normalized to 1)
  !                NB: the contribution by negative ions is neglected
  !
  !       Data supplied:
  !       CHI(5,56) Ionisation energies (eV) for first five ionisation states
  !                 of atoms from H to Ba (Z = 1 to 56)
  !
  !-

  ! Implicit typing.
  IMPLICIT NONE

  !  Subroutine Arguments (input).
  INTEGER :: K      ! Atomic number
  REAL :: T         ! Temperature (K)
  REAL :: PE        ! Electron pressure (dyn cm^-2)
  REAL :: U(6)      ! Partition functions

  !  Subroutine Arguments (output).
  REAL ::  X(6)     ! Ionization fractions

  !  Local variables.
  REAL :: F(6)             ! Ionization ratios
  REAL :: CHI(5,56)        ! Ionization potentials
  REAL :: XCHI

  INTEGER :: I, J, M

  REAL :: SUM

  REAL :: SAHA, BKDT

  !  Local constants.
  REAL, PARAMETER :: S0=0.666698474
  REAL, PARAMETER :: BKX = 1./8.6167E-05

  !  Extrrnal functions
  DOUBLE PRECISION, EXTERNAL :: TAP_IONPT

  !  Initialize ionization potentials:

  DATA  ((CHI(i,k),i=1,5),k=1,14)&
       &/13.598,  0.0  ,   0.0  ,   0.0  ,   0.0 ,&   ! H
       &24.587, 54.416,   0.0  ,   0.0  ,   0.0  ,&   ! He
       & 5.392, 75.638, 122.451,   0.0  ,   0.0  ,&   ! Li
       & 9.322, 18.211, 153.893, 217.713,   0.0  ,&   ! Be
       & 8.298, 25.154,  37.930, 259.368, 340.217,&   ! B
       &11.260, 24.383,  47.887,  64.492, 392.077,&   ! C
       &14.534, 29.601,  47.448,  77.472,  97.888,&   ! N
       &13.618, 35.116,  54.934,  77.412, 113.896,&   ! O
       &17.422, 34.970,  62.707,  87.138, 114.240,&   ! F
       &21.564, 40.962,  63.45 ,  97.11 , 126.21 ,&   ! Ne
       & 5.139, 47.286,  71.64 ,  98.91 , 138.39 ,&   ! Na
       & 7.646, 15.035,  80.143, 109.24 , 141.26 ,&   ! Mg
       & 5.986, 18.828,  28.447, 119.99 , 153.71 ,&   ! Al
       & 8.151, 16.345,  33.492,  45.141, 166.77 /    ! Si

  DATA  ((CHI(i,k),i=1,5),k=15,28)&
       &/10.486, 19.725,  30.18 ,  51.37 ,  65.023,&   ! P
       & 10.360, 23.33 ,  34.83 ,  47.30 ,  72.68 ,&   ! S
       & 12.967, 23.81 ,  39.61 ,  53.46 ,  67.8  ,&   ! Cl
       & 15.759, 27.629,  40.74 ,  59.81 ,  75.02 ,&   ! Ar
       &  4.341, 31.625,  45.72 ,  60.91 ,  82.66 ,&   ! K
       &  6.113, 11.871,  50.908,  67.10 ,  84.41 ,&   ! Ca
       &  6.54 , 12.80 ,  24.76 ,  73.47 ,  91.66 ,&   ! Sc
       &  6.82 , 13.58 ,  27.491,  43.266,  99.22 ,&   ! Ti
       &  6.74 , 14.65 ,  29.310,  46.707,  65.23 ,&   ! V
       &  6.766, 16.50 ,  30.96 ,  49.1  ,  69.3  ,&   ! Cr
       &  7.453, 15.640,  33.667,  51.2  ,  72.4  ,&   ! Mn
       &  7.870, 16.18 ,  30.651,  54.8  ,  75.0  ,&   ! Fe
       &  7.86 , 17.06 ,  33.50 ,  51.3  ,  79.5  ,&   ! Co
       &  7.635, 18.168,  35.17 ,  54.9  ,  75.5  /    ! Ni

  DATA  ((CHI(i,k),i=1,5),k=29,42)&
       &/ 7.726, 20.292,  36.83 ,  55.2  ,  79.9  ,&   ! Cu
       &  9.394, 17.964,  39.722,  59.4  ,  82.6  ,&   ! Zn
       &  5.999, 20.51 ,  30.71 ,  64.   ,  87.   ,&   ! Ga
       &  7.899, 15.934,  34.22 ,  45.71 ,  93.5  ,&   ! Ge
       &  9.81 , 18.633,  28.351,  50.13 ,  62.63 ,&   ! As
       &  9.752, 21.19 ,  30.820,  42.944,  68.3  ,&   ! Se
       & 11.814, 21.8  ,  36.   ,  47.3  ,  79.7  ,&   ! Br
       & 13.999, 24.359,  36.95 ,  52.5  ,  64.7  ,&   ! Kr
       &  4.177, 27.28 ,  40.   ,  52.6  ,  71.0  ,&   ! Rb
       &  5.695, 11.030,  43.6  ,  57.   ,  71.6  ,&   ! Sr
       &  6.38 , 12.24 ,  20.52 ,  61.8  ,  77.0  ,&   ! Y
       &  6.84 , 13.13 ,  22.99 ,  34.34 ,  81.5  ,&   ! Zr
       &  6.88 , 14.32 ,  25.04 ,  38.3  ,  50.55 ,&   ! Nb
       &  7.099, 16.15 ,  27.16 ,  46.4  ,  61.2  /    ! Mo

  DATA  ((CHI(i,k),i=1,5),k=43,56)&
       &/ 7.28 , 15.26 ,  29.54 ,  46.   ,  55.   ,&   ! Tc
       &  7.37 , 16.76 ,  28.47 ,  50.   ,  60.   ,&   ! Ru
       &  7.46 , 18.08 ,  31.06 ,  48.   ,  65.   ,&   ! Rh
       &  8.34 , 19.43 ,  32.92 ,  53.   ,  62.   ,&   ! Pd
       &  7.576, 21.49 ,  34.83 ,  56.   ,  68.   ,&   ! Ag
       &  8.993, 16.908,  37.48 ,  59.   ,  72.   ,&   ! Cd
       &  5.786, 18.869,  28.03 ,  54.4  ,  77.   ,&   ! In
       &  7.344, 14.632,  30.502,  40.734,  72.28 ,&   ! Sn
       &  8.641, 16.53 ,  25.3  ,  44.2  ,  56.   ,&   ! Sb
       &  9.009, 18.6  ,  27.96 ,  37.41 ,  58.75 ,&   ! Te
       & 10.451, 19.131,  33.   ,  42.   ,  66.   ,&   ! I
       & 12.130, 21.21 ,  32.1  ,  46.   ,  57.   ,&   ! Xe
       &  3.894, 25.1  ,  35.   ,  46.   ,  62.   ,&   ! Cs
       &  5.212, 10.004,  39.   ,  49.   ,  62.   /    ! Ba

  !.

  SAHA = S0  * T**2.5
  BKDT = BKX / T

  F = 0.

  DO I = 1,MIN(K,5)
     XCHI = TAP_IONPT (K, I-1)
     IF ( XCHI < 0.001 .AND. K < 57 ) XCHI = CHI(I,K)
     IF ( XCHI < 0.001 .AND. K > 56 ) XCHI = I + 10.* K
     ! F(I) = SAHA*U(I+1)/U(I)*EXP(-CHI(I,K)*BKDT)
     F(I) = SAHA*U(I+1)/U(I)*EXP(-XCHI*BKDT)
  END DO

  X(2) = 1.0 / (PE/F(1) + 1.0 + F(2)*(F(3)+PE)/(PE*PE))
  X(3) = X(2)*F(2)/PE
  X(4) = X(3)*F(3)/PE
  X(1) = X(2)*PE/F(1)
  X(5) = 0.
  X(6) = 0.

  !  find 5th and 6th ionization stages...
  IF (K>5) THEN
     X(5) = X(4)*F(4)/PE
     X(6) = X(5)*F(5)/PE
     SUM  = X(1) + X(2) + X(3) + X(4) + X(5) + X(6)

     !    ...by iteration
     DO I = 1,10
        X(1) = X(1)/SUM
        X(2) = X(1)*F(1)/PE
        X(3) = X(2)*F(2)/PE
        X(4) = X(3)*F(3)/PE
        X(5) = X(4)*F(4)/PE
        X(6) = X(5)*F(5)/PE
        SUM  = X(1) + X(2) + X(3) +X(4) + X(5) + X(6)
        IF (ABS(1.- SUM) .LT. 0.0001)  EXIT
     ENDDO

     IF (ABS(1.- SUM) > 0.0001) THEN
        WRITE(*,&
             &'(''OP_ION: Ionization equilibrium not converging ='',F6.4)')& 
             & SUM
        WRITE(*,*) 'atom',K,X(1),X(2),X(3),X(4),X(5),X(6)
        WRITE(*,*) 'PE', PE
     ENDIF

  ENDIF

END SUBROUTINE OP_ION
