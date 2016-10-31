!
!        Subroutine to provide ground state statistical weights
!
  SUBROUTINE GSTATW (IZ,G0,NG)

!
!     Input:  IZ      Atomic number of species for which g0 are required
!
!     Output: G0(*)   Array containing ground state statistical weights 
!             NG      Number of ionisation stages (including neutrals)
!
      implicit none

      integer :: IZ, NG
      integer :: I, J, NI
      real :: G0(NG)

      INTEGER*2 IG0(6,26)

      DATA  ((IG0(I,J),I=1,6),J=1,26)&
       &/  2,  1,  1,  1,  1,  1,    1,  2,  1,  1,  1,  1,&   !  H , He
        &  2,  1,  2,  1,  1,  1,    1,  2,  1,  2,  1,  1,&   !  Li, Be
        &  6,  1,  2,  1,  2,  1,    9,  6,  1,  2,  1,  2,&   !  B , C
        &  4,  9,  6,  1,  2,  1,    9,  4,  9,  6,  1,  2,&   !  N , O
        &  6,  9,  4,  9,  6,  1,    1,  6,  9,  4,  9,  6,&   !  F , Ne
        &  2,  1,  6,  9,  4,  9,    1,  2,  1,  6,  9,  4,&   !  Na, Mg
        &  6,  1,  2,  1,  6,  9,    9,  6,  1,  2,  1,  6,&   !  AL, Si
        &  4,  9,  6,  1,  2,  1,    9,  4,  9,  6,  1,  2,&   !  P , S
        &  6,  9,  4,  9,  6,  1,    1,  6,  9,  4,  9,  6,&   !  Cl, Ar
        &  2,  1,  6,  9,  4,  9,    1,  2,  1,  6,  9,  4,&   !  K , Ca
        & 10, 15, 10,  1,  1,  1,   21, 28, 21,  1,  1,  1,&   !  Sc, Ti
        & 28, 25, 28,  1,  1,  1,    7,  6, 25,  1,  1,  1,&   !  V , Cr
        &  6,  7,  6,  1,  1,  1,   25, 30, 25,  6,  1,  1 /  !  Mn, Fe

      IF (IZ .GT. 26) THEN
         DO I = 1,NG
           G0(I) = 1.0
         end do
      else
         NI = MIN (NG,6)
         DO I = 1,NI
           G0(I) = IG0(I,IZ)
         end do
      ENDIF


  END subroutine GSTATW
