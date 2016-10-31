      SUBROUTINE OP_IONFILL ( F_ION, U_ION, NK, NI, ND )

!+
!     
!     op_ionfill
!     
!     
!=     Fill ionisation fraction arrays throughout atmosphere
!     
!-

!     Global variables
      USE OP_IONS_GLOBAL

!     Global typing.
      IMPLICIT NONE

!     Array arguments
      INTEGER :: NK              ! dimension; number of elelemtns
      INTEGER :: NI              ! dimension: number of ions
      INTEGER :: ND              ! number of depth points
      REAL, DIMENSION(NK,NI,ND) :: F_ION ! ion fractions
      REAL, DIMENSION(NK,NI,ND) :: U_ION ! partition fractions

!.

      FH (1:2,1:ND) = F_ION( 1,1:2,1:ND)
      FHE(1:3,1:ND) = F_ION( 2,1:3,1:ND)
      FLI(1:3,1:ND) = F_ION( 3,1:3,1:ND)
      FBE(1:4,1:ND) = F_ION( 4,1:4,1:ND)
      FB (1:4,1:ND) = F_ION( 5,1:4,1:ND)
      FC (1:6,1:ND) = F_ION( 6,1:6,1:ND)
      FN (1:6,1:ND) = F_ION( 7,1:6,1:ND)
      FO (1:6,1:ND) = F_ION( 8,1:6,1:ND)
      FF (1:4,1:ND) = F_ION( 9,1:4,1:ND)
      FNA(1:4,1:ND) = F_ION(11,1:4,1:ND)
      FMG(1:6,1:ND) = F_ION(12,1:6,1:ND)
      FAL(1:4,1:ND) = F_ION(13,1:4,1:ND)
      FSI(1:6,1:ND) = F_ION(14,1:6,1:ND)
      FS (1:6,1:ND) = F_ION(16,1:6,1:ND)
      FAR(1:4,1:ND) = F_ION(18,1:4,1:ND)
      FCA(1:6,1:ND) = F_ION(20,1:6,1:ND)
      FFE(1:6,1:ND) = F_ION(26,1:6,1:ND)
      
      ZH (1:2,1:ND) = U_ION( 1,1:2,1:ND)
      ZHE(1:3,1:ND) = U_ION( 2,1:3,1:ND)
      ZLI(1:3,1:ND) = U_ION( 3,1:3,1:ND)
      ZBE(1:4,1:ND) = U_ION( 4,1:4,1:ND)
      ZB (1:4,1:ND) = U_ION( 5,1:4,1:ND)
      ZC (1:6,1:ND) = U_ION( 6,1:6,1:ND)
      ZN (1:6,1:ND) = U_ION( 7,1:6,1:ND)
      ZO (1:6,1:ND) = U_ION( 8,1:6,1:ND)
      ZF (1:4,1:ND) = U_ION( 9,1:4,1:ND)
      ZAR(1:4,1:ND) = U_ION(10,1:4,1:ND)
      ZNA(1:4,1:ND) = U_ION(11,1:4,1:ND)
      ZMG(1:6,1:ND) = U_ION(12,1:6,1:ND)
      ZAL(1:4,1:ND) = U_ION(13,1:4,1:ND)
      ZSI(1:6,1:ND) = U_ION(14,1:6,1:ND)
      ZS (1:6,1:ND) = U_ION(16,1:6,1:ND)
      ZAR(1:4,1:ND) = U_ION(18,1:4,1:ND)
      ZCA(1:6,1:ND) = U_ION(20,1:6,1:ND)
      ZFE(1:6,1:ND) = U_ION(26,1:6,1:ND)

      END SUBROUTINE OP_IONFILL
