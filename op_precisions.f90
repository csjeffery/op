MODULE OP_PRECISIONS

  !+

  !  Name:
  !     op_precisions

  !  Purpose:
  !=    Declares global kind parameters for lte-codes

  !  Module Type:
  !     Module

  !  Language:
  !     Fortran 90

  !  Authors:
  !     NTB: Natalie Behara
  !     CSJ: Simon Jeffery
  !
  !  History:
  !     21-May-2007 (CSJ): implemented as lte-precisions
  !
  !  Bugs:
  !     Needs to be moved to a location suitable for use with all lte-codes
  !
  !-

  IMPLICIT NONE

  INTEGER, PARAMETER :: LTE_SP = SELECTED_REAL_KIND(P=8,R=99)
  INTEGER, PARAMETER :: LTE_DP = SELECTED_REAL_KIND(P=14,R=99)

  INTEGER, PARAMETER :: QQ = LTE_SP
  INTEGER, PARAMETER :: RR = LTE_DP

END MODULE OP_PRECISIONS
