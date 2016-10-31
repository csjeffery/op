  real FUNCTION OP_KIEL_O1(FREQ,T)


!  O I  opacity

      use op_precisions

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FREQ

      real(kind=qq) :: FREQ1
      real :: X911

      DATA FREQ1/0./
      DATA X911 /0./

 !! External function declarations
      real, external :: SEATON
!.

      IF (FREQ.NE.FREQ1) THEN
        X911 = 0.
        IF (FREQ.GE.3.28805E15) &
          & X911 = SEATON(FREQ,3.28805D15,2.9D-18,1.D00,2.66D00)
        FREQ1 = FREQ
      ENDIF

      OP_KIEL_O1 = X911

  END function OP_KIEL_O1
