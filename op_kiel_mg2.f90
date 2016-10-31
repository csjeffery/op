  real FUNCTION OP_KIEL_MG2(FREQ,T)


!  Mg II opacity

      use op_precisions

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FREQ

      real(kind=qq) :: FREQ1
      real :: X1169, X824
      real :: TKEV, C1169

      DATA FREQ1             / 0./
      DATA X1169, X824       /0., 0./

 !! External function declarations
      real, external :: SEATON

!.
      TKEV = T * 8.6167E-5
      C1169 = 6. * EXP(-4.43/TKEV)

      IF (FREQ.NE.FREQ1) THEN 
        X1169 = 0.
        X824 = 0.
        IF (FREQ.GE.2.564306D15) THEN 
          X1169 = 5.11E-19 * (2.564306E15/FREQ)**3
          IF (FREQ.GE.3.635492E15) THEN
            X824 = SEATON(FREQ,3.635492D15,1.4D-19,4.D00,6.7D00)
          ENDIF
        ENDIF
        FREQ1 = FREQ
      ENDIF

      OP_KIEL_MG2 = (X824*2. + X1169*C1169) / 2.

  END function OP_KIEL_MG2
