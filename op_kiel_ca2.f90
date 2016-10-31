  real FUNCTION OP_KIEL_CA2(FREQ,T)

!  Ca II opacity

      use op_precisions

      implicit none

      real, intent(in) :: T
      real(kind=qq), intent(in) :: FREQ

      real(kind=qq) :: FREQ1
      real :: X1420, X1218, X1044
      real :: TKEV, C1218, C1420


 !! External function declarations
      real, external :: SEATON

      DATA FREQ1                / 0./
      DATA X1420, X1218, X1044  /0., 0., 0./
!.

      TKEV = T * 8.6167E-5
      C1218 = 10. * EXP(-1.697/TKEV)
      C1420 = 6. * EXP(-3.142/TKEV)

      IF (FREQ.NE.FREQ1) THEN
        X1420 = 0.
        X1218 = 0.
        X1044 = 0.
        IF (FREQ.GE.2.110779D15) THEN
           X1420 = SEATON(FREQ,2.110779D15,4.13D-18,3.D00,.69D00)
          IF (FREQ.GE.2.460127E15) THEN
            X1218 = 1.64E-17 * SQRT(2.460127E15/FREQ)
            IF (FREQ.GE.2.870454E15) THEN
              X1044 = 5.4E-20* (2.870454E15/FREQ)**3
            ENDIF
          ENDIF
        ENDIF
        FREQ1 = FREQ
      ENDIF

      OP_KIEL_CA2 =(X1044*2. + X1218*C1218 + X1420*C1420)/2.

  end function OP_KIEL_CA2
