SUBROUTINE OP_KIEL_RAY (WLS,JD,RAYH,RAYHE,RAYC,RAYN,RAYO)

  !+
  !
  !   Rayleigh scattering for H, He, C, N, O. 
  !   Constants for He, C, N, O
  !   from: Tarafdar \& Vardya, MNRAS 145, 171 (1969)
  !
  !-

  USE OP_PRECISIONS
  USE OP_GLOBAL, ONLY : OP_XUU

  IMPLICIT NONE

  INTEGER :: JD
  REAL, INTENT(IN) :: WLS
  REAL(KIND=QQ) :: WL2, WL4
  REAL :: RAYH,RAYHE,RAYC,RAYN,RAYO

  !.

  WL2 = WLS*WLS
  WL4 = WL2*WL2
  RAYH  = 0.0
  RAYHE = 0.0
  RAYC  = 0.0
  RAYN  = 0.0
  RAYO  = 0.0

  IF (WLS .GT. 584.) THEN 
     RAYHE = SC(5.2915e-14,0.470e+6,0.182e+12,WL2,WL4)/OP_XUU(1,2,JD)

     IF (WLS .GT. 1200.) THEN 
        RAYN = SC(2.319e-12,1.995e+6,3.362e+12,WL2,WL4)*4./OP_XUU(1,7,JD)

        IF (WLS .GT. 1215.) THEN 
           RAYH = SC(5.799e-13,2.452e+6,4.801e+12,WL2,WL4)*2./OP_XUU(1,1,JD)

           IF (WLS .GT. 1302.) THEN 
              RAYO = SC(0.7585e-12,0.997e+6,1.010e+12,WL2,WL4)*9./OP_XUU(1,8,JD)

              IF (WLS .GT. 1657.) THEN 
                 RAYC = SC(4.288e-12,2.671e+6,6.562e+12,WL2,WL4)*9./OP_XUU(1,6,JD)

              ENDIF
           ENDIF
        ENDIF
     ENDIF
  ENDIF

CONTAINS 

  REAL FUNCTION SC(a,b,c,x,y)

    IMPLICIT NONE

    REAL, INTENT(in) :: a, b, c
    REAL(kind=qq), INTENT(in)  :: x, y 

    SC = (a/y)*(1.0 + (b/x) + (c/y)) 

  END FUNCTION SC

END SUBROUTINE OP_KIEL_RAY
