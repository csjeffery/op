
      SUBROUTINE OP_LINPOL (X,Y,XX,YY,TAB,TABXY,NSPA,NZEI)

! **** LINEARE INTERPOLATION IN EINER ZWEIDIM. TABELLE ****

      use op_precisions

      implicit none

      real :: X, Y
      real, intent(in) :: YY, TAB
      real, intent(in) :: XX
      real, intent(inout) :: TABXY
      integer :: NSPA
      integer :: NZEI

      integer :: I, ISPA, IZEI
      real :: TABY1, TABY2, XR, YR

      DIMENSION         TAB(NZEI,NSPA),X(NSPA),Y(NZEI)

! ****        POSITION IN X (
      DO 1  I = 2,NSPA
        IF (X(1) .GT. X(NSPA))        GOTO 10
        IF (XX-X(I)) 2,1,1
 10        IF (XX - X(I)) 1,2,2
1     CONTINUE
      I = NSPA
2     ISPA = I

! ****        POSITION IN Y (
      DO 3  I = 2,NZEI
         IF (YY - Y(I))         4,3,3
3     CONTINUE
      I = NZEI
4     IZEI = I

! ****        INTERPOLATION (
      XR = (XX - X(ISPA-1))/(X(ISPA) - X(ISPA-1))
      YR = (YY - Y(IZEI-1))/(Y(IZEI) - Y(IZEI-1))
      TABY1 = TAB(IZEI-1,ISPA-1) +  (TAB(IZEI-1,ISPA) - &
                  & TAB(IZEI-1,ISPA-1))*XR
      TABY2 = TAB(IZEI,ISPA-1) +  (TAB(IZEI,ISPA) - TAB(IZEI,ISPA-1))* &
                 &  XR
      TABXY = TABY1 + (TABY2 - TABY1)*YR
      RETURN

      END subroutine OP_LINPOL
