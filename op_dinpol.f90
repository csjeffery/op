SUBROUTINE OP_DINPOL (X,Y,Z,XX,YY,ZZ)

  !     INTERPOLATION IN A TABLE OF TWO DIMENSIONS

  use op_precisions
  implicit none

  real :: X, Y, Z, YY
  real(kind=qq) :: XX
  real :: ZZ

  real :: A1, A2, A3, PN
  integer :: I

  DIMENSION X(2),Y(4),Z(4,2),A1(2),A2(2),A3(2),PN(2)
  !     INTERPOLATION OF THIRD DEGREE IN Y
  DO I = 1,2
     PN(I) = Z(1,I)
     A1(I) = (Z(1,I)-Z(2,I))/(Y(1)-Y(2))
     PN(I) = PN(I)+(YY-Y(1))*A1(I)
     A2(I) = (Z(2,I)-Z(3,I))/(Y(2)-Y(3))
     A1(I) = (A1(I)-A2(I))/(Y(1)-Y(3))
     PN(I) = PN(I)+(YY-Y(1))*(YY-Y(2))*A1(I)
     A3(I) = (Z(3,I)-Z(4,I))/(Y(3)-Y(4))
     A2(I) = ( A2(I)-A3(I))/(Y(2)-Y(4))
     A1(I) = (A1(I)-A2(I))/(Y(1)-Y(4))
     PN(I) = PN(I)+(YY-Y(1))*(YY-Y(2))*(YY-Y(3))*A1(I)
  ENDDO
  !     LINEAR INTERPOLATION IN X
  ZZ = PN(1)+(XX-X(1))*(PN(1)-PN(2))/(X(1)-X(2))
  RETURN

END subroutine OP_DINPOL
