      SUBROUTINE OP_FINTAB (TEMP,TABE,WL,WLS,T,ANEL,MM,LL)

!**** ZUSAMMENSTELLUNG DER FUER DIE INTERPOLATION BENOETIGTEN TABELLEN
!**** DER ABSORPTIONSKOEFFIZIENTEN

      use op_precisions

      implicit none

      real :: TEMP, TABE, WL
      real(kind=qq) :: WLS
      real, intent(in) :: T
      integer :: MM, LL
      
      real :: TEM, TAB
      real(kind=qq) :: ANEL

      integer :: I, L, M

      DIMENSION TEM(6),TAB(6,48),TABE(6,LL),TEMP(6,MM)
      DIMENSION WL(LL)

      DO 10  M = 1,MM
        IF (T .GT. TEMP(6,M))  GOTO 1
        DO 3  I = 1,6
          TEM(I) = TEMP(I,M)
          DO 3 L = 1,LL
 3          TAB(I,L) = TABE(I,L)
        CALL TABLE_OP(TEM,WL,TAB,LL,WLS,T,ANEL,-1)
          RETURN
 1        IF (T .GT. TEMP(1,M+1))         GOTO 10
        DO 4  I = 1,3
          TEM(I) = TEMP(I+3,M)
          DO 4        L = 1,LL
 4          TAB(I,L) = TABE(I+3,L)
        DO 5  I = 4,6
          TEM(I) = TEMP(I-3,M+1)
          DO 5        L = 1,LL
 5          TAB(I,L) = TABE(I-3,L)
        CALL TABLE_OP(TEM,WL,TAB,LL,WLS,T,ANEL,-1)
          RETURN
 10   CONTINUE
      RETURN

      END subroutine OP_FINTAB

