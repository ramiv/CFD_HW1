      SUBROUTINE TRIDIAG(A,B,C,D,U,N,IS,IE)
!	  Tri-diagonal matrix Solver
!
!     A, B, C, are the vectors of the diagonal and the two off-diagonals.
!     The vector D is the RHS vector, the vector U is the solution
!     vector, N is the dimension, IS is the starting point, and IE is
!     the last point.
!
      INTEGER :: N,IS,IE
      !REAL A(0:N),B(0:N),C(0:N),D(0:N),U(0:N)
      REAL A(IS:IE),B(IS:IE),C(IS:IE),D(IS:IE),U(IS:IE)
      real :: BETA
      INTEGER :: I
      DO I = IS + 1,IE
         IF(B(I-1).EQ.0.) then
            write(0,*) "Rewrite your equations "
            exit
         END IF
         BETA = A(I) / B(I-1)
         B(I) = B(I) - C(I-1) * BETA
         D(I) = D(I) - D(I-1) * BETA
      END DO
      U(IE) = D(IE) / B(IE)
      DO I = IE - 1, IS, -1
         U(I) = (D(I) - C(I) * U(I+1)) / B(I)
      END DO
      RETURN
      END
