* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE rdiapr(A,B,E,nmat,nmu)

*----------------------------------------------------------------------*
*  Calculate the full EXTENDED supermatrix product A = B E where       *
*  E is diagonal. (rdiapr = 'right diagonal product')                  *
*  The reason why the product is not limited to the integration        *
*  points is explained below Eq. (95) of de Haan et al. (1987).        *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER nmat, nmu
      DOUBLE PRECISION A(nsupMAX,nsupMAX),B(nsupMAX,nsupMAX),E(nsupMAX)
C      DIMENSION A(nsupMAX,nsupMAX),B(nsupMAX,nsupMAX),E(nsupMAX)

Cf2py intent(out) A

*----------------------------------------------------------------------*
      nsup= nmu*nmat

      DO j=1,nsup
         DO i=1,nsup
            A(i,j) = B(i,j)*E(j)
         ENDDO
      ENDDO

*----------------------------------------------------------------------*
      RETURN
      END
