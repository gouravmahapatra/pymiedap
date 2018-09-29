* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE assign(A,B,nmat,nmu)

*----------------------------------------------------------------------
*     Calculate the supermatrix assignment A=B:
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nsup,i,j,nmat,nmu

      DOUBLE PRECISION A(nsupMAX,nsupMAX),B(nsupMAX,nsupMAX)

Cf2py intent(in,out) A,B
Cf2py intent(in) nmat,nmu

*----------------------------------------------------------------------
      nsup= nmat*nmu

      DO j=1,nsup
         DO i=1,nsup
            A(i,j)= B(i,j)
         ENDDO
      ENDDO

*----------------------------------------------------------------------
      RETURN
      END
