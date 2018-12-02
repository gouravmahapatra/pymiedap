* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE trace(A,nmat,nmum,nsup,trA)

*---------------------------------------------------------------------
*  Calculate the truncated supermatrix trace(A).                    
*  The sum runs over the integration points only, see remark below
*  Eq. (124) of de Haan et al. (1987)                          
*---------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,ng,nmum,nmat,nsup
      REAL*8 trA

      REAL*8, DIMENSION(nsup,nsup) :: A !rank 2

      nsup= (nmum+1)*nmat

Cf2py intent(out) trA

*---------------------------------------------------------------------
      ng = nmum*nmat

      trA = 0.D0
      DO i=1,ng
         trA = trA + A(i,i)
      ENDDO

*---------------------------------------------------------------------
      RETURN
      END
