* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE star(A,B,nmat,nmu)

*----------------------------------------------------------------------
*  Calculate A = B* where B stands for a reflection or transmission    
*  supermatrix of a HOMOGENEOUS layer, and the star denotes            
*  illumination from below :  A = q4q3 B q3q4.                        
*  Eqs. (98)-(99) of de Haan et al. (1987)                            
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER nsup

      INCLUDE 'max_incl'

C      INTEGER nmat, nmu
C      DOUBLE PRECISION A, B
      DIMENSION A(nsupMAX,nsupMAX),B(nsupMAX,nsupMAX)

Cf2py intent(in,out) A, B
Cf2py intent(in) nmat, nmu

*----------------------------------------------------------------------
      nsup= nmu*nmat

      CALL assign(A,B,nmat,nmu)

      IF (nmat.GE.3) THEN
         DO mu=1,nmu
            ibase = (mu-1)*nmat
            DO k=3,nmat
               i= ibase+k
               DO j=1,nsup
                  A(i,j) = -A(i,j)
                  A(j,i) = -A(j,i)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

*----------------------------------------------------------------------
      RETURN
      END
