* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE transf(S,nsup,nmat)

************************************************************
*  Calculate the product D1 * SM * D1 where D1 is a 
*  supermatrix containing in each submatrix :             
*         ( 1   0   0   0 )          ( 1   0   0 )       
*         ( 0   1   1   0 )    or    ( 0   1   1 )    
*         ( 0   1  -1   0 )          ( 0   1  -1 )   
*         ( 0   0   0   1 )                        
************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
 
      INCLUDE 'max_incl'

C      INTEGER nmat, nsup
C      DOUBLE PRECISION S(nsupMAX,nsupMAX)
      DIMENSION S(nsupMAX,nsupMAX)

Cf2py intent(in) nsup, nmat
Cf2py intent(in,out) S

************************************************************
      IF ((nmat.NE.3) .AND. (nmat.NE.4)) RETURN
      DO i=2,nsup,nmat
         DO j=1,nsup
            s1 = S(i,j)
            s2 = S(i+1,j)
            S(i,j)   = s1+s2
            S(i+1,j) = s1-s2
         ENDDO
         DO j=1,nsup
            s1 = S(j,i)
            s2 = S(j,i+1)
            S(j,i)   = s1+s2
            S(j,i+1) = s1-s2
         ENDDO
      ENDDO

************************************************************
      RETURN
      END
