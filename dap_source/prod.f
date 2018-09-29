* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE prod(A,B,C,nmat,nmu,nmum)

*--------------------------------------------------------------------
*  Calculate the supermatrix product A = B * C  
*  Usually a large fraction of the execution time is spent in this
*  subroutine, especially with polarization.  
*--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER nsup,ng,nmu,nmum
C      INTEGER nsup,ng,nmu,nmum, nmat

C      DOUBLE PRECISION A(nsupMAX,nsupMAX), 
C     .          B(nsupMAX,nsupMAX), C(nsupMAX,nsupMAX)

      DIMENSION A(nsupMAX,nsupMAX), 
     .          B(nsupMAX,nsupMAX), C(nsupMAX,nsupMAX)

Cf2py intent(in,out) A, B
Cf2py intent(in) nmat, nmu, nmum

*--------------------------------------------------------------------
      nsup= nmu*nmat
      ng= nmum*nmat

      DO j=1,nsup
         DO i=1,nsup
            A(i,j)= 0.D0
            DO k=1,ng
               A(i,j)= A(i,j) + B(i,k)*C(k,j)
            ENDDO
         ENDDO
      ENDDO

*--------------------------------------------------------------------
      RETURN
      END

