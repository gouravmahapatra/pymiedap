* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE fillup(Rm,Tm,nmat,nmu,nsup)

*----------------------------------------------------------------------*
*  Fill the upper triangle of a supermatrix using symmetry relations   *
*  for a vertically homogeneous layer.                                 *
*  By the upper triangle we mean mu < mu0 : nmat by nmat submatrices   *
*  are not split.                                                      *
*  See de Haan et al. (1987): Astron. Astrophys. 183, p. 371           *
*  Eqs. (96)-(97)                                                      *
*                 R = q3 R~ q3                                         *
*                 T = q4 T~ q4                                         *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      REAL*8, DIMENSION(nsup,nsup) :: Rm, Tm
      REAL*8, DIMENSION(4) :: q3, q4

Cf2py intent(in,out) Rm, Tm
      
*-----------------------------------------------------------------------
      DO k=1,4
         q3(k)= 1.D0
         q4(k)= 1.D0
      ENDDO
      q3(3)= -1.D0
      q4(4)= -1.D0

      DO mu0=1,nmu
         jbase= (mu0-1)*nmat
         DO mu=1,mu0-1
            ibase= (mu-1)*nmat
            DO ki=1,nmat
               i= ibase+ki
               DO kj=1,nmat
                  j= jbase+kj
                  Rm(i,j)= q3(ki)*Rm(j,i)*q3(kj)
                  Tm(i,j)= q4(ki)*Tm(j,i)*q4(kj)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

*-----------------------------------------------------------------------
      RETURN
      END
