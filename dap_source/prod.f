* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE prod(A,B,C,nmat,nmu,nmum)

*--------------------------------------------------------------------
*  Calculate the supermatrix product A = B * C  
*  Usually a large fraction of the execution time is spent in this
*  subroutine, especially with polarization.
*  Edited by: Ashwyn Groot
*  Date: November 2018
*  Introduced matrix operations with f95<
*--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER nsup,ng,nmu,nmum

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Bsub,Csub !rank 2

      REAL*8, DIMENSION(nmu*nmat,nmu*nmat) :: A,B,C !rank 2

      nsup= nmu*nmat
      ng= nmum*nmat


      ALLOCATE(Bsub(nsup,ng), Csub(ng,nsup))
C     .          A(nsup,nsup), B(nsup,nsup), C(nsup,nsup))

Cf2py intent(in,out) A, B
Cf2py intent(in) nmat, nmu, nmum

*--------------------------------------------------------------------
      Bsub=B(:,1:ng)
      Csub=C(1:ng,:)
      A=MATMUL(Bsub,Csub)
*--------------------------------------------------------------------
      RETURN
      END

