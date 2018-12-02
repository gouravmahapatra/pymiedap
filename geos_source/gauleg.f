* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE gauleg(ndim,ngauss,a,b,x,w)

**********************************************************************
* Given the lower and upper limits of integration a and b, 
* and given the number of Gauss-Legendre points ngauss,
* this routine returns through array x the abscissas and through
* array w the weights of the Gauss-Legendre quadrature formula.
* Eps is the desired accuracy of the abscissas.
* This routine is documented further in:
*
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)
*   page 125 ISBN 0-521-30811-9                                   
**********************************************************************
C      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER ndim, ngauss, i
      DOUBLE PRECISION a,b,xm,xl,z,p1,p2,p3,pp,z1,pi,eps
      PARAMETER (eps= 1.d-14)

      REAL*8, DIMENSION(ndim) :: x, w

Cf2py intent(in) ndim,ngauss,a,b
Cf2py intent(out) x, w

**********************************************************************
      pi=4.D0*datan(1.D0)
      m=(ngauss+1)/2
      xm=0.5D0*(a+b)
      xl=0.5D0*(b-a)

      DO 12 i=1,m
         z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
1        CONTINUE
            p1=1.D0
            p2=0.D0
            DO j=1,ngauss
               p3= p2
               p2= p1
               p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
            ENDDO
            pp=ngauss*(z*p1-p2)/(z*z-1.d0)
            z1= z
            z= z1-p1/pp
          IF (dabs(z-z1).GT.eps) GOTO 1
          x(i)= xm-xl*z
          x(ngauss+1-i)= xm+xl*z
          w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
          w(ngauss+1-i)= w(i)
12    CONTINUE

**********************************************************************
      RETURN
      END

