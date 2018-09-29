* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE pitau(u,nmax,pi,tau)

*****************************************************************
* Calculates pi,n(u) and tau,n(u) with upward recursion         *
*                                                               *
*****************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER nmax
      DOUBLE PRECISION u, pi(nmax),tau(nmax)
C     DIMENSION pi(nmax),tau(nmax)

Cf2py intent(in) u, nmax
Cf2py intent(out) pi, tau

*     starting values:
      pi(1) = 1.D0
      pi(2) = 3.D0*u
      delta = 3.D0*u*u-1.d0
      tau(1)= u
      tau(2)= 2.D0*delta-1.d0

*     upward recursion:
      DO n=2,nmax-1
         pi(n+1) = dble(n+1)/dble(n)*delta + u*pi(n)
         delta   = u*pi(n+1) - pi(n)
         tau(n+1)= dble(n+1)*delta - pi(n)
      ENDDO

      RETURN
      END
 
