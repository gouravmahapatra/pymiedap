* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE spline(x,y,n,nMAX,y2)

Cf2py intent(in) x,y,n,nMAX
Cf2py intent(out) y2

*----------------------------------------------------------------------------
*     Spline interpolation routine from Press et al. (1986, p.88). 
*
*     Given arrays x and y of length n containing a tabulated function,
*     i.e. y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 
*     and ypn for the first derivative of the interpolating function at
*     points 1 and n respectively, this routine returns an array y2 of  
*     length n which contains the second derivatives of the interpola- 
*     ting function at the tabulated points x(i).                     
*
*     If yp1 and/or yp2 are equal to 1x10^30 or larger, the routine is 
*     signalled to set the corresponding boundary condition for a natu- 
*     ral spline, with zero second derivative on that boundary.        
*
*     n is the number of elements in x and y
*     nMAX is the maximum number of elements in x and y
*----------------------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER nMAX, n, i, k

      DOUBLE PRECISION x(nMAX),y(nMAX),y2(nMAX),u(nMAX)
      DOUBLE PRECISION p, sig, qn, un

*----------------------------------------------------------------------------
      y2(1)= 0.D0
      u(1)=  0.D0

      DO i=2,n-1
      sig= (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p= sig*y2(i-1)+2.D0
      y2(i)= (sig-1.D0)/p
      u(i)= (6.D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      ENDDO

      qn= 0.D0
      un= 0.D0
      y2(n)= (un-qn*u(n-1))/(qn*y2(n-1)+1.D0)

      DO k=n-1,1,-1
      y2(k)= y2(k)*y2(k+1)+u(k)
      ENDDO

*----------------------------------------------------------------------------
      RETURN
      END
