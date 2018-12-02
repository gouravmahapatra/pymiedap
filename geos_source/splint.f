* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE splint(xa,ya,y2a,n,x,y,nn)

Cf2py intent(in) xa, ya, y2a, n, x, nn
Cf2py intent(out) y

*----------------------------------------------------------------------------
*     Spline interpolation routine from Press et al. (1986, p.88).   
*                                                                 
*     Given the arrays xa and ya of length n, which tabulate a function
*     (with the xa(i)'s in order), and given the array y2a, which is  
*     the output from SPLINE above, and given a value of x, this     
*     routine returns a cubic-spline interpolated value y.       
*----------------------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER n, nn, klo, khi,k
      DOUBLE PRECISION xa(nn),ya(nn),y2a(nn)
      DOUBLE PRECISION x, y, a, b, h

*----------------------------------------------------------------------------
      klo=1
      khi=n

1     IF (khi-klo.GT.1) THEN
        k= (klo+khi)/2  
         IF (xa(k).GT.x) THEN
             khi= k
         ELSE
        klo= k
         ENDIF
         GOTO 1
      ENDIF

      h= xa(khi)-xa(klo)

      IF (DABS(h).LT.1.D-10) WRITE(*,10) 
      a= (xa(khi)-x)/h
      b= (x-xa(klo))/h
      y= a*ya(klo)+b*ya(khi)+
     +  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.D0

*----------------------------------------------------------------------------
10    FORMAT('ERROR in routine SPLINT: Bad XA input.')

      RETURN
      END
