* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE sizedis(idis,par1,par2,par3,weight2,r,numr,nwithr)

************************************************************************
*  Calculate the size distribution n(r) for the numr radius values     *
*  contained in array r and RETURN the results through the array nwithr*
*  The size distributions are normalized such that the integral over   *
*  all r is equal to one.                                              *
************************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER idis, numr
      DOUBLE PRECISION par1, par2, par3, weight2
C      DOUBLE PRECISION nwithr, logC, logC1, logC2
      DOUBLE PRECISION logC, logC1, logC2
      DOUBLE PRECISION r(numr),nwithr(numr)

C     DIMENSION        r(numr),nwithr(numr)

Cf2py intent(in) idis,par1,par2,par3,weight2,r,numr
Cf2py intent(out) nwithr

      pi     = dacos(-1.d0)
      root2p = dsqrt(pi+pi)

************************************************************************
      IF (idis .EQ. 0) RETURN
      GOTO (10, 20, 30, 40, 50, 60, 70, 80 ) idis
      WRITE(*,*) ' sizedis: illegal index : ',idis
      STOP 'illegal index in sizedis'

************************************************************************
*  1 TWO PARAMETER GAMMA with alpha and b given                        *
************************************************************************
   10 alpha = par1
      b     = par2
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammln(alpha1)
      DO i=1,numr
         nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
      ENDDO
      GOTO 999

************************************************************************
*  2 TWO PARAMETER GAMMA with par1= reff and par2= veff given          *
************************************************************************
   20 alpha = 1.D0/par2 - 3.D0
      b     = 1.D0/(par1*par2)
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammln(alpha1)
      DO i=1,numr
         nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
      ENDDO
      GOTO 999

************************************************************************
*  3 BIMODAL GAMMA with equal mode weights                             *
************************************************************************
   30 alpha = 1.D0/par3 - 3.D0
      b1    = 1.D0/(par1*par3)
      b2    = 1.D0/(par2*par3)
      gamlna= gammln(alpha+1.D0)
      logC1 = (alpha+1.D0)*dlog(b1)-gamlna
      logC2 = (alpha+1.D0)*dlog(b2)-gamlna
      DO i=1,numr
         nwithr(i) = (1.D0/(1.D0 + weight2)) * 
     +               ( dexp(logC1+alpha*dlog(r(i))-b1*r(i)) +
     +               weight2 * dexp(logC2+alpha*dlog(r(i))-b2*r(i)) )
      ENDDO
      GOTO 999

************************************************************************
*  4 LOG NORMAL with rg and sigma given                                *
************************************************************************
   40 flogrg = dlog(par1)
      flogsi = dabs(dlog(par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      DO i=1,numr
         nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
      ENDDO
      GOTO 999

************************************************************************
*  5 LOG NORMAL with reff and veff given                               *
************************************************************************
   50 rg     = par1/(1.D0+par2)**2.5D0
      flogrg = dlog(rg)
      flogsi = dsqrt(dlog(1.D0+par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      DO i=1,numr
         nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
      ENDDO
      GOTO 999

************************************************************************
*  6 POWER LAW                                                         *
************************************************************************
   60 alpha = par1
      rmin  = par2
      rmax  = par3
      IF (dabs(alpha+1.D0) .lt. 1.d-10) THEN
         C = 1.D0/dlog(rmax/rmin)
      ELSE
         alpha1 = alpha-1.d0
         C = alpha1 * rmax**alpha1 / ((rmax/rmin)**alpha1-1.d0)
      ENDIF
      DO i=1,numr
         IF ((r(i) .lt. rmax) .and. (r(i) .gt. rmin)) THEN
            nwithr(i) = C*r(i)**(-alpha)
         ELSE
            nwithr(i) = 0.D0
         ENDIF
      ENDDO
      GOTO 999

************************************************************************
*  7 MODIFIED GAMMA with alpha, rc and gamma given                     *
************************************************************************
   70 alpha = par1
      rc    = par2
      gamma = par3
      b     = alpha / (gamma*rc**gamma)
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammln(aperg)
      DO i=1,numr
         nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
      ENDDO
      GOTO 999

************************************************************************
*  8 MODIFIED GAMMA with alpha, b and gamma given                      *
************************************************************************
   80 alpha = par1
      b     = par2
      gamma = par3
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammln(aperg)
      DO i=1,numr
         nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
      ENDDO
      GOTO 999

************************************************************************
  999 IF (numr.LE.1) RETURN
      RETURN
      END

      FUNCTION gammln(xarg)
************************************************************************
*  Return the value of the natural logarithm of the gamma function.    *
*  The argument xarg must be real and positive.                        *
*  This function is DOcumented in :                                    *
*                                                                      *
*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
*  page 157 (ISBN 0-521-30811)                                         *
*                                                                      *
*  When the argument xarg is between zero and one, the relation (6.1.4)*
*  on page 156 of the book by Press is used.                           *
*                                         V.L. Dolman April 18 1989    *
************************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION xarg
      PARAMETER( eps = 1.d-7, one = 1.D0, two = 2.D0, half = 0.5D0
     +         , fpf = 5.5D0 )
      DIMENSION cof(6)
      DATA cof,stp/76.18009173D0,-86.50532033D0, 24.01409822D0
     +   ,-1.231739516D0, 0.120858003D-2, -0.536382D-5, 2.50662827465D0/
      pi = 4.D0*datan(1.D0)

      IF (xarg .le. 0.D0) THEN
        WRITE(*,*) ' gammln: called with negative argument xarg = ',xarg
        STOP 'function gammln called with negative value'
      ENDIF
      IF (dabs(xarg-one) .lt. eps) THEN
         WRITE(*,*) ' gammln: argument too close to one for algorithm'
         STOP ' in function gammln argument too close to one'
      ENDIF
      IF (xarg .ge. one) THEN
         xx = xarg
      ELSE
         xx = xarg+two
      ENDIF
      x = xx-one
      tmp = x+fpf
      tmp = (x+half)*dlog(tmp)-tmp
      ser = one
      DO j=1,6
         x = x+one
         ser = ser+cof(j)/x
      ENDDO
      gtmp = tmp+dlog(stp*ser)
      IF (xarg .gt. one) THEN
         gammln = gtmp
      ELSE
         pix = pi*(one-xarg)
         gammln = dlog(pix/dsin(pix))-gtmp
      ENDIF

      RETURN
      END
