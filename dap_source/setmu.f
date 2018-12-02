* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE setmu(nmug,nmu,xmu,smf)

**********************************************************************
*  Initialize the mu-values and the supermatrixfactors.         
*
*  On entry :                                          
*      nmug     : number of desired Gauss points                     
*      nmu      : number of desired Gauss points + 1
*
*  On exit :                                                         
*      smf      : array containing the supermatrix factors,        
*                 dsqrt(2*w*mu) for Gauss points, 1 for extra points.
*      xmu      : array containing the gaussian abscissas
**********************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER nmug,nmu

      REAL*8, DIMENSION(nmug+1) :: xmu, smf, wmu

Cf2py intent(in) nmug
Cf2py intent(out) nmu,xmu,smf

      nmu= nmug+1

*---------------------------------------------------------------------
*     Get the nmug Gaussian points and weights:
*---------------------------------------------------------------------
      CALL gauleg(nmu,nmug,0.D0,1.D0,xmu,wmu)

*----------------------------------------------------------------------
*     Add the extra mu value for the nadir direction:
*----------------------------------------------------------------------
      xmu(nmu)= 1.D0
      wmu(nmu)= 0.5D0

*----------------------------------------------------------------------
*     Change the Gaussian weights into supermatrix weights:
*----------------------------------------------------------------------
      smf=DSQRT(2.D0*wmu*xmu)

**********************************************************************
      RETURN
      END
