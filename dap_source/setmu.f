      SUBROUTINE setmu(nmat,nmug,iunfou,nmu,xmu,smf)

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

      INTEGER iunfou,nmug,nmat,nmu

      DOUBLE PRECISION xmu(nmuMAX),smf(nmuMAX),wmu(nmuMAX)
C      DIMENSION xmu(nmuMAX),smf(nmuMAX),wmu(nmuMAX)

Cf2py intent(in) nmat,nmug,iunfou
Cf2py intent(out) nmu,xmu,smf

*---------------------------------------------------------------------
*     Get the nmug Gaussian points and weights:
*---------------------------------------------------------------------
      CALL gauleg(nmuMAX,nmug,0.D0,1.D0,xmu,wmu)

*----------------------------------------------------------------------
*     Add the extra mu value for the nadir direction:
*----------------------------------------------------------------------
      nmu= nmug+1
      xmu(nmu)= 1.D0
      smf(nmu)= 1.D0

*----------------------------------------------------------------------
*     Change the Gaussian weights into supermatrix weights:
*----------------------------------------------------------------------
      DO i=1,nmug
         smf(i) = DSQRT(2.D0*wmu(i)*xmu(i))
      ENDDO

*----------------------------------------------------------------------
*     Write the angles to the Fourier output file,
*     and also the accuracy eps:
*----------------------------------------------------------------------
*      WRITE(iunfou,'(E12.6)') eps
      WRITE(iunfou,'(I3)') nmat
      WRITE(iunfou,'(I3)') nmu
      DO i=1,nmu
         WRITE(iunfou,'(E16.8,3X,E16.8)') xmu(i),wmu(i)
      ENDDO

**********************************************************************
      RETURN
      END
