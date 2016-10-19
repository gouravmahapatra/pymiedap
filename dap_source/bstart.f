      SUBROUTINE bstart(m,layer,coefs,ncoef,M0,xmumin,
     .                  a,b,b0,ndoubl)

*----------------------------------------------------------------------*
*  Calculate the optical thickness b0 at which doubling should be      *
*  started to obtain an error less than epsilon in the total layer     *
*  with thickness b. This is done for the m-th Fourier component.      *
*  The values of b0 and ndoubl returned are such that b = b0*2**ndoubl *
*  The algorithm is described in de Haan et al. (1987).                *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER m, layer
C      INTEGER M0(nlaysMAX), ncoef(nlaysMAX)
C      DOUBLE PRECISION coefs 
C      DOUBLE PRECISION xmumin, a, b, b0, ndoubl

      DIMENSION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX), 
     .          ncoef(nlaysMAX),M0(nlaysMAX)
      LOGICAL verbo

Cf2py intent(out) b0, ndoubl

      verbo = .false.
      IF (m.GT.M0(layer)) THEN
         PRINT *,' bstart: m = ',m,' larger than M0 = ',M0(layer)
         STOP 'in bstart m larger than M0'
      ENDIF
      IF (ncoef(layer).GT.M0(layer)) THEN
         PRINT *,' bstart: ncoef =',ncoef(layer)
     .           ,' larger than M0 =',M0(layer)
         STOP 'in bstart ncoef larger than M0'
      ENDIF

*----------------------------------------------------------------------*
*  Calculate the effective albedo am, Eq. (143)                        *
*----------------------------------------------------------------------*
      am = 0.D0
      DO k=m,M0(layer)
         amtry = a*dabs(coefs(1,1,k,layer))/(2.D0*k+1.D0)
         IF (am.LT.amtry) am = amtry
      ENDDO

*----------------------------------------------------------------------*
*  Right hand side of Eq. (142)                                        *
*----------------------------------------------------------------------*
      rhs = 4.D0*eps/(9.D0*b*am**3*dble(2*m+1))
      ndoubl = -1
      b0     = 2.D0*b
200      ndoubl = ndoubl+1
         b0     = 0.5D0*b0
         fb0mu  = b0/xmumin
         IF (fb0mu.GT.1.D0) fb0mu = 1.D0
         IF (b0.GE.(rhs/fb0mu)) GOTO 200
         IF (b0.GE.(1.D0/3.D0)) GOTO 200
      IF (verbo) PRINT *,' bstart: start optical thickness = ',b0
     .                                         ,' ndoubl=',ndoubl,m
      RETURN
      END
