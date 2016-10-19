      SUBROUTINE endfou(m,M1,M0,nlays,nextm)

*----------------------------------------------------------------------*
*  Decide if the next Fourier term (m+1) is needed, if so return       *
*  nextm = .true. else return nextm = .false.                          *
*  Because first order is treated separately without Fourier expansion *
*  it is enough to stop when only first order is needed.               *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER m, nlays, maxM
C      INTEGER M1, M0
C      DOUBLE PRECISION M1, M0
      DIMENSION M1(nlaysMAX),M0(nlaysMAX)

      LOGICAL nextm,verbo
      verbo = .false.

      nextm = .true.

Cf2py intent(out) nextm

*----------------------------------------------------------------------*
*  First test whether the user supplied bound nfou has been reached *
*----------------------------------------------------------------------*
      IF (m.GE.nfouMAX) THEN
         nextm = .false.
         IF (verbo) PRINT *,' endfou: stop Fourier series after m = '
     .                      ,m,' (nfouMAX reached)'
         GOTO 999
      ENDIF

*----------------------------------------------------------------------*
*  Test if M1 has been reached (see de Haan et al. (1987) section 7.2) *
*  so that higher Fourier terms will only have single scattering.      *
*  Of course we should take the maximum M1 over all layers.            *
*  If first order scattering is not excepted (except = .false.) then   *
*  we must sum the Fourier series all the way to M0.                   *
*----------------------------------------------------------------------*
      maxM=0
      DO layer=1,nlays
         IF (maxM.LT.M0(layer)) maxM = M0(layer)
      ENDDO

      IF (m.GE.maxM) THEN
         nextm= .false.
         IF (verbo) PRINT *,' endfou: stop Fourier series after m = '
     .                      ,m,' (M1 reached)'
         GOTO 999
      ENDIF

*----------------------------------------------------------------------*
999   RETURN
      END
