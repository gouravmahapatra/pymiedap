* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE renorm(Zmmin,Zmplus,nmu,nmat,xmu,smf)

*----------------------------------------------------------------------*
*  Renormalize the phase matrix by updating the diagonal elements      *
*  of the phase matrix.                                                *
*  This routine only makes sense if called for m=0 Fourier component.  *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER nmum,nmu, nmat
      INTEGER nmum,nmu

C      DOUBLE PRECISION Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX),
C     .          xmu(nmuMAX),smf(nmuMAX), w(nmuMAX)

      DIMENSION Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX),
     .          xmu(nmuMAX),smf(nmuMAX)
      DIMENSION w(nmuMAX)

      LOGICAL verbo
      verbo = .false.

Cf2py intent(in) nmu, nmat, xmu, smf
Cf2py intent(in,out) Zmmin, Zmplus

*----------------------------------------------------------------------*
*  Retrieve the weights from the supermatrix factors smf :             *
*----------------------------------------------------------------------*
      nmum= nmu-1

      DO i=1,nmum
         w(i) = 0.5D0*smf(i)**2/xmu(i)
      ENDDO
      fmax  =-1.0D0
      fmaxex=-1.0D0

      DO 400 j=1,nmu
*----------------------------------------------------------------------*
*         Calculate normalization integral r+t, which should be 2      *
*----------------------------------------------------------------------*
          r = 0.D0
          t = 0.D0
          jsup = (j-1)*nmat+1
          DO i=1,nmum
             isup = (i-1)*nmat+1
             r = r + Zmmin(isup,jsup)*w(i)
             t = t + Zmplus(isup,jsup)*w(i)
          ENDDO

*----------------------------------------------------------------------*
*         Update diagonal elements of phase matrix Fourier component    *
*         Only change the transmission part Zmplus, not Zmmin          *
*         For the extra mu-points the two points closest to the        *
*         forward direction are updated, each weighted by interpolation*
*----------------------------------------------------------------------*
          IF (j.LE.nmum) THEN
             fac= 1.D0 + (2.D0-r-t)/(Zmplus(jsup,jsup)*w(j))
             IF (fac.GT.fmax) fmax= fac
             DO k=1,nmat
                Zmplus(jsup+k-1,jsup+k-1)=fac*Zmplus(jsup+k-1,jsup+k-1)
             ENDDO
          ELSE
             CALL brack(xmu(j),xmu,nmum,nmuMAX,i1,i2)
             IF (i1.NE.i2) THEN
                relw1= (xmu(i2)-xmu(j))/(xmu(i2)-xmu(i1))
                relw2= (xmu(j)-xmu(i1))/(xmu(i2)-xmu(i1))
             ELSE
                relw1 = 1.0d0
                relw2 = 0.0d0
             ENDIF
             isup1 = (i1-1)*nmat+1
             isup2 = (i2-1)*nmat+1
             fac1=1.0d0 + relw1*(2.0d0-r-t)/(Zmplus(isup1,jsup)*w(i1))
             fac2=1.0d0 + relw2*(2.0d0-r-t)/(Zmplus(isup2,jsup)*w(i2))
             IF (fac1.GT.fmaxex) fmaxex = fac1
             DO k=1,nmat
              Zmplus(isup1+k-1,jsup+k-1)=fac1*Zmplus(isup1+k-1,jsup+k-1)
              Zmplus(jsup+k-1,isup1+k-1)=Zmplus(isup1+k-1,jsup+k-1)
              Zmplus(isup2+k-1,jsup+k-1)=fac2*Zmplus(isup2+k-1,jsup+k-1)
              Zmplus(jsup+k-1,isup2+k-1)=Zmplus(isup2+k-1,jsup+k-1)
             ENDDO
          ENDIF
400   CONTINUE

*---------------------------------------------------------------------
*     Print the maximum renormalization factor (verbose only):
*---------------------------------------------------------------------
      IF (verbo) THEN
        PRINT *,' renorm: max renorm. factor ',fmax,' Gauss pts'
        PRINT *,'         ',fmaxex,' (extra points)'
      ENDIF

*---------------------------------------------------------------------
      RETURN
      END
