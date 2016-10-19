      SUBROUTINE baerosols(nlays,nmat,baer,coefin,ncoefin,
     .                     basca,baabs,coefs,coefsa,ncoefsa)

************************************************************************
* DATE: July 2004
* UPDATE: May 2012
*
* AUTHOR: D. M. Stam
************************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,k,l,m,nlays,nmat
      INTEGER ncoefin(nlaysMAX)

      INTEGER ncoefsa(nlaysMAX)

      INTEGER ico
      PARAMETER (ico=57)

      DOUBLE PRECISION alb

      DOUBLE PRECISION baer(nlaysMAX),
     .                 basca(nlaysMAX),baabs(nlaysMAX),
     .                 coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 coefin(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 coefsa(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX)

Cf2py intent(in) nlays, nmat, baer, coefin, ncoefin
Cf2py intent(out) baabs, basca, coefs, coefsa, ncoefsa

*-----------------------------------------------------------------------
*     Initialisation:
*-----------------------------------------------------------------------
      DO i=1,nlaysMAX
         basca(i)= 0.D0
         baabs(i)= 0.D0
         ncoefsa(i)= 0
         DO j=1,nmatMAX
            DO k=1,nmatMAX
               DO l=0,ncoefsMAX
                  coefs(j,k,l,i)= 0.D0
                  coefsa(j,k,l,i)= 0.D0
               ENDDO
            ENDDO
         ENDDO
      ENDDO

*-----------------------------------------------------------------------
*     Loop over the layers:
*-----------------------------------------------------------------------
      DO i=1,nlays
         IF (baer(i).GT.1.D-5) THEN

*-----------------------------------------------------------------------
*           Read the expansion coefficients file:
*-----------------------------------------------------------------------
            ncoefsa(i)= ncoefin(i)

*-----------------------------------------------------------------------
*           Combine the data for the layer:
*-----------------------------------------------------------------------
            basca(i)= alb*baer(i)
            baabs(i)= (1.D0-alb)*baer(i)

            DO m=0,ncoefsa(i)
               DO j=1,nmat
                  DO k=1,nmat
                     coefsa(j,k,m,i)= coefin(j,k,m,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      ENDDO

************************************************************************
      RETURN
999   WRITE(*,*) 'baerosols: error reading expansion coefs. files'
      STOP
      END
