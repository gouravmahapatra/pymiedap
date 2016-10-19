      SUBROUTINE layer0(surfmat,smf,nmu,nmat,ebbot,Rmbot,Tmbot,Rmsbot)

*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,nsup,nmu,nmat,ibase,jbase

      DOUBLE PRECISION w

      DOUBLE PRECISION ebbot(nmuMAX),smf(nmuMAX),
     .                 Rmbot(nsupMAX,nsupMAX),Tmbot(nsupMAX,nsupMAX),
     .                 Rmsbot(nsupMAX,nsupMAX),
     .                 surfmat(nsupMAX,nsupMAX)

Cf2py intent(in) surfmat,smf,nmu,nmat
Cf2py intent(out) ebbot
Cf2py intent(in,out) Rmbot,Tmbot,Rmsbot

*-----------------------------------------------------------------------
      nsup= nmu*nmat

*-----------------------------------------------------------------------
*     Fill the reflection arrays with Lambertian reflection, and 
*     the transmission arrays with zero's:
*-----------------------------------------------------------------------
      DO i=1,nsup
         DO j=1,nsup
            Rmbot(i,j)= 0.D0
            Tmbot(i,j)= 0.D0
            Rmsbot(i,j)= 0.D0
         ENDDO
      ENDDO

      DO i=1,nmu
         ibase= (i-1)*nmat
         DO j=1,nmu
            jbase= (j-1)*nmat
            w= smf(i)*smf(j)
            Rmbot(ibase+1,jbase+1)= surfmat(ibase+1,jbase+1)
         ENDDO
      ENDDO

      CALL star(Rmsbot,Rmbot,nmat,nmu)

*-----------------------------------------------------------------------
*     Fill array ebbot (the direct transmission through the surface),
*     which equals zero, but is needed anyway:
*-----------------------------------------------------------------------
      DO i=1,nmu
         ebbot(i)= 0.D0
      ENDDO

*-----------------------------------------------------------------------
      RETURN
      END
