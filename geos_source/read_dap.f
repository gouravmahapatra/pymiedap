      SUBROUTINE READ_DAP(foufile,ngeos,alpha,thet0,theta,phi,beta,
     .                     rfou, Sv)

Cf2py intent(in) foufile,alpha,ngeos,thet0,theta,phi,beta,rfou
Cf2py intent(out) Sv

*----------------------------------------------------------------------------
*     Read a Fourier coefficients file and calculate the Stokes vector
*     for a given geometry.
*
*     It is assumed that the Stokes vector of the incoming sunlight 
*     is [1,0,0,0], with the flux measured perpendicular to the direction
*     of incidence equal to pi.
*
*     Author: Daphne M. Stam
*     Date: July 2013
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,nmat,nmugs,nfou,ngeos,ni,ki,j1,j2,m,k

      DOUBLE PRECISION fac,mu,mu0,eps,be,rf3,SvQ,SvU

      DOUBLE PRECISION xmu(nmuMAX),beta(ngeosMAX), alpha(ngeosMAX),
     .                 thet0(ngeosMAX),theta(ngeosMAX),phi(ngeosMAX),
     .                 rfou(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX),
     .                 rfm(nmatMAX),RM(nmatMAX),Bplus(4),
     .                 rfj(nmuMAX,nmatMAX),rf2(nmuMAX),rf(nmuMAX),
     .                 SvR(nmatMAX),Sv(nmatMAX,ngeosMAX)

      DOUBLE PRECISION pi,radfac
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

      DOUBLE PRECISION Sv0
      PARAMETER (Sv0=1.D0)

      CHARACTER foufile*50
 
*----------------------------------------------------------------------------
*     Open the output file:
*----------------------------------------------------------------------------
C     OPEN(unit=12,file='geos.out')
C     WRITE(12,800)
C     WRITE(12,801)

*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
      CALL rdfous(foufile,nfou,nmat,nmugs,xmu,rfou)
      eps= 1.D-10

*----------------------------------------------------------------------------
*     Initialize the Stokes vector of the whole planet:
*----------------------------------------------------------------------------
      DO k=1,nmatMAX
        DO j=1,ngeosMAX
         Sv(k,j)= 0.D0
         ENDDO
      ENDDO

*----------------------------------------------------------------------------
*     Loop over the geometries:
*----------------------------------------------------------------------------
      DO ni=1,ngeos
 
         mu0=DCOS(radfac*thet0(ni))
         mu= DCOS(radfac*theta(ni))
  
*----------------------------------------------------------------------------
*        Find the locations in array xmu where mu falls:
*----------------------------------------------------------------------------
         CALL bracks(mu,xmu,nmugs,nmuMAX,j1,j2)

*----------------------------------------------------------------------------
*        Initialize the 1st column of the reflection matrix:
*----------------------------------------------------------------------------
         DO ki=1,nmat
            RM(ki)= 0.D0
         ENDDO

*----------------------------------------------------------------------------
*        Loop over the Fourier coefficients:
*----------------------------------------------------------------------------
         DO m=0,nfou

            fac=1.D0
            IF (m.EQ.0) fac=0.5D0

            IF (j1.EQ.j2) THEN
               DO k=1,nmat
                  ki= (j1-1)*nmat+k 
                  DO i=1,nmugs
                     rf(i)= rfou(ki,i,m)
                  ENDDO
                  CALL spline(xmu,rf,nmugs,nmuMAX,rf2)
                  CALL splint(xmu,rf,rf2,nmugs,mu0,rf3,nmuMAX)
                  rfm(k)= rf3
               ENDDO
            ELSE
               DO j=1,nmugs
                  DO k=1,nmat
                     ki= (j-1)*nmat+k 
                     DO i=1,nmugs
                        rf(i)= rfou(ki,i,m)
                     ENDDO
                     CALL spline(xmu,rf,nmugs,nmuMAX,rf2)
                     CALL splint(xmu,rf,rf2,nmugs,mu0,rf3,nmuMAX)
                     rfj(j,k)= rf3
                  ENDDO
               ENDDO
               DO k=1,nmat
                  DO i=1,nmugs
                     rf(i)= rfj(i,k)
                  ENDDO
                  CALL spline(xmu,rf,nmugs,nmuMAX,rf2)
                  CALL splint(xmu,rf,rf2,nmugs,mu,rf3,nmuMAX)
                  rfm(k)= rf3
               ENDDO
            ENDIF

*----------------------------------------------------------------------------
*           Calculate the local reflection matrix (only the first column): 
*----------------------------------------------------------------------------
            Bplus(1)= DCOS(m*phi(ni)*radfac)
            Bplus(2)= DCOS(m*phi(ni)*radfac)
            Bplus(3)= DSIN(m*phi(ni)*radfac)
            Bplus(4)= DSIN(m*phi(ni)*radfac)

            DO k=1,nmat
               RM(k)= RM(k) + 2.D0*Bplus(k)*fac*rfm(k)
            ENDDO
         ENDDO

*----------------------------------------------------------------------------
*        Calculate the locally reflected Stokes vector with the matrix:
*----------------------------------------------------------------------------
         DO k=1,nmat
C           SvR(k)= 4.D0*mu0*RM(k)*apix/pi
            SvR(k)= RM(k)
            IF (DABS(RM(k)).LT.eps) SvR(k)=0.D0
         ENDDO

*----------------------------------------------------------------------------
*        Rotate Stokes elements Q and U to the actual reference plane:
*----------------------------------------------------------------------------
         IF (nmat.GT.1) THEN
            be= beta(ni)*radfac

            SvQ= DCOS(2.D0*be)*SvR(2) + DSIN(2.D0*be)*SvR(3) 
            SvU=-DSIN(2.D0*be)*SvR(2) + DCOS(2.D0*be)*SvR(3) 

            SvR(2)= SvQ
            SvR(3)= SvU
         ENDIF

*----------------------------------------------------------------------------
*        Add the new results to the Stokes vector elements:
*----------------------------------------------------------------------------
         DO k=1,nmat
            Sv(k,ni)= SvR(k)
         ENDDO
*----------------------------------------------------------------------------
*     Next geometry:
*----------------------------------------------------------------------------
      ENDDO

*----------------------------------------------------------------------------
*     Check for the accuracy:
*----------------------------------------------------------------------------
      DO k=1,nmat
         IF (DABS(Sv(k,ni)).LT.eps) Sv(k,ni)=0.D0
      ENDDO

*----------------------------------------------------------------------------
      WRITE(*,*) 'End of geometries program'
      END
