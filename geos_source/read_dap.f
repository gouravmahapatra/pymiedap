* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE READ_DAP(foufile,derivs,rfou,xmu,thet0,theta,phi,beta,
     .                  filetype,nfou,nmat,nmugs,ngeos,Sv)

Cf2py intent(in) foufile,rfou,derivs,thet0,theta,phi,beta,filetype
Cf2py intent(in) nfou,nmat,nmugs,ngeos
Cf2py intent(out) Sv

*----------------------------------------------------------------------------
*     Read a Fourier coefficients file and calculate the Stokes vector
*     for a given geometry.
*
*     It is assumed that the Stokes vector of the incoming sunlight 
*     is [1,0,0,0], with the flux measured perpendicular to the direction
*     of incidence equal to pi.
*
*     Author: Ashwyn Groot, based on read_dap.f routine of Daphne M. Stam
*     Date: November 2018
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,nmat,nmugs,nfou,ngeos,m,k,
     .          ier,filetype

      INTEGER, DIMENSION(:), ALLOCATABLE :: ki
      REAL*8, DIMENSION(:), ALLOCATABLE :: xd_1d,yd_1d
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rfm,RM,zd
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: WK

      REAL*8, DIMENSION(ngeos) :: mu,mu0,be,SvQ,SvU,xi,yi,zi,beta,thet0,
     .                  theta,phi

      REAL*8 xmu(nmugs)
      REAL*8 Bplus(4,ngeos)
      REAL*8 Sv(nmat,ngeos)

      REAL*8 derivs(0:nfou-1,nmat,3,nmugs,nmugs)
      REAL*8 rfou(nmat*nmugs,nmugs,0:nfou-1)
      REAL*8 pi,radfac,fac,eps,eps1,Sv0
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

      PARAMETER (Sv0=1.D0)

      CHARACTER(len=200) :: foufile
*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
C      IF (filetype.EQ.1) THEN
C         CALL rdfousderiv(foufile,nfou,nmat,nmugs,rfou,derivs,xmu)
C      ELSEIF (filetype.EQ.2) THEN
C         CALL rdfous(foufile,nfou,nmat,nmugs,rfou,xmu)
C      ENDIF

*----------------------------------------------------------------------------
*     Allocate new variables:
*----------------------------------------------------------------------------
      ALLOCATE(xd_1d(nmugs),yd_1d(nmugs),ki(nmugs),rfm(nmat,ngeos),
     .              RM(nmat,ngeos),zd(nmugs,nmugs),
     .              WK(3,nmugs,nmugs))

*----------------------------------------------------------------------------
*     Initialize the Stokes vector of the whole planet:
*----------------------------------------------------------------------------
      eps= 1.D-10
      eps1= 1.D-100

*----------------------------------------------------------------------------
*     Initialize mu values for all ngeos:
*----------------------------------------------------------------------------
      mu0=DCOS(radfac*thet0(:ngeos))
      mu= DCOS(radfac*theta(:ngeos))
*----------------------------------------------------------------------------
*        Initialize the 1st column of the reflection matrix:
*----------------------------------------------------------------------------
      RM=0.D0
      Sv=0.D0

      xi=mu
      yi=mu0
      xd_1d=xmu
      yd_1d=xmu
      If (filetype.EQ.1) THEN
         WHERE (abs(derivs).LT.eps1)
            derivs=0.D0
         ENDWHERE
         IF (sum(derivs).LT.eps1) THEN
            STOP 'USE OTHER FILETYPE: NO DERIVATIVES IN FILE!'
         ENDIF
      ENDIF
*----------------------------------------------------------------------------
*        Loop over the Fourier coefficients:
*----------------------------------------------------------------------------
!  NOTE: be very carefull with the meaning of nfou in this routine and
!       read_dapascii.f. nfou is used here, in the hdf5 rdfous and hdf5
!       newfou routines as the dimension of all fourier coefficients, hence
!       why we use nfou-1 here.
      DO m=0,nfou-1
C         print *, 'Progress: ',100.*(m+1)/nfou,'%'
         fac=1.D0
         IF (m.EQ.0) fac=0.5D0

         DO k=1,nmat
            ki=((/(i,i=1,nmugs,1)/)-1)*nmat+k
            zd=rfou(ki,:nmugs,m)
            WHERE (abs(zd).LT.eps1)
               zd=0.D0
            ENDWHERE
            IF (filetype.EQ.1) THEN
               WK=derivs(m,k,:,:,:)
            ENDIF
            IF (nmugs.GE.4) THEN
               call RGBI3P(1,xd_1d,yd_1d,zd,xi,yi,zi,ier,filetype,
     .                  WK,nmugs,nmugs,ngeos)
            ENDIF
            rfm(k,:)= zi
         ENDDO

*----------------------------------------------------------------------------
*           Calculate the local reflection matrix (only the first column): 
*----------------------------------------------------------------------------
         Bplus(1,:)= DCOS(m*phi(:ngeos)*radfac)
         Bplus(2,:)= DCOS(m*phi(:ngeos)*radfac)
         Bplus(3,:)= DSIN(m*phi(:ngeos)*radfac)
         Bplus(4,:)= DSIN(m*phi(:ngeos)*radfac)
         DO k=1,nmat
            RM(k,:)= RM(k,:) + 2.D0*Bplus(k,:)*fac*rfm(k,:)
         ENDDO
      ENDDO
      DEALLOCATE(xd_1d,yd_1d,ki,rfm,zd,WK)
*----------------------------------------------------------------------------
*        Calculate the locally reflected Stokes vector with the matrix:
*----------------------------------------------------------------------------
      DO k=1,nmat
C         SvR(k)= 4.D0*mu0*RM(k)*apix/pi
C         SvR(k,:)= RM(k,:)
         Sv(k,:ngeos)= RM(k,:)
         WHERE (DABS(Sv(k,:)).LT.eps)
            Sv(k,:)=0.D0
         ENDWHERE
      ENDDO
      DEALLOCATE(RM)
*----------------------------------------------------------------------------
*        Rotate Stokes elements Q and U to the actual reference plane:
*----------------------------------------------------------------------------
      IF (nmat.GT.1) THEN
         be= beta(:ngeos)*radfac

         SvQ= DCOS(2.D0*be)*Sv(2,:ngeos) + DSIN(2.D0*be)*Sv(3,:ngeos)
         SvU=-DSIN(2.D0*be)*Sv(2,:ngeos) + DCOS(2.D0*be)*Sv(3,:ngeos)

         Sv(2,:ngeos)= SvQ
         Sv(3,:ngeos)= SvU
      ENDIF
*----------------------------------------------------------------------------
*     Check for the accuracy:
*----------------------------------------------------------------------------
      DO k=1,nmat
         WHERE (DABS(Sv(k,:)).LT.eps)
            Sv(k,:)=0.D0
         ENDWHERE
      ENDDO

*----------------------------------------------------------------------------
      WRITE(*,*) 'End of geometries program'
      END
