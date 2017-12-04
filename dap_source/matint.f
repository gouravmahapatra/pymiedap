* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE matint(matrix, phi, nmat, nmu, nphi, nm, matexp) 
Cf2py intent(in) matrix, phi, nmat, nmu, nphi, nm
Cf2py intent(out) matexp

**********************************************************************
*            Expansion of a matrix in Fourier terms with 
*          integration method (see de Haan 1987, section 4)
**********************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nmat, nmu, nphi, nm
      INTEGER i,j,k,l,p,m

      INTEGER nphiMAX, nmMAX
      DOUBLE PRECISION pi
      PARAMETER (nphiMAX=1000,
     .           nmMAX=200,pi=3.14159d0)
      DOUBLE PRECISION matrix(nmat,nmat,nmu,nmu,nphi)
      DOUBLE PRECISION matexp(nmat,nmat,nmu,nmu,nmMAX)
      DOUBLE PRECISION BL(nmat,nmat)
      DOUBLE PRECISION Bplus(nmat, nmat, nphi,nmMAX)
      DOUBLE PRECISION Bmin(nmat, nmat, nphi,nmMAX)
      DOUBLE PRECISION Btot(nmat, nmat, nphi,nmMAX)
      DOUBLE PRECISION phi(nphi)


* --  INITIALISATION
      
      DO i=1,nmatMAX
         DO j=1,nmatMAX
            DO p=1,nphiMAX
                DO m=1,nmMAX
                    Bplus(i,j,p,m)= 0.D0
                    Bmin(i,j,p,m)= 0.D0
                ENDDO
            ENDDO
         ENDDO
      ENDDO

*--     Generating the B+ and B- matrices

        DO p=1,nphi
            DO m=0,nm
                !Bplus(1,1,p,m+1)= DCOS(m*phi(p))
                !Bplus(2,2,p,m+1)= DCOS(m*phi(p))
                !Bplus(3,3,p,m+1)= DSIN(m*phi(p))
                !Bplus(4,4,p,m+1)= DSIN(m*phi(p))

                !Bmin(1,1,p,m+1)= -DSIN(m*phi(p))
                !Bmin(2,2,p,m+1)= -DSIN(m*phi(p))
                !Bmin(3,3,p,m+1)= DCOS(m*phi(p))
                !Bmin(4,4,p,m+1)= DCOS(m*phi(p))

                Btot(1,1,p,m+1)= DCOS(m*phi(p)) - DSIN(m*phi(p))
                Btot(2,2,p,m+1)= DCOS(m*phi(p)) - DSIN(m*phi(p))
                Btot(3,3,p,m+1)= DSIN(m*phi(p)) + DCOS(m*phi(p))
                Btot(4,4,p,m+1)= DSIN(m*phi(p)) + DCOS(m*phi(p))
            ENDDO
        ENDDO

*-  Product of matrices
      
      DO k=1,nmu
        DO l=1,nmu
            DO p=1,nphi
                DO m=1,nm
                !Btemp(i,j) = Btot(i,j,p,m)
                !Mtemp(i,j) = matrix(i,j,k,l,m)
                write (*,*) 'BLOUP'
                BL(1:nmatMAX,1:nmatMAX) = 
     .          matmul( Btot(1:nmatMAX,1:nmatMAX,p,m), 
     .                 matrix(1:nmatMAX,1:nmatMAX,k,l,p) )

                        matexp(1:nmat,1:nmat,k,l,m) = 
     .                   matexp(1:nmat,1:nmat,k,l,m)+
     .                  BL(1:nmat,1:nmat)
     .                  *(phi(p)-phi(p-1))/(2*pi)
                ENDDO
            ENDDO
        ENDDO
      ENDDO
      write (*,*) 'BLANG!'

*--     Integration

      !DO i=1,nmat
        !DO j=1,nmat
!           DO k=1,nmu
!               DO l=1,nmu
!                   DO m=1,nmMAX
!                       DO p=2,nphi
!                       write (*,*) i,j,k,l,m,p
!                        matexp(i,j,k,l,m) = matexp(i,j,k,l,m)+
!     .                  BL(i,j,k,l,p-1)*(phi(p)-phi(p-1))/(2*pi)
!                       matexp(1:nmat,1:nmat,k,l,m) = 
!    .                   matexp(1:nmat,1:nmat,k,l,m)+
!    .                  BL(1:nmat,1:nmat,k,l,p-1)
!    .                  *(phi(p)-phi(p-1))/(2*pi)
!                       ENDDO
!                   ENDDO
!               ENDDO
!           !ENDDO
!       !ENDDO
!     ENDDO

      RETURN
      END SUBROUTINE MATINT
