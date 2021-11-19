* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE fichid(m,x,nchi,nmax,nD,psi,chi,D)

************************************************************************
*  Calculate functions psi(x)  chi(x) and D(z) where z = mx.           *
*  On entry, the following should be supplied :                        *
*      m      : complex index of refraction                            *
*      x      : sizeparameter                                          *
*      nchi   : starting index for backwards recurrence of chi         *
*      nmax   : number of chi, psi and D that must be available        *
*      nD     : starting index for backwards recurrence of D           *
*  On exit, the desired functions are returned through psi, chi and D  *
************************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER nchi, nmax, nD
      DOUBLE PRECISION x
C     DOUBLE COMPLEX D,m,z,perz,zn1
C     possible error with type of D
      DOUBLE COMPLEX m,z,perz,zn1
      DOUBLE COMPLEX D(nd)
C     DIMENSION psi(0:nchi),chi(0:nmax+1),D(nd)
      DOUBLE PRECISION psi(0:nchi),chi(0:nmax+1)

Cf2py intent(in) m,x,nchi,nmax,nD
Cf2py intent(out) psi,chi,D

      z = m*x
      perz= 1.D0/z
      perx= 1.D0/x
      sinx= dsin(x)
      cosx= dcos(x)
************************************************************************
*  (mis-) use the array psi to calculate the functions rn(x)
*  De Rooij and van der Stap Eq. (A6)
************************************************************************
      DO n=nchi-1,0,-1
         psi(n)= 1.D0 / (dble(2*n+1)/x - psi(n+1))
      ENDDO

************************************************************************
*  Calculate functions D(z) by backward recurrence
*  De Rooij and van der Stap Eq. (A11)
************************************************************************
      D(nD)= dcmplx(0.D0,0.D0)
      DO n=nD - 1,1,-1
         zn1 = dble(n+1)*perz
         D(n)= zn1 - 1.D0/(D(n+1)+zn1)
      ENDDO

************************************************************************
*  De Rooij and van der Stap Eqs. (A3) and (A1)
*  De Rooij and van der Stap Eq. (A5) where psi(n) was equal to r(n)
*  and Eq. (A2)
************************************************************************
      psi(0)= sinx
      psi1  = psi(0)*perx - cosx
      IF (dabs(psi1).GT.1.d-4) THEN
         psi(1)= psi1
         DO n=2,nmax
            psi(n)= psi(n)*psi(n-1)
         ENDDO
      ELSE
         DO n=1,nmax
            psi(n)= psi(n)*psi(n-1)
         ENDDO
      ENDIF

************************************************************************
*  De Rooij and van der Stap Eqs. (A4) and (A2)
************************************************************************
      chi(0)= cosx
      chi(1)= chi(0)*perx + sinx
      DO n=1,nmax-1
         chi(n+1)= dble(2*n+1)*chi(n)*perx - chi(n-1)
      ENDDO

      RETURN
      END

