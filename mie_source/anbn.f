      SUBROUTINE anbn(m,x,nmax,psi,chi,D,an,bn)

************************************************************************
*  Calculate the Mie coefficients an and bn.                           *
************************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER ( NDn=30000 )
      INTEGER nmax
      DOUBLE PRECISION x, psi(0:NDn), chi(0:NDn)
      DOUBLE COMPLEX m, zn, znm1, save, perm
      DOUBLE COMPLEX an(NDn), bn(NDn), D(NDn)
C     DIMENSION psi(0:NDn), chi(0:NDn)

Cf2py intent(in) m,x,nmax,psi,chi,d
Cf2py intent(out) an, bn

      perm= 1.D0/m
      perx= 1.D0/x
      xn  = 0.D0

      DO n=1,nmax
         zn   = dcmplx(psi(n),  chi(n))
         znm1 = dcmplx(psi(n-1),chi(n-1))
         xn   = dble(n)*perx
         save = D(n)*perm+xn
         an(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
         save = m*D(n)+xn
         bn(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
      ENDDO

      RETURN
      END
