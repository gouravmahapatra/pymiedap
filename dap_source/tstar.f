      SUBROUTINE tstar(Ts,T,nmat,nmu)

*----------------------------------------------------------------------*
*  Calculate the transmission supermatrix Ts for illumination from     *
*  below from the normal transmission supermatrix T by symmetry :      *
*                         Ts = q3 T~ q3                                *
*  where T~ is the transpose of T, and q3 is defined above Eq. (96)    *
*  This symmetry is also valid for vertically inhomogeneous            *
*  atmospheres. It is described in Hovenier (1970) page I-6, however,  *
*  one should be aware of the difference between Hovenier's operator   *
*  T(mu,mu0,phi-phi0) and our supermatrix Tm(mu,mu0) !!                *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'max_incl'

C      INTEGER nmat, nmu
C      DOUBLE PRECISION T, Ts
      DIMENSION T(nsupMAX,nsupMAX),Ts(nsupMAX,nsupMAX)

Cf2py intent(out) Ts

      nsup= nmu*nmat

*----------------------------------------------------------------------*
*  Transpose T and put it in Ts                                        *
*----------------------------------------------------------------------*
      DO j=1,nsup
         DO i=1,nsup
            Ts(i,j) = T(j,i)
         ENDDO
      ENDDO

*----------------------------------------------------------------------*
*  Put a minus sign in every third row and in every third column       *
*----------------------------------------------------------------------*
      IF (nmat.GE.3) THEN
         DO i=3,nsup,nmat
            DO j=1,nsup
               Ts(i,j) = -Ts(i,j)
               Ts(j,i) = -Ts(j,i)
            ENDDO
         ENDDO
      ENDIF

*----------------------------------------------------------------------*
      RETURN
      END
