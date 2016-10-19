      SUBROUTINE trace(A,trA,nmat,nmum)

*---------------------------------------------------------------------
*  Calculate the truncated supermatrix trace(A).                    
*  The sum runs over the integration points only, see remark below
*  Eq. (124) of de Haan et al. (1987)                          
*---------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,ng,nmum,nmat

      DOUBLE PRECISION trA,A(nsupMAX,nsupMAX)

Cf2py intent(out) trA

*---------------------------------------------------------------------
      ng = nmum*nmat

      trA = 0.D0
      DO i=1,ng
         trA = trA + A(i,i)
      ENDDO

*---------------------------------------------------------------------
      RETURN
      END
