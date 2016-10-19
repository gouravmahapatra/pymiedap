      SUBROUTINE addsm(C,A,B,nmat,nmu)

*-----------------------------------------------------------------------
*     Calculate the supermatrix sum C = A+B   
*
*     On entry:                            
*     A, B    : supermatrices to be added          
*     nmat    : number of elements of the Stokes vector taken into 
*               account (4 = full polarization, 3 = 3x3 approximation,
*               2 = illegal, 1 = scalar)                    
*     nmu    : total number of distinct mu points        
*
*     On exit:                                         
*     C       : supermatrix sum of A and B     
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER nmat, nmu
      DOUBLE PRECISION A(nsupMAX,nsupMAX),
     .                 B(nsupMAX,nsupMAX), 
     .                 C(nsupMAX,nsupMAX)

Cf2py intent(out) C

      nsup= nmu*nmat
      DO j=1,nsup
         DO i=1,nsup
            C(i,j)= A(i,j)+B(i,j)
         ENDDO 
      ENDDO

*-----------------------------------------------------------------------
      RETURN
      END
