      SUBROUTINE file2coefs(namein,ncoef,coefs)

****************************************************************
* PURPOSE: 
* In this program a file containing expansion 
* coefficients is read and the coefficients are
* returned
*
* DATE:
* 2017-10-05 L. Rossi
*
* AUTHOR:
* D. M. Stam, L. Rossi
****************************************************************
      IMPLICIT NONE

      INTEGER NDcoef, ncoef
      PARAMETER (NDcoef=4000)

      DOUBLE PRECISION coefs,SSA
      DIMENSION coefs(4,4,0:NDcoef)

      DOUBLE PRECISION pi
      PARAMETER (pi=3.1415926535898D0)

      CHARACTER*20 namein
      CHARACTER*23 title

Cf2py intent(in) namein
Cf2py intent(out) ncoef,coefs
****************************************************************
* 2 Read the expansion coefficients file:
****************************************************************
      OPEN (unit=10,file=namein,status='old')
         READ (10,'(A23)') title
         IF (title.NE.' EXPANSION COEFFICIENTS') THEN
         STOP 'wrong mie.sc file in readmiesc.f!'
      ENDIF
      CALL readsc(10,coefs,NDcoef,ncoef,SSA)
      CLOSE(10)

      END
