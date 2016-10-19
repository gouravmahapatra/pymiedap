      SUBROUTINE bracks(el,array,n,ND,i1,i2)

Cf2py intent(in) el, array, n, ND
Cf2py intent(out) i1, i2

*-------------------------------------------------------------------------------
* PURPOSE:
* Find a bracket around "el" in the first "n" elements of
* "array".
* 
* INPUT:
* 	el	: element to be bracketed
* 	array   : array in which bracket must be found
* 	n       : number of elements in array to be considered
* 	ND      : dimension of array
*
* OUTPUT:
*	i1	: index of left bracket element
* 	i2	: index of right bracket element
*
*
* COMMENTS:
* If "el" is outside the range of "array", "i1" and "i2" are 
* set to the appropriate extreme index value.    
* If "el" is equal to one of the elements of "array", "i2" is
* set to the appropriate index value.
*
* DATE: 
* ...
*
* AUTHOR:
* V. Dolman
*
* MODIFICATIONS:
* P. Stammes, 18 Jan. 1993:
*    adapted
* W. Wauben, 25 May 1994:
*    bug repaired in the decreasing array option
* D. Stam, July 2013:
*    removed the option for decreasing array
*-------------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION array(ND)
      DOUBLE PRECISION el

      INTEGER n, ND, i1, i2

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-8)

*-------------------------------------------------------------------------------
      IF (el.LE.(array(1)+eps)) THEN
         i1= 1
         i2= 1
      ELSEIF (el.GT.(array(n)-eps)) THEN
          i1= n
          i2= n
      ELSE
         DO i=2,n
        IF ((el.LE.(array(i)+eps)).AND.(el.GT.(array(i-1)+eps))) THEN
           i2= i
           i1= i-1
            ENDIF
         ENDDO
      ENDIF

*-------------------------------------------------------------------------------
      RETURN
      END
