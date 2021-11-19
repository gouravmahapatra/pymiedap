* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE brack(el,array,n,ND,i1,i2)

****************************************************************
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
*
* AUTHOR:
* V. Dolman
*
*
* MODIFICATIONS:
* P. Stammes, 18 Jan. 1993:
*    adapted
* W. Wauben, 25 May 1994:
*    bug repaired in the decreasing array option
****************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

C      INTEGER n, ND 
C      INTEGER i1, i2
C      DOUBLE PRECISION el, array
      DIMENSION  array(ND)

Cf2py intent(out) i1, i2

*---------------------------------------------------------------------
      IF (array(1).LE.array(n)) THEN

*---------------------------------------------------------------------
*     The array is increasing:
*---------------------------------------------------------------------
         IF (el.LE.array(1)) THEN
            i1= 1
            i2= 1
         ELSEIF (el.GT.array(n)) THEN
            i1= n
            i2= n
         ELSE
            DO i=2,n
               IF ((el.LE.array(i)).AND.(el.GT.array(i-1))) THEN
                  i2= i
                  i1= i-1
               ENDIF
            ENDDO
         ENDIF
      ELSE

*---------------------------------------------------------------------
*     The array is decreasing:
*---------------------------------------------------------------------
         IF (el.GT.array(1)) THEN
            i1= 1
            i2= 1
         ELSEIF (el.LE.array(n)) THEN
            i1= n
            i2= n
         ELSE
            DO i=n,2,-1
               IF ((el.GE.array(i)).AND.(el.LT.array(i-1))) THEN
                  i2= i
                  i1= i-1
               ENDIF
            ENDDO
         ENDIF

      ENDIF

*---------------------------------------------------------------------
      RETURN
      END

