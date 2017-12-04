* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE expbmu(b,xmu,nmu,ebmu)

*----------------------------------------------------------------------
*     Calculate dexp(-b/mu) for all mu in xmu and put the 
*     result in array ebmu:
*----------------------------------------------------------------------
      IMPLICIT NONE 

      INCLUDE 'max_incl'

      INTEGER i,nmu

      DOUBLE PRECISION b,xmu(nmuMAX),ebmu(nmuMAX)

Cf2py intent(in) b, xmu, nmu
Cf2py intent(out) ebmu

*----------------------------------------------------------------------
      DO i=1,nmu
         IF (xmu(i).GT.(b/100.D0)) THEN
            ebmu(i)= DEXP(-b/xmu(i))
         ELSE
            ebmu(i)= 0.D0
         ENDIF
      ENDDO

*----------------------------------------------------------------------
      RETURN
      END
