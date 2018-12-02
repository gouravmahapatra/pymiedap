* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE expbmu(b,xmu,nmu,ebmu)

*----------------------------------------------------------------------
*     Calculate dexp(-b/mu) for all mu in xmu and put the 
*     result in array ebmu:
*  Edited by: Ashwyn Groot                                             *
*  Date: November 2018                                                 *
*  Introduced matrix operations with f95<                              *
*----------------------------------------------------------------------
      IMPLICIT NONE 

      INCLUDE 'max_incl'

      INTEGER nmu

      REAL*8 :: b

      REAL*8, DIMENSION(nmu) :: xmu, ebmu !rank 1

Cf2py intent(in) b, xmu, nmu
Cf2py intent(out) ebmu

*----------------------------------------------------------------------
      WHERE (xmu.GT.(b/100.D0))
         ebmu= DEXP(-b/xmu)
      ELSEWHERE
         ebmu=0.D0
      ENDWHERE
*----------------------------------------------------------------------
      RETURN
      END
