* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE layer0(surfmat,smf,nmat,ebbot,Rmbot,Tmbot,
     .              Rmsbot,nmu,nsup)
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nsup,nmu,nmat,z
      INTEGER, DIMENSION(nmu) :: i, j, ibase, jbase

      REAL*8, DIMENSION(nmu) :: ebbot, smf, w

      REAL*8, DIMENSION(nsup,nsup) :: Rmbot, Tmbot, Rmsbot, surfmat

Cf2py intent(in) surfmat,smf,nmu,nmat,nsup
Cf2py intent(out) ebbot
Cf2py intent(in,out) Rmbot,Tmbot,Rmsbot

*-----------------------------------------------------------------------
*     Fill the reflection arrays with Lambertian reflection, and 
*     the transmission arrays with zero's:
*  Edited by: Ashwyn Groot
*  Date: November 2018
*  Introduced matrix operations with f95<
*-----------------------------------------------------------------------

      Rmbot=0.D0
      Tmbot=0.D0
      Rmsbot=0.D0

      i = (/(z, z=1,nmu, 1)/)
      j = i
      ibase = (i-1)*nmat
      jbase = (j-1)*nmat
      w = smf(i)*smf(j)
      Rmbot(ibase+1,jbase+1)= surfmat(ibase+1,jbase+1)

      CALL star(Rmsbot,Rmbot,nmat,nmu,nsup)

*-----------------------------------------------------------------------
*     Fill array ebbot (the direct transmission through the surface),
*     which equals zero, but is needed anyway:
*-----------------------------------------------------------------------
      ebbot=0.D0
*-----------------------------------------------------------------------
      RETURN
      END
