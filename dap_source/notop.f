* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE notop(Rmbot,Tmbot,ebtop,ebbot,nmat,nmu,nsup)

*----------------------------------------------------------------------
*     Use the adding equations when there is no scattering in 
*     the top layer:
*  Edited by: Ashwyn Groot
*  Date: November 2018
*  Introduced matrix operations with f95<
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,k,nmu,nmat,nsup

      REAL*8, DIMENSION(nsup,nsup) :: Rmbot,Tmbot,Etop,Ebot,Rmbotsub

      REAL*8, DIMENSION(nmu) :: ebtop, ebbot

Cf2py intent(in,out) Rmbot,Tmbot,ebtop,ebbot

*----------------------------------------------------------------------
      Etop=0.D0
      Ebot=0.D0
      DO i=1,nmu
         DO k=1,nmat
            Etop((i-1)*nmat+k,(i-1)*nmat+k) = ebtop(i)
            Ebot((i-1)*nmat+k,(i-1)*nmat+k) = ebbot(i)
         ENDDO
      ENDDO 

*----------------------------------------------------------------------
*     Reflection  R= E'R"E' 
*----------------------------------------------------------------------

      Rmbotsub=MATMUL(Rmbot,Etop)

      Rmbot=MATMUL(Etop,Rmbotsub)
*----------------------------------------------------------------------
*     Transmission T= T"E'  
*----------------------------------------------------------------------

      Tmbot=MATMUL(Tmbot,Etop)

*----------------------------------------------------------------------
*     Reflection star R* = R"* is trivial !                               
*----------------------------------------------------------------------
*     Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu)        
*----------------------------------------------------------------------
      ebbot=ebbot*ebtop

*----------------------------------------------------------------------
      RETURN
      END
