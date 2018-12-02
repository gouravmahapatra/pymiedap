* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE nobot(Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot,
     .                 ebtop,ebbot,nmat,nmu,nsup)

*----------------------------------------------------------------------
*     Use the adding equations when there is no scattering in the 
*     bottom layer:
*  Edited by: Ashwyn Groot
*  Date: November 2018
*  Introduced matrix operations with f95<
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,k,nmat,nmu,nsup

      REAL*8, DIMENSION(nsup,nsup) :: Rmtop,Tmtop,Rmbot,
     .          Tmbot,Rmsbot,Etop, Ebot

      REAL*8, DIMENSION(nmu) :: ebtop, ebbot

Cf2py intent(in,out) Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot,ebtop,ebbot

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
*     Reflection  R= R' 
*----------------------------------------------------------------------
      Rmbot=Rmtop

*----------------------------------------------------------------------
*     Transmission T= E"T'                                    
*----------------------------------------------------------------------
      Tmbot=MATMUL(Ebot,Tmtop)

*----------------------------------------------------------------------
*     Reflection star R*= E"R'*E"                                      
*----------------------------------------------------------------------
      CALL star(Rmsbot,Rmtop,nmat,nmu,nsup)
      Rmsbot=MATMUL(Rmsbot,Ebot)
      Rmsbot=MATMUL(Ebot,Rmsbot)
*----------------------------------------------------------------------
*     Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu) 
*----------------------------------------------------------------------
      ebbot=ebbot*ebtop

*----------------------------------------------------------------------
      RETURN
      END
