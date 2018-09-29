* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE top2bot(nmat,nmu,ebtop,ebbot,Rmtop,Tmtop,
     .                   Rmbot,Tmbot,Rmsbot)

*----------------------------------------------------------------------
*     Copy data pertaining to the top layer to the bottom layer:
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,nmat,nmu

      DOUBLE PRECISION Rmtop(nsupMAX,nsupMAX),Rmbot(nsupMAX,nsupMAX),
     .                 Tmtop(nsupMAX,nsupMAX),Tmbot(nsupMAX,nsupMAX),
     .                 Rmsbot(nsupMAX,nsupMAX)

      DOUBLE PRECISION ebtop(nmuMAX),ebbot(nmuMAX)

Cf2py intent(in) ebtop,ebbot,Rmtop,Tmtop
Cf2py intent(out) Rmbot,Tmbot,Rmsbot

*----------------------------------------------------------------------
      CALL assign(Rmbot,Rmtop,nmat,nmu)
      CALL assign(Tmbot,Tmtop,nmat,nmu)
      CALL star(Rmsbot,Rmbot,nmat,nmu)

      DO i=1,nmu
         ebbot(i)= ebtop(i)
      ENDDO

*----------------------------------------------------------------------
      RETURN
      END
