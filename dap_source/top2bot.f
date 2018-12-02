* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE top2bot(nmat,ebtop,ebbot,Rmtop,Tmtop,nmu,nsup,
     .                   Rmbot,Tmbot,Rmsbot)

*----------------------------------------------------------------------
*     Copy data pertaining to the top layer to the bottom layer:
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nmat,nmu, nsup

      REAL*8, DIMENSION(nsup,nsup) :: Rmtop,Tmtop,Rmbot,
     .          Tmbot,Rmsbot

      REAL*8, DIMENSION(nmu) :: ebtop, ebbot

Cf2py intent(in) ebtop,ebbot,Rmtop,Tmtop
Cf2py intent(out) Rmbot,Tmbot,Rmsbot

*----------------------------------------------------------------------
      Rmbot=Rmtop
      Tmbot=Tmtop

      CALL star(Rmsbot,Rmbot,nmat,nmu,nsup)

      ebbot=ebtop

*----------------------------------------------------------------------
      RETURN
      END
