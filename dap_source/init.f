* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE init(coefs,coefsm,coefsa)

************************************************************************
* DATE: December 2002
*
* AUTHOR: D. M. Stam
*
*  Edited by: Ashwyn Groot
*  Date: November 2018
*  Introduced matrix operations with f95<
************************************************************************
      INCLUDE 'max_incl'
 
      DOUBLE PRECISION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 coefsa(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 coefsm(4,4,0:ncoefsMAX)

Cf2py intent(in,out) coefs, coefsm, coefsa

*-----------------------------------------------------------------------
      coefsm=0.D0
      coefs=0.D0
      coefsa=0.D0

************************************************************************
      RETURN
      END
