* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE init(coefs,coefsm,coefsa)

************************************************************************
* DATE: December 2002
*
* AUTHOR: D. M. Stam
*
************************************************************************
      INCLUDE 'max_incl'

      INTEGER i,j,m,k
 
      DOUBLE PRECISION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 coefsa(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 coefsm(4,4,0:ncoefsMAX)

Cf2py intent(in,out) coefs, coefsm, coefsa

*-----------------------------------------------------------------------
      DO i=1,4
         DO j=1,4
            DO k=0,ncoefsMAX
               coefsm(i,j,k)= 0.D0
            ENDDO
         ENDDO
      ENDDO

      DO i=1,nmatMAX
         DO j=1,nmatMAX
            DO k=0,ncoefsMAX
               DO m=1,nlaysMAX
                  coefs(i,j,k,m)= 0.D0
                  coefsa(i,j,k,m)= 0.D0
               ENDDO
            ENDDO
         ENDDO
      ENDDO

************************************************************************
      RETURN
      END
