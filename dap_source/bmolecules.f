* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE bmolecules(wav,nlays,pres,depol,ri,mma,
     .                      grav,bmsca,bmabs,coefsm)

************************************************************************
*     DATE: July 2014
*     AUTHOR: D. M. Stam
*
************************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,k,nlays

      DOUBLE PRECISION depol,ri,grav,mma,ma,dpres,wav,dep,ssca,
     .                 rindex,w2

      DOUBLE PRECISION coefsm(4,4,0:ncoefsMAX),
     .                 pres(nlaysMAX),
     .                 bmsca(nlaysMAX),bmabs(nlaysMAX)

*     Avogadro's number (in mole^-1):
      DOUBLE PRECISION avogad,pi
      PARAMETER (avogad=6.022169D23,pi=3.1415926535898D0)

* LOSCHMIDT'S' NUMBER (in m^-3 at 273.15 K and 1013.25 hPa):
      DOUBLE PRECISION loschmidt
      PARAMETER (loschmidt=2.54743D25)

Cf2py intent(in) wav, nlays, pres,depol,ri,mma,grav
Cf2py intent(out) bmsca, bmabs, coefsm
                      
*-----------------------------------------------------------------------
*     Some parameter values:
*
*     depol: depolarization value (CO2: 0.09)
*     mma: molecular mass (in amu)
*     grav: acceleration of gravity (m/s^2)
*-----------------------------------------------------------------------
C      depol= 0.09D0
C      mma= 44.0D0
C      grav= 8.87D0

*-----------------------------------------------------------------------
*     Calculate the refractive index "ri" for dry air at STP (288.15 K):
*     Peck and Reader (1972)
*-----------------------------------------------------------------------
C     w2= 1.D0/(wav*wav)
C     IF (wav.GT.0.23D0) THEN
C        ri= 1.D0 +
C    .       1.D-8*( 5791817.D0/ (238.0185D0 - w2)) +
C    .       1.D-8*( 167909.D0/ (57.362D0 - w2))
C     ELSE
C        ri= 1.D0 + 1.D-8*8060.51D0 +
C    .       1.D-8*( 2480990.D0/ (132.274D0 - w2) ) +
C    .       1.D-8*( 17455.7D0/ (39.32957D0 - w2) )
C     ENDIF

*-----------------------------------------------------------------------
*     Calculate the molecular scattering cross section (in m^2):
*-----------------------------------------------------------------------
      dep= (6.D0+3.D0*depol)/(6.D0-7.D0*depol)
      rindex= (ri*ri-1.D0)*(ri*ri-1.D0)/((ri*ri+2.D0)*(ri*ri+2.D0))

      ssca= 24.D0*pi*pi*pi*rindex*dep/
     .      (wav*wav*wav*wav*loschmidt*loschmidt*1.D-24)

*-----------------------------------------------------------------------
*     Calculate the column number density of the molecules (in m^-2), 
*     assuming hydrostatic equilibrium:
*-----------------------------------------------------------------------
      DO i=1,nlays
         ma= mma*1.D-3/avogad
         dpres= (pres(i)-pres(i+1))*10.D0**5
         bmsca(i)= ssca*dpres/(ma*grav)
         bmabs(i)= 0.D0
      ENDDO

*-----------------------------------------------------------------------
*     Calculate the molecular scattering expansion coefficients:
*-----------------------------------------------------------------------
      DO i=1,nmatMAX
         DO j=1,nmatMAX
            DO k=0,ncoefsMAX
               coefsm(i,j,k)= 0.D0
            ENDDO
         ENDDO
      ENDDO

      coefsm(1,1,0)= 1.D0
      coefsm(1,1,2)= (1.D0-depol)/(2.D0+depol)
      coefsm(2,2,2)= (6.D0-6.D0*depol)/(2.D0+depol)
      coefsm(2,1,2)= DSQRT(1.5D0)*(2.D0-2.D0*depol)/(2.D0+depol)
      coefsm(1,2,2)= coefsm(2,1,2)

************************************************************************
      RETURN
      END
