* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE layer0(surfmat,nsurf,smf,nmu,nmat,ebbot,
     .                  Rmbot,Tmbot,Rmsbot,has_surface)

*----------------------------------------------------------------------
*     Fill the bottom-boundary arrays with the Fourier coefficient of
*     the surface reflection matrix.  Older PyMieDAP versions only used
*     the intensity-to-intensity element and only for m=0.  This version
*     accepts a full Stokes supermatrix coefficient for any Fourier term,
*     which is required for non-Lambertian rough-ocean glint.
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,nsup,nmu,nmat,nsurf

      DOUBLE PRECISION ebbot(nmuMAX),smf(nmuMAX),
     .                 Rmbot(nsupMAX,nsupMAX),Tmbot(nsupMAX,nsupMAX),
     .                 Rmsbot(nsupMAX,nsupMAX),
     .                 surfmat(nsurf,nsurf)

      LOGICAL has_surface

Cf2py intent(in) surfmat,nsurf,smf,nmu,nmat
Cf2py intent(out) ebbot,has_surface
Cf2py intent(in,out) Rmbot,Tmbot,Rmsbot

      nsup= nmu*nmat
      has_surface= .false.

      DO i=1,nsup
         DO j=1,nsup
            Rmbot(i,j)= 0.D0
            Tmbot(i,j)= 0.D0
            Rmsbot(i,j)= 0.D0
         ENDDO
      ENDDO

      DO i=1,nsup
         DO j=1,nsup
            Rmbot(i,j)= surfmat(i,j)
            IF (DABS(surfmat(i,j)).GT.eps) has_surface= .true.
         ENDDO
      ENDDO

      CALL star(Rmsbot,Rmbot,nmat,nmu)

*----------------------------------------------------------------------
*     Fill array ebbot (the direct transmission through the surface),
*     which equals zero, but is needed anyway.
*----------------------------------------------------------------------
      DO i=1,nmu
         ebbot(i)= 0.D0
      ENDDO

      RETURN
      END
