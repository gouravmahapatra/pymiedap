* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE obsmask(phi,lamb,phi0,lamb0,rot,npix,alpha,nlon,
     .                      nlat,mask)

Cf2py intent(in) phi,lamb,phi0,lamb0,rot,npix,alpha,nlon,nlat
Cf2py intent(out) mask
*----------------------------------------------------------------------------
*     Create an array with indexes of the disk pixels, applied on a
*           geographical grid.
*
*     Author: Ashwyn Groot
*     Date: December 2018
*----------------------------------------------------------------------------

*----------------------------------------------------------------------------
*     Define properties:
*----------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'max_incl'
      INTEGER  nlon,nlat,j,pixnum,npix,ngeos

      REAL*8    step,phi0,lamb0,pi,radfac,alpha,apix,rot
      REAL*8    thet0dum(ngeosMAX),thetadum(ngeosMAX),
     .          phidum(ngeosMAX), betadum(ngeosMAX),latdum(ngeosMAX),
     .          longdum(ngeosMAX),xvals(ngeosMAX), yvals(ngeosMAX),
     .          xrot(ngeosMAX),yrot(ngeosMAX)

      REAL*8    phi(nlat,nlon), lamb(nlat,nlon),
     .      x(nlat,nlon), y(nlat,nlon), cosc(nlat,nlon),
     .      y_upbound(nlat,nlon),y_botbound(nlat,nlon),
     .      x_upbound(nlat,nlon),x_botbound(nlat,nlon)

      REAL*8,DIMENSION(:,:),ALLOCATABLE :: xbounds,ybounds

      INTEGER     mask(nlat,nlon)

      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

*----------------------------------------------------------------------------
*     Prepare the mask calculations:
*----------------------------------------------------------------------------
      mask=-100
      phi=phi*radfac
      lamb=lamb*radfac
      phi0=phi0*radfac
      lamb0=lamb0*radfac
      step=2.D0/DBLE(npix)

*     Set disk pixel parameters
      call getgeos(alpha,npix,ngeos,apix,thet0dum,thetadum,phidum,
     .              betadum,latdum, longdum, xvals, yvals)

      ALLOCATE(xbounds(2,ngeos),ybounds(2,ngeos))
      xbounds(1,:)=xvals(:ngeos)-0.5*step
      xbounds(2,:)=xvals(:ngeos)+0.5*step
      ybounds(1,:)=yvals(:ngeos)-0.5*step
      ybounds(2,:)=yvals(:ngeos)+0.5*step

*     Rotate the disk pixel coordinates
      xrot=xvals*DCOS(rot*radfac)+yvals*DSIN(rot*radfac)
      yrot=-xvals*DSIN(rot*radfac)+yvals*DCOS(rot*radfac)

*     Convert lat/lon grid to y/x grid
      x=DCOS(phi)*DSIN(lamb-lamb0)
      y=DCOS(phi0)*DSIN(phi)-DSIN(phi0)*DCOS(phi)*DCOS(lamb-lamb0)
      cosc=DSIN(phi0)*DSIN(phi)+DCOS(phi0)*DCOS(phi)*DCOS(lamb-lamb0)

*     Apply negative mask for coordinates that lie outside the map projection
      WHERE (cosc.LT.0.D0)
         x=-100.D0
         y=-100.D0
      ENDWHERE

*     For no rotation or 180 degrees we apply the following:
      IF ((rot.EQ.0.D0).or.(rot.EQ.180.D0)) THEN
          IF (rot.EQ.0.D0) THEN
             pixnum=0
          ELSE IF (rot.EQ.180.D0) THEN
             pixnum=ngeos-1
          ENDIF
          DO j=1,ngeos
             WHERE ((x.GE.xbounds(1,j)).and.(x.LE.xbounds(2,j))
     .              .and.(y.GE.ybounds(1,j)).and.(y.LE.ybounds(2,j)))
                mask=pixnum
             ENDWHERE
             IF (rot.EQ.0.D0) THEN
                pixnum=pixnum+1
             ELSE IF (rot.EQ.180.D0) THEN
                pixnum=pixnum-1
             ENDIF
          ENDDO

*     For no -90,90 rotation the following is applied:
      ELSE IF ((rot.EQ.90.D0).or.(rot.EQ.-90.D0)) THEN
          IF (rot.EQ.90.D0) THEN
             pixnum=ngeos-1
          ELSE IF (rot.EQ.-90.D0) THEN
             pixnum=0
          ENDIF
          DO j=1,ngeos
             WHERE ((x.GE.ybounds(1,ngeos+1-j)).and.
     .                                      (x.LE.ybounds(2,ngeos+1-j))
     .              .and.(y.GE.xbounds(1,j)).and.
     .                                      (y.LE.xbounds(2,j)))
                mask=pixnum
             ENDWHERE
             IF (rot.EQ.90.D0) THEN
                pixnum=pixnum-1
             ELSE IF (rot.EQ.-90.D0) THEN
                pixnum=pixnum+1
             ENDIF
          ENDDO
      DEALLOCATE(xbounds,ybounds)

*     For other rotations the following is used, where we distinguish
*           positive and negative rotations:
      ELSE IF ((DABS(rot).NE.0.D0).or.(DABS(rot).NE.90.D0).or.
     .                              (DABS(rot).NE.180.D0)) THEN
          pixnum=ngeos-1
          DO j=1,ngeos
             IF (rot.GT.0.D0) THEN
                x_botbound=(x/DTAN(rot*radfac-pi/2))-
     .                      ((yvals(j)+0.5*step)/DCOS(rot*radfac))
                x_upbound=(x/DTAN(rot*radfac-pi/2))-
     .                      ((yvals(j)-0.5*step)/DCOS(rot*radfac))
                y_upbound=(x/DTAN(rot*radfac))+
     .                      ((xvals(j)+0.5*step)/DSIN(rot*radfac))
                y_botbound=(x/DTAN(rot*radfac))+
     .                      ((xvals(j)-0.5*step)/DSIN(rot*radfac))
                IF (DABS(rot).LT.90) THEN
                   WHERE ((y.LE.x_upbound).and.(y.GE.x_botbound)
     .              .and.(y.LE.y_upbound).and.(y.GE.y_botbound))
                      mask=pixnum
                   ENDWHERE
                ELSE IF (DABS(rot).GT.90) THEN
                   WHERE ((y.GE.x_upbound).and.(y.LE.x_botbound)
     .              .and.(y.LE.y_upbound).and.(y.GE.y_botbound))
                      mask=pixnum
                   ENDWHERE
                ENDIF
             ELSE IF (rot.LT.0.D0) THEN
                x_botbound=(x/DTAN(rot*radfac-pi/2))-
     .                      ((yvals(j)+0.5*step)/DCOS(rot*radfac))
                x_upbound=(x/DTAN(rot*radfac-pi/2))-
     .                      ((yvals(j)-0.5*step)/DCOS(rot*radfac))
                y_upbound=(x/DTAN(rot*radfac))+
     .                      ((xvals(j)+0.5*step)/DSIN(rot*radfac))
                y_botbound=(x/DTAN(rot*radfac))+
     .                      ((xvals(j)-0.5*step)/DSIN(rot*radfac))
                IF (DABS(rot).LT.90) THEN
                   WHERE ((y.LE.x_upbound).and.(y.GE.x_botbound)
     .              .and.(y.GE.y_upbound).and.(y.LE.y_botbound))
                      mask=pixnum
                   ENDWHERE
                ELSE IF (DABS(rot).GT.90) THEN
                   WHERE ((y.GE.x_upbound).and.(y.LE.x_botbound)
     .              .and.(y.GE.y_upbound).and.(y.LE.y_botbound))
                      mask=pixnum
                   ENDWHERE
                ENDIF
             ENDIF
             pixnum=pixnum-1
          ENDDO
      ENDIF
*     Apply a negative mask for the unvisible part of the map projection,
*           actually not necessary:
      WHERE (cosc.LT.0.D0)
         mask=-100
      ENDWHERE
*----------------------------------------------------------------------------
      RETURN
      END
