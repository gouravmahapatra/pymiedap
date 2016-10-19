      SUBROUTINE getgeos(alpha,npix,ngeos,apix,thet0,theta,phi,beta,
     .                   lat, long, xs, ys)

Cf2py intent(in) alpha, npix
Cf2py intent(out) ngeos, apix, thet0, theta, phi, beta, lat, long
Cf2py intent(out) xs, ys

*----------------------------------------------------------------------------
*     Open and read the geometries file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,ngeos,npix,nn

      DOUBLE PRECISION ph,alpha,angle,apix,bet,
     .                 calp,step,y,xlat,clat,slat,x,xlon,clon,slon,
     .                 alon,cthet,cthe0,sthet,sthe0,slst

      DOUBLE PRECISION t1,t2,t3,the,th0,xstar,zstar,
     .                 z,the2,aa1,aa2,th02,alph

      DOUBLE PRECISION thet0(ngeosMAX),theta(ngeosMAX),
     .                 phi(ngeosMAX), beta(ngeosMAX),
     .                 lat(ngeosMAX), long(ngeosMAX),
     .                 xs(ngeosMAX), ys(ngeosMAX)
 
      DOUBLE PRECISION pi,r2d,eps
      PARAMETER (pi=3.141592D0,r2d=180.D0/pi,eps=1.D-6)

*----------------------------------------------------------------------------
*     Prepare the phase angle calculations:
*----------------------------------------------------------------------------
      calp= DCOS(alpha/r2d)
      alpha = alpha/r2d
      step= 2.D0/DBLE(npix)

*----------------------------------------------------------------------------
*     The unit vector pointing to the star:
*----------------------------------------------------------------------------
      zstar= DCOS(alpha)
      xstar= DSIN(alpha)

*----------------------------------------------------------------------------
*     Loop over the y-direction (latitude):
*----------------------------------------------------------------------------
      ngeos=0
      nn=0

      DO i=1,npix

*        y is the middle of each pixel!
         y= -1.D0 + 0.5D0*step + DBLE(i-1)*step

*        xlat is the latitude of the middle of the pixel in radians!
         xlat= DASIN(y)

         slat= DSIN(xlat)

*----------------------------------------------------------------------------
*        Loop over the x-direction (longitude):
*        nn is the number of pixels across the whole disk
*----------------------------------------------------------------------------
         DO j=1,npix

*           x is the middle of each pixel!
            x= -1.D0 + 0.5D0*step + DBLE(j-1)*step

*           check if the pixel falls within the disk:
            IF ((x*x+y*y).LE.1.D0) THEN
                nn=nn+1
*----------------------------------------------------------------------------
*           Calculate the z-coordinate of the location on the sphere:
*----------------------------------------------------------------------------
                z= DSQRT(1.D0-x*x-y*y)
                xlon= DATAN(x/z)

*----------------------------------------------------------------------------
*           Calculate the local illumination angle theta0:
*----------------------------------------------------------------------------
                th0= DACOS(x*xstar+z*zstar)*r2d
                th02= DACOS(DCOS(xlat)*DCOS(alpha-xlon))*r2d
                    IF (th0.LT.90.D0) THEN
                        ngeos= ngeos+1

*----------------------------------------------------------------------------
*               Calculate the local viewing angle theta:
*----------------------------------------------------------------------------
                        the= DACOS(z)*r2d
                        the2= DACOS(DCOS(xlat)*DCOS(xlon))*r2d

*----------------------------------------------------------------------------
*               Calculate azimuthal angle phi-phi_0:
*----------------------------------------------------------------------------
                        t1= DCOS(alpha) - DCOS(the/r2d)*DCOS(th0/r2d)
                        t2= DSIN(the/r2d)*DSIN(th0/r2d)
    
                        IF (DABS(t2).LT.eps) THEN
                            ph= 0.D0
                        ELSE
                            t3= t1/t2
                            IF (t3.GT.1.D0)  t3=1.D0
                            IF (t3.LT.-1.D0) t3=-1.D0
                            ph= (pi-DACOS(t3))*r2d
                        ENDIF
                        IF (y.LT.0.D0) ph=-ph

                        theta(ngeos)= the
                        thet0(ngeos)= th0
                        phi(ngeos)= ph
                        lat(ngeos) = xlat*r2d
                        long(ngeos) = xlon*r2d
                        xs(ngeos) = x
                        ys(ngeos) = y

*----------------------------------------------------------------------------
*                 Calculate angle beta:
*----------------------------------------------------------------------------
                    cthet= DCOS(theta(ngeos)/r2d)
                    cthe0= DCOS(thet0(ngeos)/r2d)
                    sthet= DSIN(theta(ngeos)/r2d)
                    sthe0= DSIN(thet0(ngeos)/r2d)

                    IF (sthet.LE.0.D0) THEN
                        beta(ngeos)= 0.D0
                    ELSE
                        slst= slat/sthet
                        IF (slst.GT.1.D0) THEN
                            slst= 1.D0
                        ELSEIF (slst.LT.-1.D0) THEN
                            slst= -1.D0
                        ENDIF
                        bet= DASIN(slst)*r2d
                        IF (xlat.GE.0.D0) THEN
                            beta(ngeos)= bet
                            IF (xlon.LT.0.D0) beta(ngeos)= 180.D0-bet
                        ELSEIF (xlat.LT.0.D0) THEN
                            beta(ngeos)= 180.D0+bet
                            IF (xlon.LT.0.D0) beta(ngeos)= -bet
                        ENDIF
                    ENDIF

                ENDIF ! end of illumination test
            ENDIF !end of on-disk test
         ENDDO !end of x loop
      ENDDO !end of yloop

*----------------------------------------------------------------------------
*     Calculate the size of a pixel: 
*----------------------------------------------------------------------------
      apix= pi/nn
      WRITE(*,*)
      WRITE(*,*) 'apix:',apix,npix,ngeos
      WRITE(*,*)

*----------------------------------------------------------------------------
      RETURN
      END
