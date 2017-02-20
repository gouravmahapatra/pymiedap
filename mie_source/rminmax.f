      SUBROUTINE rminmax(idis,par1,par2,par3,weight2,cutoff,
     .                   rmin,rmax)

*-----------------------------------------------------------------------
*  Find the integration bounds rmin and rmax for the integration over 
*  a size distribution. These bounds are chosen such that the size   
*  distribution falls below the user specIFied cutoff. It is essential
*  that the size distribution is normalized such that the integral   
*  over all r is equal to one !                                     
*
*  This is programmed rather clumsy and will in the future be changed
*-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER idis

      DOUBLE PRECISION par1,par2,par3,cutoff,rmin,rmax,
     .                 sef,ref,rref,r0,r1,weight2

      DOUBLE PRECISION nwithr(1),r(1)

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-10)

Cf2py intent(in) idis, par1, par2, par3, weight2, cutoff
Cf2py intent(out) rmin, rmax

*-----------------------------------------------------------------------
*     No size distribution:
*-----------------------------------------------------------------------
      IF (idis.EQ.0) THEN
         rmin= par1
         rmax= par1

*-----------------------------------------------------------------------
*     Other distributions:
*-----------------------------------------------------------------------
      ELSE
         GOTO (10,20,30,40,50,60,70,80) idis

         WRITE(*,*) 'rminmax: illegal size distribution index :',idis
         STOP 

*-----------------------------------------------------------------------
*        Two parameter gamma:
*-----------------------------------------------------------------------
10       sef = 1.D0/DSQRT(par2+3.D0)
         ref = 1.D0/(sef*sef*par2)
         rref= ref
         GOTO 100

*-----------------------------------------------------------------------
*        Two parameter gamma:
*-----------------------------------------------------------------------
20       ref = par1
         sef = DSQRT(par2)
         rref= ref
         GOTO 100

*-----------------------------------------------------------------------
*        Bimodal gamma:
*-----------------------------------------------------------------------
30       sef = DSQRT(par3)
         ref = dmax1(par1,par2)+sef
         rref= dmin1(par1,par2)
         GOTO 100

*-----------------------------------------------------------------------
*        Log normal:
*-----------------------------------------------------------------------
40       sef = DSQRT(dexp(dlog(par2)**2)-1.d0)
         ref = par1*(1.D0+sef*sef)**0.4D0
         rref= ref
         GOTO 100

*-----------------------------------------------------------------------
*        Log normal:
*-----------------------------------------------------------------------
50       ref = par1
         sef = DSQRT(ref)
         rref= ref
         GOTO 100

*-----------------------------------------------------------------------
*        Power law:
*-----------------------------------------------------------------------
60       rmin= par2
         rmax= par3
         GOTO 999

*-----------------------------------------------------------------------
*        Modified gamma:
*-----------------------------------------------------------------------
70       ref = par2
         sef = 2.D0*ref
         rref=0.5D0*ref
         GOTO 100

*-----------------------------------------------------------------------
*        Modified gamma:
*-----------------------------------------------------------------------
80       ref = (par1/(par2*par3))**par3
         sef = 2.D0*ref
         rref= 0.5D0*ref

*-----------------------------------------------------------------------
*        Search for a value of r such that the size distribution is less
*        than the cutoff. Start the search at ref+sef which guarantees 
*        that such a value will be found on the TAIL of the 
*        distribution:                            
*-----------------------------------------------------------------------
100      r(1)= ref+sef
         r0  = ref
200      CALL sizedis(idis,par1,par2,par3,weight2,r,1,nwithr)

         IF (nwithr(1).GT.cutoff) THEN
            r0   = r(1)
            r(1) = 2.D0*r(1)
            GOTO 200
         ENDIF
         r1= r(1)

*-----------------------------------------------------------------------
*        Now the size distribution assumes the cutoff value somewhere
*        between r0 and r1  Use bisection to find the corresponding r:
*-----------------------------------------------------------------------
300      r(1)= 0.5D0*(r0+r1)
         CALL sizedis(idis,par1,par2,par3,weight2,r,1,nwithr)

         IF (nwithr(1).GT.cutoff) THEN
            r0= r(1)
         ELSE
            r1= r(1)
         ENDIF
         IF ((r1-r0).GT.eps) GOTO 300
         rmax = 0.5D0*(r0+r1)

*-----------------------------------------------------------------------
*        Search for a value of r on the low end of the size distribution 
*        such that the distribution falls below the cutoff. There is no 
*        guarantee that such a value exists, so use an extra test to 
*        see if the search comes very near to r=0:
*-----------------------------------------------------------------------
         r1 = rref
         r0 = 0.D0
400      r(1) = 0.5D0*r1
         CALL sizedis(idis,par1,par2,par3,weight2,r,1,nwithr)

         IF (nwithr(1).GT.cutoff) THEN
            r1 = r(1)
            IF (r1.GT.eps) GOTO 400
         ELSE
            r0 = r(1)
         ENDIF

*-----------------------------------------------------------------------
*        Possibly the size distribution goes through cutoff between r0 
*        and r1 try to find the exact value of r where this happens by  
*        bisection.                                                  
*
*        In case there is no solution, the algorithm will terminate:
*-----------------------------------------------------------------------
500      r(1) = 0.5D0*(r0+r1)
         CALL sizedis(idis,par1,par2,par3,weight2,r,1,nwithr)

         IF (nwithr(1).GT.cutoff) THEN
            r1= r(1)
         ELSE
            r0= r(1)
         ENDIF

         IF ((r1-r0).GT.eps) GOTO 500
         IF (r1.LE.eps) THEN
            rmin = 0.D0
         ELSE
            rmin = 0.5D0*(r0+r1)
         ENDIF
      ENDIF

*-----------------------------------------------------------------------
999   RETURN
      END
