      SUBROUTINE adding(outputname,a,b,coefs,ncoefs,nlays,
     .                  nmug,nmat,surfmat) 
**********************************************************************
*                A D D I N G    M E T H O D          
*          F O R    P O L A R I Z E D   L I G H T      
**********************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nlays,nmat,nmug,iunfou,nmu,iad,i,j,l,m
 
      INTEGER ncoefs(nlaysMAX),
     .        M0(nlaysMAX),M1(nlaysMAX),M2(nlaysMAX), 
     .        iadd(nlaysMAX,0:nfouMAX)

      DOUBLE PRECISION binc,binfac,bind,binfad
      DOUBLE PRECISION wav

      DOUBLE PRECISION a(nlaysMAX),b(nlaysMAX),
     .          coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX)

      DOUBLE PRECISION xmu(nmuMAX),smf(nmuMAX)

      DOUBLE PRECISION Rmbot(nsupMAX,nsupMAX),Rmtop(nsupMAX,nsupMAX),
     .                 Tmbot(nsupMAX,nsupMAX),Tmtop(nsupMAX,nsupMAX),
     .                 Rmsbot(nsupMAX,nsupMAX),
     .                 surfmat(nsupMAX,nsupMAX)

      DOUBLE PRECISION ebbot(nmuMAX),ebtop(nmuMAX)
      
      DOUBLE PRECISION Zmplus(nsupMAX,nsupMAX),Zmmin(nsupMAX,nsupMAX)

      LOGICAL nextm
      
      CHARACTER*16 outputname

Cf2py intent(in) outputname,a,b,coefs,ncoefs,nlays,nmug,nmat,surfmat

*---------------------------------------------------------------------
      binc= 1.D0
      bind= 0.D0
      iunfou= 123
*       Value set arbitrarily

      OPEN(unit=iunfou,file=outputname)

      WRITE(iunfou,330)
      WRITE(iunfou,331)
      WRITE(iunfou,332) wav
      WRITE(iunfou,333)


*---------------------------------------------------------------------
*     Initialize the mu-values and the supermatrix factors:
*---------------------------------------------------------------------
      CALL setmu(nmat,nmug,iunfou,nmu,xmu,smf)

*---------------------------------------------------------------------
*     Calculate the bounds M0, M1, and M2 on the Fourier-index:
*---------------------------------------------------------------------
      CALL setfou(coefs,ncoefs,nlays,a,b,xmu,nmug,M0,M1,M2,iadd)

*---------------------------------------------------------------------
*     Loop over the Fourier terms:
*---------------------------------------------------------------------
      m= -1
1000  m= m+1
C      WRITE(*,*) m

*---------------------------------------------------------------------
*     Calculate two binomial factors:
*---------------------------------------------------------------------
      IF (m.NE.0) THEN 
         binc= binc - binc/(2.D0*DBLE(m))
      ENDIF
      binfac= DSQRT(binc)

      IF (m.EQ.2) THEN
         bind= 1.D0/16.D0
      ELSEIF (m.GT.2) THEN
         bind= bind*(DBLE(m*m)-0.5D0*DBLE(m))/(DBLE(m*m)-4.D0)
      ENDIF
      binfad= -DSQRT(bind)

*---------------------------------------------------------------------
*     Initialization:
*---------------------------------------------------------------------
      DO i=1,nsupMAX
         DO j=1,nsupMAX
            Rmbot(i,j)= 0.D0
            Rmtop(i,j)= 0.D0
            Tmbot(i,j)= 0.D0
            Tmtop(i,j)= 0.D0
            Rmsbot(i,j)= 0.D0
         ENDDO
      ENDDO

*---------------------------------------------------------------------
*     Fill the arrays with the Lambertian surface reflection:
*---------------------------------------------------------------------
      IF (m.EQ.0) THEN
         CALL layer0(surfmat,smf,nmu,nmat,ebbot,Rmbot,Tmbot,Rmsbot)
      ENDIF

*---------------------------------------------------------------------
*     Loop over the atmospheric layers:
*---------------------------------------------------------------------
      DO l=1,nlays
         iad= iadd(l,m)

*---------------------------------------------------------------------
*        Calculate the m-th Fourier coefficient of the scat. matrix:
*---------------------------------------------------------------------
         CALL setzm(m,l,coefs,ncoefs,xmu,binfac,binfad,nmu,nmat,
     .              Zmmin,Zmplus)

         IF (m.EQ.0) CALL renorm(Zmmin,Zmplus,nmu,nmat,xmu,smf)

*---------------------------------------------------------------------
*        Calculate the m-th Fourier coefficient of the layer:
*---------------------------------------------------------------------
         CALL layerm(m,M0,M1,M2,l,xmu,smf,nmu,nmat,coefs,ncoefs,
     .               Zmplus,Zmmin,a(l),b(l),ebtop,Rmtop,Tmtop)

*---------------------------------------------------------------------
*        If this is the bottom layer, two actions may be taken:
*---------------------------------------------------------------------
         IF (l.EQ.1) THEN

*---------------------------------------------------------------------
*           If m=0, the new arrays are added to the surface arrays:
*---------------------------------------------------------------------
            IF (m.EQ.0) THEN

               CALL addlay(nmat,nmu,ebtop,ebbot,iad,Rmtop,Tmtop,
     .                  Rmbot,Tmbot,Rmsbot)

*---------------------------------------------------------------------
*           If m>0, the new arrays are shifted to the bottom arrays:
*---------------------------------------------------------------------
            ELSE
               CALL top2bot(nmat,nmu,ebtop,ebbot,Rmtop,Tmtop,
     .                      Rmbot,Tmbot,Rmsbot)
            ENDIF

*---------------------------------------------------------------------
*        If this is a higher layer, add the new to the bottom arrays:
*---------------------------------------------------------------------
         ELSE
            CALL addlay(nmat,nmu,ebtop,ebbot,iad,Rmtop,Tmtop,
     .                  Rmbot,Tmbot,Rmsbot)
         ENDIF

*---------------------------------------------------------------------
*     End of loop over atmospheric layers:
*---------------------------------------------------------------------
      ENDDO

*---------------------------------------------------------------------
*     Write the m-th Fourier coefficient to output file, 
*     and check the convergence of the Fourier series:
*---------------------------------------------------------------------
      CALL newfou(m,Rmbot,smf,iunfou,nmat,nmu,M0,nlays,nextm)

*---------------------------------------------------------------------
*     Next term of Fourier-loop:
*---------------------------------------------------------------------
      IF (m.LT.2) nextm=.true.
      IF (nextm) GOTO 1000

      CLOSE(iunfou)

*-----------------------------------------------------------------------
*     Formats:
*-----------------------------------------------------------------------
330   FORMAT('# For the correct use of this datafile, see Stam and',
     .       ' de Haan 2013')
331   FORMAT('# Please refer to Stam and de Haan 2013 when using',
     .       ' this data')
332   FORMAT('# The wavelength for the calculations:',F12.6)
333   FORMAT('#')


************************************************************************
      RETURN
      END
