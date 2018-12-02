* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE addingascii(outputname,a,b,coefsin,ncoefs,max_ncoefs,
     .                  nmug,nmat,surfmatin,lg,nlays)
**********************************************************************
*                A D D I N G    M E T H O D          
*          F O R    P O L A R I Z E D   L I G H T      
**********************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nlays,nmat,nmug,nmu,iad,l,m,nsup,max_ncoefs,
     .                  lg,iunfou
 
      INTEGER ncoefs(nlays),
     .        M0(nlays),M1(nlays),M2(nlays),
     .        iadd(nlays,0:nfouMAX)

      DOUBLE PRECISION binc,binfac,bind,binfad,wav

      REAL*8 a(nlays),b(nlays),surfmatin(nsupMAX,nsupMAX),eps1,
     .          coefsin(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX)
      REAL*8, DIMENSION(nmat,nmat,0:max_ncoefs,nlays) :: coefs

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Rmbot, Rmtop,Tmbot,
     .          Tmtop,Rmsbot,surfmat,Zmplus,Zmmin,w !rank 2

      REAL*8, DIMENSION(:), ALLOCATABLE ::ebbot,ebtop,xmu,smf !rank 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: jb,nmu_i,nmu_j

      LOGICAL nextm
      CHARACTER(len=200) :: outputname

      coefs=coefsin(1:nmat,1:nmat,0:max_ncoefs,1:nlays)
Cf2py intent(in) outputname,a,b,coefs,ncoefs,max_ncoefs,nmug,nmat
Cf2py intent(in) surfmatin,lg,filetype,nlays
*---------------------------------------------------------------------
      binc= 1.D0
      bind= 0.D0
      iunfou= 123
      eps1= 1.D-100
*       Value set arbitrarily
      OPEN(unit=iunfou,file=outputname)

      WRITE(iunfou,330)
      WRITE(iunfou,331)
      WRITE(iunfou,332) wav
      WRITE(iunfou,333)

      ALLOCATE(xmu(nmug+1),smf(nmug+1))
*---------------------------------------------------------------------
*     Initialize the mu-values and the supermatrix factors:
*---------------------------------------------------------------------
      CALL setmuascii(nmug,iunfou,nmat,nmu,xmu,smf)
*---------------------------------------------------------------------
*     Define allocatable arrays here:
*---------------------------------------------------------------------
      nsup=nmat*nmu
      ALLOCATE(Rmbot(nsup,nsup),Rmtop(nsup,nsup),Tmbot(nsup,nsup),
     .              Tmtop(nsup,nsup),Rmsbot(nsup,nsup),
     .              surfmat(nsup,nsup),Zmplus(nsup,nsup),
     .              Zmmin(nsup,nsup),ebbot(nmu),ebtop(nmu),jb(nmu),
     .              w(nmat*nmu,nmu), nmu_i(nmu*nmat),nmu_j(nmu))
      surfmat(:,:)=surfmatin(:nsup,:nsup)
*---------------------------------------------------------------------
*     Calculate the bounds M0, M1, and M2 on the Fourier-index:
*---------------------------------------------------------------------
      CALL setfouascii(coefs,ncoefs,a,b,xmu,
     .              nmug,max_ncoefs,nlays,nmat,nmu,M0,M1,M2,iadd)

*---------------------------------------------------------------------
*     Loop over the Fourier terms:
*---------------------------------------------------------------------

      m= -1
1000  m= m+1
      WRITE(*,*) 'm= ',m

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
      Rmbot= 0.D0
      Rmtop= 0.D0
      Tmbot= 0.D0
      Tmtop= 0.D0
      Rmsbot= 0.D0

*---------------------------------------------------------------------
*     Fill the arrays with the Lambertian surface reflection:
*---------------------------------------------------------------------
      IF (m.EQ.0) THEN
         CALL layer0(surfmat,smf,nmat,ebbot,Rmbot,Tmbot,Rmsbot,
     .                  nmu,nsup)
      ENDIF

*---------------------------------------------------------------------
*     Loop over the atmospheric layers:
*---------------------------------------------------------------------
      DO l=1,nlays
         iad= iadd(l,m)

*---------------------------------------------------------------------
*        Calculate the m-th Fourier coefficient of the scat. matrix:
*---------------------------------------------------------------------
         CALL setzm(m,l,coefs,ncoefs,xmu,binfac,binfad,nsup,
     .              nlays,max_ncoefs,nmu,nmat,Zmmin,Zmplus)
         IF (m.EQ.0) CALL renorm(Zmmin,Zmplus,nmat,xmu,smf,nmu,nsup)
*---------------------------------------------------------------------
*        Calculate the m-th Fourier coefficient of the layer:
*---------------------------------------------------------------------
         CALL layerm(m,M0,M1,M2,l,xmu,smf,coefs,
     .              ncoefs,Zmplus,Zmmin,a(l),b(l),ebtop,
     .              Rmtop,Tmtop,nlays,nmu,nmat,nsup,max_ncoefs)
*---------------------------------------------------------------------
*        If this is the bottom layer, two actions may be taken:
*---------------------------------------------------------------------
         IF (l.EQ.1) THEN
*---------------------------------------------------------------------
*           If m=0, the new arrays are added to the surface arrays:
*---------------------------------------------------------------------
            IF (m.EQ.0) THEN

               CALL addlay(nmat,ebtop,ebbot,iad,Rmtop,Tmtop,nmu,nsup,
     .                  Rmbot,Tmbot,Rmsbot)

*---------------------------------------------------------------------
*           If m>0, the new arrays are shifted to the bottom arrays:
*---------------------------------------------------------------------
            ELSE
               CALL top2bot(nmat,ebtop,ebbot,Rmtop,Tmtop,nmu,nsup,
     .                      Rmbot,Tmbot,Rmsbot)
            ENDIF

*---------------------------------------------------------------------
*        If this is a higher layer, add the new to the bottom arrays:
*---------------------------------------------------------------------
         ELSE
            CALL addlay(nmat,ebtop,ebbot,iad,Rmtop,Tmtop,nmu,nsup,
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
      CALL newfouascii(m,Rmbot,smf,iunfou,nmat,M0,nmu,nsup,
     .                  nlays,nextm)

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
