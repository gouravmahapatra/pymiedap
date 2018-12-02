* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE layerm(m,M0,M1,M2,l,xmu,smf,coefs,
     .                  ncoefs,Zmplus,Zmmin,a,b,ebtop,
     .                  Rmtop,Tmtop,nlays,nmu,nmat,nsup,max_ncoefs)

*----------------------------------------------------------------------
*     Calculate the m-th Fourier component of reflection and 
*     transmission of a homogeneous layer.  
*
*     When  0 <= m <= M2      use doubling
*     when  M2 < m <= M1      use first plus second order scattering
*     when  M1 < m <= M0      use first order scattering
*     when  M0 < m            use no scattering at all
*
*     Rm1top and Tm1top are used as scratch space!
*  Edited by: Ashwyn Groot
*  Date: November 2018
*  Introduced matrix operations with f95<
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER m,nmu,nmat,nsup,j,l,nlays,max_ncoefs

      INTEGER, DIMENSION(nlays) :: M0, M1, M2, ncoefs

      REAL*8, DIMENSION(nmu) :: xmu, smf, ebmu, ebtop

      REAL*8, DIMENSION(nsup,nsup) :: Zmmin,Zmplus,Rmtop,Tmtop,Rm1top,
     .              Tm1top

      REAL*8, DIMENSION(nmat,nmat,0:max_ncoefs,nlays) :: coefs

      LOGICAL verbo
      verbo = .false.

Cf2py intent(in) m,M0,M1,M2,l,xmu,smf,nmu,nmat,nsup,max_ncoefs,coefs,ncoefs
Cf2py intent(in,out) Zmplus,Zmmin,a,b,ebtop
Cf2py intent(in,out) Rmtop,Tmtop

*-----------------------------------------------------------------------
      IF (verbo) THEN
         WRITE(*,*) ' layerm: start calculating R and T for layer ',l
         WRITE(*,*) '     M0=',M0(l),' M1=',M1(l),' M2=',M2(l)
      ENDIF

      CALL expbmu(b,xmu,nmu,ebtop)

*-----------------------------------------------------------------------
*     No scattering:
*-----------------------------------------------------------------------
      IF (m.GT.M0(l)) THEN
         IF (verbo) WRITE(*,*) ' layerm: no scattering for m = ',m

         Rmtop=0.D0
         Tmtop=0.D0

         CALL expbmu(b,xmu,nmu,ebtop)

*-----------------------------------------------------------------------
*     One order of scattering:
*-----------------------------------------------------------------------
      ELSEIF (m.GT.M1(l)) THEN
         IF (verbo) WRITE(*,*) ' layerm: one order suffices for m = ',m

         CALL ord1m(xmu,smf,nmat,Zmplus,Zmmin,a,b,ebtop,nmu,nsup,
     .              Rm1top,Tm1top)

         Rmtop=Rm1top
         Tmtop=Tm1top
*-----------------------------------------------------------------------
*     Two orders of scattering:
*-----------------------------------------------------------------------
      ELSEIF (m.GT.M2(l)) THEN
         IF (verbo) WRITE(*,*) ' layerm: two orders suffice for m = ',m
         CALL ord1m(xmu,smf,nmat,Zmplus,Zmmin,a,b,ebtop,nmu,nsup,
     .              Rm1top,Tm1top)
         Rmtop=Rm1top
         Tmtop=Tm1top
         CALL ord2m(xmu,smf,nmat,Zmplus,Zmmin,a,b,ebtop,nmu,nsup,
     .              Rmtop,Tmtop)

*-----------------------------------------------------------------------
*     Doubling:
*-----------------------------------------------------------------------
      ELSE
         IF (verbo) WRITE(*,*) ' layerm: doubling needed for m = ',m

*-----------------------------------------------------------------------
*     Following call to ord1m can be commented because its output seems to be
*     unused? Loïc?
*-----------------------------------------------------------------------
         CALL ord1m(xmu,smf,nmat,Zmplus,Zmmin,a,b,ebtop,nmu,nsup,
     .              Rm1top,Tm1top)
*-----------------------------------------------------------------------

         xmumin= xmu(1)
         CALL bstart(m,l,coefs,ncoefs,M0,
     .              xmumin,a,b,nlays,max_ncoefs,nmat,b0,ndoubl)
         CALL expbmu(b0,xmu,nmu,ebmu)
         CALL ord1m(xmu,smf,nmat,Zmplus,Zmmin,a,b0,ebmu,nmu,nsup,
     .              Rmtop,Tmtop)
         CALL ord2m(xmu,smf,nmat,Zmplus,Zmmin,a,b0,ebmu,nmu,nsup,
     .              Rmtop,Tmtop)

         IF (ndoubl.GT.0) THEN
            bb = b0
            DO j=1,ndoubl
               CALL double(Rmtop,Tmtop,ebmu,nmu,nmat,nsup)
               bb= 2.D0*bb
               CALL expbmu(bb,xmu,nmu,ebmu)
            ENDDO
         ENDIF
      ENDIF

*-----------------------------------------------------------------------
      RETURN
      END
