      SUBROUTINE layerm(m,M0,M1,M2,l,xmu,smf,nmu,nmat,coefs,ncoefs,
     .                  Zmplus,Zmmin,a,b,ebtop,
     .                  Rmtop,Tmtop)

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
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER i,m,nmu,nmat,nsup,j,l

      INTEGER M0(nlaysMAX),M1(nlaysMAX),M2(nlaysMAX),
     .        ncoefs(nlaysMAX) 

      DOUBLE PRECISION xmu(nmuMAX),smf(nmuMAX)
C      DOUBLE PRECISION a, b

      DOUBLE PRECISION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX), 
     .                 Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX),
     .                 Rmtop(nsupMAX,nsupMAX), Tmtop(nsupMAX,nsupMAX),
     .                 Rm1top(nsupMAX,nsupMAX),Tm1top(nsupMAX,nsupMAX)

      DOUBLE PRECISION ebmu(nmuMAX),ebtop(nmuMAX)

      LOGICAL verbo
      verbo = .false.

Cf2py intent(in) m,M0,M1,M2,l,xmu,smf,nmu,nmat,coefs,ncoefs
Cf2py intent(in,out) Zmplus,Zmmin,a,b,ebtop
Cf2py intent(in,out) Rmtop,Tmtop

*-----------------------------------------------------------------------
      nsup= nmu*nmat

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

         DO j=1,nsup
            DO i=1,nsup
               Rmtop(i,j)= 0.D0
               Tmtop(i,j)= 0.D0
            ENDDO
         ENDDO

         CALL expbmu(b,xmu,nmu,ebtop)

*-----------------------------------------------------------------------
*     One order of scattering:
*-----------------------------------------------------------------------
      ELSEIF (m.GT.M1(l)) THEN
         IF (verbo) WRITE(*,*) ' layerm: one order suffices for m = ',m

         CALL ord1m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b,ebtop,
     .              Rm1top,Tm1top)

         CALL assign(Rmtop,Rm1top,nmat,nmu)
         CALL assign(Tmtop,Tm1top,nmat,nmu)

*-----------------------------------------------------------------------
*     Two orders of scattering:
*-----------------------------------------------------------------------
      ELSEIF (m.GT.M2(l)) THEN
         IF (verbo) WRITE(*,*) ' layerm: two orders suffice for m = ',m

         CALL ord1m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b,ebtop,
     .              Rm1top,Tm1top)

         CALL assign(Rmtop,Rm1top,nmat,nmu)
         CALL assign(Tmtop,Tm1top,nmat,nmu)

         CALL ord2m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b,ebtop,
     .              Rmtop,Tmtop)

*-----------------------------------------------------------------------
*     Doubling:
*-----------------------------------------------------------------------
      ELSE
         IF (verbo) WRITE(*,*) ' layerm: doubling needed for m = ',m

         CALL ord1m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b,ebtop,
     .              Rm1top,Tm1top)

         xmumin= xmu(1)

         CALL bstart(m,l,coefs,ncoefs,M0,xmumin,a,b,
     .               b0,ndoubl)
         CALL expbmu(b0,xmu,nmu,ebmu)
         CALL ord1m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b0,ebmu,
     .              Rmtop,Tmtop)

         CALL ord2m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b0,ebmu,
     .              Rmtop,Tmtop)

         IF (ndoubl.GT.0) THEN
            bb = b0
            DO j=1,ndoubl
               CALL double(Rmtop,Tmtop,ebmu,nmu,nmat)
               bb= 2.D0*bb
               CALL expbmu(bb,xmu,nmu,ebmu)
            ENDDO
         ENDIF

      ENDIF

*-----------------------------------------------------------------------
      RETURN
      END
