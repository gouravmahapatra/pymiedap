      SUBROUTINE scatmat(m1,m2,wav,idis,thmin,thmax,step,nsubr,ngaur,
     .                   rmin,rmax,par1,par2,par3,ratio,weight2,delta,
     .                   u,wg,F,miec,nangle)

************************************************************************
*  Calculate the scattering matrix of an ensemble of homogenous        *
*  spheres. On entry, the following must be supplied :                 *
*     m            : complex index of refraction                       *
*     wav       : wavelength                                        *
*     idis         : index of the size distribution                    *
*     nsubr         : number of subintervals for integration over r     *
*     ngaur       : number of Gauss points used per subinterval       *
*     rmin         : lower bound for integration over r                *
*     rmax         : upper bound for integration over r                *
*     par1,2,3     : parameters of the size distribution
*     ratio       : ratio of inner radius to outer radius              *
*     weight2     : weight of the second mode in case of bimodal dist  *
*     delta        : cutoff used in truncation of the Mie sum          *
*     thmin        : minimum scattering angle in degrees               *
*     thmax        : maximum scattering angle in degrees               *
*     step         : step in scattering angle in degrees               *
*  On exit, the following results are RETURNed :                       *
*     u            : cosines of scattering angles                      *
*     wg          : Gaussian weights associated with u                *
*     F            : scattering matrix for all cosines in u            *
*     miec         : array containing cross sections etc.              *
*     nangle       : the number of scattering angles                   *
************************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER idis,nsubr,ngaur,nangle,nhalf,nstop,nfi,nD,n,
     .        nfou,l,i,nmax,nfac,j,k

      INTEGER NDn,NDr
      PARAMETER (NDn=50000,NDr=10000)

      DOUBLE PRECISION wav,thmin,thmax,step,rmin,rmax,
     .                 par1,par2,par3,weight2,delta,numpar,
     .                 volume,xeff,albedo,Qext,Qsca,aux,pie,
     .                 Cscasum,Cextsum,zabs,x,y,sw,dr,fac90,
     .                 reff,G,Cext,Csca,fakt,rtox,radfac,
     .                 ratio

      DOUBLE PRECISION F(nmatMAX,nangMAX),u(nangMAX),wg(nangMAX),
     .                 miec(13),nwithr(NDr),pi(NDn),tau(NDn),
C     .                 fi(0:NDn),chi(0:NDn),D(NDn),r(NDr),w(NDr),
     .                 fi(0:NDn),chi(0:NDn),r(NDr),w(NDr),
     .                 facf(NDn),facb(NDn)

      DOUBLE COMPLEX m1,m2,ci,Splusf,Sminf,cSplusf,cSminf,Splusb,
     .               Sminb,cSplusb,cSminb

      DOUBLE COMPLEX an(NDn),bn(NDn)
      DOUBLE COMPLEX D(NDn)

      LOGICAL     symth

Cf2py intent(in) m, wav, idis, thmin, thmax, step, nsubr, ngaur, rmin
Cf2py intent(in) rmax, par1, par2, par3, weight2, delta
Cf2py intent(out) u, wg, F, miec, nangle 

*-----------------------------------------------------------------------
*     Constants:
*        nfac is the number of precomputed factors (2n+1)/(n*(n+1))
*-----------------------------------------------------------------------
      pie  = dacos(-1.d0)
      radfac= pie/180.D0
      rtox = 2.D0*pie/wav
      fakt = wav*wav/(2.D0*pie)
      nfac = 0

*-----------------------------------------------------------------------
*     Initialise:
*-----------------------------------------------------------------------
      DO j=1,nangMAX
         DO k=1,nmatMAX
            F(k,j)=0.D0
         ENDDO
      ENDDO

      Csca  = 0.D0
      Cext  = 0.D0
      numpar= 0.D0
      G     = 0.D0
      reff  = 0.D0
      nfou  = 0
      fac90 = 1.D0
      ci    = dcmplx(0.D0,1.D0)

      symth = .true.

      IF (.NOT.symth) THEN
         WRITE(*,*) ' tstsym: theta points NOT symmetrical !!'
      ENDIF

*-----------------------------------------------------------------------
*     Distinguish between distribution or not:
*-----------------------------------------------------------------------
      IF (idis.EQ.0) THEN
         w(1)  = 1.D0
         r(1)  = rmin
         nwithr(1)= 1.D0
         nsubr = 1
         ngaur = 1
         dr = 0.D0
      ELSE
         dr= (rmax-rmin)/DBLE(nsubr)
         CALL gauleg(ngaur,ngaur,(rmax-dr),rmax,r,w)
         CALL sizedis(idis,par1,par2,par3,weight2,r,ngaur,nwithr)
      ENDIF

*-----------------------------------------------------------------------
*     Start integration over radius r with largest radius:
*-----------------------------------------------------------------------
      DO 60 l=nsubr,1,-1
         write(*,*) l
         DO 50 i=ngaur,1,-1
            write(*,*) i

            sw= nwithr(i)*w(i)
            x = rtox*r(i) !x is for mantle
            y = rtox*r(i)*ratio !y is for core
            !nmax = x + 4.05D0*x**(1.D0/3.D0) + 20
            nmax= x+4.*x**0.3333+2.0
            nfi = nmax+60
            zabs = x*cdabs(m1)
            nD = zabs + 4.05D0*zabs**(1.D0/3.D0) + 70
            IF ((nD.GT.NDn).OR.(nfi.GT.NDn)) THEN
               WRITE(*,*) 'scatmat: estimated number of Mie-terms:',nD
               WRITE(*,*) '         for particle sizeparameter :',x
               WRITE(*,*) '         maximum NDn is only        : ',NDn
               STOP
            ENDIF

            write (*,*) 'm1=',m1
            write (*,*) 'm2=',m2
            CALL bhcoat(x,y,m1,m2,an,bn)

*-----------------------------------------------------------------------
*           Precompute the factor (2n+1)/(n*(n+1)) needed in Mie 
*           sum over n:
*-----------------------------------------------------------------------
            IF (nmax.GT.nfac) THEN
               DO n=nfac+1,nmax
                  facf(n)= dble(2*n+1)/dble(n*(n+1))
                  facb(n)= facf(n)
                  IF (mod(n,2).EQ.1) facb(n)= -facb(n)
               ENDDO
               nfac= nmax
            ENDIF

*-----------------------------------------------------------------------
*           Calculate extinction and scattering cross section 
*           Use the convergence criterion to determine the number of 
*           terms that will later be used in the Mie sum for the 
*           scattering matrix itself:
*-----------------------------------------------------------------------
            Cextsum= 0.D0
            Cscasum= 0.D0
            nstop= nmax
            DO n=1,nmax
               aux= (2.D0*dble(n)+1.D0) *
     .              dabs(dble(an(n)*conjg(an(n)) + bn(n)*conjg(bn(n))))
               Cscasum = Cscasum + aux
               Cextsum = Cextsum + (2.D0*n+1.D0)*dble(an(n)+bn(n))
               IF (aux.LT.delta) THEN
                  nstop= n
                  GOTO 53
               ENDIF
            ENDDO
53          nfou= nstop
            IF (nfou.GE.nmax) THEN
         WRITE(*,*) ' WARNING from scatmat : Mie sum not converged for'
     .             ,' scattering cross section'
         WRITE(*,*) '   radius r = ',r(i),' sizeparameter x = ',x
     .             ,' sizedistribution nr. ',idis
         WRITE(*,*) '   Re(m) = ',dble(m1),' Im(m) = ',dimag(m1)
         WRITE(*,*) '   a priori estimate of number of Mie terms:',nmax
         WRITE(*,*) '   term ',nmax,' for Csca was ',aux
         WRITE(*,*) '   should have been less than ',delta
         WRITE(*,*) '   the a priori estimate will be used'
            ENDIF

*-----------------------------------------------------------------------
*           Only for the first run through the loop set points in 
*           u= dcos(th):
*-----------------------------------------------------------------------
            IF ((l.EQ.nsubr).AND.(i.EQ.ngaur)) THEN

*-----------------------------------------------------------------------
*              Expansion in GSF : set Gauss points for dcos(th) 
*              Ensure accurate integrations: add two terms: 
*              nangle = 2*nfou+2
*              One should be sufficient, but total should be even!
*-----------------------------------------------------------------------
               nangle= 2*nfou+2
               IF (nangle.GT.nangMAX) THEN
                  WRITE(*,*) 'scatmat: need too many int. angles'
     .                      ,' nangle=',nangle
                  WRITE(*,*) '       maximum array size= ',nangMAX
                  STOP
               ENDIF
               CALL gauleg(nangle,nangle,-1.d0,1.D0,u,wg)
            ENDIF

*-----------------------------------------------------------------------
*           Integration for normalization of size distibution, 
*           geometrical cross section and effective radius:
*-----------------------------------------------------------------------
            numpar= numpar+sw
            G     = G     +sw*r(i)*r(i)
            reff  = reff  +sw*r(i)*r(i)*r(i)

************************************************************************
*  Start loop over scattering angles, DO only half and use symmetry    *
*  between forward and backward scattering angles                      *
*  The factor fac90 will later be used to correct for the fact that    *
*  for a symmetrical set of scattering angles with an odd number of    *
*  angles the scattering matrix is a factor 2 too big at 90 degrees    *
*  due to the way we programmed the symmetry relations                 *
************************************************************************
            IF (symth) THEN
               IF (mod(nangle,2).EQ.1) THEN
                  nhalf= (nangle+1)/2
                  fac90= 0.5D0
               ELSE
                  nhalf= nangle/2
               ENDIF

               DO j=1,nhalf
                  CALL pitau(u(j),nmax,pi,tau)
                  Splusf= dcmplx(0.D0,0.D0)
                  Sminf = dcmplx(0.D0,0.D0)
                  Splusb= dcmplx(0.D0,0.D0)
                  Sminb = dcmplx(0.D0,0.D0)

*  THIS IS THE INNERMOST LOOP !! (Mie sum)
*  can be programmed more efficiently by taking the facf multiplication
*  outside the angle loop over index j 

           DO n=1,nfou
              Splusf= Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
              Sminf = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
              Splusb= Splusb + facb(n)*(an(n)+bn(n)) * (pi(n)-tau(n))
              Sminb = Sminb  + facb(n)*(an(n)-bn(n)) * (pi(n)+tau(n))
           ENDDO
           cSplusf= conjg(Splusf)
           cSminf = conjg(Sminf )
           cSplusb= conjg(Splusb)
           cSminb = conjg(Sminb )
           k= nangle-j+1
*  the forward scattering elements
           F(1,j)= F(1,j) +    sw*(Splusf*cSplusf + Sminf *cSminf)
           F(2,j)= F(2,j) -    sw*(Sminf *cSplusf + Splusf*cSminf)
           F(3,j)= F(3,j) +    sw*(Splusf*cSplusf - Sminf *cSminf)
           F(4,j)= F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
*  the backward scattering elements
           F(1,k)= F(1,k) +    sw*(Splusb*cSplusb + Sminb *cSminb)
           F(2,k)= F(2,k) -    sw*(Sminb *cSplusb + Splusb*cSminb)
           F(3,k)= F(3,k) +    sw*(Splusb*cSplusb - Sminb *cSminb)
           F(4,k)= F(4,k) + ci*sw*(Sminb *cSplusb - Splusb*cSminb)
        ENDDO
      ELSE
*-----------------------------------------------------------------------
*       Start loop over scattering angles, do all angles:
*-----------------------------------------------------------------------
        DO j=1,nangle
           CALL pitau(u(j),nmax,pi,tau)
           Splusf= dcmplx(0.D0,0.D0)
           Sminf = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
           DO n=1,nfou
              Splusf= Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
              Sminf = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
           ENDDO
           cSplusf= conjg(Splusf)
           cSminf = conjg(Sminf )
           k = nangle-j+1
*  the forward scattering elements
           F(1,j)= F(1,j) +      sw*(Splusf*cSplusf + Sminf *cSminf)
           F(2,j)= F(2,j) -      sw*(Sminf *cSplusf + Splusf*cSminf)
           F(3,j)= F(3,j) +      sw*(Splusf*cSplusf - Sminf *cSminf)
           F(4,j)= F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
        ENDDO
      ENDIF

*-----------------------------------------------------------------------
*   Integration for cross sections, shift radius to next subinterval:
*-----------------------------------------------------------------------
            Csca= Csca + sw*Cscasum
            Cext= Cext + sw*Cextsum
            r(i)= r(i) - dr
50       CONTINUE
         IF (l.NE.1) CALL sizedis(idis,par1,par2,par3,weight2,
     .                            r,ngaur,nwithr)
60    CONTINUE

*-----------------------------------------------------------------------
*     End of integration over size distribution  
*     Some final corrections :   
*-----------------------------------------------------------------------
      DO j=1,nangle
         DO k=1,nmatMAX
            F(k,j)= F(k,j)/(2.D0*Csca)
         ENDDO
      ENDDO

      IF (symth) THEN
         DO k=1,nmatMAX
            F(k,nhalf)= fac90*F(k,nhalf)
         ENDDO
      ENDIF

      G     = pie*G
      Csca  = Csca*fakt
      Cext  = Cext*fakt
      Qsca  = Csca/G
      Qext  = Cext/G
      albedo= Csca/Cext
      volume= (4.d0/3.d0)*pie*reff
      reff  = pie*reff/G
      xeff  = rtox*reff

      miec(1) = Csca
      miec(2) = Cext
      miec(3) = Qsca
      miec(4) = Qext
      miec(5) = albedo
      miec(6) = G
      miec(7) = reff
      miec(8) = xeff
      miec(9) = numpar
      miec(10)= volume

*-----------------------------------------------------------------------
      RETURN
      END

