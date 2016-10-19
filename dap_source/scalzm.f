      SUBROUTINE scalzm(m,layer,coefs,ncoefs,xmu,binfac,
     .                  nmu,Zmmin,Zmplus)

**********************************************************************
*  Calculate the m-th Fourier component of the phase matrix
*  Zm(mu0,mu) from the expansion coefficients of the    
*  scattering matrix into generalized spherical functions. 
*  The formulae can be found in :
*                                                         
*    J.F. de Haan et al.: 1987, Astron. Astrophys. 183, 
*     pp. 371-391. 
*                                                     
*  Essentially, Eqs. (66)-(82) are used.                
*                                                       
*  NO POLARIZATION IS INCLUDED !!                     
*
*  October 2014: Daphne changed the calculation of the
*                binomial factor to allow calculations
*                up to m=1000
*                For larger values of m, a different
*                algorithm should be used!
**********************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER m, layer, nmu
C      INTEGER ncoefs(nlaysMAX)
      DOUBLE PRECISION binfac

C      DOUBLE PRECISION Plm(nmuMAX,3,2),sqlm(0:ncoefsMAX),
C     .          coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
C     .          xmu(nmuMAX),
C     .          Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX)
      DIMENSION Plm(nmuMAX,3,2),sqlm(0:ncoefsMAX),
     .          coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .          ncoefs(nlaysMAX),xmu(nmuMAX),
     .          Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX)

      LOGICAL odd

Cf2py intent(in) m,layer,coefs,ncoefs,xmu,binfac,nmu
Cf2py intent(out) Zmmin,Zmplus 

*---------------------------------------------------------------------
*  Precompute the factor dsqrt(l**2-m**2) needed in Eqs. (81)-(82):
*---------------------------------------------------------------------
      m2= m*m
      sqlm(m)= 0.D0
      DO l=m+1,ncoefs(layer)
         sqlm(l)= DSQRT(DBLE(l)**2-DBLE(m2))
      ENDDO

*---------------------------------------------------------------------
*     Initialize Plm for l=m, Eq.(77) without factor i**-m 
*     this factor drops out anyway in Eq. (66)            
*---------------------------------------------------------------------
      lold= 1
      lnew= 2
      DO i=1,nmu
         u= xmu(i)
         rootu= DSQRT(DABS(1.D0-u*u))
         Plm(i,1,lold) = 0.D0
         IF (m.NE.0) THEN
             Plm(i,1,lnew)= binfac*rootu**m
         ELSE
             Plm(i,1,lnew)= binfac
         ENDIF
      ENDDO

*---------------------------------------------------------------------
*     Initialize phase matrix with term l=m, Eq.(66):
*---------------------------------------------------------------------
      DO i=1,nmu
         SP= coefs(1,1,m,layer)*Plm(i,1,lnew)
         DO j=1,nmu
            Zmplus(i,j)= SP*Plm(j,1,lnew)
            Zmmin(i,j) = Zmplus(i,j)
         ENDDO
      ENDDO

      IF (ncoefs(layer).GT.m) THEN
*---------------------------------------------------------------------
*         Start loop over l (summation index in Eq. (66))  
*---------------------------------------------------------------------
          odd= .false.

          DO 1200 l=m+1,ncoefs(layer)
              odd = .not. odd

*---------------------------------------------------------------------
*             Do one step in recurrence for Plm, Eq.(81) 
*---------------------------------------------------------------------
              c1= DBLE(l+l-1)/sqlm(l)
              c2= sqlm(l-1)/sqlm(l)
              DO i=1,nmu
                 u= xmu(i)
                 Plm(i,1,lold)= c1*u*Plm(i,1,lnew)-c2*Plm(i,1,lold)
              ENDDO
              itmp= lnew
              lnew= lold
              lold= itmp

*---------------------------------------------------------------------
*             Add a new term to Zm, Eq.(66)              
*---------------------------------------------------------------------
              IF (odd) THEN
                 DO i=1,nmu
                    SP = coefs(1,1,l,layer)*Plm(i,1,lnew)
                    DO j=1,nmu
                       Zmplus(i,j)= Zmplus(i,j)+ SP*Plm(j,1,lnew)
                       Zmmin(i,j)= Zmmin(i,j) - SP*Plm(j,1,lnew)
                    ENDDO
                 ENDDO
              ELSE
                 DO i=1,nmu
                    SP = coefs(1,1,l,layer)*Plm(i,1,lnew)
                    DO j=1,nmu
                       Zmplus(i,j)= Zmplus(i,j) + SP*Plm(j,1,lnew)
                       Zmmin(i,j)= Zmmin(i,j) + SP*Plm(j,1,lnew)
                    ENDDO
                 ENDDO
              ENDIF
1200       CONTINUE

*---------------------------------------------------------------------
*          End of summation loop over l      
*---------------------------------------------------------------------
      ENDIF

      RETURN
      END
