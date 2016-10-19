      SUBROUTINE setzm(m,layer,coefs,ncoefs,xmu,binfac,binfad,
     .                 nmu,nmat,Zmmin,Zmplus)

**********************************************************************
*  Calculate the m-th Fourier component of the phase matrix 
*  Zm(mu0,mu) from the expansion coefficients of the
*  scattering matrix intgeneralized spherical functions. 
*  The formulae can be found in : 
*                                                
*    J.F. de Haan et al.: 1987, Astron. Astrophys. 183, 
*    pp. 371-391.
*                                                
*  Essentially, Eqs. (66)-(82) are used. The suggestion below 
*  Eq. (82) to diagonalize the matrix Plm is followed here. 
*  To this end we define the two matrices :     
*                                                 
*      ( 1     0     0     0 )      ( 1     0     0     0 )
*  D1= ( 0     1     1     0 )  D2= ( 0    0.5   0.5    0 )   
*      ( 0     1    -1     0 )      ( 0    0.5  -0.5    0 )  
*      ( 0     0     0     1 )      ( 0     0     0     1 ) 
*                                                         
*  It is clear that D1 and D2 are each other's inverse.  
*  We diagonalize the Plm given in Eq. (74) :        
*                                                        
*           D1*Plm*D2  =  diag( Plm0,  Plm-2,  Plm+2,  Plm0 ) 
*                                                
*  The sum in Eq. (66) is written :            
*                                                  
*  sum(Plm*Sl*Plm)=D1*sum(D2*Plm*D1*D2*Sl*D2*D1*Plm*D2)*D1
*                                                           
*  The matrix D2*Sl*D2 for a given l is given by :           
*                                                            
*  ( alpha1      beta1/2           beta1/ 2         0     ) 
*  ( beta1/2 (alpha2+alpha3)/4 (alpha2-alpha3)/4  beta2/2 ) 
*  ( beta1/2 (alpha2-alpha3)/4 (alpha2+alpha3)/4 -beta2/2 ) 
*  (   0        -beta2/2           beta2/2        alpha4  )
*                                                      
*  where the alpha's and beta's are given by Eqs. (68)-(73) 
**********************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER m, layer, nmu, nmat
C      INTEGER ncoefs
      DOUBLE PRECISION binfac,binfad
C      DOUBLE PRECISION coefs, xmu, Zmmin, Zmplus

      DIMENSION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .          cs(4,4,0:ncoefsMAX),
     .          ncoefs(nlaysMAX), 
     .          xmu(nmuMAX),
     .          Zmmin(nsupMAX,nsupMAX),
     .          Zmplus(nsupMAX,nsupMAX)

      DIMENSION Plm(nmuMAX,3,2),DPDpl(nsupMAX),
     .          DPDmi(nsupMAX),DSD(4,4),
     .          sqlm(0:ncoefsMAX),sql4(ncoefsMAX),
     .          rootu(nmuMAX)

      LOGICAL verbo 
      verbo = .false.

Cf2py intent(in) m,layer,coefs,ncoefs,xmu,binfac,binfad,nmu,nmat
Cf2py intent(out) Zmmin,Zmplus 

*---------------------------------------------------------------------
*     Initialize things:
*---------------------------------------------------------------------
      DO i=1,nmat
         DO j=1,nmat
            DO k=0,ncoefs(layer)
               cs(i,j,k)= coefs(i,j,k,layer)
            ENDDO
         ENDDO
      ENDDO

*---------------------------------------------------------------------
*     No polarization, use the scalar phase matrix:
*---------------------------------------------------------------------
      qroot6 = -0.25D0*DSQRT(6.D0)
      IF (nmat.EQ.1) THEN
          IF (verbo) print *,' setzm: use scalar phase matrix'
          CALL scalzm(m,layer,coefs,ncoefs,xmu,binfac,
     .                nmu,nmat,Zmmin,Zmplus)
          GOTO 999
      ENDIF

*---------------------------------------------------------------------
*     Precompute the factor DSQRT(l**2-m**2) needed in Eqs.(81)-(82) 
*     and also the factor DSQRT(l**2-4) needed in Eq. (82):
*---------------------------------------------------------------------
      m2= m*m
      DO l=m,ncoefs(layer)
         sqlm(l)= DSQRT(DABS(DBLE(l)**2-DBLE(m2)))
      ENDDO
      DO l=2,ncoefs(layer)
         sql4(l)= DSQRT(DABS(DBLE(l)**2-4.D0))
      ENDDO

*---------------------------------------------------------------------
*     Initialize phase matrix to zero:
*---------------------------------------------------------------------
      nsup= nmat*nmu
      DO j=1,nsup
         DO i=1,nsup
            Zmplus(i,j)= 0.D0
            Zmmin(i,j) = 0.D0
         ENDDO
      ENDDO

*---------------------------------------------------------------------
*     Initialize Plm0 for l=m, Eq.(77), and for l=m-1 Eq. (76):
*---------------------------------------------------------------------
      lold= 1
      lnew= 2
      DO i=1,nmu
         rootu(i)= DSQRT(DABS(1.D0-xmu(i)**2))
         IF (m.NE.0) THEN
             Plm(i,1,lnew)= binfac*rootu(i)**m
         ELSE
             Plm(i,1,lnew)= binfac
         ENDIF
         Plm(i,1,lold)= 0.D0
      ENDDO

*---------------------------------------------------------------------
*     Set Plm2 and Plm-2 to zero: initialization will be done 
*     inside loop 
*---------------------------------------------------------------------
      DO i=1,nmu
         Plm(i,2,lnew) = 0.D0
         Plm(i,2,lold) = 0.D0
         Plm(i,3,lnew) = 0.D0
         Plm(i,3,lold) = 0.D0
      ENDDO

*---------------------------------------------------------------------
*     Start loop over l (summation index in Eq. (66))             
*     Parity of Plm is (-1)**(l-m)                         
*---------------------------------------------------------------------
      parity= -1.D0
      DO 1200 l=m,ncoefs(layer)
          parity= -parity

*---------------------------------------------------------------------
*  Initialize Plm for l=max(m,2) Eqs.(78)-(80) without 
*  factor i**-m. This factor is cancelled in Eq. (66) by the
*  factor -1**m. The exception for m=2 is needed to handle 
*  u=+-1 in Eq. (80)      
*---------------------------------------------------------------------
          IF (l.EQ.max0(m,2)) THEN
              IF (m.EQ.0) THEN
                  DO i=1,nmu
                     Plm(i,2,lnew)=qroot6*rootu(i)*rootu(i)
                     Plm(i,3,lnew) = Plm(i,2,lnew)
                  ENDDO
              ELSE IF (m.EQ.1) THEN
                  DO i=1,nmu
                     u = xmu(i)
                     Plm(i,2,lnew)=-0.5D0*rootu(i)*(1.D0-u)
                     Plm(i,3,lnew)= 0.5D0*rootu(i)*(1.D0+u)
                  ENDDO
              ELSE IF (m.EQ.2) THEN
                  DO i=1,nmu
                     u = xmu(i)
                     Plm(i,2,lnew) = -0.25D0*(1.D0-u)**2
                     Plm(i,3,lnew) = -0.25D0*(1.D0+u)**2
                  ENDDO
              ELSE
                  DO i=1,nmu
                     u = xmu(i)
                     urootm = rootu(i)**(m-2)
                     Plm(i,2,lnew) = binfad*urootm*(1.D0-u)*
     .                            (1.D0-u)
                     Plm(i,3,lnew) = binfad*urootm*(1.D0+u)*
     .                            (1.D0+u)
                  ENDDO
              ENDIF
          ENDIF

*---------------------------------------------------------------------
*         Construct supervectors corresponding to the 
*         diagonal elements of the matrix D1 * Plm * D2 
*---------------------------------------------------------------------
          DO i=1,nmu
             isup = nmat*(i-1)
             DPDpl(isup+1) = Plm(i,1,lnew)
             DPDpl(isup+2) = Plm(i,2,lnew)
             DPDpl(isup+3) = Plm(i,3,lnew)
             DPDmi(isup+1) = parity*Plm(i,1,lnew)
             DPDmi(isup+2) = parity*Plm(i,3,lnew)
             DPDmi(isup+3) = parity*Plm(i,2,lnew)
          ENDDO

          IF (nmat.EQ.4) THEN
              DO i=4, nsup, 4
                 DPDpl(i) = DPDpl(i-3)
                 DPDmi(i) = DPDmi(i-3)
              ENDDO
          ENDIF

*---------------------------------------------------------------------
*         Construct the matrix D2 * S * D2            
*---------------------------------------------------------------------
          DSD(1,1) = cs(1,1,l)
          DSD(2,1) = 0.5D0*cs(1,2,l)
          DSD(2,2) = 0.25D0*(cs(2,2,l)+cs(3,3,l))
          DSD(3,2) = 0.25D0*(cs(2,2,l)-cs(3,3,l))
          DSD(3,1) = DSD(2,1)
          DSD(1,2) = DSD(2,1)
          DSD(1,3) = DSD(2,1)
          DSD(2,3) = DSD(3,2)
          DSD(3,3) = DSD(2,2)
          IF (nmat.EQ.4) THEN
              DSD(1,4)= 0.D0
              DSD(2,4)= 0.5D0*cs(3,4,l)
              DSD(3,4)= -DSD(2,4)
              DSD(4,4)= cs(4,4,l)
              DSD(4,1)= 0.D0
              DSD(4,2)= -DSD(2,4)
              DSD(4,3)= -DSD(3,4)
          ENDIF

*---------------------------------------------------------------------
*         Add a new term to the sum in Eq. (66)                 
*         The factor (-1)**m is cancelled by i**m in the Plm
*---------------------------------------------------------------------
          DO k2=1,nmat
             DO k1=1,nmat
                DO j=k2,nsup, nmat
                   SPj = DSD(k1,k2)*DPDpl(j)
                   DO i=k1,nsup, nmat
                      Zmplus(i,j)=Zmplus(i,j)+DPDpl(i)*SPj
                      Zmmin(i,j) =Zmmin(i,j)+DPDmi(i)*SPj
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

*---------------------------------------------------------------------
*         When last coefficient has been treated : 
*         skip recurrence  
*---------------------------------------------------------------------
          IF (l.EQ.ncoefs(layer) ) GOTO 1200

*---------------------------------------------------------------------
*         Do one step in recurrence for Plm0, Eq. (81)     
*---------------------------------------------------------------------
          twol1 = 2.D0*l+1.D0
          f1new = twol1/sqlm(l+1)
          f1old = sqlm(l)/sqlm(l+1)
          DO i=1,nmu
             u  = xmu(i)
             Plm(i,1,lold) = f1new*u*Plm(i,1,lnew)- 
     .                    f1old*Plm(i,1,lold)
          ENDDO

*---------------------------------------------------------------------
*         Do one step in recurrence for Plm2 and Plm-2, 
*         Eq. (82)  
*         only when they have been initialized:l >= max(m,2)
*---------------------------------------------------------------------
          IF (l.GE.max0(m,2)) THEN
              tmp   = 1.D0/(dble(l)*sql4(l+1)*sqlm(l+1))
              f2new = twol1*dble(l)*dble(l+1)*tmp
              f2newa= twol1*dble(2*m)*tmp
              f2old = dble(l+1)*sql4(l)*sqlm(l)*tmp
              DO i=1,nmu
                 u = xmu(i)
                 Plm(i,2,lold)=(f2new*u+f2newa)*
     .               Plm(i,2,lnew) - f2old*Plm(i,2,lold)
                 Plm(i,3,lold) = (f2new*u-f2newa)*
     .               Plm(i,3,lnew) - f2old*Plm(i,3,lold)
              ENDDO
          ENDIF
          itmp = lnew
          lnew = lold
          lold = itmp
1200  CONTINUE

*---------------------------------------------------------------------
*     End of summation loop over l                 
*     Calculate D1 * sum * D1                  
*---------------------------------------------------------------------
      CALL transf(Zmmin,nsup,nmat)
      CALL transf(Zmplus,nsup,nmat)

*---------------------------------------------------------------------
999   RETURN
      END
