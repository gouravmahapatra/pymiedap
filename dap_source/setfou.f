* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE setfou(coefs,ncoefs,nlays,a,b,xmu,nmug,
     .                  M0,M1,M2,iadd)

**********************************************************************
*  Calculate bounds M0, M1 and M2 on the Fourier index m such that
*  	for  0 <= m <= M2 doubling must be used,     
*  	for  M2 < m <= M1 first plus second order scattering suffice, 
*  	for  M1 < m <= M0 only first order scattering is needed.   
*  For  M0 < m there is no scattering at all !      
*
*  The criteria are described in de Haan et al. (1987) p. 38.
*
*  Also determine the options indicating for each layer which type
*  of adding procedure should be used :          
*     iadd = 1 : normal adding                   
*     iadd = 2 : top layer has no scattering        
*     iadd = 3 : bottom layer has no scattering      
*
*  This follows directly from the values of M0 for all layers.
*  The number of coefficients is truncated such that lmax <= M0,
*  see de Haan et al. (1987) first few lines of section 7.2.
**********************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,nmug,M0bot,m,nlays

      INTEGER ncoefs(nlaysMAX)

      DOUBLE PRECISION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlaysMAX),
     .                 xmu(nmuMAX),a(nlaysMAX),b(nlaysMAX)

      INTEGER M0(nlaysMAX),M1(nlaysMAX),M2(nlaysMAX),
     .        iadd(nlaysMAX,0:nfouMAX)

      DOUBLE PRECISION xmumin,qf,qff,qff2,fbmu,f3b,am,h

      LOGICAL verbo
      verbo = .FALSE.

Cf2py intent(in) coefs,ncoefs,nlays,a,b,xmu,nmug
Cf2py intent(out) M0,M1,M2,iadd

*---------------------------------------------------------------------
*     Determine the minimum mu-value to be used in criteria:
*     (array xmu contains the Gaussian abscissae)
*---------------------------------------------------------------------
      xmumin= 1.D0
      DO i=1,nmug
         IF (xmu(i).LT.xmumin) xmumin= xmu(i)
      ENDDO

*---------------------------------------------------------------------
*     Start the loop over the layers:
*---------------------------------------------------------------------
      DO i=1,nlays

*---------------------------------------------------------------------
*        There is NO scattering in the layer:
*---------------------------------------------------------------------
         IF (a(i).LT.1.D-10) THEN
            M0(i)= -1
            M1(i)= -1
            M2(i)= -1

*---------------------------------------------------------------------
*        There is scattering in the layer:
*---------------------------------------------------------------------
         ELSE
            fbmu= b(i)/xmumin
            IF (fbmu.GT.1.D0) fbmu= 1.D0
            f3b= 3.D0*b(i)
            IF (f3b.GT.1.D0) f3b= 1.D0
            qf  = 0.25D0*fbmu
            qff = 0.25D0*fbmu*f3b
            qff2= 0.25D0*fbmu*f3b**2

*---------------------------------------------------------------------
*           Start loop over Fourier index m to find M0 
*           (using the Eqs. below (140)): 
*---------------------------------------------------------------------
            m= ncoefs(i)

100         am= a(i)*coefs(1,1,m,i)/DBLE(2*m+1)
            h = qf*am*DBLE(2*m+1)
            IF ((h.LE.eps).AND.(m.GT.0)) THEN
               m= m-1
               GOTO 100
            ENDIF
            M0(i)= m

*---------------------------------------------------------------------
*           Continue loop over Fourier index m to find M1 
*           (using the Eqs. below (140)):
*---------------------------------------------------------------------
200         am= a(i)*coefs(1,1,m,i)/DBLE(2*m+1)
            h = qff*am**2*DBLE(2*m+1)
            IF ((h.LE.eps).AND.(m.GT.0)) THEN
               m= m-1
               GOTO 200
            ENDIF
            M1(i)= m

*---------------------------------------------------------------------
*           Continue loop over Fourier index m to find M2 
*           (using the Eqs. below (140)):
*---------------------------------------------------------------------
300         am= a(i)*coefs(1,1,m,i)/DBLE(2*m+1)
            h = qff2*am**3*DBLE(2*m+1)
            IF ((h.LE.eps).AND.(m.GT.0)) THEN
               m= m-1
               GOTO 300
            ENDIF
            M2(i)= m

         ENDIF

         ncoefs(i)= M0(i)

*---------------------------------------------------------------------
*        Check M0, M1, and M2:
*---------------------------------------------------------------------
         IF (verbo) print *,' setfou: layer',i,'  M0 =',M0(i)
     .                      ,'  M1 =',M1(i),'  M2 =',M2(i)

*---------------------------------------------------------------------
*     End of the loop over the layers:
*---------------------------------------------------------------------
      ENDDO 

*---------------------------------------------------------------------
*     Now determine the adding option iadd for each layer and each  
*     Fourier index m.                              
*---------------------------------------------------------------------
      DO m=0,nfouMAX
         M0bot= -1
         DO i=1,nlays
            IF ((m.GT.M0(i)) .OR. (a(i).LT.1.D-10)) THEN
               iadd(i,m)= 2
            ELSE IF (m.GT.M0bot) THEN
               iadd(i,m)= 3
            ELSE
               iadd(i,m)= 1
            ENDIF
            M0bot= MAX0(M0bot,M0(i))
         ENDDO
      ENDDO

*---------------------------------------------------------------------
*     To facilitate the reflection by the surface:
*---------------------------------------------------------------------
      IF (iadd(1,0).EQ.3) iadd(1,0)=1

**********************************************************************
      RETURN
      END
