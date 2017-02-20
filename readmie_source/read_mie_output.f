      SUBROUTINE readmieoutput(namein,nameout,thetas,Fs)

****************************************************************
* PURPOSE: 
* In this program the elements of the scattering matrix are
* calculated as a function of the scattering angle
*
* input: aerosol type
*        wavelength 
*
* output: 
*        as function of theta: F11, F21, F22, F33, F34, F44, 
*        and POL (degree of polarisation) 
*
* DATE:
* September 28, 1995
* 2015-03-16 L. Rossi
*
* AUTHOR:
* D. M. Stam, L. Rossi
****************************************************************
      IMPLICIT NONE

      INTEGER nstep,n,NDcoef, ncoef
      PARAMETER (NDcoef=3000)

C      DOUBLE PRECISION theta,wavel,comp,coefs,F,POL,SSA
      DOUBLE PRECISION theta,coefs,F,POL,SSA
      DIMENSION coefs(4,4,0:NDcoef),F(6)

      DOUBLE PRECISION thetas, Fs
      DIMENSION thetas(901), Fs(6,901)

      DOUBLE PRECISION pi
      PARAMETER (pi=3.1415926535898D0)

C      CHARACTER*8 dummy,namewavel
      CHARACTER*20 nameout,namein
C      CHARACTER*18 name
      CHARACTER*23 title

Cf2py intent(in) namein, nameout
Cf2py intent(out) thetas, Fs

****************************************************************
* 1.  Get the input:
****************************************************************
C      WRITE(*,*) 
C      WRITE(*,*) 'PROGRAM aer'
C      WRITE(*,*) 
C      WRITE(*,*) '  Give the name of the input file:'
C      READ(*,'(A20)') namein
C      WRITE(*,*) '  Give the name of the output file:'
C      READ(*,'(A20)') nameout

****************************************************************
* 3.  Open the output file and write the input away:
****************************************************************
      OPEN(unit=4,file=nameout)

300   FORMAT('# OUTPUT aer_p.f')
301   FORMAT('#')
302   FORMAT('# aerosol specie       :',2X,A1)
312   FORMAT('# wavelength (nm):',2X,F8.3)
314   FORMAT('# last wavelength (nm) :',2X,F8.3)
316   FORMAT('#  theta       F11           F22',
     .       '           F33           F44',
     .       '           F12           F34        POL')

****************************************************************
* 7.a Read the expansion coefficients file:
****************************************************************
      OPEN (unit=10,file=namein,status='old')
         READ (10,'(A23)') title
         IF (title.NE.' EXPANSION COEFFICIENTS') THEN
         STOP 'wrong mie.sc file in readmiesc.f!'
      ENDIF
      CALL readsc(10,coefs,NDcoef,ncoef,SSA)
      CLOSE(10)

****************************************************************
* 4.  Loop over the scattering angle:
****************************************************************
*      nstep= 181
      nstep= 901
*      theta=-1.0D0
      theta=-0.2D0
      DO 410 n=1,nstep
*         theta= theta + 1.D0
         thetas(n) = theta
         theta= theta + 0.2D0

****************************************************************
* 7.b Evaluate the expansion in GSF:
****************************************************************
         CALL expand(ncoef,coefs,theta,F)
            
****************************************************************
* 7.c Write the results to the output file:
****************************************************************
         POL= -F(5)/F(1)
         Fs(1,n) = F(1)
         Fs(2,n) = F(2)
         Fs(3,n) = F(3)
         Fs(4,n) = F(4)
         Fs(5,n) = F(5)
         Fs(6,n) = F(6)
         WRITE(4,700) theta,F(1),F(2),F(3),F(4),F(5),F(6),POL

700   FORMAT(F10.3,2X,6(2X,F14.8),2X,F8.6)

****************************************************************
410   ENDDO

****************************************************************
      CLOSE(4)
*      RETURN
      END

****************************************************************
      SUBROUTINE readsc(iunit,alfbet,ndl,lmax,SSA)

****************************************************************
*  This routine is an adapted version of subroutine readsc in 
*  the adding program gap, written by V.L. Dolman at the Free
*  University, Amsterdam.  
****************************************************************
*  On entry :                                                
*      iunit    unit number of the coefficient file      
*      ndl      third dimension of array alfbet (ndl >= lmax) 
*  On exit in case of success :                  
*      The coefficients and other information is read from the 
*      file with unit number iunit in the standard format.  
*      alfbet   array containing the coefficients, according to: 
*                                                             
*                 ( alpha1(l) beta2(l)     0         0      )  
*  alfbet(i,j,l)= ( beta2(l)  alpha2(l)    0         0      ) 
*                 (   0         0       alpha3(l) beta2(l)  )
*                 (   0         0      -beta2(l)  alpha4(l) ) 
*                                                  
*      lmax     maximum index l for which coefficients are 
*               non-zero 
*      eof      .false.                                          
*  On exit in case no block of coefficients was found on the 
*               file:     *
*      eof       .true.                                
*      all other variables are unchanged.  
****************************************************************
      IMPLICIT NONE

      INTEGER iunit,ndl,lmax,i,l,j

      DOUBLE PRECISION alfbet,a1,a2,a3,a4,b1,b2,SSA
      DIMENSION  alfbet(4,4,0:ndl)
      LOGICAL eof

****************************************************************
      eof=.false.

      DO 10 i=1,6
        READ(iunit,*,end=999)  
10    CONTINUE
      READ(iunit,'(32X,E25.14)') SSA
      DO 20 i=1,6
        READ(iunit,*,end=999)
20    CONTINUE

****************************************************************
*  Set the coefficients to zero :            
****************************************************************
      DO 100 l=0,ndl
         DO 100 j=1,4
            DO 100 i=1,4
               alfbet(i,j,l)=0.0
100   CONTINUE

****************************************************************
*  Loop to read in the lines containing the coefficients until
*  the termination line with l=-1 is reached.         
****************************************************************
      lmax=-1
200   READ(iunit,4,end=999) l,a1,a2,a3,a4,b1,b2
      IF (l .gt. lmax) lmax=l
      IF ((l .gt. ndl) .or. (l .lt. -1)) THEN
         print *,' readsc: illegal index read l = ',l
         print *,'         minimum is -1 maximum is ndl = ',ndl
         stop 'in readsc illegal index read from coefficient file'
      ENDIF
      IF (l .gt. -1) THEN
         alfbet(1,1,l)=a1
         alfbet(2,2,l)=a2
         alfbet(3,3,l)=a3
         alfbet(4,4,l)=a4
         alfbet(1,2,l)=b1
         alfbet(2,1,l)=b1
         alfbet(3,4,l)=b2
         alfbet(4,3,l)=-b2
         GOTO 200
      ENDIF

****************************************************************
4     FORMAT(i4,6f19.14)
      RETURN
998   eof=.true.
      RETURN
999   print *,' readsc: unexpected end of file on unit ',iunit
      STOP 'in readsc unexpected end of file encountered'
      END



****************************************************************
      SUBROUTINE expand(ncoef,coefs,scangle,F)

****************************************************************
*  This routine is an adapted version of subroutine expand in
*  the Meerhoff Mie Program 1.0, written by V.L. Dolman at the   
*  Free University, Amsterdam.              
****************************************************************
C LR : it would be good to get rid of that IMPLICIT
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      INTEGER k, lnew, lold, l
      DOUBLE PRECISION u, qroot6, fac1, fac2, fac3
      DOUBLE PRECISION sql4, sql41, twol1, tmp1, tmp2, denom
      DOUBLE PRECISION itmp, F2i, F3i

      PARAMETER (pi=3.141592653589793238462643d0,NDcoef=1000)
      DIMENSION F(6),coefs(4,4,0:NDcoef),
     .          P00(2),P02(2),P22(2),P2m2(2)

****************************************************************
*  Initialize:
****************************************************************
      u= COS(scangle*pi/180.D0)
      qroot6 = -0.25d0*dsqrt(6.d0)

****************************************************************
*  Set scattering matrix F to zero: 
****************************************************************
      DO k=1,6
         F(k)= 0.d0
      ENDDO

****************************************************************
*  Start loop over the coefficient index l     
*  first update generalized spherical functions, then calculate
*  coefs. lold and lnew are pointer-like indices used in 
*  recurrence  
****************************************************************
      lnew= 1
      lold= 2
      DO 90 l=0, ncoef
         IF (l.eq.0) THEN
****************************************************************
*           Adding paper Eq. (77) with m=0            
****************************************************************
            P00(lold) = 1.d0
            P00(lnew) = 0.d0
            P02(lold) = 0.d0
            P22(lold) = 0.d0
            P2m2(lold)= 0.d0
            P02(lnew) = 0.d0
            P22(lnew) = 0.d0
            P2m2(lnew)= 0.d0
         ELSE
            fac1 = (2.d0*dble(l)-1.d0)/dble(l)
            fac2 = (dble(l)-1.d0)/dble(l)
****************************************************************
*           Adding paper Eq. (81) with m=0            
****************************************************************
            P00(lold) = fac1*u*P00(lnew) - fac2*P00(lold)
         ENDIF
         IF (l.eq.2) THEN
****************************************************************
*           Adding paper Eq. (78)              
*           sql4 contains the factor sqrt((l+1)*(l+1)-4) 
*           needed in the recurrence Eqs. (81) and (82)  
****************************************************************
            P02(lold)  = qroot6*(1.d0-u*u)
            P22(lold)  = 0.25d0*(1.d0+u)*(1.d0+u)
            P2m2(lold) = 0.25d0*(1.d0-u)*(1.d0-u)
            P02(lnew)  = 0.d0
            P22(lnew)  = 0.d0
            P2m2(lnew) = 0.d0
            sql41 = 0.d0
         ELSEIF (l.gt.2) THEN
****************************************************************
*             Adding paper Eq. (82) with m=0 and m=2      
****************************************************************
            sql4  = sql41
            sql41 = dsqrt(dble(l*l)-4.d0)
            twol1 = 2.d0*dble(l)-1.d0
            tmp1  = twol1/sql41
            tmp2  = sql4/sql41
            denom = (dble(l)-1.d0)*(dble(l*l)-4.d0)
            fac1  = twol1*(dble(l)-1.d0)*dble(l)/denom
            fac2  = 4.d0*twol1/denom
            fac3  = dble(l)*
     .              ((dble(l)-1.d0)*(dble(l)-1.d0)-4.d0)/denom
            P02(lold) = tmp1*u*P02(lnew) - tmp2*P02(lold)
            P22(lold) =(fac1*u-fac2)*P22(lnew) - fac3*P22(lold)
            P2m2(lold)=(fac1*u+fac2)*P2m2(lnew)- fac3*P2m2(lold)
         ENDIF
****************************************************************
*    Switch indices so that lnew indicates the function with   
*    the present index value l, this mechanism prevents swapping
*    of entire arrays.                                         
****************************************************************
         itmp = lnew
         lnew = lold
         lold = itmp
****************************************************************
*        Now add  the l-th term to the scattering matrix.    
*        See adding paper Eqs. (68)-(73).                   
****************************************************************
         F(1) = F(1) + coefs(1,1,l)*P00(lnew)
         F(2) = F(2) + (coefs(2,2,l)+coefs(3,3,l))*P22(lnew)
         F(3) = F(3) + (coefs(2,2,l)-coefs(3,3,l))*P2m2(lnew)
         F(4) = F(4) + coefs(4,4,l)*P00(lnew)
         F(5) = F(5) + coefs(1,2,l)*P02(lnew)
         F(6) = F(6) + coefs(3,4,l)*P02(lnew)
90    CONTINUE
      F2i    = 0.5d0*(F(2)+F(3))
      F3i    = 0.5d0*(F(2)-F(3))
      F(2) = F2i
      F(3) = F3i
****************************************************************
*     End of loop over index l           
****************************************************************
      RETURN
      END
