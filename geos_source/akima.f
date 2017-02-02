!*   Akima spline interpolation                         *
!*       Translated to FORTRAN by,                      *
!*	G. Mahapatra, 2017 TU Delft                     *
!* Original code:                                       *
!* http://www.lfd.uci.edu/~gohlke/code/akima.py.html    *   
!* Variable description:                                *
!* xx=Point/position of point to be interpolated        *
!* yy= Output value                                     *  
!* Xin= Input X array                                   *   
!* Yin= Array corresponding to the function X           * 
!* sz= size of the Input array                          *
!* nADD= No. of fourier elements to be added            *
!********************************************************
	SUBROUTINE akima(Xin,Yin,sz,xx,yy,nMAX,nADD)

        IMPLICIT NONE  

  	integer i,n,sz,diff,pos,ids(sz),j,idssz,bins,bb,nMAX,nADD

	DOUBLE PRECISION  Xin(nMAX),Yin(nMAX),X1(sz+nADD),Y1(sz+nADD),
     .                  xx,yy,mm,mmm,mp,mpp,v1,wj,dm(sz+2+nADD),
     .                  X(sz+nADD),Y(sz+nADD),m1diff(sz+2+nADD),
     .                  XM(sz+3+nADD),m(sz-1+nADD),b(sz+nADD),
     .                  Z(sz+nADD),m1(sz+3+nADD),f1(sz+nADD),
     .                  dx(sz-1+nADD),dy(sz-1+nADD),bids(sz+nADD),
     .                  c(sz-1+nADD),d(sz-1+nADD),f12(sz+nADD),
     .                  f2(sz+nADD)

C Assigning the values into smaller sized array	
	DO i = 1,sz
	   X1(i) = Xin(i)
	   Y1(i) = Yin(i)
	ENDDO

C       Add 6 extra elements. 3 to beginning. 3 to end.
	CALL addfouel(X1,Y1,sz,3,X,Y)
	
	n = SIZE(X)
	
	CALL difference(X,n,dx)
	CALL difference(Y,n,dy)	

	DO i = 1,n-1
		m(i)=dy(i)/dx(i)
	ENDDO

	mm = 2.d0*m(1)-m(2)
	mmm = 2.d0*mm-m(1)
	mp = 2.d0*m(n-1) -m(n-2)
	mpp = 2.d0*mp-m(n-1)

C       Forming the m1 array
	m1(1) = mmm
	m1(2) = mm	
	DO i = 1,n-1
		m1(2+i) = m(i)
	ENDDO
	m1(2+sz-1+1) = mp
	m1(2+sz-1+1+1) = mpp

	CALL difference(m1,n+3,m1diff)
	DO i = 1,n+2
	   dm(i) = DABS(m1diff(i))
	ENDDO

	f1 = dm(3:n+2)
	f2 = dm(1:n)
	f12 = f1+f2

	b = m1(2:n+1)

	DO i = 1,n
	   j=i
	   b(j)=(f1(j) * m1(j + 1) + f2(j) * m1(j + 2))
	   if (b(j).eq.0.d0) then
	      b(j) = 0.d0
	   else
	      b(j) = b(j)/f12(j)
	   end if
	ENDDO

	DO i = 1,n-1
	   c(i)=(3.0*m(i)-2.0*b(i)-b(i+1))/dx(i)
	   d(i) = (b(i) + b(1+i) - 2.0 * m(i)) / (dx(i) ** 2)
	ENDDO

	CALL digitize(xx,X,n,bins)
	bb = bins
	wj = xx - X(bb)
	yy = ((wj * d(bb) + c(bb)) * wj + b(bb)) * wj + Y(bb)

	return
	END

C---------------------------------------------------------------------
C  Calculates difference for a given array
	SUBROUTINE difference(X,szAr,Xdiff)

        IMPLICIT NONE  

	integer i,szAr
	DOUBLE PRECISION Xdiff(szAr-1), X(szAr)

	DO i = 1,szAr-1
	   Xdiff(i) = X(i+1) - X(i)
	ENDDO
	
	return
	END
C---------------------------------------------------------------------
C Binning subroutine which is equal to python digitize() function
C       The array containing values is "a". The array which defines
C	the bin is "b".
C       Returns the indices of the values in the bin corr. to "a".

	SUBROUTINE digitize(a,X,sz,i)

	IMPLICIT NONE

	integer c,sz,i,szAR
	DOUBLE PRECISION a,X(sz)

	if (a.ge.X(1) .and. a.le.X(sz)) then
	   
!Find relevant table interval
	   i=1
 300	   i=i+1
	   if (a > X(i)) goto 300
	   i=i-1
	else
	   write(*,*) 'Value out of bounds'
	   STOP
	end if

	return
	END
C---------------------------------------------------------------------
!*              Add fourier elements: addfouel                       *!
!*  This subroutine adds fourier elements to the left and right      *!
!*  of the fourier vector before it is passed into the interpolation *!
!*  scheme. The elements are mirrored images corresponding to the    *!
!*  position of the number in the array.                             *!
!*  Variables: 							     *!
!*	X, Y: Input arrays to be extended			     *!
!*	sz: Size of the input array				     *!
!*	n:  No. of elements to be added to both sides		     *!
!*	Xa, Ya: New arrays               			     *!
!*								     *!
!*  Author: Gourav Mahapatra TU Delft 2017                           *!

	SUBROUTINE addfouel(X,Y,sz,n,Xa,Ya)
        
	integer sz,n,c,i,d 
	DOUBLE PRECISION X(0:sz), Y(0:sz),
     .                  Xa(sz+2*n), Ya(sz+2*n)

	c = 1+n

! Assign the values of X,Y to middle of Xa,Ya
        do i = 0,sz
		Xa(c) = X(i)
		Ya(c) = Y(i)
		c = c+1
	end do 

! Add the elements to the start and end of Xa & Ya

	c=1+n 	!Counter
	d=0 	!Counter
	do i=1,n	
		Xa(c-1-d) = -Xa(c+1+d) !The points are mirrored with changed signs to the beginning
		Ya(c-1-d) = Ya(c+1+d)
        	Xa(sz+n+1+d) = Xa(sz+n) + (Xa(sz+n) - Xa(sz+n-1-d)) ! Mirroring points and preserving distance 
		Ya(sz+n+1+d) = Ya(sz+n-1-d) !
		d=d+1
	end do
c	write(*,*) 'writing from addfouel'     
	     
	return
	end

	
	

	
