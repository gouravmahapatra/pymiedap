* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE rdfousascii(foufile,rfou,nfou,nmat,nmugs,xmu)

Cf2py intent(in) foufile
Cf2py intent(out) nfou,nmat,nmugs,xmu
Cf2py intent(in,out) rfou

*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iunf 
      PARAMETER (iunf=23)

      INTEGER m,i1,i2,i,j,ki,ibase,nfou,nmat,nmugs

      DOUBLE PRECISION xmu(nmuMAX),rfou(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX)

      CHARACTER ch*1,foufile*200

*----------------------------------------------------------------------------
*     Get the input file with the Fourier coefficients:
*----------------------------------------------------------------------------
C      WRITE(*,*) 'Give the name of the Fourier coefficients file',
C     .           ' (max. 20):'

*----------------------------------------------------------------------------
*     Open the Fourier coefficients file:
*----------------------------------------------------------------------------
      OPEN(unit=iunf,file=foufile,status='old',err=999)

*----------------------------------------------------------------------------
*     Read the header (the dap.in file contents):
*----------------------------------------------------------------------------
2     READ(iunf,'(a1)',err=997,end=998) ch
      IF (ch.EQ.'#') GOTO 2
      BACKSPACE(iunf)

*----------------------------------------------------------------------------
*     Read the accuracy, the matrix size, the number of Gauss points and the 
*     Gaussian integration points:
*----------------------------------------------------------------------------
      READ(iunf,*) nmat
      READ(iunf,*) nmugs
      DO i=1,nmugs
         READ(iunf,*) xmu(i)
      ENDDO

      nfou=0
20    DO i=1,nmugs
         ibase= (i-1)*nmat
         DO j=1,nmugs
            READ(iunf,*,END=21) m,i1,i2,
     .                    (rfou(ibase+ki,j,nfou),ki=1,nmat)
         ENDDO
      ENDDO
      nfou= nfou+1
      GOTO 20
21    CLOSE(iunf)

      GOTO 1000

*----------------------------------------------------------------------------
999   WRITE(*,*) 'Error opening Fourier file!'
      STOP
998   WRITE(*,*) 'Error: unexpected end of the Fourier file!'
      STOP
997   WRITE(*,*) 'Error reading the Fourier file!'
      STOP

*----------------------------------------------------------------------------
1000  RETURN
      END
