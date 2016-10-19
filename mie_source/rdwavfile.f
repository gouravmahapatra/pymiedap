      SUBROUTINE rdwavfile(wavfile,wavels,nr,ni,nwavels)

************************************************************************
* PURPOSE: 
*  Open and read the wavelength data file
*
* DATE: December 2002
*
* AUTHOR: D. M. Stam
*
* Version for wavelength dependent refractive index
* January 2008
*
************************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,nwavels,iunwav,j
      PARAMETER (iunwav=54)

      DOUBLE PRECISION wavels(nwavelsMAX),a1,a2,a3,
     .                 nr(nwavelsMAX),ni(nwavelsMAX)

      CHARACTER wavfile*16,ch*1

*-----------------------------------------------------------------------
      DO i=1,nwavelsMAX
         wavels(i)= 0.D0
         nr(i)= 0.D0
         ni(i)= 0.D0
      ENDDO

*-----------------------------------------------------------------------
*     Open and read the wavelength data file:
*     wavels in microns!
*-----------------------------------------------------------------------
      OPEN (file=wavfile,unit=iunwav,status='old',err=999)

1     READ(iunwav,'(a1)',err=999,end=998) ch
      IF (ch.EQ.'#') GOTO 1
      BACKSPACE(iunwav)

      i=0
10    READ(iunwav,*,err=999,end=11) j,wavels(i+1),nr(i+1),ni(i+1)
      i=i+1
      IF (i.EQ.nwavelsMAX) THEN
         WRITE(*,*) 'rdwavfile: nwavels > nwavelsMAX'
         STOP
      ENDIF
      GOTO 10
11    nwavels= i

      CLOSE (iunwav)

*-----------------------------------------------------------------------
*     Write to standard output:
*-----------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,*) 'The total number of wavelengths is ',nwavels
      WRITE(*,*)

************************************************************************
      RETURN
998   WRITE(*,*) 'rdwavfile: unexpected end of file'
      STOP
999   WRITE(*,*) 'rdwavfile: ERROR'
      STOP
      END
 
