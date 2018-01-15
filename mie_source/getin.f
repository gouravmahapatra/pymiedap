* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE getin(wavfile,delta,cutoff,thmin,thmax,step,
     .                 nsubr,ngaur,idis,par1,par2,par3,specie)

*-----------------------------------------------------------------------
*     For the version with a wavelength dependent refractive index!
*     January 2008
*
*     Read the input file:
*-----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iuni
      PARAMETER (iuni=13)

      INTEGER i,idis,ngaur,nsubr

      DOUBLE PRECISION par1,par2,par3,delta,cutoff,thmin,thmax,
     .                 step

      CHARACTER specie*1,wavfile*16

*-----------------------------------------------------------------------
*     Open the input file mie.in and read the first lines:
*-----------------------------------------------------------------------
      OPEN(iuni,file='mie.in',err=998)

      DO i=1,4
         READ(iuni,*,err=999)
      ENDDO

*-----------------------------------------------------------------------
*     Read the GENERAL parameters:
*-----------------------------------------------------------------------
      READ(iuni,*,err=999) wavfile
      READ(iuni,*,err=999) delta
      READ(iuni,*,err=999) cutoff
      READ(iuni,*,err=999) thmin
      READ(iuni,*,err=999) thmax
      READ(iuni,*,err=999) step

      DO i=1,19
         READ(iuni,*,err=999)
      ENDDO

*-----------------------------------------------------------------------
*     Check the parameters:
*-----------------------------------------------------------------------
      IF (delta.LT.0.D0) THEN
         WRITE(*,*) 'getin: truncation value Mie sum < 0!'
         STOP
      ENDIF
      IF (cutoff.LT.0.D0) THEN
         WRITE(*,*) 'getin: cutoff value size distribution < 0!'
         STOP
      ENDIF
      IF ((thmin.LT.0.D0).OR.(thmin.GT.180.D0)) THEN
         WRITE(*,*) 'getin: illegal value minimum scattering angle!'
         STOP
      ENDIF
      IF ((thmax.LT.0.D0).OR.(thmax.GT.180.D0)) THEN
         WRITE(*,*) 'getin: illegal value maximum scattering angle!'
         STOP
      ENDIF
      IF (thmax.LT.thmin) THEN
         WRITE(*,*) 'getin: maximum scattering angle < minimum angle!'
         STOP
      ENDIF
      IF ((step.LT.0.D0).OR.(step.GT.180.D0)) THEN
         WRITE(*,*) 'getin: illegal value step scattering angle!'
         STOP
      ENDIF

*-----------------------------------------------------------------------
*     Read the AEROSOL parameters:
*-----------------------------------------------------------------------
      READ(iuni,*,err=999) specie
      READ(iuni,*,err=999) idis
      READ(iuni,*,err=999) nsubr
      READ(iuni,*,err=999) ngaur
      READ(iuni,*,err=999) par1
      READ(iuni,*,err=999) par2
      READ(iuni,*,err=999) par3

*-----------------------------------------------------------------------
*     Check the parameters:
*-----------------------------------------------------------------------
      IF ((idis.LT.0).OR.(idis.GT.8)) THEN
         WRITE(*,*) 'getin: illegal value size distribution!'
         STOP
      ENDIF
      IF (nsubr.LT.0) THEN
         WRITE(*,*) 'getin: negative number of subintervals!'
         STOP
      ENDIF
      IF (ngaur.LT.0) THEN
         WRITE(*,*) 'getin: negative number of Gaussian points!'
         STOP
      ENDIF

*-----------------------------------------------------------------------
*     Close the input file and return:
*-----------------------------------------------------------------------
      CLOSE(iuni)
      RETURN

*-----------------------------------------------------------------------
998   WRITE(*,*) 'getin: error opening mie.in!'
      STOP
999   WRITE(*,*) 'getin: error reading mie.in!'
      STOP

*-----------------------------------------------------------------------
      END
