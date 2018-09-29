* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE wrout(alpha,nmat,Sv)

*----------------------------------------------------------------------------
*     Write the output to the output file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER nmat,ki

      INTEGER iunf
      PARAMETER (iunf=23)

      DOUBLE PRECISION pol,alpha,Sv(nmatMAX) 

      CHARACTER outfile*16

*----------------------------------------------------------------------------
*     Open the output file:
*----------------------------------------------------------------------------
      OPEN(unit=iunf,file='geos.out')
      WRITE(iunf,800)
      WRITE(iunf,801)

      IF (nmat.EQ.1) THEN
         WRITE(iunf,802)
      ELSEIF (nmat.EQ.3) THEN
         WRITE(iunf,803)
      ELSEIF (nmat.EQ.4) THEN
         WRITE(iunf,804)
      ENDIF

*----------------------------------------------------------------------------
*     Write the elements of the Stokes vector to the output file:
*----------------------------------------------------------------------------
      IF (nmat.EQ.1) THEN
         WRITE(iunf,900) alpha,Sv(1)
      ELSEIF (nmat.EQ.3) THEN
         pol= DSQRT(Sv(2)*Sv(2)+Sv(3)*Sv(3))/Sv(1)
         WRITE(iunf,901) alpha,(Sv(ki),ki=1,nmat),pol
      ELSEIF (nmat.EQ.4) THEN
         pol= DSQRT(Sv(2)*Sv(2)+Sv(3)*Sv(3)+Sv(4)*Sv(4))/Sv(1)
         WRITE(iunf,902) alpha,(Sv(ki),ki=1,nmat),pol
      ENDIF

*----------------------------------------------------------------------------
*     Close the output file:
*----------------------------------------------------------------------------
      CLOSE(iunf)

*----------------------------------------------------------------------------
800   FORMAT('# Output of program geos')
801   FORMAT('#')
802   FORMAT('# alpha       I')
803   FORMAT('# alpha       I',
     .       '                Q                U               P')
804   FORMAT('# alpha       I',
     .       '                Q                U',
     .       '                V               P')
900   FORMAT(F9.3,3X,E16.8)
901   FORMAT(F9.3,3X,3(E16.8,1X),2X,F8.6)
902   FORMAT(F9.3,3X,4(E16.8,1X),2X,F8.6)

*----------------------------------------------------------------------------
      RETURN
      END
