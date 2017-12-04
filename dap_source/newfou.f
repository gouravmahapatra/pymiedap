* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE newfou(m,Rmbot,smf,iunfou,nmat,nmu,M0,nlays,nextm)

*----------------------------------------------------------------------
*     Write the Fourier-coefficients to file:
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iunfou,nmat,nmu,ib,jb,i,j,m,nlays,maxM

      INTEGER M0(nlaysMAX)

      DOUBLE PRECISION w,t1,t2,t3,t4,tot,tot_eps,
     .       Rmbot(nsupMAX,nsupMAX),smf(nmuMAX)

      LOGICAL nextm
      nextm= .true.

Cf2py intent(in) Rmbot, smf, m, iunfou, nmat, nmu, M0, nlays
Cf2py intent(out) nextm

*----------------------------------------------------------------------
*     Write the supermatrix elements to file: 
*
*     Note: the supermatrix factors are removed before writing!
*-----------------------------------------------------------------------
      tot= 0.D0
      tot_eps= 0.D0

      DO i=1,nmu
         ib= (i-1)*nmat
         DO j=1,nmu
            jb= (j-1)*nmat
            w=1.D0/(smf(i)*smf(j))

            IF (nmat.EQ.1) THEN
               t1= w*Rmbot(ib+1,jb+1)
               IF (DABS(t1).LT.eps) t1= 0.D0

               WRITE(iunfou,20) m,i,j,t1
            ELSEIF (nmat.EQ.3) THEN
               t1= w*Rmbot(ib+1,jb+1)
               IF (DABS(t1).LT.eps) t1= 0.D0
               t2= w*Rmbot(ib+2,jb+1)
               IF (DABS(t2).LT.eps) t2= 0.D0
               t3= w*Rmbot(ib+3,jb+1)
               IF (DABS(t3).LT.eps) t3= 0.D0

               WRITE(iunfou,23) m,i,j,t1,t2,t3
            ELSE
               t1= w*Rmbot(ib+1,jb+1)
               IF (DABS(t1).LT.eps) t1= 0.D0
               t2= w*Rmbot(ib+2,jb+1)
               IF (DABS(t2).LT.eps) t2= 0.D0
               t3= w*Rmbot(ib+3,jb+1)
               IF (DABS(t3).LT.eps) t3= 0.D0
               t4= w*Rmbot(ib+4,jb+1)
               IF (DABS(t4).LT.eps) t4= 0.D0

               WRITE(iunfou,24) m,i,j,t1,t2,t3,t4
            ENDIF

            tot= tot + DABS(w*Rmbot(ib+1,jb+1))
            tot_eps= tot_eps + eps

         ENDDO
      ENDDO

*----------------------------------------------------------------------
*     Test whether any of the coefficients is larger than zero:
*----------------------------------------------------------------------
      IF (tot.LT.tot_eps) THEN
         nextm= .false.
         GOTO 999
      ENDIF

*----------------------------------------------------------------------
*     Test whether nfouMAX has been reached:
*----------------------------------------------------------------------
      IF (m.GE.nfouMAX) THEN
         nextm= .false.
         GOTO 999
      ENDIF

*----------------------------------------------------------------------
*     Else, we must sum the Fourier series all the way to M0:
*----------------------------------------------------------------------
      maxM=0
      DO i=1,nlays
         IF (maxM.LT.M0(i)) maxM= M0(i)
      ENDDO

      IF (m.GE.maxM) THEN
         nextm= .false.
         GOTO 999
      ENDIF

*-----------------------------------------------------------------------
20    FORMAT(I4,2X,2(I3,2X),E16.8)
23    FORMAT(I4,2X,2(I3,2X),3(E16.8,1X))
24    FORMAT(I4,2X,2(I3,2X),4(E16.8,1X))

*-----------------------------------------------------------------------
999   RETURN
      END
