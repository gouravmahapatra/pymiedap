* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE addlay(nmat,ebtop,ebbot,iad,Rmtop,Tmtop,
     .                  nmu,nsup,Rmbot,Tmbot,Rmsbot)

*----------------------------------------------------------------------
*  Use the adding method to calculate reflection and transmission  
*  of the combination of top and bottom layer. The top layer is     
*  assumed homogeneous. The resulting reflection and transmission    
*  supermatrices are returned through arrays Rmbot, Tmbot and Rmsbot. 
*
*  BEWARE:  Rmtop and Tmtop are used as scratch space !                
*
*  The option variable iad indicates the following cases : 
*     iad= 1 : normal adding                          
*     iad= 2 : no scattering in top layer  
*     iad= 3 : no scattering in bottom layer      
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER iad,nsup,nmat,nmu,nmum

      PARAMETER (maxrep=15,trmin=1.D-10)

      REAL*8, DIMENSION(nsup,nsup) :: W,X,Y,Z,Rmtop,Tmtop,Rmbot,
     .          Tmbot,Rmsbot,Etop, Ebot

      REAL*8, DIMENSION(nmu) :: ebtop, ebbot

      LOGICAL verbo
      verbo = .false.

Cf2py intent(out) Rmbot,Tmbot,Rmsbot

*----------------------------------------------------------------------
      nmum= nmu-1

      Etop=0.D0
      Ebot=0.D0

*----------------------------------------------------------------------
*     Handle the case of no scattering in the top layer:
*----------------------------------------------------------------------
      IF (iad.EQ.2) THEN
         CALL notop(Rmbot,Tmbot,ebtop,ebbot,nmat,nmu,nsup)
         GOTO 999

*----------------------------------------------------------------------
*     Handle the case of no scattering in the bottom layer:
*----------------------------------------------------------------------
      ELSEIF (iad.EQ.3) THEN
         CALL nobot(Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot,
     .              ebtop,ebbot,nmat,nmu,nsup)
         GOTO 999
      ENDIF

*----------------------------------------------------------------------
*     Use the product method to calculate the repeated reflections     
*     between the layers, see de Haan et al. (1987) Eqs. (111)-(115)  
*----------------------------------------------------------------------
*     X= R'*
      CALL star(X,Rmtop,nmat,nmu,nsup)

*     Y= R'*R"= C1
      CALL prod(Y,X,Rmbot,nmat,nmu,nmum)

*     Z= C1= S1
      Z=Y

      ir= 0
100      ir= ir+1

*        X= Cr Cr= Cr+1
         CALL prod(X,Y,Y,nmat,nmu,nmum)

*        Y= Cr+1
         Y=X

*        X= Sr Cr+1
         CALL prod(X,Z,Y,nmat,nmu,nmum)

         Z=Z+X
*                                                    Z = Sr + Sr Cr+1
         Z=Z+Y
*                                         Z = Sr + Sr Cr+1 + Cr+1 = Sr+1
         CALL trace(Y,nmat,nmum,nsup,trC)
*                                                    trC = trace(Cr+1)
         IF (verbo) PRINT *,' addlay: r = ',ir,' trace = ',trC

         IF ((trC.GT.trmin) .AND. (ir.LE.maxrep)) GOTO 100

*----------------------------------------------------------------------*
      IF (ir.GT.maxrep) THEN
         PRINT *,' addlay: WARNING repeated reflections did not'
     .           ,' converge after',maxrep,' steps'
         PRINT *,'         proceed anyway !'
      ENDIF

*----------------------------------------------------------------------
*     Now Z contains the matrix Q Eq. (115)                               
*     Use the adding Eqs. (85)-(91)                                      
*----------------------------------------------------------------------
      DO i=1,nmu
         DO k=1,nmat
            Etop((i-1)*nmat+k,(i-1)*nmat+k) = ebtop(i)
            Ebot((i-1)*nmat+k,(i-1)*nmat+k) = ebbot(i)
         ENDDO
      ENDDO

*     X = QT'
      CALL prod(X,Z,Tmtop,nmat,nmu,nmum)

*     W = QE'
      W=MATMUL(Z,Etop)

*     X = QT' + QE'
      X=X+W

*     X = T' + QT' +QE' = D
      X=X+Tmtop
*----------------------------------------------------------------------
*     W = R"D
      CALL prod(W,Rmbot,X,nmat,nmu,nmum)

*     Y = R"E'
      Y=MATMUL(Rmbot,Etop)

*     W = R"E' + R"D' = U
      W=W+Y
*----------------------------------------------------------------------
*     Y = T"D
      CALL prod(Y,Tmbot,X,nmat,nmu,nmum)

*     X = E"D
      X=MATMUL(Ebot,X)

*     Y = T"D + E"D
      Y=Y+X

*     X = T"E'
      X=MATMUL(Tmbot,Etop)

*     Y = T"D + E"D + T"E' =Ttot
      Y=Y+X

*----------------------------------------------------------------------
*     Tmtop = T'*
      CALL star(Tmtop,Tmtop,nmat,nmu,nsup)

*     X = T'*U
      CALL prod(X,Tmtop,W,nmat,nmu,nmum)

*     Rmbot = R' + T'*U
      Rmbot=Rmtop+X

*     W = E'U
      W=MATMUL(Etop,W)

*     Rmbot = R' + T'*U + E'U = Rtot
      Rmbot=Rmbot+W

*----------------------------------------------------------------------
*     Rmtop= R'*
      CALL star(Rmtop,Rmtop,nmat,nmu,nsup)

*     X= QR'*
      CALL prod(X,Z,Rmtop,nmat,nmu,nmum)

*     W= T"*
      CALL tstar(Tmbot,nmat,nsup,W)

*     Z= QR'*T"*
      CALL prod(Z,X,Tmbot,nmat,nmu,nmum)

*     Z= QR'*T"*
      CALL prod(Z,X,W,nmat,nmu,nmum)

*     X= QR'*E"
      X=MATMUL(X,Ebot)

*     X= QR'*E" + QR'*T"*
      X=X+Z

*     W= T"*
      CALL tstar(Tmbot,nmat,nsup,W)

*     Z= R'*T"*
      CALL prod(Z,Rmtop,W,nmat,nmu,nmum)

*     X= QR'*E" + QR'*T"* + R'*T"*
      X=X+Z

*     Z= R'*E"
      Z=MATMUL(Rmbot,Ebot)

*     X= QR'*E" + QR'*T"* + R'*T"* + R'*E" = U*
      X=X+Z

*----------------------------------------------------------------------
*     Z= T"U*
      CALL prod(Z,Tmbot,X,nmat,nmu,nmum)

*     Rmsbot= R"* + T"U*
      Rmsbot=Rmsbot+Z

*     X= E"U*
      X=MATMUL(Ebot,X)

*     Rmsbot= R"* + T"U* + E"U*= Rtot*
      Rmsbot=Rmsbot+X

*     Tmbot= T"D + E"D + T"E'= Ttot
      Tmbot=Y

*----------------------------------------------------------------------
*     Direct transmission:
*       exp(-b/mu)= exp(-bbot/mu)*exp(-btop/mu)
*----------------------------------------------------------------------
      ebbot=ebbot*ebtop

*----------------------------------------------------------------------
  999 RETURN
      END
