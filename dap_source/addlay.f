      SUBROUTINE addlay(nmat,nmu,ebtop,ebbot,iad,
     .                  Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot)

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


C      DOUBLE PRECISION ebtop, ebbot, Rmtop, Tmtop, Rmbot, 
C     .                 Tmbot, Rmsbot
      DIMENSION W(nsupMAX,nsupMAX),X(nsupMAX,nsupMAX),
     .          Y(nsupMAX,nsupMAX),Z(nsupMAX,nsupMAX),
     .          Rmtop(nsupMAX,nsupMAX), Tmtop(nsupMAX,nsupMAX),
     .          Rmbot(nsupMAX,nsupMAX), Tmbot(nsupMAX,nsupMAX),
     .          Rmsbot(nsupMAX,nsupMAX),
     .          Etop(nsupMAX),Ebot(nsupMAX),
     .          ebtop(nmuMAX),ebbot(nmuMAX)

      LOGICAL verbo
      verbo = .false.

Cf2py intent(out) Rmbot,Tmbot,Rmsbot

*----------------------------------------------------------------------
      nmum= nmu-1
      nsup= nmat*nmu

      DO i=1,nsup
         Etop(i)= 0.D0
         Ebot(i)= 0.D0
      ENDDO

*----------------------------------------------------------------------
*     Handle the case of no scattering in the top layer:
*----------------------------------------------------------------------
      IF (iad.EQ.2) THEN
         CALL notop(Rmbot,Tmbot,ebtop,ebbot,nmu,nmat)
         GOTO 999

*----------------------------------------------------------------------
*     Handle the case of no scattering in the bottom layer:
*----------------------------------------------------------------------
      ELSEIF (iad.EQ.3) THEN
         CALL nobot(Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot,
     .              ebtop,ebbot,nmu,nmat)
         GOTO 999
      ENDIF

*----------------------------------------------------------------------
*     Use the product method to calculate the repeated reflections     
*     between the layers, see de Haan et al. (1987) Eqs. (111)-(115)  
*----------------------------------------------------------------------
*     X= R'*
      CALL star(X,Rmtop,nmat,nmu)

*     Y= R'*R"= C1
      CALL prod(Y,X,Rmbot,nmat,nmu,nmum)

*     Z= C1= S1
      CALL assign(Z,Y,nmat,nmu)

      ir= 0
100      ir= ir+1

*        X= Cr Cr= Cr+1
         CALL prod(X,Y,Y,nmat,nmu,nmum)

*        Y= Cr+1
         CALL assign(Y,X,nmat,nmu)

*        X= Sr Cr+1
         CALL prod(X,Z,Y,nmat,nmu,nmum)

         CALL addsm(Z,Z,X,nmat,nmu)
*                                                    Z = Sr + Sr Cr+1
         CALL addsm(Z,Z,Y,nmat,nmu)
*                                         Z = Sr + Sr Cr+1 + Cr+1 = Sr+1
         CALL trace(Y,trC,nmat,nmum)
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
            Etop((i-1)*nmat+k) = ebtop(i)
            Ebot((i-1)*nmat+k) = ebbot(i)
         ENDDO
      ENDDO

*     X = QT'
      CALL prod(X,Z,Tmtop,nmat,nmu,nmum)

*     W = QE'
      CALL rdiapr(W,Z,Etop,nmat,nmu)

*     X = QT' + QE'
      CALL addsm(X,X,W,nmat,nmu)

*     X = T' + QT' +QE' = D
      CALL addsm(X,X,Tmtop,nmat,nmu)

*----------------------------------------------------------------------
*     W = R"D
      CALL prod(W,Rmbot,X,nmat,nmu,nmum)

*     Y = R"E'
      CALL rdiapr(Y,Rmbot,Etop,nmat,nmu)

*     W = R"E' + R"D' = U
      CALL addsm(W,W,Y,nmat,nmu)

*----------------------------------------------------------------------
*     Y = T"D
      CALL prod(Y,Tmbot,X,nmat,nmu,nmum)

*     X = E"D
      CALL ldiapr(X,Ebot,X,nmat,nmu)

*     Y = T"D + E"D
      CALL addsm(Y,Y,X,nmat,nmu)

*     X = T"E'
      CALL rdiapr(X,Tmbot,Etop,nmat,nmu)

*     Y = T"D + E"D + T"E' =Ttot
      CALL addsm(Y,Y,X,nmat,nmu)

*----------------------------------------------------------------------
*     Tmtop = T'*
      CALL star(Tmtop,Tmtop,nmat,nmu)

*     X = T'*U
      CALL prod(X,Tmtop,W,nmat,nmu,nmum)

*     Rmbot = R' + T'*U
      CALL addsm(Rmbot,Rmtop,X,nmat,nmu)

*     W = E'U
      CALL ldiapr(W,Etop,W,nmat,nmu)

*     Rmbot = R' + T'*U + E'U = Rtot
      CALL addsm(Rmbot,Rmbot,W,nmat,nmu)

*----------------------------------------------------------------------
*     Rmtop= R'*
      CALL star(Rmtop,Rmtop,nmat,nmu)

*     X= QR'*
      CALL prod(X,Z,Rmtop,nmat,nmu,nmum)

*     W= T"*
      CALL tstar(W,Tmbot,nmat,nmu)

*     Z= QR'*T"*
      CALL prod(Z,X,Tmbot,nmat,nmu,nmum)

*     Z= QR'*T"*
      CALL prod(Z,X,W,nmat,nmu,nmum)

*     X= QR'*E"
      CALL rdiapr(X,X,Ebot,nmat,nmu)

*     X= QR'*E" + QR'*T"*
      CALL addsm(X,X,Z,nmat,nmu)

*     W= T"*
      CALL tstar(W,Tmbot,nmat,nmu)

*     Z= R'*T"*
      CALL prod(Z,Rmtop,W,nmat,nmu,nmum)

*     X= QR'*E" + QR'*T"* + R'*T"*
      CALL addsm(X,X,Z,nmat,nmu)

*     Z= R'*E"
      CALL rdiapr(Z,Rmtop,Ebot,nmat,nmu)

*     X= QR'*E" + QR'*T"* + R'*T"* + R'*E" = U*
      CALL addsm(X,X,Z,nmat,nmu)

*----------------------------------------------------------------------
*     Z= T"U*
      CALL prod(Z,Tmbot,X,nmat,nmu,nmum)

*     Rmsbot= R"* + T"U*
      CALL addsm(Rmsbot,Rmsbot,Z,nmat,nmu)

*     X= E"U*
      CALL ldiapr(X,Ebot,X,nmat,nmu)

*     Rmsbot= R"* + T"U* + E"U*= Rtot*
      CALL addsm(Rmsbot,Rmsbot,X,nmat,nmu)

*     Tmbot= T"D + E"D + T"E'= Ttot
      CALL assign(Tmbot,Y,nmat,nmu)

*----------------------------------------------------------------------
*     Direct transmission:
*       exp(-b/mu)= exp(-bbot/mu)*exp(-btop/mu)
*----------------------------------------------------------------------
      DO i=1,nmu
         ebbot(i)= ebbot(i)*ebtop(i)
      ENDDO

*----------------------------------------------------------------------
  999 RETURN
      END
