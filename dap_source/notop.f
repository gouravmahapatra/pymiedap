      SUBROUTINE notop(Rmbot,Tmbot,ebtop,ebbot,nmu,nmat)

*----------------------------------------------------------------------
*     Use the adding equations when there is no scattering in 
*     the top layer:
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,k,nmu,nmat

      DOUBLE PRECISION Rmbot(nsupMAX,nsupMAX),Tmbot(nsupMAX,nsupMAX),
     .                 ebtop(nmuMAX),ebbot(nmuMAX),
     .                 Etop(nsupMAX),Ebot(nsupMAX)

Cf2py intent(in,out) Rmbot,Tmbot,ebtop,ebbot

*----------------------------------------------------------------------
      DO i=1,nmu
         DO k=1,nmat
            Etop((i-1)*nmat+k) = ebtop(i)
            Ebot((i-1)*nmat+k) = ebbot(i)
         ENDDO
      ENDDO 

*----------------------------------------------------------------------
*     Reflection  R= E'R"E' 
*----------------------------------------------------------------------
      CALL rdiapr(Rmbot,Rmbot,Etop,nmat,nmu)
      CALL ldiapr(Rmbot,Etop,Rmbot,nmat,nmu)

*----------------------------------------------------------------------
*     Transmission T= T"E'  
*----------------------------------------------------------------------
      CALL rdiapr(Tmbot,Tmbot,Etop,nmat,nmu)

*----------------------------------------------------------------------
*     Reflection star R* = R"* is trivial !                               
*----------------------------------------------------------------------
*     Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu)        
*----------------------------------------------------------------------
      DO i=1,nmu
         ebbot(i)= ebbot(i)*ebtop(i)
      ENDDO

*----------------------------------------------------------------------
      RETURN
      END
