      SUBROUTINE nobot(Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot,
     .                 ebtop,ebbot,nmu,nmat)

*----------------------------------------------------------------------
*     Use the adding equations when there is no scattering in the 
*     bottom layer:
*----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,k,nmat,nmu,nsup

      DOUBLE PRECISION Rmtop(nsupMAX,nsupMAX),Tmtop(nsupMAX,nsupMAX),
     .                 Rmbot(nsupMAX,nsupMAX),Tmbot(nsupMAX,nsupMAX),
     .                 Rmsbot(nsupMAX,nsupMAX),
     .                 ebtop(nmuMAX),ebbot(nmuMAX),
     .                 Etop(nsupMAX),Ebot(nsupMAX)

Cf2py intent(in,out) Rmtop,Tmtop,Rmbot,Tmbot,Rmsbot,ebtop,ebbot

*----------------------------------------------------------------------
      nsup= nmat*nmu

      DO i=1,nmu
         DO k=1,nmat
            Etop((i-1)*nmat+k) = ebtop(i)
            Ebot((i-1)*nmat+k) = ebbot(i)
         ENDDO
      ENDDO

*----------------------------------------------------------------------
*     Reflection  R= R' 
*----------------------------------------------------------------------
      CALL assign(Rmbot,Rmtop,nmat,nmu)

*----------------------------------------------------------------------
*     Transmission T= E"T'                                    
*----------------------------------------------------------------------
      CALL ldiapr(Tmbot,Ebot,Tmtop,nmat,nmu)

*----------------------------------------------------------------------
*     Reflection star R*= E"R'*E"                                      
*----------------------------------------------------------------------
      CALL star(Rmsbot,Rmtop,nmat,nmu)
      CALL rdiapr(Rmsbot,Rmsbot,Ebot,nmat,nmu)
      CALL ldiapr(Rmsbot,Ebot,Rmsbot,nmat,nmu)

*----------------------------------------------------------------------
*     Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu) 
*----------------------------------------------------------------------
      DO i=1,nmu
         ebbot(i)= ebbot(i)*ebtop(i)
      ENDDO

*----------------------------------------------------------------------
      RETURN
      END
