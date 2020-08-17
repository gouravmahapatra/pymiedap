* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE devel(ncoefs,nangle,u,wg,F,coefs)

************************************************************************
*  Calculate the expansion coefficients of the scattering matrix in    *
*  generalized spherical functions by numerical integration over the   *
*  scattering angle.                                                   *
************************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER ncoefs, nangle
      DOUBLE PRECISION u(nangMAX),wg(nangMAX),F(6,nangMAX),
     .                 coefs(nmatMAX,nmatMAX,0:ncoefsMAX),
     .                 P00(nangMAX,2),P02(nangMAX,2),
     .                 P22(nangMAX,2),P2m2(nangMAX,2)

Cf2py intent(in) ncoefs, nangle, u, wg, F
Cf2py intent(out) coefs

************************************************************************
*  Initialization                                                      *
************************************************************************
      qroot6 = -0.25D0*dsqrt(6.D0)

      DO j=0,ncoefsMAX
         DO i=1,nmatMAX
            DO ii=1,nmatMAX
               coefs(ii,i,j)=0.D0
            ENDDO
         ENDDO
      ENDDO

************************************************************************
*  Multiply the scattering matrix F with the weights w for all angles  *
*  We DO this here because otherwise it should be DOne for each l      *
************************************************************************
      DO k=1,6
         DO i=1,nangle
            F(k,i)= wg(i)*F(k,i)
         ENDDO
      ENDDO

************************************************************************
*  Start loop over the coefficient index l                             *
*  first update generalized spherical functions, then calculate coefs. *
*  lold and lnew are pointer-like indices used in recurrence           *
************************************************************************
      lnew= 1
      lold= 2

      DO 70 l=0,ncoefs
         IF (l.EQ.0) THEN

************************************************************************
*           Adding paper Eq. (77) with m=0                           *
************************************************************************
            DO i=1,nangle
               P00(i,lold) = 1.D0
               P00(i,lnew) = 0.D0
               P02(i,lold) = 0.D0
               P22(i,lold) = 0.D0
               P2m2(i,lold)= 0.D0
               P02(i,lnew) = 0.D0
               P22(i,lnew) = 0.D0
               P2m2(i,lnew)= 0.D0
            ENDDO
         ELSE
            fac1 = (2.D0*l-1.d0)/dble(l)
            fac2 = dble(l-1.d0)/dble(l)

************************************************************************
*           Adding paper Eq. (81) with m=0                           *
************************************************************************
            DO i=1,nangle
               P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
            ENDDO
         ENDIF
         IF (l.EQ.2) THEN
************************************************************************
*           Adding paper Eqs. (78) and (80)  with m=2                *
*           sql4 contains the factor dsqrt(l*l-4) needed in          *
*           the recurrence Eqs. (81) and (82)                        *
************************************************************************
            DO i=1,nangle
               P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
               P22(i,lold) = 0.25D0*(1.D0+u(i))*(1.D0+u(i))
               P2m2(i,lold)= 0.25D0*(1.D0-u(i))*(1.D0-u(i))
               P02(i,lnew) = 0.D0
               P22(i,lnew) = 0.D0
               P2m2(i,lnew)= 0.D0
            ENDDO
            sql41= 0.D0
         ELSE IF (l.GT.2) THEN
************************************************************************
*           Adding paper Eq. (82) with m=0 and m=2                   *
************************************************************************
            sql4 = sql41
            sql41= dsqrt(dble(l*l)-4.d0)
            twol1= 2.D0*dble(l)-1.d0
            tmp1 = twol1/sql41
            tmp2 = sql4/sql41
            denom= (dble(l)-1.d0)*(dble(l*l)-4.d0)
            fac1 = twol1*(dble(l)-1.d0)*dble(l)/denom
            fac2 = 4.D0*twol1/denom
            fac3 = dble(l)*((dble(l)-1.d0)*(dble(l)-1.d0)-4.d0)/denom
            DO i=1,nangle
               P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
               P22(i,lold) = (fac1*u(i)-fac2)*P22(i,lnew)
     +                                             - fac3*P22(i,lold)
               P2m2(i,lold)= (fac1*u(i)+fac2)*P2m2(i,lnew)
     +                                            - fac3*P2m2(i,lold)
            ENDDO
         ENDIF
************************************************************************
*        Switch indices so that lnew indicates the function with      *
*        the present index value l, this mechanism prevents swapping  *
*        of entire arrays.                                            *
************************************************************************
         itmp = lnew
         lnew = lold
         lold = itmp
************************************************************************
*        Now calculate the coefficients by integration over angle     *
*        See de Haan et al. (1987) Eqs. (68)-(73).                    *
*        Remember for Mie scattering : F11 = F22 and F33 = F44        *
************************************************************************
         alfap= 0.D0
         alfam= 0.D0
         DO i=1,nangle
            coefs(1,1,l) = coefs(1,1,l) + P00(i,lnew)*F(1,i)
            alfap = alfap + P22(i,lnew)*(F(2,i)+F(3,i))
            alfam = alfam + P2m2(i,lnew)*(F(2,i)-F(3,i))
            coefs(4,4,l) = coefs(4,4,l) + P00(i,lnew)*F(4,i)
            coefs(1,2,l) = coefs(1,2,l) + P02(i,lnew)*F(5,i)
            coefs(3,4,l) = coefs(3,4,l) + P02(i,lnew)*F(6,i)
         ENDDO
************************************************************************
*         Multiply with trivial factors like 0.5D0*(2*l+1)             *
************************************************************************
          fl = dble(l)+0.5D0
          coefs(1,1,l) =  fl*coefs(1,1,l)
          coefs(2,2,l) =  fl*0.5D0*(alfap+alfam)
          coefs(3,3,l) =  fl*0.5D0*(alfap-alfam)
          coefs(4,4,l) =  fl*coefs(4,4,l)
          coefs(1,2,l) =  fl*coefs(1,2,l)
          coefs(3,4,l) =  fl*coefs(3,4,l)
          coefs(2,1,l) =     coefs(1,2,l)
          coefs(4,3,l) =    -coefs(3,4,l)
   70 CONTINUE
************************************************************************
*     End of loop over index l                                         *
************************************************************************

************************************************************************
*     Remove the weight factor from the scattering matrix              *
************************************************************************
      DO k=1,6
         DO i=1,nangle
            F(k,i) = F(k,i)/wg(i)
         ENDDO
      ENDDO

*-----------------------------------------------------------------------
      RETURN
      END
