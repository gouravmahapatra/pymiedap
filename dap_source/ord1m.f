* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE ord1m(xmu,smf,nmu,nmat,Zmplus,Zmmin,a,b,ebmu,
     .                 Rm1,Tm1)

*----------------------------------------------------------------------
*     Calculate the first order scattering contribution to the m-th
*     Fourier component of reflection and transmission of a homogeneous 
*     layer. 
*
*     Formulae can be found in de Haan et al. (1987) Eqs. (A-21)-(A-23)                                                
*     The resulting supermatrices are returned through Rm1 and Tm1. 
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

C      INTEGER nmu,nmat

C      DOUBLE PRECISION a, b
C      DOUBLE PRECISION xmu(nmuMAX),smf(nmuMAX),ebmu(nmuMAX),
C     .          Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX),
C     .          Rm1(nsupMAX,nsupMAX),Tm1(nsupMAX,nsupMAX)

      DIMENSION xmu(nmuMAX),smf(nmuMAX),ebmu(nmuMAX),
     .          Zmmin(nsupMAX,nsupMAX),Zmplus(nsupMAX,nsupMAX),
     .          Rm1(nsupMAX,nsupMAX),Tm1(nsupMAX,nsupMAX)

Cf2py intent(out) Rm1,Tm1

*-----------------------------------------------------------------------
      quarta = 0.25D0*a

      DO 300 i=1,nmu
         xmui = xmu(i)
         ei   = ebmu(i)
         awi  = quarta*smf(i)
         im   = (i-1)*nmat
         DO j=1,nmu
            xmuj= xmu(j)
            ej  = ebmu(j)
            awij= awi*smf(j)
            jm  = (j-1)*nmat
            eiej= 1.D0-ei*ej
            IF (eiej.GT.1.D-3) THEN
               h = xmui+xmuj
               IF (h.GT.1.D-10) h = 1.D0/h
                  hR = awij*h*eiej
                  IF (dabs(xmui-xmuj).GT.1.D-10) THEN
                     hT = awij*(ei-ej)/(xmui-xmuj)
                  ELSE
                     h = 0.D0
                     IF (xmui.GT.1.D-10) h = b*ei/(xmui*xmui)
                     hT = awij*h
                  ENDIF
               ELSE

*-----------------------------------------------------------------------
*              Use Taylor series to avoid loss of accuracy when b<<1:
*-----------------------------------------------------------------------
                  bperi  = b/xmui
                  bperj  = b/xmuj
                  bperij = bperi+bperj
                  h      = 1.D0-0.5D0*bperij*(1.D0-bperij/3.D0)
                  y      = awij*bperi/xmuj
                  hR     = y*h
                  hT     = y*(h-bperi*bperj/6.D0)
              ENDIF
              DO k=1,nmat
                 ik= im+k
                 DO l=1,nmat
                    jl = jm+l
                    Rm1(ik,jl)= hR*Zmmin(ik,jl)
                    Tm1(ik,jl)= hT*Zmplus(ik,jl)
                 ENDDO
              ENDDO
           ENDDO
300   CONTINUE

*-----------------------------------------------------------------------
      RETURN
      END
