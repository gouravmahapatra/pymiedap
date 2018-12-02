* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE ord2m(xmu,smf,nmat,
     .                 Zmplus,Zmmin,a,b,ebmu,nmu,nsup,Rm,Tm)

*----------------------------------------------------------------------*
*  Calculate second order scattering contribution to the m-th Fourier  *
*  component of reflection and transmission of a homogeneous layer.    *
*  Formulae can be found in                                            *
*      Hovenier, (1971) Astron. Astrophys. 13, p. 7-29                 *
*      Eqs. (A-26)-(A-38)                                              *
*  The second order scattering is added to whatever is in Rm and Tm.   *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER nmu,nmum,nsup

      REAL*8, DIMENSION(nmu) :: xmu, smf, ebmu !rank 1

      REAL*8, DIMENSION(nsup,nsup) :: Zmmin, Zmplus, Rm, Tm !rank 2

      DIMENSION yik(4,4)

Cf2py intent(out) Rm, Tm

************************************************************************
      nmum= nmu-1

      aa16=a*a/16.D0
      DO kl=1,nmat
         DO il=1,nmat
            yik(il,kl) = 1.D0
            IF ((kl.GE.3) .AND. (il.LE.2)) yik(il,kl) = -1.D0
            IF ((il.GE.3) .AND. (kl.LE.2)) yik(il,kl) = -1.D0
         ENDDO
      ENDDO

      DO 1100 i=1,nmu
         xi = xmu(i)
         ei = ebmu(i)
         asmfi = aa16*smf(i)
         im = (i-1)*nmat
         iz = 0
         IF (xi.GT.1.D-10) iz = 1
         bpi = 0.D0
         IF (iz.EQ.1D0) bpi = b/xi
         beigxi = bpi*ei
         DO 1000 j = 1,i
            xj = xmu(j)
            ej = ebmu(j)
            asmfij = asmfi*smf(j)
            jm = (j-1)*nmat
            jz = 0
            IF (xj.GT.1.D-10) jz = 1

*----------------------------------------------------------------------*
*  When both mu(i) and mu(j) are almost zero : set reflection          *
*  and transmission to zero and go to next mu(j) if there is one.      *
*----------------------------------------------------------------------*
            IF ((jz.EQ.0).AND.(iz.EQ.0)) THEN
               DO k=1,nmat
                  ik = im+k
                  DO l=1,nmat
                     jl = jm+l
                     Rm(ik,jl) = 0.D0
                     Tm(ik,jl) = 0.D0
                  ENDDO
               ENDDO 
            ENDIF

            bpj = 0.D0
            IF (jz.EQ.1) bpj = b/xj
            bij = bpi+bpj
            xipxj = xi+xj
            ximxj = xi-xj
            xjximxj = 0.D0
            IF (i.NE.j) xjximxj = xj/ximxj
            eimej = ei-ej
            eiej = 1.D0-ei*ej
*----------------------------------------------------------------------*
*  Taylor series to avoid loss of accuracy when b << 1.                *
*----------------------------------------------------------------------*
            IF (eiej.LT.1.D-3) THEN
               eiej = bij*(1.D0-bij/2.D0*(1.D0-bij/3.D0))
               eimej = bpj*(1.D0-bpj/2.D0*(1.D0-bpj/3.D0))-
     .                 bpi*(1.D0-bpi/2.D0*(1.D0-bpi/3.D0))
            ENDIF
            e1 = xjximxj*eimej
            g1 = xj/xipxj*eiej
*----------------------------------------------------------------------*
*  Start integration over mu'                                          *
*----------------------------------------------------------------------*
            DO k=1,nmum
               xk = xmu(k)
               ek = ebmu(k)
               km = (k-1)*nmat
               bpk = b/xk
               bik = bpi+bpk
               bjk = bpj+bpk
*----------------------------------------------------------------------*
*  Calculate the functions e, f, g and h in Hovenier (1971)            *
*----------------------------------------------------------------------*
               IF ((bpi.LT.1.D-3) .AND. (bpj.LT.1.D-3) 
     .                            .AND. (bpk.LT.1.D-3)) THEN
*----------------------------------------------------------------------*
*         Use Taylor series to avoid loss of accuracy when b << 1.     *
*----------------------------------------------------------------------*
                  IF (iz.EQ.0) THEN
                     z = b/xk/xj
                     e = 0.D0
                     g = z*(1.D0-bjk/2.D0*(1.D0-bjk/3.D0))
                     f = g-z*z*b/6.D0
                     h = 0.D0
                  ELSEIF (jz.EQ.0) THEN
                     z = b/xk/xi
                     e = 0.D0
                     h = z*(1.D0-bik/2.D0*(1.D0-bik/3.D0))
                     f = h-z*z*b/6.D0
                     g = 0.D0
                  ELSE
                     z = bpi*bpj*bpk/b/2.D0
                     e = z*(1.D0-(bpk+2.D0*(bpi+bpj))/3.D0)
                     f = z*(1.D0-(bpk+bpi+bpj)/3.D0)
                     g = z*(1.D0-(bpk+bpi+2.D0*bpj)/3.D0)
                     h = z*(1.D0-(bpk+bpj+2.D0*bpi)/3.D0)
                  ENDIF
               ELSE
*----------------------------------------------------------------------*
*         No Taylor series is needed.                                  *
*----------------------------------------------------------------------*
                  xipxk = xi+xk
                  xjpxk = xj+xk
                  ximxk = xi-xk
                  xjmxk = xj-xk
                  eimek = ei-ek
                  ejmek = ej-ek
                  eiek = 1.D0-ei*ek
                  ejek = 1.D0-ej*ek
*----------------------------------------------------------------------*
*         Use Taylor series to avoid loss of accuracy when b << 1.     *
*----------------------------------------------------------------------*
                  IF (eiek.LT.1.D-3) THEN
                     eiek = bik*(1.D0-bik/2.D0*(1.D0-bik/3.D0))
                     eimek = bpk*(1.D0-bpk/2.D0*(1.D0-bpk/3.D0))-
     .                       bpi*(1.D0-bpi/2.D0*(1.D0-bpi/3.D0))
                  ENDIF
                  IF (ejek.LT.1.D-3) THEN
                     ejek = bjk*(1.D0-bjk/2.D0*(1.D0-bjk/3.D0))
                     ejmek = bpk*(1.D0-bpk/2.D0*(1.D0-bpk/3.D0))-
     .                       bpj*(1.D0-bpj/2.D0*(1.D0-bpj/3.D0))
                  ENDIF
                  IF (i.EQ.j) THEN
                     IF (i.NE.k) THEN
                        e = (b/xj*ej-xk/xipxk*ej*ejek)/xjpxk
                        f = (b/xj*ej-xk/xjmxk*ejmek)/xjmxk
                        g = (g1-xk/ximxk*ej*eimek)/xjpxk
                        h = (g1-xk/xipxk*eiek)/xjmxk
                     ELSE
                        e = (b/xj*ej-xk/xipxk*ej*ejek)/xjpxk
                        f = b*b/2.D0/xj/xj/xj*ej
                        g = (g1-b/xi*ei*ej)/xipxj
                        h = g
                     ENDIF
                  ELSEIF (i.EQ.k) THEN
                     e = (e1-xk/xipxk*ej*eiek)/xjpxk
                     f = (beigxi-e1)/xj*xjximxj
                     g = (g1-beigxi*ej)/xipxj
                     h = (g1-xk/xipxk*eiek)/xjmxk
                  ELSEIF (j.EQ.k) THEN
                     e = (e1-xk/xipxk*ej*eiek)/xjpxk
                     f = (xi/xj*e1-b/xj*ej)/xj*xjximxj
                     g = (g1-xk/ximxk*ej*eimek)/xjpxk
                     h = (xi/xipxj*eiej-b/xj*ei*ej)/xipxj
                  ELSE
                     e = (e1-xk/xipxk*ej*eiek)/xjpxk
                     f = (e1-xk/ximxk*eimek)/xjmxk
                     g = (g1-xk/ximxk*ej*eimek)/xjpxk
                     h = (g1-xk/xipxk*eiek)/xjmxk
                  ENDIF
               ENDIF
*----------------------------------------------------------------------*
*  End of calculation of functions e, f, g and h                       *
*----------------------------------------------------------------------*
               y = smf(k)*smf(k)*asmfij/xk
               DO il=1,nmat
                  iil= im+il
                  DO jl=1,nmat
                     jjl= jm+jl
                     DO kl=1,nmat
                        kkl = km+kl
                        Zpnik = Zmplus(iil,kkl)
                        Zpnkj = Zmplus(kkl,jjl)
                        Zmnik = Zmmin(iil,kkl)
                        Zmnkj = Zmmin(kkl,jjl)
                        Zmsik = yik(il,kl)*Zmnik
                        Zpsik = yik(il,kl)*Zpnik
                        Rm(iil,jjl) = Rm(iil,jjl) +
     .                                y*(Zpsik*Zmnkj*g+Zmnik*Zpnkj*h)
                        Tm(iil,jjl) = Tm(iil,jjl) +
     .                                y*(Zmsik*Zmnkj*e+Zpnik*Zpnkj*f)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
*----------------------------------------------------------------------*
*  End of integration over mu'                                         *
*----------------------------------------------------------------------*
1000     CONTINUE

*----------------------------------------------------------------------*
*  End of loop over column index j                                     *
*----------------------------------------------------------------------*
1100  CONTINUE

*----------------------------------------------------------------------*
*  End of loop over row index i                                        *
*  Fill upper triangle of the supermatrices using symmetry relations.  *
*----------------------------------------------------------------------*
      CALL fillup(Rm,Tm,nmat,nmu,nsup)

*-----------------------------------------------------------------------
      RETURN
      END
