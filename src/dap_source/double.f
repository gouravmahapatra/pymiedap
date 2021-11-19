* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE double(Rm,Tm,ebmu,nmu,nmat)

*----------------------------------------------------------------------*
*  Calculate the m-th Fourier term of reflection ans transmission of   *
*  a homogeneous layer from the reflection and transmission of a layer *
*  with only half the optical thickness.                               *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER nmu,nmum,nmat

      PARAMETER(maxrep=15,trmin=1.D-8)

C      DOUBLE PRECISION X(nsupMAX,nsupMAX),Y(nsupMAX,nsupMAX),
C     .          Z(nsupMAX,nsupMAX),Rm(nsupMAX,nsupMAX),
C     .          Tm(nsupMAX,nsupMAX),ebmu(nmuMAX),
C     .          E(nsupMAX)


      DIMENSION X(nsupMAX,nsupMAX),Y(nsupMAX,nsupMAX),
     .          Z(nsupMAX,nsupMAX),Rm(nsupMAX,nsupMAX),
     .          Tm(nsupMAX,nsupMAX),ebmu(nmuMAX),
     .          E(nsupMAX)

      LOGICAL verbo
      verbo = .false.

Cf2py intent(in,out) Rm,Tm,ebmu

*----------------------------------------------------------------------*
      nmum= nmu-1

*----------------------------------------------------------------------*
*  Use the product method to calculate the repeated reflections        *
*  between the layers, se de Haan et al. (1987) Eqs. (111)-(115)       *
*----------------------------------------------------------------------*
      CALL star(X,Rm,nmat,nmu)
*                                                    X = R*
      CALL prod(Y,X,Rm,nmat,nmu,nmum)
*                                                    Y = R*R = C1
      CALL assign(Z,Y,nmat,nmu)
*                                                    Z = C1 = S1
      ir= 0
  100     ir= ir+1
          CALL prod(X,Y,Y,nmat,nmu,nmum)
*                                                    X = Cr Cr = Cr+1
          CALL assign(Y,X,nmat,nmu)
*                                                    Y = Cr+1
          CALL prod(X,Z,Y,nmat,nmu,nmum)
*                                                    X = Sr Cr+1
          CALL addsm(Z,Z,X,nmat,nmu)
*                                                    Z = Sr + Sr Cr+1
          CALL addsm(Z,Z,Y,nmat,nmu)
*                                         Z = Sr + Sr Cr+1 + Cr+1 = Sr+1
          CALL trace(Y,trC,nmat,nmum)
*                                                    trC = trace(Cr+1)
          IF (verbo) PRINT *,' double: r = ',ir,' trace = ',trC
          IF ((trC.GT.trmin) .AND. (ir.LE.maxrep)) GOTO 100

*-----------------------------------------------------------------------
      IF (ir.GT.maxrep) THEN
         PRINT *,' double: WARNING repeated reflections did not'
     .           ,' converge after',maxrep,' steps'
         PRINT *,'         proceed anyway !'
      ENDIF

*----------------------------------------------------------------------*
*  Now Z contains the matrix Q Eq. (115)                               *
*  Use the adding Eqs. (85)-(91) with identical layers                 *
*----------------------------------------------------------------------*
      DO i=1,nmu
         DO k=1,nmat
            E((i-1)*nmat+k)= ebmu(i)
         ENDDO
      ENDDO

      CALL prod(X,Z,Tm,nmat,nmu,nmum)
*                                              X = QT
      CALL rdiapr(Z,Z,E,nmat,nmu)
*                                              Z = QE
      CALL addsm(X,X,Z,nmat,nmu)
*                                              X = QT + QE
      CALL addsm(X,X,Tm,nmat,nmu)
*                                              X = T + QT +QE = D
*-----------------------------------------------------------------------
      CALL prod(Z,Rm,X,nmat,nmu,nmum)
*                                              Z = RD
      CALL rdiapr(Y,Rm,E,nmat,nmu)
*                                              Y = RE
      CALL addsm(Z,Z,Y,nmat,nmu)
*                                              Z = RE + RD = U
*-----------------------------------------------------------------------
      CALL prod(Y,Tm,X,nmat,nmu,nmum)
*                                              Y = TD
      CALL ldiapr(X,E,X,nmat,nmu)
*                                              X = ED
      CALL addsm(Y,Y,X,nmat,nmu)
*                                              Y = TD + ED
      CALL rdiapr(X,Tm,E,nmat,nmu)
*                                              X = TE
      CALL addsm(Y,Y,X,nmat,nmu)
*                                              Y = TD + ED + TE = Ttot
*-----------------------------------------------------------------------
      CALL star(Tm,Tm,nmat,nmu)
*                                              Tm = T*
      CALL prod(X,Tm,Z,nmat,nmu,nmum)
*                                              X = T*U
      CALL addsm(Rm,Rm,X,nmat,nmu)
*                                              Rm = R + T*U
      CALL ldiapr(Z,E,Z,nmat,nmu)
*                                              Z = EU
      CALL addsm(Rm,Rm,Z,nmat,nmu)
*                                              Rm = R + T*U + EU = Rtot

*-----------------------------------------------------------------------
*     Tm= Ttot :
*-----------------------------------------------------------------------
      CALL assign(Tm,Y,nmat,nmu)

*-----------------------------------------------------------------------
      RETURN
      END
