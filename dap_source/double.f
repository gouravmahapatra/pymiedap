* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE double(Rm,Tm,ebmu,nmu,nmat,nsup)

*----------------------------------------------------------------------*
*  Calculate the m-th Fourier term of reflection ans transmission of   *
*  a homogeneous layer from the reflection and transmission of a layer *
*  with only half the optical thickness.                               *
*  Edited by: Ashwyn Groot                                             *
*  Date: November 2018                                                 *
*  Introduced matrix operations with f95<                              *
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INCLUDE 'max_incl'

      INTEGER nmu,nmum,nmat,nsup

      PARAMETER(maxrep=15,trmin=1.D-8)

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: X, Y, Z, E !rank 2
      REAL*8, DIMENSION(nsup,nsup) :: Rm, Tm !rank 2
      REAL*8, DIMENSION(nmu) :: ebmu !rank 1

      LOGICAL verbo
      verbo = .false.
      ALLOCATE( X(nsup,nsup), Y(nsup,nsup), Z(nsup,nsup),
     .          E(nsup,nsup))

Cf2py intent(in,out) Rm,Tm,ebmu

*----------------------------------------------------------------------*
      nmum= nmu-1

*----------------------------------------------------------------------*
*  Use the product method to calculate the repeated reflections        *
*  between the layers, se de Haan et al. (1987) Eqs. (111)-(115)       *
*----------------------------------------------------------------------*
      CALL star(X,Rm,nmat,nmu,nsup)
*                                                    X = R*
      CALL prod(Y,X,Rm,nmat,nmu,nmum)
*                                                    Y = R*R = C1
      Z=Y
*                                                    Z = C1 = S1
      ir= 0
  100     ir= ir+1
          CALL prod(X,Y,Y,nmat,nmu,nmum)
*                                                    X = Cr Cr = Cr+1
          Y=X
*                                                    Y = Cr+1
          CALL prod(X,Z,Y,nmat,nmu,nmum)
*                                                    X = Sr Cr+1
          Z=Z+X
*                                                    Z = Sr + Sr Cr+1
          Z=Z+Y
*                                         Z = Sr + Sr Cr+1 + Cr+1 = Sr+1
          CALL trace(Y,nmat,nmum,nsup,trC)
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
      E=0.D0
      DO i=1,nmu
         DO k=1,nmat
            E((i-1)*nmat+k,(i-1)*nmat+k)= ebmu(i)
         ENDDO
      ENDDO

      CALL prod(X,Z,Tm,nmat,nmu,nmum)
*                                              X = QT
      Z=MATMUL(Z,E)
*                                              Z = QE
      X=X+Z
*                                              X = QT + QE
      X=X+Tm
*                                              X = T + QT +QE = D
*-----------------------------------------------------------------------
      CALL prod(Z,Rm,X,nmat,nmu,nmum)
*                                              Z = RD
      Y=MATMUL(Rm,E)
*                                              Y = RE
      Z=Z+Y
*                                              Z = RE + RD = U
*-----------------------------------------------------------------------
      CALL prod(Y,Tm,X,nmat,nmu,nmum)
*                                              Y = TD
      X=MATMUL(E,X)
*                                              X = ED
      Y=Y+X
*                                              Y = TD + ED
      X=MATMUL(Tm,E)
*                                              X = TE
      Y=Y+X
*                                              Y = TD + ED + TE = Ttot
*-----------------------------------------------------------------------
      CALL star(Tm,Tm,nmat,nmu,nsup)
*                                              Tm = T*
      CALL prod(X,Tm,Z,nmat,nmu,nmum)
*                                              X = T*U
      Rm=Rm+X
*                                              Rm = R + T*U
      Z=MATMUL(E,Z)
*                                              Z = EU
      Rm=Rm+Z
*                                              Rm = R + T*U + EU = Rtot

*-----------------------------------------------------------------------
*     Tm= Ttot :
*-----------------------------------------------------------------------
      Tm=Y

*-----------------------------------------------------------------------
      RETURN
      END
