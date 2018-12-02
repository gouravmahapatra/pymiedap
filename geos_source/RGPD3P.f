* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.
      SUBROUTINE RGPD3P(XD,YD,ZD,PDD,NXD,NYD)

Cf2py intent(in) XD,YD,ZD,NXD,NYD
Cf2py intent(in,out) PDD

*----------------------------------------------------------------------------
*     Implementation of RGPD3P in pymiedap.
*
*     Implemented by: Ashwyn Groot
*     Date: November 2018
*----------------------------------------------------------------------------
*
* Partial derivatives of a bivariate function on a rectangular grid
* (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08
*
* This subroutine estimates three partial derivatives, zx, zy, and
* zxy, of a bivariate function, z(x,y), on a rectangular grid in
* the x-y plane.  It is based on the revised Akima method that has
* the accuracy of a bicubic polynomial.
*
* The input arguments are
*   NXD = number of the input-grid data points in the x
*         coordinate (must be 2 or greater),
*   NYD = number of the input-grid data points in the y
*         coordinate (must be 2 or greater),
*   XD  = array of dimension NXD containing the x coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   YD  = array of dimension NYD containing the y coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   ZD  = two-dimensional array of dimension NXD*NYD
*         containing the z(x,y) values at the input-grid data
*         points.
*
* The output argument is
*   PDD = three-dimensional array of dimension 3*NXD*NYD,
*         where the estimated zx, zy, and zxy values at the
*         input-grid data points are to be stored.
*
*
* Specification statements
*     .. Scalar Arguments ..
      INTEGER          NXD,NYD
*     ..
*     .. Array Arguments ..
      REAL*8           PDD(3,NXD,NYD),XD(NXD),YD(NYD),ZD(NXD,NYD)
*     ..
*     .. Local Scalars ..
      REAL*8           B00,B00X,B00Y,B01,B10,B11,CX1,CX2,CX3,CY1,CY2,
     +                 CY3,DISF,DNM,DZ00,DZ01,DZ02,DZ03,DZ10,DZ11,DZ12,
     +                 DZ13,DZ20,DZ21,DZ22,DZ23,DZ30,DZ31,DZ32,DZ33,
     +                 DZX10,DZX20,DZX30,DZXY11,DZXY12,DZXY13,DZXY21,
     +                 DZXY22,DZXY23,DZXY31,DZXY32,DZXY33,DZY01,DZY02,
     +                 DZY03,EPSLN,PEZX,PEZXY,PEZY,SMPEF,SMPEI,SMWTF,
     +                 SMWTI,SX,SXX,SXXY,SXXYY,SXY,SXYY,SXYZ,SXZ,SY,SYY,
     +                 SYZ,SZ,VOLF,WT,X0,X1,X2,X3,XX1,XX2,XX3,Y0,Y1,Y2,
     +                 Y3,Z00,Z01,Z02,Z03,Z10,Z11,Z12,Z13,Z20,Z21,Z22,
     +                 Z23,Z30,Z31,Z32,Z33,ZXDI,ZXYDI,ZYDI,ZZ0,ZZ1,ZZ2
      INTEGER          IPEX,IPEY,IX0,IX1,IX2,IX3,IY0,IY1,IY2,IY3,JPEXY,
     +                 JXY,NX0,NY0
*     ..
*     .. Local Arrays ..
      REAL*8           B00XA(4),B00YA(4),B01A(4),B10A(4),CXA(3,4),
     +                 CYA(3,4),SXA(4),SXXA(4),SYA(4),SYYA(4),XA(3,4),
     +                 YA(3,4),Z0IA(3,4),ZI0A(3,4)
      INTEGER          IDLT(3,4)
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MAX
*     ..
*     .. Statement Functions ..
      REAL*8           Z2F,Z3F
*     ..
* Data statements 
      DATA             ((IDLT(JXY,JPEXY),JPEXY=1,4),JXY=1,3)/-3,-2,-1,1,
     +                 -2,-1,1,2,-1,1,2,3/
*     ..
* Statement Function definitions 
      Z2F(XX1,XX2,ZZ0,ZZ1) = (ZZ1-ZZ0)*XX2/XX1 + ZZ0
      Z3F(XX1,XX2,XX3,ZZ0,ZZ1,ZZ2) = ((ZZ2-ZZ0)* (XX3-XX1)/XX2-
     +                               (ZZ1-ZZ0)* (XX3-XX2)/XX1)*
     +                               (XX3/ (XX2-XX1)) + ZZ0
*     ..
* Calculation
* Initial setting of some local variables
      NX0 = MAX(4,NXD)
      NY0 = MAX(4,NYD)
* Double DO-loop with respect to the input grid points
      DO 60 IY0 = 1,NYD
          DO 50 IX0 = 1,NXD
              X0 = XD(IX0)
              Y0 = YD(IY0)
              Z00 = ZD(IX0,IY0)
* Part 1.  Estimation of ZXDI
* Initial setting
              SMPEF = 0.0
              SMWTF = 0.0
              SMPEI = 0.0
              SMWTI = 0.0
* DO-loop with respect to the primary estimate
              DO 10 IPEX = 1,4
* Selects necessary grid points in the x direction.
                  IX1 = IX0 + IDLT(1,IPEX)
                  IX2 = IX0 + IDLT(2,IPEX)
                  IX3 = IX0 + IDLT(3,IPEX)
                  IF ((IX1.LT.1) .OR. (IX2.LT.1) .OR. (IX3.LT.1) .OR.
     +                (IX1.GT.NX0) .OR. (IX2.GT.NX0) .OR.
     +                (IX3.GT.NX0)) GO TO 10
* Selects and/or supplements the x and z values.
                  X1 = XD(IX1) - X0
                  Z10 = ZD(IX1,IY0)
                  IF (NXD.GE.4) THEN
                      X2 = XD(IX2) - X0
                      X3 = XD(IX3) - X0
                      Z20 = ZD(IX2,IY0)
                      Z30 = ZD(IX3,IY0)
                  ELSE IF (NXD.EQ.3) THEN
                      X2 = XD(IX2) - X0
                      Z20 = ZD(IX2,IY0)
                      X3 = 2*XD(3) - XD(2) - X0
                      Z30 = Z3F(X1,X2,X3,Z00,Z10,Z20)
                  ELSE IF (NXD.EQ.2) THEN
                      X2 = 2*XD(2) - XD(1) - X0
                      Z20 = Z2F(X1,X2,Z00,Z10)
                      X3 = 2*XD(1) - XD(2) - X0
                      Z30 = Z2F(X1,X3,Z00,Z10)
                  END IF
                  DZX10 = (Z10-Z00)/X1
                  DZX20 = (Z20-Z00)/X2
                  DZX30 = (Z30-Z00)/X3
* Calculates the primary estimate of partial derivative zx as
* the coefficient of the bicubic polynomial.
                  CX1 = X2*X3/ ((X1-X2)* (X1-X3))
                  CX2 = X3*X1/ ((X2-X3)* (X2-X1))
                  CX3 = X1*X2/ ((X3-X1)* (X3-X2))
                  PEZX = CX1*DZX10 + CX2*DZX20 + CX3*DZX30
* Calculates the volatility factor and distance factor in the x
* direction for the primary estimate of zx.
                  SX = X1 + X2 + X3
                  SZ = Z00 + Z10 + Z20 + Z30
                  SXX = X1*X1 + X2*X2 + X3*X3
                  SXZ = X1*Z10 + X2*Z20 + X3*Z30
                  DNM = 4.0*SXX - SX*SX
                  B00 = (SXX*SZ-SX*SXZ)/DNM
                  B10 = (4.0*SXZ-SX*SZ)/DNM
                  DZ00 = Z00 - B00
                  DZ10 = Z10 - (B00+B10*X1)
                  DZ20 = Z20 - (B00+B10*X2)
                  DZ30 = Z30 - (B00+B10*X3)
                  VOLF = DZ00**2 + DZ10**2 + DZ20**2 + DZ30**2
                  DISF = SXX
* Calculates the EPSLN value, which is used to decide whether or
* not the volatility factor is essentially zero.
                  EPSLN = (Z00**2+Z10**2+Z20**2+Z30**2)*1.0E-12
* Accumulates the weighted primary estimates of zx and their
* weights.
                  IF (VOLF.GT.EPSLN) THEN
* - For a finite weight.
                      WT = 1.0/ (VOLF*DISF)
                      SMPEF = SMPEF + WT*PEZX
                      SMWTF = SMWTF + WT
                  ELSE
* - For an infinite weight.
                      SMPEI = SMPEI + PEZX
                      SMWTI = SMWTI + 1.0
                  END IF
* Saves the necessary values for estimating zxy
                  XA(1,IPEX) = X1
                  XA(2,IPEX) = X2
                  XA(3,IPEX) = X3
                  ZI0A(1,IPEX) = Z10
                  ZI0A(2,IPEX) = Z20
                  ZI0A(3,IPEX) = Z30
                  CXA(1,IPEX) = CX1
                  CXA(2,IPEX) = CX2
                  CXA(3,IPEX) = CX3
                  SXA(IPEX) = SX
                  SXXA(IPEX) = SXX
                  B00XA(IPEX) = B00
                  B10A(IPEX) = B10
   10         CONTINUE
* Calculates the final estimate of zx.
              IF (SMWTI.LT.0.5) THEN
* - When no infinite weights exist.
                  ZXDI = SMPEF/SMWTF
              ELSE
* - When infinite weights exist.
                  ZXDI = SMPEI/SMWTI
              END IF
* End of Part 1.
* Part 2.  Estimation of ZYDI
* Initial setting
              SMPEF = 0.0
              SMWTF = 0.0
              SMPEI = 0.0
              SMWTI = 0.0
* DO-loop with respect to the primary estimate
              DO 20 IPEY = 1,4
* Selects necessary grid points in the y direction.
                  IY1 = IY0 + IDLT(1,IPEY)
                  IY2 = IY0 + IDLT(2,IPEY)
                  IY3 = IY0 + IDLT(3,IPEY)
                  IF ((IY1.LT.1) .OR. (IY2.LT.1) .OR. (IY3.LT.1) .OR.
     +                (IY1.GT.NY0) .OR. (IY2.GT.NY0) .OR.
     +                (IY3.GT.NY0)) GO TO 20
* Selects and/or supplements the y and z values.
                  Y1 = YD(IY1) - Y0
                  Z01 = ZD(IX0,IY1)
                  IF (NYD.GE.4) THEN
                      Y2 = YD(IY2) - Y0
                      Y3 = YD(IY3) - Y0
                      Z02 = ZD(IX0,IY2)
                      Z03 = ZD(IX0,IY3)
                  ELSE IF (NYD.EQ.3) THEN
                      Y2 = YD(IY2) - Y0
                      Z02 = ZD(IX0,IY2)
                      Y3 = 2*YD(3) - YD(2) - Y0
                      Z03 = Z3F(Y1,Y2,Y3,Z00,Z01,Z02)
                  ELSE IF (NYD.EQ.2) THEN
                      Y2 = 2*YD(2) - YD(1) - Y0
                      Z02 = Z2F(Y1,Y2,Z00,Z01)
                      Y3 = 2*YD(1) - YD(2) - Y0
                      Z03 = Z2F(Y1,Y3,Z00,Z01)
                  END IF
                  DZY01 = (Z01-Z00)/Y1
                  DZY02 = (Z02-Z00)/Y2
                  DZY03 = (Z03-Z00)/Y3
* Calculates the primary estimate of partial derivative zy as
* the coefficient of the bicubic polynomial.
                  CY1 = Y2*Y3/ ((Y1-Y2)* (Y1-Y3))
                  CY2 = Y3*Y1/ ((Y2-Y3)* (Y2-Y1))
                  CY3 = Y1*Y2/ ((Y3-Y1)* (Y3-Y2))
                  PEZY = CY1*DZY01 + CY2*DZY02 + CY3*DZY03
* Calculates the volatility factor and distance factor in the y
* direction for the primary estimate of zy.
                  SY = Y1 + Y2 + Y3
                  SZ = Z00 + Z01 + Z02 + Z03
                  SYY = Y1*Y1 + Y2*Y2 + Y3*Y3
                  SYZ = Y1*Z01 + Y2*Z02 + Y3*Z03
                  DNM = 4.0*SYY - SY*SY
                  B00 = (SYY*SZ-SY*SYZ)/DNM
                  B01 = (4.0*SYZ-SY*SZ)/DNM
                  DZ00 = Z00 - B00
                  DZ01 = Z01 - (B00+B01*Y1)
                  DZ02 = Z02 - (B00+B01*Y2)
                  DZ03 = Z03 - (B00+B01*Y3)
                  VOLF = DZ00**2 + DZ01**2 + DZ02**2 + DZ03**2
                  DISF = SYY
* Calculates the EPSLN value, which is used to decide whether or
* not the volatility factor is essentially zero.
                  EPSLN = (Z00**2+Z01**2+Z02**2+Z03**2)*1.0E-12
* Accumulates the weighted primary estimates of zy and their
* weights.
                  IF (VOLF.GT.EPSLN) THEN
* - For a finite weight.
                      WT = 1.0/ (VOLF*DISF)
                      SMPEF = SMPEF + WT*PEZY
                      SMWTF = SMWTF + WT
                  ELSE
* - For an infinite weight.
                      SMPEI = SMPEI + PEZY
                      SMWTI = SMWTI + 1.0
                  END IF
* Saves the necessary values for estimating zxy
                  YA(1,IPEY) = Y1
                  YA(2,IPEY) = Y2
                  YA(3,IPEY) = Y3
                  Z0IA(1,IPEY) = Z01
                  Z0IA(2,IPEY) = Z02
                  Z0IA(3,IPEY) = Z03
                  CYA(1,IPEY) = CY1
                  CYA(2,IPEY) = CY2
                  CYA(3,IPEY) = CY3
                  SYA(IPEY) = SY
                  SYYA(IPEY) = SYY
                  B00YA(IPEY) = B00
                  B01A(IPEY) = B01
   20         CONTINUE
* Calculates the final estimate of zy.
              IF (SMWTI.LT.0.5) THEN
* - When no infinite weights exist.
                  ZYDI = SMPEF/SMWTF
              ELSE
* - When infinite weights exist.
                  ZYDI = SMPEI/SMWTI
              END IF
* End of Part 2.
* Part 3.  Estimation of ZXYDI
* Initial setting
              SMPEF = 0.0
              SMWTF = 0.0
              SMPEI = 0.0
              SMWTI = 0.0
* Outer DO-loops with respect to the primary estimates in the x
* direction
              DO 40 IPEX = 1,4
                  IX1 = IX0 + IDLT(1,IPEX)
                  IX2 = IX0 + IDLT(2,IPEX)
                  IX3 = IX0 + IDLT(3,IPEX)
                  IF ((IX1.LT.1) .OR. (IX2.LT.1) .OR. (IX3.LT.1) .OR.
     +                (IX1.GT.NX0) .OR. (IX2.GT.NX0) .OR.
     +                (IX3.GT.NX0)) GO TO 40
* Retrieves the necessary values for estimating zxy in the x
* direction.
                  X1 = XA(1,IPEX)
                  X2 = XA(2,IPEX)
                  X3 = XA(3,IPEX)
                  Z10 = ZI0A(1,IPEX)
                  Z20 = ZI0A(2,IPEX)
                  Z30 = ZI0A(3,IPEX)
                  CX1 = CXA(1,IPEX)
                  CX2 = CXA(2,IPEX)
                  CX3 = CXA(3,IPEX)
                  SX = SXA(IPEX)
                  SXX = SXXA(IPEX)
                  B00X = B00XA(IPEX)
                  B10 = B10A(IPEX)
* Inner DO-loops with respect to the primary estimates in the y
* direction
                  DO 30 IPEY = 1,4
                      IY1 = IY0 + IDLT(1,IPEY)
                      IY2 = IY0 + IDLT(2,IPEY)
                      IY3 = IY0 + IDLT(3,IPEY)
                      IF ((IY1.LT.1) .OR. (IY2.LT.1) .OR.
     +                    (IY3.LT.1) .OR. (IY1.GT.NY0) .OR.
     +                    (IY2.GT.NY0) .OR. (IY3.GT.NY0)) GO TO 30
* Retrieves the necessary values for estimating zxy in the y
* direction.
                      Y1 = YA(1,IPEY)
                      Y2 = YA(2,IPEY)
                      Y3 = YA(3,IPEY)
                      Z01 = Z0IA(1,IPEY)
                      Z02 = Z0IA(2,IPEY)
                      Z03 = Z0IA(3,IPEY)
                      CY1 = CYA(1,IPEY)
                      CY2 = CYA(2,IPEY)
                      CY3 = CYA(3,IPEY)
                      SY = SYA(IPEY)
                      SYY = SYYA(IPEY)
                      B00Y = B00YA(IPEY)
                      B01 = B01A(IPEY)
* Selects and/or supplements the z values.
                      IF (NYD.GE.4) THEN
                          Z11 = ZD(IX1,IY1)
                          Z12 = ZD(IX1,IY2)
                          Z13 = ZD(IX1,IY3)
                          IF (NXD.GE.4) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z23 = ZD(IX2,IY3)
                              Z31 = ZD(IX3,IY1)
                              Z32 = ZD(IX3,IY2)
                              Z33 = ZD(IX3,IY3)
                          ELSE IF (NXD.EQ.3) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z23 = ZD(IX2,IY3)
                              Z31 = Z3F(X1,X2,X3,Z01,Z11,Z21)
                              Z32 = Z3F(X1,X2,X3,Z02,Z12,Z22)
                              Z33 = Z3F(X1,X2,X3,Z03,Z13,Z23)
                          ELSE IF (NXD.EQ.2) THEN
                              Z21 = Z2F(X1,X2,Z01,Z11)
                              Z22 = Z2F(X1,X2,Z02,Z12)
                              Z23 = Z2F(X1,X2,Z03,Z13)
                              Z31 = Z2F(X1,X3,Z01,Z11)
                              Z32 = Z2F(X1,X3,Z02,Z12)
                              Z33 = Z2F(X1,X3,Z03,Z13)
                          END IF
                      ELSE IF (NYD.EQ.3) THEN
                          Z11 = ZD(IX1,IY1)
                          Z12 = ZD(IX1,IY2)
                          Z13 = Z3F(Y1,Y2,Y3,Z10,Z11,Z12)
                          IF (NXD.GE.4) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z31 = ZD(IX3,IY1)
                              Z32 = ZD(IX3,IY2)
                          ELSE IF (NXD.EQ.3) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z31 = Z3F(X1,X2,X3,Z01,Z11,Z21)
                              Z32 = Z3F(X1,X2,X3,Z02,Z12,Z22)
                          ELSE IF (NXD.EQ.2) THEN
                              Z21 = Z2F(X1,X2,Z01,Z11)
                              Z22 = Z2F(X1,X2,Z02,Z12)
                              Z31 = Z2F(X1,X3,Z01,Z11)
                              Z32 = Z2F(X1,X3,Z02,Z12)
                          END IF
                          Z23 = Z3F(Y1,Y2,Y3,Z20,Z21,Z22)
                          Z33 = Z3F(Y1,Y2,Y3,Z30,Z31,Z32)
                      ELSE IF (NYD.EQ.2) THEN
                          Z11 = ZD(IX1,IY1)
                          Z12 = Z2F(Y1,Y2,Z10,Z11)
                          Z13 = Z2F(Y1,Y3,Z10,Z11)
                          IF (NXD.GE.4) THEN
                              Z21 = ZD(IX2,IY1)
                              Z31 = ZD(IX3,IY1)
                          ELSE IF (NXD.EQ.3) THEN
                              Z21 = ZD(IX2,IY1)
                              Z31 = Z3F(X1,X2,X3,Z01,Z11,Z21)
                          ELSE IF (NXD.EQ.2) THEN
                              Z21 = Z2F(X1,X2,Z01,Z11)
                              Z31 = Z2F(X1,X3,Z01,Z11)
                          END IF
                          Z22 = Z2F(Y1,Y2,Z20,Z21)
                          Z23 = Z2F(Y1,Y3,Z20,Z21)
                          Z32 = Z2F(Y1,Y2,Z30,Z31)
                          Z33 = Z2F(Y1,Y3,Z30,Z31)
                      END IF
* Calculates the primary estimate of partial derivative zxy as
* the coefficient of the bicubic polynomial.
                      DZXY11 = (Z11-Z10-Z01+Z00)/ (X1*Y1)
                      DZXY12 = (Z12-Z10-Z02+Z00)/ (X1*Y2)
                      DZXY13 = (Z13-Z10-Z03+Z00)/ (X1*Y3)
                      DZXY21 = (Z21-Z20-Z01+Z00)/ (X2*Y1)
                      DZXY22 = (Z22-Z20-Z02+Z00)/ (X2*Y2)
                      DZXY23 = (Z23-Z20-Z03+Z00)/ (X2*Y3)
                      DZXY31 = (Z31-Z30-Z01+Z00)/ (X3*Y1)
                      DZXY32 = (Z32-Z30-Z02+Z00)/ (X3*Y2)
                      DZXY33 = (Z33-Z30-Z03+Z00)/ (X3*Y3)
                      PEZXY = CX1* (CY1*DZXY11+CY2*DZXY12+CY3*DZXY13) +
     +                        CX2* (CY1*DZXY21+CY2*DZXY22+CY3*DZXY23) +
     +                        CX3* (CY1*DZXY31+CY2*DZXY32+CY3*DZXY33)
* Calculates the volatility factor and distance factor in the x
* and y directions for the primary estimate of zxy.
                      B00 = (B00X+B00Y)/2.0
                      SXY = SX*SY
                      SXXY = SXX*SY
                      SXYY = SX*SYY
                      SXXYY = SXX*SYY
                      SXYZ = X1* (Y1*Z11+Y2*Z12+Y3*Z13) +
     +                       X2* (Y1*Z21+Y2*Z22+Y3*Z23) +
     +                       X3* (Y1*Z31+Y2*Z32+Y3*Z33)
                      B11 = (SXYZ-B00*SXY-B10*SXXY-B01*SXYY)/SXXYY
                      DZ00 = Z00 - B00
                      DZ01 = Z01 - (B00+B01*Y1)
                      DZ02 = Z02 - (B00+B01*Y2)
                      DZ03 = Z03 - (B00+B01*Y3)
                      DZ10 = Z10 - (B00+B10*X1)
                      DZ11 = Z11 - (B00+B01*Y1+X1* (B10+B11*Y1))
                      DZ12 = Z12 - (B00+B01*Y2+X1* (B10+B11*Y2))
                      DZ13 = Z13 - (B00+B01*Y3+X1* (B10+B11*Y3))
                      DZ20 = Z20 - (B00+B10*X2)
                      DZ21 = Z21 - (B00+B01*Y1+X2* (B10+B11*Y1))
                      DZ22 = Z22 - (B00+B01*Y2+X2* (B10+B11*Y2))
                      DZ23 = Z23 - (B00+B01*Y3+X2* (B10+B11*Y3))
                      DZ30 = Z30 - (B00+B10*X3)
                      DZ31 = Z31 - (B00+B01*Y1+X3* (B10+B11*Y1))
                      DZ32 = Z32 - (B00+B01*Y2+X3* (B10+B11*Y2))
                      DZ33 = Z33 - (B00+B01*Y3+X3* (B10+B11*Y3))
                      VOLF = DZ00**2 + DZ01**2 + DZ02**2 + DZ03**2 +
     +                       DZ10**2 + DZ11**2 + DZ12**2 + DZ13**2 +
     +                       DZ20**2 + DZ21**2 + DZ22**2 + DZ23**2 +
     +                       DZ30**2 + DZ31**2 + DZ32**2 + DZ33**2
                      DISF = SXX*SYY
* Calculates EPSLN.
                      EPSLN = (Z00**2+Z01**2+Z02**2+Z03**2+Z10**2+
     +                        Z11**2+Z12**2+Z13**2+Z20**2+Z21**2+Z22**2+
     +                        Z23**2+Z30**2+Z31**2+Z32**2+Z33**2)*
     +                        1.0E-12
* Accumulates the weighted primary estimates of zxy and their
* weights.
                      IF (VOLF.GT.EPSLN) THEN
* - For a finite weight.
                          WT = 1.0/ (VOLF*DISF)
                          SMPEF = SMPEF + WT*PEZXY
                          SMWTF = SMWTF + WT
                      ELSE
* - For an infinite weight.
                          SMPEI = SMPEI + PEZXY
                          SMWTI = SMWTI + 1.0
                      END IF
   30             CONTINUE
   40         CONTINUE
* Calculates the final estimate of zxy.
              IF (SMWTI.LT.0.5) THEN
* - When no infinite weights exist.
                  ZXYDI = SMPEF/SMWTF
              ELSE
* - When infinite weights exist.
                  ZXYDI = SMPEI/SMWTI
              END IF
* End of Part 3
              PDD(1,IX0,IY0) = ZXDI
              PDD(2,IX0,IY0) = ZYDI
              PDD(3,IX0,IY0) = ZXYDI
   50     CONTINUE
   60 CONTINUE
      RETURN
      END
