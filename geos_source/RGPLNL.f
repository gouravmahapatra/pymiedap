* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.
      SUBROUTINE RGPLNL(XD,YD,ZD,PDD,
     +                  XI,YI,INXI,INYI,ZI,NXD,NYD,NIP)

Cf2py intent(in) XD,YD,ZD,PDD,XI,YI,INXI,INYI,NXD,NYD,NIP
Cf2py intent(out) ZI

*----------------------------------------------------------------------------
*     Implementation of RGPLNL in pymiedap.
*
*     Implemented by: Ashwyn Groot
*     Date: November 2018
*----------------------------------------------------------------------------
*
* Polynomials for rectangular-grid bivariate interpolation and
* surface fitting
* (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08
*
* This subroutine determines a polynomial in x and y for a rectangle
* of the input grid in the x-y plane and calculates the z value for
* the desired points by evaluating the polynomial for rectangular-
* grid bivariate interpolation and surface fitting.
*
* The input arguments are
*   NXD  = number of the input-grid data points in the x
*          coordinate (must be 2 or greater),
*   NYD  = number of the input-grid data points in the y
*          coordinate (must be 2 or greater),
*   XD   = array of dimension NXD containing the x coordinates
*          of the input-grid data points (must be in a
*          monotonic increasing order),
*   YD   = array of dimension NYD containing the y coordinates
*          of the input-grid data points (must be in a
*          monotonic increasing order),
*   ZD   = two-dimensional array of dimension NXD*NYD
*          containing the z(x,y) values at the input-grid data
*          points,
*   PDD  = three-dimensional array of dimension 3*NXD*NYD
*          containing the estimated zx, zy, and zxy values
*          at the input-grid data points,
*   NIP  = number of the output points at which interpolation
*          is to be performed,
*   XI   = array of dimension NIP containing the x coordinates
*          of the output points,
*   YI   = array of dimension NIP containing the y coordinates
*          of the output points,
*   INXI = integer array of dimension NIP containing the
*          interval numbers of the input grid intervals in the
*          x direction where the x coordinates of the output
*          points lie,
*   INYI = integer array of dimension NIP containing the
*          interval numbers of the input grid intervals in the
*          y direction where the y coordinates of the output
*          points lie.
*
* The output argument is
*   ZI   = array of dimension NIP, where the interpolated z
*          values at the output points are to be stored.
*
*
* Specification statements
*     .. Scalar Arguments ..
      INTEGER          NIP,NXD,NYD
*     ..
*     .. Array Arguments ..
      REAL*8           PDD(3,NXD,NYD),XD(NXD),XI(NIP),YD(NYD),YI(NIP),
     +                 ZD(NXD,NYD),ZI(NIP)
      INTEGER          INXI(NIP),INYI(NIP)
*     ..
*     .. Local Scalars ..
      REAL*8           A,B,C,D,DX,DXSQ,DY,DYSQ,P00,P01,P02,P03,P10,P11,
     +                 P12,P13,P20,P21,P22,P23,P30,P31,P32,P33,Q0,Q1,Q2,
     +                 Q3,U,V,X0,XII,Y0,YII,Z00,Z01,Z0DX,Z0DY,Z10,Z11,
     +                 Z1DX,Z1DY,ZDXDY,ZII,ZX00,ZX01,ZX0DY,ZX10,ZX11,
     +                 ZX1DY,ZXY00,ZXY01,ZXY10,ZXY11,ZY00,ZY01,ZY0DX,
     +                 ZY10,ZY11,ZY1DX
      INTEGER          IIP,IXD0,IXD1,IXDI,IXDIPV,IYD0,IYD1,IYDI,IYDIPV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MAX
*     ..
* Calculation
* Outermost DO-loop with respect to the output point
      DO 10 IIP = 1,NIP
          XII = XI(IIP)
          YII = YI(IIP)
          IF (IIP.EQ.1) THEN
              IXDIPV = -1
              IYDIPV = -1
          ELSE
              IXDIPV = IXDI
              IYDIPV = IYDI
          END IF
          IXDI = INXI(IIP)
          IYDI = INYI(IIP)
* Retrieves the z and partial derivative values at the origin of
* the coordinate for the rectangle.
          IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
              IXD0 = MAX(1,IXDI)
              IYD0 = MAX(1,IYDI)
              X0 = XD(IXD0)
              Y0 = YD(IYD0)
              Z00 = ZD(IXD0,IYD0)
              ZX00 = PDD(1,IXD0,IYD0)
              ZY00 = PDD(2,IXD0,IYD0)
              ZXY00 = PDD(3,IXD0,IYD0)
          END IF
* Case 1.  When the rectangle is inside the data area in both the
* x and y directions.
          IF ((IXDI.GT.0.AND.IXDI.LT.NXD) .AND.
     +        (IYDI.GT.0.AND.IYDI.LT.NYD)) THEN
* Retrieves the z and partial derivative values at the other three
* vertexes of the rectangle.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  IXD1 = IXD0 + 1
                  DX = XD(IXD1) - X0
                  DXSQ = DX*DX
                  IYD1 = IYD0 + 1
                  DY = YD(IYD1) - Y0
                  DYSQ = DY*DY
                  Z10 = ZD(IXD1,IYD0)
                  Z01 = ZD(IXD0,IYD1)
                  Z11 = ZD(IXD1,IYD1)
                  ZX10 = PDD(1,IXD1,IYD0)
                  ZX01 = PDD(1,IXD0,IYD1)
                  ZX11 = PDD(1,IXD1,IYD1)
                  ZY10 = PDD(2,IXD1,IYD0)
                  ZY01 = PDD(2,IXD0,IYD1)
                  ZY11 = PDD(2,IXD1,IYD1)
                  ZXY10 = PDD(3,IXD1,IYD0)
                  ZXY01 = PDD(3,IXD0,IYD1)
                  ZXY11 = PDD(3,IXD1,IYD1)
* Calculates the polynomial coefficients.
                  Z0DX = (Z10-Z00)/DX
                  Z1DX = (Z11-Z01)/DX
                  Z0DY = (Z01-Z00)/DY
                  Z1DY = (Z11-Z10)/DY
                  ZX0DY = (ZX01-ZX00)/DY
                  ZX1DY = (ZX11-ZX10)/DY
                  ZY0DX = (ZY10-ZY00)/DX
                  ZY1DX = (ZY11-ZY01)/DX
                  ZDXDY = (Z1DY-Z0DY)/DX
                  A = ZDXDY - ZX0DY - ZY0DX + ZXY00
                  B = ZX1DY - ZX0DY - ZXY10 + ZXY00
                  C = ZY1DX - ZY0DX - ZXY01 + ZXY00
                  D = ZXY11 - ZXY10 - ZXY01 + ZXY00
                  P00 = Z00
                  P01 = ZY00
                  P02 = (2.0* (Z0DY-ZY00)+Z0DY-ZY01)/DY
                  P03 = (-2.0*Z0DY+ZY01+ZY00)/DYSQ
                  P10 = ZX00
                  P11 = ZXY00
                  P12 = (2.0* (ZX0DY-ZXY00)+ZX0DY-ZXY01)/DY
                  P13 = (-2.0*ZX0DY+ZXY01+ZXY00)/DYSQ
                  P20 = (2.0* (Z0DX-ZX00)+Z0DX-ZX10)/DX
                  P21 = (2.0* (ZY0DX-ZXY00)+ZY0DX-ZXY10)/DX
                  P22 = (3.0* (3.0*A-B-C)+D)/ (DX*DY)
                  P23 = (-6.0*A+2.0*B+3.0*C-D)/ (DX*DYSQ)
                  P30 = (-2.0*Z0DX+ZX10+ZX00)/DXSQ
                  P31 = (-2.0*ZY0DX+ZXY10+ZXY00)/DXSQ
                  P32 = (-6.0*A+3.0*B+2.0*C-D)/ (DXSQ*DY)
                  P33 = (2.0* (2.0*A-B-C)+D)/ (DXSQ*DYSQ)
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V* (P01+V* (P02+V*P03))
              Q1 = P10 + V* (P11+V* (P12+V*P13))
              Q2 = P20 + V* (P21+V* (P22+V*P23))
              Q3 = P30 + V* (P31+V* (P32+V*P33))
              ZII = Q0 + U* (Q1+U* (Q2+U*Q3))
* End of Case 1
* Case 2.  When the rectangle is inside the data area in the x
* direction but outside in the y direction.
          ELSE IF ((IXDI.GT.0.AND.IXDI.LT.NXD) .AND.
     +             (IYDI.LE.0.OR.IYDI.GE.NYD)) THEN
* Retrieves the z and partial derivative values at the other
* vertex of the semi-infinite rectangle.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  IXD1 = IXD0 + 1
                  DX = XD(IXD1) - X0
                  DXSQ = DX*DX
                  Z10 = ZD(IXD1,IYD0)
                  ZX10 = PDD(1,IXD1,IYD0)
                  ZY10 = PDD(2,IXD1,IYD0)
                  ZXY10 = PDD(3,IXD1,IYD0)
* Calculates the polynomial coefficients.
                  Z0DX = (Z10-Z00)/DX
                  ZY0DX = (ZY10-ZY00)/DX
                  P00 = Z00
                  P01 = ZY00
                  P10 = ZX00
                  P11 = ZXY00
                  P20 = (2.0* (Z0DX-ZX00)+Z0DX-ZX10)/DX
                  P21 = (2.0* (ZY0DX-ZXY00)+ZY0DX-ZXY10)/DX
                  P30 = (-2.0*Z0DX+ZX10+ZX00)/DXSQ
                  P31 = (-2.0*ZY0DX+ZXY10+ZXY00)/DXSQ
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V*P01
              Q1 = P10 + V*P11
              Q2 = P20 + V*P21
              Q3 = P30 + V*P31
              ZII = Q0 + U* (Q1+U* (Q2+U*Q3))
* End of Case 2
* Case 3.  When the rectangle is outside the data area in the x
* direction but inside in the y direction.
          ELSE IF ((IXDI.LE.0.OR.IXDI.GE.NXD) .AND.
     +             (IYDI.GT.0.AND.IYDI.LT.NYD)) THEN
* Retrieves the z and partial derivative values at the other
* vertex of the semi-infinite rectangle.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  IYD1 = IYD0 + 1
                  DY = YD(IYD1) - Y0
                  DYSQ = DY*DY
                  Z01 = ZD(IXD0,IYD1)
                  ZX01 = PDD(1,IXD0,IYD1)
                  ZY01 = PDD(2,IXD0,IYD1)
                  ZXY01 = PDD(3,IXD0,IYD1)
* Calculates the polynomial coefficients.
                  Z0DY = (Z01-Z00)/DY
                  ZX0DY = (ZX01-ZX00)/DY
                  P00 = Z00
                  P01 = ZY00
                  P02 = (2.0* (Z0DY-ZY00)+Z0DY-ZY01)/DY
                  P03 = (-2.0*Z0DY+ZY01+ZY00)/DYSQ
                  P10 = ZX00
                  P11 = ZXY00
                  P12 = (2.0* (ZX0DY-ZXY00)+ZX0DY-ZXY01)/DY
                  P13 = (-2.0*ZX0DY+ZXY01+ZXY00)/DYSQ
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V* (P01+V* (P02+V*P03))
              Q1 = P10 + V* (P11+V* (P12+V*P13))
              ZII = Q0 + U*Q1
* End of Case 3
* Case 4.  When the rectangle is outside the data area in both the
* x and y direction.
          ELSE IF ((IXDI.LE.0.OR.IXDI.GE.NXD) .AND.
     +             (IYDI.LE.0.OR.IYDI.GE.NYD)) THEN
* Calculates the polynomial coefficients.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  P00 = Z00
                  P01 = ZY00
                  P10 = ZX00
                  P11 = ZXY00
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V*P01
              Q1 = P10 + V*P11
              ZII = Q0 + U*Q1
          END IF
* End of Case 4
          ZI(IIP) = ZII
   10 CONTINUE
      RETURN
      END
