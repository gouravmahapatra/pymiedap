* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.
      SUBROUTINE RGLCTN(XD,YD,XI,YI,INXI,INYI,NXD,NYD,NIP)

Cf2py intent(in) XD,YD,XI,YI,NXD,NYD,NIP
Cf2py intent(out) INXI,INYI

*----------------------------------------------------------------------------
*     Implementation of RGLCTN in pymiedap.
*
*     Implemented by: Ashwyn Groot
*     Date: November 2018
*----------------------------------------------------------------------------
*
* Location of the desired points in a rectangular grid
* (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08
*
* This subroutine locates the desired points in a rectangular grid
* in the x-y plane.
*
* The grid lines can be unevenly spaced.
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
*   NIP  = number of the output points to be located (must be
*          1 or greater),
*   XI   = array of dimension NIP containing the x coordinates
*          of the output points to be located,
*   YI   = array of dimension NIP containing the y coordinates
*          of the output points to be located.
*
* The output arguments are
*   INXI = integer array of dimension NIP where the interval
*          numbers of the XI array elements are to be stored,
*   INYI = integer array of dimension NIP where the interval
*          numbers of the YI array elements are to be stored.
* The interval numbers are between 0 and NXD and between 0 and NYD,
* respectively.
*
*
* Specification statements
*     .. Scalar Arguments ..
      INTEGER          NIP,NXD,NYD
*     ..
*     .. Array Arguments ..
      REAL*8           XD(NXD),XI(NIP),YD(NYD),YI(NIP)
      INTEGER          INXI(NIP),INYI(NIP)
*     ..
*     .. Local Scalars ..
      REAL*8           XII,YII
      INTEGER          IIP,IMD,IMN,IMX,IXD,IYD,NINTX,NINTY
*     ..
* DO-loop with respect to IIP, which is the point number of the
* output point
      DO 30 IIP = 1,NIP
          XII = XI(IIP)
          YII = YI(IIP)
* Checks if the x coordinate of the IIPth output point, XII, is
* in a new interval.  (NINTX is the new-interval flag.)
          IF (IIP.EQ.1) THEN
              NINTX = 1
          ELSE
              NINTX = 0
              IF (IXD.EQ.0) THEN
                  IF (XII.GT.XD(1)) NINTX = 1
              ELSE IF (IXD.LT.NXD) THEN
                  IF ((XII.LT.XD(IXD)) .OR.
     +                (XII.GT.XD(IXD+1))) NINTX = 1
              ELSE
                  IF (XII.LT.XD(NXD)) NINTX = 1
              END IF
          END IF
* Locates the output point by binary search if XII is in a new
* interval.  Determines IXD for which XII lies between XD(IXD)
* and XD(IXD+1).
          IF (NINTX.EQ.1) THEN
              IF (XII.LE.XD(1)) THEN
                  IXD = 0
              ELSE IF (XII.LT.XD(NXD)) THEN
                  IMN = 1
                  IMX = NXD
                  IMD = (IMN+IMX)/2
   10             IF (XII.GE.XD(IMD)) THEN
                      IMN = IMD
                  ELSE
                      IMX = IMD
                  END IF
                  IMD = (IMN+IMX)/2
                  IF (IMD.GT.IMN) GO TO 10
                  IXD = IMD
              ELSE
                  IXD = NXD
              END IF
          END IF
          INXI(IIP) = IXD
* Checks if the y coordinate of the IIPth output point, YII, is
* in a new interval.  (NINTY is the new-interval flag.)
          IF (IIP.EQ.1) THEN
              NINTY = 1
          ELSE
              NINTY = 0
              IF (IYD.EQ.0) THEN
                  IF (YII.GT.YD(1)) NINTY = 1
              ELSE IF (IYD.LT.NYD) THEN
                  IF ((YII.LT.YD(IYD)) .OR.
     +                (YII.GT.YD(IYD+1))) NINTY = 1
              ELSE
                  IF (YII.LT.YD(NYD)) NINTY = 1
              END IF
          END IF
* Locates the output point by binary search if YII is in a new
* interval.  Determines IYD for which YII lies between YD(IYD)
* and YD(IYD+1).
          IF (NINTY.EQ.1) THEN
              IF (YII.LE.YD(1)) THEN
                  IYD = 0
              ELSE IF (YII.LT.YD(NYD)) THEN
                  IMN = 1
                  IMX = NYD
                  IMD = (IMN+IMX)/2
   20             IF (YII.GE.YD(IMD)) THEN
                      IMN = IMD
                  ELSE
                      IMX = IMD
                  END IF
                  IMD = (IMN+IMX)/2
                  IF (IMD.GT.IMN) GO TO 20
                  IYD = IMD
              ELSE
                  IYD = NYD
              END IF
          END IF
          INYI(IIP) = IYD
   30 CONTINUE
      RETURN
      END
