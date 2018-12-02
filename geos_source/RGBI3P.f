* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.
      SUBROUTINE RGBI3P(MD,XD,YD,ZD,XI,YI,
     +            ZI,IER,filetype,WK,NXD,NYD,NIP)

Cf2py intent(in) MD,XD,YD,ZD,XI,YI,filetype,WK,NXD,NYD,NIP
Cf2py intent(out) ZI,IER

*----------------------------------------------------------------------------
*     Implementation of RGBI3P in pymiedap.
*
*     Editer: Ashwyn Groot
*     Date: November 2018
*----------------------------------------------------------------------------
*
* Rectangular-grid bivariate interpolation
* (a master subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08
*
* This subroutine performs interpolation of a bivariate function,
* z(x,y), on a rectangular grid in the x-y plane.  It is based on
* the revised Akima method.
*
* In this subroutine, the interpolating function is a piecewise
* function composed of a set of bicubic (bivariate third-degree)
* polynomials, each applicable to a rectangle of the input grid
* in the x-y plane.  Each polynomial is determined locally.
*
* This subroutine has the accuracy of a bicubic polynomial, i.e.,
* it interpolates accurately when all data points lie on a
* surface of a bicubic polynomial.
*
* The grid lines can be unevenly spaced.
*
* The input arguments are
*   MD  = mode of computation
*       = 1 for new XD, YD, or ZD data (default)
*       = 2 for old XD, YD, and ZD data,
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
*         points,
*   NIP = number of the output points at which interpolation
*         of the z value is desired (must be 1 or greater),
*   XI  = array of dimension NIP containing the x coordinates
*         of the output points,
*   YI  = array of dimension NIP containing the y coordinates
*         of the output points.
*
* The output arguments are
*   ZI  = array of dimension NIP where the interpolated z
*         values at the output points are to be stored,
*   IER = error flag
*       = 0 for no errors
*       = 1 for NXD = 1 or less
*       = 2 for NYD = 1 or less
*       = 3 for identical XD values or
*               XD values out of sequence
*       = 4 for identical YD values or
*               YD values out of sequence
*       = 5 for NIP = 0 or less.
*
* The other argument is
*   WK  = three dimensional array of dimension 3*NXD*NYD used
*         internally as a work area.
*
* The very fisrt call to this subroutine and the call with a new
* XD, YD, and ZD array must be made with MD=1.  The call with MD=2
* must be preceded by another call with the same XD, YD, and ZD
* arrays.  Between the call with MD=2 and its preceding call, the
* WK array must not be disturbed.
*
* The constant in the PARAMETER statement below is
*   NIPIMX = maximum number of output points to be processed
*            at a time.
* The constant value has been selected empirically.
*
* This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.

* Specification statements
*     .. Parameters ..
      INTEGER          NIPIMX
      PARAMETER        (NIPIMX=51)

*     .. Scalar Arguments ..
      INTEGER          IER,MD,NIP,NXD,NYD

*     .. Array Arguments ..
      REAL*8           WK(3,NXD,NYD),XD(NXD),XI(NIP),YD(NYD),YI(NIP),
     +                 ZD(NXD,NYD),ZI(NIP)

*     .. Local Scalars ..
      INTEGER          IIP,IX,IY,NIPI

*     .. Local Arrays ..
      INTEGER          INXI(NIPIMX),INYI(NIPIMX)

*     .. External Subroutines ..
C      EXTERNAL         RGLCTN,RGPD3P,RGPLNL

*     .. Intrinsic Functions ..
      INTRINSIC        MIN

* Preliminary processing
* Error check
      IF (NXD.LE.1) GO TO 40
      IF (NYD.LE.1) GO TO 50
      DO 10 IX = 2,NXD
          IF (XD(IX).LE.XD(IX-1)) GO TO 60
   10 CONTINUE
      DO 20 IY = 2,NYD
          IF (YD(IY).LE.YD(IY-1)) GO TO 70
   20 CONTINUE
      IF (NIP.LE.0) GO TO 80
      IER = 0
* Calculation
* Estimates partial derivatives at all input-grid data points
* (for MD=1).
      IF (filetype.GT.1) THEN
         IF (MD.NE.2) THEN
            WK=0.D0
            CALL RGPD3P(XD,YD,ZD,WK,NXD,NYD)
*          CALL RGPD3P(NXD,NYD,XD,YD,ZD, PDD)
         END IF
      ENDIF
* DO-loop with respect to the output point
* Processes NIPIMX output points, at most, at a time.
      DO 30 IIP = 1,NIP,NIPIMX
          NIPI = MIN(NIP- (IIP-1),NIPIMX)
* Locates the output points.
          CALL RGLCTN(XD,YD,XI(IIP),YI(IIP),
     +            INXI,INYI,NXD,NYD,NIPI)
*         CALL RGLCTN(NXD,NYD,XD,YD,NIP,XI,YI, INXI,INYI)
* Calculates the z values at the output points.
          CALL RGPLNL(XD,YD,ZD,WK,XI(IIP),
     +            YI(IIP),INXI,INYI,ZI(IIP),NXD,NYD,NIPI)
*         CALL RGPLNL(NXD,NYD,XD,YD,ZD,PDD,NIP,XI,YI,INXI,INYI, ZI)
   30 CONTINUE
      RETURN

* Error exit
   40 WRITE (*,FMT=9000)
      IER = 1
      GO TO 90
   50 WRITE (*,FMT=9010)
      IER = 2
      GO TO 90
   60 WRITE (*,FMT=9020) IX,XD(IX)
      IER = 3
      GO TO 90
   70 WRITE (*,FMT=9030) IY,YD(IY)
      IER = 4
      GO TO 90
   80 WRITE (*,FMT=9040)
      IER = 5
   90 WRITE (*,FMT=9050) NXD,NYD,NIP
      RETURN
* Format statements for error messages
 9000 FORMAT (1X,/,'*** RGBI3P Error 1: NXD = 1 or less')
 9010 FORMAT (1X,/,'*** RGBI3P Error 2: NYD = 1 or less')
 9020 FORMAT (1X,/,'*** RGBI3P Error 3: Identical XD values or',
     +       ' XD values out of sequence',/,'    IX =',I6,',  XD(IX) =',
     +       E11.3)
 9030 FORMAT (1X,/,'*** RGBI3P Error 4: Identical YD values or',
     +       ' YD values out of sequence',/,'    IY =',I6,',  YD(IY) =',
     +       E11.3)
 9040 FORMAT (1X,/,'*** RGBI3P Error 5: NIP = 0 or less')
 9050 FORMAT ('    NXD =',I5,',  NYD =',I5,',  NIP =',I5,/)
      END
