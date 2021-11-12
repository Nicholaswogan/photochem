      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew , ierr)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
*=  ascending order                                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
*=  Y    - REAL vector of length LD, y-values                            (IO)=*
*=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
*=         program                                                           =*
*=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
*=         N < LD.  On exit, N is incremented by 1.                          =*
*=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
*=  YNEW - REAL, y-value of point to be added                             (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

C calling parameters

      INTEGER ld, n
      REAL*8 x(ld), y(ld)
      REAL*8 xnew, ynew
      INTEGER ierr

C local variables

      INTEGER insert
      INTEGER i

C-----------------------------------------------------------------------

* initialize error flag

      ierr = 0

* check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
         WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE(0,*) '                        All elements used.'
         ierr = 1
      ENDIF

      insert = 1
      i = 2

* check, whether x is already sorted.
* also, use this loop to find the point at which xnew needs to be inserted
* into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           WRITE(0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//
     >                'in ascending order!', i,x(i),x(i-1)
           ierr = 2
        ELSE
           IF (xnew .GT. x(i)) insert = i + 1
        ENDIF
        i = i+1
        GOTO 10
      ENDIF

* if <xnew,ynew> needs to be appended at the end, just do so,
* otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN

         x(n+1) = xnew
         y(n+1) = ynew

      ELSE

* shift all existing points one index up

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

* insert new point

         x(insert) = xnew
         y(insert) = ynew

      ENDIF

* increase total number of elements in x, y

      n = n+1

      END

      SUBROUTINE inter2(ng,xg,yg,n,x,y,ierr)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on single, discrete points onto a set of target     =*
*=  bins.                                                                    =*
*=  The original input data are given on single, discrete points of an       =*
*=  arbitrary grid and are being linearly interpolated onto a specified set  =*
*=  of target bins.  In general, this is the case for most of the weighting  =*
*=  functions (action spectra, molecular cross section, and quantum yield    =*
*=  data), which have to be matched onto the specified wavelength intervals. =*
*=  The average value in each target bin is found by averaging the trapezoi- =*
*=  dal area underneath the input data curve (constructed by linearly connec-=*
*=  ting the discrete input values).                                         =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data set does not span the range of the target grid, an error      =*
*=  message is printed and the execution is stopped, as extrapolation of the =*
*=  data is not permitted.                                                   =*
*=  If the input data does not encompass the target grid, use ADDPNT to      =*
*=  expand the input array.                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
*=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
*=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
*=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
*=        bin i (i = 1..NG-1)                                                =*
*=  N   - INTEGER, number of points in input grid                         (I)=*
*=  X   - REAL, grid on which input data are defined                      (I)=*
*=  Y   - REAL, input y-data                                              (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER ng, n
      REAL*8 x(n), y(n), xg(ng)
      INTEGER ierr
* output:
      REAL*8 yg(ng)

* local:
      REAL*8 area, xgl, xgu
      REAL*8 darea, slope
      REAL*8 a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart

      ierr = 0
*_______________________________________________________________________

*  test for correct ordering of data, by increasing value of x

      DO 10, i = 2, n
         IF (x(i) .LE. x(i-1)) THEN
            ierr = 1
            WRITE(*,*)'data not sorted', i,x(i),x(i-1)
            RETURN
         ENDIF
   10 CONTINUE


      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
           ierr = 2
          WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
          RETURN
        ENDIF
      ENDDO




* check for xg-values outside the x-range

      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
          WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span '//
     >               'grid.  '
          WRITE(0,*) '                        Use ADDPNT to '//
     >               'expand data and re-run.'
          return
      ENDIF

*  find the integral of each grid interval and use this to
*  calculate the average y value for the interval
*  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1
      ngintv = ng - 1
      DO 50, i = 1,ngintv

* initalize:

            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)

*  discard data before the first grid interval and after the
*  last grid interval
*  for internal grid intervals, start calculating area by interpolating
*  between the last point which lies in the previous interval and the
*  first point inside the current interval

            k = jstart

            IF (k .LE. n-1) THEN

*  if both points are before the first grid, go to the next point
   30         CONTINUE
                IF (x(k+1) .LE. xgl) THEN
                   jstart = k - 1
                   k = k+1
                   IF (k .LE. n-1) GO TO 30
                ENDIF

*  if the last point is beyond the end of the grid, complete and go to the next
*  grid
   40         CONTINUE
                 IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN

                    jstart = k-1

* compute x-coordinates of increment

                    a1 = MAX(x(k),xgl)
                    a2 = MIN(x(k+1),xgu)

*  if points coincide, contribution is zero

                    IF (x(k+1).EQ.x(k)) THEN
                       darea = 0.e0
                    ELSE
                       slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                       b1 = y(k) + slope*(a1 - x(k))
                       b2 = y(k) + slope*(a2 - x(k))
                       darea = (a2 - a1)*(b2 + b1)/2.
c                       print *,a2,a1,k,y(k),slope,b2,b1,darea
                    ENDIF


*  find the area under the trapezoid from a1 to a2

                    area = area + darea

* go to next point

                    k = k+1
                    GO TO 40

                ENDIF

            ENDIF

*  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)


   50 CONTINUE
*_______________________________________________________________________

      RETURN
      END
