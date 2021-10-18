!*==ADDPNT_.spg  processed by SPAG 6.72Dc at 17:56 on  5 Aug 2021
      SUBROUTINE ADDPNT(X,Y,Ld,N,Xnew,Ynew,Ierr)
 
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!=  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!=  Y    - REAL vector of length LD, y-values                            (IO)=*
!=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
!=         program                                                           =*
!=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
!=         N < LD.  On exit, N is incremented by 1.                          =*
!=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!=  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*
 
      IMPLICIT NONE
!*--ADDPNT_21
 
! calling parameters
 
      INTEGER Ld , N
      REAL*8 X(Ld) , Y(Ld)
      REAL*8 Xnew , Ynew
      INTEGER Ierr
 
! local variables
 
      INTEGER insert
      INTEGER i
 
!-----------------------------------------------------------------------
 
! initialize error flag
 
      Ierr = 0
 
! check n<ld to make sure x will hold another point
 
      IF ( N.GE.Ld ) THEN
         WRITE (0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE (0,*) '                        All elements used.'
         Ierr = 1
      ENDIF
 
      insert = 1
      i = 2
 
! check, whether x is already sorted.
! also, use this loop to find the point at which xnew needs to be inserted
! into vector x, if x is sorted.
 
 100  IF ( i.LT.N ) THEN
         IF ( X(i).LT.X(i-1) ) THEN
            WRITE (0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//     &
                       &'in ascending order!' , i , X(i) , X(i-1)
            Ierr = 2
         ELSE
            IF ( Xnew.GT.X(i) ) insert = i + 1
         ENDIF
         i = i + 1
         GOTO 100
      ENDIF
 
! if <xnew,ynew> needs to be appended at the end, just do so,
! otherwise, insert <xnew,ynew> at position INSERT
 
      IF ( Xnew.GT.X(N) ) THEN
 
         X(N+1) = Xnew
         Y(N+1) = Ynew
 
      ELSE
 
! shift all existing points one index up
 
         DO i = N , insert , -1
            X(i+1) = X(i)
            Y(i+1) = Y(i)
         ENDDO
 
! insert new point
 
         X(insert) = Xnew
         Y(insert) = Ynew
 
      ENDIF
 
! increase total number of elements in x, y
 
      N = N + 1
 
      END
!*==INTER2_.spg  processed by SPAG 6.72Dc at 17:56 on  5 Aug 2021
 
      SUBROUTINE INTER2(Ng,Xg,Yg,N,X,Y,Ierr)
 
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points onto a set of target     =*
!=  bins.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!=  of target bins.  In general, this is the case for most of the weighting  =*
!=  functions (action spectra, molecular cross section, and quantum yield    =*
!=  data), which have to be matched onto the specified wavelength intervals. =*
!=  The average value in each target bin is found by averaging the trapezoi- =*
!=  dal area underneath the input data curve (constructed by linearly connec-=*
!=  ting the discrete input values).                                         =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data set does not span the range of the target grid, an error      =*
!=  message is printed and the execution is stopped, as extrapolation of the =*
!=  data is not permitted.                                                   =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - INTEGER, number of points in input grid                         (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*
 
      IMPLICIT NONE
!*--INTER2_132
 
! input:
      INTEGER Ng , N
      REAL*8 X(N) , Y(N) , Xg(Ng)
      INTEGER Ierr
! output:
      REAL*8 Yg(Ng)
 
! local:
      REAL*8 area , xgl , xgu
      REAL*8 darea , slope
      REAL*8 a1 , a2 , b1 , b2
      INTEGER ngintv
      INTEGER i , k , jstart
 
      Ierr = 0
!_______________________________________________________________________
 
!  test for correct ordering of data, by increasing value of x
 
      DO i = 2 , N
         IF ( X(i).LE.X(i-1) ) THEN
            Ierr = 1
            WRITE (*,*) 'data not sorted' , i , X(i) , X(i-1)
            RETURN
         ENDIF
      ENDDO
 
 
      DO i = 2 , Ng
         IF ( Xg(i).LE.Xg(i-1) ) THEN
            Ierr = 2
            WRITE (0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
            RETURN
         ENDIF
      ENDDO
 
 
 
 
! check for xg-values outside the x-range
 
      IF ( (X(1).GT.Xg(1)) .OR. (X(N).LT.Xg(Ng)) ) THEN
         WRITE (0,*) '>>> ERROR (inter2) <<<  Data do not span '//      &
                    &'grid.  '
         WRITE (0,*) '                        Use ADDPNT to '//         &
                    &'expand data and re-run.'
         RETURN
      ENDIF
 
!  find the integral of each grid interval and use this to
!  calculate the average y value for the interval
!  xgl and xgu are the lower and upper limits of the grid interval
 
      jstart = 1
      ngintv = Ng - 1
      DO i = 1 , ngintv
 
! initalize:
 
         area = 0.0
         xgl = Xg(i)
         xgu = Xg(i+1)
 
!  discard data before the first grid interval and after the
!  last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval
 
         k = jstart
 
         IF ( k.LE.N-1 ) THEN
 
!  if both points are before the first grid, go to the next point
 20         IF ( X(k+1).LE.xgl ) THEN
               jstart = k - 1
               k = k + 1
               IF ( k.LE.N-1 ) GOTO 20
            ENDIF
 
!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
 40         IF ( (k.LE.N-1) .AND. (X(k).LT.xgu) ) THEN
 
               jstart = k - 1
 
! compute x-coordinates of increment
 
               a1 = MAX(X(k),xgl)
               a2 = MIN(X(k+1),xgu)
 
!  if points coincide, contribution is zero
 
               IF ( X(k+1).EQ.X(k) ) THEN
                  darea = 0.E0
               ELSE
                  slope = (Y(k+1)-Y(k))/(X(k+1)-X(k))
                  b1 = Y(k) + slope*(a1-X(k))
                  b2 = Y(k) + slope*(a2-X(k))
                  darea = (a2-a1)*(b2+b1)/2.
!                       print *,a2,a1,k,y(k),slope,b2,b1,darea
               ENDIF
 
 
!  find the area under the trapezoid from a1 to a2
 
               area = area + darea
 
! go to next point
 
               k = k + 1
               GOTO 40
 
            ENDIF
 
         ENDIF
 
!  calculate the average y after summing the areas in the interval
         Yg(i) = area/(xgu-xgl)
 
 
      ENDDO
!_______________________________________________________________________
 
      END
