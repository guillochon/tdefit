! Support file for redlaws, released under GPL license [http://adsabs.harvard.edu/abs/2004AJ....128.2144M]

      subroutine splt_p (x_o, y_o, nxy, xsorted, sigma, x, y, yp)

! (out)LINE:
! Subroutine to compute first derivatives for the
! "Tense" cubic spline interpolator adapted from IDL.  
! Normally set sigma to 1.0.
! Converted to Fortran by Marshall McCall, 99 June, using as
! a reference the original Fortran source code from FITPACK
! from which the IDL routine was derived.
! Note subtle differences between FITPACK and IDL code, particularly
! in setting initial approximation for endpoint slopes for nxy=2.
! Also, the routine does not allow sigma to be zero (formulae
! in that extreme were omitted).
! Adjustments were made to check if input arrays are sorted.

! 2004 Jun 10: Set dimensions to be determined by argument list
! 2003 May 9:  Fixed bug for situation where only 2 points
! Checked by reproducing reddening curve of Fitzpatrick.

! $Id: spline.pro,v 1.6 1998/01/15 18:44:11 scottm Exp $

 
! INPUTS:
!	x_o:    The abcissa vector, in original order.
!	y_o:    The vector of ordinate values corresponding to x_o.
!       nxy:    The number of elements in x_o and y_o.
!	sigma:  The amount of "tension" that is applied to the curve. 
!               The default value is 1.0. If sigma is close to 0, 
!               (e.g., .01), then effectively there is a cubic spline 
!               fit. If sigma is large, (e.g., greater than 10), then 
!               the fit will be like a polynomial interpolation.

! (out)PUTS:
!	x:      The abscissa vector, in ascending order.
!	y:      The vector of ordinate values corresponding to x.
!	yp:     The vector of second derivatives at each x.
 
! RESTRICTIONS:
!	Abcissa values, once sorted, must be monotonically increasing.
 
      logical xsorted
      integer nxy
      integer i, ibak, indx(nxy)
      real*8 c1, c2, c3, dcosh, deln, delnm1, delnn, dels, &
             delx1, delx12, delx2, diag1, diag2, diagin,  &
             dsinh, dx1, dx2, &
             sigma, sigmap, sinhin, sinhs,  &
             slpp1, slppn, spdiag
      real*8 temp(nxy), x(nxy), x_o(nxy),  &
             y(nxy), y_o(nxy), yp(nxy)

!     Check that there is more than one point in data arrays

      if (nxy .le. 1) then
         write (*,*) 'Error in splt_p:  Too few points to interpolate'
         return
      end if

!     Check that the data array is sorted into ascending order
!     If not, sort it, and load x and y.

      if (.not. xsorted) then
         call hpsort (x_o, nxy, indx, x)
         do i=1,nxy
            y(i) = y_o(indx(i))
         end do
      else
         do i=1,nxy
            x(i) = x_o(i)
            y(i) = y_o(i)
         end do
      end if

!     Approximate slopes at endpoints of x.
!     Do so by setting the tension to be zero.
!     If nxy=2, then set approximations to the slopes to zero.

!     slpp1 = slope at x(1)
!     slppn = slope at x(nxy)

      delx1 = x(2) - x(1)
      if (delx1 .eq. 0.) then
         write (*,4000) x(1), x(2)
 4000    format('In splt_p, x(1) is same as x(2):  ',  &
                 d14.7, 1x, d14.7)
         return
      end if
 
      if (nxy .le. 2) then
         yp(1) = 0.
         yp(2) = 0.
      else
         dx1 = (y(2) - y(1)) / delx1
         delx2 = x(3) - x(2)
         if (delx2 .eq. 0.) then
            write (*,5000) x(2), x(3)
 5000       format ('In splt_p, x(2) is same as x(3):  ',  &
                    d14.7, 1x, d14.7)
            return
         end if
         delx12 = x(3) - x(1)
         c1 = -(delx12 + delx1) / (delx12 * delx1)
         c2 = delx12 / (delx1 * delx2)
         c3 = -delx1 / (delx12 * delx2)
         slpp1 = c1 * y(1) + c2 *y(2) + c3 * y(3)
         deln = x(nxy) - x(nxy-1)
         delnm1 = x(nxy-1) - x(nxy-2)
         delnn = x(nxy) - x(nxy-2)
         c1 = (delnn + deln) / (delnn * deln)
         c2 = -delnn / (deln * delnm1)
         c3 = deln / (delnn * delnm1)
         slppn = c3 * y(nxy-2) + c2 * y(nxy-1) + c1 * y(nxy)

!        Revise slopes at endpoints
!        First denormalize tension factor
!        sigmap = denormalized value of sigma
 
         sigmap = sigma * (nxy - 1) / (x(nxy) - x(1))
         dels = sigmap * delx1
         sinhs = dsinh(dels)
         sinhin = 1. / (delx1 * sinhs)
         diag1 = sinhin * (dels * dcosh(dels) - sinhs)
         diagin = 1. / diag1
         yp(1) = diagin * (dx1 - slpp1)
         spdiag = sinhin * (sinhs - dels)
         temp(1) = diagin * spdiag

!        Now compute slopes in between endpoints
 
         do i=2,nxy-1
            delx2 = x(i+1) - x(i)
            if (delx2 .eq. 0.) then
               write (*,6000) x(i), x(i+1)
 6000        format ('In splt_p, x not strictly increasing:  ',  &
                     d14.7, 1x, d14.7)
               return
            end if
            dx2=(y(i+1) - y(i)) / delx2
            dels = sigmap * delx2
            sinhs = dsinh(dels)
            sinhin = 1. / (delx2 * sinhs)
            diag2 = sinhin * (dels * dcosh(dels) - sinhs)
            diagin = 1. / (diag1 + diag2 - spdiag * temp(i-1))
            yp(i) = diagin * (dx2 - dx1 - spdiag * yp(i-1))
            spdiag = sinhin * (sinhs - dels)
            temp(i) = diagin * spdiag
            dx1 = dx2
            diag1 = diag2
         end do
 
         diagin = 1. / (diag1 - spdiag * temp(nxy-1))
         yp(nxy) = diagin * (slppn - dx1 - spdiag * yp(nxy-1))

         do i=2,nxy
            ibak = nxy + 1 - i
            yp(ibak) = yp(ibak) - temp(ibak) * yp(ibak+1)
         end do

      end if

      return
      end
