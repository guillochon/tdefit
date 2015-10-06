      subroutine splt (x, y, yp, nxy, t_o, nt, tsorted,  &
                       sigma, spl_o)

! OUTLINE:
! "Tense" cubic spline interpolator, adapted from IDL.  
! Normally set sigma to 1.0.
! Converted to Fortran by Marshall McCall, 99 June, using as
! a reference the original Fortran source code from FITPACK
! from which the IDL routine was derived.
! Note subtle differences between FITPACK and IDL code, particularly
! in setting initial approximation for endpoint slopes for nxy=2.
! Also, the routine does not allow sigma to be zero (formulae
! in that extreme were omitted).
! Adjustments were made to check if input arrays are sorted.

! Extrapolations are not allowed.  If it is required, the
! result is set to the value of the nearest endpoint.

! Checked by reproducing reddening curve of Fitzpatrick.

! 2005 Aug 12:  Fixed error in dimensioning of subs
! 2004 Jun 10:  Set dimensions to be determined by argument list

! $Id: spline.pro,v 1.6 1998/01/15 18:44:11 scottm Exp $

 
! INPUTS:
!	x:    The abcissa vector, sorted into ascending order.
!	y:    The vector of ordinate values corresponding to x.
!	yp:   The vector of second derivatives at each x.
!       nxy:    The number of elements in x, y, and yp.
!	t_o:    The vector of abcissae values for which the ordinate is 
!		desired. 
!       nt:     The number of elements of t_o.
!       tsorted:  A Boolean describing whether or not t has been
!                 sorted into a monotonicallly ascending order.
!	sigma:  The amount of "tension" that is applied to the curve. 
!               The default value is 1.0. If sigma is close to 0, 
!               (e.g., .01), then effectively there is a cubic spline 
!               fit. If sigma is large, (e.g., greater than 10), then 
!               the fit will be like a polynomial interpolation.
 
! OUTPUTS:
!	spline returns a vector of interpolated ordinates.
!	spl_o(i) = value of the function at t_o(i).
 
! RESTRICTIONS:
!       Must run splt_p first to get x, y, and yp from x_o and y_o.
 

      logical tsorted
      integer nt, nxy
      integer eflag, i, indt(nt), &
              j, jhi, jlo, k, khi, klo,  &
              nt_inter, sub, sub1, subs(nt)
      real*8 del1, del2, dels, dels1, dels2, ds1ds2,  &
             dsinh, sigma, sigmap, sinhd1, sinhd2, sinhs
      real*8 spl(nt), spl_o(nt), &
             t(nt), t_o(nt), &
             x(nxy), y(nxy), yp(nxy)

!     Check that there is more than one point in data arrays

      if (nxy .le. 1) then
         write (*,*) 'Error in splt:  Too few points to interpolate'
         return
      end if

!     Check that the points to interpolate are sorted into 
!     ascending order.  If not, sort them.

      if (.not. tsorted) then
         call hpsort (t_o, nt, indt, t)
      else
         do i=1,nt
            t(i) = t_o(i)
         end do
      end if

!     Check if extrapolations are required.  Set spl to nearest
!     endpoint wherever they are.  If start and end points to
!     interpolate are identical to endpoints of interpolation array,
!     set their values to the endpoints too.  This is necessary
!     to avoid an endpoint problem when subscripts of next highest
!     value in interpolation array are assigned.

!     nt_inter = number of points to interpolate

      jlo = 1
      jhi = nt
      eflag = 0
      do j=1,nt
         if (x(1) .ge. t(j)) then
            eflag = eflag + 1
            spl(j) = y(1)
            subs(j) = 1
            jlo = j + 1
         else
            go to 2
         end if
      end do
    2 do j=nt,1,-1
         if (x(nxy) .le. t(j)) then
            eflag = eflag + 1
            spl(j) = y(nxy)
            subs(j) = nxy + 1
            jhi = j - 1
         else
            go to 4
         end if
      end do
    4 nt_inter = jhi - jlo + 1
      if (nt_inter .lt. 1) go to 98

!     Find subscript khi where x(subs) > t(jlo) > x(subs-1).
!     This is done by the method of bisection.

      klo = 1
      khi = nxy
    6 if (khi - klo .gt. 1) then
         k = (khi + klo) / 2
         if (x(k) .gt. t(jlo)) then
           khi = k
         else
           klo = k
         end if
         go to 6
      end if
      subs(jlo) = khi
 
!     Find subscripts where x(subs) > t(j) > x(subs-1)
!     Remember that both x and t are sorted into ascending order
 
      if (nt_inter .eq. 1) go to 10
      j = jlo + 1
      do i=khi,nxy
    8    if (x(i) .gt. t(j)) then
            subs(j) = i
            j = j + 1
            if (j .gt. jhi) go to 10 
            go to 8
         end if
      end do

!     interpolate

   10 sigmap = sigma * (nxy - 1) / (x(nxy) - x(1))
      do i=jlo,jhi
         sub = subs(i)
         sub1 = subs(i) - 1
         del1 = t(i) - x(sub1)
         del2 = x(sub) - t(i)
         dels = x(sub) - x(sub1)
         dels1 = sigmap * del1
         sinhd1 = dsinh(dels1)
         dels2 = sigmap * del2
         sinhd2 = dsinh(dels2)
         ds1ds2 = dels1 + dels2
         sinhs = dsinh(ds1ds2)
         spl(i) = (yp(sub) * sinhd1 + yp(sub1) * sinhd2) / sinhs + &
               ((y(sub) - yp(sub)) * del1 +  &
               (y(sub1) - yp(sub1)) * del2) / dels
      end do

!     Put answers in input order

   98 if (.not. tsorted) then
         do i=1,nt
            spl_o(indt(i)) = spl(i)
         end do
      else
         do i=1,nt
            spl_o(i) = spl(i)
         end do
      end if

      return
      end
