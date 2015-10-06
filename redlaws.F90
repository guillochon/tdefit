      subroutine law_fmr (x_o, nx, xsorted, R, a1_o)

!     2005 Aug 22:  Added law_odo (i.e., O'Donnell 1994)
!     2004 Jul 29:  Revised output formats in error messages
!     2004 JUN 18:  Constants converted to double precision
!     2003 AUG 11:  Made ! first character in error messages.

!----------------------------------------------------------------
!  COMPUTE R-DEPENDENT EXTINCTION CURVE FROM IR THROUGH OPTICAL
!
!  Uses Fitzpatrick/Massa fitting function in UV and cubic 
!  spline interpolation in the optical/IR region.
!
!  Approach discussed in Fitzpatrick (1999)
!
!  Given wave numbers x_o (in inverse microns) 
!  and the ratio of total to
!  selective extinction "R", this subroutine computes the
!  extinction coefficients a1 at each x_o, namely the ratio of
!  A_lambda / E(B-V).  Values of A_lambda are monochromatic,
!  but the E(B-V) suffers from bandpass effects.  The E(B-V)
!  applies to a 30,000 K star for which E(B-V)=0.5.
!
!----------------------------------------------------------------

!     INPUT:
!     x_o =  wave number, in inverse microns 
!     nx = number of elements in x
!     xsorted = Boolean defining whether or not x_o is sorted into
!               ascending order
!             = true, if x_o is already sorted
!             = false, if not
!     R = ratio of total-to-selective extinction

!     OUTPUT:
!     a1_o = extinction coefficient


!     The following specify array dimensions are needed for 
!     those variables in synext6 which aren't passed in subroutine calls

!     nfilpmax = maximum number of filters which can be picked
!     nflxmax = maximum number of input points in each object spectrum
!     nresmax = maximum number of input points in each filter response curve
!     nsubmax = maximum number of response curve samples
!     ntaumax = maximum number of values of tau1
!     nzmax = maximum number of values of z

      use tdefit_data, only: fitz_xspl, fitz_nspl, fitz_nspltot
      use constants, only: fitz_xcutuv, fitz_x0, fitz_gamma, fitz_c3, fitz_c4
      include 'synext6.inc'

!     Now the declarations

      logical sorted, xsorted
      integer flag, i, iuv, j, nopir, nuv, nx
      integer indx(nflxmax)
      real*8 a1(nflxmax), a1_o(nflxmax), c1, c2, &
             coespl(3), &
             R, tension, x(nflxmax), R2, &
             xspl_s(fitz_nspltot), xtouse, x_o(nflxmax), &
             xopir(nflxmax), xuv(nflxmax), &
             yopir(nflxmax), yspl(fitz_nspltot), yspl_s(fitz_nspltot), &
             ysplp_s(fitz_nspltot), yuv(nflxmax)

      flag = 0

!     Sort x_o if necessary.  Place in x.

      if (xsorted) then
         do i=1,nx
            indx(i) = i
            x(i) = x_o(i)
         end do
      else
         call hpsort (x_o, nx, indx, x)
      end if
      sorted = .true.

!     Now classify input wave numbers, based upon whether they
!     are red or blue of a cutoff xcutuv.

!     nopir = counter for optical/IR input wave numbers
!     nuv = counter for UV wave numbers (starts at 2, for
!           2700 and 2600 A, but can change if any input
!           wave numbers exceed xcutuv.
!     xcutuv = boundary between UV and optical/IR
!     xuv = values of x in UV
!     xopir = values of x in optical/IR

      nopir = 0
      nuv = 0

      do i=1,nx
         if (x(i) .lt. 0.) then
            x(i) = 0.0d0
            flag = 1
         end if
         if (x(i) .gt. 11.0) then
            x(i) = 11.0d0
            flag = 1
         end if
         if (x(i) .ge. fitz_xcutuv) then
            nuv = nuv + 1
            xuv(nuv) = x(i)
         else
            nopir = nopir + 1
            xopir(nopir) = x(i)
         end if
      end do
      if (flag .eq. 1) then
         write (*,1000) (10000.0d0 / x_o(indx(nx))), &
                        (10000.0d0 / x_o(indx(1)))
 1000    format(t1, &
                '! Input line ', i3, /, &
                t1, '! Warning:  Extinction coefficients are ', &
                'required from ', /, &
                t1, '! ', f10.1, ' to ', f10.1, ' Angstroms', /, &
                t1, '! At least one wavelength lies ', &
                'outside the bounds of applicability ', &
                /, t1, '! of the reddening law fitzpatrick_99', &
                /, t1,  '! Wavelengths must be longer than ', &
                '909.1 Angstroms.', &
                /, t1, '! Outside the bounds, extinction ', &
                'coefficients are set to be the same as for ', &
                /, t1, '! nearest wavelength encompassed by law.', &
                /, t1, '!')
      endif

!     Append the two UV spline anchor points to the end of xuv

      xuv(nuv+1) = fitz_xspl(8)
      xuv(nuv+2) = fitz_xspl(9)

! ............................................
! COMPUTE UV VALUES OF A(lambda)/E(B-V)
! USING FM FITTING FUNCTION AND R-DEPENDENT 
! COEFFICIENTS
! ............................................

!     Note that if nuv is zero, there are still 2 points
!     in the UV for which extinction coefficients must
!     be computed, namely the last two spline anchor points
!     required to evaluate extinction coefficients in the
!     optical/IR.

!     Now, the coefficients which correlate with R

      c2    = -0.824d0 + 4.717d0 / R
      c1    =  2.030d0 - 3.007d0 * c2

!     Now, compute the extinction coefficients in the UV,
!     remembering that the last two elements of xuv are the
!     spline anchor points.
!     yuv = extinction coefficient at xuv
 
      do iuv=1,nuv+2
         yuv(iuv) = c1  + c2 * xuv(iuv)
         yuv(iuv) = yuv(iuv) + fitz_c3 * xuv(iuv)**2 &
                    / ((xuv(iuv)**2-fitz_x0**2)**2 &
                    + (xuv(iuv)*fitz_gamma)**2)
         if (xuv(iuv) .gt. 5.9) then
            xtouse = xuv(iuv)
         else
            xtouse = 5.9d0
         end if
         yuv(iuv) = yuv(iuv) + fitz_c4 * (0.5392d0 * (xtouse - 5.9d0)**2 &
                   + 0.05644d0 * (xtouse - 5.9d0)**3)
         yuv(iuv) = yuv(iuv) + R
      end do
      
! .....................................................
! COMPUTE OPTICAL/IR VALUES OF A(lambda)/E(B-V)
! USING S CUBIC SPLINE ANCHORED IN UV, OPTICAL, and IR
! (too many significant figures here, but...)
! .....................................................

      if (nopir .ne. 0) then

!        Set the extinction coefficients at the spline anchor points
!        yspl = extinction coefficient at xspl

         data coespl /0.0d0, 0.26469d0, 0.82925d0/
         do i=1,3
            yspl(i)   = coespl(i) * R / 3.1d0
         end do
         R2 = R*R
         yspl(4)  = -4.22809d-01 + 1.00270d0*R + 2.13572d-04*R2
         yspl(5)  = -5.13540d-02 + 1.00216d0*R - 7.35778d-05*R2
         yspl(6)  =  7.00127d-01 + 1.00184d0*R - 3.32598d-05*R2
         yspl(7)  =  1.19456d0   + 1.01707d0*R - 5.46959d-03*R2 &
                           + 7.97809d-04*R*R*R - 4.45636d-05*R*R*R*R
         yspl(8)  = yuv(nuv+1)
         yspl(9)  = yuv(nuv+2)

!        Interpolate extinction coefficients for all xopir
!        Use a tension of 1.0 for the spline, as Fitzpatrick must have
!        yopir = extinction coeffient at xopir
  
         tension = 1.0d0
         call splt_p (fitz_xspl, yspl, fitz_nspl, sorted, tension,  &
                      xspl_s, yspl_s, ysplp_s)
         call splt (xspl_s, yspl_s, ysplp_s, fitz_nspl,  &
                    xopir, nopir, sorted, tension, yopir)
      end if
 
! ........................................................
! PUT ALL THE ANSWERS TOGETHER
! ........................................................

      if (nopir .ne. 0 .and. nuv .eq. 0) then
         do i=1,nx
            a1(i) = yopir(i)
         end do
      else if (nopir .eq. 0 .and. nuv .ne. 0) then
         do i=1,nx
            a1(i) = yuv(i)
         end do
      else if (nopir .ne. 0 .and. nuv .ne. 0) then
         do i=1,nopir
            a1(i) = yopir(i)
         end do
         do i=1,nuv
            j = nopir + i
            a1(j) = yuv(i)
         end do
      end if

!
!     Put a1 into correct order.
!

      do i=1,nx
         a1_o(indx(i)) = a1(i)
      end do

      return
      end


      subroutine law_ccm (x_o, nx, xsorted, R, a1_o)
!
!     Reddening law of Cardelli, Clayton, and Mathis 1989
!
!     A_lambda / A_V = a(x) + b(x) / R
!     A_lambda / E(B-V) = a(x) * R + b(x)
!
!     INPUT:
!     x_o =  wave number, in inverse microns 
!     nx = number of elements in x
!     xsorted = Boolean defining whether or not x_o is sorted into
!               ascending order
!             = true, if x_o is already sorted
!             = false, if not
!     R = ratio of total-to-selective extinction for OB stars
!         without regard to bandpass effects brought on by
!         the level of extinction

!     OUTPUT:
!     a1_o = extinction coefficient A_lambda/E(B-V)
 
!     Let
!     
!     a = a(x)
!     b = b(x)
!     x = 1 / lambda in microns^(-1)
!     a1_o = A_lambda / E(B-V)


!     The following specify array dimensions are needed for 
!     those variables in synext6 which aren't passed in subroutine calls

!     nfilpmax = maximum number of filters which can be picked
!     nflxmax = maximum number of input points in each object spectrum
!     nresmax = maximum number of input points in each filter response curve
!     nsubmax = maximum number of response curve samples
!     ntaumax = maximum number of values of tau1
!     nzmax = maximum number of values of z

      include 'synext6.inc'

!     Now the declarations

      logical sorted, xsorted
      integer flag, i, indx(nflxmax), nx
      real*8 a, a1_o(nflxmax), b, fa, fb, R,  &
             x, x_o(nflxmax), x_s(nflxmax), y

      flag = 0

!     Sort x_o if necessary.  Place in x.

      if (xsorted) then
         do i=1,nx
            indx(i) = i
            x_s(i) = x_o(i)
         end do
      else
         call hpsort (x_o, nx, indx, x_s)
      end if
      sorted = .true.

      do i=1,nx
         if (x_s(i) .lt. 0.) then
            x_s(i) = 0.0d0
            flag = 1
         else if (x_s(i) .gt. 10.0) then
            x_s(i) = 10.0d0
            flag = 1
         end if
         x = x_s(i)
         y = x - 1.82d0
!
!     First, 0.3 < x < 1.1 microns
!
         if (x .le. 1.1) then
            a = 0.574d0 * x**1.61d0
            b = -0.527d0 * x**1.61d0
         else
            if (x .le. 3.3) then
               a = 1 + 0.17699d0 * y - 0.50447d0 * y**2  &
                     - 0.02427d0 * y**3 + 0.72085d0 * y**4 &
                     + 0.01979d0 * y**5  &
                     - 0.77530d0 * y**6 + 0.32999d0 * y**7
               b = 1.41338d0 * y + 2.28305d0 * y**2 + 1.07233d0 * y**3 &
                     - 5.38434d0 * y**4 &
                     - 0.62251d0 * y**5 + 5.30260d0 * y**6  &
                     - 2.09002d0 * y**7
            else
               if (x .le. 8.0) then
                  if (x .lt. 5.9) then
                     fa = 0.0d0
                     fb = 0.0d0
                  else
                     fa = -0.04473d0 * (x - 5.9d0)**2  &
                          - 0.009779d0 * (x - 5.9d0)**3
                     fb = 0.2130d0 * (x - 5.9d0)**2  &
                          + 0.1207d0 * (x - 5.9d0)**3
                  end if
                     a = 1.752d0 - 0.316d0 * x  &
                         - 0.104d0 / ((x - 4.67d0)**2 + 0.341d0) + fa
                     b = -3.090d0 + 1.825d0 * x  &
                         + 1.206d0 / ((x - 4.62d0)**2 + 0.263d0) + fb
               else
 
!     note here that equations are specifically valid only out 
!     to 10 microns^(-1)

                  a = -1.073d0 - 0.628d0 * (x - 8)  &
                        + 0.137d0 * (x - 8)**2 - 0.070d0 * (x - 8)**3
                  b = 13.670d0 + 4.257d0 * (x - 8)  &
                      - 0.420d0 * (x - 8)**2 + 0.374d0 * (x - 8)**3
               end if
            end if
         end if

         a1_o(indx(i)) = a * R + b

      end do

      if (flag .eq. 1) then
         write (*,1000) (10000.0d0 / x_o(indx(nx))),  &
                        (10000.0d0 / x_o(indx(1)))
 1000    format(t1, &
                '! Input line ', i3, /, &
                t1, '! Warning:  Extinction coefficients are ', &
                'required from ', /, &
                t1, '! ', f10.1, ' to ', f10.1, ' Angstroms', /, &
                t1, '! At least one wavelength lies ', &
                'outside the bounds of applicability ', &
                /, t1, '! of the reddening law cardelli_89', &
                /, t1,  '! Wavelengths must be longer than ', &
                '1000.0 Angstroms.', &
                /, t1, '! Outside the bounds, extinction ', &
                'coefficients are set to be the same as for ', &
                /, t1, '! nearest wavelength encompassed by law.', &
                /, t1, '!')
      end if

      return
      end


      subroutine law_odo (x_o, nx, xsorted, R, a1_o)
!
!     Reddening law of O'Donnell 1994.  This is the same
!     as that of Cardelli, Clayton, and Mathis (1989), except
!     that formulae for the optical have been modified.
!
!     A_lambda / A_V = a(x) + b(x) / R
!     A_lambda / E(B-V) = a(x) * R + b(x)
!
!     INPUT:
!     x_o =  wave number, in inverse microns 
!     nx = number of elements in x
!     xsorted = Boolean defining whether or not x_o is sorted into
!               ascending order
!             = true, if x_o is already sorted
!             = false, if not
!     R = ratio of total-to-selective extinction for OB stars
!         without regard to bandpass effects brought on by
!         the level of extinction

!     OUTPUT:
!     a1_o = extinction coefficient A_lambda/E(B-V)
 
!     Let
!     
!     a = a(x)
!     b = b(x)
!     x = 1 / lambda in microns^(-1)
!     a1_o = A_lambda / E(B-V)


!     The following specify array dimensions are needed for 
!     those variables in synext6 which aren't passed in subroutine calls

!     nfilpmax = maximum number of filters which can be picked
!     nflxmax = maximum number of input points in each object spectrum
!     nresmax = maximum number of input points in each filter response curve
!     nsubmax = maximum number of response curve samples
!     ntaumax = maximum number of values of tau1
!     nzmax = maximum number of values of z

      include 'synext6.inc'

!     Now the declarations

      logical sorted, xsorted
      integer flag, i, indx(nflxmax), nx
      real*8 a, a1_o(nflxmax), b, fa, fb, R,  &
             x, x_o(nflxmax), x_s(nflxmax), y

      flag = 0

!     Sort x_o if necessary.  Place in x.

      if (xsorted) then
         do i=1,nx
            indx(i) = i
            x_s(i) = x_o(i)
         end do
      else
         call hpsort (x_o, nx, indx, x_s)
      end if
      sorted = .true.

      do i=1,nx
         if (x_s(i) .lt. 0.) then
            x_s(i) = 0.0d0
            flag = 1
         else if (x_s(i) .gt. 10.0) then
            x_s(i) = 10.0d0
            flag = 1
         end if
         x = x_s(i)
         y = x - 1.82d0
!
!     First, 0.3 < x < 1.1 microns
!
         if (x .le. 1.1) then
            a = 0.574d0 * x**1.61d0
            b = -0.527d0 * x**1.61d0
         else
            if (x .le. 3.3) then
               a = 1 + 0.104d0 * y - 0.609d0 * y**2  &
                     + 0.701d0 * y**3 + 1.137d0 * y**4 &
                     - 1.718d0 * y**5  &
                     - 0.827d0 * y**6 + 1.647d0 * y**7 &
                     - 0.505d0 * y**8
               b = 1.952d0 * y + 2.908d0 * y**2 - 3.989d0 * y**3 &
                     - 7.985d0 * y**4 &
                     + 11.102d0 * y**5 + 5.491d0 * y**6  &
                     - 10.805d0 * y**7 + 3.347 * y**8
            else
               if (x .le. 8.0) then
                  if (x .lt. 5.9) then
                     fa = 0.0d0
                     fb = 0.0d0
                  else
                     fa = -0.04473d0 * (x - 5.9d0)**2  &
                          - 0.009779d0 * (x - 5.9d0)**3
                     fb = 0.2130d0 * (x - 5.9d0)**2  &
                          + 0.1207d0 * (x - 5.9d0)**3
                  end if
                     a = 1.752d0 - 0.316d0 * x  &
                         - 0.104d0 / ((x - 4.67d0)**2 + 0.341d0) + fa
                     b = -3.090d0 + 1.825d0 * x  &
                         + 1.206d0 / ((x - 4.62d0)**2 + 0.263d0) + fb
               else
 
!     note here that equations are specifically valid only out 
!     to 10 microns^(-1)

                  a = -1.073d0 - 0.628d0 * (x - 8)  &
                        + 0.137d0 * (x - 8)**2 - 0.070d0 * (x - 8)**3
                  b = 13.670d0 + 4.257d0 * (x - 8)  &
                      - 0.420d0 * (x - 8)**2 + 0.374d0 * (x - 8)**3
               end if
            end if
         end if

         a1_o(indx(i)) = a * R + b

      end do

      if (flag .eq. 1) then
         write (*,1000) (10000.0d0 / x_o(indx(nx))),  &
                        (10000.0d0 / x_o(indx(1)))
 1000    format(t1, &
                '! Input line ', i3, /, &
                t1, '! Warning:  Extinction coefficients are ', &
                'required from ', /, &
                t1, '! ', f10.1, ' to ', f10.1, ' Angstroms', /, &
                t1, '! At least one wavelength lies ', &
                'outside the bounds of applicability ', &
                /, t1, '! of the reddening law odonnell_94', &
                /, t1,  '! Wavelengths must be longer than ', &
                '1000.0 Angstroms.', &
                /, t1, '! Outside the bounds, extinction ', &
                'coefficients are set to be the same as for ', &
                /, t1, '! nearest wavelength encompassed by law.', &
                /, t1, '!')
      end if

      return
      end


      subroutine law_cal (x_o, nx, xsorted, a1_o)

!     Reddening law of Calzetti (1997)

!     This function computes the monochromatic extinction A_lambda
!     relative to E(B-V) for a set of wavelengths.
!     The reddening law is defined by formulae derived by
!     Calzetti (1997) by examining UV continuum slopes of
!     star-forming galaxies (especially BCD's) as a function of
!     the differential extinction implied by the ratio of Halpha to
!     Hbeta.  Thus, the law applies to the CONTINUUM of
!     star-forming galaxies.  The extinction of line emission
!     appears to be significantly higher than the extinction of
!     the continuum.  The shape of the reddening law is fixed.

!     What is done here is to compute A_lambda
!     relative to E(B-V) directly from the formulae, i.e. without
!     modification.  Note that no
!     corrections for effective wavelength shifts, as
!     would arise from the kinds of stars used to define the law,
!     are applied here.  The output here
!     is exactly what would arise from the Calzetti's formulae
!     without adjustments.

!     INPUT:
!     x_o =  wave number, in inverse microns 
!     nx = number of elements in x
!     xsorted = Boolean defining whether or not x_o is sorted into
!               ascending order
!             = true, if x_o is already sorted
!             = false, if not

!     OUTPUT:
!     a1_o = extinction coefficient A_lambda/E(B-V)
!     flag = flag to indicate extrapolation in wavelength.
 
!     Let

!     x = 1 / lambda in microns^(-1)
!     a1_o = A_lambda / E(B-V)

!     The following specify array dimensions are needed for 
!     those variables in synext6 which aren't passed in subroutine calls

!     nfilpmax = maximum number of filters which can be picked
!     nflxmax = maximum number of input points in each object spectrum
!     nresmax = maximum number of input points in each filter 
!               response curve
!     nsubmax = maximum number of response curve samples 
!     ntaumax = maximum number of optical depths
!     nzmax = maximum number of redshifts

      include 'synext6.inc'

!     Now the declarations

      logical xsorted
      integer flag, i, indx(nflxmax), nx
      real*8 a1_o(nflxmax), &
             x, xmax, xmid, x_o(nflxmax), x_s(nflxmax)

      flag = 0

!     Sort x_o if necessary.  Place in x.

      if (xsorted) then
         do i=1,nx
            indx(i) = i
            x_s(i) = x_o(i)
         end do
      else
         call hpsort (x_o, nx, indx, x_s)
      end if

!     Work out extinction coefficients.

      xmid = 1.0d0 / 0.63d0
      xmax = 1.0d0 / 0.12d0
      do i=1,nx

!     First determine if wavelength is out of bounds.

         if (x_s(i) .gt. xmax) then
            x_s(i) = xmax
            flag = 1
         else if (x_s(i) .lt. 1.0) then
            x_s(i) = 1.0d0
            flag = 1
         end if
         x = x_s(i)

!        Now compute coefficients on basis of x

         if (x .le. xmid) then
             a1_o(i) = ((1.86d0 - 0.48d0 * x) * x - 0.1d0) * x + 1.73d0
         else
             a1_o(i) = 2.656d0 * (-2.156d0 + 1.509d0 * x  &
                       - 0.198d0 * x**2 + 0.011d0 * x**3) + 4.88d0
         end if
         
      end do

!     Flag errors

      if (flag .eq. 1) then
         write (*,1000) (10000.0d0 / x_o(indx(nx))),  &
                        (10000.0d0 / x_o(indx(1)))
 1000    format(t1, &
                '! Input line ', i3, /, &
                t1, '! Warning:  Extinction coefficients are ', &
                'required from', /, &
                t1, '! ', f10.1, ' to ', f10.1, ' Angstroms', /, &
                t1, '! At least one wavelength lies ', &
                'outside the bounds of applicability ', &
                /, t1, '! of the reddening law calzetti_97', &
                /, t1, '! Wavelengths must lie between ', &
                '1200.0 and 10000.0 Angstroms.', &
                /, t1, '! Outside the bounds, extinction ', &
                'coefficients are set to be the same as for ', &
                /, t1, '! nearest wavelength encompassed by law.', &
                /, t1, '!')
      end if

      return
      end
