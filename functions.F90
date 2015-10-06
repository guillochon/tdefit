function sigma(mh) result(sig)
    use constants
    real, intent(in) :: mh
    real :: sig

    sig = 2.43184d5*(mh*imsun)**0.2358490566037736d0 !Gultekin 2009
end function

function rs(mstar) result(rad) !Tout 1996
    use constants
    real, intent(in) :: mstar
    real :: rad

    rad = rsun*(1.715359d0*mstar**2.5d0 + 6.597788d0*mstar**6.5d0 + 10.08855d0*mstar**11 + 1.012495d0*mstar**19 + 0.07490166d0*mstar**19.5d0)/&
          (0.01077422d0 + 3.082234d0*mstar**2 + 17.84778d0*mstar**8.5d0 + mstar**18.5d0 + 0.00022582d0*mstar**19.5d0)
end function

function kroupa(mstar) result(prob)
    use constants
    use tdefit_data
    real, intent(in) :: mstar
    real :: prob

    if (mstar .lt. kc1) then
        prob = mstar**(-0.3d0)
    elseif (mstar .ge. kc1 .and. mstar .lt. kc2) then
        prob = kc1*mstar**(-1.3d0)
    else
        prob = kc1*kc2*mstar**(-2.3d0)
    endif
end function

function chabrier(mstar) result(prob)
    real, intent(in) :: mstar
    real :: prob

    if (mstar .gt. 1.d0) then
        prob = 4.43d-2*mstar**(-1.3d0)
    else
        prob = 0.158*dexp(-0.5d0*((dlog10(mstar) - dlog10(0.079))/0.69d0)**2)
    endif
end function

!function wang(a, m) result(prob) !From section 5.2 of Wang & Merritt 2004.
!    use tdefit_data
!    use constants
!#include "tdefit.fpp"
!    real, intent(in) :: a
!    integer, intent(in) :: m
!
!    real :: mh, lambda, rt, r0, ne, tre, prob
!    
!    mh = model_pars(m,MOD_BHMASS)
!    rt = rsun*(mh/msun)**one_th/min_beta
!    if (a .le. rt) then
!        prob = 0.d0
!    else
!        lambda = mh/msun
!        r0 = rt/a
!        ne = mcamax(m)**(gamr - 3.d0)*mh*(3.d0 - gamr)*a**(2.d0 - gamr)/4.d0/pi
!        tre = sqrt2*(G*mh/a)**three_halfs/pi/G**2.d0/msun/ne/dlog(lambda)
!        prob = ne/tre/dlog(1.d0/r0)
!    endif
!end function
!
!function tal(a, mstar, rstar, m) result(prob) !From Alexander & Hopman 2003.
!    use tdefit_data
!    use constants
!#include "tdefit.fpp"
!    real, intent(in) :: a, mstar, rstar
!    integer, intent(in) :: m
!
!    real :: mh, lambda, rt, prob, r
!    
!    mh = model_pars(m,MOD_BHMASS)
!    lambda = mh/msun
!    rt = rstar*(mh/mstar)**one_th/min_beta
!    r = 2.d0*a - rt
!
!    prob = (dlog(lambda)/dlog((2.d0*a - rt)/rt) - 1.d0/4.d0)*a**(7.d0/2.d0 - 2.d0*gamr)
!end function
