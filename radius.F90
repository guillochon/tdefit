! From Tout 1996, input is in solar masses, output in solar radii.

function radius(m) result(r)
    use tdefit_data, only: use_solar_radius, trial_object_type, cur_event
    use constants
    real, intent(in) :: m
    real :: r, mi

    if (use_solar_radius) then
        r = 1.d0
        return
    endif

    mi = m

    if (trial_object_type(cur_event) .eq. 0) then
        if (mi .gt. 0.1d0) then
            r = (1.71535900d0*mi**2.5d0 + 6.59778800d0*mi**6.5d0 + 10.0885500d0*mi**11.d0 + &
                 1.01249500d0*mi**19.d0 + 0.07490166d0*mi**19.5d0)/&
                (0.01077422d0 + 3.08223400d0*mi**2.d0 + 17.8477800d0*mi**8.5d0 + &
                 mi**18.5d0 + 0.00022582d0*mi**19.5d0)
        else
            if (mi .lt. 2.d-3) mi = 2.d-3
            mi = dlog10(mi) + 3.d0
            r = 1.d1**((-127.32993911050335*dsqrt(mi) + 367.70087766013705*mi - 398.36566829235505*mi**1.5 + &
                        191.90631559845060*mi**2 - 34.68499944805778*mi**2.5)/ &
                       (132.28529792732346*dsqrt(mi) - 383.55084153652444*mi + 417.23446941374976*mi**1.5 - &
                        201.82646889377915*mi**2 + 36.63019728351757*mi**2.5))
        endif
    else
        if (mi .ge. 1.433) then
            r = 1.d-3
            return
        endif
        r = 0.0112d0*dsqrt((1.433d0/mi)**two_th - (mi/1.433e0)**two_th)
    endif
end function
