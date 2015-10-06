function obs_bb_func(lnu) result(flux)
    ! This function returns the observed spectrum, including the filter and extinction.
#include "tdefit.fpp"
    !use, intrinsic :: ieee_arithmetic
    use constants
    use tdefit_data
    use tdefit_interface, ONLY: filterintfunc, alambdaz, is_abnormal


    real, intent(in) :: lnu
    real :: nu, flux, nusrc

    ! Mancone and Gonzalez 2006
    ! Blanton & Brinkmann 2003
    ! This function can also be used for calculating AB magnitudes, as the divided nu cancels the extra nu introduced
    ! by integrating over log space.

    bbcalls = bbcalls + 1
    nu = 1.d1**lnu
    nusrc = nu*bb1pz

    flux = x_const*nusrc/bbtemp
    if (flux .gt. lhuge) then
        flux = 0.d0
    else
        flux = nusrc**3*filterintfunc(nu)*alambdaz(nu,bbz,bbnh,bbnhsrc)*flux_const/&
            (dexp(flux) - 1.d0)
    endif

    if (bbmultbynu) then
        flux = flux*nu
    endif

    flux = flux*l10

#ifdef DEBUG
    !if (.not. ieee_is_normal(flux)) then
    if (is_abnormal(flux)) then
        print *, 'flux abnormal in obs_bb_func', nu, bbtemp
        call exit(0)
    endif
#endif
end function

function arv_bb_func(lnu) result(flux)
#include "tdefit.fpp"
    ! This function returns the observed spectrum, including the filter and extinction.
    !use, intrinsic :: ieee_arithmetic
    use constants
    use tdefit_data
    use tdefit_interface, ONLY: alambdaz, is_abnormal


    real, intent(in) :: lnu
    real :: nu, flux, nusrc

    ! Mancone and Gonzalez 2006
    ! Blanton & Brinkmann 2003
    ! This function can also be used for calculating AB magnitudes, as the divided nu cancels the extra nu introduced
    ! by integrating over log space.

    bbcalls = bbcalls + 1
    nu = 1.d1**lnu
    nusrc = nu*bb1pz
    flux = x_const*nusrc/bbtemp
    if (flux .gt. lhuge) then
        flux = 0.d0
    else
        flux = nusrc**3*alambdaz(nu,bbz,bbnh,bbnhsrc)*flux_const/&
            (dexp(flux) - 1.d0)
    endif

    if (bbmultbynu) then
        flux = flux*nu
    endif

    flux = flux*l10

#ifdef DEBUG
    !if (.not. ieee_is_normal(flux)) then
    if (is_abnormal(flux)) then
        print *, 'flux abnormal in arv_bb_func', nu, bbtemp
        call exit(0)
    endif
#endif
end function

function src_bb_func(lnu) result(flux)
#include "tdefit.fpp"
    ! This function returns the emitted spectrum, with no filter or extinction.
    !use, intrinsic :: ieee_arithmetic
    use constants
    use tdefit_data


    real, intent(in) :: lnu
    real :: nu, flux

    bbcalls = bbcalls + 1
    nu = 10.d0**lnu
    flux = x_const*nu/bbtemp
    if (flux .gt. lhuge) then
        flux = 0.d0
    else
        flux = flux_const * nu**3 / (dexp(flux) - 1.d0) ! Times l10, moved to bbflux for efficiency
    endif

    if (bbmultbynu) then
        flux = flux*nu
    endif

    flux = flux*l10

#ifdef DEBUG
    !if (.not. ieee_is_normal(flux)) then
    if (is_abnormal(flux)) then
        print *, 'flux abnormal in arv_bb_func', nu, bbtemp
        call exit(0)
    endif
#endif
end function

function filtnormfunc(lnu) result(flux)
#include "tdefit.fpp"
    !use, intrinsic :: ieee_arithmetic
    use constants
    use tdefit_interface, ONLY: filterintfunc, is_abnormal

    real, intent(in) :: lnu
    real :: flux

    flux = l10*filterintfunc(1.d1**lnu)
#ifdef DEBUG
    !if (.not. ieee_is_normal(flux)) then
    if (is_abnormal(flux)) then
        print *, 'flux abnormal in filtnormfunc', lnu
        call exit(0)
    endif
#endif
end function
