function cosmo_ez(z) result(ez)
    use constants
    real, intent(in) :: z
    real :: ez

    ez = dsqrt(omega_m*(1.d0 + z)**3 + omega_k*(1.d0 + z)**2 + omega_l)
end function

function cosmo_dc(z) result(dc)
    use tdefit_interface, ONLY: cosmo_ez
    use constants
    real, parameter :: tol = 1.d-6
    real, intent(in) :: z
    real :: dc, err
    integer :: neval, ierr

    call qag(cosmo_ez,0.d0,z,tol,tol,3,dc,err,neval,ierr)
    dc = dc*cosmo_dH
end function

function cosmo_dm(z) result(dm)
    use tdefit_interface, ONLY: cosmo_dc
    use constants
    real, intent(in) :: z
    real :: dm

    if (omega_k .gt. 0.d0) then
        dm = cosmo_dH/dsqrt(omega_k)*dsinh(dsqrt(omega_k)*cosmo_dc(z)/cosmo_dH)
    elseif (omega_k .lt. 0.d0) then
        dm = cosmo_dH/dsqrt(dabs(omega_k))*dsin(dsqrt(dabs(omega_k))*cosmo_dc(z)/cosmo_dH)
    else
        dm = cosmo_dc(z)
    endif
end function

function cosmo_da(z) result(da)
    use tdefit_interface, ONLY: cosmo_dm
    real, intent(in) :: z
    real :: da

    da = cosmo_dm(z) / (1.d0 + z)
end function

function cosmo_dl(z) result(dl)
    use tdefit_interface, ONLY: cosmo_dm
    real, intent(in) :: z
    real :: dl

    dl = (1.d0 + z) * cosmo_dm(z)
end function
