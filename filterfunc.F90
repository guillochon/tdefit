function filterfunc(nu) result(frac)
    use tdefit_data
    use tdefit_interface, ONLY: bisect

    integer :: i
    real, INTENT(IN) :: nu
    real :: frac

    i = bisect(filt_resp(bandi,1:filt_len(bandi),1), nu)
!#ifdef DEBUG
    if (i .eq. 0 .or. i .eq. filt_len(bandi) + 1) then
        print *, nu, filt_resp(bandi,1,1)
        print *, "filterfunc: nu out of filter range"
        call exit(0)
    endif
!#endif
    frac = filt_resp(bandi,i,2) + dfilt_resp(bandi,i)*(nu - filt_resp(bandi,i,1))
end function

function filterintfunc(nu) result(frac)
    use tdefit_data
    use constants

    integer :: i
    real, INTENT(IN) :: nu
    real :: frac!, frac2

    if (bandi .eq. -1) then
        frac = 1.d0
        return
    endif

    i = min(max(ceiling(filt_const(bandi)*(nu - filtint(bandi,1,1))),1),filtintlen-1)
    !i = ceiling(filt_const(bandi)*(nu - filtint(bandi,1,1)))
#ifdef DEBUG
    if (i .eq. 0 .or. i .eq. filtintlen + 1) then
        print *, nu, filtint(bandi,1,1), filtint(bandi,filtintlen,1)
        print *, "filterintfunc: nu out of filter range"
        call exit(0)
    endif
#endif
    frac = filtint(bandi,i,2) + filtint(bandi,i,3)*(nu - filtint(bandi,i,1))

    !frac = (nu - filtint(bandi,i,1))/(filtint(bandi,i+1,1) - filtint(bandi,i,1))
    !frac2 = 0.5d0*(1.d0 - dcos(frac*pi))
    !frac = filtint(bandi,i,2)*(1.d0 - frac2) + filtint(bandi,i+1,2)*frac2
end function
