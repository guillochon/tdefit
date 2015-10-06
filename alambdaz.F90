function alambdaz(nu, z, nh, nhsrc) result(avfrac)
    use tdefit_interface, ONLY: alambda
    use tdefit_data, ONLY: trial_source_rv, local_rv, bbband, no_source_extinction, cur_event


    real, intent(in) :: nu, z, nh, nhsrc
    real :: avfrac, avsrcfrac, dummy

    if (bbband .eq. 'Lb') then
        avfrac = 1.d0
    elseif (no_source_extinction .or. nhsrc .eq. 0.d0) then
        call alambda(nu, nh, local_rv, avfrac)
        avfrac = 10.d0**(-0.4d0*avfrac)
    else
        call alambda(nu, nh, local_rv, avfrac)
        call alambda(nu*(1.d0 + z), nhsrc, trial_source_rv(cur_event), avsrcfrac)
        avfrac = 10.d0**(-0.4d0*(avfrac+avsrcfrac))
    endif
end function
