function ang_frac(r) result(frac)
    use constants
    use tdefit_data

    real, intent(in) :: r
    real :: frac, t_visc

    frac = 1.d0
    if (.not. trial_full_disk_coverage(cur_event)) then
        if (r .lt. 2.d0*trial_rp(cur_event) + &
            (8.d0*((dftime - first_accretion_time)/pi)**2*trial_gmh(cur_event))**one_th) then
            frac = min(dsqrt(trial_rp(cur_event)/r), 1.d0)
        else
            frac = 0.d0
        endif
        if (.not. fixed_angle_frac .and. r .gt. 2.d0*trial_rp(cur_event)) then
            t_visc = twopi*dsqrt(r**3.d0/trial_gmh(cur_event))/trial_alphhr(cur_event)
            frac = max(0.d0,min(frac, (dftime - (first_accretion_time + &
                   pi*dsqrt(0.125d0*(r - 2.d0*trial_rp(cur_event))**3/&
                   trial_gmh(cur_event))))/t_visc))
        endif
    endif
end function
