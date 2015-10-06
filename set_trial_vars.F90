function set_trial_vars(x) result(penalty)
    use tdefit_data
    use tdefit_interface, only: radius, set_var


#include "tdefit.fpp"

    real, dimension(:), intent(in) :: x
    real, dimension(size(x)) :: xlimited
    real :: penalty, min_allowed_beta, max_allowed_beta, req_aspin, trial_val, penconst, val
    logical :: z_search

    integer :: i, j

    penalty = 1.d0

    xlimited = x

    if (minval(x) .lt. 0.d0 .or. maxval(x) .gt. 1.d0) then
        if (.not. use_soft_penalties .and. (hard_penalties .eq. 2 .or. nstep .gt. nanneal)) then
#ifdef PRINT_PENALTY_REASONS
            print *, 'Penalty: Selected point out of domain.'
            print *, minval(x), maxval(x)
#endif
            penalty = huge(1.d0)
            trial_penalty = penalty
            return
        else
            do i = 1, size(x)
                if (x(i) .lt. 0.d0) then
                    penalty = penalty*dexp(-x(i)*mag_penalty)
                    xlimited(i) = 0.d0
                elseif (x(i) .gt. 1.d0) then
                    penalty = penalty*dexp((x(i) - 1.d0)*mag_penalty)
                    xlimited(i) = 1.d0
                endif
            enddo
        endif
    endif

    j = 0
    do i = 1, nvars
        if (var_events(i) .ne. 0 .and. var_events(i) .ne. cur_event) cycle
        if (var_locks(i)) then
            val = locked_trial_vars(i)
        else
            val = xlimited(i)
        endif
        if (var_types(i) .eq. 0) then
            trial_val = min_search(i) + val*(max_search(i) - min_search(i))
        elseif (var_types(i) .eq. 1) then
            trial_val = min_search(i)*10.d0**(val*(dlog10(max_search(i)) - dlog10(min_search(i))))
        elseif (var_types(i) .eq. 2) then
            trial_val = val*(max_search(i) - min_search(i)) + min_search(i)
        elseif (var_types(i) .eq. 3) then
            trial_val = val
        else
            print *, '[set_trial_vars] Invalid variable type.'
            call exit(0)
        endif 
        call set_var(var_names(i), trial_val)
    enddo

    ! Special circumstances apply to beta, whose range depends on trial_model.
    min_allowed_beta = min_model_beta(trial_model(cur_event)+1)
    max_allowed_beta = max_model_beta(trial_model(cur_event)+1)

    z_search = .false.
    max_aspin = 0.998d0 ! Maximum spin possible, only change if max_aspin is set.
    do i = 1, nvars
        if (var_events(i) .ne. 0 .and. var_events(i) .ne. cur_event) cycle
        if (var_locks(i)) then
            val = locked_trial_vars(i)
        else
            val = xlimited(i)
        endif
        if (var_names(i) .eq. "aspin") then
            max_aspin = max_search(i)
        elseif (var_names(i) .eq. "beta") then
            trial_val = min_allowed_beta + val*(max_allowed_beta - min_allowed_beta) 
            call set_var(var_names(i), trial_val)
        endif
        if (var_names(i) .eq. "z") then
            z_search = .true.
        endif
    enddo

    if (.not. z_search) trial_z(cur_event) = event_claimed_z(cur_event)

    call set_derived_trial_vars

    if (.not. use_soft_penalties .and. (hard_penalties .eq. 2 .or. nstep .gt. nanneal)) then
        if (trial_rp(cur_event) .lt. trial_r_ibco(cur_event)) then
#ifdef PRINT_PENALTY_REASONS
            print *, 'Penalty: rp < r_ibco'
            print *, 'rp:', trial_rp(cur_event), 'r_ibco:', trial_r_ibco(cur_event)
#endif
            penalty = huge(1.d0)    
            trial_penalty = penalty
            return
        endif
!        if (trial_beta(cur_event) .lt. min_allowed_beta .or. trial_beta(cur_event) .gt. max_allowed_beta) then
!#ifdef PRINT_PENALTY_REASONS
!            print *, 'Penalty: Beta out of range'
!            print *, 'Beta:', trial_beta(cur_event), min_allowed_beta, max_allowed_beta
!
!#endif
!            penalty = huge(1.d0)    
!            trial_penalty = penalty
!            return
!        endif
        if (disallow_unphysical_models) then
            if (trial_object_type(cur_event) .eq. 0) then
                if ((trial_model(cur_event) .eq. 0 .and. trial_ms0(cur_event) .lt. 0.1d0) .or. &
                    (trial_model(cur_event) .eq. 1 .and. trial_ms0(cur_event) .gt. 1.0d0 .and. trial_ms0(cur_event) .lt. 20.d0)) then
                    penalty = huge(1.d0)
                    trial_penalty = penalty
                    return
                endif
            else
                if ((trial_ms0(cur_event) .gt. 1.433d0 .or. trial_ms0(cur_event) .lt. 1.d-2) .or. &
                   (((trial_model(cur_event) .eq. 0 .and. trial_ms0(cur_event) .lt. 1.d0) .or. &
                    (trial_model(cur_event) .eq. 1 .and. trial_ms0(cur_event) .ge. 1.d0)))) then
                    penalty = huge(1.d0)
                    trial_penalty = penalty
                    return
                endif
            endif
        endif
    else
        if (use_soft_penalties .or. hard_penalties .eq. 0) then
            penconst = mag_penalty
        else
            penconst = mag_penalty*dlog10(huge(1.d0))*(dble(nstep)/nanneal)**(mctemp*(1.d0 - dble(nstep)/nanneal))
        endif

        if (disallow_unphysical_models) then
            if (trial_object_type(cur_event) .eq. 0) then
                if (trial_model(cur_event) .eq. 0 .and. trial_ms0(cur_event) .lt. 0.1d0) then
                    penalty = penalty*10.d0**((0.1d0/trial_ms0(cur_event) - 1.d0)*penconst)
                endif
                if (trial_model(cur_event) .eq. 1 .and. trial_ms0(cur_event) .gt. 1.0d0 .and. trial_model(cur_event) .lt. 20.d0) then
                    penalty = penalty*10.d0**((min(trial_ms0(cur_event)/1.0d0, 20.d0/trial_ms0(cur_event)) - 1.d0)*penconst)
                endif
            else
                if (trial_ms0(cur_event) .lt. 1.d-2) then
                    penalty = penalty*10.d0**((1.d-2/trial_ms0(cur_event) - 1.d0)*penconst)
                endif
                if (trial_ms0(cur_event) .gt. 1.433d0) then
                    penalty = penalty*10.d0**((trial_ms0(cur_event)/1.433d0 - 1.d0)*penconst)
                endif
                if (trial_model(cur_event) .eq. 0 .and. trial_ms0(cur_event) .lt. 1.d0) then
                    penalty = penalty*10.d0**((1.d0/trial_ms0(cur_event) - 1.d0)*penconst)
                endif
                if (trial_model(cur_event) .eq. 1 .and. trial_ms0(cur_event) .ge. 1.d0) then
                    penalty = penalty*10.d0**((trial_ms0(cur_event)/1.d0 - 1.d0)*penconst)
                endif
            endif
        endif

        if (trial_beta(cur_event) .gt. max_allowed_beta) then
            penalty = penalty*10.d0**((trial_beta(cur_event)/max_allowed_beta - 1.d0)*penconst)
            trial_beta(cur_event) = max_allowed_beta
            trial_rp(cur_event) = (trial_mh(cur_event)/trial_ms(cur_event))**one_th*trial_rs(cur_event)/trial_beta(cur_event)
        endif

        !if (trial_rp(cur_event) .lt. trial_rg(cur_event)/2.d0) then
        !    penalty = 10.d0**((trial_rg(cur_event)/2.d0/trial_rp(cur_event) - 1.d0)*penconst)
        !
        !    trial_beta(cur_event) = trial_rp(cur_event)/trial_rg(cur_event)*2.d0
        !    trial_aspin(cur_event) = max_aspin
        !    trial_rp(cur_event) = (trial_mh(cur_event)/trial_ms(cur_event))**one_th*trial_rs(cur_event)/trial_beta(cur_event)
        !endif

        if (trial_beta(cur_event) .lt. min_allowed_beta) then
            penalty = penalty*10.d0**((min_allowed_beta/trial_beta(cur_event) - 1.d0)*penconst)
            trial_beta(cur_event) = min_allowed_beta
            trial_rp(cur_event) = (trial_mh(cur_event)/trial_ms(cur_event))**one_th*trial_rs(cur_event)/trial_beta(cur_event)
        endif

        if (trial_rp(cur_event) .lt. trial_r_ibco(cur_event)) then
            req_aspin = 2.d0*dsqrt(trial_rp(cur_event)*(trial_r_ibco(cur_event) - trial_rp(cur_event)))/(trial_r_ibco(cur_event))
            if (trial_aspin(cur_event) .lt. req_aspin) then
                penalty = penalty*10.d0**(((1.d0 - trial_aspin(cur_event)) / (1.d0 - req_aspin))*penconst)
                trial_aspin(cur_event) = req_aspin
            endif
            if (max_aspin .lt. trial_aspin(cur_event)) then
                penalty = penalty*10.d0**(((1.d0 - max_aspin) / (1.d0 - trial_aspin(cur_event)))*penconst)
                trial_aspin(cur_event) = max_aspin
            endif
            trial_rp(cur_event) = trial_r_ibco(cur_event)*min_r_fac
        endif
    endif

    trial_penalty = penalty
end function

subroutine set_derived_trial_vars
    use tdefit_data
    use tdefit_interface, only: radius, cosmo_dl

    real :: z1, z2

    real :: spin_acos

    trial_ms0(cur_event) = trial_ms(cur_event)*imsun
    trial_mh0(cur_event) = trial_mh(cur_event)*imsun
    trial_rs(cur_event) = rsun*radius(trial_ms0(cur_event))*trial_rsc(cur_event)
    trial_rs0(cur_event) = trial_rs(cur_event)*irsun
    trial_gmh(cur_event) = G*trial_mh(cur_event)
    trial_rg(cur_event) = 2.d0*trial_gmh(cur_event)/c2
    trial_rp(cur_event) = (trial_mh(cur_event)/trial_ms(cur_event))**one_th*trial_rs(cur_event)/trial_beta(cur_event)

    ! Bardeen & Teukolsky 1972
    z1 = ((1.d0-trial_aspin(cur_event)**2)**one_th)
    z1 = z1*(((1.d0+trial_aspin(cur_event))**one_th)+((1.d0-trial_aspin(cur_event))**one_th))
    z1 = 1.d0+z1
    z2 = dsqrt(3.d0*trial_aspin(cur_event)*trial_aspin(cur_event)+z1*z1)
    trial_bh_rms(cur_event) = 3.d0+z2-dsqrt((3.d0-z1)*(3.d0+z1+2.d0*z2))

    ! Prograde
    trial_bh_rbs(cur_event) = 2.d0-trial_aspin(cur_event)+2.d0*dsqrt(1.d0-trial_aspin(cur_event))
    ! Retrograde
    !trial_bh_rbs(cur_event) = 2.d0+trial_aspin(cur_event)+2.d0*dsqrt(1.d0+trial_aspin(cur_event))

    trial_eps_edd(cur_event) = 1.d0-dsqrt(1.d0-2.d0/(3.d0*trial_bh_rms(cur_event)))

    trial_r_isco(cur_event) = trial_bh_rms(cur_event)*trial_rg(cur_event)/2.d0*min_r_fac
    trial_r_ibco(cur_event) = trial_bh_rbs(cur_event)*trial_rg(cur_event)/2.d0*min_r_fac

    if (disallow_unphysical_models) then
        if (trial_object_type(cur_event) .eq. 0) then
            trial_mu_e(cur_event) = 1.18d0
        elseif (trial_object_type(cur_event) .eq. 1) then
            trial_mu_e(cur_event) = 2.d0
        endif
    endif

    trial_ledd(cur_event) = edd_const*trial_mh(cur_event)*trial_mu_e(cur_event)

    trial_dl(cur_event) = cosmo_dl(trial_z(cur_event))
    trial_1pz(cur_event) = trial_z(cur_event) + 1.d0

    trial_variance2(cur_event) = trial_variance(cur_event)**2
    trial_variability2(cur_event) = trial_variability(cur_event)**2

    ! These constants are from Done 2011, moved to here to avoid having to
    ! calculate them with each call to disk_temp.
    trial_yms(cur_event)=dsqrt(trial_bh_rms(cur_event))
    spin_acos=dacos(trial_aspin(cur_event))
    trial_y1(cur_event)=2.0*dcos(one_th*(spin_acos-pi))
    trial_y2(cur_event)=2.0*dcos(one_th*(spin_acos+pi))
    trial_y3(cur_event)=-2.0*dcos(one_th*spin_acos)
end subroutine
