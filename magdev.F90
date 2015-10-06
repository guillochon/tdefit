function annealydev(x) result(dev)
    use tdefit_interface, only: magdev, likelihood

    real, dimension(:), intent(in) :: x
    real :: dev

    dev = magdev(x, .true.)
    dev = dev - likelihood()
end function

function ydev(x) result(dev)
    use tdefit_interface, only: magdev, likelihood
    use tdefit_data, only: print_likelihood

    real, dimension(:), intent(in) :: x
    real :: dev

    dev = magdev(x, .true.)
    dev = likelihood() - dev
end function

function magdev(x, max_likelihood) result(dev)
    use tdefit_data
    use tdefit_interface, only: bandmag, dmdt, radius, set_trial_vars, get_sim_index, set_event

#include "tdefit.fpp"

    real, dimension(:), intent(in) :: x
    logical, intent(in), optional :: max_likelihood
    real :: dev, dev_const, tscl, penconst, dummy, ldev
    integer :: min_md_loc, max_md_loc, i, e
    logical :: lnmult

    dfreject = .false.

#ifdef DEBUG
    if (size(x) .ne. nvars) then
        call tdefit_print('Wrong size input to magdev')
        call exit(0)
    endif
#endif

    lnmult = .true.
    if (present(max_likelihood)) then
        if (.not. max_likelihood) lnmult = .false.
    else
        lnmult = .false.
    endif

    dev = 0.0d0
    magfail = .false.

    if (discard_failed_integrals) then
        dfcnt = 0
        dffailcnt = 0
        bbcnt = 0
        bbfailcnt = 0
    endif

    if (use_soft_penalties .or. hard_penalties .eq. 0 .or. (hard_penalties .eq. 1 .and. nstep .le. nanneal)) then
        if (use_soft_penalties .or. hard_penalties .eq. 0) then
            penconst = 0.d0
        else
            penconst = dlog10(huge(1.d0))*(dble(nstep)/nanneal)**(mctemp*(1.d0 - dble(nstep)/nanneal))
        endif

        dev_const = mag_penalty*cur_npts*10.d0**penconst
    endif

    do e = 1, event_n
        call set_event(e)

        dummy = set_trial_vars(x)

        if (trial_penalty .eq. huge(1.d0)) then
#ifdef PRINT_PENALTY_REASONS
            print *, 'Penalty: Bad trial vars.'
#endif
            dev = huge(1.d0)
            dfreject = .true.
            return
        endif

        if (.not. use_soft_penalties .and. (hard_penalties .eq. 2 .or. nstep .gt. nanneal)) then
            if (trial_beta(cur_event) .lt. min_model_beta(trial_model(cur_event)+1) .or. trial_beta(cur_event) .gt. max_model_beta(trial_model(cur_event)+1)) then
#ifdef PRINT_PENALTY_REASONS
                print *, 'Penalty: Beta out of range.'
#endif
                dev = huge(1.d0)
                dfreject = .true.
                return
            endif
        endif

        trial_fbs(:cur_npts,cur_event) = 0.d0
        trial_mdots(:cur_npts,cur_event) = 0.d0
        trial_times(:cur_npts,cur_event) = event_times(:cur_npts,cur_event)/(1.d0 + trial_z(cur_event))
        call dmdt(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_fbs(:cur_npts,cur_event), &
                  .false., trial_menv(:cur_npts,cur_event))
        call dmdt(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_mdots(:cur_npts,cur_event), &
                  .true., trial_menv(:cur_npts,cur_event))
        !print *, 'magdev 2'
        call bandmag(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_fbs(:cur_npts,cur_event), &
                     trial_mdots(:cur_npts,cur_event), event_bands(:cur_npts,cur_event), &
                     trial_mags(:cur_npts,cur_event), event_penalties(:cur_npts,cur_event))

        !if (discard_failed_integrals .and. (dffailcnt .gt. 0 .or. bbfailcnt .gt. 0)) then
        if (discard_failed_integrals .and. dffailcnt .gt. 0) then
#ifdef PRINT_PENALTY_REASONS
            print *, 'Penalty: Failed df integral.'
#endif
            dev = huge(1.d0)
            dfreject = .true.
            return
        endif

        if (use_soft_penalties .or. hard_penalties .eq. 0 .or. (hard_penalties .eq. 1 .and. nstep .le. nanneal)) then
            if (penalize_early_time) then
                min_md_loc = minloc(trial_fbs(:cur_npts,cur_event), 1, mask = event_penalties(:cur_npts,cur_event) .eq. 0)
                tscl = dabs(maxval(trial_times(:cur_npts,cur_event)) - minval(trial_times(:cur_npts,cur_event)))
                if (min_md_loc .eq. 0) then
                    min_md_loc = maxloc(trial_fbs(:cur_npts,cur_event), 1, mask = event_penalties(:cur_npts,cur_event) .eq. 2)
                    if (count(event_types(:cur_npts,cur_event) .eq. 0 .and. &
                        event_penalties(:cur_npts,cur_event) .eq. 2) .eq. 1) then
                        dev = huge(1.d0)
                        dfreject = .true.
                        return
                    endif
                    if (min_md_loc .eq. 0) then
                        dev = huge(1.d0)
                        dfreject = .true.
                        return
                    endif
                endif
            endif

            ! It might need to be sqrt of errors for lnmult part!
            do i = 1, cur_npts
                if (event_types(i,cur_event) .eq. 0 .and. event_penalties(i,cur_event) .eq. 0) then
                    dev = dev + &
                        0.5d0*event_weights(i,cur_event)*(trial_mags(i,cur_event) - event_ABs(i,cur_event))**2 / &
                        (event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event))
                    if (lnmult) dev = dev + &
                        0.5d0*event_weights(i,cur_event)*dlog(event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event))
                endif

                if (penalize_early_time) then
                    if (event_types(i,cur_event) .eq. 0 .and. event_penalties(i,cur_event) .eq. 1) then
                        dev = dev + &
                            dev_const * (1.d0 + dabs((trial_times(i,cur_event) + trial_toff(cur_event))/tscl))
                    endif
                    if (event_types(i,cur_event) .eq. 0 .and. event_penalties(i,cur_event) .eq. 2) then
                        dev = dev + &
                            dev_const * (1.d0 + dabs(dlog10(trial_fbs(min_md_loc,cur_event)) - dlog10(trial_fbs(i,cur_event))))
                    endif
                endif

                if (lnmult .and. event_types(i,cur_event) .eq. 1 .and. event_penalties(i,cur_event) .eq. 0) then
                    ldev = 1.d0 + erf((trial_mags(i,cur_event) - event_ABs(i,cur_event))/&
                           (dsqrt(2.d0*(event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event)))))
                    if (ldev .eq. 0.d0) then
                        dev = huge(1.d0)
                        dfreject = .true.
                        return
                    else
                        dev = dev - event_weights(i,cur_event)*dlog(0.5d0*ldev)
                    endif
                endif
            enddo
#ifdef DEBUG
            if (dev .ne. dev .or. dev .ge. huge(1.d0) .or. dev .le. 0.d0) then
                do i = 1,cur_npts
                    print *, event_weights(i,cur_event)
                    print *, event_ABs(i,cur_event)
                    print *, event_errs(i,cur_event)
                    print *, trial_mags(i,cur_event)
                    print *, trial_fbs(i,cur_event)
                enddo
                call exit(0)
            endif
#endif
        else
            if (penalize_early_time) then
                if (any(event_penalties(:cur_npts,cur_event) .ne. 0 .and. event_types(:cur_npts,cur_event) .eq. 0)) then
                    dev = huge(1.d0)
                    dfreject = .true.
#ifdef PRINT_PENALTY_REASONS
                    print *, 'Penalty: Either no dm/dt or magnitude in model where flux is observed.'
                    print *, 'Number of violating points:', count(event_penalties(:cur_npts,cur_event) .eq. 1 .and. &
                        event_types(:cur_npts,cur_event) .eq. 0), &
                        count(event_penalties(:cur_npts,cur_event) .eq. 2 .and. event_types(:cur_npts,cur_event) .eq. 0)
#endif
                    return
                endif
            else
                max_md_loc = maxloc(trial_fbs(:cur_npts,cur_event), 1, mask = event_penalties(:cur_npts,cur_event) .eq. 0)
                if (any(event_penalties(max_md_loc:cur_npts,cur_event) .ne. 0 .and. &
                    event_types(max_md_loc:cur_npts,cur_event) .eq. 0)) then
                    dev = huge(1.d0)
                    dfreject = .true.
#ifdef PRINT_PENALTY_REASONS
                    print *, 'Penalty: Either no dm/dt or magnitude in model where flux is observed.'
                    print *, 'Number of violating points:', count(event_penalties(:cur_npts,cur_event) .eq. 1 .and. &
                        event_types(:cur_npts,cur_event) .eq. 0), &
                        count(event_penalties(:cur_npts,cur_event) .eq. 2 .and. event_types(:cur_npts,cur_event) .eq. 0)
#endif
                    return
                endif
            endif
            do i = 1, cur_npts
                if (event_types(i,cur_event) .eq. 0 .and. event_penalties(i,cur_event) .eq. 0) then
                    dev = dev + &
                        0.5d0*event_weights(i,cur_event)*(trial_mags(i,cur_event) - event_ABs(i,cur_event))**2 / &
                        (event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event))
                    if (lnmult) dev = dev + &
                        0.5d0*event_weights(i,cur_event)*dlog(event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event))
                endif
                if (lnmult .and. event_types(i,cur_event) .eq. 1 .and. event_penalties(i,cur_event) .eq. 0 .and. &
                    event_ABs(i,cur_event) .gt. trial_mags(i,cur_event)) then
                    ldev = 1.d0 + erf((trial_mags(i,cur_event) - event_ABs(i,cur_event))/&
                           (dsqrt(2.d0*(event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event)))))
                    if (ldev .eq. 0.d0) then
                        dev = huge(1.d0)
                        dfreject = .true.
                        return
                    else
                        dev = dev - event_weights(i,cur_event)*dlog(0.5d0*ldev)
                    endif
                endif
            enddo
        endif

#ifdef DEBUG
        if (dev .ne. dev .or. dev .eq. huge(1.d0) .or. dev .le. 0.d0) then
            print *, 'failed at 7a'
            call exit(0)
        endif
#endif

        if (dev_print) then
            do i = 1, cur_npts
                if (event_types(i,cur_event) .eq. 0 .or. event_ABs(i,cur_event) .gt. trial_mags(i,cur_event)) then
                    print *, dev, event_ABs(i,cur_event), trial_mags(i,cur_event), trial_fbs(i,cur_event), &
                        event_types(i,cur_event), event_penalties(i,cur_event), &
                        0.5d0*event_weights(i,cur_event)*((trial_mags(i,cur_event) - event_ABs(i,cur_event))**2 / &
                        (event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event)) + &
                        count([lnmult])*0.5d0*(dlog(event_errs(i,cur_event) + trial_variance2(cur_event) + trial_variability2(cur_event))))
                endif
            enddo
        endif
        if (any(event_types(:cur_npts,cur_event) .eq. 0 .and. penalize_early_time .and. &
            event_penalties(:cur_npts,cur_event) .ne. 0)) then
            magfail = .true.
        endif
#ifdef DEBUG
        if (dev .ne. dev .or. dev .ge. 1.d300 .or. dev .le. 0.d0) then
            print *, trial_penalty, cur_npts
            print *, 'failed at 7'
            call exit(0)
        endif
#endif
    enddo

    dev = trial_penalty * dev
end function
