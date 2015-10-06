!This file is part of TDEFit.

!TDEFit is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TDEFit is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with TDEFit.  If not, see <http://www.gnu.org/licenses/>.

function likelihood() result(l)
    use tdefit_data
    use tdefit_interface, only: kroupa, chabrier, get_sim_index, dmdt, bandmag, &
                                bbflux, src_bb_func, set_event
    use constants


#include "tdefit.fpp"

    real :: l
    real :: velouter, blr_dist, lmh, lmh2, lmh3
    real :: relativity_corr, blr_err, rblr, int_scat, var_mean, var_var, var_norm

    real, dimension(1) :: blr_fbs, blr_mdots, blr_mags, blr_routs, blr_rphots
    integer, dimension(1) :: blr_penalties

    character*2 :: this_band
    integer :: bi, ei, i, e

    if (trial_penalty .eq. huge(1.d0)) then
        l = -lhuge
        return
    endif

    l = 1.d0

    do e = 1, event_n
        call set_event(e)
        if (lf_ms .ne. 0) then
            if (trial_object_type(cur_event) .eq. 0) then
                if (lf_ms .eq. 1) then
                    l = l * kroupa(trial_ms0(cur_event))
                elseif (lf_ms .eq. 2) then
                    l = l * chabrier(trial_ms0(cur_event))
                endif
            else
                l = l * dexp(-0.5d0*((trial_ms0(cur_event) - 0.575d0)/0.05d0)**2) ! Eyeballed from DeGennaro 2007
            endif
        endif

        if (lf_beta .eq. 1) then
            l = l / trial_beta(cur_event)**2
        endif

        if (lf_band_sigma .gt. 0.d0) then
            l = l * dexp(-(trial_offset_X1(cur_event)**2 + trial_offset_GN(cur_event)**2 + trial_offset_Pg(cur_event)**2 + &
                           trial_offset_Pr(cur_event)**2 + trial_offset_Pi(cur_event)**2 + trial_offset_Pz(cur_event)**2)/lf_band_sigma**2)
        endif

        if (lf_variability .eq. 1) then
            lmh = dlog10(trial_mh0(cur_event))
            lmh2 = lmh*lmh
            lmh3 = lmh2*lmh
            var_mean = (-2.47699 - 15.4362*lmh + 4.44035*lmh2 - 0.314841*lmh3)/(1. - 0.747169*lmh - 1.06062*lmh2 + 0.00563869*lmh3)
            var_var  = (5.09238 + 6.93489*lmh - 2.25845*lmh2 + 0.166751*lmh3)/(1. + 3.02952*lmh - 0.608163*lmh2 + 0.0527849*lmh3)
            var_norm = (1. + 6.8571*lmh - 1.98696*lmh2 + 0.141459*lmh3)/(3.57951 + 3.21078*lmh - 1.11082*lmh2 + 0.08313*lmh3)
            l = l * dexp(-0.5d0*(trial_variability(cur_event) - var_mean)**2/var_var**2)*var_norm
        endif

        ! Custom for PS1-10jh
        if (lf_mh .eq. 1) then
            warnings = "Warning: Special conditions being used for PS1-10jh, please disable this code in general."
            ! From Gezari 2012
            if (trial_mh0(cur_event) .gt. 4.d6) then
                l = l * dexp(-(trial_mh0(cur_event) - 4.d6)**2/(4.d6**2))
            else
                l = l * dexp(-(trial_mh0(cur_event) - 4.d6)**2/(2.d6**2))
            endif
        endif

        if (lf_heii_dispersion .eq. 1 .or. lf_halpha_dispersion .eq. 1) then
            do i = 1, event_blrpts(cur_event)
                if (event_blr_bands(i,cur_event) .eq. 'Ha' .and. lf_halpha_dispersion .eq. 0) cycle
                if (event_blr_bands(i,cur_event) .eq. 'He' .and. lf_heii_dispersion .eq. 0) cycle

                this_band = '51'
                call dmdt(event_blr_times(i:i,cur_event) + trial_toff(cur_event), blr_fbs, &
                          .false., trial_menv(:cur_npts,cur_event))
                call dmdt(event_blr_times(i:i,cur_event) + trial_toff(cur_event), blr_mdots, &
                          .true., trial_menv(:cur_npts,cur_event))
                call bandmag(event_blr_times(i:i,cur_event) + trial_toff(cur_event), blr_fbs, &
                             blr_mdots, [this_band], blr_mags, blr_penalties, blr_routs, blr_rphots)
                if (any(blr_penalties .ne. 0)) then
                    l = -lhuge
                    return
                endif

                blr_mags = dlog10(blr_mags)

                if (event_blr_bands(i,cur_event) .eq. 'Ha') then
                    blr_dist = -23.5351 + 0.572604*blr_mags(1)! + 1.d0
                else
                    blr_dist = -24.3247 + 0.572604*blr_mags(1)! - 1.d0
                endif
                blr_dist = 10.d0**blr_dist*c*day

                if (trial_blr_model(cur_event) .eq. 0) then
                    rblr = blr_routs(1)
                elseif (trial_blr_model(cur_event) .eq. 1) then
                    rblr = blr_rphots(1)
                else
                    rblr = max(blr_rphots(1),blr_routs(1))
                endif

                ! Check if material is out far enough or out too far.
                if (event_blr_exists(i,cur_event) .eqv. .false.) then
                    !int_scat = 0.011569 ! Scatter from measurement errors only
                    !int_scat = 0.32412  ! Includes instrinsic scatter
                    !if (rblr .gt. blr_dist) then
                    !    l = -lhuge
                    !    return
                    !endif
                    int_scat = 0.0520734d0 ! Max likelihood result
                    blr_err = dabs(int_scat+1.0084*dsqrt(3.5061-0.15952*blr_mags(1)+0.0018148*blr_mags(1)**2))
                    l = l * (1.d0 - 0.5d0*erfc((dlog10(blr_dist) - dlog10(rblr))/(sqrt2*blr_err)))
                    if (print_extra) then
                        print *, '[Extra]: No exists ', event_blr_bands(i,cur_event), &
                            (event_blr_times(i,cur_event) + trial_toff(cur_event))/day, &
                            (1.d0 - 0.5d0*erfc((dlog10(blr_dist) - dlog10(rblr))/(sqrt2*blr_err)))
                    endif
                else
                    !int_scat = 0.0054674 ! Scatter from measurement errors only
                    !int_scat = 0.059088 ! Includes instrinsic scatter
                    !if (rblr .lt. blr_dist) then
                    !    l = -lhuge
                    !    return
                    !endif
                    int_scat = 0.0358824 ! Max likelihood result
                    blr_err = dabs(int_scat+1.0084*dsqrt(3.5061-0.15952*blr_mags(1)+0.0018148*blr_mags(1)**2))
                    l = l * (1.d0 - 0.5d0*erfc((dlog10(rblr) - dlog10(blr_dist))/(sqrt2*blr_err)))
                    if (print_extra) then
                        print *, '[Extra]: Exists ', event_blr_bands(i,cur_event), &
                            (event_blr_times(i,cur_event) + trial_toff(cur_event))/day, &
                            (1.d0 - 0.5d0*erfc((dlog10(rblr) - dlog10(blr_dist))/(sqrt2*blr_err)))
                    endif
                endif

                ! Check if line velocity is correct.
                if (event_blr_exists(i,cur_event) .and. event_blr_vels(i,cur_event) .ne. 0.d0) then
                    ! Velocity must correspond to between Keplerian and escape
                    !velouter = max(dsqrt(trial_gmh(cur_event)/blr_dist), (twopi*trial_gmh(cur_event)/event_blr_times(i,cur_event))**one_th)
                    velouter = dsqrt(trial_gmh(cur_event)/blr_dist)

                    !print *, 'velouter', dlog10(velouter), dlog10(trial_mh(cur_event))
                    if (event_blr_vels(i,cur_event) .lt. velouter) then
                        if (print_extra) then
                            print *, '[Extra]: Measured velocity too small compared to model.', dabs(event_blr_vels(i,cur_event) - velouter)/7.d7, &
                                dexp(-0.5d0*((event_blr_vels(i,cur_event) - velouter)/7.d7)**2)
                        endif
                        l = l * dexp(-0.5d0*((event_blr_vels(i,cur_event) - velouter)/7.d7)**2)
                        !l = -lhuge
                        !return
                    endif
                    if (event_blr_vels(i,cur_event) .gt. sqrt2*velouter) then
                        if (print_extra) then
                            print *, '[Extra]: Measured velocity too large compared to model.', dabs(event_blr_vels(i,cur_event) - sqrt2*velouter)/7.d7, &
                                dexp(-0.5d0*((event_blr_vels(i,cur_event) - sqrt2*velouter)/7.d7)**2)
                        endif
                        l = l * dexp(-0.5d0*((event_blr_vels(i,cur_event) - sqrt2*velouter)/7.d7)**2)
                        !l = -lhuge
                        !return
                    endif
                endif
            enddo
        endif

        ! Added just for Dougie
        !if (trial_beta(cur_event) .lt. model_beta_destroy(trial_model(cur_event)+1)) then
        !    l = l * (1.d0 - 0.5d0*erfc((orb_period - 4.d0*yr)/(sqrt2*yr)))
        !    warnings = 'Period constraint for Dougie imposed'
        !endif

        if (print_extra) then
            if (trial_beta(cur_event) .gt. model_beta_destroy(trial_model(cur_event)+1)) then
                print *, 'Fully destroyed, no orbital period'
            elseif (trial_ecor(cur_event) .gt. 0.d0) then
                print *, 'Orbital period:', orb_period/yr
            else
                print *, 'Partial disruption, orbital energy not currently calculated.'
            endif
        endif
    enddo

    if (l .eq. 0.d0) then
        l = -lhuge
    else
        l = dlog(l)
    endif
end function
