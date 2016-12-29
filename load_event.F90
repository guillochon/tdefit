!This file is part of TDEFit.

!TDEFit is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TDEFit is distributed in the hope that it will be useful,
!but WITH(out) ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with TDEFit.  If not, see <http://www.gnu.org/licenses/>.

subroutine load_events
    use tdefit_data

    integer :: fn, e, i, event_max_nbest_bands

    open(unit = fn, file = "event_list.dat", status='old', action='read')
    read(fn, *) event_n
    allocate(event_blrpts(event_n))
    allocate(event_claimed_z(event_n))
    allocate(event_fnames(event_n))
    allocate(event_min_aspin(event_n))
    allocate(event_nbest_bands(event_n))
    allocate(event_nh(event_n))
    allocate(event_nhcorr(event_n))
    allocate(event_npts(event_n))
    allocate(event_restframe(event_n))
    allocate(ll(event_n))
    do e = 1, event_n
        read(fn, *) event_fnames(e)
    enddo
    close(fn)

    event_npts = 0
    event_blrpts = 0
    do e = 1, event_n
        call load_event(e,.true.)
    enddo

    event_max_npts = maxval(event_npts)
    event_max_blrpts = maxval(event_blrpts)

    allocate(event_bands(event_max_npts,event_n))
    allocate(event_time_units(event_max_npts,event_n))
    allocate(event_times(event_max_npts,event_n))
    allocate(event_ABs(event_max_npts,event_n))
    allocate(event_errs(event_max_npts,event_n))
    allocate(event_devs(event_max_npts,event_n))
    allocate(event_weights(event_max_npts,event_n))
    allocate(event_penalties(event_max_npts,event_n))
    allocate(event_types(event_max_npts,event_n))
    allocate(event_blr_time_units(event_max_blrpts,event_n))
    allocate(event_blr_times(event_max_blrpts,event_n))
    allocate(event_blr_vels(event_max_blrpts,event_n))
    allocate(event_blr_bands(event_max_blrpts,event_n))
    allocate(event_blr_exists(event_max_blrpts,event_n))

    allocate(trial_1pz(event_n))
    allocate(trial_alphhr(event_n))
    allocate(trial_aspin(event_n))
    allocate(trial_beta(event_n))
    allocate(trial_bh_rbs(event_n))
    allocate(trial_bh_rms(event_n))
    allocate(trial_blr_model(event_n))
    allocate(trial_log_above_edd(event_n))
    allocate(trial_cap_at_edd(event_n))
    allocate(trial_dl(event_n))
    allocate(trial_ecor(event_n))
    allocate(trial_eps_edd(event_n))
    allocate(trial_exp_1(event_n))
    allocate(trial_exp_2(event_n))
    allocate(trial_exp_3(event_n))
    allocate(trial_exp_4(event_n))
    allocate(trial_fcor(event_n))
    allocate(trial_fout(event_n))
    allocate(trial_full_disk_coverage(event_n))
    allocate(trial_gmh(event_n))
    allocate(trial_ledd(event_n))
    allocate(trial_magoff(event_n))
    allocate(trial_mdot_floor(event_n))
    allocate(trial_mh(event_n))
    allocate(trial_mh0(event_n))
    allocate(trial_model(event_n))
    allocate(trial_ms(event_n))
    allocate(trial_ms0(event_n))
    allocate(trial_mu_e(event_n))
    allocate(trial_nhsrc(event_n))
    allocate(trial_object_type(event_n))
    allocate(trial_offset_GN(event_n))
    allocate(trial_offset_Pg(event_n))
    allocate(trial_offset_Pi(event_n))
    allocate(trial_offset_Pr(event_n))
    allocate(trial_offset_Pz(event_n))
    allocate(trial_offset_RO(event_n))
    allocate(trial_offset_U1(event_n))
    allocate(trial_offset_U2(event_n))
    allocate(trial_offset_Ub(event_n))
    allocate(trial_offset_Um(event_n))
    allocate(trial_offset_Uu(event_n))
    allocate(trial_offset_Uv(event_n))
    allocate(trial_offset_X1(event_n))
    allocate(trial_offset_X2(event_n))
    allocate(trial_offset_bV(event_n))
    allocate(trial_offset_bI(event_n))
    allocate(trial_opacity(event_n))
    allocate(trial_outflow_frac(event_n))
    allocate(trial_outflow_model(event_n))
    allocate(trial_phi(event_n))
    allocate(trial_r_ibco(event_n))
    allocate(trial_r_isco(event_n))
    allocate(trial_reprocess_temp(event_n))
    allocate(trial_rg(event_n))
    allocate(trial_rin(event_n))
    allocate(trial_rout(event_n))
    allocate(trial_rp(event_n))
    allocate(trial_rphot(event_n))
    allocate(trial_rs(event_n))
    allocate(trial_rs0(event_n))
    allocate(trial_rsc(event_n))
    allocate(trial_simple_bb(event_n))
    allocate(trial_source_rv(event_n))
    allocate(trial_temp_mult(event_n))
    allocate(trial_temperature_model(event_n))
    allocate(trial_time_dep_rin(event_n))
    allocate(trial_time_dep_rout(event_n))
    allocate(trial_tlimit(event_n))
    allocate(trial_toff(event_n))
    allocate(trial_use_fcor(event_n))
    allocate(trial_variability(event_n))
    allocate(trial_variability2(event_n))
    allocate(trial_variance(event_n))
    allocate(trial_variance2(event_n))
    allocate(trial_viscous_time(event_n))
    allocate(trial_y1(event_n))
    allocate(trial_y2(event_n))
    allocate(trial_y3(event_n))
    allocate(trial_yms(event_n))
    allocate(trial_z(event_n))

    allocate(trial_times(event_max_npts,event_n))
    allocate(trial_fbs(event_max_npts,event_n))
    allocate(trial_mdots(event_max_npts,event_n))
    allocate(trial_menv(event_max_npts,event_n))
    allocate(trial_rins(event_max_npts,event_n))
    allocate(trial_routs(event_max_npts,event_n))
    allocate(trial_rphots(event_max_npts,event_n))
    allocate(trial_mags(event_max_npts,event_n))

    nvars = 0
    do e = 1, event_n
        call load_event(e,.false.)
        call set_event(e)
        call load_defaults(1)
        call load_user_vars(1)
        cur => ll(e)%p
        event_nbest_bands(e) = nextra_bands
        do while (associated(cur))
            event_nbest_bands(e) = event_nbest_bands(e) + 1
            cur => cur%next
        enddo
    enddo

    event_max_nbest_bands = maxval(event_nbest_bands)
    allocate(event_best_bands(event_max_nbest_bands,event_n))

    do e = 1, event_n
        cur => ll(e)%p
        i = 0
        do while (associated(cur))
            i = i + 1
            event_best_bands(i,e) = cur%band
            cur => cur%next
        enddo
    enddo

    do e = 1, event_n
        do i = 1, nextra_bands
            event_best_bands(event_nbest_bands(e) - nextra_bands + i,e) = extra_bands(i)
        enddo
    enddo

    if (sum(event_npts) .le. nvars - 1) then
        print *, 'Warning: Number of measurement points must be +2 larger than ' // &
                 'number of degrees of freedom for proper chi-square measurement.'
    endif
end subroutine

subroutine load_event(e, prepare)
#include "tdefit.fpp"
    use tdefit_data
    use constants
    use tdefit_interface, only: bbflux, obs_bb_func, get_band_type, tdefit_print
    use tdefit_util, only: sort2


    integer, intent(in)                     :: e
    logical, intent(in)                     :: prepare

    integer                                 :: phoi, blri, nbuf, stat, begk, endk
    integer                                 :: loc, locspace, loctab, loccomma, loc2
    integer                                 :: locbeginquote, locendquote
    integer                                 :: fn, i, j, bi, ei, band_count, version
    real                                    :: flux, first_time
    logical                                 :: next_band, begin_quote
    character*1                             :: band_type
    character*2                             :: cur_band
    character*50                            :: var_name, str_value
    character*500                           :: buffer
    character*2, dimension(:), allocatable  :: temp_bands
    real, dimension(:), allocatable         :: temp_times

    fn = 11

    if (prepare) then
        if (my_pe .eq. 0) print *, "Preparing to load event '" // trim(event_fnames(e)) // "'..."
    else
        if (my_pe .eq. 0) print *, "Loading event '" // trim(event_fnames(e)) // "'..."
    endif
    open(unit = fn, file = trim(event_path) // trim(event_fnames(e)) // ".dat", status='old', action='read')

    ! Read through file twice, first is to simply count number of photometric and BLR points
    stat = 0
    phoi = 0
    blri = 0
    linedo: do
        buffer = ""
        read(fn, '(A500)', iostat=stat, size=nbuf, advance='no') buffer
        if (is_iostat_eor(stat)) then
            j = 0
            loc = 1
            locspace = 1
            loctab = 1
            loccomma = 1
            locbeginquote = 1
            locendquote = 1

            do while(loc .le. nbuf .and. &
                     (locspace .ne. 0 .or. loctab .ne. 0 .or. &
                      loccomma .ne. 0 .or. locbeginquote .ne. 0))

                locspace = index(trim(buffer(loc:)), " ")
                loctab = index(trim(buffer(loc:)), char(9))
                loccomma = index(trim(buffer(loc:)), ",")
                if (locspace .eq. 1 .or. loctab .eq. 1 .or. loccomma .eq. 1) then
                    loc = loc + 1
                    cycle
                endif
                locbeginquote = index(trim(buffer(loc:)), '"')
                j = j + 1
                begin_quote = .false.
                if (locbeginquote .eq. 1) then
                    begin_quote = .true.
                    loc = loc + 1
                    begk = loc
                    if (loc .le. nbuf) then
                        locendquote = index(trim(buffer(loc:)), '"')
                    else
                        locendquote = 0
                    endif
                    if (locendquote .eq. 0) then
                        call tdefit_print("Error, unclosed quote found in event file")
                        call exit(0)
                    elseif (locendquote .eq. 1) then
                        call tdefit_print("Empty value detected in event file")
                        call exit(0)
                    endif
                    endk = loc + locendquote - 2
                    loc = loc + locendquote + 1
                endif
                if (locspace .eq. 0 .and. loctab .eq. 0 .and. loccomma .eq. 0 .and. locbeginquote .eq. 0) then
                    begk = loc
                    endk = nbuf
                elseif (.not. begin_quote) then
                    loc2 = minval( (/ locspace, loctab, loccomma /), &
                        mask = (/ locspace, loctab, loccomma /) .ne. 0)
                    begk = loc
                    endk = loc + loc2 - 2
                    loc = loc + loc2
                endif

                if (j .eq. 1) then
                    var_name = ''
                    read(buffer(begk:endk),*) var_name
                    select case (trim(var_name))
                        case ('photometry', 'broad_line', 'redshift', 'nh', &
                              'nhcorr', 'restframe')
                        case default
                            cycle linedo
                    end select
                else
                    select case (trim(var_name))
                        case ('photometry')
                            if (prepare) then
                                if (j .eq. 2) event_npts(e) = event_npts(e) + 1
                                cycle linedo
                            endif
                            select case (j)
                                case (2)
                                    phoi = phoi + 1
                                    read(buffer(begk:endk),*) event_time_units(phoi,e)
                                case (3)
                                    read(buffer(begk:endk),*) event_times(phoi,e)
                                case (4)
                                    read(buffer(begk:endk),*) event_bands(phoi,e)
                                case (5)
                                    read(buffer(begk:endk),*) event_ABs(phoi,e)
                                case (6)
                                    read(buffer(begk:endk),*) event_errs(phoi,e)
                                case (7)
                                    read(buffer(begk:endk),*) event_types(phoi,e)
                                case default
                                    if (my_pe .eq. 0) then
                                        print *, 'Too many columns for photometric point! Aborting.'
                                        call exit(0)
                                    endif
                            end select
                        case ('broad_line')
                            if (prepare) then
                                if (j .eq. 2) event_blrpts(e) = event_blrpts(e) + 1
                                cycle linedo
                            endif
                            select case (j)
                                case (2)
                                    blri = blri + 1
                                    read(buffer(begk:endk),*) event_blr_time_units(blri,e)
                                case (3)
                                    read(buffer(begk:endk),*) event_blr_times(blri,e)
                                case (4)
                                    read(buffer(begk:endk),*) event_blr_vels(blri,e)
                                case (5)
                                    read(buffer(begk:endk),*) event_blr_bands(blri,e)
                                case (6)
                                    read(buffer(begk:endk),*) event_blr_exists(blri,e)
                                case default
                                    if (my_pe .eq. 0) then
                                        print *, 'Too many columns for broad line point! Aborting.'
                                        call exit(0)
                                    endif
                            end select
                        case ('redshift')
                            read(buffer(begk:endk),*) event_claimed_z(e)
                        case ('nh')
                            read(buffer(begk:endk),*) event_nh(e)
                        case ('nhcorr')
                            read(buffer(begk:endk),*) event_nhcorr(e)
                        case ('restframe')
                            read(buffer(begk:endk),*) event_restframe(e)
                        case default
                            cycle linedo
                    end select
                endif
            enddo
            if (.not. prepare .and. my_pe .eq. 0) then
                write(*,'(A15,X)',advance='no') var_name
                select case (trim(var_name))
                    case ('photometry')
                        write(*,'(A3,X,E10.3,X,A2,X,E10.3,X,E10.3,X,I1)') &
                            event_time_units(phoi,e), event_times(phoi,e), &
                            event_bands(phoi,e), event_ABs(phoi,e), event_errs(phoi,e), &
                            event_types(phoi,e)
                    case ('redshift')
                        write(*,'(E10.3)') event_claimed_z(e)
                    case ('nh')
                        write(*,'(E10.3)') event_nh(e)
                    case ('nhcorr')
                        write(*,'(L1)') event_nhcorr(e)
                    case ('restframe')
                        write(*,'(L1)') event_restframe(e)
                end select
            endif
        elseif (is_iostat_end(stat)) then
            exit
        else
            cycle
        endif
    enddo linedo

    close(fn)

    if (prepare) return

    ! Sanity check some inputs
    if (event_claimed_z(e) .le. 0.d0) then
        if (my_pe .eq. 0) print *, 'Invalid redshift specified, aborting.'
        call exit(0)
    endif
    if (event_nh(e) .lt. 0.d0) then
        if (my_pe .eq. 0) print *, 'Invalid nh specified, aborting.'
        call exit(0)
    endif
    if (event_npts(e) .le. 0) then
        if (my_pe .eq. 0) print *, 'Must specify at least one photometric point.'
        call exit(0)
    endif

    do i = 1, event_npts(e)
        band_type = get_band_type(event_bands(i,e))
        if (band_type == 'X') then
            event_ABs(i,e) = -event_ABs(i,e)*mag_fac
        endif
    enddo

    do i = 1, event_npts(e)
        if (event_time_units(i,e) == 'MJD' .or. trim(event_time_units(i,e)) == 'JD') then
            event_times(i,e) = event_times(i,e)*day
            if (trim(event_time_units(i,e)) == 'JD') then
                event_times(i,e) = event_times(i,e) - 2400000.5*day
            endif
        elseif (event_time_units(i,e) == 'yrs') then
            event_times(i,e) = event_times(i,e)*yr
        else
            print *, 'Invalid time unit specified for photometric point in event file, aborting.'
            call exit(0)
        endif
    enddo

    do i = 1, event_blrpts(e)
        if (event_blr_time_units(i,e) == 'MJD' .or. trim(event_blr_time_units(i,e)) == 'JD') then
            event_blr_times(i,e) = event_blr_times(i,e)*day
            if (trim(event_blr_time_units(i,e)) == 'JD') then
                event_blr_times(i,e) = event_blr_times(i,e) - 2400000.5*day
            endif
        elseif (event_blr_time_units(i,e) == 'yrs') then
            event_blr_times(i,e) = event_blr_times(i,e)*yr
        else
            print *, 'Invalid time unit specified for broad line point in event file, aborting.'
            call exit(0)
        endif
    enddo

    first_time = minval(event_times(:,e))
    event_times(:,e) = event_times(:,e) - first_time
    event_blr_times(:,e) = event_blr_times(:,e) - first_time

    if (event_restframe(e)) then
        event_times(:,e) = event_times(:,e)*(1.d0 + event_claimed_z(e)) ! Redshift stretches events in time, need to remove to get actual time-evolution.
        event_ABs(:,e) = event_ABs(:,e) + mag_fac*dlog(1.d0 + event_claimed_z(e))
    endif
    event_penalties(:,e) = 0
    event_weights(:,e) = 1.d0

    allocate(temp_bands(event_npts(e)))

    temp_bands = event_bands(:event_npts(e),e)
    call sort2(temp_bands, event_times(:event_npts(e),e), 2)
    temp_bands = event_bands(:event_npts(e),e)
    call sort2(temp_bands, event_ABs(:event_npts(e),e), 2)
    temp_bands = event_bands(:event_npts(e),e)
    call sort2(temp_bands, event_errs(:event_npts(e),e), 2)
    temp_bands = event_bands(:event_npts(e),e)
    call sort2(temp_bands, event_types(:event_npts(e),e), 2)
    event_bands(:event_npts(e),e) = temp_bands

    deallocate(temp_bands)

    bi = 1
    ei = 1
    band_count = 0
    cur_band = event_bands(1,e)

    allocate(ll(e)%p)
    cur => ll(e)%p
    cur%band = cur_band

    do i = 2, event_npts(e)+1
        next_band = .false.
        if (i .le. event_npts(e)) then
            if (event_bands(i,e) .ne. cur_band) then
                ei = i - 1
                next_band = .true.
            endif
        else
            ei = i - 1
            next_band = .true.
        endif
        if (next_band) then
            band_count = band_count + 1
            allocate(temp_times(ei - bi + 1))
            temp_times = event_times(bi:ei,e)
            call sort2(temp_times, event_errs(bi:ei,e))
            temp_times = event_times(bi:ei,e)
            call sort2(temp_times, event_ABs(bi:ei,e))
            temp_times = event_times(bi:ei,e)
            call sort2(temp_times, event_types(bi:ei,e))
            temp_times = event_times(bi:ei,e)
            call sort2(temp_times, event_weights(bi:ei,e))
            event_times(bi:ei,e) = temp_times
            deallocate(temp_times)
            do j = bi, ei
                if (event_errs(j,e) .eq. 0.0d0) event_errs(j,e) = upp_lim_err
                if (time_weighted) then
                    if (j .eq. bi) then
                        event_weights(j,e) = (event_times(j+1,e) - event_times(j,e))
                    elseif (j .eq. ei) then
                        event_weights(j,e) = (event_times(j,e) - event_times(j-1,e))
                    else
                        event_weights(j,e) = 0.5d0 * (event_times(j+1,e) - event_times(j-1,e))
                    endif
                endif
            enddo
            if (i .le. event_npts(e)) then
                cur_band = event_bands(i,e)
                allocate(cur%next)
                cur => cur%next
                cur%band = cur_band
                bi = i
                ei = i
            endif
        endif
    enddo

    if (time_weighted) event_weights(:,e) = event_weights(:,e) / sum(event_weights(:,e))

    allocate(temp_times(event_npts(e)))
    temp_times = event_times(:event_npts(e),e)
    call sort2(temp_times, event_bands(:event_npts(e),e), 2)
    temp_times = event_times(:event_npts(e),e)
    call sort2(temp_times, event_ABs(:event_npts(e),e))
    temp_times = event_times(:event_npts(e),e)
    call sort2(temp_times, event_errs(:event_npts(e),e))
    temp_times = event_times(:event_npts(e),e)
    call sort2(temp_times, event_types(:event_npts(e),e))
    temp_times = event_times(:event_npts(e),e)
    call sort2(temp_times, event_weights(:event_npts(e),e))
    event_times(:event_npts(e),e) = temp_times
    deallocate(temp_times)

    cur => ll(e)%p
    do while (associated(cur))
        if (remove_extinction_corr .and. event_nhcorr(e)) then
            ! The light curves are already corrected for extinction, this "uncorrects" the magnitudes by presuming the flux in the band is in the Rayleigh Jeans limit.
            flux = bbflux(obs_bb_func, cur%band, 1.d6, 0.d0, event_nh(e), 0.d0)
            flux = bbflux(obs_bb_func, cur%band, 1.d6, 0.d0, 0.d0, 0.d0)/flux
            do i = 1, event_npts(e)
                if (event_bands(i,e) .ne. cur%band) cycle
                event_ABs(i,e) = event_ABs(i,e) + mag_fac*dlog10(flux)
            enddo
        endif
        cur => cur%next
    enddo

    event_errs(:event_npts(e),e) = event_errs(:,e)**2
end subroutine
