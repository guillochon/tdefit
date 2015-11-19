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

subroutine load_event(e)
#include "tdefit.fpp"
    use tdefit_data
    use constants
    use tdefit_interface, only: bbflux, obs_bb_func
    use tdefit_util, only: sort2


    integer, intent(in)                     :: e
    integer                                 :: fn, i, j, bi, ei, band_count, version
    real                                    :: flux, dummy, first_time
    logical                                 :: next_band
    character*2                             :: cur_band
    character*2, dimension(:), allocatable  :: temp_bands
    character*3                             :: time_unit
    real, dimension(:), allocatable         :: temp_times

    fn = 11

    if (my_pe .eq. 0) print *, "Loading event '" // trim(event_fnames(e)) // "'..."
    open(unit = fn, file = trim(event_path) // trim(event_fnames(e)) // ".dat", status='old', action='read')
    read(fn, *) version, dummy, dummy, dummy, time_unit
    read(fn, *) event_claimed_z(e)
    read(fn, *) event_nh(e)
    read(fn, *) event_nhcorr(e), event_restframe(e)
    do i = 1, event_npts(e)
        read(fn, *) event_bands(i,e), event_times(i,e), event_ABs(i,e), event_errs(i,e), event_types(i,e)
    enddo
    do i = 1, event_blrpts(e)
        read(fn, *) event_blr_times(i,e), event_blr_vels(i,e), event_blr_bands(i,e), event_blr_exists(i,e)
    enddo
    close(fn)

    first_time = minval(event_times(:,e))
    event_times(:,e) = event_times(:,e) - first_time
    event_blr_times(:,e) = event_blr_times(:,e) - first_time

    if (time_unit == "MJD") then
        event_times(:,e) = event_times(:,e)*day
        event_blr_times(:,e) = event_blr_times(:,e)*day
    else
        event_times(:,e) = event_times(:,e)*yr
        event_blr_times(:,e) = event_blr_times(:,e)*yr
    endif
    if (event_restframe(e) .eq. 1) then
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

    if (band_count .ne. event_nbest_bands(e) - nextra_bands) then
        write (*, *), 'Error, actual band count does not match header.'
        call exit(0)
    endif

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
        if (remove_extinction_corr .and. event_nhcorr(e) .eq. 1) then
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

    ! Loop through band list
    cur => ll(e)%p
    i = 0
    do while (associated(cur))
        i = i + 1
        event_best_bands(i,e) = cur%band
        cur => cur%next
    enddo

    do i = 1, nextra_bands
        event_best_bands(event_nbest_bands(e) - nextra_bands + i,e) = extra_bands(i)
    enddo

    ! The error bars don't include intrinsic variability, this updates the error bars to more appropriately reflect the uncertainty.
    !where (event_types .eq. 0)
    !    event_errs = event_errs + mag_fac*dlog10(1.d0 + variability)
    !endwhere
    event_errs(:event_npts(e),e) = event_errs(:,e)**2
end subroutine
