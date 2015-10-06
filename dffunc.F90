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

function obs_df_func(lr) result(flux)
    use tdefit_interface, only: dffunc, obs_bb_func


    real, intent(in) :: lr
    real :: flux

    flux = dffunc(obs_bb_func, lr)
end function

function arv_df_func(lr) result(flux)
    use tdefit_interface, only: dffunc, arv_bb_func


    real, intent(in) :: lr
    real :: flux

    flux = dffunc(arv_bb_func, lr)
end function

function src_df_func(lr) result(flux)
    use tdefit_interface, only: dffunc, src_bb_func


    real, intent(in) :: lr
    real :: flux

    flux = dffunc(src_bb_func, lr)
end function

function dffunc(bbfunc, lr) result(flux)
#include "tdefit.fpp"
    use constants
    use tdefit_data
    use tdefit_interface, only: disk_temp, bbflux, bbsed, ang_frac


    real, intent(in) :: lr
    real, external :: bbfunc

    real :: flux, temp, fcor, r, frac, afrac

    r = 10.d0**lr

    flux = 0.d0

    !if (r .gt. trial_rp(cur_event)) then
    !    temp = (dfmd*2.d0*G*trial_ms(cur_event)/trial_rs(cur_event)/(pi_sigma_b*trial_rp(cur_event)*trial_rs(cur_event)))**(0.25d0)!*&
    !        !(r/trial_rp(cur_event))**(-0.75d0)
    !else
        temp = disk_temp(dfmd, r)*trial_temp_mult(cur_event)*df_temp_mult
    !endif

    if (temp .ne. temp) then
        print *, 'nan temp', r
        call exit(0)
    endif

    afrac = ang_frac(r)

    fcor = 1.d0

    if (trial_use_fcor(cur_event)) then
        if (temp .lt. trial_tlimit(cur_event)) then
            fcor = fcor*trial_tlimit(cur_event)/temp
        elseif (temp .gt. 3.d4) then
            fcor = fcor*(min(temp, 1.d5)/3.d4)**(0.82d0*trial_fcor(cur_event))
        endif
    endif

    flux = bbflux(bbfunc, dfband, temp*fcor, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
    flux = flux / fcor**4

    if (dfreprocess .and. r .le. df_rphot) then
        flux = min(1.d0, max(1.d0 - dfcovering, 0.d0))*flux
    endif

    if (make_sed) then
        zone_sed_table = bbsed(bbfunc, temp*fcor, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event)) / fcor**4
        if (dfreprocess .and. r .lt. df_rphot) then
            zone_sed_table = max(1.d0 - dfcovering, 0.d0)*zone_sed_table
        endif
    endif

    ! Extra r factor below and l10 comes because we are integrating over log r.
    flux = fourpi * flux * r * r * afrac * l10

    if (record_max_disk .and. flux .gt. dfmaxdiskflux) then
        dfmaxdiskflux = flux
        dfmaxdisktemp = temp
    endif

    if (make_sed) then
        zone_sed_table = fourpi * r * r * afrac * zone_sed_table
    endif
end function
