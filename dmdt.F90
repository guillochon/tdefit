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

subroutine get_sim_index(betafrac, bi, ei)
    use tdefit_data
    use tdefit_interface, ONLY: bisect


#include "tdefit.fpp"
    real, intent(out) :: betafrac
    integer, intent(out) :: bi, ei
    integer :: beta_bi, beta_ei

    beta_bi = model_index(trial_model(cur_event)+1)
    if (trial_model(cur_event) .eq. nmodels - 1) then
        beta_ei = nruns
    else
        beta_ei = model_index(trial_model(cur_event)+2)-1
    endif

    bi = beta_bi + bisect(edat(beta_bi:beta_ei,E_BETA), trial_beta(cur_event)) - 1
    ei = bi + 1

    betafrac = (trial_beta(cur_event) - edat(bi,E_BETA)) / (edat(ei,E_BETA) - edat(bi,E_BETA))
end subroutine

subroutine dmdt(tdes, dm, add_delay, im, rhom, mode, ades)
    use constants
    use tdefit_data, only: first_accretion_time, trial_ms0, trial_rs0, trial_ecor, &
                    trial_mh, trial_gmh, trial_r_ibco, d_emin, d_emid, d_emax, &
                    trial_rp, trial_beta, edat, ddat_ncols, trial_alphhr, orb_period, &
                    model_beta_destroy, trial_model, viscous_dmdt, ddat, &
                    cur_event, trial_fout, trial_viscous_time, maxdmdttime
    use tdefit_interface, only: interp_flash_output, get_sim_index


#include "tdefit.fpp"

    real, dimension(:), intent(in) :: tdes
    real, dimension(size(tdes)), intent(out) :: dm
    logical, intent(in) :: add_delay
    real, dimension(size(tdes)), intent(out), optional :: im, rhom
    real, dimension(:), allocatable :: temporary
    integer, intent(in), optional :: mode
    integer, parameter :: intl = 100
    real, intent(in), optional :: ades

    real :: betafrac, emin, emid, emax, epscor, edenom, orb_ener_correct, &
            time_corr, relativity_corr, sc_mh, circular_time, tnot, peaktime, &
            tmidint, dtint
    real, dimension(size(tdes)) :: ratio, md1, md2, es, newt, imd1, imd2, rhomd1, rhomd2
    real, dimension(size(tdes),intl) :: md1int, md2int, mdint, tint, es1int, es2int
    !real, dimension(:), allocatable :: kernel, temp_mdot
    integer :: i, j, z, bi, ei, begi, endi, cbegi

    dm = 0.d0

    if (present(im)) then
        im = 0.d0
    endif

    if (present(rhom)) then
        rhom = 0.d0
    endif

    ! Adjusts trial_beta(cur_event) to account for stars becoming less centrally concentrated. Should be replaced in the future by 5/3-5/3 parameter study.

    !trial_beta(cur_event) = min(trial_beta(cur_event), max_sim_beta)

    !if (trial_beta(cur_event) .lt. max(min_sim_beta, beta_noloss)) then
    !    print *, 'trial_beta(cur_event)', trial_beta(cur_event), min_sim_beta, beta_noloss
    !endif
    !if (trial_beta(cur_event) .lt. max(min_sim_beta, beta_noloss)) return

    call get_sim_index(betafrac, bi, ei)

    !Some truncation errors can arise, this checks for them
    epscor = 0.d0
    if (present(ades)) then
        epscor = epscor - 0.5d0*trial_gmh(cur_event)/ades
    endif

    ! From Kesden 2012
    relativity_corr = 1.d0/dsqrt(1.d0 - 0.5d0*trial_r_ibco(cur_event)/trial_rp(cur_event))

    ! For debugging
    !print *, 'warning betafrac set to 0.d0 for debugging'
    !betafrac = 0.0d0
    
    ! emin should be set to the same small number for all runs so light curve can extend to infinite time.
    emin = (betafrac*(d_emin(ei) - d_emin(bi)) + d_emin(bi))
    emid = (betafrac*(d_emid(ei) - d_emid(bi)) + d_emid(bi))
    emax = (betafrac*(d_emax(ei) - d_emax(bi)) + d_emax(bi))
    if (emin .gt. 0.d0) then
        print *, 'No bound material! Probably an error...'
        print *, d_emin(ei), d_emin(bi), d_emax(ei), d_emax(bi), bi
        print *, emin, emax, betafrac, trial_beta(cur_event)
        call exit(0)
    endif

    sc_mh = edat(bi,E_MPERT)

    ! Use unscaled time for accessing dm/de table.

    ! The way trial_ecor(cur_event) is put in here is probably not correct!!!!!!
    !newt = pi_G*sc_mh*isqrt2/((sqrt2*tdes/(pi_G*sc_mh))**(-two_th) - &
    !    G*trial_mh(cur_event)/trial_rp(cur_event)**2*rsun*trial_ecor(cur_event))**three_halfs
    !newt = tdes

    !This still seems to be glitchy, not sure why...
    orb_ener_correct = -trial_ecor(cur_event)*G*trial_beta(cur_event)**2*sc_mh/(rsun*(sc_mh*imsun)**two_th)

    if (trial_ecor(cur_event) .gt. 0.d0 .and. trial_beta(cur_event) .lt. model_beta_destroy(trial_model(cur_event)+1)) then
        orb_period = pi_G*isqrt2/dabs(orb_ener_correct)**three_halfs*sc_mh
    else
        orb_period = huge(1.d0)
    endif

    time_corr = dsqrt(trial_mh(cur_event)/sc_mh)/trial_ms0(cur_event)*trial_rs0(cur_event)**three_halfs/relativity_corr**three_halfs

    first_accretion_time = pi_G*isqrt2*sc_mh/(-emin - orb_ener_correct)**three_halfs
        
    ! Based on Ulmer 1999
    !first_accretion_time = max(first_accretion_time,&
    !    pi_G*isqrt2*sc_mh*rsun**three_halfs/&
    !    (G*(sc_mh*msun**2)**one_th*trial_beta(cur_event)**2 - rsun*orb_ener_correct)**three_halfs)

    if (add_delay .and. viscous_dmdt) then
        peaktime = ((betafrac*(maxdmdttime(ei) - maxdmdttime(bi)) + maxdmdttime(bi)))*&
            dsqrt(trial_mh(cur_event)/edat(bi,E_MPERT))/trial_ms0(cur_event)*trial_rs0(cur_event)**three_halfs*&
            (1.d0/dsqrt(1.d0 - 0.5d0*trial_r_ibco(cur_event)/trial_rp(cur_event)))**three_halfs

        if (trial_viscous_time(cur_event) .ne. 0.d0) then
            circular_time = trial_viscous_time(cur_event)
        else
            circular_time = first_accretion_time*(1.d0/trial_alphhr(cur_event) - 1.d0)
                        !*trial_fout(cur_event)**three_halfs
        endif

        time_corr = time_corr * (1.d0 + circular_time/peaktime)
    else
        circular_time = 0.d0
    endif

    first_accretion_time = first_accretion_time*time_corr
    newt = tdes/time_corr

    begi = 1
    endi = size(tdes)

    edenom = -emin

    do i = begi, endi
        if (newt(i) .le. 0.d0) then
            begi = i+1
            cycle
        endif

        ratio(i) = (-((pi_G*isqrt2*sc_mh/newt(i))**two_th) - emin - orb_ener_correct)/edenom

        if (ratio(i) .le. 0.d0) then
            begi = i+1
            cycle
        endif

        if (ratio(i) .ge. 1.d0 .or. newt(i) .gt. orb_period) then
            endi = i-1
            exit
        endif
    enddo

    if (begi .gt. endi) return

    if (viscous_dmdt) then
        do i = begi, endi
            do j = 1, intl
                tint(i,j) = (pi_G*isqrt2*sc_mh)/(-((j-1.)/(intl-1.)*ratio(i)*edenom - emin - orb_ener_correct))**1.5d0
                es1int(i,j) = d_emin(bi) + (j-1.)/(intl-1.)*ratio(i)*(-d_emin(bi))
                es2int(i,j) = d_emin(ei) + (j-1.)/(intl-1.)*ratio(i)*(-d_emin(ei))
            enddo
            call interp_flash_output(DDAT_ARR, bi, begi, es1int(i,:), md1int(i,:))
            call interp_flash_output(DDAT_ARR, ei, begi, es2int(i,:), md2int(i,:))
        enddo
    endif

    es(begi:endi) = d_emin(bi) + ratio(begi:endi)*(-d_emin(bi))
    call interp_flash_output(DDAT_ARR, bi, begi, es(begi:endi), md1(begi:endi))
    if (present(im)) then
        call interp_flash_output(IDDAT_ARR, bi, begi, es(begi:endi), imd1(begi:endi))
    endif
    if (present(rhom)) then
        call interp_flash_output(RHODDAT_ARR, bi, begi, es(begi:endi), rhomd1(begi:endi))
    endif

    es(begi:endi) = d_emin(ei) + ratio(begi:endi)*(-d_emin(ei))
    call interp_flash_output(DDAT_ARR, ei, begi, es(begi:endi), md2(begi:endi))
    if (present(im)) then
        call interp_flash_output(IDDAT_ARR, ei, begi, es(begi:endi), imd2(begi:endi))
    endif
    if (present(rhom)) then
        call interp_flash_output(RHODDAT_ARR, ei, begi, es(begi:endi), rhomd2(begi:endi))
    endif

    !local_kerw = ceiling(final_kerw*(endi-begi+1))
    !allocate(kernel(6*int(local_kerw)+1))
    !do i = (size(kernel)-1)/2 + 1, size(kernel)
    !    kernel(i) = dexp(-((i - ((size(kernel)-1)/2 + 1))**2.d0/2.d0/local_kerw**2.d0))
    !enddo
    !kernel(1:(size(kernel)-1)/2) = kernel(size(kernel):(size(kernel)-1)/2+2:-1)

    if (viscous_dmdt) then
        mdint(begi:endi,:) = one_th*(twopi*G*sc_mh)**two_th/&
            tint(begi:endi,:)**five_th*dexp(betafrac*&
            (dlog(md2int(begi:endi,:))-dlog(md1int(begi:endi,:))) + dlog(md1int(begi:endi,:)))*&
            trial_ms0(cur_event)/time_corr

        do i = begi, endi
            do j = 2, intl
                tmidint = 0.5d0*(tint(i,j)+tint(i,j-1))
                dtint = tint(i,j)-tint(i,j-1)
                dm(i) = dm(i) + dexp((tmid-tdes(i))/circular_time)*mdint(i,j)*dtint
                dm(i) = dm(i)/circular_time
            enddo
        enddo
    else
        dm(begi:endi) = one_th*(twopi*G*sc_mh)**two_th/&
            newt(begi:endi)**five_th*dexp(betafrac*&
            (dlog(md2(begi:endi))-dlog(md1(begi:endi))) + dlog(md1(begi:endi)))*&
            trial_ms0(cur_event)*newt(begi:endi)/tdes(begi:endi)
    endif

    ! This bit of code would modify the fallback curve to include viscous
    ! effects. This changes the slope in dm/de to e^-2/3 if e < e_visc.
    !do i = begi, endi

    !enddo

    if (present(im)) then
        ! Normalized to one, not actual integrated density.
        ! calculation in future.
        im(begi:endi) = dexp(betafrac*(dlog(imd2(begi:endi))-dlog(imd1(begi:endi))) + &
                        dlog(imd1(begi:endi)))
    endif
    if (present(rhom)) then
        ! Normalized to one, not actual density.
        ! calculation in future.
        rhom(begi:endi) = dexp(betafrac*(dlog(rhomd2(begi:endi))-dlog(rhomd1(begi:endi))) + &
                          dlog(rhomd1(begi:endi)))
    endif

    if (present(mode)) then
        if (mode .eq. 2) then
            allocate(temporary(size(dm)))
            temporary = 0.d0
            do j = 2, size(dm)
                temporary(j) = temporary(j-1) + 0.5d0*(newt(j) - newt(j-1))*(dm(j) + dm(j-1))
            enddo
            dm = temporary
            deallocate(temporary)
        endif
    endif

    ! Convert back to the actual event scale
    if (trial_ecor(cur_event) .gt. 0.d0 .and. trial_beta(cur_event) .lt. model_beta_destroy(trial_model(cur_event)+1)) then
        orb_period = orb_period/((trial_mh(cur_event)/sc_mh)**one_th/trial_rs0(cur_event))**three_halfs
    endif

    ! Smooth the merged curve a bit
    !allocate(temp_mdot(begi:endi))
    !temp_mdot = 0.d0
    !do i = begi, endi
    !    kersum = 0.d0
    !    do j = 1, size(kernel)
    !        if (i+j-(size(kernel)-1)/2 .lt. begi) cycle
    !        if (i+j-(size(kernel)-1)/2 .gt. endi) cycle
    !        temp_mdot(i) = temp_mdot(i) + dm(i+j-(size(kernel)-1)/2)*kernel(j)
    !        kersum = kersum + kernel(j)
    !    enddo
    !    temp_mdot(i) = temp_mdot(i) / kersum
    !enddo
    !dm(begi:endi) = temp_mdot(begi:endi)
    !deallocate(temp_mdot)
end subroutine

function im_root(t)
    use tdefit_data, only: trial_fout, dfintval, cur_event
    use tdefit_interface, only: dmdt

    real, intent(in) :: t
    real :: im_root
    real, dimension(1) :: dm, im

    call dmdt([10.d0**t], dm, .false., im)

    im_root = trial_fout(cur_event)*(dfintval - im(1)) - 1.d0
end function
