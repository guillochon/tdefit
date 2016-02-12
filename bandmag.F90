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

subroutine bandmag(times, fbs, mdots, bands, mags, penalties, routs, rphots)
    use tdefit_data
    use constants
    use tdefit_interface, ONLY: ftoABmag, obs_df_func, src_df_func, disk_temp_root, bbflux, trapezoid, &
                                bbsed, obs_bb_func, src_bb_func, annulus_intercept, disk_temp, &
                                arv_bb_func, arv_df_func, get_sim_index, bisect, dmdt, &
                                ftomag, integrate_df, get_band_type
    use tdefit_util, ONLY: sort

#include "tdefit.fpp"

    real, dimension(:), intent(in) :: times
    real, dimension(:), intent(in) :: fbs, mdots
    character*2, dimension(:), intent(in) :: bands
    real, dimension(size(fbs)), intent(out) :: mags
    integer, dimension(size(fbs)), intent(out) :: penalties
    real, dimension(size(fbs)), intent(out), optional :: routs, rphots

    integer :: i, j, k, neval, ierr, n_zones, bi, ei
    real :: penalty, err, tworp, mag, error, zmin, zmax, reprocessed_lum, &
            total_lum, medd, circ_lum, outflowmd
    real :: dummy, coeff, flux_interior, teff, teffdisk, reflect_interior, cur_r
    real :: lmaxmdot, rconst, peaktime, betafrac, t1, t2, ratio, incl_corr, rcirc
    real, dimension(1) :: input, dm, im1, im2
    real, dimension(5) :: zone_locs
    character*1, dimension(size(bands)) :: band_types
    logical :: check
    
    mags = 0.d0
    penalty = 1.d0
    penalties = 0
    make_sed = .false.

    teff = 0.d0

    if (present(routs)) routs = 0.d0
    if (present(rphots)) rphots = 0.d0

    tworp = 2.d0*trial_rp(cur_event)

    if (produce_seds) then
        if (allocated(sed_table)) deallocate(sed_table)
        allocate(sed_table(2, nbest_times, sed_nsteps))
        sed_table = 0.d0
        sed_index = 0
    endif

    lmaxmdot = maxval(fbs)

    ! Relativity term from Kesden 2012
    call get_sim_index(betafrac, bi, ei)

    peaktime = ((betafrac*(maxdmdttime(ei) - maxdmdttime(bi)) + maxdmdttime(bi)))*&
        dsqrt(trial_mh(cur_event)/edat(bi,E_MPERT))/trial_ms0(cur_event)*trial_rs0(cur_event)**three_halfs*&
        (1.d0/dsqrt(1.d0 - 0.5d0*trial_r_ibco(cur_event)/trial_rp(cur_event)))**three_halfs
        
    if (viscous_dmdt) then
        peaktime = peaktime/trial_alphhr(cur_event)
    endif
            
    rconst = tworp + (8.d0*((peaktime - first_accretion_time)/pi)**2*trial_gmh(cur_event))**one_th
    rconst = min(trial_rout(cur_event)*tworp,rconst)

    ! This assumes that disk thickness scales with alphhr
    incl_corr = min(trial_alphhr(cur_event), 1.d0) + (1.d0 - min(trial_alphhr(cur_event), 1.d0))*dcos(trial_phi(cur_event))

    dffb = -huge(1.d0)
    dfmd = -huge(1.d0)
    do j = 1, size(fbs)
        dftime = times(j)
        dffb = max(0.d0, fbs(j))
        dfmd = max(0.d0, mdots(j) + trial_mdot_floor(cur_event))

        !!!! NOTE: NEED TO PASS ENVELOPE TO THIS function, INSTEAD OF RELYING ON GLOBALS
        dfenv = trial_menv(j,cur_event)

        if (dfmd .eq. lmaxmdot) then
            record_max_disk = .true.
            dfmaxdisktemp = 0.d0
            dfmaxdiskflux = 0.d0
        else
            record_max_disk = .false.
        endif

        if (trial_time_dep_rin(cur_event)) then
            dfri = max(tworp**3 - 4.5d0*trial_gmh(cur_event)*(max(dftime - first_accretion_time, 0.d0)*trial_alphhr(cur_event))**2, trial_r_isco(cur_event)**3)**one_th
        else
            dfri = trial_r_isco(cur_event)
        endif

        if (trial_time_dep_rout(cur_event)) then
            dfro = tworp + (8.d0*(max(dftime - first_accretion_time, 0.d0)/pi)**2*trial_gmh(cur_event))**one_th
            dfro = min(trial_rout(cur_event)*tworp,dfro)
        else
            dfro = trial_rout(cur_event)*tworp
        endif

        if (include_disk) then
            zone_locs(1:2) = [dfri, dfro]
            n_zones = 1

            !if (separate_outer_zone) then
            !    disk_t = trial_tlimit(cur_event)
            !    input = disk_temp_root([dfro])
            !    if (input(1) .lt. 0.d0) then
            !        input = (dfro + dfri)/2.d0
            !        call newt(disk_temp_root,input,1.e-4,check,error)
            !        if (.not. check .and. input(1) .gt. dfri .and. input(1) .lt. dfro) then
            !            zone_locs(n_zones + 2) = input(1)
            !            n_zones = n_zones + 1
            !        endif
            !    endif
            !endif

            call sort(zone_locs(1:n_zones+1))
        endif

        ! Calculate the total luminosity of the disk
        total_lum = 0.d0
        circ_lum = 0.d0
        df_temp_mult = 1.d0
        if (include_disk) then
            make_sed = .false.
            dfreprocess = .false.
            dfband = 'Lb'

            do i = 1, n_zones
                ierr = 0
                dfcnt = dfcnt + 1
                zmin = dlog10(zone_locs(i))
                zmax = dlog10(zone_locs(i+1))

                call integrate_df(src_df_func, zmin, zmax, mag, df_int_divs**2)

                total_lum = total_lum + mag
            enddo
            dfreprocess = .true.
            ! Compare total luminosity to Eddington luminosity, and adjust
            ! temperature if total is exceeded.

            if (total_lum .gt. trial_ledd(cur_event)) then
                df_temp_mult = (trial_ledd(cur_event)/total_lum)**0.25d0
            endif
        endif

        ! Using most bound debris as reference size
        rcirc = (trial_gmh(cur_event)*(first_accretion_time/twopi)**2)**one_th
        ! Using Tsvi's model as an upper limit
        !rcirc = 6.2e14*trial_ms0(cur_event)**(one_th-0.2d0)*(trial_mh0(cur_event)/1.d1**6.5d0)**two_th
        rcirc = max(min(rcirc,trial_rp(cur_event)/trial_alphhr(cur_event)), trial_rp(cur_event))
        if (include_circ) then
            circ_lum = G*trial_mh(cur_event)*dffb/rcirc*trial_fout(cur_event)
        endif

        dfcovering = 0.d0

        if (reprocess_model .ne. RM_CLOUDS .and. include_disk) then
            dftemp = 0.d0
        else
            medd = trial_ledd(cur_event)/(trial_eps_edd(cur_event)*c2)
            if (wind_phot) then
                outflowmd = dfmd
            else
                outflowmd = dffb
            endif

            if (trial_outflow_model(cur_event) .eq. 0) then
                ratio = max(0.d0,outflowmd-medd)/medd
            elseif (trial_outflow_model(cur_event) .eq. 1) then
                ratio = min(outflowmd,medd)/medd
            else
                ratio = outflowmd/medd
                !ratio = dffb/(trial_ledd(cur_event)/(c2*trial_eps_edd(cur_event)))
                !if (print_extra) then
                !    print *, 'ratio', total_lum, dffb*trial_eps_edd(cur_event)*c2, trial_ledd(cur_event)
                !endif
            endif

            !teffdisk = (trial_eps_edd(cur_event)*dfmd*c2*(trial_r_isco(cur_event)/dfri)**3/(fourpi_sigma_b*4.d0*dfri**2))**0.25d0

            ! Cap at bound material
            !if (teffdisk .eq. 0.d0) then
            !    df_rphot = 0.d0
            !else
            !    ! Expression calculated from ionizing-fraction.nb
            !    !df_rphot = ratio*10.d0**(-186.7218922*dexp(1.376393384*dlog10(teffdisk) - &
            !    !               0.4728966226*dlog10(teffdisk)**2))
            !    df_rphot = ratio
            !endif
            !if (df_rphot .eq. 0.d0) then
            !    df_rphot = dfri
            !else
                if (wind_phot) then
                    df_rphot = max(trial_r_isco(cur_event), &
                               trial_rphot(cur_event)*rconst*(ratio**trial_exp_2(cur_event)), rcirc)
                elseif (circ_phot) then
                    df_rphot = min(max(trial_r_isco(cur_event), rcirc), dfro)
                else
                    df_rphot = min(max(trial_r_isco(cur_event), &
                               trial_rphot(cur_event)*rconst*(ratio**trial_exp_2(cur_event)), rcirc), dfro)
                endif
            !endif
        endif

        if (present(routs)) routs(j) = dfro
        if (present(rphots)) rphots(j) = df_rphot

        if (dffb .eq. 0.d0 .or. dftime .le. first_accretion_time) then
            penalties(j) = 1
            cycle
        endif

        if (dfro .lt. dfri) then
            penalties(j) = 1
            cycle
        endif

        if (reprocess_model .ne. RM_NONE .and. .not. circ_phot) then
            call dmdt([dftime], dm, .false., im2)

            dfcovering = trial_opacity(cur_event)!*max(0.d0,im2(1))

            ! If the light escape time is comparable to the time, exclude this
            ! solution, as it would delay the light curve.
            !if (dfcovering**2 .gt. dftime/(df_rphot*ic)) then
            !    penalties(j) = 1
            !endif

            !dfcovering = 1.d0 - dexp(-dfcovering)

            if (df_rphot .lt. dfri) then
                penalties(j) = 1
                cycle
            endif

            ! If the photosphere is larger than the light travel time, exclude this
            ! solution, as it would imply super-luminal motion
            if ((df_rphot - dfri) .gt. c*(dftime - first_accretion_time)) then
                penalties(j) = 1
                cycle
            endif

        endif

        ! Add photosphere to zones.
        if (include_disk) then
            select case (reprocess_model)
                case (RM_CLOUDS)
                    if (df_rphot .gt. minval(zone_locs(1:n_zones+1)) .and. &
                        df_rphot .lt. maxval(zone_locs(1:n_zones+1))) then
                        zone_locs(n_zones + 2) = df_rphot
                        n_zones = n_zones + 1
                        call sort(zone_locs(1:n_zones+1))
                    endif
            end select
        endif

        ! Calculate the amount of flux that should be reprocessed.
        reprocessed_lum = 0.d0
        if (include_disk) then
            make_sed = .false.
            dfreprocess = .false.
            dfband = 'Lb'

            do i = 1, n_zones
                ierr = 0
                dfcnt = dfcnt + 1
                zmin = dlog10(zone_locs(i))
                zmax = dlog10(zone_locs(i+1))

                call integrate_df(src_df_func, zmin, zmax, mag, df_int_divs**2)

                reprocessed_lum = reprocessed_lum + mag

                if (zone_locs(i+1) .ge. df_rphot) exit
            enddo
            reprocessed_lum = reprocessed_lum * dfcovering
            dfreprocess = .true.
        endif

        !if (print_extra) then
        !    print *, 'Reprocessed luminosity from disk', reprocessed_lum, dfcovering
        !endif

        ! Currently, the extra energy radiated by the relativistic outflow
        ! doesn't contribute directly to the SED in any way, only what is
        ! reprocessed.
        !if (trial_fout(cur_event) .gt. 0.d0) then

        ! Wind_phot and circular luminosity currently mutually exclusive
        if (circ_phot) then
            reprocessed_lum = circ_lum
        elseif (.not. wind_phot) then
            reprocessed_lum = reprocessed_lum + circ_lum
        endif
        !endif

        if (circ_phot .and. wind_phot) then
            print *, 'Error: circ_phot and wind_phot mutually exclusive'
            call exit(0)
        endif


        if (wind_phot .and. trial_cap_at_edd(cur_event)) then
            !if (print_extra .and. max(total_lum - trial_ledd(cur_event),0.d0) .gt. 0.d0) then
            !    print *, 'Reprocessed luminosity from wind', &
            !        trial_fout(cur_event)*max(total_lum - trial_ledd(cur_event),0.d0)
            !    print *, total_lum, trial_ledd(cur_event)
            !endif

            ! This is for only using the Eddington excess for the extra
            ! luminosity
            !print *, 'Error: need to split fout as it has multiple meanings.'
            !call exit(0)

            reprocessed_lum = reprocessed_lum + &
                trial_fout(cur_event)*max(total_lum - trial_ledd(cur_event),0.d0)
        endif
        ! End reprocessed calculation

        !write(*, '(7G17.5)') trial_opacity(cur_event), im2(1), trial_ms(cur_event), df_rphot, rconst, dfcovering, reprocessed_lum

        dfband = bands(j)

        if (produce_seds .and. dfband .eq. 'Lb') then
            sed_index = sed_index + 1
        endif

        if (include_disk) then
            select case (reprocess_model)
                case (RM_CLOUDS)
                    teff = (reprocessed_lum/(fourpi_sigma_b*df_rphot**2))**0.25d0
                    if (teff .ne. 0.d0) then
                        if (restframe_mode) then
                            mags(j) = mags(j) + fourpi*df_rphot**2*&
                                bbflux(src_bb_func, dfband, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                        else
                            mags(j) = mags(j) + fourpi*df_rphot**2*&
                                bbflux(obs_bb_func, dfband, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                        endif
                    endif
            end select
        else
            teff = (trial_eps_edd(cur_event)*dffb*c2*(trial_r_isco(cur_event)/dfri)**3/(fourpi_sigma_b*df_rphot**2))**0.25d0

            if (restframe_mode) then
                if (dfcovering .gt. 0.d0) then
                    mags(j) = mags(j) + fourpi*dfcovering*incl_corr*df_rphot**2*bbflux(src_bb_func, dfband, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                endif
            else
                if (dfcovering .gt. 0.d0) then
                    mags(j) = mags(j) + fourpi*dfcovering*incl_corr*df_rphot**2*bbflux(obs_bb_func, dfband, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                endif
            endif
        endif

        if (include_disk) then
            do i = 1, n_zones
                ierr = 0
                dfcnt = dfcnt + 1
                zmin = dlog10(zone_locs(i))
                zmax = dlog10(zone_locs(i+1))

                make_sed = .false.

                if (restframe_mode) then
                    call integrate_df(src_df_func, zmin, zmax, mag)
                else
                    call integrate_df(obs_df_func, zmin, zmax, mag)
                endif

                if (produce_seds .and. dfband .eq. 'Lb') then
                    ! Simple trapezoid, but for an array
                    coeff = (zmax - zmin)/(dble(sed_divs) - 1.d0) 
                    make_sed = .true.
                    do k = 1, sed_divs
                        ! zone_sed_table is being generated here through data module
                        dummy = src_df_func(zmin + coeff*(dble(k) - 1.d0))
                        if (k .eq. 1 .or. k .eq. sed_divs) then
                            sed_table(1,sed_index,:) = sed_table(1,sed_index,:) + 0.5d0*zone_sed_table*coeff
                        else
                            sed_table(1,sed_index,:) = sed_table(1,sed_index,:) + zone_sed_table*coeff
                        endif
                        dummy = arv_df_func(zmin + coeff*(dble(k) - 1.d0))
                        if (k .eq. 1 .or. k .eq. sed_divs) then
                            sed_table(2,sed_index,:) = sed_table(2,sed_index,:) + 0.5d0*zone_sed_table*coeff
                        else
                            sed_table(2,sed_index,:) = sed_table(2,sed_index,:) + zone_sed_table*coeff
                        endif
                    enddo
                endif

                if (restframe_mode) then
                    mags(j) = mags(j) + mag
                else
                    mags(j) = mags(j) + mag*incl_corr
                endif
            enddo
        endif

        if (produce_seds .and. dfband .eq. 'Lb') then
            if (include_disk) sed_table(2,sed_index,:) = sed_table(2,sed_index,:)*incl_corr
            if (.not. include_disk) then
                if (dfcovering .gt. 0.d0) then
                    sed_table(1,sed_index,:) = sed_table(1,sed_index,:) + &
                        fourpi*dfcovering*df_rphot**2*bbsed(src_bb_func, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                    sed_table(2,sed_index,:) = sed_table(2,sed_index,:) + &
                        fourpi*dfcovering*df_rphot**2*bbsed(arv_bb_func, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                endif
            elseif (reprocess_model .eq. RM_CLOUDS) then
                sed_table(1,sed_index,:) = sed_table(1,sed_index,:) + &
                    fourpi*df_rphot**2*&
                    bbsed(src_bb_func, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
                sed_table(2,sed_index,:) = sed_table(2,sed_index,:) + &
                    fourpi*df_rphot**2*&
                    bbsed(arv_bb_func, teff, trial_z(cur_event), event_nh(cur_event), trial_nhsrc(cur_event))
            endif
        endif

    enddo

    where (mags .le. 0.d0 .and. fbs .ne. 0.d0)
        mags = 0.d0
        penalties = 2
    endwhere

    if (restframe_mode) return
        
    do i = 1, size(bands)
        band_types(i) = get_band_type(bands(i))
    enddo

    where (band_types .eq. 'O')
        mags = penalty*ftoABmag(mags)
    endwhere

    where (band_types .eq. 'l' .or. band_types .eq. 'X')
        mags = penalty*ftomag(mags)
    endwhere

    where ((mags .ge. huge(1.d0) .or. mags .ne. mags) .and. fbs .ne. 0.d0)
        penalties = 2
    endwhere

    where (bands .eq. 'X1')
        mags = mags + trial_offset_X1(cur_event)
    endwhere
    where (bands .eq. 'X2')
        mags = mags + trial_offset_X2(cur_event)
    endwhere
    where (bands .eq. 'GN')
        mags = mags + trial_offset_GN(cur_event)
    endwhere
    where (bands .eq. 'Pg')
        mags = mags + trial_offset_Pg(cur_event)
    endwhere
    where (bands .eq. 'Pr')
        mags = mags + trial_offset_Pr(cur_event)
    endwhere
    where (bands .eq. 'Pi')
        mags = mags + trial_offset_Pi(cur_event)
    endwhere
    where (bands .eq. 'Pz')
        mags = mags + trial_offset_Pz(cur_event)
    endwhere
    where (bands .eq. 'U1')
        mags = mags + trial_offset_U1(cur_event)
    endwhere
    where (bands .eq. 'U2')
        mags = mags + trial_offset_U2(cur_event)
    endwhere
    where (bands .eq. 'RO')
        mags = mags + trial_offset_RO(cur_event)
    endwhere
    where (bands .eq. 'Ub')
        mags = mags + trial_offset_Ub(cur_event)
    endwhere
    where (bands .eq. 'Um')
        mags = mags + trial_offset_Um(cur_event)
    endwhere
    where (bands .eq. 'Uu')
        mags = mags + trial_offset_Uu(cur_event)
    endwhere
    where (bands .eq. 'Uv')
        mags = mags + trial_offset_Uv(cur_event)
    endwhere
    where (bands .eq. 'bV')
        mags = mags + trial_offset_bV(cur_event)
    endwhere

    if (print_max_disk) then
        print *, 'dfmaxdisktemp', dfmaxdisktemp
    endif
end subroutine
