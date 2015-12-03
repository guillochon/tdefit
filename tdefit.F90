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

program tdefit
#include "tdefit.fpp"
    use tdefit_data
    use tdefit_interface, ONLY: dmdt, radius, magdev, bandmag, load_event, ydev, &
                                set_trial_vars, md_arr, annealydev, &
                                set_derived_trial_vars, get_var, print_trial_vars, acor, &
                                draw_random_walker, likelihood, load_defaults, &
                                load_user_vars
    use tdefit_util, only: amoeba_anneal
#ifndef MPIH
    use mpi
#endif


#ifdef MPIH
#include "mpif.h"
#endif

    integer :: ierr, ret_c, i, j, k, l, m, n, ndiscrete, e, modelcnt, rcount, nit, iter, iiter, &
               fn, wp, nwalkers_in, nvars_in, nburn_in, var_type_in, annealloc, mm, &
               pref_fail, w1, w2, nrejections, gnrejections, wchainl, wchainr, wchainn, &
               gdfcnt, gdffailcnt, gbbcnt, gbbfailcnt, gbbcalls, ninitfails, nstationary, ranj, &
               nreplacements, gnreplacements, event_max_npts, event_max_blrpts, &
               event_max_nbest_bands

    !real, dimension(n_nhs) :: nhsrcs
    real, allocatable, dimension(:) :: praw, chi2raw, yraw, pcombinedraw, psend, chi2best_p, ybest_p, &
                                                   yanneal, pbanneal, pbbanneal, new_p, chi2send, ysend, &
                                                   atime, var_dbles
    real, allocatable, dimension(:,:) :: y, chi2, ystep, chi2step, walkerdraw, panneal
    real, allocatable, dimension(:,:,:) :: p, walkerout, pcombined, yall, chi2all, yburn, chi2burn, pstep, walkerin
    real, allocatable, dimension(:,:,:,:) :: pall, pburn
    real :: chi2best, ybest, dummy, dummy2, penalty, pref_dev, zstretch, q, accept_frac, &
                        chi2old, yold, chi2best_y, ybest_chi2, lzscale, min_search_in, max_search_in, ybanneal, ybbanneal, &
                        temptr, anneal_y_best_initial, new_chi2, new_y, &
                        lmin_trial_time, lmax_trial_time, ybeststep, replace_frac
    real, allocatable, dimension(:) :: fat_most_probable, fat_best_fit
    integer, allocatable :: seed_arr(:)
    integer :: seed_size, version
    integer, dimension(2) :: bestloc
    integer*8, dimension(1) :: mpi_int_in, mpi_int_out
    integer, allocatable, dimension(:,:) :: scount

    character(len=128) :: str
    integer, parameter :: cur_version = 1

    logical :: redraw, force_accept, file_exists, replace

    fn = 11
    verbose = .false.

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_pe, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_pe, ierr)

    if (red_blue) then
        if (mod(n_pe, 2) .ne. 0) then
            write(*, *), 'Warning: TDEFit must be run with an even number of processors in red/blue mode, switching off red/blue.'
        endif
    endif

    call random_seed(size = seed_size)
    allocate(seed_arr(seed_size))
    call random_seed(get = seed_arr)
    seed_arr = (my_pe+1)*seed_arr
    call random_seed(put = seed_arr)

    call init                                                                                   

    if (mod(nwalkers, nchains) .ne. 0) then
        print*, 'Error: Number of walkers must be evenly divisible by number of chains'
        call exit(0)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    open(unit = fn, file = "event_list.dat", status='old', action='read')
    read(fn, *) event_n
    allocate(event_fnames(event_n))
    allocate(event_npts(event_n))
    allocate(event_nbest_bands(event_n))
    allocate(event_blrpts(event_n))
    allocate(event_nh(event_n))
    allocate(event_claimed_z(event_n))
    allocate(event_min_aspin(event_n))
    allocate(event_nhcorr(event_n))
    allocate(event_restframe(event_n))
    allocate(trial_mh(event_n))
    allocate(trial_ms(event_n))
    allocate(trial_rs(event_n))
    allocate(trial_beta(event_n))
    allocate(trial_aspin(event_n))
    allocate(trial_nhsrc(event_n))
    allocate(trial_phi(event_n))
    allocate(trial_bh_rms(event_n))
    allocate(trial_rp(event_n))
    allocate(trial_toff(event_n))
    allocate(trial_rout(event_n))
    allocate(trial_rg(event_n))
    allocate(trial_gmh(event_n))
    allocate(trial_ms0(event_n))
    allocate(trial_mh0(event_n))
    allocate(trial_rs0(event_n))
    allocate(trial_rsc(event_n))
    allocate(trial_alphhr(event_n))
    allocate(trial_eps_edd(event_n))
    allocate(trial_source_rv(event_n))
    allocate(trial_magoff(event_n))
    allocate(trial_fcor(event_n))
    allocate(trial_tlimit(event_n))
    allocate(trial_ecor(event_n))
    allocate(trial_r_isco(event_n))
    allocate(trial_reprocess_temp(event_n))
    allocate(trial_offset_X1(event_n))
    allocate(trial_offset_X2(event_n))
    allocate(trial_offset_GN(event_n))
    allocate(trial_offset_Pg(event_n))
    allocate(trial_offset_Pr(event_n))
    allocate(trial_offset_Pi(event_n))
    allocate(trial_offset_Pz(event_n))
    allocate(trial_offset_U1(event_n))
    allocate(trial_offset_U2(event_n))
    allocate(trial_offset_RO(event_n))
    allocate(trial_offset_Ub(event_n))
    allocate(trial_offset_Um(event_n))
    allocate(trial_offset_Uu(event_n))
    allocate(trial_offset_Uv(event_n))
    allocate(trial_offset_bV(event_n))
    allocate(trial_fout(event_n))
    allocate(trial_mu_e(event_n))
    allocate(trial_ledd(event_n))
    allocate(trial_outflow_frac(event_n))
    allocate(trial_bh_rbs(event_n))
    allocate(trial_r_ibco(event_n))
    allocate(trial_temp_mult(event_n))
    allocate(trial_z(event_n))
    allocate(trial_1pz(event_n))
    allocate(trial_dl(event_n))
    allocate(trial_exp_1(event_n))
    allocate(trial_exp_2(event_n))
    allocate(trial_exp_3(event_n))
    allocate(trial_exp_4(event_n))
    allocate(trial_variability(event_n))
    allocate(trial_variance(event_n))
    allocate(trial_viscous_time(event_n))
    allocate(trial_variability2(event_n))
    allocate(trial_variance2(event_n))
    allocate(trial_yms(event_n))
    allocate(trial_y1(event_n))
    allocate(trial_y2(event_n))
    allocate(trial_y3(event_n))
    allocate(trial_opacity(event_n))
    allocate(trial_rphot(event_n))
    allocate(trial_model(event_n))
    allocate(trial_outflow_model(event_n))
    allocate(trial_object_type(event_n))
    allocate(trial_temperature_model(event_n))
    allocate(trial_blr_model(event_n))
    allocate(trial_time_dep_rin(event_n))
    allocate(trial_time_dep_rout(event_n))
    allocate(trial_simple_bb(event_n))
    allocate(trial_use_fcor(event_n))
    allocate(trial_cap_at_edd(event_n))
    allocate(trial_full_disk_coverage(event_n))
    allocate(ll(event_n))
    do e = 1, event_n
        read(fn, *) event_fnames(e)
    enddo
    close(fn)
    allocate(fat_best_fit(event_n))
    allocate(fat_most_probable(event_n))

    call tdefit_print('Loading events')

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

    allocate(trial_times(event_max_npts,event_n))
    allocate(trial_fbs(event_max_npts,event_n))
    allocate(trial_mdots(event_max_npts,event_n))
    allocate(trial_menv(event_max_npts,event_n))
    allocate(trial_routs(event_max_npts,event_n))
    allocate(trial_rphots(event_max_npts,event_n))
    allocate(trial_mags(event_max_npts,event_n))

    nvars = 0
    do e = 1, event_n
        call load_event(e,.false.)
        call set_event(e)
        call load_defaults(1)
        call load_user_vars(1)
    enddo

    event_max_nbest_bands = maxval(event_nbest_bands)
    allocate(event_best_bands(event_max_nbest_bands,event_n))

    do e = 1, event_n
        do i = 1, nextra_bands
            event_best_bands(event_nbest_bands(e) - nextra_bands + i,e) = extra_bands(i)
        enddo
    enddo

    if (sum(event_npts) .le. nvars - 1) then
        print *, 'Warning: Number of measurement points must be +2 larger than ' // &
                 'number of degrees of freedom for proper chi-square measurement.'
    endif

    reduced_chi2_const = 1.d0 / max(sum(event_npts) - nvars - 1 - count(var_types .eq. 2), 1)

    allocate(y(nmcsteps,nwalkers))
    allocate(chi2(nmcsteps,nwalkers))
    allocate(p(nmcsteps,nwalkers,nvars))
    allocate(pcombinedraw(n_pe*nwalkers*nvars))
    allocate(pcombined(nwalkers,nvars,n_pe))
    allocate(scount(count(var_types .eq. 2), &
        nint(maxval(max_search - min_search, var_types .eq. 2)) + 1))
    allocate(panneal(nvars+1,nvars))
    allocate(pbanneal(nvars))
    allocate(pbbanneal(nvars))
    allocate(yanneal(nvars+1))
    allocate(locked_trial_vars(nvars))
    allocate(var_locks(nvars))
    allocate(new_p(nvars))
    allocate(atime(nvars))
    allocate(var_dbles(size(all_var_types)))

    var_locks = .false.

    if (my_pe .eq. 0) then
        allocate(pstep(nwalkers, nvars, n_pe))
        allocate(chi2step(nwalkers, n_pe))
        allocate(ystep(nwalkers, n_pe))
        allocate(chi2best_p(nvars))
        allocate(ybest_p(nvars))
        if (nvars .eq. 0) then
            write(*, *), 'No variables are set to be minimized over, immediately producing output.'
        endif
    endif

    call tdefit_print('Initialized arrays.')

    if (preferred_mode) then
        open(unit = fn, file = trim(event_fnames(e)) // "_preferred.dat", status='old', action='read')
        read(fn,*) pref_dev, trial_toff(cur_event), trial_nhsrc(cur_event), trial_rout(cur_event), trial_aspin(cur_event), trial_phi(cur_event), trial_beta(cur_event), &
            trial_mh(cur_event), trial_ms(cur_event), trial_rsc(cur_event), trial_alphhr(cur_event)
        pref_fail = 0
        close(fn)

        call set_derived_trial_vars
        write(*, *), pref_dev, trial_toff(cur_event), trial_nhsrc(cur_event), trial_rout(cur_event), trial_aspin(cur_event), trial_phi(cur_event), trial_beta(cur_event), &
            trial_mh(cur_event), trial_ms(cur_event), trial_rsc(cur_event), trial_alphhr(cur_event)
    else
        modelcnt = 0

        call tdefit_print('Beginning annealing process...')

        modelcnt = modelcnt + 1

        if (my_pe .eq. 0) then
            allocate(praw(n_pe*nwalkers*nvars))
            allocate(chi2raw(n_pe*nwalkers))
            allocate(yraw(n_pe*nwalkers))
        else
            allocate(praw(1))
            allocate(chi2raw(1))
            allocate(yraw(1))
        endif
        allocate(psend(nwalkers*nvars))
        allocate(chi2send(nwalkers))
        allocate(ysend(nwalkers))

        accept_frac = 1.d0

        chi2best = huge(1.d0)
        ybest = -huge(1.d0)
        ybeststep = -huge(1.d0)

        if (nvars .eq. 0) then
            chi2best = magdev(p(1,1,:))
            ybest = ydev(p(1,1,:))
            ybest_chi2 = chi2best
            chi2best_y = ybest
        elseif (initial_walker_dist .eq. IW_STATISTICAL .or. &
            initial_walker_dist .eq. IW_PROBABLE) then
            if (initial_walker_dist .eq. IW_STATISTICAL) then
                call tdefit_print('Drawing initial walkers from walker distribution.')
            elseif (initial_walker_dist .eq. IW_PROBABLE) then
                call tdefit_print('Drawing initial walkers based on priors.')
            endif
            open(unit = fn, file = trim(output_path) // 'burned_walkers' // '.dat', status = 'unknown', action = 'read', form = 'unformatted')
            read(fn) nburn_in, nwalkers_in, nvars_in
            if (nburn_in .lt. 0) then
                write(*, *), 'Error, nburn_in must be positive when using statistical or probable'
                call exit(0)
            endif
            if (nvars .ne. nvars_in) then
                if (my_pe .eq. 0) write(*, *), 'Error, number of variables to minimize over inconsistent.'
                call exit(0)
            endif
            do i = 1, nvars
                read(fn) min_search_in, max_search_in, var_type_in
                if (min_search_in .ne. min_search(i) .or. max_search_in .ne. max_search(i) .or. var_type_in .ne. var_types(i)) then
                    if (my_pe .eq. 0) write(*, *), 'Error, search domain inconsistent.'
                    call exit(0)
                endif
            enddo

            allocate(walkerin(nburn_in,nwalkers_in,nvars_in+2))
            allocate(walkerdraw(nwalkers_in*nburn_in,nvars_in+1))

            read(fn) (((walkerin(i,j,k),k=1,nvars_in+2),j=1,nwalkers_in),i=1,nburn_in)
            close(fn)

            if (initial_walker_dist .eq. IW_PROBABLE) then
                do i = 1, nburn_in
                    do j = 1, nwalkers_in
                        walkerdraw((i-1)*nwalkers_in + j,1:) = walkerin(i,j,2:)
                    enddo
                enddo
                i = maxloc(walkerdraw(:,1),1)
                print *, walkerdraw(i,1)
                p(1,1,:) = walkerdraw(i,2:)
                chi2(1,1) = magdev(p(1,1,:))
                y(1,1) = ydev(p(1,1,:))
                do w1 = 2, nwalkers
                    p(1,w1,:) = p(1,1,:)
                    chi2(1,w1) = chi2(1,1)
                    y(1,w1) = y(1,1)
                enddo
                print *, likelihood()
                print *, y(1,1), chi2(1,1)
            else
                walkerdraw(1,1) = dexp(walkerin(1,2,1))
                do i = 1, nburn_in
                    do j = 1, nwalkers_in
                        if (i .eq. 1 .and. j .eq. 1) cycle
                        walkerdraw((i-1)*nwalkers_in + j,1) = walkerdraw((i-1)*nwalkers_in + j - 1,1) + dexp(walkerin(i,j,2))
                        walkerdraw((i-1)*nwalkers_in + j,2:) = walkerin(i,j,3:)
                    enddo
                enddo

                walkerdraw(:,1) = walkerdraw(:,1)/walkerdraw(nwalkers_in*nburn_in,1)
                walkerdraw(:,1) = walkerdraw(:,1) - walkerdraw(1,1)

                do w1 = 1, nwalkers
                    call random_number(dummy)
                    do i = nwalkers_in*nburn_in, 1, -1
                        if (dummy .lt. walkerdraw(i,1)) cycle
                        p(1,w1,:) = walkerdraw(i,2:)
                        chi2(1,w1) = magdev(p(1,w1,:))
                        y(1,w1) = ydev(p(1,w1,:))
                        exit
                    enddo
                enddo
            endif

            deallocate(walkerdraw)
            deallocate(walkerin)
        else
            if (initial_walker_dist .eq. IW_REGULAR) then
                call tdefit_print('Drawing initial walkers on regular grid.')
            endif
            do w1 = 1, nwalkers
#ifdef REJECT_INITIAL_FAILS
                ninitfails = -1
                dfreject = .true.
                do while (dfreject)
                    ninitfails = ninitfails + 1
#endif
                    do i = 1, nvars
                        ! Only works for one variable for now
                        if (initial_walker_dist .eq. IW_REGULAR) then
                            if (nvars .ne. 1) then
                                call tdefit_print('Error: IW_REGULAR only supports one free variable')
                                call exit(0)
                            endif
                            p(1,w1,i) = dble(w1-1)/dble(nwalkers-1)
                        else
                            call random_number(dummy)
                            p(1,w1,i) = dummy
                        endif
                        if (var_types(i) .eq. 2 .or. var_types(i) .eq. 3) then
                            p(1,w1,i) = dble(floor(p(1,w1,i)*(max_search(i) - min_search(i) + 1.d0)))/(max_search(i) - min_search(i))
                        endif
                    enddo
                    chi2(1,w1) = magdev(p(1,w1,:))

                    if (initial_walker_dist .eq. IW_REGULAR) then
                        call tdefit_print('IW_REGULAR permits initial fails.')
                        exit
                    endif
#ifdef DEBUG
                    if (dfreject) write(*, *), 'ninitfails', ninitfails, my_pe
#endif
#ifdef REJECT_INITIAL_FAILS
                enddo
#endif
                !if (nanneal .eq. 0) then
                    y(1,w1) = ydev(p(1,w1,:))
                !else
                !    y(1,w1) = (likelihood() - chi2(1,w1))/annealt0
                !endif
                !if (chi2(1,w1) .gt. 2.d0 + dble(nvars)) then
                !    y(1,w1) = likelihood()*gcf(0.5d0*dble(nvars), 0.5d0*chi2(1,w1), dummy)
                !else
                !    y(1,w1) = likelihood()*gser(0.5d0*dble(nvars), 0.5d0*chi2(1,w1), dummy)
                !endif
            enddo

#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering [1a]'
#endif
            rcount = nwalkers*nvars
            psend = reshape(p(1,:,:), [rcount])
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(psend, rcount, MPI_DOUBLE_PRECISION, praw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering finished [1a]'
            if (my_pe .eq. 0) write(*, *), 'Gathering [1b]'
#endif
            rcount = nwalkers
            chi2send = reshape(chi2(1,:), [rcount])
            ysend = reshape(y(1,:), [rcount])
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(chi2send, rcount, MPI_DOUBLE_PRECISION, chi2raw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(ysend, rcount, MPI_DOUBLE_PRECISION, yraw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering finished [1b]'
#endif
            if (my_pe .eq. 0) then
                pstep = reshape(praw, [nwalkers, nvars, n_pe])
                chi2step = reshape(chi2raw, [nwalkers, n_pe])
                ystep = reshape(yraw, [nwalkers, n_pe])
                ybeststep = maxval(ystep)

                bestloc = minloc(chi2step)
                chi2best_p = pstep(bestloc(1),:,bestloc(2))
                chi2best_y = ystep(bestloc(1),bestloc(2))
                chi2best = minval(chi2step)

                bestloc = maxloc(ystep)
                ybest_p = pstep(bestloc(1),:,bestloc(2))
                ybest_chi2 = chi2step(bestloc(1),bestloc(2))
                ybest = ybeststep
            endif
        endif

        ndiscrete = 0
        do i = 1, nvars
            if (var_types(i) .eq. 2) then
                ndiscrete = ndiscrete + 1
                do j = nint(min_search(i)) + 1, nint(max_search(i)) + 1
#ifdef VERBOSE_MPI
                    if (my_pe .eq. 0) write(*, *), 'Reducing [1]'
#endif
                    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                    call MPI_ALLREDUCE(count(p(1,:,i) .eq. (dble(j-1) - min_search(i))/(max_search(i) - min_search(i))), &
                        scount(ndiscrete,j-nint(min_search(i))), 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
                    if (my_pe .eq. 0) write(*, *), 'Reducing finished [1]'
#endif
                enddo
            endif
        enddo

        do nstep = 2, nmcsteps
            bbcalls = 0
            if (.not. discard_failed_integrals) then
                dfcnt = 0
                dffailcnt = 0
                bbcnt = 0
                bbfailcnt = 0
            endif

            nrejections = 0
            nreplacements = 0

            if (my_pe .eq. 0) write(*, *), 'Step: ', nstep

            rcount = nwalkers*nvars
            psend = reshape(p(nstep-1,:,:), [rcount])

            ! Stop if a file named "stop" is present.
            inquire(file="stop", exist=file_exists)
            if (file_exists) then
                nmcsteps = nstep - 1
                nburn = min(nburn, nmcsteps)
                burn_in = nmcsteps - nburn + 1
                exit
            endif
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering [2]'
#endif
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_ALLGATHER(psend, rcount, MPI_DOUBLE_PRECISION, &
                pcombinedraw, rcount, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            write(*, *), 'Gathering finished [2]', my_pe
#endif
            pcombined = reshape(pcombinedraw, [nwalkers, nvars, n_pe])

            do w1 = 1, nwalkers
                ! Randomly redraw walkers that have been stationary for a long time, or are significantly worse
                ! than the all-time best. Only do this during annealing phase.
                wchainl = int(dble(w1-1)/(nwalkers/nchains)) * nwalkers/nchains + 1
                wchainr = wchainl + nwalkers/nchains - 1
                wchainn = floor((dble(wchainl) - 0.5d0)/nchains)

                if (nanneal .eq. 0) then
                    mctemp = 1.d0
                else
                    ! Single temp chains
                    mctemp = max(annealt0*(1.d0/annealt0)**(dble(nstep)/nanneal/anneal_hot_frac), 1.d0)
                    ! Multi temp chains
                    !mctemp = 10.d0**((dble(wchainn) - 1.d0)/(dble(nchains) - 1.d0)*&
                    !    max(dlog10(annealt0*(1.d0/annealt0)**(dble(nstep)/nanneal)), 0.d0))
                endif

                redraw = .false.
                force_accept = .false.

                !if (redraw_stationary) then
                !    if (nstep .gt. stationary_redraw .and. nstep .le. nanneal) then
                !        if (y(nstep-1,w1) .eq. y(nstep-stationary_redraw,w1)) then
                !            !write(*, *), 'Redrew a stationary walker.'
                !            if (redraw_using_best) then
                !                annealloc = maxloc(y(nstep-1,:),1)
                !                p(nstep-1,w1,:) = p(nstep-1,annealloc,:)
                !                chi2(nstep-1,w1) = chi2(nstep-1,annealloc)
                !                y(nstep-1,w1) = y(nstep-1,annealloc)
                !                !force_accept = .true.
                !            else
                !                do i = 1, nvars
                !                    call random_number(dummy)
                !                    p(nstep,w1,i) = dummy
                !                    if (var_types(i) .eq. 2 .or. var_types(i) .eq. 3) then
                !                        p(nstep,w1,i) = dble(floor(p(nstep,w1,i)*&
                !                            (max_search(i) - min_search(i) + 1.d0)))/(max_search(i) - min_search(i))
                !                    endif
                !                enddo
                !                redraw = .true.
                !            endif
                !        endif
                !    endif
                !endif

                !if (redraw_bad_walkers) then
                !    if ((ybest - y(nstep-1,w1))/mctemp .gt. dlog(dble(nwalkers*n_pe)) + 2.d0) then
                !        call draw_random_walker(w1, wchainl, wchainr, w2, wp)
                !        print *, 'Redrew bad walker.'
                !        p(nstep-1,w1,:) = pcombined(w2,:,wp)
                !    endif
                !endif

                !if (ybeststep .eq. y(nstep-1,w1)) then
                !    print *, 'Found best walker.'
                !    redraw = .true.
                !    force_accept = .true.
                !endif

                ! Select a random walker
                if (.not. redraw) then
#ifndef REJECT_OUT_OF_BOUNDS
                    p(nstep,w1,:) = -1.d0
                    do while (any(p(nstep,w1,:) .gt. 1.d0) .or. any(p(nstep,w1,:) .lt. 0.d0))
#endif
                        call draw_random_walker(w1, wchainl, wchainr, w2, wp)
                        !do nstationary = 0, nstep-3
                        !    if (nstep .le. 2) exit
                        !    if (y(nstep-2-nstationary,w1) .ne. y(nstep-1,w1)) exit
                        !enddo

                        !if (nstep .lt. nanneal) then
                        !    lzscale = (zscale - 1.d0)*dexp(-dble(nstationary)/scale_exp_fac) + 1.d0
                        !else
                            lzscale = zscale
                        !endif

                        ! Below is from "Markov Chain Monte Carlo and Gibbs Sampling" by B. Walsh 2004
                        call random_number(dummy)

                        ! Exclude large reverse stretches
                        !dummy = dummy*min(1.d0, (1.d0 - lzscale + sqrt2*dsqrt(lzscale - 2.d0*lzscale**2 + lzscale**3))/&
                        !    (1.d0 - 2.d0*lzscale + lzscale**2.d0))

                        !zstretch = 1.d0/zscale + 0.25d0*dummy*(4.d0/dsqrt(zscale) + dummy)
                        zstretch = (dummy*(lzscale - 1.d0) + 1.d0)**2.d0/lzscale

                        replace = .false.
                        if (nstep .le. nanneal) then
                            if (redraw_bad_walkers) then
                                if ((ybeststep - y(nstep-1,w1))/mctemp .gt. dlog(dble(nwalkers*n_pe)) + dlog(lzscale)*(dble(nvars) - 1.d0)) then
                                    !print *, 'Replacing bad walker.'
                                    replace = .true.
                                endif
                            endif

                            if (redraw_stationary .and. nstep .gt. stationary_redraw) then
                                if (y(nstep-1,w1) .eq. y(nstep-stationary_redraw,w1)) then
                                    !print *, 'Replacing stationary walker'
                                    replace = .true.
                                    !force_accept = .true.
                                endif
                            endif
                        endif

                        ! Draw new point
                        do i = 1, nvars
                            if (var_types(i) .ne. 2 .and. var_types(i) .ne. 3) then
                                !p(nstep,w1,i) = pother(w2,i) + zstretch*(p(nstep-1,w1,i) - pother(w2,i))
                                if (replace) then
                                    ! Trying "replacement" instead of "stretch"
                                    p(nstep,w1,i) = p(nstep-1,w1,i) + zstretch*(pcombined(w2,i,wp) - p(nstep-1,w1,i))
                                else
                                    p(nstep,w1,i) = pcombined(w2,i,wp) + zstretch*(p(nstep-1,w1,i) - pcombined(w2,i,wp))
                                endif
                            else
                                ! Occaisionally select another model. This is to guarantee that there are always
                                ! some walkers present in all models.
                                !if (p(nstep-1,w1,i) .eq. pother(w2,i)) then
                                if (p(nstep-1,w1,i) .eq. pcombined(w2,i,wp)) then
                                    call random_number(dummy2)
                                    if (dummy2 .le. max(mutation_min, &
                                        max(0.d0,dble(nanneal - nstep)/nanneal)*mutation_prob/ndiscrete)) then
                                        call random_number(dummy2)
                                        p(nstep,w1,i) = floor(dummy2*(1.d0 + max_search(i) - min_search(i)))/&
                                            (max_search(i) - min_search(i))
                                    else
                                        p(nstep,w1,i) = p(nstep-1,w1,i)
                                    endif
                                else
                                    if (replace) then
                                        p(nstep,w1,i) = pcombined(w2,i,wp)
                                    else
                                        call random_number(dummy2)
                                        if (dummy2 .lt. 1.d0/lzscale) then
                                        ! Compare to 50% for CDF
                                        !if (zstretch .gt. 0.25d0*(1.d0 + lzscale)**2/lzscale) then
                                            p(nstep,w1,i) = p(nstep-1,w1,i)
                                        else
                                            !p(nstep,w1,i) = pother(w2,i)
                                            p(nstep,w1,i) = pcombined(w2,i,wp)
                                        endif
                                    endif
                                endif
                            endif
                        enddo

                        if (replace) then
                            nreplacements = nreplacements + 1
                        endif
#ifndef REJECT_OUT_OF_BOUNDS    
                    enddo
#endif
                endif

#ifdef REJECT_OUT_OF_BOUNDS
                if (all(p(nstep,w1,:) .le. 1.d0) .and. all(p(nstep,w1,:) .ge. 0.d0)) then
#endif
                    if (hard_penalties .eq. 1 .and. nstep .le. nanneal) then
                        chi2old = magdev(p(nstep-1,w1,:))
                        yold = ydev(p(nstep-1,w1,:))
                    else
                        yold = y(nstep-1,w1)
                    endif
                    chi2(nstep,w1) = magdev(p(nstep,w1,:))
                    y(nstep,w1) = ydev(p(nstep,w1,:))

                    ! accept_frac from Feroz and Hobson 2008
                    q = dlog(zstretch)*(dble(nvars) - 1.d0) + (y(nstep,w1) - yold)/mctemp
                    if (boost_accepts) then
                        if (nstep .le. nanneal) q = q - &
                            min(0.d0, 20.d0*(dble(nanneal - nstep)/nanneal)*(accept_frac - 0.2d0))
                    endif
                    q = min(0.d0, q)

                    ! Small chance to accept bad locations during annealing phase
                    !q = max(-2.d0/max(0.d0, dble(nanneal - nstep) / nanneal), q)

                    call random_number(dummy)
                    !write(*, '(A,5ES12.5)'), 'q comps', dlog(zstretch)*(dble(nvars) - 1.d0), y(nstep,w1), yold, (y(nstep,w1) - yold)/mctemp, dlog(dummy)
                    if (.not. force_accept .and. (trial_penalty .ge. huge(1.d0) .or. q .lt. dlog(dummy))) then
                        !call random_number(dummy)
                        !if (-dlog(mctemp) .gt. dlog(dummy)) then
                            nrejections = nrejections + 1
                            p(nstep,w1,:) = p(nstep-1,w1,:)
                            y(nstep,w1) = y(nstep-1,w1)
                            chi2(nstep,w1) = chi2(nstep-1,w1)
                        !endif
                    endif
#ifdef REJECT_OUT_OF_BOUNDS
                else
                    nrejections = nrejections + 1
                    p(nstep,w1,:) = p(nstep-1,w1,:)
                    y(nstep,w1) = y(nstep-1,w1)
                    chi2(nstep,w1) = chi2(nstep-1,w1)
                endif
#endif

            enddo

#ifdef VERBOSE_MPI
            write(*, *), 'Gathering [5]', my_pe
#endif
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            mpi_int_in = nrejections
            call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
            gnrejections = mpi_int_out(1)

            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            mpi_int_in = nreplacements
            call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
            gnreplacements = mpi_int_out(1)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering finished [3]'
#endif

            accept_frac = 1.d0 - dble(gnrejections)/n_pe/nwalkers
            replace_frac = dble(gnreplacements)/n_pe/nwalkers

            if (my_pe .eq. 0) then
                write(*, *), 'Acceptance fraction: ', accept_frac
                write(*, *), 'Replacement fraction: ', replace_frac
            endif

            if (.not. discard_failed_integrals) then
#ifdef VERBOSE_MPI
                if (my_pe .eq. 0) write(*, *), 'Reducing [2]'
#endif
                !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                mpi_int_in = dfcnt
                call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
                gdfcnt = mpi_int_out(1)
#ifdef VERBOSE_MPI
                if (my_pe .eq. 0) write(*, *), 'Reducing [3]'
#endif
                !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                mpi_int_in = dffailcnt
                call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
                gdffailcnt = mpi_int_out(1)
#ifdef VERBOSE_MPI
                if (my_pe .eq. 0) write(*, *), 'Reducing [4]'
#endif
                !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                mpi_int_in = bbcnt
                call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
                gbbcnt = mpi_int_out(1)
#ifdef VERBOSE_MPI
                if (my_pe .eq. 0) write(*, *), 'Reducing [5]'
#endif
                !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                mpi_int_in = bbfailcnt
                call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
                gbbfailcnt = mpi_int_out(1)
#ifdef VERBOSE_MPI
                if (my_pe .eq. 0) write(*, *), 'Reducing finished [5]'
#endif

                !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                mpi_int_in = bbcalls
                call MPI_ALLREDUCE(mpi_int_in, mpi_int_out, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
                gbbcalls = mpi_int_out(1)

                if (my_pe .eq. 0) then
                    if (gdfcnt .ne. 0) &
                        write(*, *), 'df fail fraction:  ', dble(gdffailcnt)/gdfcnt
                    if (gbbcnt .ne. 0) &
                        write(*, *), 'bb fail fraction:  ', dble(gbbfailcnt)/gbbcnt
                    write(*, *), 'bb ints performed: ', dble(gbbcnt)
                    write(*, *), 'bb calls per int:  ', dble(gbbcalls)/gbbcnt
                endif
            endif

            mm = 0
            do i = 1, nvars
                if (var_types(i) .eq. 2) then
                    mm = mm + 1
                    do j = nint(min_search(i)) + 1, nint(max_search(i)) + 1
#ifdef VERBOSE_MPI
                        if (my_pe .eq. 0) write(*, *), 'Reducing [6]'
#endif
                        !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                        call MPI_ALLREDUCE(count(p(nstep,:,i) .eq. (dble(j-1) - min_search(i))/(max_search(i) - min_search(i))), &
                            scount(mm,j-nint(min_search(i))), 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
                        if (my_pe .eq. 0) write(*, *), 'Reducing finished [6]'
#endif
                    enddo
                endif
            enddo

            if (amoeba_best .and. mod(nstep-1,amoeba_step) .eq. 0 .and. (nstep - 1) .le. nanneal) then
                !use_soft_penalties = .false.
                if (my_pe .eq. 0) write(*, *), 'Annealing best walkers...'

                var_locks = .false.
                do i = 1, nvars
                    if (var_types(i) .eq. 2 .or. var_types(i) .eq. 3) then
                        var_locks(i) = .true.
                    endif
                enddo

                do j = 1, nchains
                    wchainl = (j - 1) * nwalkers/nchains + 1
                    wchainr = wchainl + nwalkers/nchains - 1

                    call random_number(dummy)
                    ranj = ceiling(dummy*nwalkers/nchains)

                    annealloc = wchainl + maxloc(y(nstep,wchainl:wchainr),1) - 1
                    !annealloc = wchainl + ranj - 1

                    do i = 1, nvars
                        if (var_types(i) .eq. 2 .or. var_types(i) .eq. 3) then
                            locked_trial_vars(i) = p(nstep,annealloc,i)
                        endif
                    enddo

                    panneal = 0.d0
                    panneal(1,:) = p(nstep,annealloc,:)
                    do m = 2, size(panneal, 1)
                        do n = 1, size(panneal, 2)
                            if (var_types(n) .ne. 2 .and. var_types(n) .ne. 3) then
                                call random_number(dummy)
                                panneal(m,n) = max(0.d0,min(1.d0,(1.d0 - amoeba_spread) * &
                                    p(nstep,annealloc,n) + 2.d0*amoeba_spread*dummy))
                            else
                                panneal(m,n) = p(nstep,annealloc,n)
                            endif
                        enddo
                    enddo

                    !use_soft_penalties = .true.

                    do m = 1, size(panneal, 1)
                        yanneal(m) = annealydev(panneal(m,:))
                    enddo

                    anneal_y_best_initial = yanneal(1)

                    ybanneal = huge(1.d0)
                    ybbanneal = huge(1.d0)
                    temptr = mctemp
                    
                    nit = 0
                    iiter = 2*nvars
                    do m = 1, 10
                        iter = iiter
                        temptr = (0.1d0/mctemp)**(1.d0/(nvars))*temptr
                        call amoeba_anneal(panneal, yanneal, pbanneal, ybanneal, amoeba_tol, annealydev, iter, temptr)
                        nit = nit + iiter - iter
                        if (ybanneal .lt. ybbanneal) then
                            ybbanneal = ybanneal
                            pbbanneal = pbanneal
#ifdef PRINT_ANNEALING_INFO
                            write(str,'(I2)') size(pbanneal)
                            write(*,'(1x,i6,i6,e10.3,' // trim(str) //'f11.5,e20.7)')&
                                m,nit,temptr,(pbanneal(n),n=1,size(pbanneal)),-ybanneal
#endif
                        endif
                        if (iter .gt. 0) then
                            exit
                        endif
                    enddo

                    !use_soft_penalties = .false.

                    if (ybbanneal .gt. anneal_y_best_initial) then
                        if (print_annealing) then
                            write (*, '(A,I4,A,I4,A)') 'Annealing on proc: ', my_pe, ', chain: ', j, ' could not improve solution.'
                        endif
                    else
                        new_p = pbbanneal
                        do i = 1, nvars
                            if (var_locks(i)) then
                                new_p(i) = locked_trial_vars(i)
                            endif
                        enddo
                        new_chi2 = magdev(new_p)
                        new_y = ydev(new_p)
                        if (new_y .gt. y(nstep,annealloc) .and. trial_penalty .ne. huge(1.d0)) then
                            if (print_annealing) then
                                write (*, '(A,I4,A,I4,A,ES10.3,A,ES10.3)') 'Annealing on proc: ', my_pe, &
                                    ', chain: ', j, ' complete, old: ', y(nstep,annealloc), ', new: ', new_y
                            endif

                            p(nstep,annealloc,:) = new_p
                            chi2(nstep,annealloc) = new_chi2
                            y(nstep,annealloc) = new_y
                        else
                            if (print_annealing) then
                                write (*, '(A,I4,A,I4,A,ES10.3,A,ES10.3)') 'Annealing on proc: ', my_pe, &
                                    ', chain: ', j, ' did not improve solution because of soft penalties in amoeba.'
                            endif
                        endif
                    endif
                enddo
                var_locks = .false.
                use_soft_penalties = .false.
            endif

            rcount = nwalkers*nvars
            psend = reshape(p(nstep,:,:), [rcount])
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(psend, rcount, MPI_DOUBLE_PRECISION, praw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering [4]'
#endif
            rcount = nwalkers
            chi2send = reshape(chi2(nstep,:), [rcount])
            ysend = reshape(y(nstep,:), [rcount])
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(chi2send, rcount, MPI_DOUBLE_PRECISION, chi2raw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(ysend, rcount, MPI_DOUBLE_PRECISION, yraw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering finished [4]'
#endif

            if (my_pe .eq. 0) then
                pstep = reshape(praw, [nwalkers, nvars, n_pe])
                chi2step = reshape(chi2raw, [nwalkers, n_pe])
                ystep = reshape(yraw, [nwalkers, n_pe])
                ybeststep = maxval(ystep)
                if (hard_penalties .eq. 1 .and. nstep .le. nanneal) then
                    chi2best = huge(1.d0)
                    ybest = -huge(1.d0)
                endif
                if (minval(chi2step) .lt. chi2best) then
                    bestloc = minloc(chi2step)
                    chi2best_p = pstep(bestloc(1),:,bestloc(2))
                    chi2best_y = ystep(bestloc(1),bestloc(2))
                    chi2best = minval(chi2step)
                endif
                if (ybeststep .gt. ybest) then
                    bestloc = maxloc(ystep)
                    ybest_p = pstep(bestloc(1),:,bestloc(2))
                    ybest_chi2 = chi2step(bestloc(1),bestloc(2))
                    ybest = ybeststep
                endif

                print_extra = .true.
                dummy = ydev(ybest_p)
                print_extra = .false.

                if (hard_penalties .eq. 1 .and. nstep .le. nanneal) then
                    write(*, '(A40,ES12.5,A2,ES12.5,A1)'), "Current best (annealing phase): ", chi2best * reduced_chi2_const, &
                        " (", chi2best_y, ")"
                    write(*, '(A40,ES12.5,A2,ES12.5,A1)'), "Current most-probable (annealing phase): ", ybest, &
                        " (", ybest_chi2 * reduced_chi2_const, ")"
                else
                    write(*, '(A25,ES12.5,A2,ES12.5,A1)'), "All-time best: ", chi2best * reduced_chi2_const, &
                        " (", chi2best_y, ")"
                    write(*, '(A25,ES12.5,A2,ES12.5,A1)'), "All-time most-probable: ", ybest, &
                        " (", ybest_chi2 * reduced_chi2_const, ")"
                endif
                mm = 0
                do i = 1, nvars
                    if (var_types(i) .eq. 2) then
                        mm = mm + 1
                        if (mm .eq. 1) write(*, *), 'Model counts:'
                        if (min_search(i) .ne. max_search(i)) then
                            write(*,'(2X,A15)',advance='no'), var_names(i)
                            do j = nint(min_search(i)) + 1, nint(max_search(i)) + 1
                                write(*, '(I6)', advance='no'), scount(mm,j-nint(min_search(i)))
                            enddo
                            write(*,*) ''
                        endif
                    endif
                enddo
            endif

            call MPI_BCAST(ybeststep, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(ybest, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

#ifdef CLEAR_BETWEEN_STEPS
            call system('clear')
#endif
            if (my_pe .eq. 0) then
                print_extra = .true.
                do e = 1, event_n
                    call set_event(e)
                    call print_trial_vars()
                enddo
                print_extra = .false.
                write (*,'(A)'), warnings
            endif
        enddo

        deallocate(praw)
        deallocate(chi2raw)
        deallocate(yraw)
        deallocate(psend)
        deallocate(chi2send)
        deallocate(ysend)

        if (calc_acor .or. (dump_all_walkers .and. nvars .ne. 0)) then
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering [5]'
#endif
            if (my_pe .eq. 0) then
                allocate(praw(n_pe*nmcsteps*nwalkers*nvars))
                allocate(chi2raw(n_pe*nmcsteps*nwalkers))
                allocate(yraw(n_pe*nmcsteps*nwalkers))
                allocate(pall(nmcsteps, nwalkers, nvars, n_pe))
                allocate(chi2all(nmcsteps, nwalkers, n_pe))
                allocate(yall(nmcsteps, nwalkers, n_pe))
            else
                allocate(praw(1))
                allocate(chi2raw(1))
                allocate(yraw(1))
            endif
            rcount = nmcsteps*nwalkers*nvars
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(reshape(p, [rcount]), rcount, MPI_DOUBLE_PRECISION, praw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            rcount = nmcsteps*nwalkers
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(reshape(chi2, [rcount]), rcount, MPI_DOUBLE_PRECISION, chi2raw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(reshape(y, [rcount]), rcount, MPI_DOUBLE_PRECISION, yraw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering finished [5]'
#endif
            if (my_pe .eq. 0) then
                pall = reshape(praw, [nmcsteps, nwalkers, nvars, n_pe])
                chi2all = reshape(chi2raw, [nmcsteps, nwalkers, n_pe])
                yall = reshape(yraw, [nmcsteps, nwalkers, n_pe])
            endif

            deallocate(praw)
            deallocate(chi2raw)
            deallocate(yraw)
        endif

        if (dump_burned_walkers) then
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering [5]'
#endif
            if (my_pe .eq. 0) then
                allocate(praw(n_pe*nburn*nwalkers*nvars))
                allocate(chi2raw(n_pe*nburn*nwalkers))
                allocate(yraw(n_pe*nburn*nwalkers))
                allocate(pburn(nburn, nwalkers, nvars, n_pe))
                allocate(chi2burn(nburn, nwalkers, n_pe))
                allocate(yburn(nburn, nwalkers, n_pe))
            else
                allocate(praw(1))
                allocate(chi2raw(1))
                allocate(yraw(1))
            endif

            allocate(psend(nburn*nwalkers*nvars))
            allocate(chi2send(nburn*nwalkers))
            allocate(ysend(nburn*nwalkers))

            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            rcount = nburn*nwalkers*nvars
            psend = reshape(p(burn_in:burn_in+nburn-1,:,:), [rcount])
            call MPI_GATHER(psend, rcount, MPI_DOUBLE_PRECISION, praw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            rcount = nburn*nwalkers
            chi2send = reshape(chi2(burn_in:burn_in+nburn-1,:), [rcount])
            ysend = reshape(y(burn_in:burn_in+nburn-1,:), [rcount])
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(chi2send, rcount, MPI_DOUBLE_PRECISION, chi2raw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHER(ysend, rcount, MPI_DOUBLE_PRECISION, yraw, &
                rcount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef VERBOSE_MPI
            if (my_pe .eq. 0) write(*, *), 'Gathering finished [5]'
#endif
            if (my_pe .eq. 0) then
                pburn = reshape(praw, [nburn, nwalkers, nvars, n_pe])
                chi2burn = reshape(chi2raw, [nburn, nwalkers, n_pe])
                yburn = reshape(yraw, [nburn, nwalkers, n_pe])
            endif

            deallocate(psend)
            deallocate(chi2send)
            deallocate(ysend)
        endif

        if (my_pe .eq. 0) then
            penalty = set_trial_vars(chi2best_p)
            write(*, *), 'Best fit found: '
            write(*, '(A7,G14.4E3)'),    'Chi2: ', chi2best * reduced_chi2_const
            write(*, '(A7,G14.4E3)'),    'y: ', chi2best_y

            call print_trial_vars()

            penalty = set_trial_vars(ybest_p)
            write(*, *), 'Most-probable fit found: '
            write(*, '(A7,G14.4E3)'),    'Chi2: ', ybest_chi2 * reduced_chi2_const
            write(*, '(A7,G14.4E3)'),    'y: ', ybest

            call print_trial_vars()

            write(*,*) ''
            write(*, '(G11.4,A)'), 100.d0*modelcnt, ' % complete'
            write(*,*) ''
        endif
    endif

    call tdefit_print('Writing to disk')
    
    if (my_pe .eq. 0) then
        do e = 1, event_n
            ! Write out the best fit, first half
            if (preferred_mode) then
                open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_preferred_best_fit' // '.dat', status = 'unknown')
            else
                penalty = set_trial_vars(chi2best_p)
                open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_best_fit' // '.dat', status = 'unknown')
            endif

            call set_event(e)

            trial_fbs(:cur_npts,cur_event) = 0.d0
            trial_mdots(:cur_npts,cur_event) = 0.d0
            trial_mags(:cur_npts,cur_event) = 0.d0
            trial_routs(:cur_npts,cur_event) = 0.d0
            trial_rphots(:cur_npts,cur_event) = 0.d0
            trial_times(:cur_npts,cur_event) = event_times(:cur_npts,cur_event) / trial_1pz(cur_event)
            call dmdt(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_fbs(:cur_npts,cur_event), &
                      .false., trial_menv(:cur_npts,cur_event))
            call dmdt(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_mdots(:cur_npts,cur_event), &
                      .true., trial_menv(:cur_npts,cur_event))
            call bandmag(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_fbs(:cur_npts,cur_event), &
                         trial_mdots(:cur_npts,cur_event), event_bands(:cur_npts,cur_event), &
                         trial_mags(:cur_npts,cur_event), event_penalties(:cur_npts,cur_event), trial_routs(:cur_npts,cur_event), trial_rphots(:cur_npts,cur_event))
            fat_best_fit(cur_event) = first_accretion_time
            if (output_restframe) then
                trial_mags(:cur_npts,cur_event) = trial_mags(:cur_npts,cur_event) - mag_fac*dlog10(trial_1pz(cur_event))
            endif
            
            call write_vars(fn)
            write(fn,'(A)') ''
            write(fn,'(2(G18.8,X))') chi2best * reduced_chi2_const, chi2best_y

            cur => ll(cur_event)%p
            do while (associated(cur))
                write(fn,'(A,1X)',advance='no'), cur%band
                cur => cur%next
            enddo
            write(fn,'(A)') ''

            write(str,'(I9)') cur_npts
            if (output_restframe) then
                write(fn,'('//trim(str)//'G17.8E3)') (event_times(i,cur_event)/trial_1pz(cur_event),i=1,cur_npts)
            else
                write(fn,'('//trim(str)//'G17.8E3)') (event_times(i,cur_event),i=1,cur_npts)
            endif
            write(fn,'('//trim(str)//'G17.8E3)') (trial_fbs(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mdots(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'(A,X))') (event_bands(i,cur_event),i=1,cur_npts)
            if (output_restframe) then
                write(fn,'('//trim(str)//'G17.8E3)') (event_ABs(i,cur_event) - mag_fac*dlog10(trial_1pz(cur_event)),&
                    i=1,cur_npts)
            else
                write(fn,'('//trim(str)//'G17.8E3)') (event_ABs(i,cur_event),i=1,cur_npts)
            endif
            write(fn,'('//trim(str)//'I2)') (event_types(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (dsqrt(event_errs(i,cur_event)),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mags(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_routs(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_rphots(i,cur_event),i=1,cur_npts)

            close(fn)

            ! Now write out the most probable fit, first half

            penalty = set_trial_vars(ybest_p)
            open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_most_probable' // '.dat', status = 'unknown')

            trial_fbs(:cur_npts,cur_event) = 0.d0
            trial_mdots(:cur_npts,cur_event) = 0.d0
            trial_mags(:cur_npts,cur_event) = 0.d0
            trial_routs(:cur_npts,cur_event) = 0.d0
            trial_rphots(:cur_npts,cur_event) = 0.d0
            trial_times(:cur_npts,cur_event) = event_times(:cur_npts,cur_event) / trial_1pz(cur_event)
            call dmdt(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_fbs(:cur_npts,cur_event), &
                      .false., trial_menv(:cur_npts,cur_event))
            call dmdt(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_mdots(:cur_npts,cur_event), &
                      .true., trial_menv(:cur_npts,cur_event))
            call bandmag(trial_times(:cur_npts,cur_event) + trial_toff(cur_event), trial_fbs(:cur_npts,cur_event), &
                         trial_mdots(:cur_npts,cur_event), event_bands(:cur_npts,cur_event), &
                         trial_mags(:cur_npts,cur_event), event_penalties(:cur_npts,cur_event), trial_routs(:cur_npts,cur_event), trial_rphots(:cur_npts,cur_event))
            fat_most_probable(cur_event) = first_accretion_time
            if (output_restframe) then
                trial_mags(:cur_npts,cur_event) = trial_mags(:cur_npts,cur_event) - mag_fac*dlog10(trial_1pz(cur_event))
            endif

            call write_vars(fn)
            write(fn,'(A)') ''
            write(fn,'(2(G18.8,1X))') ybest_chi2 * reduced_chi2_const, ybest

            cur => ll(cur_event)%p
            do while (associated(cur))
                write(fn,'(A,1X)',advance='no'), cur%band
                cur => cur%next
            enddo
            write(fn,'(A)') ''

            write(str,'(I9)') cur_npts
            if (output_restframe) then
                write(fn,'('//trim(str)//'G17.8E3)') (event_times(i,cur_event)/trial_1pz(cur_event),i=1,cur_npts)
            else
                write(fn,'('//trim(str)//'G17.8E3)') (event_times(i,cur_event),i=1,cur_npts)
            endif
            write(fn,'('//trim(str)//'G17.8E3)') (trial_fbs(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mdots(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'(A,X))') (event_bands(i,cur_event),i=1,cur_npts)
            if (output_restframe) then
                write(fn,'('//trim(str)//'G17.8E3)') (event_ABs(i,cur_event) - mag_fac*dlog10(trial_1pz(cur_event)),&
                    i=1,cur_npts)
            else
                write(fn,'('//trim(str)//'G17.8E3)') (event_ABs(i,cur_event),i=1,cur_npts)
            endif
            write(fn,'('//trim(str)//'I2)') (event_types(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (dsqrt(event_errs(i,cur_event)),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mags(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_routs(i,cur_event),i=1,cur_npts)
            write(fn,'('//trim(str)//'G17.8E3)') (trial_rphots(i,cur_event),i=1,cur_npts)

            close(fn)
        enddo

        ! Allocate combined walker array
        if (calc_acor .or. (dump_all_walkers .and. nvars .ne. 0)) then
            print *, 'Allocating walker arrays...'
            print *, nmcsteps*nwalkers*n_pe*(nvars+2)

            allocate(walkerout(nmcsteps,nwalkers*n_pe,nvars+2))

            print *, 'Filling walker arrays...'
            do j = 1, nwalkers
                do k = 1, n_pe
                    walkerout(:,(k-1)*nwalkers+j,1) = chi2all(:,j,k)
                    walkerout(:,(k-1)*nwalkers+j,2) = yall(:,j,k)
                    walkerout(:,(k-1)*nwalkers+j,3:) = pall(:,j,:,k)
                enddo
            enddo

            print *, 'Calculating auto-correlation time...'
            call acor(walkerout(nanneal:,:,3:), atime)
            print *, 'atimes', atime
        endif

        ! Write out the walker positions.
        if (dump_all_walkers .and. nvars .ne. 0) then
            write(*, *), 'Dumping all walkers to file...'
            open(unit = fn, file = trim(output_path) // 'walkers' // '.dat', status = 'unknown', form = 'unformatted')
            write(fn) nmcsteps, n_pe*nwalkers, nvars
            do i = 1, nvars
                write(fn) var_names(i), min_search(i), max_search(i), var_types(i)
            enddo
            write(fn) (((walkerout(i,j,k),k=1,nvars+2),j=1,nwalkers*n_pe),i=1,nmcsteps)
            close(fn)
        endif

        if (calc_acor .or. (dump_all_walkers .and. nvars .ne. 0)) then
            deallocate(walkerout)
        endif

        if (dump_burned_walkers .or. dump_burned_ensembles) then
            allocate(walkerout(nburn,nwalkers*n_pe,nvars+2))
            do i = 1, nburn
                do j = 1, nwalkers
                    do k = 1, n_pe
                        walkerout(i,(k-1)*nwalkers+j,1) = chi2burn(i,j,k) * &
                            reduced_chi2_const
                        walkerout(i,(k-1)*nwalkers+j,2) = yburn(i,j,k)
                        do l = 1, nvars
                            walkerout(i,(k-1)*nwalkers+j,l+2) = pburn(i,j,l,k)
                        enddo
                    enddo
                enddo
            enddo
        endif

        if (dump_burned_walkers) then
            write(*, *), 'Dumping burned walkers to file...'
            open(unit = fn, file = trim(output_path) // 'burned_walkers' // '.dat', status = 'unknown', form = 'unformatted')
            write(fn) nburn, n_pe*nwalkers, nvars
            do i = 1, nvars
                write(fn) var_names(i), min_search(i), max_search(i), var_types(i)
            enddo
            write(fn) (((walkerout(i,j,k),k=1,nvars+2),j=1,nwalkers*n_pe),i=1,nburn)
            close(fn)

            if (.not. dump_burned_ensembles) then
                deallocate(walkerout)
            endif
        endif

        deallocate(event_bands)
        deallocate(event_times)
        deallocate(event_penalties)
        deallocate(trial_fbs)
        deallocate(trial_mdots)
        deallocate(trial_menv)
        deallocate(trial_routs)
        deallocate(trial_rphots)
        deallocate(trial_mags)
        deallocate(trial_times)

        ! Always produce SEDs for the smooth fits
        produce_seds = .true.

        ! Now write out a smooth curve for the best fit
        penalty = set_trial_vars(chi2best_p)

        write(*, *), 'Producing smooth light curves for best fit...'

        do e = 1, event_n
            call set_event(e)

            ! Second index is e:e because functions assume first dimension refers
            ! to event ID.
            allocate(event_bands(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(event_times(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(event_penalties(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_fbs(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_mdots(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_menv(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_routs(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_rphots(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_mags(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_times(nbest_times*event_nbest_bands(cur_event),e:e))

            open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_best_fit' // '.dat', status = 'unknown', position = 'append')

            ! For constructing nicely sampled light curves when we're done.
            lmin_trial_time = dlog10(max(min_best_time, fat_best_fit(cur_event))) 
            lmax_trial_time = dlog10(max_best_time)
            do i = 1, nbest_times*event_nbest_bands(cur_event)
                event_bands(i,cur_event) = event_best_bands(mod(i - 1, event_nbest_bands(cur_event)) + 1,cur_event)
                if (mod(i - 1, event_nbest_bands(cur_event)) .ne. 0) then
                    trial_times(i,cur_event) = trial_times(i-1,cur_event)
                else
                    if (i .eq. 1) then
                        trial_times(i,cur_event) = lmin_trial_time
                    else
                        trial_times(i,cur_event) = trial_times(i-1,cur_event) + &
                            (lmax_trial_time - lmin_trial_time) / (nbest_times - 1)
                    endif
                endif
            enddo
            trial_times(:,cur_event) = 10.d0**trial_times(:,cur_event)
            event_times(:,cur_event) = (trial_times(:,cur_event) - trial_toff(cur_event))
            if (.not. output_restframe) event_times(:,cur_event) = event_times(:,cur_event)*trial_1pz(cur_event)

            trial_fbs(:,cur_event) = 0.d0
            trial_mdots(:,cur_event) = 0.d0
            trial_routs(:,cur_event) = 0.d0
            trial_rphots(:,cur_event) = 0.d0
            trial_mags(:,cur_event) = 0.d0
            call dmdt(trial_times(:,cur_event), trial_fbs(:,cur_event), .false., trial_menv(:,cur_event))
            call dmdt(trial_times(:,cur_event), trial_mdots(:,cur_event), .true., trial_menv(:,cur_event))
            call bandmag(trial_times(:,cur_event), trial_fbs(:,cur_event), trial_mdots(:,cur_event), &
                         event_bands(:,cur_event), trial_mags(:,cur_event), event_penalties(:,cur_event), &
                         trial_routs(:,cur_event), trial_rphots(:,cur_event))
            if (output_restframe) then
                trial_mags(:,cur_event) = trial_mags(:,cur_event) - mag_fac*dlog10(trial_1pz(cur_event))
            endif

            write(str,'(I9)') nbest_times*event_nbest_bands(cur_event)
            write(fn,'('//trim(str)//'(A,X))') (event_bands(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (event_times(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_fbs(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mdots(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mags(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_routs(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_rphots(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))

            close(fn)

            write(*, *), 'Writing best fit SEDs'
            open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_best_fit_seds' // '.dat', status = 'unknown', form = 'unformatted')
            write(fn) nbest_times, sed_nsteps, sed_min, sed_max, sed_step
            write(fn) (((sed_table(i,j,k),k=1,sed_nsteps),j=1,nbest_times),i=1,2)
            close(fn)

            deallocate(event_bands)
            deallocate(event_times)
            deallocate(event_penalties)
            deallocate(trial_fbs)
            deallocate(trial_mdots)
            deallocate(trial_menv)
            deallocate(trial_routs)
            deallocate(trial_rphots)
            deallocate(trial_mags)
            deallocate(trial_times)
        enddo

        ! Now write out a smooth curve for the most-probable fit
        penalty = set_trial_vars(ybest_p)

        write(*, *), 'Producing smooth light curves for most probable fit...'

        ! For constructing nicely sampled light curves when we're done.
        do e = 1, event_n
            call set_event(e)

            allocate(event_bands(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(event_times(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(event_penalties(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_fbs(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_mdots(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_menv(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_routs(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_rphots(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_mags(nbest_times*event_nbest_bands(cur_event),e:e))
            allocate(trial_times(nbest_times*event_nbest_bands(cur_event),e:e))

            open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_most_probable' // '.dat', status = 'unknown', position = 'append')

            lmin_trial_time = dlog10(max(min_best_time, fat_most_probable(cur_event))) 
            lmax_trial_time = dlog10(max_best_time)
            do i = 1, nbest_times*event_nbest_bands(cur_event)
                event_bands(i,cur_event) = event_best_bands(mod(i - 1, event_nbest_bands(cur_event)) + 1,cur_event)
                if (mod(i - 1, event_nbest_bands(cur_event)) .ne. 0) then
                    trial_times(i,cur_event) = trial_times(i-1,cur_event)
                else
                    if (i .eq. 1) then
                        trial_times(i,cur_event) = lmin_trial_time
                    else
                        trial_times(i,cur_event) = trial_times(i-1,cur_event) + &
                            (lmax_trial_time - lmin_trial_time) / (nbest_times - 1)
                    endif
                endif
            enddo
            trial_times(:,cur_event) = 10.d0**trial_times(:,cur_event)
            event_times(:,cur_event) = trial_times(:,cur_event) - trial_toff(cur_event)
            if (.not. output_restframe) event_times(:,cur_event) = event_times(:,cur_event)*trial_1pz(cur_event)

            trial_fbs(:,cur_event) = 0.d0
            trial_mdots(:,cur_event) = 0.d0
            trial_routs(:,cur_event) = 0.d0
            trial_rphots(:,cur_event) = 0.d0
            trial_mags(:,cur_event) = 0.d0
            call dmdt(trial_times(:,cur_event), trial_fbs(:,cur_event), .false., trial_menv(:,cur_event))
            call dmdt(trial_times(:,cur_event), trial_mdots(:,cur_event), .true., trial_menv(:,cur_event))
            call bandmag(trial_times(:,cur_event), trial_fbs(:,cur_event), trial_mdots(:,cur_event), &
                         event_bands(:,cur_event), trial_mags(:,cur_event), &
                         event_penalties(:,cur_event), trial_routs(:,cur_event), trial_rphots(:,cur_event))
            if (output_restframe) then
                trial_mags(:,cur_event) = trial_mags(:,cur_event) - mag_fac*dlog10(trial_1pz(cur_event))
            endif

            write(str,'(I9)') nbest_times*event_nbest_bands(cur_event)
            write(fn,'('//trim(str)//'(A,X))') (event_bands(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (event_times(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_fbs(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mdots(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_mags(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_routs(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))
            write(fn,'('//trim(str)//'G17.8E3)') (trial_rphots(i,cur_event),i=1,nbest_times*event_nbest_bands(cur_event))

            close(fn)

            write(*, *), 'Writing most probable SEDs'
            open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_most_probable_seds' // '.dat', status = 'unknown', form = 'unformatted')
            write(fn) nbest_times, sed_nsteps, sed_min, sed_max, sed_step
            write(fn) (((sed_table(i,j,k),k=1,sed_nsteps),j=1,nbest_times),i=1,2)
            close(fn)

            deallocate(event_bands)
            deallocate(event_times)
            deallocate(event_penalties)
            deallocate(trial_fbs)
            deallocate(trial_mdots)
            deallocate(trial_menv)
            deallocate(trial_routs)
            deallocate(trial_rphots)
            deallocate(trial_mags)
            deallocate(trial_times)
        enddo

        ! Now write out an ensemble of walker fits
        if (dump_burned_ensembles) then
            write (*, *), 'Writing ensemble of walker fits.'

            if (initial_walker_dist .eq. IW_REGULAR) then
                call tdefit_print('When initial_walker_dist is set to IW_REGULAR, ensemble includes all parameter combos.')
                k = 0
                l = nwalkers*n_pe
                do j = 1, l
                    penalty = set_trial_vars(walkerout(nburn,j,3:))
                    if (penalty .eq. huge(1.d0)) cycle
                    k = k + 1
                enddo
            else
                k = nensemble
                l = k
            endif

            do e = 1, event_n
                call set_event(e)

                allocate(event_bands(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(event_times(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(event_penalties(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_fbs(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_mdots(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_menv(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_routs(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_rphots(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_mags(nbest_times*event_nbest_bands(cur_event),e:e))
                allocate(trial_times(nbest_times*event_nbest_bands(cur_event),e:e))

                open(unit = fn, file = trim(output_path) // trim(event_fnames(e)) // '_ensemble_fits' // '.dat', status = 'unknown', form='unformatted')

                write(fn), k, nbest_times

                lmin_trial_time = dlog10(max(min_best_time, fat_most_probable(cur_event))) 
                lmax_trial_time = dlog10(max_best_time)
                do j = 1, l
                    do i = 1, nbest_times*event_nbest_bands(cur_event)
                        event_bands(i,cur_event) = event_best_bands(mod(i - 1, event_nbest_bands(cur_event)) + 1,cur_event)
                        if (mod(i - 1, event_nbest_bands(cur_event)) .ne. 0) then
                            trial_times(i,cur_event) = trial_times(i-1,cur_event)
                        else
                            if (i .eq. 1) then
                                trial_times(i,cur_event) = lmin_trial_time
                            else
                                trial_times(i,cur_event) = trial_times(i-1,cur_event) + &
                                    (lmax_trial_time - lmin_trial_time) / (nbest_times - 1)
                            endif
                        endif
                    enddo
                    print *, 'Outputting realization #', j
                    if (initial_walker_dist .eq. IW_REGULAR) then
                        ranj = j
                    else
                        call random_number(dummy)
                        ranj = ceiling(dummy*nwalkers*n_pe)
                    endif

                    penalty = set_trial_vars(walkerout(nburn,ranj,3:))

                    trial_times(:,cur_event) = 10.d0**trial_times(:,cur_event)
                    event_times(:,cur_event) = trial_times(:,cur_event) - trial_toff(cur_event)
                    if (.not. output_restframe) event_times(:,cur_event) = event_times(:,cur_event)*trial_1pz(cur_event)

                    call print_trial_vars()
                    if (penalty .eq. huge(1.d0)) cycle

                    trial_fbs(:,cur_event) = 0.d0
                    trial_mdots(:,cur_event) = 0.d0
                    trial_mags(:,cur_event) = 0.d0
                    trial_routs(:,cur_event) = 0.d0
                    trial_rphots(:,cur_event) = 0.d0
                    call dmdt(trial_times(:,cur_event), trial_fbs(:,cur_event), .false., trial_menv(:,cur_event))
                    call dmdt(trial_times(:,cur_event), trial_mdots(:,cur_event), .true., trial_menv(:,cur_event))
                    call bandmag(trial_times(:,cur_event), trial_fbs(:,cur_event), trial_mdots(:,cur_event), &
                                 event_bands(:,cur_event), trial_mags(:,cur_event), &
                                 event_penalties(:,cur_event), trial_routs(:,cur_event), trial_rphots(:,cur_event))
                    if (output_restframe) then
                        trial_mags(:,cur_event) = trial_mags(:,cur_event) - mag_fac*dlog10(trial_1pz(cur_event))
                    endif

                    do i = 1, size(all_var_names)
                        var_dbles(i) = get_var(all_var_names(i))
                    enddo
                    write(fn) var_dbles(:)
                    write(fn) event_times(:,cur_event)
                    write(fn) trial_fbs(:,cur_event)
                    write(fn) trial_mdots(:,cur_event)
                    write(fn) trial_mags(:,cur_event)
                    write(fn) trial_routs(:,cur_event)
                    write(fn) trial_rphots(:,cur_event)
                enddo

                close(fn)

                deallocate(event_bands)
                deallocate(event_times)
                deallocate(event_penalties)
                deallocate(trial_fbs)
                deallocate(trial_mdots)
                deallocate(trial_menv)
                deallocate(trial_routs)
                deallocate(trial_rphots)
                deallocate(trial_mags)
                deallocate(trial_times)
            enddo

            deallocate(walkerout)
        endif
        
    endif

    ! Clear band list
    do e = 1, event_n
        call set_event(e)

        cur => ll(cur_event)%p
        do while (associated(cur))
            ll(cur_event)%p => cur%next
            deallocate(cur)
            cur => ll(cur_event)%p
        enddo
    enddo

    deallocate(event_fnames)
    deallocate(event_npts)
    deallocate(event_blrpts)
    deallocate(event_nh)
    deallocate(event_claimed_z)
    deallocate(event_min_aspin)
    deallocate(event_nhcorr)
    deallocate(event_restframe)
    deallocate(event_time_units)
    deallocate(event_blr_time_units)
    deallocate(event_blr_times)
    deallocate(event_blr_vels)
    deallocate(event_blr_bands)
    deallocate(event_blr_exists)
    deallocate(trial_mh)
    deallocate(trial_ms)
    deallocate(trial_rs)
    deallocate(trial_beta)
    deallocate(trial_aspin)
    deallocate(trial_nhsrc)
    deallocate(trial_phi)
    deallocate(trial_bh_rms)
    deallocate(trial_rp)
    deallocate(trial_toff)
    deallocate(trial_rout)
    deallocate(trial_rg)
    deallocate(trial_gmh)
    deallocate(trial_ms0)
    deallocate(trial_mh0)
    deallocate(trial_rs0)
    deallocate(trial_rsc)
    deallocate(trial_alphhr)
    deallocate(trial_eps_edd)
    deallocate(trial_source_rv)
    deallocate(trial_magoff)
    deallocate(trial_fcor)
    deallocate(trial_tlimit)
    deallocate(trial_ecor)
    deallocate(trial_r_isco)
    deallocate(trial_reprocess_temp)
    deallocate(trial_offset_X1)
    deallocate(trial_offset_X2)
    deallocate(trial_offset_GN)
    deallocate(trial_offset_Pg)
    deallocate(trial_offset_Pr)
    deallocate(trial_offset_Pi)
    deallocate(trial_offset_Pz)
    deallocate(trial_offset_U1)
    deallocate(trial_offset_U2)
    deallocate(trial_offset_RO)
    deallocate(trial_offset_Ub)
    deallocate(trial_offset_Um)
    deallocate(trial_offset_Uu)
    deallocate(trial_offset_Uv)
    deallocate(trial_offset_bV)
    deallocate(trial_fout)
    deallocate(trial_mu_e)
    deallocate(trial_ledd)
    deallocate(trial_outflow_frac)
    deallocate(trial_bh_rbs)
    deallocate(trial_r_ibco)
    deallocate(trial_temp_mult)
    deallocate(trial_z)
    deallocate(trial_1pz)
    deallocate(trial_dl)
    deallocate(trial_exp_1)
    deallocate(trial_exp_2)
    deallocate(trial_exp_3)
    deallocate(trial_exp_4)
    deallocate(trial_variability)
    deallocate(trial_variance)
    deallocate(trial_viscous_time)
    deallocate(trial_variability2)
    deallocate(trial_variance2)
    deallocate(trial_yms)
    deallocate(trial_y1)
    deallocate(trial_y2)
    deallocate(trial_y3)
    deallocate(trial_opacity)
    deallocate(trial_rphot)
    deallocate(trial_model)
    deallocate(trial_outflow_model)
    deallocate(trial_object_type)
    deallocate(trial_temperature_model)
    deallocate(trial_blr_model)
    deallocate(trial_time_dep_rin)
    deallocate(trial_time_dep_rout)
    deallocate(trial_simple_bb)
    deallocate(trial_use_fcor)
    deallocate(trial_cap_at_edd)
    deallocate(trial_full_disk_coverage)
    deallocate(event_best_bands)
    deallocate(event_ABs)
    deallocate(event_errs)
    deallocate(event_devs)
    deallocate(event_weights)
    deallocate(event_types)
    deallocate(ll)

    deallocate(filt_names)
    deallocate(filt_files)
    deallocate(filt_len)
    deallocate(filt_resp)
    deallocate(dfilt_resp)
    deallocate(filt_lmin)
    deallocate(filt_lmax)

    call MPI_FINALIZE(ret_c)

    if (my_pe .eq. 0) write(*, *), 'Done!'
end program
