subroutine init
#include "tdefit.fpp"
#include "version.fpp"
    use constants
    use tdefit_interface, only: filterfunc, filtnormfunc, interp_flash_output, &
                                load_defaults, load_user_vars
    use tdefit_data
#ifndef MPIH
    use mpi
#endif
    use adapt_quad, only: qxgs
    use tdefit_util, only: least_sq


#ifdef MPIH
#include "mpif.h"
#endif

    character :: dummy
    character*200 :: dir_path, bin_path, str

    integer :: i, j, k, s, ierr, mark, imark, fn, fn2, cur_model, i_zero, &
               i_roche_begin, i_roche_end, i_early_begin, i_early_end
    real :: kersum, beta_destroy, early_e
    logical :: exists
    integer, dimension(13) :: statb1, statb2

    real :: epscor, epsscl, roche_ener, maxdmde, log_offset, min_neg_e, max_pos_e
    real, dimension(:,:), allocatable :: temp_ddat, ddat_temp

    real :: err, dx, fita, fitb, siga, sigb, fitchi2, fitq, denom, window, arg
    integer :: neval, fcnt, maxdmdeloc, maxtauloc

    fn = 11
    fn2 = 12

    !First read the paths file to set file read/write locations.
    open(unit = fn, file = path_file, status='old', action='read')
    read(fn,"(A)") exec_path 
    read(fn,"(A)") output_path 
    read(fn,"(A)") binary_path
    read(fn,"(A)") input_path 
    read(fn,"(A)") event_path 
    close(fn)

    call load_defaults(0)

    call load_user_vars(0)

    call check_options

    if (burn_in .le. 0) then
        burn_in = nmcsteps + burn_in
    endif

    nburn = nmcsteps - burn_in + 1

    open(unit = fn, file = trim(input_path) // "run_info.dat", status='old', action='read')
    read(fn, *) nmodels
    allocate(model_index(nmodels))
    allocate(model_beta_destroy(nmodels))
    do i = 1, nmodels
        read(fn, *) model_index(i), model_beta_destroy(i)
    enddo
    nruns = sum(model_index)
    do i = nmodels, 2, -1
        model_index(i) = model_index(i-1)
    enddo
    model_index(1) = 1
    do i = 2, nmodels
        model_index(i) = model_index(i-1) + model_index(i)
    enddo

    allocate(dir_names(1:nruns))
    allocate(sim_dele(1:nruns))
    allocate(sim_mfinal(1:nruns))
    do i = 1, nruns
        read(fn,'(A)') dir_names(i)
    enddo
    close(fn)

    allocate(mdat(1:nruns, 1:maxnrows, 1:MDAT_NCOLS))
    allocate(odat(1:nruns, 1:maxnrows, 1:ODAT_NCOLS))
    allocate(ddat(1:nruns, 1:2, 1:maxnrows))
    allocate(edat(1:nruns, 1:NUM_EDAT))
    allocate(mdat_row(1:MDAT_NCOLS))
    allocate(odat_row(1:ODAT_NCOLS))
    allocate(mdat_nrows(1:nruns))
    allocate(odat_nrows(1:nruns))
    allocate(ddat_ncols(1:nruns))
    allocate(ddat_bnd(1:nruns))
    allocate(ddat_time(1:nruns))
    allocate(d_emin(1:nruns))
    allocate(d_emid(1:nruns))
    allocate(d_emax(1:nruns))

    sim_dele = 0.d0
    mdat_nrows = 0
    odat_nrows = 0
    ddat_ncols = 0
    ddat_bnd = 0.d0
    ddat_time = 0.d0

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call tdefit_print('Reading simulation data files...')
    do i = 1, nruns
#ifdef VERBOSE_INIT
        call tdefit_print('Reading ' // trim(dir_names(i)))
#endif

        write(str,'(I9)') i
        str = adjustl(str)
        dir_path = trim(input_path) // trim(dir_names(i))
        bin_path = trim(binary_path) // 'extras_' // trim(str) // '.bin'

#ifdef VERBOSE_INIT
        call tdefit_print('Reading extras')
#endif
        inquire(file = bin_path, exist = exists)
        !if (exists) then
        !    call stat(trim(dir_path) // '/extras.bin', statb1, ierr)
        !    call stat(trim(dir_path) // '/extras.dat', statb2, ierr)
        !    if (statb1(fn) .lt. statb2(fn)) then
        !        print *, 'dat file newer, regenerating bin'
        !        exists = .false.
        !    endif
        !endif
        if (force_reload_data .or. .not. exists) then
            if (my_pe .eq. 0) then
                call tdefit_print('Binary file not found, converting...')
                open(unit = fn, file = trim(dir_path) // '/extras.dat', status='old', action='read')
                do j = 1, NUM_EDAT
                    read(fn, *) edat(i,j)
                enddo
                close(fn)

                open(unit = fn, file = bin_path, status = 'unknown', form = 'unformatted')
                write(fn) (edat(i,j),j=1,NUM_EDAT)
                call flush(fn)
                close(fn)
            endif
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do while(.not. exists)
            call sleep(1)
            inquire(file = bin_path, exist = exists)
        enddo

        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) (edat(i,j),j=1,NUM_EDAT)
        close(fn)

#ifdef VERBOSE_INIT
        call tdefit_print('Reading multipoly')
#endif
        bin_path = trim(binary_path) // 'multipoly_' // trim(str) // '.bin'
        inquire(file = bin_path, exist = exists)
        !if (exists) then
        !    call stat(trim(dir_path) // '/multipoly.bin', statb1, ierr)
        !    call stat(trim(dir_path) // '/multipoly.dat', statb2, ierr)
        !    if (statb1(fn) .lt. statb2(fn)) then
        !        print *, 'dat file newer, regenerating bin'
        !        exists = .false.
        !    endif
        !endif
        if (force_reload_data .or. .not. exists) then
            if (my_pe .eq. 0) then
                call tdefit_print('Binary file not found, converting...')
                open(unit = fn, file = trim(dir_path) // '/multipoly.dat', status='old', action='read')
                read(fn, fmt='(A1)') dummy !Skip header line
                j = 1
                do
                    read(fn, *, iostat=ierr) mdat_row
                    if (ierr /= 0) then
                        mdat_nrows(i) = j - 1

                        exit
                    endif
                    mdat(i,j,:) = mdat_row
                    j = j + 1
                enddo
                close(fn)

                open(unit = fn, file = bin_path, status = 'unknown', form = 'unformatted')
                write(fn) mdat_nrows(i)
                do j = 1, mdat_nrows(i)
                    write(fn) (mdat(i,j,k),k=1,MDAT_NCOLS)
                enddo
                call flush(fn)
                close(fn)
            endif
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do while(.not. exists)
            call sleep(1)
            inquire(file = bin_path, exist = exists)
        enddo

        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) mdat_nrows(i)
        close(fn)

#ifdef VERBOSE_INIT
        call tdefit_print('Reading orbit')
#endif
        bin_path = trim(binary_path) // 'orbit_' // trim(str) // '.bin'
        inquire(file = bin_path, exist = exists)
        !if (exists) then
        !    call stat(trim(dir_path) // '/orbit.bin', statb1, ierr)
        !    call stat(trim(dir_path) // '/orbit.dat', statb2, ierr)
        !    if (statb1(fn) .lt. statb2(fn)) then
        !        print *, 'dat file newer, regenerating bin'
        !        exists = .false.
        !    endif
        !endif
        if (force_reload_data .or. .not. exists) then
            if (my_pe .eq. 0) then
                call tdefit_print('Creating binary file')
                open(unit = fn, file = trim(dir_path) // '/orbit.dat', status='old', action='read')
                read(fn, fmt='(A1)') dummy !Skip header line
                j = 1
                do
                    read(fn, *, iostat=ierr) odat_row
                    if (ierr /= 0) then
                        odat_nrows(i) = j - 1
                        exit
                    endif
                    odat(i,j,:) = odat_row
                    j = j + 1
                enddo
                close(fn)

                open(unit = fn, file = bin_path, status = 'unknown', form = 'unformatted')
                write(fn) odat_nrows(i)
                do j = 1, odat_nrows(i)
                    write(fn) (odat(i,j,k),k=1,ODAT_NCOLS)
                enddo
                call flush(fn)
                close(fn)
            endif
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do while(.not. exists)
            call sleep(1)
            inquire(file = bin_path, exist = exists)
        enddo

        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) odat_nrows(i)
        close(fn)

        if (mdat_nrows(i) .ne. odat_nrows(i)) then
            print *, mdat_nrows(i), odat_nrows(i)
            call tdefit_print('ERROR: Number of rows in multipoly.dat and orbit.dat do not match!')
            call exit(0)
        endif

#ifdef VERBOSE_INIT
        call tdefit_print('Reading dmde')
#endif
        bin_path = trim(binary_path) // 'dmde_' // trim(str) // '.bin'
        inquire(file = bin_path, exist = exists)
        !if (exists) then
        !    call stat(trim(dir_path) // '/dmde.bin', statb1, ierr)
        !    call stat(trim(dir_path) // '/dmde.dat', statb2, ierr)
        !    if (statb1(fn) .lt. statb2(fn)) then
        !        print *, 'dat file newer, regenerating bin'
        !        exists = .false.
        !    endif
        !endif
        if (force_reload_data .or. .not. exists) then
            if (my_pe .eq. 0) then
                call tdefit_print('Binary file not found, converting...')
                open(unit = fn, file = trim(dir_path) // '/dmde.dat', status='old', action='read')
                read(fn, *) ddat_ncols(i)
                read(fn, *) ddat_bnd(i)
                read(fn, *) ddat_time(i)
                do j = 1, 2
                    read(fn, *) ddat(i,j,1:ddat_ncols(i))
                enddo
                close(fn)

                open(unit = fn, file = bin_path, status = 'unknown', form = 'unformatted')
                write(fn) ddat_ncols(i)
                write(fn) ddat_bnd(i)
                write(fn) ddat_time(i)
                do j = 1, 2
                    write(fn) (ddat(i,j,k),k=1,ddat_ncols(i))
                enddo
                call flush(fn)
                close(fn)
            endif
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do while(.not. exists)
            call sleep(1)
            inquire(file = bin_path, exist = exists)
        enddo

        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) ddat_ncols(i)
        close(fn)
    enddo

    nrows = maxval(mdat_nrows)
    deallocate(mdat, odat, ddat)
    allocate(mdat(1:nruns, 1:nrows, 1:MDAT_NCOLS))
    allocate(odat(1:nruns, 1:nrows, 1:ODAT_NCOLS))
    allocate(ddat(1:nruns, 1:2, 1:(maxval(ddat_ncols)+n_early_bins)))
    allocate(iddat(1:nruns, 1:2, 1:(maxval(ddat_ncols)+n_early_bins)))
    allocate(rhoddat(1:nruns, 1:2, 1:(maxval(ddat_ncols)+n_early_bins)))
    allocate(ddat_coeff(1:nruns, 1:(maxval(ddat_ncols)+n_early_bins)))
    allocate(iddat_coeff(1:nruns, 1:(maxval(ddat_ncols)+n_early_bins)))
    allocate(rhoddat_coeff(1:nruns, 1:(maxval(ddat_ncols)+n_early_bins)))

    do i = 1, nruns
        dir_path = trim(input_path) // trim(dir_names(i))

        write(str,'(I9)') i
        str = adjustl(str)
        bin_path = trim(binary_path) // 'multipoly_' // trim(str) // '.bin'
        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) mdat_nrows(i)
        do j = 1, mdat_nrows(i)
            read(fn) (mdat(i,j,k),k=1,MDAT_NCOLS)
        enddo
        close(fn)

        bin_path = trim(binary_path) // 'orbit_' // trim(str) // '.bin'
        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) odat_nrows(i)
        do j = 1, odat_nrows(i)
            read(fn) (odat(i,j,k),k=1,ODAT_NCOLS)
        enddo
        close(fn)

        bin_path = trim(binary_path) // 'dmde_' // trim(str) // '.bin'
        open(unit = fn, file = bin_path, form = 'unformatted', status='old', action='read')
        read(fn) ddat_ncols(i)
        read(fn) ddat_bnd(i)
        read(fn) ddat_time(i)
        do j = 1, 2
            read(fn) (ddat(i,j,k),k=1,ddat_ncols(i))
        enddo
        close(fn)
    enddo

    !allocate(mu_tot(nruns, nrows))
    !allocate(mu_bnd(nruns, nrows))
    allocate(vec1(nruns, nrows, 6))
    allocate(vec2(nruns, nrows, 6))
    allocate(bndcom(nruns, nrows, 6))
    !allocate(mpolecom(nruns, nrows, 6))
    allocate(totcom(nruns, nrows, 6))
    !allocate(barycenter(nruns, nrows, 6))
    !allocate(v1(nruns, nrows))
    !allocate(v2(nruns, nrows))
    allocate(rad(nruns, nrows))
    allocate(vel(nruns, nrows))
    allocate(semimaj(nruns, nrows))
    allocate(maxdmdttime(nruns))
    allocate(maxtautime(nruns))

    do i = 1, nrows
        !mu_tot(:,i) = mdat(:,i,MTOT_VAR) / (mdat(:,i,MTOT_VAR) + edat(:,E_MPERT))
        !mu_bnd(:,i) = odat(:,i,MBND_VAR) / (odat(:,i,MBND_VAR) + edat(:,E_MPERT))
        vec1(:,i,:) = odat(:,i,VEC1X_VAR:VEC1VZ_VAR)
        vec2(:,i,:) = odat(:,i,VEC2X_VAR:VEC2VZ_VAR)
        bndcom(:,i,:) = odat(:,i,BNDCOMX_VAR:BNDCOMVZ_VAR)
        !mpolecom(:,i,:) = odat(:,i,MPOLECOMX_VAR:MPOLECOMVZ_VAR)
        totcom(:,i,:) = odat(:,i,TOTCOMX_VAR:TOTCOMVZ_VAR)
        do j = 1, nruns
            !barycenter(j,i,:) = (1.d0 - mu_tot(j,i))*vec1(j,i,:) + mu_tot(j,i)*vec2(j,i,:) + totcom(j,i,:)
            !bndcom(j,i,:) = vec2(j,i,:) + bndcom(j,i,:) - barycenter(j,i,:)
            !mpolecom(j,i,:) = vec2(j,i,:) + mpolecom(j,i,:) - barycenter(j,i,:)
            !totcom(j,i,:) = totcom(j,i,:) - barycenter(j,i,:)
            !vec1(j,i,:) = vec1(j,i,:) - barycenter(j,i,:)
            !vec2(j,i,:) = vec2(j,i,:) - barycenter(j,i,:)
            rad(j,i) = sqrt(sum((vec2(j,i,1:3) + bndcom(j,i,1:3) - vec1(j,i,1:3) - totcom(j,i,1:3))**2))
            vel(j,i) = sqrt(sum((vec2(j,i,4:6) + bndcom(j,i,4:6) - vec1(j,i,4:6) - totcom(j,i,4:6))**2))
            !v1(j,i) = sqrt(sum(vec1(j,i,4:6)**2))
            !v2(j,i) = sqrt(sum(vec2(j,i,4:6)**2))
            ! semimaj should technically be M_PERT + M_BOUND, but M_BOUND << M_PERT and thus can be ignored.
            semimaj(j,i) = G*(edat(j,E_MPERT) + odat(j,i,MBND_VAR))*rad(j,i)/&
                (2.d0*G*(edat(j,E_MPERT) + odat(j,i,MBND_VAR)) - rad(j,i)*vel(j,i)**2)
        enddo
    enddo

    cur_model = 0
    do i = 1, nruns
        if (cur_model .ne. nmodels) then
            if (i .eq. model_index(cur_model+1)) cur_model = cur_model + 1
        endif
        beta_destroy = model_beta_destroy(cur_model)

        if (edat(i,E_BETA) .le. beta_destroy) then
            imark = -1
            !print *, maxval(1.d0 - mdat(i,1:odat_nrows(i),MTOT_VAR)/mdat(i,1,MTOT_VAR))
            do j = 1, odat_nrows(i)
                if (odat(i,j,TIME_VAR) .ge. dele_time) then
                    imark = j
                    exit
                endif
                if (1.d0 - mdat(i,j,MTOT_VAR)/mdat(i,1,MTOT_VAR) .gt. dele_dmtot) then
                    imark = j
                    exit
                endif
            enddo
            if (imark .eq. -1) then
                call tdefit_print('dele_time and dele_dmtot criteria never met!')
                imark = odat_nrows(i)
                print *, 1.d0 - mdat(i,imark,MTOT_VAR)/mdat(i,1,MTOT_VAR)
                print *, odat(i,odat_nrows(i),TIME_VAR)
                call exit(0)
            endif
            sim_dele(i) = G*edat(i,E_MPERT)/2.d0*(1.d0/semimaj(i,1) - 1.d0/semimaj(i,imark))*msun/odat(i,1,MBND_VAR)
        endif

        if (edat(i,E_BETA) .gt. beta_destroy) then
            sim_mfinal(i) = 0.d0
        else
            sim_mfinal(i) = min(odat(i,odat_nrows(i),MBND_VAR),odat(i,1,MBND_VAR))*msun/odat(i,1,MBND_VAR)
        endif
        ddat(i,1,1:ddat_ncols(i)) = ddat(i,1,ddat_ncols(i):1:-1)
        ddat(i,2,1:ddat_ncols(i)) = ddat(i,2,ddat_ncols(i):1:-1)

        !Adjust data based on desired mbh, etc.
        epscor = -G*edat(i,E_MPERT)/(2.d0*edat(i,E_ROBJ)*(edat(i,E_MPERT)&
                 /odat(i,1,MBND_VAR))**one_th)*(1.d0 - edat(i,E_ECC))*edat(i,E_BETA) - sim_dele(i)
        epsscl = (edat(i,E_ROBJ)/rsun) * (msun/odat(i,1,MBND_VAR))**one_th

        ddat(i,1,1:ddat_ncols(i)) = epsscl*(ddat(i,1,1:ddat_ncols(i)) + epscor)

        !Get rid of small bins
        allocate(ddat_temp(ddat_ncols(i),2))
        ddat_temp = 0.d0
        k = 0
        do j = 1, ddat_ncols(i)
            if (ddat(i,2,j) .le. min_dm) then
                cycle
            endif
            k = k + 1
            ddat_temp(k,:) = ddat(i,:,j)
        enddo
        ddat(i,:,:) = 0
        do j = 1, k
            ddat(i,:,j) = ddat_temp(j,:)
        enddo
        ddat_ncols(i) = k
        deallocate(ddat_temp)

        !Divide by de to get dm/de
        do j = 1, ddat_ncols(i)
            ddat(i,2,j) = ddat(i,2,j)/(ddat(i,1,min(j+1,ddat_ncols(i))) - ddat(i,1,min(j,ddat_ncols(i)-1)))
        enddo

        !Remesh dmde into equally-spaced log bins before smoothing
        i_zero = -1
        do j = 1, ddat_ncols(i)
            if (ddat(i,1,j) .ge. 0.d0) then
                i_zero = j - 1
                exit
            endif
        enddo
        if (i_zero .le. 0) then
            print *, 'Error, zero-energy not found.'
            print *, my_pe, maxval(-ddat(i,1,1:ddat_ncols(i))), ddat_ncols(i)
            call exit(0)
        endif

        allocate(ddat_temp(ddat_ncols(i),2))
        min_neg_e = dlog10(abs(minval(ddat(i,1,1:i_zero))))
        max_pos_e = dlog10(abs(maxval(ddat(i,1,i_zero+1:ddat_ncols(i)))))
        do j = 1, i_zero
            ddat_temp(j,1) = min_neg_e + (j-1.d0)*(min_abs_e - min_neg_e)/dble(i_zero-1)
        enddo
        ddat_temp(1:i_zero,1) = -1.d1**ddat_temp(1:i_zero,1)
        do j = i_zero+1, ddat_ncols(i)
            ddat_temp(j,1) = min_abs_e + (j-i_zero-1.d0)*(max_pos_e - min_abs_e)/dble(ddat_ncols(i)-i_zero-1)
        enddo
        ddat_temp(i_zero+1:ddat_ncols(i),1) = 1.d1**ddat_temp(i_zero+1:ddat_ncols(i),1)
        call interp_flash_output(DDAT_ARR, i, 1, ddat_temp(:,1), ddat_temp(:,2))
        do j = 1, ddat_ncols(i)
            ddat(i,:,j) = ddat_temp(j,:)
        enddo
        deallocate(ddat_temp)

        !Smooth dmde here
        do s = 4, 0, -1
            if (my_pe .eq. 0) print *, 'Smoothing dmde', i, s
            allocate(ddat_temp(ddat_ncols(i),2))
            ddat_temp = 0.d0
            do j = 1, ddat_ncols(i)
                ddat_temp(j,1) = sign(dlog10(dabs(ddat(i,1,j))),ddat(i,1,j))
                ddat_temp(j,2) = dlog10(ddat(i,2,j))
            enddo
            imark = -1
            do j = 1, ddat_ncols(i)
                if (ddat_temp(j,1) .gt. 0.d0) then
                    log_offset = -ddat_temp(j,1)
                    imark = j
                endif
            enddo
            if (imark .le. 0) then
                print *, 'Error, zero energy not found.'
                call exit(0)
            endif
            log_offset = -ddat_temp(imark-1,1) + ddat_temp(imark,1)
            ddat_temp(imark:,1) = ddat_temp(imark:,1) - log_offset
            do j = 1, ddat_ncols(i)
                !ddat(i,2,j) = sum(ddat_temp(:,2)*dexp(-((ddat_temp(:,1) - ddat_temp(j,1))/kerw)**2)) / &
                !              sum(dexp(-((ddat_temp(:,1) - ddat_temp(j,1))/kerw)**2))
                denom = 0.d0
                ddat(i,2,j) = 0.d0
                do k = 1, ddat_ncols(i)
                    arg = (ddat_temp(k,1) - ddat_temp(j,1))/kerw*2.d0**s
                    if (arg .lt. -1.d0 .or. arg .gt. 1.d0) cycle
                    window = dsin(0.5*pi*(1.d0 + arg))
                    ddat(i,2,j) = ddat(i,2,j) + ddat_temp(k,2)*window
                    denom = denom + window
                enddo
                ddat(i,2,j) = ddat(i,2,j) / denom
            enddo
            ddat(i,2,:) = 10.d0**ddat(i,2,:)
            deallocate(ddat_temp)
        enddo

        imark = -1
        !print *, maxval(1.d0 - mdat(i,1:odat_nrows(i),MTOT_VAR)/mdat(i,1,MTOT_VAR))
        do j = 1, odat_nrows(i)
            if (odat(i,j,TIME_VAR) .ge. ddat_time(i)) then
                imark = j
                exit
            endif
        enddo
        if (imark .eq. -1) then
            imark = odat_nrows(i)
            if (my_pe .eq. 0) then
                print *, "Warning: dmdt is not contained within orbit.dat's time range."
                print *, 'max orbit.dat time: ', odat(i,odat_nrows(i),TIME_VAR)
                print *, 'dmdt.dat time:      ', ddat_time(i)
            endif
            !call exit(0)
        endif

        i_zero = -1
        do j = 1, ddat_ncols(i)
            if (ddat(i,1,j) .ge. 0.d0) then
                i_zero = j - 1
                exit
            endif
        enddo
        if (i_zero .le. 0) then
            print *, 'Error, zero-energy not found.'
            print *, my_pe, maxval(-ddat(i,1,1:ddat_ncols(i))), ddat_ncols(i)
            call exit(0)
        endif

        if (exclude_roche_zone) then
            !roche_ener = roche_cut*(odat(i,imark,MBND_VAR)/3.d0/edat(i,E_MPERT))**one_th*G*edat(i,E_MPERT)/rad(i,imark)
            !roche_ener = roche_cut*G*edat(i,E_MPERT)/rad(i,imark)*&
            !    (max(odat(i,imark,MBND_VAR)/(3.d0*edat(i,E_MPERT)), 1.d-8))**one_th
            roche_ener = roche_cut

            i_roche_begin = -1
            i_roche_end = -1
            do j = 1, ddat_ncols(i)
                if (i_roche_begin .eq. -1 .and. ddat(i,1,j) .ge. -roche_ener) then
                    i_roche_begin = j
                endif
                if (i_roche_begin .ne. -1 .and. &
                    ddat(i,1,j) .ge. -1.d1**(dlog10(roche_ener) - roche_range)) then
                    i_roche_end = j
                    exit
                endif
            enddo
            if (i_roche_begin .eq. -1 .or. i_roche_end .eq. -1) then
                print *, 'Error, roche mark not found.'
                print *, my_pe, maxval(-ddat(i,1,1:ddat_ncols(i))), roche_ener, ddat_ncols(i)
                call exit(0)
            endif
            call least_sq(dlog10(-ddat(i,1,i_roche_begin:i_roche_end)), &
                          dlog10(ddat(i,2,i_roche_begin:i_roche_end)), &
                fita, fitb, siga, sigb, fitchi2, fitq)
            if (my_pe .eq. 0) then
                print *, 'Roche energy cutoff', roche_ener
                print *, 'Fitted roche energy slope/beta', -fitb, edat(i,E_BETA), fitchi2
            endif
            ddat(i,2,i_roche_begin:i_zero) = 1.d1**(dlog10(ddat(i,2,i_roche_begin-1)) + fitb*&
                (dlog10(-ddat(i,1,i_roche_begin:i_zero))-dlog10(-ddat(i,1,i_roche_begin-1))))
        endif

        ! Apply cuts to early-time rise
        maxdmde = dlog10(maxval((-ddat(i,1,1:i_zero))**2.5d0*ddat(i,2,1:i_zero)))
        maxdmdeloc = maxloc((-ddat(i,1,1:i_zero))**2.5d0*ddat(i,2,1:i_zero), 1)
        ! Save maxdmdttime for use in eliminating models
        maxdmdttime(i) = pi/sqrt2*G*edat(i,E_MPERT)*dabs(ddat(i,1,maxdmdeloc))**(-1.5d0)
        i_early_begin = 1
        i_early_end = maxdmdeloc 
        do j = maxdmdeloc, 1, -1
            if (dlog10((-ddat(i,1,j))**2.5d0*ddat(i,2,j)) .lt. maxdmde - early_cut_const) then
                i_early_begin = j
                exit
            endif
        enddo
        early_e = dlog10(-ddat(i,1,i_early_begin)) - early_range
        do j = i_early_begin, maxdmdeloc
            if (dlog10(-ddat(i,1,j)) .lt. early_e) then
                i_early_end = j
                exit
            endif
        enddo

        i_early_begin = min(i_early_begin, i_early_end - 2)

        if (i_early_begin .eq. -1 .or. i_early_end .eq. -1) then
            print *, 'Error, unable to bracket early-time rise', &
                     i_early_begin, i_early_end, early_e, maxdmdeloc
            call exit(0)
        endif

        call least_sq(dlog10(-ddat(i,1,i_early_begin:i_early_end)), &
                      dlog10(ddat(i,2,i_early_begin:i_early_end)), &
                      fita, fitb, siga, sigb, fitchi2, fitq)
        !ddat(i,:,1:ddat_ncols(i)-imark+1) = ddat(i,:,imark:ddat_ncols(i))
        !ddat_ncols(i) = ddat_ncols(i) - imark + 1
        ddat(i,:,n_early_bins+1:ddat_ncols(i)-i_early_begin+1+n_early_bins) = &
            ddat(i,:,i_early_begin:ddat_ncols(i))
        ddat_ncols(i) = ddat_ncols(i) - i_early_begin + 1 + n_early_bins

        do j = n_early_bins, 1, -1
            ddat(i,1,j) = ddat(i,1,j+1) + (ddat(i,1,j+1) - ddat(i,1,j+2))
        enddo
        !ddat(i,2,1:imark+n_early_bins-1) = 1.d1**(fita + fitb*dlog10(-ddat(i,1,1:imark+n_early_bins-1)))
        do j = n_early_bins, 1, -1
            ddat(i,2,j) = 1.d1**(dlog10(ddat(i,2,j+1)) + &
                          fitb*(dlog10(-ddat(i,1,j))-dlog10(-ddat(i,1,j+1))))
        enddo

        ! Now generate the integrated dm/dt arrays
        iddat(i,1,:) = ddat(i,1,:)
        iddat(i,2,1) = ddat(i,2,1)
        rhoddat(i,1,:) = ddat(i,1,:)
        rhoddat(i,2,:) = ddat(i,2,:)*dsqrt(dabs(ddat(i,2,:)))
        rhoddat(i,2,:) = rhoddat(i,2,:)/maxval(rhoddat(i,2,:))
        do j = 2, ddat_ncols(i)
            iddat(i,2,j) = iddat(i,2,j-1) + ddat(i,2,j)!*dsqrt(dabs(ddat(i,1,j)))
        enddo
        iddat(i,2,1:ddat_ncols(i)) = iddat(i,2,1:ddat_ncols(i))/maxval(iddat(i,2,1:ddat_ncols(i)))

        maxtauloc = maxloc(ddat(i,2,1:ddat_ncols(i))*dsqrt(dabs(ddat(i,1,1:ddat_ncols(i)))), 1, &
                           mask = ddat(i,1,1:ddat_ncols(i)) .lt. 0.d0)
        maxtautime(i) = pi/sqrt2*G*edat(i,E_MPERT)*dabs(ddat(i,1,maxtauloc))**(-1.5d0)

        ! Coefficients to make interpolation faster.
        do j = 1, ddat_ncols(i)
            ddat_coeff(i,j) = (ddat(i,2,j+1) - ddat(i,2,j))/(ddat(i,1,j+1) - ddat(i,1,j))
            iddat_coeff(i,j) = (iddat(i,2,j+1) - iddat(i,2,j))/(iddat(i,1,j+1) - iddat(i,1,j))
            rhoddat_coeff(i,j) = (rhoddat(i,2,j+1) - rhoddat(i,2,j))/(rhoddat(i,1,j+1) - rhoddat(i,1,j))
        enddo
    enddo

    do i = 1, nruns
        d_emin(i) = ddat(i,1,1)
        d_emid(i) = ddat(i,1,maxloc((dabs(ddat(i,1,1:ddat_ncols(i))))**2.5d0*ddat(i,2,1:ddat_ncols(i)),1, &
                        mask = ddat(i,1,1:ddat_ncols(i)) .lt. 0.d0))
        d_emax(i) = ddat(i,1,ddat_ncols(i))
    enddo

    if (my_pe .eq. 0) then
        open(unit = fn, file = 'dmde_library.dat')
        do i = 1, nruns
            write(str,'(I9)') ddat_ncols(i)
            write(fn,*) ddat_ncols(i), edat(i,E_BETA)
            write(fn,'('//trim(str)//'G17.8E3)') (ddat(i,1,j),j=1,ddat_ncols(i))
            write(fn,'('//trim(str)//'G17.8E3)') (ddat(i,2,j),j=1,ddat_ncols(i))
        enddo
        close(fn)
    endif

    min_sim_beta = minval(edat(:,E_BETA))
    max_sim_beta = maxval(edat(:,E_BETA))

    allocate(min_model_beta(nmodels))
    allocate(max_model_beta(nmodels))
    do i = 1, nmodels
        if (i .lt. nmodels) then
            min_model_beta(i) = minval(edat(model_index(i):model_index(i+1)-1,E_BETA))
            max_model_beta(i) = maxval(edat(model_index(i):model_index(i+1)-1,E_BETA))
        else
            min_model_beta(i) = minval(edat(model_index(i):,E_BETA))
            max_model_beta(i) = maxval(edat(model_index(i):,E_BETA))
        endif
    enddo

    ! Deallocate some unneeded arrays
    deallocate(rad)
    deallocate(vel)
    deallocate(vec1)
    deallocate(vec2)
    deallocate(bndcom)
    deallocate(totcom)
    deallocate(semimaj)

    ! odat and mdat are not needed anymore
    deallocate(odat)
    deallocate(mdat)

    ! Initialize sed structures
    sed_min = dlog10(1.00001*0.3/lamb_const) ! Set by alambda
    sed_max = dlog10(0.5*max_x_xray/lamb_const) ! Set by alambda
    sed_step = 0.02d0
    sed_nsteps = nint((sed_max - sed_min) / sed_step) + 1
    sed_divs = 100
    allocate(zone_sed_table(sed_nsteps))

    ! Initialize filters and flux table
    call tdefit_print('Initializing filters')
    open(unit = fn, file = trim(exec_path) // "filters/filterlist.dat", status='old', action='read')
    read(fn, *) num_filt
    allocate(filt_names(num_filt))
    allocate(filt_files(num_filt))
    allocate(filt_len(num_filt))
    do i = 1, num_filt
        read(fn,'(A,A)') filt_names(i), filt_files(i)
        open(unit = fn2, file = trim(adjustl(exec_path)) // "filters/" // &
            trim(adjustl(filt_files(i))), status='old', action='read')
        read(fn2, *) filt_len(i)
        close(fn2)
    enddo
    close(fn)

    allocate(filt_resp(num_filt,maxval(filt_len),2))
    allocate(dfilt_resp(num_filt,maxval(filt_len)))
    allocate(filtnorm(num_filt))
    allocate(filt_lmin(num_filt))
    allocate(filt_lmax(num_filt))
    do i = 1, num_filt
        open(unit = fn2, file = trim(adjustl(exec_path)) // "filters/" // &
            trim(adjustl(filt_files(i))), status='old', action='read')
        read(fn2, *) dummy
        do j = 1, filt_len(i)
            read(fn2, *) (filt_resp(i,j,k), k=1,2)
        enddo
        close(fn2)
        filt_resp(i,:filt_len(i),1) = c/filt_resp(i,:filt_len(i),1)/angstrom
    enddo

    do i = 1, num_filt
        filt_resp(i,:filt_len(i),:) = filt_resp(i,filt_len(i):1:-1,:)
        do j = 1, filt_len(i)-1
            dfilt_resp(i,j) = (filt_resp(i,j+1,2) - filt_resp(i,j,2))/(filt_resp(i,j+1,1) - filt_resp(i,j,1))
        enddo
    enddo

    allocate(filtint(num_filt,filtintlen,3))
    allocate(filt_const(num_filt))
    do j = 1, num_filt
        bbband = filt_names(j)
        bandi = j
        i = 1
        do while (filt_resp(j,i,2) .le. 0.d0)
            i = i + 1
        enddo
        filtint(j,1,1) = filt_resp(j,max(i-1,1),1)
        i = filt_len(j)
        do while (filt_resp(j,i,2) .le. 0.d0)
            i = i - 1
        enddo
        filtint(j,filtintlen,1) = filt_resp(j,min(i+1,filt_len(j)),1)
        do i = 1, filtintlen
            filtint(j,i,1) = dble(i-1)/(dble(filtintlen-1))*(filtint(j,filtintlen,1) - filtint(j,1,1)) + filtint(j,1,1)
            filtint(j,i,2) = max(0.d0, filterfunc(filtint(j,i,1)))
        enddo
        do i = 1, filtintlen-1
            filtint(j,i,3) = (filtint(j,i+1,2) - filtint(j,i,2))/(filtint(j,i+1,1) - filtint(j,i,1))
        enddo
        filt_const(j) = dble(filtintlen)/(filtint(j,filtintlen,1) - filtint(j,1,1))
        filt_lmin(j) = dlog10(filtint(j,1,1))
        filt_lmax(j) = dlog10(filtint(j,filtintlen,1))
    enddo

    do i = 1, num_filt
        bbband = filt_names(i)
        bandi = i
        call qxgs(filtnormfunc,filt_lmin(i),filt_lmax(i),0.d0,1.d-2*bb_int_tol,filtnorm(i),err,ierr,bb_int_subdiv,neval)
    enddo

    !allocate(avint(avintlen,5))
    !avint(1,1) = 0.3/lamb_const
    !avint(avintlen,1) = 88.67744852d0/lamb_const
    !do i = 1, avintlen
    !    avint(i,1) = dble(i-1)/(dble(avintlen-1))*(avint(avintlen,1) - avint(1,1)) + avint(1,1)
    !    call alambda(avint(i,1), 0.d0, 1.d0, err, avint(i,2), avint(i,3))
    !enddo
    !do i = 1, avintlen-1
    !    avint(i,4) = (avint(i+1,2) - avint(i,2))/(avint(i+1,1) - avint(i,1))
    !    avint(i,5) = (avint(i+1,3) - avint(i,3))/(avint(i+1,1) - avint(i,1))
    !enddo
    !avconst = dble(avintlen)/(avint(avintlen,1) - avint(1,1))

    call tdefit_print('Initializing reddening arrays.')
!   Compute wave numbers of spline anchor points
!   fitz_lamspl = wavelengths of spline anchor points
!            (first is at infinity)
!   fitz_xspl = wave numbers of spline anchor points, in ascending order
!   fitz_nspl = number of anchor points

    fitz_xspl(1) = 0.0d0
    do i=1,8
       fitz_xspl(i+1) = 10000.0d0 / fitz_lamspl(i) 
    end do
    fitz_nspl = 9

    call tdefit_print('Finished initialization.')
    print *, 'Version ', VERSION, ', my_pe', my_pe
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine

