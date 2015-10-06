function bbflux(bbfunc, band, T, z, nh, nhsrc) result(flux)
    use tdefit_interface, ONLY: alambdaz, filterintfunc, trapezoid, is_abnormal
    use adapt_quad, ONLY: qxgs
    use tdefit_data


    character*2, intent(in) :: band
    real, intent(in) :: T, z, nh, nhsrc
    real, external :: bbfunc

    real :: flux, fluxtemp, err, minlnu, maxlnu!, midlnu, dummy, numin, numax
    integer :: neval, ierr, i, qng_div

    ierr = 0
    qng_div = 1

    flux = 0.d0
    if (T .le. 1.d0 .or. T .gt. 1.d12 .or. is_abnormal(T)) then
        return
    endif

    bbband = band
    bbtemp = T
    bbnh = nh
    bbnhsrc = nhsrc
    bbz = z
    bb1pz = 1.d0 + z
    bbpenalty = .false.

    if (band .eq. 'Lb') then
        flux = sigma_b*bbtemp**4
        return
    else
        bandi = 0
        do i = 1, num_filt
            if (band .eq. filt_names(i)) then
                bandi = i
                exit
            endif
        enddo
    endif
    if (bandi .eq. 0) then
        print *, "Error, requested band name not found:", band
        call exit(0)
    endif

    maxlnu = min(dlog10(bbtemp*wein_max/(x_const*bb1pz)),filt_lmax(bandi)/filt_limit_mult,&
                 dlog10(max_x_xray/lamb_const) - dlog10(bb1pz))
    minlnu = filt_lmin(bandi)*filt_limit_mult

    if (maxlnu .lt. minlnu) then
        flux = 0.d0
        return
    endif

    select case (band)
        case ('Hl', 'HL', '51')
            bbmultbynu = .true.
        case default
            bbmultbynu = .false.
    endselect

    if (bb_int_method .eq. 1) then
        bbcnt = bbcnt + 1
        call qag(bbfunc,minlnu,maxlnu,0.d0,bb_int_tol,bb_int_mode,flux,err,neval,ierr)
        if (ierr .eq. 1) then
            print *, "Error, maximum number of steps executed in qag [bbflux].", flux, neval
            call exit(0)
        endif
    elseif (bb_int_method .eq. 2 .or. bb_int_method .eq. 4) then
        !if (midlnu .lt. maxlnu .and. midlnu .gt. minlnu) then
        !    call qng(bbfunc,minlnu,midlnu,0.d0,bb_int_tol,dummy,err,neval,ierr)
        !    if (ierr .ne. 0) then
        !        bbfailcnt = bbfailcnt + 1
        !    endif
        !    bbcnt = bbcnt + 1
        !    call qng(bbfunc,midlnu,maxlnu,0.d0,bb_int_tol,flux,err,neval,ierr)
        !    if (ierr .ne. 0) then
        !        bbfailcnt = bbfailcnt + 1
        !    endif
        !    flux = flux + dummy
        !else
        !flux = 0.d0
        !numin = minlnu
        !numax = maxlnu
        !do while (numax .le. maxlnu)
        !    bbcnt = bbcnt + 1
        !    ierr = 0
        !    call qng(bbfunc,numin,numax,0.d0,bb_int_tol,fluxtemp,err,neval,ierr)
        !    if (ierr .ne. 0) then
        !        !print *, bbfailcnt, bbband, bbtemp, fluxtemp, err, err/fluxtemp, numin, numax, filterintfunc(10.d0**numin), &
        !        !    filterintfunc(10.d0**numax), alambdaz(10.d0**numin, bbz, bbnh, bbnhsrc), &
        !        !    alambdaz(10.d0**numax, bbz, bbnh, bbnhsrc)
        !        !call exit(0)
        !        bbfailcnt = bbfailcnt + 1
        !        qng_div = 2*qng_div
        !        if (qng_div .gt. bb_int_qng_div) exit
        !        numax = numin + 0.5d0*(numax - numin)
        !        cycle
        !    endif
        !    flux = flux + fluxtemp
        !    numin = numax
        !    numax = numin + (maxlnu - minlnu)/dble(qng_div)
        !enddo

        ierr = 0
        call qng(bbfunc,minlnu,maxlnu,0.d0,bb_int_tol,fluxtemp,err,neval,ierr)
        if (ierr .ne. 0) then
            bbfailcnt = bbfailcnt + 1
        endif
        !endif
        !if (ierr .eq. 1) then
        !    print *, "Error, maximum number of steps executed in qng [bbflux]."
        !    call exit(0)
        !endif
    endif

    if (bb_int_method .eq. 3 .or. (bb_int_method .eq. 4 .and. ierr .ne. 0)) then
        bbcnt = bbcnt + 1
        !if (midlnu .lt. maxlnu .and. midlnu .gt. minlnu) then
        !    call qxgs(bbfunc,minlnu,midlnu,0.d0,bb_int_tol,dummy,err,ierr,bb_int_subdiv,neval)
        !    if (ierr .ne. 0) then
        !        bbfailcnt = bbfailcnt + 1
        !    endif
        !    call qxgs(bbfunc,midlnu,maxlnu,0.d0,bb_int_tol,flux,err,ierr,bb_int_subdiv,neval)
        !    if (ierr .ne. 0) then
        !        bbfailcnt = bbfailcnt + 1
        !    endif
        !    flux = flux + dummy
        !else
            ierr = 0
            call qxgs(bbfunc,minlnu,maxlnu,0.d0,bb_int_tol,flux,err,ierr,bb_int_subdiv,neval)
            if (ierr .ne. 0) then
                bbfailcnt = bbfailcnt + 1
                print *, "Warning, maximum number of sub-divisions created in qxgs [bbflux].", flux, neval, ierr
                !call exit(0)
            endif
        !endif
    endif

    if (bb_int_method .eq. 5) then
        bbcnt = bbcnt + 1
        call trapezoid(bbfunc, minlnu, maxlnu, bb_int_divs, flux)
    endif
    !if (verbose) then
    !    print *, 'bb int evals: ', neval
    !endif

    if (bbpenalty) then
        flux = 0.d0
        return
    endif

    select case (band)
        case ('X1', 'X2')
            flux = flux/h
        case ('Hl', 'HL')
        case default
            ! Necessary to calculate AB magnitudes
            flux = flux/filtnorm(bandi)
    endselect

    if (flux .ne. flux) then
        print *, 'bbflux: nan flux', T, bandi, filtnorm(bandi)
        call exit(0)
    endif
end function
