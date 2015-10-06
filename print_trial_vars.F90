subroutine print_trial_vars()
    use tdefit_data
    use tdefit_interface, ONLY: get_var

    integer :: i, j

    j = 0
    do i = 1, nvars
        if (var_events(i) .ne. 0 .and. var_events(i) .ne. cur_event) cycle
        j = j + 1
        if (cur_event .gt. 1 .and. i .gt. 1) then
            if (var_names(i) .ne. var_names(i-1)) cycle
        endif
        if (var_types(i) .ne. 1) then
            write(*, '(A25)', advance='no'), '   ' // trim(var_names(i)) // ': '
            write(*, '(ES12.5)', advance='no'), get_var(var_names(i))
        else
            write(*, '(A25)', advance='no'), '   log ' // trim(var_names(i)) // ': '
            write(*, '(ES12.5)', advance='no'), dlog10(get_var(var_names(i)))
        endif
        if (j .eq. 5) then
            write(*,*)
            j = 0
        endif
    enddo
    write(*,*) ''
end subroutine
