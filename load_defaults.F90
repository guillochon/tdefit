subroutine load_defaults(mode)
    use tdefit_data, only: exec_path, my_pe, all_var_names, all_var_types
    use tdefit_interface, only: set_var, set_string

    integer, intent(in) :: mode

    character*50 :: var_name, str_value
    real :: dbl_value
    logical :: log_value
    integer :: fn, stat, int_value, var_type, var_trial, n
    character*50, dimension(200) :: temp_var_names

    fn = 14
    stat = 0
    n = 0
    open(unit = fn, file = trim(exec_path) // "defaults.par", status='old', action='read')
    do while (.not. is_iostat_end(stat))
        read(fn, *, iostat=stat) var_name, var_type
        if (is_iostat_end(stat)) exit
        n = n + 1
        temp_var_names(n) = var_name
        backspace(fn)
        select case (var_type)
            case (0)
                read(fn, *, iostat=stat) var_name, var_type, var_trial, log_value
                if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                if (log_value) then
                    call set_var(var_name, 1.d0)
                else
                    call set_var(var_name, 0.d0)
                endif
                if (my_pe .eq. 0) write(*,'(A,L)'), var_name // ': ', log_value
            case (1)
                read(fn, *, iostat=stat) var_name, var_type, var_trial, int_value
                if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                call set_var(var_name, dble(int_value))
                if (my_pe .eq. 0) write(*,'(A,I10)'), var_name // ': ', int_value
            case (2)
                read(fn, *, iostat=stat) var_name, var_type, var_trial, dbl_value
                if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                call set_var(var_name, dbl_value)
                if (my_pe .eq. 0) write(*,'(A,ES10.3)'), var_name // ': ', dbl_value
            case (3)
                read(fn, *, iostat=stat) var_name, var_type, var_trial, str_value
                if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                call set_string(var_name, str_value)
                if (my_pe .eq. 0) write(*,'(A,A)'), var_name // ': ' // adjustl(trim(str_value))
        end select
    enddo

    if (mode .eq. 0) then
        allocate(all_var_names(n))
        allocate(all_var_types(n))
        all_var_names = temp_var_names(1:n)
        all_var_types = -1
    endif

    close(fn)
end subroutine
