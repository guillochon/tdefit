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

subroutine load_user_vars(mode)
    use tdefit_data
    use tdefit_interface, only: set_var, set_string, logic2dbl

    integer, intent(in) :: mode

    character*50 :: var_name, str_value
    character*200 :: buffer
    real :: dbl_value, dbl_mvalue
    logical :: log_value, log_mvalue, is_log, is_search_var
    integer :: fn, stat, int_value, int_mvalue, var_kind, ngrid, i, &
               n, nbuf, j, loc, loc2, begk, endk, m, var_trial, nb, ne, o, &
               locspace, loctab, loccomma

    fn = 15
    open(unit = fn, file = "tdefit.par", status='old', action='read')
    ! Read through file twice, first is to simply count number of minimizing variables
    do i = 1, 2
        stat = 0
        n = 0
        linedo: do
            read(fn, '(A200)', iostat=stat, size=nbuf, advance='no') buffer
            if (is_iostat_eor(stat)) then
                j = 0
                loc = 1
                is_search_var = .false.
                ngrid = 1
                is_log = .false.
                locspace = 1
                loctab = 1
                loccomma = 1
                do while(locspace .ne. 0 .or. loctab .ne. 0 .or. loccomma .ne. 0)
                    locspace = index(trim(buffer(loc:)), " ")
                    loctab = index(trim(buffer(loc:)), char(9))
                    loccomma = index(trim(buffer(loc:)), ",")
                    if (locspace .eq. 1 .or. loctab .eq. 1 .or. loccomma .eq. 1) then
                        loc = loc + 1
                        cycle
                    endif
                    if (j .eq. 0) then
                        if (buffer(loc:loc) .eq. '#') cycle linedo
                    endif
                    j = j + 1
                    if (locspace .eq. 0 .and. loctab .eq. 0 .and. loccomma .eq. 0) then
                        begk = loc
                        endk = nbuf
                    else
                        loc2 = minval( (/ locspace, loctab, loccomma /), &
                            mask = (/ locspace, loctab, loccomma /) .ne. 0)
                        begk = loc
                        endk = loc + loc2 - 2
                        loc = loc + loc2
                    endif

                    select case (j)
                        case (1)
                            read(buffer(begk:endk),*) var_name
                        case (2)
                            read(buffer(begk:endk),*) var_kind
                        case (3)
                            read(buffer(begk:endk),*) var_trial
                        case (4)
                            if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                            if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                            select case (var_kind)
                                case (0)
                                    read(buffer(begk:endk),*) log_value
                                    dbl_value = logic2dbl(log_value)
                                case (1)
                                    read(buffer(begk:endk),*) int_value
                                    dbl_value = dble(int_value)
                                case (2)
                                    read(buffer(begk:endk),*) dbl_value
                                case (3)
                                    read(buffer(begk:endk),*) str_value
                            end select
                        case (5)
                            if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                            if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                            is_search_var = .true.
                            select case (var_kind)
                                case (0)
                                    read(buffer(begk:endk),*) log_mvalue
                                    dbl_mvalue = logic2dbl(log_mvalue)
                                case (1)
                                    read(buffer(begk:endk),*) int_mvalue
                                    dbl_mvalue = dble(int_mvalue)
                                case (2)
                                    read(buffer(begk:endk),*) dbl_mvalue
                                case (3)
                                    print *, 'Error: String variable cannot be search variable.', var_name
                                    call exit(0)
                            end select
                        case (6)
                            if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                            if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                            read(buffer(begk:endk),*) ngrid
                        case (7)
                            if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                            if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                            read(buffer(begk:endk),*) is_log
                    end select
                enddo
                if (is_search_var) then
                    if (var_trial .eq. 2) then
                        nb = n + 1
                        ne = n + event_n
                    else
                        nb = n + 1
                        ne = n + 1
                    endif
                    n = ne
                    if (i .eq. 2) then
                        do o = nb, ne
                            var_names(o) = var_name
                            if (var_trial .eq. 2) then
                                var_events(o) = o - nb + 1
                            endif
                            min_search(o) = dbl_value
                            max_search(o) = dbl_mvalue
                            select case (var_kind)
                                case (0)
                                    var_types(o) = 3
                                case (1)
                                    var_types(o) = 2
                                case (2)
                                    if (is_log) then
                                        var_types(o) = 1
                                    else
                                        var_types(o) = 0
                                    endif
                            end select
                            if (o .eq. nb) then
                                do m = 1, size(all_var_types)
                                    if (trim(all_var_names(m)) .eq. trim(var_name)) then
                                        all_var_types(m) = var_types(o)
                                        exit
                                    endif
                                enddo
                                if (my_pe .eq. 0) then
                                    select case (var_kind)
                                        case (0)
                                            write(*,'(A20,L10,X,L10)') var_name, log_value, log_mvalue
                                        case (1)
                                            write(*,'(A20,I10,X,I10)') var_name, int_value, int_mvalue
                                        case (2)
                                            write(*,'(A20,ES10.3,X,ES10.3)') var_name, dbl_value, dbl_mvalue
                                    end select
                                endif
                            endif
                        enddo
                    endif
                else
                    if (i .eq. 2) then
                        if (mode .eq. 0 .and. var_trial .ne. 0) cycle
                        if (mode .ne. 0 .and. var_trial .eq. 0) cycle
                        select case (var_kind)
                            case (0)
                                call set_var(var_name, dbl_value)
                                if (my_pe .eq. 0) write(*,'(A20,L10)') var_name, log_value
                            case (1)
                                call set_var(var_name, dbl_value)
                                if (my_pe .eq. 0) write(*,'(A20,I10)') var_name, int_value
                            case (2)
                                call set_var(var_name, dbl_value)
                                if (my_pe .eq. 0) write(*,'(A20,E10.3)') var_name, dbl_value
                            case (3)
                                call set_string(var_name, str_value)
                                if (my_pe .eq. 0) write(*,'(A20,A)') var_name, str_value
                        end select
                    endif
                endif
            elseif (is_iostat_end(stat)) then
                exit
            endif
        enddo linedo
        if (i .eq. 1) then
            rewind(fn)
            ! This is currently hacked to work for one event only, still need to
            ! figure out how to set free parameter arrays for multiple events.
            if (mode .eq. 1 .and. cur_event .eq. 1) then
                nvars = n
                allocate(var_names(nvars))
                allocate(var_types(nvars))
                allocate(var_events(nvars))
                allocate(min_search(nvars))
                allocate(max_search(nvars))

                var_events = 0
            endif
        endif
    enddo
    close(fn)
end subroutine
