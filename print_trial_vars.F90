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
