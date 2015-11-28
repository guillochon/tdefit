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

subroutine write_vars(fn)
#include "tdefit.fpp"
#include "version.fpp"
    use tdefit_data, only: all_var_names, all_var_types
    use tdefit_interface, only: get_var

    integer, intent(in) :: fn
    integer :: i
    do i = 1, size(all_var_names)
        write(fn,'(A,1X)',advance='no') trim(all_var_names(i))
    enddo
    write(fn,'(A)') ''
    do i = 1, size(all_var_names)
        write(fn,'(G18.8,1X)',advance='no') get_var(all_var_names(i))
    enddo
    write(fn,'(A)') ''
    do i = 1, size(all_var_names)
        if (all_var_types(i) .ne. -1) then
            write(fn,'(I1,1X)',advance='no') 1
        else
            write(fn,'(I1,1X)',advance='no') 0
        endif
    enddo
    write(fn,'(A)') ''
    do i = 1, size(all_var_names)
        write(fn,'(I2,1X)',advance='no') all_var_types(i)
    enddo
end subroutine
