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

subroutine tdefit_print(message)
    use tdefit_data

    character*(*), intent(in) :: message

    if (print_children) then
        write(*,'(A12,I2,A2,A)'), '[Processor #', my_pe, '] ', trim(message)
    elseif (my_pe .eq. 0) then
        write(*,'(A)'), trim(message)
    endif
end subroutine
