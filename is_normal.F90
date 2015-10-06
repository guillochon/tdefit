!This file is part of TDEFit.

!TDEFit is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TDEFit is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with TDEFit.  If not, see <http://www.gnu.org/licenses/>.

logical function is_normal(a)
    real, intent(in) :: a

    if (a .ne. a) then
        is_normal = .false.
    else
        is_normal = .true.
    endif
end function

logical function is_abnormal(a)
    use tdefit_interface, only: is_normal
    real, intent(in) :: a

    if (is_normal(a)) then
        is_abnormal = .false.
    else
        is_abnormal = .true.
    endif
end function
