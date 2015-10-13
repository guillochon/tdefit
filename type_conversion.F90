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

real function logic2dbl(a)
    logical, intent(in) :: a
  
    if (a .eqv. .true.) then
        logic2dbl = 1.d0
    else
        logic2dbl = 0.d0
    end if
end function logic2dbl

logical function int2logic(a)
    integer, intent(in) :: a

    if (a .eq. 1) then
        int2logic = .true.
    else
        int2logic = .false.
    endif
end function

logical function dbl2logic(a)
    real, intent(in) :: a
    integer :: tempi

    tempi = nint(a)


    if (tempi .eq. 1) then
        dbl2logic = .true.
    else
        dbl2logic = .false.
    endif
 end function
