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

function bisect(arr, val, retj) result(i)
    use tdefit_data, only : bisect_lim

    real, intent(in), dimension(:) :: arr
    real, intent(in) :: val
    logical, intent(in), optional :: retj
    integer :: i, j, k

    i=1   
    j=size(arr)
    if (val <= arr(1)) then
        i=1
        return
    elseif (val >= arr(j)) then
        i=j-1
        return
    endif

    do 
        k=(i+j)/2   
        if (val < arr(k)) then 
            j=k  
        else
            i=k
        end if
        if (i+1 >= j) exit
    end do

    if (present(retj)) then
        if (retj) i = j
    endif
end function
