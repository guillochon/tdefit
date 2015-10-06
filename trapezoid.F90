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

recursive subroutine trapezoid(func,minx,maxx,div,val)
#include "tdefit.fpp"
    !use, intrinsic :: ieee_arithmetic
    use tdefit_interface, only: is_abnormal

    real, external    :: func
    real, intent(in)  :: minx, maxx
    integer, intent(in)           :: div
    real, intent(out) :: val

    integer :: i
    real :: coeff, temporary

    val = 0.d0
    coeff = (maxx - minx)/(dble(div) - 1.d0)
    do i = 1, div
        temporary = func(minx + coeff*(dble(i)-1.d0))
#ifdef DEBUG
        !if (.not. ieee_is_normal(temporary)) then
        if (is_abnormal(temporary)) then
            print *, 'val abnormal', coeff, maxx, minx, div, temporary, i
            call exit(0)
        endif
#endif
        if (i .eq. 1 .or. i .eq. div) then
            val = val + 0.5d0*temporary
        else
            val = val + temporary
        endif
    enddo
    val = val * coeff
end subroutine
