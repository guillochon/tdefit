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

subroutine alambda(nu, nh, rv, al)
    use constants
    use tdefit_interface, only : bisect
    use tdefit_data

    real, intent(in) :: nu, nh, rv
    real, intent(out) :: al
    real :: x, fa, fb, y, y2, y3, y4, y5, y6, y7
    real, dimension(1) :: xarr, alarr
    integer :: i

    !if (nh .le. 0) then
    !    al = 0.0
    !    return
    !endif

    x = lamb_const*nu

    if (x .ge. 0.d0 .and. x .lt. x_lyman) then
        xarr = x
        alarr = al
        call law_odo(xarr, 1, .true., rv, alarr)
        !call law_cal(x, 1, .true., al)
        al = nhconst*nh*(alarr(1)/rv)
        !print *, x, al, alambdaold
    elseif (x .ge. x_lyman .and. x .le. max_x_xray) then
        y = max(h*nu/keV, 0.03d0)
        i = bisect(mm83(:,1),y)
        al = 1.d-24*(mm83(i,2) + mm83(i,3)*y + mm83(i,4)*y**2)/y**3
        ! For less than 0.03 keV, just assume cross-section scales as E^-3.
        ! http://ned.ipac.caltech.edu/level5/Madau6/Madau1_2.html
        if (x .lt. min_x_xray) then
            al = al*(min_x_xray/x)**3
        endif
        al = mag_fac*il10*nh*al
    else
        print *, 'alambda: nu out of range, x = ', x
        call exit(0)
    endif

    if (al .lt. 0.d0) then
#ifdef PRINT_PENALTY_REASONS
        print *, 'Error: alambda less than 0!'
#endif
        bbpenalty = .true.
    endif
end subroutine

function avintfunc(nu, rv) result(al)
    use tdefit_data, only: avint, avintlen, avconst, bbband
    use tdefit_interface, only: get_band_type
    use constants
    !use tdefit_interface, ONLY : bisect


    real, intent(in) :: nu, rv
    real :: al, nuval
    integer :: i

    !if (nh .le. 0) then
    !    al = 0.0
    !    return
    !endif

    !i = bisect(avint(:,1), nu)
    i = min(max(ceiling(avconst*(nu - avint(1,1))),1),avintlen)
    !i = ceiling(avconst*(nu - avint(1,1)))

    nuval = nu - avint(i,1)

    if (get_band_type(bbband) == 'X') then
        al = mag_fac*il10*(nuval*avint(i,4) + avint(i,2))
    else
        al = nhconst*(nuval*avint(i,4) + avint(i,2) + &
                      (nuval*avint(i,5) + avint(i,3))/rv)
    endif
end function
