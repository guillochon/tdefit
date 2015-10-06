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

subroutine draw_random_walker(w1, wchainl, wchainr, w2, wp)
    use tdefit_data

    integer, intent(in) :: w1, wchainl, wchainr
    integer, intent(out) :: w2, wp

    real :: dummy
    w2 = w1
    wp = my_pe + 1
    do while (w2 .eq. w1 .and. wp .eq. my_pe + 1)
        call random_number(dummy)
        if (red_blue .and. mod(n_pe,2) .eq. 0) then
            if (wp .le. n_pe/2) then
                wp = n_pe/2 + int(dummy*(n_pe/2)) + 1
            else
                wp = int(dummy*(n_pe/2)) + 1
            endif
        else
            wp = int(dummy*n_pe) + 1
        endif

        call random_number(dummy)
        if (mod(nstep, mixstep) .eq. 0 .or. nstep .gt. nanneal) then
            w2 = int(dummy*nwalkers) + 1
        else
            w2 = int(dummy*(wchainr-wchainl+1)) + wchainl
        endif
    enddo
end subroutine
