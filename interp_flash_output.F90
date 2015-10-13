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

subroutine interp_flash_output(arr, run, ind, x, val)
    use constants
    use tdefit_data
    use tdefit_interface, ONLY: bisect

#include "tdefit.fpp"

    integer, intent(in) :: arr, run, ind
    real, dimension(:), intent(in) :: x
    real, dimension(size(x)), intent(out) :: val
    real, dimension(:), pointer :: xarr, yarr

    integer :: i, j, l, n, ierr, nx
    
    nx = size(x)

    select case (arr)
        case (MDAT_ARR)
            n = mdat_nrows(run)
            xarr => mdat(run,1:n,TIME_VAR)
            yarr => mdat(run,1:n,ind)
        case (ODAT_ARR)
            n = odat_nrows(run)
            xarr => odat(run,1:n,TIME_VAR)
            yarr => odat(run,1:n,ind)
        case (DDAT_ARR)
            n = ddat_ncols(run)
            xarr => ddat(run,1,1:n)
            yarr => ddat(run,2,1:n)
        case (IDDAT_ARR)
            n = ddat_ncols(run)
            xarr => iddat(run,1,1:n)
            yarr => iddat(run,2,1:n)
        case (MFINAL_ARR)
            n = nruns
            xarr => edat(1:n,E_BETA)
            yarr => sim_mfinal(1:n)
        case (DELE_ARR)
            n = nruns
            xarr => edat(1:n,E_BETA)
            yarr => sim_dele(1:n)
        case (PT_ARR)
            n = ptneta
            xarr => ptdat(1:n,1)
            yarr => ptdat(1:n,2)
        case default
            print *, "Unrecognized array selected."
            call exit(ierr)
    end select

    !if (arr .eq. DDAT_ARR) then
    !    do l = 1, nx
    !        i = floor((x(l)-ddat(run,1,1))/(ddat(run,1,n)-ddat(run,1,1))*n)+1 
    !        val(l) = ddat(run,2,i) + (x(l) - ddat(run,1,i))*ddat_coeff(run,i)
    !    enddo
    !elseif (arr .eq. IDDAT_ARR) then
    !    do l = 1, nx
    !        i = floor((x(l)-iddat(run,1,1))/(iddat(run,1,n)-iddat(run,1,1))*n)+1 
    !        val(l) = iddat(run,2,i) + (x(l) - iddat(run,1,i))*iddat_coeff(run,i)
    !    enddo
    !else
        do l = 1, nx
#ifdef DEBUG
            if (x(l) .lt. xarr(1)) then
                print *, 'Warning: x less than min xind'
                print *, arr, run, ind, x(l), xarr(1)
            elseif (x(l) .gt. xarr(n)) then
                print *, 'Warning: x greater than max xind'
                print *, arr, run, ind, x(l), xarr(n)
            endif
#endif

            i = bisect(xarr, x(l))
            if (i .eq. 0) then
                val(l) = 0.d0
            else
                j = i + 1
                val(l) = yarr(i) + (x(l) - xarr(i))/(xarr(j) - xarr(i))*(yarr(j) - yarr(i))
            endif
        enddo
    !endif

end subroutine
