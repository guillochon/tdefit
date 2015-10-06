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

function bbsed(bbfunc, T, z, nh, nhsrc) result(sed)
    use tdefit_data


#include "tdefit.fpp"

    real, intent(in) :: T, z, nh, nhsrc
    real, external :: bbfunc
    real, dimension(sed_nsteps) :: sed
    integer :: i

    bbmultbynu = .false.
    bbtemp = T
    bbz = z
    bb1pz = 1.d0 + z
    bbnh = nh
    bbnhsrc = nhsrc
    bbpenalty = .false.

    do i = 1, sed_nsteps
        sed(i) = bbfunc(sed_min + (i-1)*sed_step)
    enddo

    if (bbpenalty) then
        print *, 'Error generating SED: Penalty in bbfunc.'
        call exit(0)
    endif
end function
