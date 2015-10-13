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

function alambdaz(nu, z, nh, nhsrc) result(avfrac)
    use tdefit_interface, ONLY: alambda
    use tdefit_data, ONLY: trial_source_rv, local_rv, bbband, no_source_extinction, cur_event


    real, intent(in) :: nu, z, nh, nhsrc
    real :: avfrac, avsrcfrac, dummy

    if (bbband .eq. 'Lb') then
        avfrac = 1.d0
    elseif (no_source_extinction .or. nhsrc .eq. 0.d0) then
        call alambda(nu, nh, local_rv, avfrac)
        avfrac = 10.d0**(-0.4d0*avfrac)
    else
        call alambda(nu, nh, local_rv, avfrac)
        call alambda(nu*(1.d0 + z), nhsrc, trial_source_rv(cur_event), avsrcfrac)
        avfrac = 10.d0**(-0.4d0*(avfrac+avsrcfrac))
    endif
end function
