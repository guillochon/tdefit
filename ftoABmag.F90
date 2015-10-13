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

function ftoABmag(f) result(ab)
    use constants
    use tdefit_data, only: trial_phi, trial_magoff, trial_dl, cur_event
    use tdefit_interface, only: cosmo_dl


    real, intent(in), dimension(:) :: f
    real, dimension(size(f)) :: ab

    ab = f/(fourpi*trial_dl(cur_event)**2)
    where (ab .le. 0.d0)
        ab = huge(1.d0)
    elsewhere
        ab = -48.6d0 - mag_fac*(dlog10(ab)) + trial_magoff(cur_event)
    endwhere
end function

function ftomag(f) result(ab)
    use constants
    use tdefit_data, only: trial_phi, trial_magoff, trial_dl, cur_event
    use tdefit_interface, only: cosmo_dl


    real, intent(in), dimension(:) :: f
    real, dimension(size(f)) :: ab

    ab = f/(fourpi*trial_dl(cur_event)**2)
    where (ab .le. 0.d0)
        ab = huge(1.d0)
    elsewhere
        ab = -mag_fac*(dlog10(ab)) + trial_magoff(cur_event)
    endwhere
end function
