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

function ang_frac(r) result(frac)
    use constants
    use tdefit_data

    real, intent(in) :: r
    real :: frac, t_visc

    frac = 1.d0
    if (.not. trial_full_disk_coverage(cur_event)) then
        if (r .lt. 2.d0*trial_rp(cur_event) + &
            (8.d0*((dftime - first_accretion_time)/pi)**2*trial_gmh(cur_event))**one_th) then
            frac = min(dsqrt(trial_rp(cur_event)/r), 1.d0)
        else
            frac = 0.d0
        endif
        if (.not. fixed_angle_frac .and. r .gt. 2.d0*trial_rp(cur_event)) then
            t_visc = twopi*dsqrt(r**3.d0/trial_gmh(cur_event))/trial_alphhr(cur_event)
            frac = max(0.d0,min(frac, (dftime - (first_accretion_time + &
                   pi*dsqrt(0.125d0*(r - 2.d0*trial_rp(cur_event))**3/&
                   trial_gmh(cur_event))))/t_visc))
        endif
    endif
end function
