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

function annulus_intercept(lr) result(f)
    use tdefit_data
    use constants
    use tdefit_interface, only: ang_frac


    real, intent(in) :: lr
    real :: f, r

    r = 1.e1**lr
    f = ang_frac(r)
    f = trial_fout(cur_event)*f*l10*dexp(-0.5d0*((r - ai_mu)/ai_sig)**2)
end function
