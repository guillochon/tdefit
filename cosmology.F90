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

function cosmo_ez(z) result(ez)
    use constants
    real, intent(in) :: z
    real :: ez

    ez = dsqrt(omega_m*(1.d0 + z)**3 + omega_k*(1.d0 + z)**2 + omega_l)
end function

function cosmo_dc(z) result(dc)
    use tdefit_interface, ONLY: cosmo_ez
    use constants
    real, parameter :: tol = 1.d-6
    real, intent(in) :: z
    real :: dc, err
    integer :: neval, ierr

    call qag(cosmo_ez,0.d0,z,tol,tol,3,dc,err,neval,ierr)
    dc = dc*cosmo_dH
end function

function cosmo_dm(z) result(dm)
    use tdefit_interface, ONLY: cosmo_dc
    use constants
    real, intent(in) :: z
    real :: dm

    if (omega_k .gt. 0.d0) then
        dm = cosmo_dH/dsqrt(omega_k)*dsinh(dsqrt(omega_k)*cosmo_dc(z)/cosmo_dH)
    elseif (omega_k .lt. 0.d0) then
        dm = cosmo_dH/dsqrt(dabs(omega_k))*dsin(dsqrt(dabs(omega_k))*cosmo_dc(z)/cosmo_dH)
    else
        dm = cosmo_dc(z)
    endif
end function

function cosmo_da(z) result(da)
    use tdefit_interface, ONLY: cosmo_dm
    real, intent(in) :: z
    real :: da

    da = cosmo_dm(z) / (1.d0 + z)
end function

function cosmo_dl(z) result(dl)
    use tdefit_interface, ONLY: cosmo_dm
    real, intent(in) :: z
    real :: dl

    dl = (1.d0 + z) * cosmo_dm(z)
end function
