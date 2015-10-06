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

subroutine integrate_df(func, zmin, zmax, mag, divsin)
    use tdefit_data
    use tdefit_interface, ONLY: trapezoid
    use adapt_quad, ONLY: qxgs

    real, external    :: func
    real, intent(in) :: zmin, zmax
    real, intent(out) :: mag
    integer, intent(in), optional :: divsin

    integer :: neval, ierr, divs
    real :: err

    if (present(divsin)) then
        divs = divsin
    else
        divs = df_int_divs
    endif

    if (df_int_method .eq. 1) then
        call qag(func,zmin,zmax,0.d0,df_int_tol,df_int_mode,mag,err,neval,ierr)
        if (ierr .eq. 1) then
            print *, "Error, maximum number of steps executed in qag [bandmag, 1]."
            call exit(0)
        endif
    elseif (df_int_method .eq. 2 .or. df_int_method .eq. 4) then
        call qng(func,zmin,zmax,0.d0,df_int_tol,mag,err,neval,ierr)
        if (ierr .ne. 0 .and. df_int_method .ne. 4) then
            dffailcnt = dffailcnt + 1
        endif
    endif

    if (df_int_method .eq. 3 .or. (df_int_method .eq. 4 .and. ierr .ne. 0)) then
        call qxgs(func,zmin,zmax,0.d0,df_int_tol,mag,err,ierr,df_int_subdiv,neval)
        if (ierr .ne. 0) then
            if (print_integral_warnings) print *, "Warning, maximum number of sub-divisions created in qxgs [bandmag, 1]."
            dffailcnt = dffailcnt + 1
        endif
    endif

    if (df_int_method .eq. 5 .or. ierr .ne. 0) then
        call trapezoid(func, zmin, zmax, df_int_divs, mag)
    endif
end subroutine
