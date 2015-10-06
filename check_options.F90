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

subroutine check_options
    use tdefit_data
    
    if (my_pe .eq. 0) then
        if (.not. dump_burned_walkers .and. dump_burned_ensembles) then
            print *, 'Error: Dump burned walkers must be enabled if dump burned ensembles is enabled'
            call exit(0)
        endif
    endif
end subroutine
