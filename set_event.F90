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

subroutine set_event(event)
#include "tdefit.fpp"
    use tdefit_interface, only: tdefit_print
    use tdefit_data

    integer, intent(in) :: event

    if (event .le. 0 .or. event .gt. event_n) then
        call tdefit_print('Error: Invalid event ID specified.')
    endif
    cur_event = event
    cur_npts = event_npts(event)
    cur_blrpts = event_blrpts(event)
end subroutine
