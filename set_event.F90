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
