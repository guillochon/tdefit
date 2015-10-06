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
