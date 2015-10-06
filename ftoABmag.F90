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
