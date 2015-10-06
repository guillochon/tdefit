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
