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

subroutine set_var(var, val)
    use tdefit_data
    use constants, only: halfpi
    use tdefit_interface, only: dbl2logic


    character*50, intent(in) :: var
    real, intent(in) :: val

    select case (trim(var))
        case ("toff")
            trial_toff(cur_event) = val
        case ("nhsrc")
            trial_nhsrc(cur_event) = val
        case ("rout")
            trial_rout(cur_event) = val
        case ("aspin")
            trial_aspin(cur_event) = val
        case ("phi")
            trial_phi(cur_event) = val*halfpi
        case ("beta")
            trial_beta(cur_event) = val
        case ("mh")
            trial_mh(cur_event) = val*msun
        case ("ms")
            trial_ms(cur_event) = val*msun
        case ("rsc")
            trial_rsc(cur_event) = val
        case ("alphhr")
            trial_alphhr(cur_event) = val
        case ("source_rv")
            trial_source_rv(cur_event) = val
        case ("magoff")
            trial_magoff(cur_event) = val
        case ("fcor")
            trial_fcor(cur_event) = val
        case ("model")
            trial_model(cur_event) = nint(val)
        case ("outflow_model")
            trial_outflow_model(cur_event) = nint(val)
        case ("temperature_model")
            trial_temperature_model(cur_event) = nint(val)
        case ("blr_model")
            trial_blr_model(cur_event) = nint(val)
        case ("object_type")
            trial_object_type(cur_event) = nint(val)
        case ("tlimit")
            trial_tlimit(cur_event) = val
        case ("ecor")
            trial_ecor(cur_event) = val
        case ("time_dep_rin")
            trial_time_dep_rin(cur_event) = dbl2logic(val)
        case ("time_dep_rout")
            trial_time_dep_rout(cur_event) = dbl2logic(val)
        case ("simple_bb")
            trial_simple_bb(cur_event) = dbl2logic(val)
        case ("use_fcor")
            trial_use_fcor(cur_event) = dbl2logic(val)
        case ("cap_at_edd")
            trial_cap_at_edd(cur_event) = dbl2logic(val)
        case ("full_disk_coverage")
            trial_full_disk_coverage(cur_event) = dbl2logic(val)
        case ("reprocess_temp")
            trial_reprocess_temp(cur_event) = val
        case ("rphot")
            trial_rphot(cur_event) = val
        case ("opacity")
            trial_opacity(cur_event) = val
        case ("fout")
            trial_fout(cur_event) = val
        case ("outflow_frac")
            trial_outflow_frac(cur_event) = val
        case ("mu_e")
            trial_mu_e(cur_event) = val
        case ("offset_X1")
            trial_offset_X1(cur_event) = val
        case ("offset_X2")
            trial_offset_X2(cur_event) = val
        case ("offset_GN")
            trial_offset_GN(cur_event) = val
        case ("offset_Pg")
            trial_offset_Pg(cur_event) = val
        case ("offset_Pr")
            trial_offset_Pr(cur_event) = val
        case ("offset_Pi")
            trial_offset_Pi(cur_event) = val
        case ("offset_Pz")
            trial_offset_Pz(cur_event) = val
        case ("offset_U1")
            trial_offset_U1(cur_event) = val
        case ("offset_U2")
            trial_offset_U2(cur_event) = val
        case ("offset_RO")
            trial_offset_RO(cur_event) = val
        case ("offset_Ub")
            trial_offset_Ub(cur_event) = val
        case ("offset_Um")
            trial_offset_Um(cur_event) = val
        case ("offset_Uu")
            trial_offset_Uu(cur_event) = val
        case ("offset_Uv")
            trial_offset_Uv(cur_event) = val
        case ("offset_bV")
            trial_offset_bV(cur_event) = val
        case ("temp_mult")
            trial_temp_mult(cur_event) = val
        case ("z")
            trial_z(cur_event) = val
        case ("exp_1")
            trial_exp_1(cur_event) = val
        case ("exp_2")
            trial_exp_2(cur_event) = val
        case ("exp_3")
            trial_exp_3(cur_event) = val
        case ("exp_4")
            trial_exp_4(cur_event) = val
        case ("variability")
            trial_variability(cur_event) = val
        case ("variance")
            trial_variance(cur_event) = val
        case ("viscous_time")
            trial_viscous_time(cur_event) = val
        !Non-trial vars
        case ("zscale")
            zscale = val
        case ("mag_penalty")
            mag_penalty = val
        case ("wind_phot")
            wind_phot = dbl2logic(val)
        case ("circ_phot")
            circ_phot = dbl2logic(val)
        case ("include_circ")
            include_circ = dbl2logic(val)
        case ("viscous_dmdt")
            viscous_dmdt = dbl2logic(val)
        case ("remove_extinction_corr")
            remove_extinction_corr = dbl2logic(val)
        case ("dump_all_walkers")
            dump_all_walkers = dbl2logic(val)
        case ("dump_burned_walkers")
            dump_burned_walkers = dbl2logic(val)
        case ("dump_burned_ensembles")
            dump_burned_ensembles = dbl2logic(val)
        case ("calc_acor")
            calc_acor = dbl2logic(val)
        case ("discard_failed_integrals")
            discard_failed_integrals = dbl2logic(val)
        case ("output_test")
            output_test = dbl2logic(val)
        case ("print_integral_warnings")
            print_integral_warnings = dbl2logic(val)
        case ("nensemble")
            nensemble = nint(val)
        case ("nmcsteps")
            nmcsteps = nint(val)
        case ("burn_in")
            burn_in = nint(val)
        case ("nanneal")
            nanneal = nint(val)
        case ("annealt0")
            annealt0 = val
        case ("anneal_hot_frac")
            anneal_hot_frac = val
        case ("mutation_prob")
            mutation_prob = val
        case ("mutation_min")
            mutation_min = val
        case ("roche_cut")
            roche_cut = val
        case ("roche_range")
            roche_range = val
        case ("early_range")
            early_range = val
        case ("min_best_time")
            min_best_time = val*yr
        case ("max_best_time")
            max_best_time = val*yr
        case ("nbest_times")
            nbest_times = nint(val)
        case ("stationary_redraw")
            stationary_redraw = nint(val)
        case ("nwalkers")
            nwalkers = nint(val)
        case ("nchains")
            nchains = nint(val)
        case ("mixstep")
            mixstep = nint(val)
        case ("bb_int_method")
            bb_int_method = nint(val)
        case ("df_int_method")
            df_int_method = nint(val)
        case ("use_solar_radius")
            use_solar_radius = dbl2logic(val)
        case ("include_disk")
            include_disk = dbl2logic(val)
        case ("redraw_stationary")
            redraw_stationary = dbl2logic(val)
        case ("redraw_bad_walkers")
            redraw_bad_walkers = dbl2logic(val)
        case ("amoeba_best")
            amoeba_best = dbl2logic(val)
        case ("amoeba_step")
            amoeba_step = nint(val)
        case ("redraw_using_best")
            redraw_using_best = dbl2logic(val)
        case ("exclude_roche_zone")
            exclude_roche_zone = dbl2logic(val)
        case ("lf_ms")
            lf_ms = nint(val)
        case ("lf_mh")
            lf_mh = nint(val)
        case ("lf_beta")
            lf_beta = nint(val)
        case ("lf_heii_dispersion")
            lf_heii_dispersion = nint(val)
        case ("lf_halpha_dispersion")
            lf_halpha_dispersion = nint(val)
        case ("lf_variability")
            lf_variability = nint(val)
        case ("bb_int_divs")
            bb_int_divs = nint(val)
        case ("df_int_divs")
            df_int_divs = nint(val)
        case ("extra_reflect")
            extra_reflect = dbl2logic(val)
        case ("fixed_angle_frac")
            fixed_angle_frac = dbl2logic(val)
        case ("disallow_unphysical_models")
            disallow_unphysical_models = dbl2logic(val)
        case ("annulus_width")
            annulus_width = val
        case ("kerw")
            kerw = val
        case ("amoeba_spread")
            amoeba_spread = val
        case ("amoeba_tol")
            amoeba_tol = val
        case ("bb_int_tol")
            bb_int_tol = val
        case ("df_int_tol")
            df_int_tol = val
        case default
            print *, 'Error [set_var]: Invalid variable name specified --- ', var
            call exit(0)
    end select
end subroutine

subroutine set_string(var, val)
#include "tdefit.fpp"
    use tdefit_data


    character*50, intent(in) :: var
    character*50, intent(in) :: val

    select case(var)
        case ('initial_walker_dist')
            select case(val)
                case ('random')
                    initial_walker_dist = IW_RANDOM
                case ('file_statistical')
                    initial_walker_dist = IW_STATISTICAL
                case ('file_most_probable')
                    initial_walker_dist = IW_PROBABLE
                case ('regular')
                    initial_walker_dist = IW_REGULAR
                case default
                    print *, 'Error: Invalid initial walker distribution selected.'
            end select
        case ('reprocess_model')
            select case(val)
                case ('none')
                    reprocess_model = RM_NONE
                case ('clumpy_inner')
                    reprocess_model = RM_CLUMPY_INNER
                case ('clumpy_outer')
                    reprocess_model = RM_CLUMPY_OUTER
                case ('annulus')
                    reprocess_model = RM_ANNULUS
                case ('clouds')
                    reprocess_model = RM_CLOUDS
                case ('circular')
                    reprocess_model = RM_CIRCULAR
                case default
                    print *, 'Error: Invalid reprocessing mode selected'
            end select
        case default
            print *, 'Error [set_string]: Invalid variable name specified --- ', var
            call exit(0)
    end select
end subroutine
