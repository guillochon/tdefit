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

function get_var(var) result(val)
    use tdefit_data

    character*50, intent(in) :: var
    real :: val

    select case (trim(var))
        case ("toff")
            val = trial_toff(cur_event)
        case ("nhsrc")
            val = trial_nhsrc(cur_event)
        case ("rin")
            val = trial_rin(cur_event)
        case ("rout")
            val = trial_rout(cur_event)
        case ("aspin")
            val = trial_aspin(cur_event)
        case ("phi")
            val = trial_phi(cur_event)/halfpi
        case ("beta")
            val = trial_beta(cur_event)
        case ("mh")
            val = trial_mh(cur_event)*imsun
        case ("ms")
            val = trial_ms(cur_event)*imsun
        case ("rsc")
            val = trial_rsc(cur_event)
        case ("alphhr")
            val = trial_alphhr(cur_event)
        case ("source_rv")
            val = trial_source_rv(cur_event)
        case ("magoff")
            val = trial_magoff(cur_event)
        case ("fcor")
            val = trial_fcor(cur_event)
        case ("model")
            val = dble(trial_model(cur_event))
        case ("outflow_model")
            val = dble(trial_outflow_model(cur_event))
        case ("temperature_model")
            val = dble(trial_temperature_model(cur_event))
        case ("blr_model")
            val = dble(trial_blr_model(cur_event))
        case ("object_type")
            val = dble(trial_object_type(cur_event))
        case ("tlimit")
            val = trial_tlimit(cur_event)
        case ("ecor")
            val = trial_ecor(cur_event)
        case ("time_dep_rin")
            if (trial_time_dep_rin(cur_event)) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("time_dep_rout")
            if (trial_time_dep_rout(cur_event)) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("simple_bb")
            if (trial_simple_bb(cur_event)) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("use_fcor")
            if (trial_use_fcor(cur_event)) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("cap_at_edd")
            if (trial_cap_at_edd(cur_event)) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("full_disk_coverage")
            if (trial_full_disk_coverage(cur_event)) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("reprocess_temp")
            val = trial_reprocess_temp(cur_event)
        case ("rphot")
            val = trial_rphot(cur_event)
        case ("opacity")
            val = trial_opacity(cur_event)
        case ("fout")
            val = trial_fout(cur_event)
        case ("outflow_frac")
            val = trial_outflow_frac(cur_event)
        case ("mu_e")
            val = trial_mu_e(cur_event)
        case ("offset_X1")
            val = trial_offset_X1(cur_event)
        case ("offset_X2")
            val = trial_offset_X2(cur_event)
        case ("offset_GN")
            val = trial_offset_GN(cur_event)
        case ("offset_Pg")
            val = trial_offset_Pg(cur_event)
        case ("offset_Pr")
            val = trial_offset_Pr(cur_event)
        case ("offset_Pi")
            val = trial_offset_Pi(cur_event)
        case ("offset_Pz")
            val = trial_offset_Pz(cur_event)
        case ("offset_U1")
            val = trial_offset_U1(cur_event)
        case ("offset_U2")
            val = trial_offset_U2(cur_event)
        case ("offset_RO")
            val = trial_offset_RO(cur_event)
        case ("offset_Ub")
            val = trial_offset_Ub(cur_event)
        case ("offset_Um")
            val = trial_offset_Um(cur_event)
        case ("offset_Uu")
            val = trial_offset_Uu(cur_event)
        case ("offset_Uv")
            val = trial_offset_Uv(cur_event)
        case ("offset_bV")
            val = trial_offset_bV(cur_event)
        case ("offset_bI")
            val = trial_offset_bI(cur_event)
        case ("temp_mult")
            val = trial_temp_mult(cur_event)
        case ("z")
            val = trial_z(cur_event)
        case ("exp_1")
            val = trial_exp_1(cur_event)
        case ("exp_2")
            val = trial_exp_2(cur_event)
        case ("exp_3")
            val = trial_exp_3(cur_event)
        case ("exp_4")
            val = trial_exp_4(cur_event)
        case ("variability")
            val = trial_variability(cur_event)
        case ("variance")
            val = trial_variance(cur_event)
        case ("viscous_time")
            val = trial_viscous_time(cur_event)
        case ("mdot_floor")
            val = trial_mdot_floor(cur_event)
        ! Non-trial vars
        !   Logical vars
        case ("wind_phot")
            if (wind_phot) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("circ_phot")
            if (circ_phot) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("include_circ")
            if (include_circ) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("viscous_dmdt")
            if (viscous_dmdt) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("remove_extinction_corr")
            if (remove_extinction_corr) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("dump_all_walkers")
            if (dump_all_walkers) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("dump_burned_walkers")
            if (dump_burned_walkers) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("dump_burned_ensembles")
            if (dump_burned_ensembles) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("calc_acor")
            if (calc_acor) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("discard_failed_integrals")
            if (discard_failed_integrals) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("output_test")
            if (output_test) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("print_integral_warnings")
            if (print_integral_warnings) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("use_solar_radius")
            if (use_solar_radius) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("include_disk")
            if (include_disk) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("redraw_stationary")
            if (redraw_stationary) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("redraw_bad_walkers")
            if (redraw_bad_walkers) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("amoeba_best")
            if (amoeba_best) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("redraw_using_best")
            if (redraw_using_best) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("exclude_roche_zone")
            if (exclude_roche_zone) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("extra_reflect")
            if (extra_reflect) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("fixed_angle_frac")
            if (fixed_angle_frac) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("disallow_unphysical_models")
            if (disallow_unphysical_models) then
                val = 1.d0
            else
                val = 0.d0
            endif
        case ("zscale")
            val = zscale 
        case ("mag_penalty")
            val = mag_penalty 
        case ("nensemble")
            val = nensemble 
        case ("nmcsteps")
            val = nmcsteps 
        case ("burn_in")
            val = burn_in 
        case ("nanneal")
            val = nanneal 
        case ("annealt0")
            val = annealt0 
        case ("anneal_hot_frac")
            val = anneal_hot_frac 
        case ("mutation_prob")
            val = mutation_prob 
        case ("mutation_min")
            val = mutation_min 
        case ("roche_cut")
            val = roche_cut 
        case ("roche_range")
            val = roche_range 
        case ("early_range")
            val = early_range 
        case ("min_best_time")
            val = min_best_time 
        case ("max_best_time")
            val = max_best_time 
        case ("nbest_times")
            val = nbest_times 
        case ("stationary_redraw")
            val = stationary_redraw 
        case ("nwalkers")
            val = nwalkers 
        case ("nchains")
            val = nchains 
        case ("mixstep")
            val = mixstep 
        case ("bb_int_method")
            val = bb_int_method 
        case ("df_int_method")
            val = df_int_method 
        case ("amoeba_step")
            val = amoeba_step 
        case ("lf_ms")
            val = lf_ms 
        case ("lf_mh")
            val = lf_mh 
        case ("lf_beta")
            val = lf_beta 
        case ("lf_heii_dispersion")
            val = lf_heii_dispersion 
        case ("lf_halpha_dispersion")
            val = lf_halpha_dispersion 
        case ("lf_variability")
            val = lf_variability 
        case ("bb_int_divs")
            val = bb_int_divs 
        case ("df_int_divs")
            val = df_int_divs 
        case ("annulus_width")
            val = annulus_width 
        case ("kerw")
            val = kerw 
        case ("amoeba_spread")
            val = amoeba_spread 
        case ("amoeba_tol")
            val = amoeba_tol 
        case ("bb_int_tol")
            val = bb_int_tol 
        case ("df_int_tol")
            val = df_int_tol 
        ! String variables, but stored as integers after loading
        case ("initial_walker_dist")
            val = initial_walker_dist
        case ("reprocess_model")
            val = reprocess_model
        case default
            print *, 'Error [get_var]: Invalid variable name specified --- ', var
            call exit(0)
    end select
end function
