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

module tdefit_data
    use constants
    use types

#include "tdefit.fpp"

    character*100, parameter :: path_file = "paths.par"
    real, parameter :: local_rv = 3.1d0 !3.1 for MW

    logical, save :: bbmultbynu, bbpenalty
    real, save :: bbtemp, bbnh, bbnhsrc, bbz, bb1pz
    character*2, save :: bbband

    logical, save :: print_max_disk = .false.
    logical, save :: record_max_disk = .false.
    logical, save :: print_likelihood = .true.
    logical, save :: print_annealing = .false.
    real, save :: dffb, dfmd, dftime, dfenv, dfcovering, dftemp, dfri, dfro, &
                  dfmaxdisktemp, dfmaxdiskflux, dfintval, dfabsorb
    character*2, save :: dfband

    integer, save :: num_filt
    real, allocatable, dimension(:,:,:), save :: filt_resp
    real, allocatable, dimension(:,:), save :: dfilt_resp
    real, allocatable, dimension(:), save :: filtnorm, filt_lmin, filt_lmax
    integer, save, allocatable, dimension(:) :: filt_len
    character*2, save, allocatable, dimension(:) :: filt_names
    character*100, save, allocatable, dimension(:) :: filt_files

    ! Variables to store event data.
    integer, save :: event_n
    integer, save :: cur_event, cur_npts, cur_blrpts
    integer, save, allocatable, dimension(:) :: event_npts, event_blrpts, event_nbest_bands
    real, save, allocatable, dimension(:) :: event_nh, event_claimed_z, event_min_aspin
    character*100, save, allocatable, dimension(:) :: event_fnames
    character*2, save, allocatable, dimension(:,:) :: event_bands, event_blr_bands, event_best_bands
    real, save, allocatable, dimension(:,:) :: event_times, event_ABs, event_devs, event_errs, event_weights, &
                                               event_blr_times, event_blr_vels
    integer, save, allocatable, dimension(:,:) :: event_types, event_penalties
    logical, save, allocatable, dimension(:,:) :: event_blr_exists
    integer, save, allocatable, dimension(:) :: event_nhcorr, event_restframe

    real, save :: reduced_chi2_const

    ! Variables that are generated for each parameter combination.
    real, save :: trial_penalty
    real, save, allocatable, dimension(:) :: trial_mh, trial_ms, trial_rs, trial_beta, trial_aspin, trial_nhsrc, trial_phi, &
                                             trial_bh_rms, trial_rp, trial_toff, trial_rout, trial_rg, trial_gmh, trial_ms0, &
                                             trial_mh0, trial_rs0, trial_rsc, trial_alphhr, trial_eps_edd, trial_source_rv, &
                                             trial_magoff, trial_fcor, trial_tlimit, trial_ecor, trial_r_isco, trial_reprocess_temp, &
                                             trial_offset_X1, trial_offset_X2, &
                                             trial_offset_GN, trial_offset_Pg, trial_offset_Pr, trial_offset_Pi, trial_offset_Pz, &
                                             trial_offset_U1, trial_offset_U2, trial_offset_RO, trial_offset_Ub, &
                                             trial_offset_Um, trial_offset_Uu, trial_offset_Uv, trial_offset_bV, &
                                             trial_fout, trial_mu_e, trial_ledd, trial_outflow_frac, trial_bh_rbs, trial_r_ibco, &
                                             trial_temp_mult, trial_z, trial_dl, trial_exp_1, trial_exp_2, trial_exp_3, trial_exp_4, &
                                             trial_variability, trial_variance, trial_variability2, trial_variance2, &
                                             trial_opacity, trial_rphot, trial_yms, trial_y1, trial_y2, trial_y3, &
                                             trial_1pz, trial_viscous_time
    integer, save, allocatable, dimension(:) :: trial_model
    integer, save, allocatable, dimension(:) :: trial_outflow_model, trial_object_type, trial_temperature_model, trial_blr_model
    logical, save, allocatable, dimension(:) :: trial_time_dep_rin
    logical, save, allocatable, dimension(:) :: trial_time_dep_rout
    logical, save, allocatable, dimension(:) :: trial_simple_bb
    logical, save, allocatable, dimension(:) :: trial_use_fcor
    logical, save, allocatable, dimension(:) :: trial_cap_at_edd
    logical, save, allocatable, dimension(:) :: trial_full_disk_coverage

    real, save, allocatable, dimension(:,:) :: trial_mags, trial_fbs, trial_times, &
                                               trial_menv, trial_routs, trial_rphots, &
                                               trial_mdots

    !integer, parameter :: nTs = 40, nnhs = 4, nzs = 20
    !real, parameter :: minT = 1.d3, maxT = 1.d8
    !real, parameter :: minnh = 1.d19, maxnh = 1.d23
    !real, parameter :: minz = 0.0, maxz = 0.37

    !real :: Ts(nTs), nhs(nnhs), zs(nzs)
    !real :: lfluxes_table(nbands,nTs,nzs,nnhs,nnhs)
    !real :: fluxes_table(nbands,nTs,nzs,nnhs,nnhs)

!!!!! USER parameterS.
    ! The following parameters are typically edited by the user before each run.

    ! Output/execution preferences.
    logical, parameter :: force_reload_data = .false.
    logical, parameter :: print_children = .false.
    logical, parameter :: dev_print = .false.
    logical, parameter :: preferred_mode = .false.  !Will read in a file describing the best model for each event, 
                                                   !generating a well-sampled set of light curves for each.
    logical, parameter :: output_restframe = .true.
    logical :: remove_extinction_corr
    logical :: dump_all_walkers
    logical :: dump_burned_walkers
    logical :: dump_burned_ensembles
    logical :: calc_acor
    logical :: discard_failed_integrals
    logical :: output_test
    logical :: print_integral_warnings
    logical :: print_extra = .false.
    integer :: initial_walker_dist

    ! Extra bands to output
    character*2, dimension(2), parameter :: extra_bands = ['Lb', '51']
    integer, parameter                   :: nextra_bands = size(extra_bands)

    ! Adjustable physical parameters.
    logical, parameter          :: time_weighted = .false.        ! Default: false
    logical, parameter          :: penalize_early_time = .true.   ! Default: true
    logical, parameter          :: no_source_extinction = .false. ! Default: false
    logical, parameter          :: include_recomb = .true.        ! Default: false (CURRENTLY NOT USED)
    logical, parameter          :: boost_accepts = .false.
    logical, save :: wind_phot
    logical, save :: circ_phot
    logical, save :: include_circ
    logical, save :: viscous_dmdt
    logical, save :: use_solar_radius
    logical, save :: include_disk
    logical, save :: redraw_stationary
    logical, save :: redraw_bad_walkers
    logical, save :: redraw_using_best
    logical, save :: extra_reflect
    logical, save :: disallow_unphysical_models
    logical, save :: fixed_angle_frac
    real, save :: annulus_width, df_rphot, df_temp_mult
    real, save :: df_intercept_frac, ai_mu, ai_sig
    integer, save :: reprocess_model
    integer, save :: annulus_divs = 1000

    logical, save :: dfreject = .false.

    logical, save :: amoeba_best
    real, save :: amoeba_spread, amoeba_tol
    integer, save :: amoeba_step

    real, parameter :: upp_lim_err = 0.001d0 !Assumed errors on upper limits with no defined errors, in magnitudes.
    !real, parameter :: tidal_ms_0 = 1.0d0, tidal_ms_1 = 1.1d0 !This set of parameters is used to mimic degen. core disruptions.

    ! Parameters relating to imported simulation data.
    real, parameter :: dele_time = 2.5d5 !Time after pericenter to sample change in energy
    real, parameter :: dele_dmtot = 1.0d-8 !Record the change in energy when this much mass has left the box.
                                            !The change in energy starts to be unreliable when mass leaves the box.
    real, save :: kerw        !In number of bins
    real, save :: roche_cut, roche_range, early_range
    real, parameter :: early_cut_const = 2.0d0 ! Number of orders of magnitude to cutoff at early times from dm/de.
    real, parameter :: min_dm = 5.d24
    real, parameter :: min_abs_e = 11.d0 ! Minimum log value in E before smoothing
    integer, parameter :: n_early_bins = 2000
    integer, parameter :: dmdt_viscl = 200 ! When integrating for viscous model, use this many bins.
    
    ! Likelihood function options
    integer, parameter :: hard_penalties = 2 ! 0: Soft penalties, 1: Gradual hardening while annealing, 2: Always hard
    integer, save :: lf_ms
    integer, save :: lf_mh
    integer, save :: lf_beta
    integer, save :: lf_heii_dispersion
    integer, save :: lf_halpha_dispersion
    integer, save :: lf_variability
    real, parameter :: lf_band_sigma = 0.d0

    ! Flux model parameters
    ! Integration method
    integer :: bb_int_method
    integer :: df_int_method
    ! Used for method 1
    integer, parameter :: bb_int_mode = 1
    integer, parameter :: df_int_mode = 1
    ! Used for method 2
    integer, parameter :: bb_int_qng_div = 1
    ! Used for method 3
    integer, parameter :: bb_int_subdiv = 100000
    integer, parameter :: df_int_subdiv = 100000
    ! Used for method 5
    integer, save :: bb_int_divs
    integer, save :: df_int_divs
    real, parameter :: filt_limit_mult = 1.d0 + 1.d-6
    real :: bb_int_tol
    real :: df_int_tol
    real, parameter :: wein_max = dlog(huge(1.d0)) !Warning: This value seems to cause issues with the integrator
    !real, parameter :: wein_max = 30.d0
    logical, parameter          :: separate_outer_zone = .false.
    logical, save               :: exclude_roche_zone

    ! Parameters for outputting fit ensemble.
    integer, save :: nensemble

    !Parameters for outputting the best-fitting light curves.
    integer, save :: nbest_times
    real, save :: min_best_time, max_best_time

    ! MCMC parameters
    logical, parameter :: red_blue = .true.
    integer, save :: nwalkers !Number of walkers per processor
    integer, save :: nmcsteps
    integer, save :: burn_in, nburn
    integer, save :: nanneal
    integer, save :: nchains
    integer, save :: mixstep
    integer, save :: stationary_redraw
    real, save :: zscale
    real, save :: mag_penalty
    real, save :: anneal_hot_frac
    real, save :: mutation_prob !Probability of random mutation if both pts have same model
                                            !during annealing phase
    real, save :: mutation_min  !Minimum probability of mutation, all phases
    real, parameter :: scale_exp_fac = 10.d0

    ! MCMC shared variables
    real, save :: annealt0
    real, save :: annealvar

    ! Set which variables to maximize over and which to generate a table from.
    integer :: nvars
    integer, save, allocatable, dimension(:) :: var_types, all_var_types, var_events
    character*50, save, allocatable, dimension(:) :: var_names, all_var_names
    real, save, allocatable, dimension(:) :: min_search, max_search

!!!!! end USER parameterS

!!!!! STATIC parameterS.
    ! These parameters should rarely change from run to run.

    integer, parameter :: maxnrows = 100000
    integer, parameter :: avintlen = 1000
    integer, parameter :: filtintlen = 1000
    integer, parameter :: bisect_lim = 100
    real, parameter :: min_r_fac = 1.0 + 1.0d-3
!!!!! end STATIC parameterS.

!!!!! SHARED VARIABLES.
    ! These variables are used internally.

    ! Variables for reading in simulation data
    character*100, save                                :: exec_path, output_path, input_path, event_path, binary_path
    real, save                                         :: min_sim_beta, max_sim_beta !Are set to min/max sim beta in init.
    real, save                                         :: bhdisk_fcor, eps_edd, max_aspin, avconst, mctemp
    real, save                                         :: first_accretion_time, full_reflect_r, orb_period
    integer, save                                      :: nruns, nrows, maxnt, ptneta, my_pe, n_pe
    integer, save                                      :: bandi, nstep
    integer, save                                      :: dfcnt, dffailcnt, bbcnt, bbfailcnt, nmodels, bbcalls
    logical, save                                      :: magfail, verbose, restframe_mode = .false.
    logical, save                                      :: dfreprocess
    type(bnptr), dimension(:), allocatable             :: ll
    type(bandnode), pointer                            :: cur
    character*50, save, dimension(:), allocatable      :: dir_names
    logical, save, dimension(:), allocatable           :: var_locks
    real, save, dimension(:), allocatable              :: ddat_bnd, ddat_time, model_beta_destroy, min_model_beta, max_model_beta
    real, save, dimension(:), allocatable              :: mdat_row, odat_row, d_emin, d_emid, d_emax, mctmin, mctmax, filt_const
    real, save, dimension(:), allocatable              :: tmaxcap, locked_trial_vars, maxdmdttime, maxtautime
    real, save, dimension(:), target, allocatable      :: sim_dele, sim_mfinal
    real, save, dimension(:,:), allocatable            :: avint, ddat_coeff, iddat_coeff, rhoddat_coeff
    real, save, dimension(:,:), target, allocatable    :: ptdat, edat
    real, save, dimension(:,:), allocatable            :: mu_tot, mu_bnd, v1, v2, rad, vel, semimaj
    real, save, dimension(:,:,:), target, allocatable  :: mdat, odat, ddat, iddat, rhoddat, filtint
    real, save, dimension(:,:,:), allocatable          :: vec1, vec2, bndcom, mpolecom, totcom, barycenter
    integer, dimension(:), save, allocatable           :: mdat_nrows, odat_nrows, ddat_ncols

    ! Annealing variables
    integer, save, allocatable, dimension(:)            :: model_index
    logical, save                                       :: use_soft_penalties = .false.

    ! Root parameters
    real, save                              :: disk_t

    ! SED parameters
    logical, save                                         :: produce_seds = .false., make_sed
    real, save, allocatable, dimension(:,:,:) :: sed_table
    real, save, allocatable, dimension(:)     :: zone_sed_table
    real, save                                :: sed_min, sed_max, sed_step
    integer, save                                         :: sed_nsteps
    integer, save                                         :: sed_index
    integer, save                                         :: sed_divs

    ! Fitzpatrick constants
    integer, parameter :: fitz_nspltot = 9
    real fitz_lamspl(8)
    data fitz_lamspl /26500.0d0, 12200.0d0, 6000.0d0, 5470.0d0, &
                      4670.0d0, 4110.0d0, 2700.0d0, 2600.0d0/
    real, save ::    fitz_xspl(fitz_nspltot) 
    integer, save :: fitz_nspl

    ! Temporary storage parameters
    real, save, dimension(:), allocatable   :: flux_save
    character*100, save                                 :: warnings = ''
!!!!! end SHARED VARIABLES
end module
