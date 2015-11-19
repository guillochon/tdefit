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

module tdefit_interface
    interface
        subroutine integrate_df(func, zmin, zmax, mag, divsin)
            use tdefit_data
            use adapt_quad, ONLY: qxgs

            real, external    :: func
            real, intent(in) :: zmin, zmax
            real, intent(out) :: mag
            integer, intent(in), optional :: divsin
        end subroutine

        function dffunc(bbfunc, lr) result(flux)       
            real, intent(in) :: lr
            real, external :: bbfunc
            real :: flux
        end function

        function obs_df_func(lr) result(flux)       
            real, intent(in) :: lr
            real :: flux
        end function

        function arv_df_func(lr) result(flux)       
            real, intent(in) :: lr
            real :: flux
        end function

        function src_df_func(lr) result(flux)       
            real, intent(in) :: lr
            real :: flux
        end function

        function obs_bb_func(lnu) result(flux)       
            real, intent(in) :: lnu
            real :: flux
        end function

        function arv_bb_func(lnu) result(flux)       
            real, intent(in) :: lnu
            real :: flux
        end function

        function src_bb_func(lnu) result(flux)       
            real, intent(in) :: lnu
            real :: flux
        end function

        function filtnormfunc(lnu) result(flux)       
            real, intent(in) :: lnu
            real :: flux
        end function

        function ang_frac(r) result(frac)
            real, intent(in) :: r
            real :: frac
        end function

        function annulus_intercept(lr) result(f)
            real, intent(in) :: lr
            real :: f
        end function

        function bisect(arr, val, retj) result(i)
            real, intent(in), dimension(:) :: arr
            real, intent(in) :: val
            logical, intent(in), optional :: retj
            integer :: i
        end function
        subroutine load_event(e)
            integer, intent(in)     :: e
        end subroutine

        function filterfunc(nu) result(frac)
            real, intent(in) :: nu
            real :: frac
        end function

        function filterintfunc(nu) result(frac)
            real, intent(in) :: nu
            real :: frac
        end function

        function alambdaz(nu, z, nh, nhsrc) result(al)
            real, intent(in) :: nu, z, nh, nhsrc
            real :: al
        end function

        subroutine alambda(nu, nh, rv, al)
            real, intent(in) :: nu, nh, rv
            real, intent(out) :: al
        end subroutine

        function avintfunc(nu, rv) result(al)
            real, intent(in) :: nu, rv
            real :: al
        end function

        function bbflux(bbfunc, band, T, z, nh, nhsrc) result(flux)
            character*2, intent(in) :: band
            real, intent(in) :: T, z, nh, nhsrc
            real, external :: bbfunc
            real :: flux
        end function

        function bbsed(bbfunc, T, z, nh, nhsrc) result(sed)
            use tdefit_data
            real, intent(in) :: T, z, nh, nhsrc
            real, external :: bbfunc
            real, dimension(sed_nsteps) :: sed
        end function

        function cosmo_ez(z) result(ez)
            real, intent(in) :: z
            real :: ez
        end function

        function cosmo_dc(z) result(dc)
            real, intent(in) :: z
            real :: dc
        end function

        function cosmo_dm(z) result(dm)
            real, intent(in) :: z
            real :: dm
        end function

        function cosmo_da(z) result(da)
            real, intent(in) :: z
            real :: da
        end function

        function cosmo_dl(z) result(dl)
            real, intent(in) :: z
            real :: dl
        end function

        function set_trial_vars(x) result(penalty)
            real, dimension(:), intent(in) :: x
            real, dimension(size(x)) :: xlimited
            real :: penalty
        end function

        subroutine print_trial_vars
        end subroutine

        subroutine set_derived_trial_vars
        end subroutine

        subroutine set_string(var, val)
            character*50, intent(in) :: var
            character*50, intent(in) :: val
        end subroutine

        subroutine set_var(var, val)
            character*50, intent(in) :: var
            real, intent(in) :: val
        end subroutine

        function get_var(var) result(val)
            character*50, intent(in) :: var
            real :: val
        end function

        subroutine md_arr(o,i)
            integer, intent(in) :: o
            integer, dimension(:), intent(out) :: i
        end subroutine

        subroutine get_sim_index(betafrac, bi, ei)
            real, intent(out) :: betafrac
            integer, intent(out) :: bi, ei
        end subroutine

        subroutine dmdt(tdes, dm, add_delay, im, rhom, mode, ades)
            real, dimension(:), intent(in) :: tdes
            real, dimension(size(tdes)), intent(out) :: dm
            logical, intent(in) :: add_delay
            real, dimension(size(tdes)), intent(out), optional :: im, rhom
            integer, intent(in), optional :: mode
            real, intent(in), optional :: ades
        end subroutine

        function im_root(t)
            real, intent(in) :: t
            real :: im_root
        end function

        subroutine interp_flash_output(arr, run, var, x, val)
            integer, intent(in) :: arr, run, var
            real, dimension(:), target, intent(in) :: x
            real, dimension(size(x)), intent(out) :: val
        end subroutine

        function disk_temp(mdot, r) result(temp)
            real, intent(in) :: mdot, r
            real :: temp
        end function

        function disk_temp_root(r)
            real, dimension(:), intent(in) :: r
            real, dimension(size(r)) :: disk_temp_root
        end function

        function ftoABmag(f) result(ab)
            real, intent(in), dimension(:) :: f
            real, dimension(size(f)) :: ab
        end function

        function ftomag(f) result(ab)
            real, intent(in), dimension(:) :: f
            real, dimension(size(f)) :: ab
        end function

        subroutine bandmag(times, fbs, mdots, bands, mags, penalties, routs, rphots)
            real, dimension(:), intent(in) :: times
            real, dimension(:), intent(in) :: fbs
            real, dimension(:), intent(in) :: mdots
            character*2, dimension(:), intent(in) :: bands
            real, dimension(size(mdots)), intent(out) :: mags
            integer, dimension(size(mdots)), intent(out) :: penalties
            real, dimension(size(mdots)), intent(out), optional :: routs, rphots
        end subroutine

        function radius(m) result(r)
            real, intent(in) :: m
            real :: r
        end function

        function ydev(x) result(dev)
            real, dimension(:), intent(in) :: x
            real :: dev
        end function

        function annealydev(x) result(dev)
            real, dimension(:), intent(in) :: x
            real :: dev
        end function

        function magdev(x, max_likelihood) result(dev)
            real, dimension(:), intent(in) :: x
            logical, intent(in), optional :: max_likelihood
            real :: dev
        end function

        function kroupa(mstar) result(prob)
            real, intent(in) :: mstar
            real :: prob
        end function

        function chabrier(mstar) result(prob)
            real, intent(in) :: mstar
            real :: prob
        end function

        function likelihood() result(l)
            real :: l
        end function

        recursive subroutine trapezoid(func,minx,maxx,div,val)
            real, external    :: func
            real, intent(in)  :: minx, maxx
            integer, intent(in)           :: div
            real, intent(out) :: val
        end subroutine

        real function logic2dbl(a)
            logical, intent(in) :: a
        end function

        logical function int2logic(a)
            integer, intent(in) :: a
        end function

        logical function dbl2logic(a)
            real, intent(in) :: a
        end function

        logical function is_normal(a)
            real, intent(in) :: a
        end function

        logical function is_abnormal(a)
            real, intent(in) :: a
        end function

        subroutine acor(walkers, atime)
            real, dimension(:,:,:), intent(in) :: walkers
            real, dimension(size(walkers, 3)), intent(out) :: atime
        end subroutine

        subroutine draw_random_walker(w1, wchainl, wchainr, w2, wp)
            integer, intent(in) :: w1, wchainl, wchainr
            integer, intent(out) :: w2, wp
        end subroutine

        subroutine set_event(event)
            integer, intent(in) :: event
        end subroutine

        subroutine tdefit_print(message)
            character*(*), intent(in) :: message
        end subroutine

        subroutine load_defaults(mode)
            integer, intent(in) :: mode
        end subroutine

        subroutine load_user_vars(mode)
            integer, intent(in) :: mode
        end subroutine
        
        subroutine write_vars(fn)
            integer, intent(in) :: fn
        end subroutine

        function get_band_type(band) result(bt)
            character*2, intent(in) :: band
            character*1, bt
        end function
    end interface
end module
