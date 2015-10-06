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

module constants
    
    real, parameter :: G = 6.67428d-8
    real, parameter :: msun = 1.98892d33
    real, parameter :: rsun = 6.955d10
    real, parameter :: c = 2.99792458d10
    real, parameter :: pc = 3.08568025d18
    real, parameter :: day = 86400.d0
    real, parameter :: yr = 3.1556926d7
    real, parameter :: h = 6.626068d-27
    real, parameter :: kb = 1.3806503d-16
    real, parameter :: mp = 1.672621777d-24
    real, parameter :: eV = 1.60217646d-12
    real, parameter :: keV = 1.d3*eV
    real, parameter :: micron = 1.d-4
    real, parameter :: angstrom = 1.d-8
    real, parameter :: one_th = 1.d0/3.d0, two_th = 2.d0/3.d0, four_th = 4.d0/3.d0, five_th = 5.d0/3.d0, &
                                   three_halfs = 3.d0/2.d0, five_halfs = 5.d0/2.d0, one_fourth = 1.d0/4.d0
    real, parameter :: pi      = 3.141592653589793238462643383279502884197
    real, parameter :: twopi   = 6.283185307179586476925286766559005768394
    real, parameter :: halfpi  = 1.57079632679489661923132169163975144209858
    real, parameter :: threepi = 3.d0*pi
    real, parameter :: fourpi  = 4.d0*pi
    real, parameter :: eightpi = 8.d0*pi
    real, parameter :: euler   = 0.5772156649015328606065120900824024310422
    real, parameter :: sqrt2   = 1.41421356237309504880168872420969807856967
    real, parameter :: isqrt2  = 1.d0/sqrt2
    real, parameter :: l10 = dlog(10.d0)
    real, parameter :: il10 = 1.d0/l10
    real, parameter :: pi_G = pi*G
    real, parameter :: sigma_b = 5.670373d-5
    real, parameter :: pi_sigma_b = pi*sigma_b
    real, parameter :: fourpi_sigma_b = fourpi*sigma_b
    real, parameter :: eightpi_sigma_b = eightpi*sigma_b
    real, parameter :: sigma_t = 6.65245854533d-25
    real, parameter :: c2 = c**2
    real, parameter :: flux_const = 2.0d0*h/c2*pi
    real, parameter :: x_const = h/kb
    real, parameter :: x_lyman = 10.9737316d0
    real, parameter :: fitz_xcutuv = 10000.0d0 / 2700.0d0
    real, parameter :: fitz_x0    = 4.596d0
    real, parameter :: fitz_gamma = 0.99d0
    real, parameter :: fitz_c3    = 3.23d0
    real, parameter :: fitz_c4    = 0.41d0
    real, parameter :: lamb_const = micron/c
    real, parameter :: edd_const = fourpi*G*mp*c/sigma_t
    real, parameter :: max_fcor = (10.d0/3.d0)**0.82d0
    real, parameter :: fcor_const = (1.d0/3.d4)**0.82d0
    real, parameter :: lhuge = dlog(huge(1.d0))
    real, parameter :: min_x_xray = 0.03d0*micron*kev/(h*c)
    real, parameter :: max_x_xray = 20.d0*micron*kev/(h*c)
    real, parameter :: mag_fac = 1.d2**0.2d0

    ! Inverted constants, for speed
    real, parameter :: imsun = 1.d0/msun
    real, parameter :: irsun = 1.d0/rsun
    real, parameter :: ic = 1.d0/c

    ! Cosmology parameters
    real :: omega_k = 0.d0 !This is not a parameter to avoid divide by 0 compile errors.
    real, parameter :: omega_m = 0.27d0
    real, parameter :: omega_l = 0.73d0
    real, parameter :: hubb_h = 0.71d0
    real, parameter :: cosmo_dH = 9.26d27 / hubb_h

    ! IMF parameters
    real, parameter :: kc1 = 0.08d0, kc2 = 0.5d0

    ! Reddening parameters
    real, parameter :: nhconst = 1.d0/2.21d21 ! Guver & Ozel 2009

    ! Morrison and McCammon 1983 data
    real, dimension(15, 4) :: mm83 = reshape( &
                                                 [ 0.03,0.0,0.0,0.0, &
                                                   0.1,17.3,608.1,-2150.0, &
                                                   0.284,34.6,267.9,-476.1, &
                                                   0.4,78.1,18.8,4.3, &
                                                   0.532,71.4,66.8,-51.4, &
                                                   0.707,95.5,145.8,-61.1, &
                                                   0.867,308.9,-380.6,294.0, &
                                                   1.303,120.6,169.3,-47.7, &
                                                   1.84,141.3,146.8,-31.5, &
                                                   2.471,202.7,104.7,-17.0, &
                                                   3.21,342.7,18.7,0.0, &
                                                   4.038,352.2,18.7,0.0, &
                                                   7.111,433.9,-2.4,0.75, &
                                                   8.331,629.0,30.9,0.0, &
                                                   10.0,701.2,25.2,0.0 &
                                                 ], [15, 4])
end module constants
