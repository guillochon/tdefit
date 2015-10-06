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

function disk_temp(mdot, r) result(temp)
    use constants
    use tdefit_data

#include "tdefit.fpp"

    real, intent(in) :: mdot, r
    real :: temp!, f, ledd

    real r0
    real cc,bb,y,part1,part2,part3

    !f = 1.d0 - dsqrt(3.d0*trial_rg(cur_event)/r)
    !ledd = 1.3d38*trial_mh(cur_event)/msun

    ! Simple model
    !temp = (3.d0*trial_gmh(cur_event)*mdot*f/(eightpi_sigma_b*r**3.d0))**one_fourth

    ! Strubbe 2009
    !temp = (3.d0*trial_gmh(cur_event)*mdot*f/(eightpi_sigma_b*r**3.d0)/&
    !    (0.5d0 + dsqrt(0.25d0 + 1.5d0*f*(mdot*c2/ledd)**2.d0/(r/rg)**2.d0)))**one_fourth

    ! Done 2011 (http://heasarc.nasa.gov/xanadu/xspec/models/optxagn.html)
    r0 = 2.d0 * r / trial_rg(cur_event)
    y=dsqrt(r0)

    part3=3.0*((trial_y3(cur_event)-trial_aspin(cur_event))**2)*dlog((y-trial_y3(cur_event))/(trial_yms(cur_event)-trial_y3(cur_event)))
    part3=part3/(y*trial_y3(cur_event)*(trial_y3(cur_event)-trial_y1(cur_event))*(trial_y3(cur_event)-trial_y2(cur_event)))
    part2=3.0*((trial_y2(cur_event)-trial_aspin(cur_event))**2)*dlog((y-trial_y2(cur_event))/(trial_yms(cur_event)-trial_y2(cur_event)))
    part2=part2/(y*trial_y2(cur_event)*(trial_y2(cur_event)-trial_y1(cur_event))*(trial_y2(cur_event)-trial_y3(cur_event)))
    part1=3.0*((trial_y1(cur_event)-trial_aspin(cur_event))**2)*dlog((y-trial_y1(cur_event))/(trial_yms(cur_event)-trial_y1(cur_event)))
    part1=part1/(y*trial_y1(cur_event)*(trial_y1(cur_event)-trial_y2(cur_event))*(trial_y1(cur_event)-trial_y3(cur_event)))

    cc=1.0-trial_yms(cur_event)/y-(3.0*trial_aspin(cur_event)/(2.0*y))*dlog(y/trial_yms(cur_event))-part1-part2-part3
    bb=1.0-3.0/r0+2.0*trial_aspin(cur_event)/(r0**1.5)
      
    temp=3.0*trial_gmh(cur_event)*mdot
    temp=temp/(eightpi_sigma_b*r**3)
    temp=(temp*cc/bb)**0.25

    ! Added to test "super-virial" disk temperature implied by teardrop model
    !temp=temp*(2.99d0/1.99d0)*(trial_mh(cur_event)/(1.d6*msun))**(-1.d0/12.d0)
end function

function disk_temp_root(r)
    use tdefit_data, only : dfmd, disk_t, trial_r_isco, trial_rg, trial_bh_rms, &
                            cur_event
    use tdefit_interface, only : disk_temp
    real, dimension(:), intent(in) :: r 
    real, dimension(size(r)) :: disk_temp_root 
    real :: rad

    rad = max(r(1), trial_r_isco(cur_event))

    disk_temp_root = disk_temp(dfmd, rad) * ((rad - r(1)) / rad + 1.d0) - disk_t
end function
