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

subroutine acor(walkers, atime)
    use tdefit_data

#include "tdefit.fpp"

    real, dimension(:,:,:), intent(in) :: walkers
    real, dimension(size(walkers, 3)), intent(out) :: atime

    real, dimension(:,:,:), allocatable :: cov
    real, dimension(:,:), allocatable :: mean, avcov
    integer :: i, j, k, s, nw, nv, ns
    integer :: ml
    
    ml = min(10, size(walkers, 1)-1)

    ns = size(walkers, 1)
    nw = size(walkers, 2)
    nv = size(walkers, 3)

    allocate(mean(nw,nv))
    allocate(cov(ml+1,nw,nv))
    allocate(avcov(ml+1,nv))

    ! Follows Foreman-Mackey 2012 and Akeret 2012
    mean = sum(walkers, 1)/dble(ns)

    cov = 0.d0

    do i = 1, ns - ml
    do j = 1, nw
    do k = 1, nv
    do s = 0, ml
        cov(s+1,j,k) = cov(s+1,j,k) + (walkers(i,j,k)-mean(j,k))*(walkers(i+s,j,k)-mean(j,k))/dble(ns-ml)
    enddo
    enddo
    enddo
    enddo

    !do i = 1, ns
    !do j = 1, nw
    !do k = 1, nv
    !    cov(i,j,k) = sum((walkers(:ns-i+1,j,k)-mean(j,k))*(walkers(i:,j,k)-mean(j,k)))/dble(ns)
    !enddo
    !enddo
    !enddo

    avcov = sum(cov, 2)/dble(nw)

    atime = 1.d0 + 2.d0*sum(avcov,1)/avcov(1,:)

    deallocate(mean,cov,avcov)
end subroutine
