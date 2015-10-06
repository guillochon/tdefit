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

subroutine sort(arr)
    USE tdefit_util, ONLY : swap,tdefit_error
    real, dimension(:), intent(inout) :: arr
    integer, parameter :: NN=15, NSTACK=50
    real :: a
    integer :: n,k,i,j,jstack,l,r
    integer, dimension(NSTACK) :: istack
    n=size(arr)
    jstack=0
    l=1
    r=n
    do
        if (r-l < NN) then
            do j=l+1,r
                a=arr(j)
                do i=j-1,l,-1
                    if (arr(i) <= a) exit
                    arr(i+1)=arr(i)
                end do
                arr(i+1)=a
            end do
            if (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
        else
            k=(l+r)/2
            call swap(arr(k),arr(l+1))
            call swap(arr(l),arr(r),arr(l)>arr(r))
            call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
            call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
            i=l+1
            j=r
            a=arr(l+1)
            do
                do
                    i=i+1
                    if (arr(i) >= a) exit
                end do
                do
                    j=j-1
                    if (arr(j) <= a) exit
                end do
                if (j < i) exit
                call swap(arr(i),arr(j))
            end do
            arr(l+1)=arr(j)
            arr(j)=a
            jstack=jstack+2
            if (jstack > NSTACK) call tdefit_error('sort: NSTACK too small')
            if (r-i+1 >= j-l) then
                istack(jstack)=r
                istack(jstack-1)=i
                r=j-1
            else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
            end if
        end if
    end do
end subroutine sort
