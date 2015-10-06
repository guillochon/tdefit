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

SUBROUTINE indexing_re(arr,index)
    USE tdefit_util, ONLY : arth,assert_eq,tdefit_error,swap
    real, DIMENSION(:), INTENT(IN) :: arr
    integer, DIMENSION(:), INTENT(OUT) :: index
    integer, PARAMETER :: NN=15, NSTACK=50
    real :: a
    integer :: n,k,i,j,indext,jstack,l,r
    integer, DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(arr),'indexing_re')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
        if (r-l < NN) then
            do j=l+1,r
                indext=index(j)
                a=arr(indext)
                do i=j-1,l,-1
                    if (arr(index(i)) <= a) exit
                    index(i+1)=index(i)
                end do
                index(i+1)=indext
            end do
            if (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
        else
            k=(l+r)/2
            call swap(index(k),index(l+1))
            call icomp_xchg(index(l),index(r))
            call icomp_xchg(index(l+1),index(r))
            call icomp_xchg(index(l),index(l+1))
            i=l+1
            j=r
            indext=index(l+1)
            a=arr(indext)
            do
                do
                    i=i+1
                    if (arr(index(i)) >= a) exit
                end do
                do
                    j=j-1
                    if (arr(index(j)) <= a) exit
                end do
                if (j < i) exit
                call swap(index(i),index(j))
            end do
            index(l+1)=index(j)
            index(j)=indext
            jstack=jstack+2
            if (jstack > NSTACK) call tdefit_error('indexing: NSTACK too small')
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
    CONTAINS

    SUBROUTINE icomp_xchg(i,j)
    integer, INTENT(INOUT) :: i,j
    integer :: swp
    if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
    end if
    END SUBROUTINE icomp_xchg
    END SUBROUTINE indexing_re

    SUBROUTINE indexing_i4b(iarr,index)
    USE tdefit_util, ONLY : arth,assert_eq,tdefit_error,swap
    integer, DIMENSION(:), INTENT(IN) :: iarr
    integer, DIMENSION(:), INTENT(OUT) :: index
    integer, PARAMETER :: NN=15, NSTACK=50
    integer :: a
    integer :: n,k,i,j,indext,jstack,l,r
    integer, DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(iarr),'indexing_re')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
        if (r-l < NN) then
            do j=l+1,r
                indext=index(j)
                a=iarr(indext)
                do i=j-1,1,-1
                    if (iarr(index(i)) <= a) exit
                    index(i+1)=index(i)
                end do
                index(i+1)=indext
            end do
            if (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
        else
            k=(l+r)/2
            call swap(index(k),index(l+1))
            call icomp_xchg(index(l),index(r))
            call icomp_xchg(index(l+1),index(r))
            call icomp_xchg(index(l),index(l+1))
            i=l+1
            j=r
            indext=index(l+1)
            a=iarr(indext)
            do
                do
                    i=i+1
                    if (iarr(index(i)) >= a) exit
                end do
                do
                    j=j-1
                    if (iarr(index(j)) <= a) exit
                end do
                if (j < i) exit
                call swap(index(i),index(j))
            end do
            index(l+1)=index(j)
            index(j)=indext
            jstack=jstack+2
            if (jstack > NSTACK) call tdefit_error('indexing: NSTACK too small')
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
    CONTAINS

    SUBROUTINE icomp_xchg(i,j)
    integer, INTENT(INOUT) :: i,j
    integer :: swp
    if (iarr(j) < iarr(i)) then
        swp=i
        i=j
        j=swp
    end if
    END SUBROUTINE icomp_xchg
END SUBROUTINE indexing_i4b
