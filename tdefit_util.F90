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

module tdefit_util
    use constants
    integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
    interface
        subroutine least_sq(x,y,a,b,siga,sigb,chi2,q,sig)
        real, dimension(:), intent(in) :: x,y
        real, intent(out) :: a, b, siga, sigb, chi2, q
        real, dimension(:), optional, intent(in) :: sig
        end subroutine least_sq
    end interface
    interface gamma_inc
        function gamma_inc(p, x)
            real ( kind = 8 ) p
            real ( kind = 8 ) x
            real ( kind = 8 ) gamma_inc
        end function gamma_inc
    end interface
    interface assert
        module procedure assert1,assert2,assert3,assert4,assert_v
    end interface
    interface assert_eq
        module procedure assert_eq2,assert_eq3,assert_eq4,assert_eqn
    end interface
    interface swap
        module procedure swap_i,swap_r,swap_rv,&
            swap_z,swap_zv,swap_zm, &
            masked_swap_rs,masked_swap_rv,masked_swap_rm
    end interface
    interface arth
        module procedure arth_r, arth_i
    end interface
    interface
        subroutine amoeba_anneal(p,y,pb,yb,ftol,func,iter,temptr)
        integer, intent(inout) :: iter
        real, intent(inout) :: yb
        real, intent(in) :: ftol,temptr
        real, dimension(:), intent(inout) :: y,pb
        real, dimension(:,:), intent(inout) :: p
        interface
            function func(x)
            real, dimension(:), intent(in) :: x
            real :: func
            end function func
        end interface
        end subroutine amoeba_anneal
    end interface
    interface indexing
        subroutine indexing_re(arr,index)
        real, dimension(:), intent(in) :: arr
        integer, dimension(:), intent(out) :: index
        end subroutine indexing_re
        subroutine indexing_i4b(iarr,index)
        integer, dimension(:), intent(in) :: iarr
        integer, dimension(:), intent(out) :: index
        end subroutine indexing_i4b
    end interface
    interface
        subroutine sort(arr)
        real, dimension(:), intent(inout) :: arr
        end subroutine sort
    end interface
    interface sort2
        subroutine sort2_rr(arr,slave)
            real, dimension(:), intent(inout) :: arr,slave
        end subroutine sort2_rr
        subroutine sort2_rc(arr,slave,length)
            integer, intent(in) :: length
            real, dimension(:), intent(inout) :: arr
            character(len=length), dimension(:), intent(inout) :: slave
        end subroutine sort2_rc
        subroutine sort2_ri(arr,slave)
            real, dimension(:), intent(inout) :: arr
            integer, dimension(:), intent(inout) :: slave
        end subroutine sort2_ri
        subroutine sort2_cr(arr,slave,length)
            integer, intent(in) :: length
            character(len=length), dimension(:), intent(inout) :: arr
            real, dimension(:), intent(inout) :: slave
        end subroutine sort2_cr
        subroutine sort2_ci(arr,slave,length)
            integer, intent(in) :: length
            character(len=length), dimension(:), intent(inout) :: arr
            integer, dimension(:), intent(inout) :: slave
        end subroutine sort2_ci
        subroutine sort2_cc(arr,slave)
            character, dimension(:), intent(inout) :: arr, slave
        end subroutine sort2_cc
    end interface
CONTAINS
    subroutine tdefit_error(string)
        character(len=*), intent(in) :: string
        print *, string
        call exit(0)
    end subroutine tdefit_error
    subroutine assert1(n1,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1
    if (.not. n1) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        stop 'program terminated by assert1'
    end if
    end subroutine assert1

    subroutine assert2(n1,n2,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1,n2
    if (.not. (n1 .and. n2)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        stop 'program terminated by assert2'
    end if
    end subroutine assert2

    subroutine assert3(n1,n2,n3,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        stop 'program terminated by assert3'
    end if
    end subroutine assert3

    subroutine assert4(n1,n2,n3,n4,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        stop 'program terminated by assert4'
    end if
    end subroutine assert4

    subroutine assert_v(n,string)
    character(len=*), intent(in) :: string
    logical, dimension(:), intent(in) :: n
    if (.not. all(n)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        stop 'program terminated by assert_v'
    end if
    end subroutine assert_v

    function assert_eq2(n1,n2,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq2
    if (n1 == n2) then
        assert_eq2=n1
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        stop 'program terminated by assert_eq2'
    end if
    end function assert_eq2

    function assert_eq3(n1,n2,n3,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2,n3
    integer :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
        assert_eq3=n1
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        stop 'program terminated by assert_eq3'
    end if
    end function assert_eq3

    function assert_eq4(n1,n2,n3,n4,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2,n3,n4
    integer :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
        assert_eq4=n1
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        stop 'program terminated by assert_eq4'
    end if
    end function assert_eq4

    function assert_eqn(nn,string)
    character(len=*), intent(in) :: string
    integer, dimension(:), intent(in) :: nn
    integer :: assert_eqn
    if (all(nn(2:) == nn(1))) then
        assert_eqn=nn(1)
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        stop 'program terminated by assert_eqn'
    end if
    end function assert_eqn

    subroutine swap_i(a,b)
    integer, intent(inout) :: a,b
    integer :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_i

    subroutine swap_r(a,b)
    real, intent(inout) :: a,b
    real :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_r

    subroutine swap_rv(a,b)
    real, dimension(:), intent(inout) :: a,b
    real, dimension(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_rv

    subroutine swap_z(a,b)
    complex, intent(inout) :: a,b
    complex :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_z

    subroutine swap_zv(a,b)
    complex, dimension(:), intent(inout) :: a,b
    complex, dimension(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_zv

    subroutine swap_zm(a,b)
    complex, dimension(:,:), intent(inout) :: a,b
    complex, dimension(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_zm

    subroutine masked_swap_rs(a,b,mask)
    real, intent(inout) :: a,b
    logical, intent(in) :: mask
    real :: swp
    if (mask) then
        swp=a
        a=b
        b=swp
    end if
    end subroutine masked_swap_rs

    subroutine masked_swap_rv(a,b,mask)
    real, dimension(:), intent(inout) :: a,b
    logical, dimension(:), intent(in) :: mask
    real, dimension(size(a)) :: swp
    where (mask)
        swp=a
        a=b
        b=swp
    end where
    end subroutine masked_swap_rv

    subroutine masked_swap_rm(a,b,mask)
    real, dimension(:,:), intent(inout) :: a,b
    logical, dimension(:,:), intent(in) :: mask
    real, dimension(size(a,1),size(a,2)) :: swp
    where (mask)
        swp=a
        a=b
        b=swp
    end where
    end subroutine masked_swap_rm

    function arth_r(first,increment,n)
    real, intent(in) :: first,increment
    integer, intent(in) :: n
    real, dimension(n) :: arth_r
    integer :: k,k2
    real :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth_r(k)=arth_r(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth_r(k)=arth_r(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
    end function arth_r

    function arth_i(first,increment,n)
    integer, intent(in) :: first,increment,n
    integer, dimension(n) :: arth_i
    integer :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth_i(k)=arth_i(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
    end function arth_i
end module tdefit_util
