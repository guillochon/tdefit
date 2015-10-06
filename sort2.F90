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

subroutine sort2_rr(arr,slave)
    USE tdefit_util, ONLY : assert_eq, indexing
    real, dimension(:), intent(inout) :: arr,slave
    integer :: ndum
    integer, dimension(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(arr,index)
    arr=arr(index)
    slave=slave(index)
end subroutine sort2_rr

subroutine sort2_rc(arr,slave,length)
    USE tdefit_util, ONLY : assert_eq, indexing
    integer, intent(in) :: length
    real, dimension(:), intent(inout) :: arr
    character(len=length), dimension(:), intent(inout) :: slave
    integer :: ndum
    integer, dimension(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2_rc')
    call indexing(arr,index)
    arr=arr(index)
    slave=slave(index)
end subroutine sort2_rc

subroutine sort2_ri(arr,slave)
    USE tdefit_util, ONLY : assert_eq, indexing
    real, dimension(:), intent(inout) :: arr
    integer, dimension(:), intent(inout) :: slave
    integer :: ndum
    integer, dimension(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(arr,index)
    arr=arr(index)
    slave=slave(index)
end subroutine sort2_ri

subroutine sort2_cr(arr,slave,length)
    USE tdefit_util, ONLY : assert_eq, indexing
    integer, intent(in) :: length
    character(len=length), dimension(:), intent(inout) :: arr
    real, dimension(:), intent(inout) :: slave
    integer :: ndum, i, j
    integer, dimension(size(arr)) :: index, iarr
    character(len=length) :: temp_c
    iarr = 0
    do i = 1, size(arr)
        do j = 1, length
            temp_c = arr(i)
            iarr(i) = iarr(i) + 10**(3*(j-1))*iachar(temp_c(j:j))
        enddo
    enddo
    ndum=assert_eq(size(arr),size(slave),'sort2_cr')
    call indexing(iarr,index)
    arr=arr(index)
    slave=slave(index)
end subroutine sort2_cr

subroutine sort2_ci(arr,slave,length)
    USE tdefit_util, ONLY : assert_eq, indexing
    integer, intent(in) :: length
    character(len=length), dimension(:), intent(inout) :: arr
    integer, dimension(:), intent(inout) :: slave
    integer :: ndum, i, j
    integer, dimension(size(arr)) :: index, iarr
    character(len=length) :: temp_c
    iarr = 0
    do i = 1, size(arr)
        do j = 1, length
            temp_c = arr(i)
            iarr(i) = iarr(i) + 10**(3*(j-1))*iachar(temp_c(j:j))
        enddo
    enddo
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(iarr,index)
    arr=arr(index)
    slave=slave(index)
end subroutine sort2_ci

subroutine sort2_cc(arr,slave)
    USE tdefit_util, ONLY : assert_eq, indexing
    character, dimension(:), intent(inout) :: arr, slave
    integer :: ndum, i
    integer, dimension(size(arr)) :: index, iarr
    do i = 1, size(arr)
        iarr(i) = ichar(arr(i))
    enddo
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(iarr,index)
    arr=arr(index)
    slave=slave(index)
end subroutine sort2_cc
