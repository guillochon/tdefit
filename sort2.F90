SUBROUTINE sort2_rr(arr,slave)
    USE tdefit_util, ONLY : assert_eq, indexing
    real, DIMENSION(:), INTENT(INOUT) :: arr,slave
    integer :: ndum
    integer, DIMENSION(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(arr,index)
    arr=arr(index)
    slave=slave(index)
END SUBROUTINE sort2_rr

SUBROUTINE sort2_rc(arr,slave,length)
    USE tdefit_util, ONLY : assert_eq, indexing
    INTEGER, INTENT(IN) :: length
    real, DIMENSION(:), INTENT(INOUT) :: arr
    CHARACTER(len=length), DIMENSION(:), INTENT(INOUT) :: slave
    integer :: ndum
    integer, DIMENSION(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2_rc')
    call indexing(arr,index)
    arr=arr(index)
    slave=slave(index)
END SUBROUTINE sort2_rc

SUBROUTINE sort2_ri(arr,slave)
    USE tdefit_util, ONLY : assert_eq, indexing
    real, DIMENSION(:), INTENT(INOUT) :: arr
    integer, DIMENSION(:), INTENT(INOUT) :: slave
    integer :: ndum
    integer, DIMENSION(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(arr,index)
    arr=arr(index)
    slave=slave(index)
END SUBROUTINE sort2_ri

SUBROUTINE sort2_cr(arr,slave,length)
    USE tdefit_util, ONLY : assert_eq, indexing
    INTEGER, INTENT(IN) :: length
    CHARACTER(len=length), DIMENSION(:), INTENT(INOUT) :: arr
    real, DIMENSION(:), INTENT(INOUT) :: slave
    integer :: ndum, i, j
    integer, DIMENSION(size(arr)) :: index, iarr
    CHARACTER(len=length) :: temp_c
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
END SUBROUTINE sort2_cr

SUBROUTINE sort2_ci(arr,slave,length)
    USE tdefit_util, ONLY : assert_eq, indexing
    INTEGER, INTENT(IN) :: length
    CHARACTER(len=length), DIMENSION(:), INTENT(INOUT) :: arr
    integer, DIMENSION(:), INTENT(INOUT) :: slave
    integer :: ndum, i, j
    integer, DIMENSION(size(arr)) :: index, iarr
    CHARACTER(len=length) :: temp_c
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
END SUBROUTINE sort2_ci

SUBROUTINE sort2_cc(arr,slave)
    USE tdefit_util, ONLY : assert_eq, indexing
    CHARACTER, DIMENSION(:), INTENT(INOUT) :: arr, slave
    integer :: ndum, i
    integer, DIMENSION(size(arr)) :: index, iarr
    do i = 1, size(arr)
        iarr(i) = ichar(arr(i))
    enddo
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexing(iarr,index)
    arr=arr(index)
    slave=slave(index)
END SUBROUTINE sort2_cc
