MODULE tdefit_util
    use constants
    integer, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    INTERFACE
        SUBROUTINE least_sq(x,y,a,b,siga,sigb,chi2,q,sig)
        real, DIMENSION(:), INTENT(IN) :: x,y
        real, INTENT(OUT) :: a, b, siga, sigb, chi2, q
        real, DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
        END SUBROUTINE least_sq
    END INTERFACE
    INTERFACE gamma_inc
        FUNCTION gamma_inc(p, x)
            real ( kind = 8 ) p
            real ( kind = 8 ) x
            real ( kind = 8 ) gamma_inc
        END FUNCTION gamma_inc
    END INTERFACE
    INTERFACE assert
        MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
    END INTERFACE
    INTERFACE assert_eq
        MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
    END INTERFACE
    INTERFACE swap
        MODULE PROCEDURE swap_i,swap_r,swap_rv,&
            swap_z,swap_zv,swap_zm, &
            masked_swap_rs,masked_swap_rv,masked_swap_rm
    END INTERFACE
    INTERFACE arth
        MODULE PROCEDURE arth_r, arth_i
    END INTERFACE
    INTERFACE
        SUBROUTINE amoeba_anneal(p,y,pb,yb,ftol,func,iter,temptr)
        integer, INTENT(INOUT) :: iter
        real, INTENT(INOUT) :: yb
        real, INTENT(IN) :: ftol,temptr
        real, DIMENSION(:), INTENT(INOUT) :: y,pb
        real, DIMENSION(:,:), INTENT(INOUT) :: p
        INTERFACE
            FUNCTION func(x)
            real, DIMENSION(:), INTENT(IN) :: x
            real :: func
            END FUNCTION func
        END INTERFACE
        END SUBROUTINE amoeba_anneal
    END INTERFACE
    INTERFACE indexing
        SUBROUTINE indexing_re(arr,index)
        real, DIMENSION(:), INTENT(IN) :: arr
        integer, DIMENSION(:), INTENT(OUT) :: index
        END SUBROUTINE indexing_re
        SUBROUTINE indexing_i4b(iarr,index)
        integer, DIMENSION(:), INTENT(IN) :: iarr
        integer, DIMENSION(:), INTENT(OUT) :: index
        END SUBROUTINE indexing_i4b
    END INTERFACE
    INTERFACE
        SUBROUTINE sort(arr)
        real, DIMENSION(:), INTENT(INOUT) :: arr
        END SUBROUTINE sort
    END INTERFACE
    INTERFACE sort2
        SUBROUTINE sort2_rr(arr,slave)
            real, DIMENSION(:), INTENT(INOUT) :: arr,slave
        END SUBROUTINE sort2_rr
        SUBROUTINE sort2_rc(arr,slave,length)
            INTEGER, INTENT(IN) :: length
            real, DIMENSION(:), INTENT(INOUT) :: arr
            CHARACTER(len=length), DIMENSION(:), INTENT(INOUT) :: slave
        END SUBROUTINE sort2_rc
        SUBROUTINE sort2_ri(arr,slave)
            real, DIMENSION(:), INTENT(INOUT) :: arr
            integer, DIMENSION(:), INTENT(INOUT) :: slave
        END SUBROUTINE sort2_ri
        SUBROUTINE sort2_cr(arr,slave,length)
            INTEGER, INTENT(IN) :: length
            CHARACTER(len=length), DIMENSION(:), INTENT(INOUT) :: arr
            real, DIMENSION(:), INTENT(INOUT) :: slave
        END SUBROUTINE sort2_cr
        SUBROUTINE sort2_ci(arr,slave,length)
            INTEGER, INTENT(IN) :: length
            CHARACTER(len=length), DIMENSION(:), INTENT(INOUT) :: arr
            integer, DIMENSION(:), INTENT(INOUT) :: slave
        END SUBROUTINE sort2_ci
        SUBROUTINE sort2_cc(arr,slave)
            CHARACTER, DIMENSION(:), INTENT(INOUT) :: arr, slave
        END SUBROUTINE sort2_cc
    END INTERFACE
CONTAINS
    subroutine tdefit_error(string)
        character(len=*), intent(in) :: string
        print *, string
        call exit(0)
    end subroutine tdefit_error
    SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        STOP 'program terminated by assert1'
    end if
    END SUBROUTINE assert1

    SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        STOP 'program terminated by assert2'
    end if
    END SUBROUTINE assert2

    SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        STOP 'program terminated by assert3'
    end if
    END SUBROUTINE assert3

    SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        STOP 'program terminated by assert4'
    end if
    END SUBROUTINE assert4

    SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
        write (*,*) 'error: an assertion failed with this tag:', &
            string
        STOP 'program terminated by assert_v'
    end if
    END SUBROUTINE assert_v

    FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
        assert_eq2=n1
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq2'
    end if
    END FUNCTION assert_eq2

    FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
        assert_eq3=n1
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq3'
    end if
    END FUNCTION assert_eq3

    FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
        assert_eq4=n1
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq4'
    end if
    END FUNCTION assert_eq4

    FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
        assert_eqn=nn(1)
    else
        write (*,*) 'error: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eqn'
    end if
    END FUNCTION assert_eqn

    SUBROUTINE swap_i(a,b)
    integer, INTENT(INOUT) :: a,b
    integer :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_i

    SUBROUTINE swap_r(a,b)
    real, INTENT(INOUT) :: a,b
    real :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_r

    SUBROUTINE swap_rv(a,b)
    real, DIMENSION(:), INTENT(INOUT) :: a,b
    real, DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_rv

    SUBROUTINE swap_z(a,b)
    complex, INTENT(INOUT) :: a,b
    complex :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_z

    SUBROUTINE swap_zv(a,b)
    complex, DIMENSION(:), INTENT(INOUT) :: a,b
    complex, DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_zv

    SUBROUTINE swap_zm(a,b)
    complex, DIMENSION(:,:), INTENT(INOUT) :: a,b
    complex, DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_zm

    SUBROUTINE masked_swap_rs(a,b,mask)
    real, INTENT(INOUT) :: a,b
    logical, INTENT(IN) :: mask
    real :: swp
    if (mask) then
        swp=a
        a=b
        b=swp
    end if
    END SUBROUTINE masked_swap_rs

    SUBROUTINE masked_swap_rv(a,b,mask)
    real, DIMENSION(:), INTENT(INOUT) :: a,b
    logical, DIMENSION(:), INTENT(IN) :: mask
    real, DIMENSION(size(a)) :: swp
    where (mask)
        swp=a
        a=b
        b=swp
    end where
    END SUBROUTINE masked_swap_rv

    SUBROUTINE masked_swap_rm(a,b,mask)
    real, DIMENSION(:,:), INTENT(INOUT) :: a,b
    logical, DIMENSION(:,:), INTENT(IN) :: mask
    real, DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
        swp=a
        a=b
        b=swp
    end where
    END SUBROUTINE masked_swap_rm

    FUNCTION arth_r(first,increment,n)
    real, INTENT(IN) :: first,increment
    integer, INTENT(IN) :: n
    real, DIMENSION(n) :: arth_r
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
    END FUNCTION arth_r

    FUNCTION arth_i(first,increment,n)
    integer, INTENT(IN) :: first,increment,n
    integer, DIMENSION(n) :: arth_i
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
    END FUNCTION arth_i
END MODULE tdefit_util
