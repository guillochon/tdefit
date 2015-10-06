recursive subroutine trapezoid(func,minx,maxx,div,val)
#include "tdefit.fpp"
    !use, intrinsic :: ieee_arithmetic
    use tdefit_interface, only: is_abnormal

    real, external    :: func
    real, intent(IN)  :: minx, maxx
    integer, intent(IN)           :: div
    real, intent(OUT) :: val

    integer :: i
    real :: coeff, temporary

    val = 0.d0
    coeff = (maxx - minx)/(dble(div) - 1.d0)
    do i = 1, div
        temporary = func(minx + coeff*(dble(i)-1.d0))
#ifdef DEBUG
        !if (.not. ieee_is_normal(temporary)) then
        if (is_abnormal(temporary)) then
            print *, 'val abnormal', coeff, maxx, minx, div, temporary, i
            call exit(0)
        endif
#endif
        if (i .eq. 1 .or. i .eq. div) then
            val = val + 0.5d0*temporary
        else
            val = val + temporary
        endif
    enddo
    val = val * coeff
end subroutine
