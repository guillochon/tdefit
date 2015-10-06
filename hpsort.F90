subroutine hpsort (x, n, indx, xasc)
      
    use tdefit_util, only: indexing

!   Heapsort

    integer n
    integer i, indx(n)
    real*8 x(n), xasc(n)

    call indexing (x, indx) 
    do i=1,n
       xasc(i) = x(indx(i))
    end do
    return
end
