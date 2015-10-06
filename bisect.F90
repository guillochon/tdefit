function bisect(arr, val, retj) result(i)
    use tdefit_data, only : bisect_lim

    real, intent(in), dimension(:) :: arr
    real, intent(in) :: val
    logical, intent(in), optional :: retj
    integer :: i, j, k

    i=1   
    j=size(arr)
    if (val <= arr(1)) then
        i=1
        return
    elseif (val >= arr(j)) then
        i=j-1
        return
    endif

    !i = count(arr<=val)
    !j = i+1
    !if (j .le. bisect_lim) then
    !    !i = count(arr<=val)
    !    !j = i+1
    !    do k = 2, j
    !        if (val < arr(k)) then
    !            i=k-1
    !            j=k
    !            exit
    !        endif
    !    enddo
    !else
        do 
            k=(i+j)/2   
            if (val < arr(k)) then 
                j=k  
            else
                i=k
            end if
            if (i+1 >= j) exit
        end do
    !endif

    if (present(retj)) then
        if (retj) i = j
    endif
end function
