real function logic2dbl(a)
    logical, intent(in) :: a
  
    if (a .eqv. .true.) then
        logic2dbl = 1.d0
    else
        logic2dbl = 0.d0
    end if
end function logic2dbl

logical function int2logic(a)
    integer, intent(in) :: a

    if (a .eq. 1) then
        int2logic = .true.
    else
        int2logic = .false.
    endif
end function

logical function dbl2logic(a)
    real, intent(in) :: a
    integer :: tempi

    tempi = nint(a)


    if (tempi .eq. 1) then
        dbl2logic = .true.
    else
        dbl2logic = .false.
    endif
 end function
