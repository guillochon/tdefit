logical function is_normal(a)
    real, intent(in) :: a

    if (a .ne. a) then
        is_normal = .false.
    else
        is_normal = .true.
    endif
end function

logical function is_abnormal(a)
    use tdefit_interface, only: is_normal
    real, intent(in) :: a

    if (is_normal(a)) then
        is_abnormal = .false.
    else
        is_abnormal = .true.
    endif
end function
