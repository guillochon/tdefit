subroutine tdefit_print(message)
    use tdefit_data

    character*(*), intent(in) :: message

    if (print_children) then
        write(*,'(A12,I2,A2,A)'), '[Processor #', my_pe, '] ', trim(message)
    elseif (my_pe .eq. 0) then
        write(*,'(A)'), trim(message)
    endif
end subroutine
