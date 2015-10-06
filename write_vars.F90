subroutine write_vars(fn)
#include "tdefit.fpp"
#include "version.fpp"
    use tdefit_data, only: all_var_names, all_var_types
    use tdefit_interface, only: get_var

    integer, intent(in) :: fn
    integer :: i
    do i = 1, size(all_var_names)
        write(fn,'(A,1X)',advance='no') trim(all_var_names(i))
    enddo
    write(fn,'(A)') ''
    do i = 1, size(all_var_names)
        write(fn,'(G18.8,1X)',advance='no') get_var(all_var_names(i))
    enddo
    write(fn,'(A)') ''
    do i = 1, size(all_var_names)
        if (all_var_types(i) .ne. -1) then
            write(fn,'(I1,1X)',advance='no') 1
        else
            write(fn,'(I1,1X)',advance='no') 0
        endif
    enddo
    write(fn,'(A)') ''
    do i = 1, size(all_var_names)
        write(fn,'(I2,1X)',advance='no') all_var_types(i)
    enddo
end subroutine
