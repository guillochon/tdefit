subroutine check_options
    use tdefit_data
    
    if (my_pe .eq. 0) then
        if (.not. dump_burned_walkers .and. dump_burned_ensembles) then
            print *, 'Error: Dump burned walkers must be enabled if dump burned ensembles is enabled'
            call exit(0)
        endif
    endif
end subroutine
