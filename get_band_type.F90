function get_band_type(band) result(bt)
    character*2, intent(in) :: band
    character*1 :: bt

    select case(band)
        case('X1', 'X2', 'Xs')
            bt = 'X'
        case('Lb')
            bt = 'L'
        case('Hl', 'HL', '51')
            bt = 'l'
        case default
            bt = 'O'
    end select
end function
