module types
    type bandnode
        type(bandnode), pointer :: next => null()
        character*2 :: band
    end type bandnode
    type bnptr
        type(bandnode), pointer :: p
    end type bnptr
end module types
