!     Array dimensions for synext6

!     nfilpmax = maximum number of filters which can be picked
!     nflxmax = maximum number of input points in each object spectrum
!     nresmax = maximum number of input points in each filter 
!               response curve
!     nsubmax = maximum number of response curve samples 
!     ntaumax = maximum number of optical depths
!     nzmax = maximum number of redshifts

      integer nfilpmax, nflxmax, nresmax, nsubmax, ntaumax, nzmax
      parameter (nfilpmax=100, nflxmax=1, nresmax=2100, &
                 nsubmax=10000, ntaumax=101, nzmax=202)
