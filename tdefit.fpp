!#define DEBUG
!#define VERBOSE_MPI
#define REJECT_OUT_OF_BOUNDS
#define REJECT_INITIAL_FAILS
#define MPIH
!#define PRINT_PENALTY_REASONS
#define PRINT_ANNEALING_INFO
!#define BLR_EXCLUDE
#define VERBOSE_INIT

!common variables
#define TIME_VAR         1

!multipoly.dat variables
#define MDAT_NCOLS       22
#define MTOT_VAR         2
#define KTOT_VAR         7
#define EITOT_VAR        8
#define GRTOT_VAR        9
#define JTOTX_VAR        10
#define JTOTY_VAR        11
#define JTOTZ_VAR        12
#define KOBJ_VAR         13
#define EIOBJ_VAR        14
#define GROBJ_VAR        15
#define JOBJX_VAR        16
#define JOBJY_VAR        17
#define JOBJZ_VAR        18
#define MACC_VAR         19
#define JACCX_VAR        20
#define JACCY_VAR        21
#define JACCZ_VAR        22

!orbit.dat variables
#define ODAT_NCOLS       51
#define VEC1X_VAR        2
#define VEC1Y_VAR        3
#define VEC1Z_VAR        4
#define VEC1VX_VAR       5
#define VEC1VY_VAR       6
#define VEC1VZ_VAR       7
#define VEC2X_VAR        8
#define VEC2Y_VAR        9
#define VEC2Z_VAR        10
#define VEC2VX_VAR       11
#define VEC2VY_VAR       12
#define VEC2VZ_VAR       13
#define BNDCOMX_VAR      14
#define BNDCOMY_VAR      15
#define BNDCOMZ_VAR      16
#define BNDCOMVX_VAR     17
#define BNDCOMVY_VAR     18
#define BNDCOMVZ_VAR     19
#define TOTCOMX_VAR      20
#define TOTCOMY_VAR      21
#define TOTCOMZ_VAR      22
#define TOTCOMVX_VAR     23
#define TOTCOMVY_VAR     24
#define TOTCOMVZ_VAR     25
#define PEAKCOMX_VAR     26
#define PEAKCOMY_VAR     27
#define PEAKCOMZ_VAR     28
#define PEAKCOMVX_VAR    29
#define PEAKCOMVY_VAR    30
#define PEAKCOMVZ_VAR    31
#define MPOLECOMX_VAR    32
#define MPOLECOMY_VAR    33
#define MPOLECOMZ_VAR    34
#define MPOLECOMVX_VAR   35
#define MPOLECOMVY_VAR   36
#define MPOLECOMVZ_VAR   37
#define RAD_VAR          38
#define VEL_VAR          39
#define SEMIMAJ_VAR      41
#define ECC_VAR          42
#define NU_VAR           44
#define MBND_VAR         49
#define GRENER_VAR       50
#define ENER_VAR         51

#define MDAT_ARR         1
#define ODAT_ARR         2
#define DDAT_ARR         3
#define IDDAT_ARR        4
#define RHODDAT_ARR      5
#define D_EMIN_ARR       6
#define D_EMAX_ARR       7
#define MFINAL_ARR       8
#define DELE_ARR         9
#define PT_ARR           10

#define NUM_EDAT         9
#define E_BETA           1
#define E_PTIME          3
#define E_ROBJ           4
#define E_MOBJ           5
#define E_ECC            6
#define E_MPERT          7

#define NUM_REPMODELS    6
#define RM_NONE          1
#define RM_CLUMPY_INNER  2
#define RM_CLUMPY_OUTER  3
#define RM_ANNULUS       4
#define RM_REFLECT       5
#define RM_CLOUDS        6
#define RM_CIRCULAR      7

#define NUM_IWDISTS      4
#define IW_RANDOM        1
#define IW_STATISTICAL   2
#define IW_PROBABLE      3
#define IW_REGULAR       4
