TDEFit -- Originally coded by James Guillochon (http://astrocrash.net)

Software includes a few pieces of code originally written by others:

* The [York Extinction Solver](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/staging/proc/tmp/www/YorkExtinctionSolver/fortran/), [McCall 2004](http://adsabs.harvard.edu/abs/2004AJ....128.2144M), which provides functions for a number of common reddening laws.

* [Quadpack](https://en.wikipedia.org/wiki/QUADPACK), for numerical integration, in the public domain.

To install, clone the repository and then compile the code using make. Code only requires a Fortran compiler that supports Fortran 95, and has been tested with both ifort and gfortran.

```
#!csh

hg clone ssh://hg@bitbucket.org/Guillochon/tdefit
make -j tdefit
```

Code requires two files sets of data not included in this repository: A directory of event data and a directory of dm/dt data. Both are available from astrocrash.net.