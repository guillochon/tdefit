#TDEFit - Software for fitting tidal disruption event light curves
Originally coded by James Guillochon (http://astrocrash.net).

##Installation instructions

To install, clone the repository and then compile the code using make. Code only requires a Fortran compiler that supports Fortran 95, and has been tested with both ifort and gfortran. Code requires two sets of data as inputs that are not included in this repository because of their size or proprietary nature: A directory of event data (available from a separate repository) and a directory of dm/dt data.

```
#!csh

hg clone ssh://hg@bitbucket.org/Guillochon/tdefit
cd tdefit
make -j tdefit
wget http://astrocrash.net/files/tdefit-dmdts.tar.gz
tar -xzf tdefit-dmdts.tar.gz
hg clone ssh://hg@bitbucket.org/Guillochon/TDE_events
```

##Getting started with an example

Within the TDEFit folder is an `example` folder, which will run a fit to ASASSN-14li. To run this example, change to the example directory and run TDEFit.

```
cd example
mpirun -np 4 ../tdefit
```

The first time TDEFit runs, it will spend a few minutes generating binary versions of the ASCII dm/dt files, once those files are created it will read from them for subsequent runs (so long as the `paths.par` file points to the same binary data location), which is much faster.

##Credits

Software includes a few pieces of code originally written by others:

* The [York Extinction Solver](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/staging/proc/tmp/www/YorkExtinctionSolver/fortran/), [McCall 2004](http://adsabs.harvard.edu/abs/2004AJ....128.2144M), which provides functions for a number of common reddening laws, GPL.

* [Quadpack](https://en.wikipedia.org/wiki/QUADPACK), for numerical integration, in the public domain.

* Some functions from [PROB](https://people.sc.fsu.edu/~jburkardt/f_src/prob/prob.html), a Fortran90 probability library written by John Burkardt, LGPL.

* [qxgs](http://jblevins.org/mirror/amiller/), a one-dimensional integrator written by Alan Miller, public domain.