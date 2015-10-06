TDEFit -- Originally coded by James Guillochon (http://astrocrash.net)

To install, clone the repository and then compile the code using make. Code only requires a Fortran compiler that supports Fortran 95, and has been tested with both ifort and gfortran.

```
#!csh

hg clone ssh://hg@bitbucket.org/Guillochon/tdefit
make -j tdefit
```

Code requires two files sets of data not included in this repository: A directory of event data and a directory of dm/dt data. Both are available from astrocrash.net.