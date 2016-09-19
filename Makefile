.PHONY: all info clean version.fpp

#Prefer ifort over gcc
checkifort := $(shell command -v mpiifort 2>/dev/null)
ifdef checkifort
	FC = mpiifort
else
	FC = mpif90
endif

compiler = $(shell basename $(shell ${FC} -show | cut -d " " -f1))

tip=$(shell git log -n 1  --pretty=format:"%H")

ifeq (${compiler},gfortran)
	GFORT = yes
endif
ifeq (${compiler},gfortran44)
	GFORT = yes
endif

ifdef GFORT
	NSO = -fno-strict-overflow -fwrapv
	FCFLAGS = -O3 -m64 -ffree-form -ffree-line-length-none -fimplicit-none -fdefault-double-8 -fdefault-real-8
#	FCFLAGS = -g -O0 -m64 -ffree-form -ffree-line-length-none -fimplicit-none -fdefault-double-8 -fdefault-real-8 \
	-fbacktrace -fbounds-check
else
	NSO =
#	FCFLAGS = -u -r8 -i4 -O3 -fpp -unroll -mieee-fp
#	FCFLAGS = -u -r8 -i4 -O3 -fpp -m64 -xAVX -axAVX #-profile-functions -profile-loops=all #-mieee-fp
	FCFLAGS = -u -r8 -i4 -O3 -fpp -xSSE3 -unroll -mieee-fp
#	FCFLAGS = -u -r8 -i4 -O3 -traceback -fpp -xSSE3 -unroll -mieee-fp
#	FCFLAGS = -u -r8 -i4 -O3 -pg -xSSE -fpp -unroll -mieee-fp
#	FCFLAGS = -u -r8 -i4 -O0 -g -traceback -check all -fpp -warn all #-mieee-fp
endif

SRCS=$(patsubst %.F90, %.o, $(wildcard *.F90))

NOFPE=amoeba_anneal.o

#$(SRCS): FPE:=-fpe0
#$(NOFPE): FPE:=
$(SRCS): FNO:=
$(FNO): FNO:=$(NSO)

.DEFAULT_GOAL=all

%.o: %.F90 tdefit.fpp Makefile
	$(FC) $(FCFLAGS) $(FPE) $(FNO) -c $<

tdefit_print.o: tdefit_data.o

radius.o: tdefit_data.o

tdefit_data.o: constants.o types.o

tdefit_interface.o: tdefit_data.o qxgs.o

likelihood.o: tdefit_interface.o tdefit_data.o dmdt.o bandmag.o set_event.o constants.o

load_user_vars.o: tdefit_interface.o tdefit_data.o

check_options.o: tdefit_data.o

least_sq.o: tdefit_util.o gamma_inc.o

functions.o: constants.o tdefit_interface.o tdefit_data.o

qxgs.o: constants_nswc.o

tdefit_util.o: constants.o

cosmology.o: constants.o quadpack.o tdefit_interface.o

interp_flash_output.o: bisect.o tdefit_interface.o tdefit_data.o

bisect.o: tdefit_data.o

set_trial_vars.o: tdefit_interface.o tdefit_data.o

set_var.o: constants.o tdefit_data.o type_conversion.o tdefit_interface.o

get_var.o: constants.o tdefit_data.o

ftoABmag.o: tdefit_interface.o tdefit_data.o cosmology.o

dmdt.o: tdefit_interface.o tdefit_data.o

sort.o: tdefit_util.o

sort2.o: indexing.o

indexing.o: tdefit_util.o

disk_temp.o: constants.o tdefit_interface.o tdefit_data.o

bandmag.o: tdefit_interface.o dffunc.o disk_temp.o bbflux.o \
		   tdefit_data.o annulus_intercept.o integrate_df.o dmdt.o

integrate_df.o: qxgs.o trapezoid.o quadpack.o

amoeba_anneal.o: tdefit_util.o

tdefit.o: init.o radius.o magdev.o dmdt.o bandmag.o sort2.o set_trial_vars.o \
          init_search_grid.o tdefit_data.o init_search_grid.o load_user_vars.o load_event.o \
		  amoeba_anneal.o check_options.o set_event.o tdefit_memory.o \
          print_trial_vars.o acor.o draw_random_walker.o

draw_random_walker.o: tdefit_data.o tdefit_interface.o

annulus_intercept.o: tdefit_data.o tdefit_interface.o ang_frac.o

print_trial_vars.o: tdefit_data.o get_var.o tdefit_interface.o

ang_frac.o: tdefit_data.o constants.o tdefit_interface.o

load_event.o: constants.o tdefit_data.o tdefit_interface.o

magdev.o: bandmag.o set_trial_vars.o tdefit_data.o set_event.o

tdefit_memory.o: tdefit_data.o

init.o: constants.o tdefit_interface.o alambda.o qxgs.o least_sq.o tdefit_data.o load_defaults.o tdefit_util.o

load_defaults.o: tdefit_interface.o tdefit_data.o

init_search_grid.o: constants.o tdefit_data.o

bbflux.o: tdefit_interface.o quadpack.o bbfunc.o qxgs.o trapezoid.o tdefit_data.o

trapezoid.o: tdefit_interface.o

bbsed.o: tdefit_interface.o tdefit_data.o

bbfunc.o: tdefit_interface.o filterfunc.o alambdaz.o tdefit_data.o

dffunc.o: tdefit_interface.o constants.o disk_temp.o bbflux.o tdefit_data.o bbsed.o ang_frac.o

alambdaz.o: alambda.o tdefit_interface.o tdefit_data.o

alambda.o: constants.o tdefit_interface.o tdefit_data.o bisect.o redlaws.o

redlaws.o: hpsort.o splt.o splt_p.o constants.o tdefit_data.o

hpsort.o: tdefit_util.o

acor.o: tdefit_data.o

filterfunc.o: tdefit_interface.o bisect.o tdefit_data.o

is_normal.o: tdefit_interface.o

set_event.o: tdefit_interface.o

write_vars.o: tdefit_interface.o

all: tdefit

version.fpp:
	@echo '#define VERSION "$(tip)"' > $@

%.mod: %.h

tdefit: version.fpp $(SRCS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod *__genmod* tdefit version.fpp

