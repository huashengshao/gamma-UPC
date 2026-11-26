# gammaUPC library.
# In MadGraph (default): FC, FFLAGS, libext and CREATELIB from ../../make_opts, the library in ../../../lib/,
# with the MadGraph entry photonpdfsquare (v2.0).
# Standalone (drivers, tests): make STANDALONE=1: gfortran, -O2 unless FFLAGS is set, the library here, no
# photonpdfsquare (it needs MadGraph's pdf.inc).
ifeq ($(STANDALONE),1)
libext=a
define CREATELIB
$(AR) cru $(1) $(2)
ranlib $(1)
endef
FC=gfortran
FFLAGS?=-O2
LIBDIR=./
MGOBJ=
else
LIBDIR=../../../lib/
include ../../make_opts
MGOBJ=photonpdfsquare.o
endif
LIBRARY=libgammaUPC.$(libext)

EXE=test

# make INELASTIC=1: the inelastic photon (structure functions at high Q2) with LHAPDF; otherwise stubs that stop
ifeq ($(INELASTIC),1)
LHAPDFOBJ=
LHAPDFLIBS=$(shell lhapdf-config --libs) -lc++
else
LHAPDFOBJ=lhapdf_stub.o
LHAPDFLIBS=
endif

# optional, never the default: make OPENMP=1 computes the rows of the inelastic photon in parallel (programs then
# link with -fopenmp); without it the code is serial (the OpenMP directives are comments)
ifeq ($(OPENMP),1)
override FFLAGS+=-fopenmp
endif

photonflux=tbessk.o tbessj.o elliptic_integral.o nielsen_generalized_polylog.o photoabsorption.o bessel_integration.o alpha_running.o ElasticFlux.o PhotonPhotonFlux.o OpticalGlauber_Geometry.o nintlib.o interpolation.o fmin_module.o structure_functions.o InelasticFlux.o $(LHAPDFOBJ) $(MGOBJ)

all: $(LIBDIR)$(LIBRARY)

$(LIBDIR)$(LIBRARY): $(photonflux)
	$(call CREATELIB, $@, $^)

$(EXE): test.o $(photonflux)
	$(FC) $(FFLAGS) -o $(EXE) test.o $(photonflux) $(LHAPDFLIBS)

# Each module is a prerequisite through its .o: one 'gfortran -c X.f90' writes both X.o and x.mod. A rule
# 'X.o X.mod: X.f90' declares two independent targets, so 'make -j' ran that compilation twice at the same time and
# the two collided on x.mod. gfortran also lowercases the .mod name, which a case-sensitive file system never
# matches to X.mod. (v2.0)
test.o: test.f90 PhotonPhotonFlux.o OpticalGlauber_Geometry.o
	$(FC) $(FFLAGS) -c $<

PhotonPhotonFlux.o: PhotonPhotonFlux.f90 ElasticFlux.o InelasticFlux.o OpticalGlauber_Geometry.o nintlib.o interpolation.o fmin_module.o photoabsorption.o
	$(FC) $(FFLAGS) -c $<

ElasticFlux.o: ElasticFlux.f90 bessel_integration.o alpha_running.o OpticalGlauber_Geometry.o nintlib.o interpolation.o fmin_module.o
	$(FC) $(FFLAGS) -c $<

bessel_integration.o: bessel_integration.f90
	$(FC) $(FFLAGS) -c $<

structure_functions.o: structure_functions.f90 ElasticFlux.o
	$(FC) $(FFLAGS) -c $<

InelasticFlux.o: InelasticFlux.f90 ElasticFlux.o OpticalGlauber_Geometry.o alpha_running.o structure_functions.o bessel_integration.o interpolation.o
	$(FC) $(FFLAGS) -c $<

lhapdf_stub.o: lhapdf_stub.f90
	$(FC) $(FFLAGS) -c $<

alpha_running.o: alpha_running.f90 nielsen_generalized_polylog.o
	$(FC) $(FFLAGS) -c $<

photoabsorption.o: photoabsorption.f90
	$(FC) $(FFLAGS) -c $<

photonpdfsquare.o: photonpdfsquare.f PhotonPhotonFlux.o InelasticFlux.o
	$(FC) $(FFLAGS) -c $<

tbessk.o: tbessk.f90
	$(FC) $(FFLAGS) -c $<

tbessj.o: tbessj.f90
	$(FC) $(FFLAGS) -c $<

elliptic_integral.o: elliptic_integral.f90
	$(FC) $(FFLAGS) -c $<

nielsen_generalized_polylog.o: nielsen_generalized_polylog.f90
	$(FC) $(FFLAGS) -c $<

interpolation.o: interpolation.f90
	$(FC) $(FFLAGS) -c $<

fmin_module.o: fmin_module.f90
	$(FC) $(FFLAGS) -c $<

OpticalGlauber_Geometry.o: OpticalGlauber_Geometry.f90 nintlib.o interpolation.o nielsen_generalized_polylog.o
	$(FC) $(FFLAGS) -c $<

nintlib.o: nintlib.f90
	$(FC) $(FFLAGS) -c $<

clean_all:
	$(RM) *.o *~ $(BIN)$(EXE) *.mod $(LIBDIR)$(LIBRARY)

clean:
	$(RM) *.o *~ *.mod $(BIN)$(EXE) $(LIBDIR)$(LIBRARY)
