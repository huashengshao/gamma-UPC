
libext=a
define CREATELIB
$(AR) cru $(1) $(2)
ranlib $(1)
endef

LIBRARY       = libgammaUPC.$(libext)
LIBDIR        = ./

FC=gfortran
# default optimisation; FFLAGS set in the environment or on the command line (e.g. by MadGraph) takes precedence
FFLAGS?=-O2

EXE=test

# make INELASTIC=1: the inelastic photon (structure functions at high Q2) with LHAPDF; otherwise stubs that stop
ifeq ($(INELASTIC),1)
LHAPDFOBJ=
LHAPDFLIBS=$(shell lhapdf-config --libs) -lc++
else
LHAPDFOBJ=lhapdf_stub.o
LHAPDFLIBS=
endif

# optional, never the default: make OPENMP=1 computes the rows of the inelastic photon (Inelastic_fast_building) in
# parallel (programs then link with -fopenmp); without it the code is serial (the OpenMP directives are comments)
ifeq ($(OPENMP),1)
override FFLAGS+=-fopenmp
endif

photonflux=tbessk.o tbessj.o elliptic_integral.o nielsen_generalized_polylog.o photoabsorption.o bessel_integration.o alpha_running.o ElasticFlux.o PhotonPhotonFlux.o ElasticPhotonPhotonFlux.o OpticalGlauber_Geometry.o nintlib.o interpolation.o fmin_module.o structure_functions.o InelasticFlux.o $(LHAPDFOBJ)

all: $(LIBDIR)$(LIBRARY)

$(LIBDIR)$(LIBRARY): $(photonflux)
	$(call CREATELIB, $@, $^)

$(EXE): test.o $(photonflux)
	$(FC) $(FFLAGS) -o $(EXE) test.o $(photonflux) $(LHAPDFLIBS)

test.o: test.f90 ElasticPhotonPhotonFlux.mod OpticalGlauber_Geometry.mod
	$(FC) $(FFLAGS) -c $<

ElasticPhotonPhotonFlux.o ElasticPhotonPhotonFlux.mod: ElasticPhotonPhotonFlux.f90 PhotonPhotonFlux.mod
	$(FC) $(FFLAGS) -c $<

PhotonPhotonFlux.o PhotonPhotonFlux.mod: PhotonPhotonFlux.f90 ElasticFlux.mod inelasticflux.mod OpticalGlauber_Geometry.mod nintlib.mod interpolation.mod fmin_module.mod photoabsorption.mod
	$(FC) $(FFLAGS) -c $<

ElasticFlux.o ElasticFlux.mod: ElasticFlux.f90 proton_A1_formfactors.inc alpha_running_poles.inc bessel_integration.mod alpha_running.mod OpticalGlauber_Geometry.mod nintlib.mod interpolation.mod fmin_module.mod
	$(FC) $(FFLAGS) -c $<

bessel_integration.o bessel_integration.mod: bessel_integration.f90
	$(FC) $(FFLAGS) -c $<

structure_functions.o structure_functions.mod: structure_functions.f90 clas_spline.inc interpolation.mod ElasticFlux.mod
	$(FC) $(FFLAGS) -c $<

InelasticFlux.o inelasticflux.mod: InelasticFlux.f90 ElasticFlux.mod alpha_running.mod structure_functions.mod bessel_integration.mod interpolation.mod
	$(FC) $(FFLAGS) -c $<

lhapdf_stub.o: lhapdf_stub.f90
	$(FC) $(FFLAGS) -c $<

alpha_running.o alpha_running.mod: alpha_running.f90 nielsen_generalized_polylog.mod
	$(FC) $(FFLAGS) -c $<

photoabsorption.o photoabsorption.mod: photoabsorption.f90
	$(FC) $(FFLAGS) -c $<

tbessk.o: tbessk.f90
	$(FC) $(FFLAGS) -c $<

tbessj.o: tbessj.f90
	$(FC) $(FFLAGS) -c $<

elliptic_integral.o: elliptic_integral.f90
	$(FC) $(FFLAGS) -c $<

nielsen_generalized_polylog.o nielsen_generalized_polylog.mod: nielsen_generalized_polylog.f90
	$(FC) $(FFLAGS) -c $<

interpolation.o interpolation.mod: interpolation.f90
	$(FC) $(FFLAGS) -c $<

fmin_module.o fmin_module.mod: fmin_module.f90
	$(FC) $(FFLAGS) -c $<

OpticalGlauber_Geometry.o OpticalGlauber_Geometry.mod: OpticalGlauber_Geometry.f90 nintlib.mod interpolation.mod nielsen_generalized_polylog.mod
	$(FC) $(FFLAGS) -c $<

nintlib.o nintlib.mod: nintlib.f90
	$(FC) $(FFLAGS) -c $<

clean_all:
	$(RM) *.o *~ $(BIN)$(EXE) *.mod $(LIBDIR)$(LIBRARY)

clean:
	$(RM) *.o *~ *.mod $(BIN)$(EXE) $(LIBDIR)$(LIBRARY)
