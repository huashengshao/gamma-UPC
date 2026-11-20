# Benchmark drivers for gamma-UPC (not part of the library build).
#   make LIBDIR=<dir with libgammaUPC.a, *.mod, run90.inc> OUT=<build dir> [FFLAGS=-O2] [NEW=1]
LIBDIR ?= ..
OUT ?= .
FC = gfortran
FFLAGS ?= -O2
ifeq ($(NEW),1)
CPPDEFS = -DGUPC_NEW
endif

all: $(OUT)/bench_ab

$(OUT)/bench_ab: bench_ab.F90 $(LIBDIR)/libgammaUPC.a
	mkdir -p $(OUT)
	$(FC) $(FFLAGS) $(CPPDEFS) -I$(LIBDIR) -J$(OUT) -o $@ bench_ab.F90 $(LIBDIR)/libgammaUPC.a

clean:
	rm -f $(OUT)/bench_ab $(OUT)/*.mod

.PHONY: all clean
