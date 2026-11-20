# Unit tests of the gamma-UPC numerical utilities (milestone M1).
# Compiles the library sources it needs itself (independent of the main build):
#   make          build and run the tests, writes results_utils.json
#   make ref      regenerate bessel_ref.dat, besselkn_ref.dat (python3 + mpmath)
#   make clean
# FFLAGS=-O0 etc. may be given on the command line.
FC = gfortran
FFLAGS = -O2
SRC = ../..
OUT = build

OBJS = $(OUT)/tbessk.o $(OUT)/nintlib.o $(OUT)/interpolation.o $(OUT)/tbessk_legacy.o

all: run

$(OUT):
	mkdir -p $(OUT)

$(OUT)/tbessk.o: $(SRC)/tbessk.f90 | $(OUT)
	$(FC) $(FFLAGS) -c $< -o $@

$(OUT)/nintlib.o: $(SRC)/nintlib.f90 | $(OUT)
	$(FC) $(FFLAGS) -J$(OUT) -c $< -o $@

$(OUT)/interpolation.o: $(SRC)/interpolation.f90 | $(OUT)
	$(FC) $(FFLAGS) -J$(OUT) -c $< -o $@

$(OUT)/tbessk_legacy.o: tbessk_legacy.f90 | $(OUT)
	$(FC) $(FFLAGS) -c $< -o $@

$(OUT)/test_utils: test_utils.f90 $(OBJS)
	$(FC) $(FFLAGS) -J$(OUT) -o $@ test_utils.f90 $(OBJS)

run: $(OUT)/test_utils bessel_ref.dat besselkn_ref.dat
	./$(OUT)/test_utils

ref:
	python3 gen_refs.py

clean:
	rm -rf $(OUT)

.PHONY: all run ref clean
