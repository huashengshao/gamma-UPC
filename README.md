# gamma-UPC

A Fortran library for exclusive photon–photon (γγ) processes in ultraperipheral collisions (UPCs) of protons and nuclei, with realistic photon fluxes (version 2.0).

By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN).

**If you use gamma-UPC, please cite** H.-S. Shao and D. d'Enterria, *gamma-UPC: automated generation of exclusive photon–photon processes in ultraperipheral collisions of protons and nuclei with realistic photon fluxes*, [arXiv:2207.03012](https://arxiv.org/abs/2207.03012), JHEP 09 (2022) 248.

gamma-UPC computes the effective γγ luminosities that event generators such as MadGraph5_aMC@NLO and HELAC-Onia fold with partonic cross sections. Its main features:

- **Collisions:** pp, pA/Ap and AB (nucleus–nucleus).
- **Photon fluxes:** the electric-dipole form factor (EDFF), or the charge form factor (ChFF) of the proton (dipole) and of the nuclei (Woods–Saxon).
- **Survival probability:** the probability of no hadronic interaction, from the optical Glauber model (Woods–Saxon or hard sphere), with optional closed-form MC-Glauber fits.
- **Forward-neutron tagging** (0n, 1n…4n, Xn) in AB and pA collisions.
- **Mean impact parameter** ⟨Δb⟩ of the photon emitters.

## Contents

| File | Content |
|---|---|
| `ElasticPhotonPhotonFlux.f90` | module `ElasticPhotonPhotonFlux`: photon fluxes, γγ fluxes and luminosities, ⟨Δb⟩, neutron tagging, the user switches |
| `OpticalGlauber_Geometry.f90` | module `OpticalGlauber_Geometry`: Woods–Saxon nuclei, thickness functions, T̂_AB, σ_inel(NN), MC-Glauber fits |
| `photoabsorption.f90` | photo-absorption and photo-neutron cross sections for neutron tagging |
| `nintlib.f90`, `interpolation.f90`, `fmin_module.f90` | quadrature (Gauss–Legendre, adaptive Gauss–Kronrod), Lagrange interpolation, minimisation |
| `tbessk.f90`, `tbessj.f90`, `elliptic_integral.f90`, `nielsen_generalized_polylog.f90` | special functions |
| `run90.inc` | common blocks for the beams (`ebeam`, `nuclearA_beam1/2`, `nuclearZ_beam1/2`, …), shared with the host program |
| `banner.inc`, `generate_T*hat_WoodsSaxon_grid.inc` | banner and include files of the library |
| `test.f90` | example program |
| `benchmark/` | benchmark and validation drivers (not part of the library build) |

## Building

You need gfortran.

```bash
make
```

This builds the static library `libgammaUPC.a` and the module files in the current directory.

```bash
make test
```

This builds the example program `test` from `test.f90`.

Compiler flags go in `FFLAGS`, for example `make FFLAGS="-O2"` or `make FFLAGS="-O -fPIC"`. The makefile builds without optimisation by default. Link your program against `libgammaUPC.a`, and point the compiler at this directory with `-I` for the `.mod` files and `run90.inc`.

## Quick start

```fortran
PROGRAM example
  USE OpticalGlauber_Geometry
  USE ElasticPhotonPhotonFlux
  IMPLICIT NONE
  include 'run90.inc'
  REAL(KIND(1d0)) :: W, L, flux
  ! Pb+Pb at sqrt(s_NN) = 5.02 TeV: ebeam is the energy per nucleon (GeV)
  ebeam(1) = 2510d0
  ebeam(2) = 2510d0
  nuclearA_beam1 = 208
  nuclearZ_beam1 = 82
  nuclearA_beam2 = 208
  nuclearZ_beam2 = 82
  USE_CHARGEFORMFACTOR4PHOTON = .TRUE.   ! ChFF (default .FALSE.: EDFF)
  ! photon-photon flux d^2N/dx1dx2 at x1 = x2 = 1e-3, x = E_gamma/E_beam
  flux = PhotonPhotonFlux_AB_WoodsSaxon(1d-3, 1d-3)
  ! luminosity dN/dW (GeV^-1) at W = 10 GeV: icoll = 3 (AB), iprofile = 1 (Woods-Saxon)
  L = dLgammagammadW_UPC(10d0, 3, 1)
  PRINT *, flux, L
END PROGRAM example
```

```bash
gfortran -I. example.f90 libgammaUPC.a -o example
```

### Conventions

- `ebeam(i)` is the energy per nucleon of beam i in GeV; for a proton beam it is the proton energy.
- A proton beam has `nuclearA_beam = 0` and `nuclearZ_beam = 0`. In pA either beam may be the proton.
- x = E_γ / E_beam. The `PhotonPhotonFlux_*` functions return d²N_γγ/dx₁dx₂.
- `Lgammagamma_UPC(W, icoll, iprofile)` returns dN_γγ/dτ at τ = W²/s.
- `dLgammagammadW_UPC(W, icoll, iprofile)` returns dN_γγ/dW in GeV⁻¹.
- `DeltaB_UPC_at_W(W, icoll, iprofile)` returns ⟨Δb⟩ in fm.
- `icoll`: 1 = pp, 2 = pA/Ap, 3 = AB.
- `iprofile`: 0 = no hadronic suppression (P_NOHAD = 1), 1 = Woods–Saxon, 2 = hard sphere.
- The collision-specific functions follow the same pattern: `PhotonPhotonFlux_{pp,pA_WoodsSaxon,pA_hardsphere,AB_WoodsSaxon,AB_hardsphere}(x1, x2 [, FORCEPNOHAD1])`. The plain names read from a grid built on the first call; the `_eval` variants evaluate directly. The `PhotonPhotonDeltaB_*` functions give the ⟨Δb⟩ numerator.
- Woods–Saxon parameters are built in for 42 nuclei, from H2 to Pb208 (`GetNuclearInfo`). Neutron tagging needs photo-absorption data and is available for Pb208, Au197 and O16 beams.

## Options

Set these public module variables before the first call. The grids are built on the first call; the new engine rebuilds them when its configuration changes.

| Variable | Default | Meaning |
|---|---|---|
| `USE_CHARGEFORMFACTOR4PHOTON` | `.FALSE.` | photon flux from the charge form factor (ChFF) instead of the electric-dipole one (EDFF) |
| `neutron_tagging(2)` | `(/-2,-2/)` | forward-neutron tag of each beam: −2 any (0n+Xn), −1 Xn, 0 0n, 1–4 for 1n–4n; ignored for a proton beam |
| `neutron_conj_sum` | `.FALSE.` | in AB with different tags, also add the conjugate tag combination |
| `use_MC_Glauber` | `.FALSE.` | survival probability from closed-form MC-Glauber fits (e.g. Pb–Pb) instead of the optical Glauber model |
| `IEPSILON_EDFF` | `0` | EDFF lower cut in b: 0 → R_A, 1 → R_A + a |
| `proton_FF_correction` | `.FALSE.` | with ChFF, corrections beyond the proton dipole form factor (pp and pA) |
| `GENERATE_PhotonPhotonFlux_GRID` | `.TRUE.` | tabulate the γγ flux on an (x₁, x₂) grid (otherwise evaluate directly) |
| `GENERATE_DELTAB_GRID` | `.FALSE.` | tabulate the ⟨Δb⟩ numerator as well |
| `alphaem_elasticphoton` | 1/137 | α_em of the photon fluxes |
| `NEUTRON_XSIGMA` | `0` | shift the photo-neutron cross sections by this many standard deviations |
| `PhotonPhotonFlux_Algorithm` | `2` | 2: kernel-factorised engine (below); 1: legacy integration |
| `FFE_Legacy_ChFF_Definitions` | `.FALSE.` | engine only: use the legacy ChFF photon density (transverse only below 3R, point-like beyond and for x < 10⁻⁴, b ≥ 0.1R) |
| `FFE_NodeScale` | `1` | engine only: >1 makes the quadrature denser, for convergence tests |
| `FFE_Verbose` | `0` | engine only: >0 prints stage timings and diagnostics |
| `FFE_Tag_Legacy` | `.FALSE.` | engine only: fill the tagging table from the legacy `Pbreak_*` functions |
| `FFE_Grid_Cache` | `.FALSE.` | engine only: share the final grids between processes through files (below) |
| `FFE_Grid_Cache_Dir` | `'./'` | engine only: directory of the grid cache files |

## Grid generation: the kernel-factorised engine

The pp, pA/Ap and AB (Woods–Saxon) fluxes and ⟨Δb⟩ are computed by default with a kernel-factorised engine (`PhotonPhotonFlux_Algorithm = 2`). The hard-sphere variants, and pp and pA with `proton_FF_correction`, keep the legacy integration.

**How it works.** The γγ flux is

  d²N/dx₁dx₂ = (2π α² Z₁² Z₂² / x₁x₂) ∫∫ d ln b₁ d ln b₂ W₁(x₁; b₁) W₂(x₂; b₂) Φ(b₁, b₂),

with W = b²N and Φ = ∫dφ P(b₁₂). The kernel Φ does not depend on x. It is integrated once per beam setup into a matrix K over Lagrange bases of fixed ln b grids. The whole grid is then Z = W₁ K W₂ᵀ, and any point off the grid costs about a millisecond. The quadrature follows the edges of the survival probability and has no truncation in b; far from the nuclei Φ comes from a moment expansion. The ChFF photon density is computed exactly in position space. The grid is stored as ln(x₁x₂Z) and read with 8 × 8 Lagrange interpolation.

**Speed and precision** (CPU time of the first call, gfortran -O2):

| Grid | legacy | engine |
|---|---|---|
| Pb–Pb 5.02 TeV, ChFF | 528 s | 1.7 s |
| Pb–Pb 5.02 TeV, EDFF | 51 s | 0.43 s |
| Pb–Pb ChFF with neutron tagging | 600–990 s | 4 s |
| p–Pb 8.16 TeV, ChFF | 457 s | 1.2 s |
| p–p 13 TeV, ChFF | 16 s | 0.44 s |

Against an independent brute-force reference, the engine agrees to ≤ 10⁻⁷ for 32 configurations, pp included. The legacy integration deviates from the same reference by 0.1% to 22% (about 0.1% for pp). Grid interpolation adds ≤ 4 × 10⁻⁷.

**Grid cache.** With `FFE_Grid_Cache = .TRUE.`, the first process that needs a grid writes it to `FFE_Grid_Cache_Dir`, as `gammaUPC_grid_<key>_k<kind>.bin` (about 0.6 MB). The key is a hash of the whole configuration, and other processes with the same configuration read the file instead of rebuilding: under a millisecond instead of 0.3–6 s. While one process builds, the others wait for its file; a lock directory left by a killed process is ignored after 300 s. Files are written under a temporary name and renamed, checked against the configuration and a checksum when read, and rebuilt when they do not match. MadGraph turns the cache on for its jobs. Delete the files after changing the gamma-UPC source code.

Set `PhotonPhotonFlux_Algorithm = 1` to reproduce the legacy numbers. That path includes three bug fixes: the ⟨Δb⟩ grids, pA ChFF with the nucleus as beam 1, and the b = 0 row of the ChFF density grids.

## Benchmarks and validation (`benchmark/`)

These programs use only the public API and link against `libgammaUPC.a`.

- `bench_ab.F90` reads a namelist `/config/` (beams, flags, points file) and writes timings, flux values at chosen (x₁, x₂) points, dL/dW and ⟨Δb⟩(W) as JSON. Build it with `make -C benchmark LIBDIR=.. OUT=<dir> NEW=1`; `NEW=1` enables the engine-only switches.
- `ref_ab.F90` is the independent reference: nested adaptive Gauss–Kronrod in (ln b₁, ln b₂, φ), with the ChFF density from the k-space integrals.
- `bench_blocks.F90` and `test_ffe.F90` check the building blocks of the engine.
- `unit/` holds unit tests of the quadrature, interpolation and Bessel routines. Run them with `make -C benchmark/unit`; the mpmath reference values are included.
