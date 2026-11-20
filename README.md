# gamma-UPC

A Fortran library for exclusive photon–photon (γγ) processes in ultraperipheral collisions (UPCs) of protons and nuclei, with realistic photon fluxes (version 1.8).

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
| `nintlib.f90`, `interpolation.f90`, `fmin_module.f90` | multidimensional quadrature (`ROMBERG_ND` and others), polynomial and spline interpolation, minimisation |
| `tbessk.f90`, `tbessj.f90`, `elliptic_integral.f90`, `nielsen_generalized_polylog.f90` | special functions |
| `run90.inc` | common blocks for the beams (`ebeam`, `nuclearA_beam1/2`, `nuclearZ_beam1/2`, …), shared with the host program |
| `banner.inc`, `generate_T*hat_WoodsSaxon_grid.inc` | banner and include files of the library |
| `test.f90` | example program |

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

Set these public module variables before the first call: the grids are built on the first call and kept for the rest of the run.

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

## Grid generation

With `GENERATE_PhotonPhotonFlux_GRID = .TRUE.` (the default), the first call builds the γγ flux on an 81 × 81 grid in (log₁₀ x₁, log₁₀ x₂). For Pb–Pb at 5.02 TeV this takes about 50 s with EDFF and about 9 min with ChFF (CPU, gfortran -O2), and about twice as long without optimisation. A kernel-factorised engine that builds the AB and pA grids tens to hundreds of times faster (×310 for Pb–Pb ChFF), and more precisely, is under review in [pull request #1](https://github.com/huashengshao/gamma-UPC/pull/1) (branch `faster-gammaUPC`).

## Known issues

These are fixed on `faster-gammaUPC` and will reach `main` with pull request #1.

- With `GENERATE_DELTAB_GRID = .TRUE.`, the ⟨Δb⟩ grids (AB, pA, pp) return an uninitialised value and index out of bounds. The default `.FALSE.` is not affected.
- pA with ChFF in the A+p beam order (nucleus as beam 1) applies the proton and nucleus photon densities to the wrong sides.
- The b = 0 row of the ChFF photon-density grids is left uninitialised. It is used only below b = 0.005 R.
