MODULE PhotonPhotonFlux
  USE ElasticFlux
  USE OpticalGlauber_Geometry
  USE NINTLIB ! for multiple dimensional integrations
  USE interpolation
  USE InelasticFlux,ONLY:PhotonNumberDensity_inelastic_proton,PhotonFlux_proton_inelastic_nob,inel_fast_prepare_rows,&
       inel_rhohat,inel_rhohat_cached,elastic_rhohat,inel_rhohat_key,nucleus_rhohat,nucleus_rhohat_cached,&
       elastic_rhohat_cached,vec_cache_get,vec_cache_put
  IMPLICIT NONE
  PRIVATE
  PUBLIC::alphaem_elasticphoton,alpha_running4photon,GENERATE_PROTON_FLUX_GRID
  ! For UPCs photon-photon collisions
  LOGICAL,PUBLIC::USE_CHARGEFORMFACTOR4PHOTON=.FALSE.
  ! when USE_CHARGEFORMFACTOR4PHOTON=.TRUE., the proton photon density without the small-x approximation
  ! (PhotonNumberDensity_ChargeFormFactor_proton_exact) when USE_EXACT_ELASTIC4PHOTON=.TRUE.
  LOGICAL,PUBLIC::USE_EXACT_ELASTIC4PHOTON=.FALSE.
  ! final state of the protons (pp, and pA with the Woods-Saxon profile), P = P_NOHAD: 0 el-el P (elastic only,
  ! default); 1 (el-inel + inel-el) P; 2 (el-el + el-inel + inel-el) (1-P) + inel-inel; 3 all without survival factor
  ! (0+1+2 = 3). In pA only the proton has an inelastic photon (the nucleus coherent only): 1 inel-el P, 2 (el-el +
  ! inel-el) (1-P), 3 (el + inel) el. The inelastic photon (InelasticFlux) at the scale mu = W_gammagamma
  ! (mu2 = x1 x2 s_NN); the elastic one as set by the switches above
  INTEGER,PUBLIC::photon_final_state=0
  REAL(KIND(1d0)),PRIVATE::mu2_photon=0d0 ! the scale of the inelastic photons in PhotonPhotonFlux_pp_fxn (GeV^2)
  ! alpha_norm_mu(mu2_photon) of ElasticFlux (the MU2 argument of the elastic proton functions): with alpha_running4photon
  ! each elastic photon (proton and nucleus) of the photon-photon fluxes is normalised to alpha(mu2) as the inelastic
  ! one (LUXqed, MSbar at mu = W), so the hard process uses alpha(mu2); 1 with fixed alpha
  REAL(KIND(1d0)),PRIVATE::ael_photon=1d0
  ! pp and pA with the exact proton distribution (USE_EXACT_ELASTIC4PHOTON): the proton b integrals start at EXACT_BLOW
  ! times the lower limit of the small-x approximation (0.1 Rproton), as below it lies up to 9% of the flux at x = 0.5;
  ! the part left out scales as b^4 (about 1e-5)
  REAL(KIND(1d0)),PARAMETER,PRIVATE::EXACT_BLOW=0.1d0
  PUBLIC::PhotonPhotonFlux_pp,PhotonPhotonFlux_pp_eval
  ! (M11) pp, photon_final_state 1 and 2, exact elastic photon: the (1-P) parts in momentum space,
  !   int d^2b1 d^2b2 N1 N2 (1-P)(b1-b2) = (1/2 pi) int q dq (1-P)^(q) N1^(q) N2^(q),
  ! (1-P)^ of P = (1-exp(-b^2/(2 b0)))^2 in closed form, N^ the two-dimensional Fourier transforms of the elastic
  ! (elastic_rhohat) and inelastic (inel_rhohat, at QS_NC1 Chebyshev nodes in ln q on [1e-6, QS_QM] and QS_NC2 in q on
  ! [QS_QM, 1.6] GeV, interpolated: exp(-B Delta^2) of rho /= delta) densities; q nodes:
  ! Gauss-Legendre 6 on 2 panels per decade of q in [1e-6, 1.6] GeV. Otherwise the legacy b-space integration.
  ! (M15) also pA with the Woods-Saxon profile (pA_qspace_*)
  LOGICAL,PUBLIC::pp_final_state_qspace=.TRUE.
  ! (M11) the photon-photon grid of pp states 1-3: nodes computed when first needed (the same values); (M15) also pA
  LOGICAL,PUBLIC::pp_grid_lazy=.TRUE.
  ! (M16) the grids of states 1-3: 40 and 20 nodes per decade for x > 0.1 and 0.01 < x < 0.1 (.FALSE.: 10, as v2.0)
  LOGICAL,PUBLIC::pp_grid_fine=.TRUE.
  INTEGER,PARAMETER,PRIVATE::QS_NQ=78,QS_NC1=12,QS_NC2=12,QS_NC=QS_NC1+QS_NC2,QS_MAXE=4096,QS_MAXI=16384
  REAL(KIND(1d0)),PARAMETER,PRIVATE::QS_LO=-13.815510557964274d0,QS_HI=0.47000362924573558d0 ! ln 1e-6, ln 1.6
  REAL(KIND(1d0)),PARAMETER,PRIVATE::QS_QM=0.05d0 ! the inelastic transforms: Chebyshev in ln q below, in q above
  REAL(KIND(1d0)),DIMENSION(QS_NQ),PRIVATE,SAVE::qs_q,qs_w ! nodes; weights q^2 dln q (1-P)^/(2 pi)
  REAL(KIND(1d0)),DIMENSION(QS_NC),PRIVATE,SAVE::qs_qc ! Chebyshev nodes of the inelastic transforms
  REAL(KIND(1d0)),DIMENSION(QS_NQ,QS_NC),PRIVATE,SAVE::qs_ci ! interpolation matrix
  REAL(KIND(1d0)),PRIVATE,SAVE::qs_b0=-1d0
  INTEGER,PRIVATE,SAVE::qs_ne=0,qs_ni=0,qs_ekey=-1,qs_ikey=-1
  REAL(KIND(1d0)),DIMENSION(QS_MAXE),PRIVATE,SAVE::qs_ex
  REAL(KIND(1d0)),DIMENSION(QS_NQ,QS_MAXE),PRIVATE,SAVE::qs_er
  REAL(KIND(1d0)),DIMENSION(QS_MAXI),PRIVATE,SAVE::qs_ix,qs_imu2
  REAL(KIND(1d0)),DIMENSION(QS_NC,QS_MAXI),PRIVATE,SAVE::qs_ir
  ! (M15) pA (p+A, A+p: Woods-Saxon nucleus with ChFF, exact elastic proton) final states 1, 2 in momentum space, as pp:
  ! L = (1/2 pi) int q dq (1-P_pA)^ rho_p^ rho_A^. q nodes: Gauss-Legendre QA_NG on panels of width 1.1 in ln q from
  ! 1e-6 GeV to QA_Q1R/R_A, of width QA_DQR/R_A in q from there to 1.6 GeV ((1-P_pA)^ and the nucleus transforms
  ! oscillate with period ~2 pi/R_A); weights q dq (1-P_pA)^/(2 pi). The inelastic transforms at their own Chebyshev
  ! nodes (QA_NC1 in ln q on [1e-6, QS_QM], QA_NC2 in q on [QS_QM, 1.6] GeV: twice pp's, as (1-P_pA)^ weights q ~
  ! 0.01-0.3 GeV more; with 12 + 12 the states are off by up to 9e-6 (dfjk), with 24 + 24 by <= 3.5e-6, the noise of
  ! the fast transforms: drivers/m15_pA_test against exact transforms on 32 + 24 nodes), interpolated. Stores of the proton elastic and the nucleus transforms at the q nodes and at
  ! q = 0 (node qa_nq+1, weight 0: the collinear fluxes, int d^2b N), and of the inelastic transforms per x, mu2
  REAL(KIND(1d0)),PARAMETER,PRIVATE::QA_Q1R=0.5d0,QA_DQR=1d0
  INTEGER,PARAMETER,PRIVATE::QA_NG=6,QA_NC1=24,QA_NC2=24,QA_NC=QA_NC1+QA_NC2
  INTEGER,PARAMETER,PRIVATE::QA_MAXE=2048,QA_MAXN=2048,QA_MAXI=16384
  INTEGER,PRIVATE,SAVE::qa_nq=0,qa_ne=0,qa_nn=0,qa_ekey=-1,qa_ni=0,qa_ikey=-1
  REAL(KIND(1d0)),DIMENSION(QA_NC),PRIVATE,SAVE::qa_qc
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::qa_ix,qa_imu2
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PRIVATE,SAVE::qa_ir
  REAL(KIND(1d0)),DIMENSION(6),PRIVATE,SAVE::qa_key=-1d0 ! R_A, w, a (fm), A, sigma_NN (fm^2), use_MC_Glauber
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::qa_q,qa_w,qa_ex,qa_nx
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PRIVATE,SAVE::qa_ci,qa_er,qa_nr
  ! (M15) the exact collinear flux of the nucleus: ln rho_A^(q=0; x) on QA_CPD nodes per decade of x in
  ! [QA_CXLO, QA_CXHI], 8-point Lagrange in ln x (direct outside)
  INTEGER,PARAMETER,PRIVATE::QA_CPD=40
  REAL(KIND(1d0)),PARAMETER,PRIVATE::QA_CXLO=1d-10,QA_CXHI=0.5d0
  INTEGER,PRIVATE,SAVE::qa_cn=0
  REAL(KIND(1d0)),PRIVATE,SAVE::qa_ch=0d0
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::qa_ct
  PUBLIC::PhotonFlux_proton_EDFF_nob,PhotonFlux_proton_ChFF_nob,PhotonFlux_proton_ChFF_exact_nob
  PUBLIC::PhotonPhotonFlux_pA_hardsphere,PhotonPhotonFlux_pA_hardsphere_eval
  PUBLIC::PhotonPhotonFlux_pA_WoodsSaxon,PhotonPhotonFlux_pA_WoodsSaxon_eval
  PUBLIC::DeltaB_pp_at_x1x2,PhotonPhotonDeltaB_pp,PhotonPhotonDeltaB_pp_eval
  PUBLIC::DeltaB_pA_at_x1x2,PhotonPhotonDeltaB_pA_WoodsSaxon,PhotonPhotonDeltaB_pA_WoodsSaxon_eval
  PUBLIC::PhotonFlux_nucleus_EDFF_nob,PhotonFlux_nucleus_ChFF_nob
  PUBLIC::PhotonPhotonFlux_AB_hardsphere,PhotonPhotonFlux_AB_hardsphere_eval
  PUBLIC::PhotonPhotonFlux_AB_WoodsSaxon,PhotonPhotonFlux_AB_WoodsSaxon_eval
  PUBLIC::DeltaB_AB_at_x1x2,PhotonPhotonDeltaB_AB_WoodsSaxon,PhotonPhotonDeltaB_AB_WoodsSaxon_eval
  PUBLIC::Lgammagamma_UPC,dLgammagammadW_UPC,DeltaB_UPC_at_W
  PUBLIC::PhotonNumberDensity ! EDFF
  ! ChFF
  PUBLIC::PhotonNumberDensity_ChargeFormFactor_WS,PhotonNumberDensity_ChargeFormFactor_proton
  PUBLIC::PhotonNumberDensity_ChargeFormFactor_proton_exact
  ! Forward neutron tagging
  ! tagging neutrons for two beams
  ! -2 (0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)
  ! For the proton beam, it will be ignored
  INTEGER,DIMENSION(2),PUBLIC::neutron_tagging=(/-2,-2/)
  ! If it is .TRUE. in A+B collisions and neutron_tagging(1)=!=neutron_tagging(2),
  ! it includes (neutron_tagging(1),neutron_tagging(2))+(neutron_tagging(2),neutron_tagging(1))
  LOGICAL,PUBLIC::neutron_conj_sum=.FALSE.
  PUBLIC::Pbreak_AB_WoodsSaxon,Pbreak_pA_WoodsSaxon
  PUBLIC::PbreakXn_LO_AB_WS,PbreakXn_LO_pA
  PUBLIC::Pbreak1n_LO_AB_WS,Pbreak1n_LO_pA
  PUBLIC::Pbreak2n_LO_AB_WS,Pbreak2n_LO_pA
  PUBLIC::Pbreak3n_LO_AB_WS,Pbreak3n_LO_pA
  PUBLIC::Pbreak4n_LO_AB_WS,Pbreak4n_LO_pA
  ! For generating virtualities of initial photons
  PUBLIC::generate_Q2_epa_proton_iWW,generate_Q2_epa_ion_ChFF
  PUBLIC::InitialMomentumReshuffle_PhotonPhoton
  PUBLIC::BOOSTTOECM,BOOSTFROMECM
  PUBLIC::GetASymbol,GetNuclearInfo
  INTEGER,PARAMETER,PRIVATE::SUB_FACTOR=2
  LOGICAL,PRIVATE,SAVE::print_banner=.FALSE.
  LOGICAL,PUBLIC::use_MC_Glauber=.FALSE.
  INTEGER,PUBLIC::proton_FF_correction=0 ! when USE_CHARGEFORMFACTOR4PHOTON=.TRUE., 
                                         ! the proton form factors: 0 dipole, 1 dde1 (the correction beyond the dipole
                                         ! of the former proton_FF_correction=.TRUE.), 2 luxqed (A1 fit)
                                         ! it will impact both pp and pA UPCs
  INTEGER,PUBLIC::IEPSILON_EDFF=0  ! 0: bmin = RA; 1: bmin = RA+a
  LOGICAL,PRIVATE::CALC_DELTABQ=.FALSE. ! whether we want to calculate DeltaB (.TRUE.) or not (.FALSE.)
  LOGICAL,PUBLIC::GENERATE_DELTAB_GRID=.FALSE. ! if true, it will always first generate the grid
  LOGICAL,PUBLIC::GENERATE_PhotonPhotonFlux_GRID=.TRUE. ! if true, it will always first generate the grid for PhotonPhotonFlux
  ! A, Z of the beams (A = Z = 0: proton) and energy in GeV per nucleon, from run90.inc (MadGraph: ebeamMG5 of the
  ! whole beam particle, nb_proton, nb_neutron), set by the callers (v2.0)
  INTEGER,PRIVATE,SAVE::nuclearA_beam1,nuclearA_beam2,nuclearZ_beam1,nuclearZ_beam2
  REAL(KIND(1d0)),DIMENSION(2),PRIVATE,SAVE::ebeam_PN
  ! minimal energies (GeV) to be integrated in neutron tagging (in target frame)
  ! Emin_Xn, Emin_1n, Emin_2n, Emin_3n, Emin_4n
  REAL(KIND(1d0)),PARAMETER,PRIVATE::Emin_Xn=5d-3, Emin_1n=5d-3, Emin_2n=1d-2
  REAL(KIND(1d0)),PARAMETER,PRIVATE::Emin_3n=1.5d-2, Emin_4n=2d-2
  ! maximal energies (GeV) to be integrated in neutron tagging (in target frame)
  ! if it is 0, it means we take all possible energies (xgamma up to unity)
  ! Emax_Xn, Emax_1n, Emax_2n, Emax_3n, Emax_4n
  REAL(KIND(1d0)),PARAMETER,PRIVATE::Emax_Xn=0d0, Emax_1n=2d-1, Emax_2n=2d-1
  REAL(KIND(1d0)),PARAMETER,PRIVATE::Emax_3n=2d-1, Emax_4n=2d-1
  ! it will calculate cross sections by increasing all parameters by XSIGMA sigma.
  REAL(KIND(1d0)),PUBLIC::NEUTRON_XSIGMA=0d0
  ! ---------------------------------------------------------------------------
  ! faster-gammaUPC engine (FFE): kernel-factorised photon-photon fluxes.
  ! The flux integral (2pi alpha^2 Z1^2 Z2^2/(x1 x2)) int dln b1 dln b2 W1 W2 Phi(b1,b2),
  ! W = b^2 N, Phi = int_0^{2pi} dphi P(b12), is evaluated through an x-independent matrix
  ! K (Phi integrated against Lagrange basis functions of fixed 1D grids in ln b):
  ! Z = W1 K W2^T. See ../faster-gammaUPC/PLAN.md and docs/benchmark.pdf.
  ! ---------------------------------------------------------------------------
  ! AB and pA Woods-Saxon fluxes and <Delta b> (not pp, hard sphere, or pA with proton_FF_correction):
  INTEGER,PUBLIC::PhotonPhotonFlux_Algorithm=2 ! 2 (default): kernel-factorised engine; 1: legacy integration
  LOGICAL,PUBLIC::FFE_Legacy_ChFF_Definitions=.FALSE. ! emulate the legacy ChFF photon density
                                               ! (transverse only for b<3R, EPA beyond, b>=0.1R)
  REAL(KIND(1d0)),PUBLIC::FFE_NodeScale=1d0    ! >1: denser quadrature (convergence tests)
  INTEGER,PUBLIC::FFE_Verbose=0                ! >0: print timings and diagnostics
  LOGICAL,PUBLIC::FFE_Tag_Legacy=.FALSE.       ! neutron tagging from the legacy Pbreak functions
  ! Grid cache: the final grids (with the 1D tables of the PNOHAD=1 guard) are written to and read
  ! from files in FFE_Grid_Cache_Dir, one per configuration, so that the processes of one run
  ! (e.g. the MadEvent jobs) build each grid once. Off by default.
  LOGICAL,PUBLIC::FFE_Grid_Cache=.FALSE.
  CHARACTER(len=512),PUBLIC::FFE_Grid_Cache_Dir='./'
  PUBLIC::FFE_Reset,FFE_Flux_AB,FFE_DeltaB_AB,FFE_Info,FFE_Test_Value
  TYPE ffe_beam_t
     LOGICAL::proton=.FALSE.,chff=.FALSE.
     REAL(KIND(1d0))::anum=0d0,znum=0d0,Rfm=0d0,afm=0d0,w=0d0 ! A, Z, Woods-Saxon parameters (fm)
     REAL(KIND(1d0))::R=0d0,a=0d0                         ! the same in GeV^-1
     REAL(KIND(1d0))::ebeam=0d0,gamma=1d0,mass=0d0        ! E per nucleon, Lorentz factor, x*mass*b = b*E/gamma
     REAL(KIND(1d0))::bmin=0d0                            ! lower cut on b (GeV^-1)
     REAL(KIND(1d0))::rho0=0d0                            ! normalisation of rho (GeV^3)
  END TYPE ffe_beam_t
  TYPE ffe_vgrid_t
     INTEGER::n=0,nseg=0
     REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::v,b,qw     ! nodes v = ln b, b, int L_j dv
     INTEGER,DIMENSION(:),ALLOCATABLE::slo,shi            ! node range of each uniform segment
     REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::sv0,sh     ! first node and spacing of each segment
  END TYPE ffe_vgrid_t
  TYPE ffe_chff_t                                         ! x-independent position-space data
     INTEGER::ns=0,nz=0,nmom=0,npan=0
     REAL(KIND(1d0))::smax=0d0
     REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::pa,pb      ! s-panels
     REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::s,ws       ! s-nodes and weights (ns)
     INTEGER,DIMENSION(:),ALLOCATABLE::zoff               ! z-nodes of s-node i: zoff(i)+1..zoff(i+1)
     REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::z,wrho     ! z-node, weight*2*rho
     REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::mom      ! (0:nmom, ns) moments 2 int z^2k rho dz
     INTEGER,DIMENSION(:),ALLOCATABLE::jpan               ! table node j < smax -> panel ending at b_j
  END TYPE ffe_chff_t
  TYPE(ffe_beam_t),DIMENSION(2),PRIVATE,SAVE::ffe_bm
  TYPE(ffe_vgrid_t),DIMENSION(2),PRIVATE,SAVE::ffe_vg
  TYPE(ffe_chff_t),DIMENSION(2),PRIVATE,SAVE::ffe_cf
  INTEGER,PRIVATE,SAVE::ffe_icoll=0                       ! 2: pA (beam 1 = proton), 3: AB
  LOGICAL,PRIVATE,SAVE::ffe_swap=.FALSE.                  ! pA with the nucleus as beam 1
  REAL(KIND(1d0)),DIMENSION(40),PRIVATE,SAVE::ffe_sig=0d0 ! configuration signature
  LOGICAL,PRIVATE,SAVE::ffe_ready=.FALSE.
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_alpha=0d0,ffe_sigNN=0d0
  ! survival probability: ln chi on a uniform b grid (fm), P = exp(-chi), P = 1 beyond bcut
  INTEGER,PRIVATE,SAVE::ffe_nP=0
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_hb=0d0,ffe_bcut_fm=0d0
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::ffe_lnchi
  ! neutron-tagging probability on a uniform ln b grid (b in GeV^-1)
  LOGICAL,PRIVATE,SAVE::ffe_tag=.FALSE.
  INTEGER,PRIVATE,SAVE::ffe_nT=0
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_Tv0=0d0,ffe_Th=0d0
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::ffe_PtagT
  INTEGER,PRIVATE,SAVE::ffe_nTF=0
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_ThF=0d0
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::ffe_PtagF
  ! kernel breakpoints (GeV^-1)
  INTEGER,PRIVATE,SAVE::ffe_nbeta=0
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::ffe_beta
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_bcut=0d0
  LOGICAL,PRIVATE,SAVE::ffe_pjump=.FALSE.   ! P jumps to 1 at bcut (closed-form MC-Glauber fits)
  ! kernel matrices: 1 flux, 2 <Delta b> numerator (fm), 3 <Delta b> numerator with PNOHAD=1 (fm)
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_K
  LOGICAL,DIMENSION(3),PRIVATE,SAVE::ffe_Kready=.FALSE.
  INTEGER,PRIVATE,SAVE::ffe_nodes=0,ffe_nguard=0
  ! grids in u = -log10 x: piecewise uniform, 80 nodes per decade for x > 1e-2 (the ChFF flux has a
  ! shoulder at the first zero of the form factor), 20 below; node 1 is x = 1, node ffe_MX+1 is xmin
  TYPE(ffe_vgrid_t),PRIVATE,SAVE::ffe_xg
  INTEGER,PRIVATE,SAVE::ffe_MX=0
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_xmin=1d-8
  LOGICAL,PRIVATE,SAVE::ffe_tables=.FALSE.   ! v-grids, kernel and density tables of the configuration
  LOGICAL,PRIVATE,SAVE::ffe_lnNready=.FALSE. ! the 1D tables ffe_lnN
  INTEGER,PARAMETER,PRIVATE::ffe_cache_version=1  ! bump when the grids of a configuration change
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_W  ! (MX+1, nv, beam)
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_ZL ! ln(x1 x2 Z), kinds 1..3
  LOGICAL,DIMENSION(3),PRIVATE,SAVE::ffe_Zready=.FALSE.
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_lnN ! ln(int W dv) per grid x and beam
  REAL(KIND(1d0)),DIMENSION(6),PRIVATE,SAVE::ffe_time=0d0         ! setup, W, P, K, grid, total (CPU s)
  ! fine P table and far-zone moment tables of the kernel
  INTEGER,PARAMETER,PRIVATE::FFE_NMOM=26
  INTEGER,PRIVATE,SAVE::ffe_nf=0,ffe_ng=0,ffe_nw=0
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_b0pp=0d0       ! pp: slope b0 (GeV^-2) of P = (1-exp(-b^2/(2 b0)))^2
  REAL(KIND(1d0)),PRIVATE,SAVE::ffe_hf=0d0,ffe_hg=0d0,ffe_hw=0d0
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_ht   ! tagged inner moments H_n(w)
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PRIVATE,SAVE::ffe_Pf
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_gt
  ! near-zone node pool (v1, v2, c1..c3) used while the s-adaptation runs
  INTEGER,PRIVATE,SAVE::ffe_npool=0
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_pool
  ! patch near x = 1 for the proton photon of the inelastic branch (ffe_ourW), which vanishes at x = 1 (ln Z
  ! singular at the node x = 1 of the x grid): nodes uniform in s = -log10(1-x), s = FFE_SH0 + (k-1) FFE_HH,
  ! k = 1..FFE_NH; used for FFE_XHLO <= x <= FFE_XHHI (8-point Lagrange in s), direct evaluation above
  INTEGER,PARAMETER,PRIVATE::FFE_NH=76
  REAL(KIND(1d0)),PARAMETER,PRIVATE::FFE_SH0=0.5d0,FFE_HH=0.05d0,FFE_XHLO=0.85d0,FFE_XHHI=0.9999d0
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_WH   ! (FFE_NH, nv, beam) W rows at the patch x
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_lnNH   ! (FFE_NH, beam) ln(int W dv)
  REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE,PRIVATE,SAVE::ffe_ZH ! (FFE_NH, MX+1, kind, side) ln(x1 x2 Z)
  LOGICAL,PRIVATE,SAVE::ffe_WHready=.FALSE.
  LOGICAL,DIMENSION(3),PRIVATE,SAVE::ffe_ZHready=.FALSE.
CONTAINS
  FUNCTION PNOHAD_pp(bx,by,b0)
    ! the probability of no hardonic interaction at impact b=(bx,by)
    ! for pp collisions
    ! typical value of b0=19.8 GeV-2 at the LHC
    ! fitted by DdE is b0=9.7511+0.222796*log(s/GeV**2)+0.0179103*log(s/GeV**2)**2 GeV-2 from 10^1 to 10^5 GeV dsqrt(s)
    ! a new one (see 2207.03012) b0=9.81+0.211*log(s/GeV**2)+0.0185*log(s/GeV**2)**2 GeV-2
    ! bx and by should be in unit of GeV-1
    IMPLICIT NONE
    REAL(KIND(1d0))::PNOHAD_pp
    REAL(KIND(1d0)),INTENT(IN)::bx,by,b0
    REAL(KIND(1d0))::b2,gammasb,exponent
    b2=bx**2+by**2
    exponent=b2/2d0/b0
    IF(exponent.GT.500d0)THEN
       PNOHAD_pp=1d0
    ELSE
       gammasb=DEXP(-exponent)
       PNOHAD_pp=DABS(1d0-gammasb)**2
    ENDIF
    RETURN
  END FUNCTION PNOHAD_pp

  FUNCTION PNOHAD_AB_hardsphere(bx,by,ABAB,RR,sigmaNN)
    ! the probability of no hardonic interaction at impact b=(bx,by)
    ! for AB collisions
    ! bx, by, RR should be in unit of fm
    ! sigmaNN should be in unit of fm^2
    IMPLICIT NONE
    REAL(KIND(1d0))::PNOHAD_AB_hardsphere
    REAL(KIND(1d0)),INTENT(IN)::bx,by,ABAB
    REAL(KIND(1d0)),DIMENSION(2),INTENT(IN)::RR
    REAL(KIND(1d0)),INTENT(IN)::sigmaNN ! inelastic NN xs (in unit of fm^2)
    REAL(KIND(1d0))::TAB
    TAB=TABhat_hardsphere_grid(bx,by,RR)*ABAB
    IF(TAB*sigmaNN.GT.500d0)THEN
       PNOHAD_AB_hardsphere=0d0
    ELSE
       PNOHAD_AB_hardsphere=DEXP(-TAB*sigmaNN)
    ENDIF
    RETURN
  END FUNCTION PNOHAD_AB_hardsphere

  FUNCTION PNOHAD_AB_WoodsSaxon(bx,by,RR,w,aa,A,Z,sigmaNN)
    ! the probability of no hardonic interaction at impact b=(bx,by)
    ! for AB collisions
    ! bx, by, RR should be in unit of fm
    ! sigmaNN should be in unit of fm^2
    IMPLICIT NONE
    REAL(KIND(1d0))::PNOHAD_AB_WoodsSaxon
    REAL(KIND(1d0)),INTENT(IN)::bx,by
    REAL(KIND(1d0)),DIMENSION(2),INTENT(IN)::RR,w,aa,A,Z
    REAL(KIND(1d0)),INTENT(IN)::sigmaNN ! inelastic NN xs (in unit of fm^2)
    REAL(KIND(1d0))::TAB,ABAB
    ABAB=A(1)*A(2)
    IF(.NOT.use_MC_Glauber)THEN
       TAB=TABhat_WoodsSaxon_grid(bx,by,RR,w,aa,A)*ABAB
    ELSE
       TAB=TAB_MC_Glauber(bx,by,A,Z)
    ENDIF
    IF(TAB*sigmaNN.GT.500d0)THEN
       PNOHAD_AB_WoodsSaxon=0d0
    ELSE
       PNOHAD_AB_WoodsSaxon=DEXP(-TAB*sigmaNN)
    ENDIF
    RETURN
  END FUNCTION PNOHAD_AB_WoodsSaxon

  FUNCTION PNOHAD_pA_hardsphere(bx,by,RR,A,sigmaNN)
    ! the probability of no hardonic interaction at impact b=(bx,by)
    ! for pA or Ap collisions
    ! bx, by, RR should be in unit of fm
    ! sigmaNN should be in unit of fm^2
    IMPLICIT NONE
    REAL(KIND(1d0))::PNOHAD_pA_hardsphere
    REAL(KIND(1d0)),INTENT(IN)::bx,by,RR,A
    REAL(KIND(1d0)),INTENT(IN)::sigmaNN ! inelastic NN xs (in unit of fm^2)
    REAL(KIND(1d0))::TA
    TA=TAhat_hardsphere(bx,by,RR)*A
    IF(TA*sigmaNN.GT.500d0)THEN
       PNOHAD_pA_hardsphere=0d0
    ELSE
       PNOHAD_pA_hardsphere=DEXP(-TA*sigmaNN)
    ENDIF
    RETURN
  END FUNCTION PNOHAD_pA_hardsphere

  FUNCTION PNOHAD_pA_WoodsSaxon(bx,by,RR,w,aa,A,sigmaNN)
    ! the probability of no hardonic interaction at impact b=(bx,by)
    ! for pA or Ap collisions
    ! bx, by, RR should be in unit of fm
    ! sigmaNN should be in unit of fm^2
    IMPLICIT NONE
    REAL(KIND(1d0))::PNOHAD_pA_WoodsSaxon
    REAL(KIND(1d0)),INTENT(IN)::bx,by,RR,w,aa,A
    REAL(KIND(1d0)),INTENT(IN)::sigmaNN ! inelastic NN xs (in unit of fm^2)
    REAL(KIND(1d0))::TA
    IF(.NOT.use_MC_Glauber)THEN
       TA=TAhat_WoodsSaxon(bx,by,RR,w,aa,A,1,.FALSE.)*A
    ELSE
       TA=TA_MC_Glauber(bx,by,A)
    ENDIF
    IF(TA*sigmaNN.GT.500d0)THEN
       PNOHAD_pA_WoodsSaxon=0d0
    ELSE
       PNOHAD_pA_WoodsSaxon=DEXP(-TA*sigmaNN)
    ENDIF
    RETURN
  END FUNCTION PNOHAD_pA_WoodsSaxon

  FUNCTION PhotonPhotonFlux_pp(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_pp
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER::MX_save=0,MY_save=0
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD,ZD_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    INTEGER::gstate=0,gsel
    SAVE gstate
    REAL(KIND(1d0)),PARAMETER::ZUNSET=-HUGE(1d0) ! (M11) a node of the grid not computed yet
    LOGICAL::lazy
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_pp=0d0
       RETURN
    ENDIF
    IF(gstate.NE.photon_final_state)THEN
       ! a new grid for another final state
       init=0
       MX_save=0
       MY_save=0
       gstate=photon_final_state
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(ffe_supported())THEN
       ! kernel-factorised engine (the correction beyond the proton dipole form factor
       ! is only available in the legacy integration)
       PhotonPhotonFlux_pp=ffe_api(1,x1,x2,1,force_pnohad1,GENERATE_PhotonPhotonFlux_GRID)
       RETURN
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin).AND..NOT.force_pnohad1&
         .AND.GENERATE_PhotonPhotonFlux_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=0 ! (M16) nseg nodes per decade; states 1-3 with pp_grid_fine: more in the two decades below x = 1 (nsegd)
       DO I=0,log10xmin+1,-1
          MX=MX+nsegd(I)
       ENDDO
       MY=MX
       !IF(log10xmin.NE.log10xmin_before.or.init.EQ.0)THEN
       ! try to deallocate first
       IF(ALLOCATED(XD_1D))THEN
          DEALLOCATE(XD_1D)
       ENDIF
       ALLOCATE(XD_1D(MX+1))
       IF(ALLOCATED(YD_1D))THEN
          DEALLOCATE(YD_1D)
       ENDIF
       ALLOCATE(YD_1D(MY+1))
       IF(ALLOCATED(ZD))THEN
          WRITE(*,*)"INFO: in PhotonPhotonFlux_pp, the xmin of the grid has been updated to",xmin
          ! then deallocate it
          DEALLOCATE(ZD)
       ENDIF
       ALLOCATE(ZD(MX+1,MY+1))
       !ENDIF
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nsegd(I)
             log10x1=-1d0/DBLE(nsegd(I))*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x1
             YD_1D(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP 1
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       ! (M11) states 1-3 (pp_grid_lazy): the nodes are computed when a lookup first needs them (fill_block)
       lazy=pp_grid_lazy.AND.photon_final_state.NE.0
       IF(.NOT.lazy.AND.(photon_final_state.EQ.1.OR.photon_final_state.EQ.2))THEN
          IF(qspace_ok())THEN
             ! (M11) the transforms of all nodes and scales of the grid
             CALL pp_qspace_prepare(MX+1,10d0**XD_1D,MY+1,10d0**YD_1D)
          ELSE
             ! the inelastic densities at all nodes and scales of the grid, one row per x
             CALL pp_inelastic_prepare(MX+1,10d0**XD_1D,MY+1,10d0**YD_1D)
          ENDIF
       ENDIF
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! use the ones that have been calculated before
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1.AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_save(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             IF(lazy)THEN
                ZD(I,J)=ZUNSET
                CYCLE
             ENDIF
             ZD(I,J)=PhotonPhotonFlux_pp_eval(xx1,xx2)
             IF(MX_save.EQ.0.AND.MY_save.EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   WRITE(*,*)"INFO: generate grid of photon-photon flux in pp (will take tens of seconds)"
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF(ALLOCATED(ZD_save))THEN
          DEALLOCATE(ZD_save)
       ENDIF
       ALLOCATE(ZD_save(MX+1,MY+1))
       ! save the values that have been calculated before
       MX_save=MX
       MY_save=MY
       DO I=1,MX+1
          DO J=1,MY+1
             ZD_save(I,J)=ZD(I,J)
          ENDDO
       ENDDO
       init=1
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(GENERATE_PhotonPhotonFlux_GRID)THEN
          XI(1)=DLOG10(x1)
          YI(1)=DLOG10(x2)
          ! CALL lagrange_interp_2d(MX,MY,XD_1D,YD_1D,ZD,1,XI,YI,ZI)
          IF(XI(1).GE.0d0)THEN
             ilog10x1=-1
          ELSE
             ilog10x1=FLOOR(XI(1))
          ENDIF
          IF(YI(1).GE.0d0)THEN
             ilog10x2=-1
          ELSE
             ilog10x2=FLOOR(YI(1))
          ENDIF
          !K=nseg*(-ilog10x1-1)
          !DO I=1,nseg+1
          !   XD2_1D(I)=XD_1D(K+I)
          !ENDDO
          !L=nseg*(-ilog10x2-1)
          !DO I=1,nseg+1
          !   YD2_1D(I)=YD_1D(L+I)
          !ENDDO
          !DO I=1,nseg+1
          !   DO J=1,nseg+1
          !      ZD2(I,J)=ZD(K+I,L+J)
          !   ENDDO
          !ENDDO
          !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
          DO I=1,MX_save+1 ! (M16) the first node <= log10 x1 (nodes decreasing; as the per-decade search before)
             IF(XD_1D(I).LE.XI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save+1)THEN
             K=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             K=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save+1)THEN
             K=MX_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #1 !"
             STOP 1
          ENDIF
          DO I=1,MY_save+1
             IF(YD_1D(I).LE.YI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save+1)THEN
             L=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             L=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save+1)THEN
             L=MY_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #2 !"
             STOP 1
          ENDIF
          IF(pp_grid_lazy.AND.photon_final_state.NE.0)CALL fill_block(K,L)
          DO I=1,n_interp
             XD2_1D(I)=XD_1D(K+I)
             YD2_1D(I)=YD_1D(L+I)
          ENDDO
          DO I=1,n_interp
             DO J=1,n_interp
                ZD2(I,J)=ZD(K+I,L+J)
             ENDDO
          ENDDO
          CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
               YD2_1D,ZD2,1,XI,YI,ZI)
       ELSE
          ZI(1)=PhotonPhotonFlux_pp_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    IF((photon_final_state.EQ.1.OR.photon_final_state.EQ.2).AND..NOT.force_pnohad1)THEN
       ! the reference only checks the value (zero, NaN, sign): the elastic one (the collinear inelastic fluxes are
       ! numerical, 2 per call)
       gsel=photon_final_state
       photon_final_state=0
       IF(ffe_supported().AND.GENERATE_PhotonPhotonFlux_GRID)THEN
          ! (M11) from FFE's tables of int W dv, as the PNOHAD=1 flux of state 0 with the grid (the direct one
          ! computed two W rows per point of the x integral)
          pnohadval=ffe_api(1,x1,x2,1,.TRUE.,.TRUE.)
       ELSE
          pnohadval=PhotonPhotonFlux_pp_eval(x1,x2,.TRUE.)
       ENDIF
       photon_final_state=gsel
    ELSE
       pnohadval=PhotonPhotonFlux_pp_eval(x1,x2,.TRUE.)
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pp=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.&
            (DABS(ZI(1)/pnohadval).GT.1D2.AND.photon_final_state.NE.1.AND.photon_final_state.NE.2))THEN
          PhotonPhotonFlux_pp=pnohadval
       ELSE
          PhotonPhotonFlux_pp=ZI(1)
       ENDIF
    ELSE
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pp=0d0
       ELSE
          PhotonPhotonFlux_pp=pnohadval
       ENDIF
    ENDIF
    RETURN
  CONTAINS
    INTEGER FUNCTION nsegd(id)
      ! (M16) nodes in the decade [id-1, id] of log10 x: nseg; states 1-3 with pp_grid_fine 40 (x > 0.1) and 20
      ! (0.01 < x < 0.1), where the grid lookup error of the steep fall towards x -> 1 was 1e-4 to 1e-2 of dL/dW
      INTEGER,INTENT(IN)::id
      nsegd=nseg
      IF(pp_grid_fine.AND.photon_final_state.NE.0)THEN
         IF(id.EQ.0)nsegd=40
         IF(id.EQ.-1)nsegd=20
      ENDIF
    END FUNCTION nsegd
    SUBROUTINE fill_block(K0,L0)
      ! (M11) the nodes of the stencil K0+1..K0+n_interp, L0+1..L0+n_interp not computed yet: the transforms (or the
      ! b-space rows) of their x and scales first, then PhotonPhotonFlux_pp_eval (the values the full grid has)
      INTEGER,INTENT(IN)::K0,L0
      REAL(KIND(1d0)),DIMENSION(n_interp)::xa,xb
      INTEGER::ii,jj
      IF(ALL(ZD(K0+1:K0+n_interp,L0+1:L0+n_interp).NE.ZUNSET))RETURN
      xa=10d0**XD_1D(K0+1:K0+n_interp)
      xb=10d0**YD_1D(L0+1:L0+n_interp)
      IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
         IF(qspace_ok())THEN
            ! the transforms of each x of the stencil at all scales of its grid row and column down to the
            ! smallest scale of the stencil (one call per x shares the scale-independent work; the small scales,
            ! the costly ones, only when a lookup needs them)
            CALL pp_qspace_prepare_grid(n_interp,xa,MY+1,10d0**YD_1D,n_interp,xb,MX+1,10d0**XD_1D,&
                 photon_mu2(xa(n_interp),xb(n_interp))*(1d0-1d-12))
         ELSE
            CALL pp_inelastic_prepare(n_interp,xa,n_interp,xb)
         ENDIF
      ENDIF
      DO jj=1,n_interp
         DO ii=1,n_interp
            IF(ZD(K0+ii,L0+jj).NE.ZUNSET)CYCLE
            ZD(K0+ii,L0+jj)=PhotonPhotonFlux_pp_eval(xa(ii),xb(jj))
            ZD_save(K0+ii,L0+jj)=ZD(K0+ii,L0+jj)
         ENDDO
      ENDDO
    END SUBROUTINE fill_block
  END FUNCTION PhotonPhotonFlux_pp

  FUNCTION PhotonPhotonFlux_pp_eval(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_pp_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::b0_common ! in unit of GeV-2
    COMMON/PhotonPhoton_pp/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,b0_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact
    SAVE bfact
    REAL(KIND(1d0)),PARAMETER::bupper=2d0
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    REAL(KIND(1d0))::cmenergy,el1,el2
    INTEGER,DIMENSION(3)::subn
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_pp_eval=0d0
       RETURN
    ENDIF
    IF(ffe_supported())THEN
       force_pnohad1=.FALSE.
       IF(PRESENT(FORCEPNOHAD1))force_pnohad1=FORCEPNOHAD1
       PhotonPhotonFlux_pp_eval=ffe_api(1,x1,x2,1,force_pnohad1,.FALSE.)
       RETURN
    ENDIF
    CALC_DELTABQ=.FALSE.
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       gamma1_common=ebeam_PN(1)/mproton
       gamma2_common=ebeam_PN(2)/mproton
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! get b0 from the DdE fit
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2)) ! in unit of GeV
       !b0_common=9.7511D0+0.222796D0*DLOG(cmenergy**2)&
       !     +0.0179103D0*DLOG(cmenergy**2)**2 ! in unit of GeV-2
       b0_common=9.81D0+0.211D0*DLOG(cmenergy**2)&
             +0.0185D0*DLOG(cmenergy**2)**2 ! in unit of GeV-2
       ! two Z are 1
       Z1=1d0
       Z2=1d0
       ! 0.1973 is from fm to GeV-1
       bfact=Rproton/GeVm12fm*mproton
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          bfact=bfact*LOWER_BFactor_Limit
       ENDIF
       bbx(1)=bupper
       bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we should increase the number of segments
          sub_num(1)=sub_num(1)*SUB_FACTOR
          sub_num(2)=sub_num(2)*SUB_FACTOR
       ENDIF
       sub_num(3)=10
       init=1
    ENDIF
    x1_common=x1
    x2_common=x2
    E1_common=ebeam_PN(1)*x1
    E2_common=ebeam_PN(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    CALL set_photon_scale(x1,x2)
    IF(force_pnohad1.OR.photon_final_state.EQ.3)THEN
       ! we only use PNOHAD=1
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          el1=PhotonFlux_proton_nob(x1_common,gamma1_common)
          el2=PhotonFlux_proton_nob(x2_common,gamma2_common)
       ELSE
          el1=PhotonFlux_proton_EDFF_nob(x1_common,gamma1_common)
          el2=PhotonFlux_proton_EDFF_nob(x2_common,gamma2_common)
       ENDIF
       PhotonPhotonFlux_pp_eval=final_state_nob(el1*ael_photon,el2*ael_photon,x1_common,x2_common,2)
       PhotonPhotonFlux_pp_eval=MAX(PhotonPhotonFlux_pp_eval,0d0)
       RETURN
    ENDIF
    IF((photon_final_state.EQ.1.OR.photon_final_state.EQ.2).AND.qspace_ok())THEN
       ! (M11) momentum space: 1: (el-inel + inel-el) P = collinear products - (1-P) parts; 2: (1-P) (el-el + el-inel +
       ! inel-el) + inel-inel
       CALL pp_qspace_L(x1,x2,el1,el2) ! el1 = L(el,el), el2 = L(el,inel) + L(inel,el) (alpha(0) elastic)
       IF(photon_final_state.EQ.1)THEN
          integral=-alpha**2/(x1*x2)*ael_photon*el2
          el1=PhotonFlux_proton_nob(x1_common,gamma1_common)
          el2=PhotonFlux_proton_nob(x2_common,gamma2_common)
       ELSE
          integral=alpha**2/(x1*x2)*(ael_photon**2*el1+ael_photon*el2)
          el1=0d0; el2=0d0
       ENDIF
       PhotonPhotonFlux_pp_eval=integral+final_state_nob(el1*ael_photon,el2*ael_photon,x1_common,x2_common,2)
       RETURN
    ENDIF
    IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
       CALL pp_inelastic_prepare(1,(/x1/),1,(/x2/))
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact*x1)
    ! b2*x2*mproton = Exp(bA(2))
    aax(2)=DLOG(bfact*x2)
    IF(USE_CHARGEFORMFACTOR4PHOTON.AND.USE_EXACT_ELASTIC4PHOTON)THEN
       aax(1)=aax(1)+DLOG(EXACT_BLOW)
       aax(2)=aax(2)+DLOG(EXACT_BLOW)
    ENDIF
    subn=sub_num
    IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
       ! down to the b of the inelastic photons at large k_T
       CALL pp_inelastic_range(x1,aax(1),subn(1))
       CALL pp_inelastic_range(x2,aax(2),subn(2))
    ENDIF
    CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax,bbx,3,subn,1,1d-5,&
         integral,ind,eval_num)
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0.AND.photon_final_state.NE.1)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          bbx(1)=bupper*2d0**(iter)
          bbx(2)=bupper*2d0**(iter)
          CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax,bbx,3,subn,1,1d-5,&
               integral,ind,eval_num)
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       bbx(2)=bupper
    ENDIF
    IF(integral.LT.0d0.AND.photon_final_state.NE.1)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon flux at (x1,x2)=",x1_common,x2_common
          WRITE(*,*)"WARNING: use PNOHAD=1 approx. instead (most probably need to increase bupper)"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for approximation, we just use dipole
          PhotonPhotonFlux_pp_eval=PhotonFlux_proton_ChFF_nob(x1_common,gamma1_common)
          PhotonPhotonFlux_pp_eval=PhotonPhotonFlux_pp_eval*&
               PhotonFlux_proton_ChFF_nob(x2_common,gamma2_common)
       ELSE
          PhotonPhotonFlux_pp_eval=PhotonFlux_proton_EDFF_nob(x1_common,gamma1_common)
          PhotonPhotonFlux_pp_eval=PhotonPhotonFlux_pp_eval*&
               PhotonFlux_proton_EDFF_nob(x2_common,gamma2_common)
       ENDIF
       PhotonPhotonFlux_pp_eval=MAX(PhotonPhotonFlux_pp_eval,0d0)*ael_photon**2
    ELSE
       PhotonPhotonFlux_pp_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
       IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
          ! 1: (el-inel + inel-el) P = (el-inel + inel-el) - (el-inel + inel-el) (1-P), the first from the collinear
          ! fluxes (the b integral of the inelastic photon rings with the sharp k_T < mu cut), the second, local in
          ! b12, integrated above; 2: + inel-inel (P=1)
          el1=0d0
          el2=0d0
          IF(photon_final_state.EQ.1)THEN
             IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
                el1=PhotonFlux_proton_nob(x1_common,gamma1_common)
                el2=PhotonFlux_proton_nob(x2_common,gamma2_common)
             ELSE
                el1=PhotonFlux_proton_EDFF_nob(x1_common,gamma1_common)
                el2=PhotonFlux_proton_EDFF_nob(x2_common,gamma2_common)
             ENDIF
          ENDIF
          PhotonPhotonFlux_pp_eval=PhotonPhotonFlux_pp_eval+&
               final_state_nob(el1*ael_photon,el2*ael_photon,x1_common,x2_common,2)
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_pp_eval

  FUNCTION DeltaB_pp_at_x1x2(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::DeltaB_pp_at_x1x2
    REAL(KIND(1d0))::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    LOGICAL::force_pnohad1
    REAL(KIND(1d0))::num,den
    LOGICAL::use_grid_bu
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    ! numerator
    ! PhotonPhotonFlux*Delta b in unit of fm
    num=PhotonPhotonDeltaB_pp(x1,x2,force_pnohad1)
    ! denominator
    ! PhotonPhotonFlux
    use_grid_bu=GENERATE_PhotonPhotonFlux_GRID
    ! Let us make it to be coherent with DELTAB grid case
    GENERATE_PhotonPhotonFlux_GRID=GENERATE_DELTAB_GRID
    den=PhotonPhotonFlux_pp(x1,x2,force_pnohad1)
    GENERATE_PhotonPhotonFlux_GRID=use_grid_bu ! recover the original one
    IF(den.LE.0d0.OR.num.EQ.0d0)THEN
       DeltaB_pp_at_x1x2=0d0
    ELSE
       DeltaB_pp_at_x1x2=num/den
    ENDIF
    RETURN
  END FUNCTION DeltaB_pp_at_x1x2

  FUNCTION PhotonPhotonDeltaB_pp(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonDeltaB_pp
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER,DIMENSION(2)::MX_save=(/0,0/),MY_save=(/0,0/)
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::ZD
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD1_save,ZD2_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD1_save,ZD2_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    LOGICAL,DIMENSION(2)::gridready=(/.FALSE.,.FALSE./)
    SAVE gridready
    INTEGER::igrid
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonDeltaB_pp=0d0
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(ffe_supported())THEN
       PhotonPhotonDeltaB_pp=ffe_api(1,x1,x2,2,force_pnohad1,GENERATE_DELTAB_GRID)
       RETURN
    ENDIF
    IF(force_pnohad1)THEN
       igrid=2
    ELSE
       igrid=1
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin.OR..NOT.gridready(igrid)).AND.GENERATE_DELTAB_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=nseg*(-log10xmin)
       MY=MX
       !IF(log10xmin.NE.log10xmin_before.or.init.EQ.0)THEN
       ! try to deallocate first
       IF(init.EQ.0.OR.log10xmin_before.GT.log10xmin)THEN
          IF(ALLOCATED(XD_1D))THEN
             DEALLOCATE(XD_1D)
          ENDIF
          ALLOCATE(XD_1D(MX+1))
          IF(ALLOCATED(YD_1D))THEN
             DEALLOCATE(YD_1D)
          ENDIF
          ALLOCATE(YD_1D(MY+1))
          IF(ALLOCATED(ZD))THEN
             WRITE(*,*)"INFO: in PhotonPhotonDeltaB_pp, the xmin of the grid has been updated to ",xmin
             ! save the values that have been calculated before
             IF(gridready(1))THEN
                ! the saved grid has the size before the xmin update
                MX_save(1)=nseg*(-log10xmin_before)
                MY_save(1)=MX_save(1)
                IF(ALLOCATED(ZD1_save))THEN
                   DEALLOCATE(ZD1_save)
                ENDIF
                ALLOCATE(ZD1_save(MX_save(1)+1,MY_save(1)+1))
                DO I=1,MX_save(1)+1
                   DO J=1,MY_save(1)+1
                      ZD1_save(I,J)=ZD(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(gridready(2))THEN
                ! the saved grid has the size before the xmin update
                MX_save(2)=nseg*(-log10xmin_before)
                MY_save(2)=MX_save(2)
                IF(ALLOCATED(ZD2_save))THEN
                   DEALLOCATE(ZD2_save)
                ENDIF
                ALLOCATE(ZD2_save(MX_save(2)+1,MY_save(2)+1))
                DO I=1,MX_save(2)+1
                   DO J=1,MY_save(2)+1
                      ZD2_save(I,J)=ZD(2,I,J)
                   ENDDO
                ENDDO
             ENDIF
             DEALLOCATE(ZD)
          ENDIF
          ALLOCATE(ZD(2,MX+1,MY+1))
          K=0
          DO I=0,log10xmin+1,-1
             DO J=1,nseg
                log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
                K=K+1
                XD_1D(K)=log10x1
                YD_1D(K)=log10x1
             ENDDO
          ENDDO
          IF(K.NE.MX)THEN
             WRITE(*,*)"ERROR: K != MX"
             STOP 1
          ENDIF
          XD_1D(MX+1)=DBLE(log10xmin)
          YD_1D(MY+1)=DBLE(log10xmin)
          gridready(1:2)=.FALSE.
       ENDIF
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! save the values that have been calculated before
             IF(MX_save(igrid).GT.0.AND.MY_save(igrid).GT.0.AND.I.LE.MX_save(igrid)+1&
                  .AND.J.LE.MY_save(igrid)+1)THEN
                IF(igrid.EQ.1)THEN
                   ZD(igrid,I,J)=ZD1_save(I,J)
                ELSEIF(igrid.EQ.2)THEN
                   ZD(igrid,I,J)=ZD2_save(I,J)
                ELSE
                   WRITE(*,*)"ERROR: do not known igrid=",igrid
                ENDIF
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             ZD(igrid,I,J)=PhotonPhotonDeltaB_pp_eval(xx1,xx2,force_pnohad1)
             IF(MX_save(igrid).EQ.0.AND.MY_save(igrid).EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   IF(igrid.EQ.1)THEN
                      WRITE(*,*)"INFO: generate grid of photon-photon <delta b> in pp (will take tens of seconds)"
                   ELSE
                      WRITE(*,*)"INFO: generate grid of photon-photon <delta b> in pp with PNOHAD=1 (will take tens of seconds)"
                   ENDIF
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       gridready(igrid)=.TRUE.
       init=1
    ENDIF
    IF(GENERATE_DELTAB_GRID)THEN
       XI(1)=DLOG10(x1)
       YI(1)=DLOG10(x2)
       IF(XI(1).GE.0d0)THEN
          ilog10x1=-1
       ELSE
          ilog10x1=FLOOR(XI(1))
       ENDIF
       IF(YI(1).GE.0d0)THEN
          ilog10x2=-1
       ELSE
          ilog10x2=FLOOR(YI(1))
       ENDIF
       !K=nseg*(-ilog10x1-1)
       !DO I=1,nseg+1
       !   XD2_1D(I)=XD_1D(K+I)
       !ENDDO
       !L=nseg*(-ilog10x2-1)
       !DO I=1,nseg+1
       !   YD2_1D(I)=YD_1D(L+I)
       !ENDDO
       !DO I=1,nseg+1
       !   DO J=1,nseg+1
       !      ZD2(I,J)=ZD(igrid,K+I,L+J)
       !   ENDDO
       !ENDDO
       !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
       K=nseg*(-ilog10x1-1)
       DO I=K+1,K+nseg+1
          if(XD_1D(I).LE.XI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX+1)THEN
          K=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MX+1)THEN
          K=MX+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #1 !"
          STOP 1
       ENDIF
       L=nseg*(-ilog10x2-1)
       DO I=L+1,L+nseg+1
          if(YD_1D(I).LE.YI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY+1)THEN
          L=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          L=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MY+1)THEN
          L=MY+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #2 !"
          STOP 1
       ENDIF
       DO I=1,n_interp
          XD2_1D(I)=XD_1D(K+I)
          YD2_1D(I)=YD_1D(L+I)
       ENDDO
       DO I=1,n_interp
          DO J=1,n_interp
             ZD2(I,J)=ZD(igrid,K+I,L+J)
          ENDDO
       ENDDO
       CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
             YD2_1D,ZD2,1,XI,YI,ZI)
       IF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0)THEN
          PhotonPhotonDeltaB_pp=0d0
       ELSE
          PhotonPhotonDeltaB_pp=ZI(1)
       ENDIF
    ELSE
       PhotonPhotonDeltaB_pp=PhotonPhotonDeltaB_pp_eval(x1,x2,force_pnohad1)
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonDeltaB_pp

  FUNCTION PhotonPhotonDeltaB_pp_eval(x1,x2,FORCEPNOHAD1)
    ! this is PhotonPhotonFlux*Delta b, where Delta b is in unit of fm
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonDeltaB_pp_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::b0_common ! in unit of GeV-2
    COMMON/PhotonPhoton_pp/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,b0_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact
    SAVE bfact
    REAL(KIND(1d0)),PARAMETER::bupper=2d0
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    REAL(KIND(1d0))::cmenergy
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonDeltaB_pp_eval=0d0
       RETURN
    ENDIF
    IF(ffe_supported())THEN
       force_pnohad1=.FALSE.
       IF(PRESENT(FORCEPNOHAD1))force_pnohad1=FORCEPNOHAD1
       PhotonPhotonDeltaB_pp_eval=ffe_api(1,x1,x2,2,force_pnohad1,.FALSE.)
       RETURN
    ENDIF
    CALC_DELTABQ=.TRUE.
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       gamma1_common=ebeam_PN(1)/mproton
       gamma2_common=ebeam_PN(2)/mproton
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! get b0 from the DdE fit
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2)) ! in unit of GeV
       !b0_common=9.7511D0+0.222796D0*DLOG(cmenergy**2)&
       !     +0.0179103D0*DLOG(cmenergy**2)**2 ! in unit of GeV-2
       b0_common=9.81D0+0.211D0*DLOG(cmenergy**2)&
             +0.0185D0*DLOG(cmenergy**2)**2 ! in unit of GeV-2
       ! two Z are 1
       Z1=1d0
       Z2=1d0
       ! 0.1973 is from fm to GeV-1
       bfact=Rproton/GeVm12fm*mproton
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          bfact=bfact*LOWER_BFactor_Limit
       ENDIF
       bbx(1)=bupper
       bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we should increase the number of segments
          sub_num(1)=sub_num(1)*SUB_FACTOR
          sub_num(2)=sub_num(2)*SUB_FACTOR
       ENDIF
       sub_num(3)=10
       init=1
    ENDIF
    x1_common=x1
    x2_common=x2
    E1_common=ebeam_PN(1)*x1
    E2_common=ebeam_PN(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    CALL set_photon_scale(x1,x2)
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact*x1)
    ! b2*x2*mproton = Exp(bA(2))
    aax(2)=DLOG(bfact*x2)
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax(1:2),bbx(1:2),2,sub_num(1:2),1,1d-5,&
            integral,ind,eval_num)
    ELSE
       CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax,bbx,3,sub_num,1,1d-5,&
            integral,ind,eval_num)
    ENDIF
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          bbx(1)=bupper*2d0**(iter)
          bbx(2)=bupper*2d0**(iter)
          IF(force_pnohad1)THEN
             ! we only use PNOHAD=1
             CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax(1:2),bbx(1:2),2,sub_num(1:2),1,1d-5,&
                  integral,ind,eval_num)
          ELSE
             CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax,bbx,3,sub_num,1,1d-5,&
                  integral,ind,eval_num)
          ENDIF
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       bbx(2)=bupper
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon <delta b> at (x1,x2)=",x1_common,x2_common
          WRITE(*,*)"WARNING: most probably need to increase bupper"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       PhotonPhotonDeltaB_pp_eval=0d0
    ELSE
       PhotonPhotonDeltaB_pp_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonDeltaB_pp_eval

  FUNCTION PhotonPhotonFlux_pp_fxn(dim_num,bA)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_pp_fxn
    INTEGER,INTENT(IN)::dim_num ! should be 3 or 2
    ! 1/0.1973d0 from fm to GeV-1 for b 
    ! x1*b1*mproton=Exp(bA(1))
    ! x2*b2*mproton=Exp(bA(2))
    ! bA(3) = theta_{12}
    REAL(KIND(1d0)),DIMENSION(dim_num),INTENT(IN)::bA
    REAL(KIND(1d0))::b1,b2,b12,costh,pnohad
    REAL(KIND(1d0))::Ngamma1,Ngamma2,Ninel1,Ninel2,blo_el
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! as PhotonPhotonFlux_pp_eval
    !REAL(KIND(1d0)),PARAMETER::b0=19.8d0 ! in unit of GeV-2
    REAL(KIND(1d0))::b0_common ! in unit of GeV-2
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    COMMON/PhotonPhoton_pp/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,b0_common
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),EXTERNAL::elliptic_em
    REAL(KIND(1d0))::arg1,arg2
    IF(dim_num.NE.3.AND.dim_num.NE.2)THEN
       WRITE(*,*)"ERROR: PhotonPhotonFlux_pp_fxn is not a three or two dimensional function"
       STOP 1
    ENDIF
    ! in unit of GeV-1
    ! x1*b1*mproton=Exp(bA(1))
    b1=DEXP(bA(1))/x1_common/mproton
    ! x2*b2*mproton=Exp(bA(2))
    b2=DEXP(bA(2))/x2_common/mproton
    IF(dim_num.EQ.3)THEN
       costh=DCOS(bA(3))
       b12=DSQRT(b1**2+b2**2-2d0*b1*b2*costh)
       pnohad=PNOHAD_pp(b12,0d0,b0_common)
       ! hadronic interaction (1 and 2: the parts with 1-P, see PhotonPhotonFlux_pp_eval)
       IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)pnohad=1d0-pnohad
       IF(CALC_DELTABQ)THEN
          pnohad=pnohad*b12*GeVm12fm ! convert from GeV-1 to fm
       ENDIF
    ELSE
       ! dim_num = 2
       ! it must be a case with pnohad=1
       ! and CALC_DELTABQ=.TRUE.
       IF(.NOT.CALC_DELTABQ)THEN
          WRITE(*,*)"ERROR: cannot reach here with CALC_DELTABQ=.FALSE."
       ENDIF
       ! pnohad is Integrate[Sqrt[b1^2 + b2^2 - 2 b1*b2*Cos[phi]], {phi, 0, 2 Pi}]
       arg2=4d0*b1*b2/(b1+b2)**2
       IF(b1.NE.b2.AND.arg2.LT.1d0)THEN
          arg1=-4d0*b1*b2/(b1-b2)**2
          ! elliptic_em(m) is EllipticE[m] with m<1
          pnohad=DABS(b1-b2)*elliptic_em(arg1)
          pnohad=pnohad+(b1+b2)*elliptic_em(arg2)
       ELSE
          pnohad=4d0*b1
       ENDIF
       pnohad=pnohad*2d0*GeVm12fm ! convert from GeV-1 to fm
    ENDIF
    IF(pnohad.LE.0d0)THEN
       PhotonPhotonFlux_pp_fxn=0d0
       RETURN
    ENDIF
    blo_el=0d0
    IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
       ! the elastic photons on the b range of photon_final_state=0 (the lower limit of PhotonPhotonFlux_pp_eval);
       ! the range is extended for the inelastic photons only (pp_inelastic_range)
       blo_el=Rproton/GeVm12fm
       IF(USE_CHARGEFORMFACTOR4PHOTON)blo_el=blo_el*LOWER_BFactor_Limit
       IF(USE_CHARGEFORMFACTOR4PHOTON.AND.USE_EXACT_ELASTIC4PHOTON)blo_el=blo_el*EXACT_BLOW
    ENDIF
    Ngamma1=0d0
    Ngamma2=0d0
    IF(.NOT.USE_CHARGEFORMFACTOR4PHOTON)THEN
       ! EDFF
       IF(b1.GE.blo_el)Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
       IF(b2.GE.blo_el)Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    ELSE
       ! ChFF
       IF(b1.GE.blo_el)Ngamma1=PhotonNumberDensity_proton(b1,E1_common,gamma1_common,1)
       IF(b2.GE.blo_el)Ngamma2=PhotonNumberDensity_proton(b2,E2_common,gamma2_common,2)
    ENDIF
    ! the elastic photons normalised to alpha(mu2) with running alpha (set_photon_scale; the nucleus photon has no
    ! running alpha(Q2), Q2 <~ 1/R_A**2, only this normalisation)
    Ngamma1=Ngamma1*ael_photon
    Ngamma2=Ngamma2*ael_photon
    IF(photon_final_state.EQ.0)THEN
       PhotonPhotonFlux_pp_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
    ELSE
       Ninel1=PhotonNumberDensity_inelastic_beam(b1,E1_common,gamma1_common,1)
       Ninel2=PhotonNumberDensity_inelastic_beam(b2,E2_common,gamma2_common,2)
       IF(photon_final_state.EQ.1)THEN
          ! -(1-P) (el-inel + inel-el); (el-inel + inel-el) is added in PhotonPhotonFlux_pp_eval
          PhotonPhotonFlux_pp_fxn=-b1**2*b2**2*pnohad*(Ngamma1*Ninel2+Ninel1*Ngamma2)
       ELSE
          ! 2: (1-P) (el-el + el-inel + inel-el); inel-inel (P=1) is added in PhotonPhotonFlux_pp_eval
          PhotonPhotonFlux_pp_fxn=b1**2*b2**2*pnohad*(Ngamma1*Ngamma2+Ngamma1*Ninel2+Ninel1*Ngamma2)
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_pp_fxn

  FUNCTION PhotonPhotonFlux_pA_hardsphere(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_pA_hardsphere
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER::MX_save=0,MY_save=0
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD,ZD_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_pA_hardsphere=0d0
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin).AND..NOT.force_pnohad1&
         .AND.GENERATE_PhotonPhotonFlux_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=nseg*(-log10xmin)
       MY=MX
       ! try to deallocate first
       IF(ALLOCATED(XD_1D))THEN
          DEALLOCATE(XD_1D)
       ENDIF
       ALLOCATE(XD_1D(MX+1))
       IF(ALLOCATED(YD_1D))THEN
          DEALLOCATE(YD_1D)
       ENDIF
       ALLOCATE(YD_1D(MY+1))
       IF(ALLOCATED(ZD))THEN
          ! then deallocate it
          DEALLOCATE(ZD)
       ENDIF
       ALLOCATE(ZD(MX+1,MY+1))
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nseg
             log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x1
             YD_1D(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP 1
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1&
                  .AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_save(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             ZD(I,J)=PhotonPhotonFlux_pA_hardsphere_eval(xx1,xx2)
             IF(MX_save.EQ.0.AND.MY_save.EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   WRITE(*,*)"INFO: generate grid of photon-photon flux in pA or Ap (will take tens of seconds)"
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF(ALLOCATED(ZD_save))THEN
          DEALLOCATE(ZD_save)
       ENDIF
       ALLOCATE(ZD_save(MX+1,MY+1))
       ! save the values that have been calculated before
       MX_save=MX
       MY_save=MY
       DO I=1,MX+1
          DO J=1,MY+1
             ZD_save(I,J)=ZD(I,J)
          ENDDO
       ENDDO
       init=1
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(GENERATE_PhotonPhotonFlux_GRID)THEN
          XI(1)=DLOG10(x1)
          YI(1)=DLOG10(x2)
          IF(XI(1).GE.0d0)THEN
             ilog10x1=-1
          ELSE
             ilog10x1=FLOOR(XI(1))
          ENDIF
          IF(YI(1).GE.0d0)THEN
             ilog10x2=-1
          ELSE
             ilog10x2=FLOOR(YI(1))
          ENDIF
          !K=nseg*(-ilog10x1-1)
          !DO I=1,nseg+1
          !   XD2_1D(I)=XD_1D(K+I)
          !ENDDO
          !L=nseg*(-ilog10x2-1)
          !DO I=1,nseg+1
          !   YD2_1D(I)=YD_1D(L+I)
          !ENDDO
          !DO I=1,nseg+1
          !   DO J=1,nseg+1
          !      ZD2(I,J)=ZD(K+I,L+J)
          !   ENDDO
          !ENDDO
          !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
          K=nseg*(-ilog10x1-1)
          DO I=K+1,K+nseg+1
             if(XD_1D(I).LE.XI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save+1)THEN
             K=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             K=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save+1)THEN
             K=MX_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #1 !"
             STOP 1
          ENDIF
          L=nseg*(-ilog10x2-1)
          DO I=L+1,L+nseg+1
             if(YD_1D(I).LE.YI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save+1)THEN
             L=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             L=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save+1)THEN
             L=MY_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #2 !"
             STOP 1
          ENDIF
          DO I=1,n_interp
             XD2_1D(I)=XD_1D(K+I)
             YD2_1D(I)=YD_1D(L+I)
          ENDDO
          DO I=1,n_interp
             DO J=1,n_interp
                ZD2(I,J)=ZD(K+I,L+J)
             ENDDO
          ENDDO
          CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
               YD2_1D,ZD2,1,XI,YI,ZI)
       ELSE
          ZI(1)=PhotonPhotonFlux_pA_hardsphere_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    pnohadval=PhotonPhotonFlux_pA_hardsphere_eval(x1,x2,.TRUE.)
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pA_hardsphere=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.DABS(ZI(1)/pnohadval).GT.1D2)THEN
          PhotonPhotonFlux_pA_hardsphere=pnohadval
       ELSE
          PhotonPhotonFlux_pA_hardsphere=ZI(1)
       ENDIF
    ELSE
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pA_hardsphere=0d0
       ELSE
          PhotonPhotonFlux_pA_hardsphere=pnohadval
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_pA_hardsphere

  FUNCTION PhotonPhotonFlux_pA_hardsphere_eval(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_pA_hardsphere_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, only evaluate with PNOHAD=1
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0))::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    COMMON/PhotonPhoton_pA_HS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact1,bfact2
    SAVE bfact1,bfact2
    REAL(KIND(1d0)),PARAMETER::bupper=3d0
    REAL(KIND(1d0))::aaVal,wVal
    CHARACTER(len=7)::Aname
    REAL(KIND(1d0))::cmenergy
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_pA_hardsphere_eval=0d0
       RETURN
    ENDIF
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_common=ebeam_PN(2)/mproton
          gamma2_common=ebeam_PN(1)/mN
       ELSE
          gamma1_common=ebeam_PN(1)/mproton
          gamma2_common=ebeam_PN(2)/mN
       ENDIF
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! proton Z is 1
       Z1=1d0
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       IF(nuclearA_beam1.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       ELSEIF(nuclearA_beam2.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
       ELSE
          WRITE(*,*)"ERROR: please set nuclearA_beam1 or nuclearA_beam2 nonzero first !"
          STOP 1
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in p+"//TRIM(Aname)//" collisions"
       CALL GetNuclearInfo(Aname,A_common,Z2,RA_common,aaval,wval)
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
       sigNN_inel_common=sigma_inelastic(cmenergy)
       sigNN_inel_common=sigNN_inel_common*0.1d0 ! from mb to fm^2
       ! 0.1973 is from fm to GeV-1
       bfact1=Rproton/GeVm12fm*mproton
       bfact2=RA_common/GeVm12fm*mN
       bbx(1)=bupper
       bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       sub_num(3)=10
       init=1
    ENDIF
    IF(nuclearA_beam1.NE.0)THEN
       ! swap two beams
       x1_common=x2
       x2_common=x1
       E1_common=ebeam_PN(2)*x2
       E2_common=ebeam_PN(1)*x1
    ELSE
       x1_common=x1
       x2_common=x2
       E1_common=ebeam_PN(1)*x1
       E2_common=ebeam_PN(2)*x2
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       PhotonPhotonFlux_pA_hardsphere_eval=PhotonFlux_proton_EDFF_nob(x1_common,&
            gamma1_common)
       PhotonPhotonFlux_pA_hardsphere_eval=PhotonPhotonFlux_pA_hardsphere_eval*&
            PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
       PhotonPhotonFlux_pA_hardsphere_eval=MAX(PhotonPhotonFlux_pA_hardsphere_eval,0d0)
       RETURN
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact1*x1_common)
    ! b2*x2*mN = Exp(bA(2))
    aax(2)=DLOG(bfact2*x2_common)
    CALL ROMBERG_ND(PhotonPhotonFlux_pA_hardsphere_fxn,aax,bbx,3,sub_num,1,1d-5,&
         integral,ind,eval_num)
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          bbx(1)=bupper*2d0**(iter)
          bbx(2)=bupper*2d0**(iter)
          CALL ROMBERG_ND(PhotonPhotonFlux_pA_hardsphere_fxn,aax,bbx,3,sub_num,1,1d-5,&
               integral,ind,eval_num)
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       bbx(2)=bupper
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon flux at (x1,x2)=",x1_common,x2_common
          WRITE(*,*)"WARNING: use PNOHAD=1 approx. instead (most probably need to increase bupper)"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       PhotonPhotonFlux_pA_hardsphere_eval=PhotonFlux_proton_EDFF_nob(x1_common,&
            gamma1_common)
       PhotonPhotonFlux_pA_hardsphere_eval=PhotonPhotonFlux_pA_hardsphere_eval*&
            PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
       PhotonPhotonFlux_pA_hardsphere_eval=MAX(PhotonPhotonFlux_pA_hardsphere_eval,0d0)
    ELSE
       PhotonPhotonFlux_pA_hardsphere_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_pA_hardsphere_eval

  FUNCTION PhotonPhotonFlux_pA_hardsphere_fxn(dim_num,bA)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_pA_hardsphere_fxn
    INTEGER,INTENT(IN)::dim_num ! should be 3
    ! 1/0.1973d0 from fm to GeV-1 for b
    ! x1*b1*mproton=Exp(bA(1))
    ! x2*b2*mproton=Exp(bA(2))
    ! bA(3) = theta_{12}
    REAL(KIND(1d0)),DIMENSION(dim_num),INTENT(IN)::bA
    REAL(KIND(1d0))::b1,b2,b12,costh,pnohad
    REAL(KIND(1d0))::Ngamma1,Ngamma2
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0))::RA_common, A_common      ! radius of nuclei and atom number of nuclei 
    COMMON/PhotonPhoton_pA_HS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    IF(dim_num.NE.3)THEN
       WRITE(*,*)"ERROR: PhotonPhotonFlux_pA_hardsphere_fxn is not a three dimensional function"
       STOP 1
    ENDIF
    costh=DCOS(bA(3))
    ! in unit of GeV-1
    ! x1*b1*mproton=Exp(bA(1))
    b1=DEXP(bA(1))/x1_common/mproton
    ! x2*b2*mN=Exp(bA(2))
    b2=DEXP(bA(2))/x2_common/mN
    b12=DSQRT(b1**2+b2**2-2d0*b1*b2*costh)
    b12=b12*GeVm12fm ! from GeV-1 to fm
    pnohad=PNOHAD_pA_hardsphere(b12,0d0,RA_common,A_common,sigNN_inel_common)
    IF(pnohad.LE.0d0)THEN
       PhotonPhotonFlux_pA_hardsphere_fxn=0d0
       RETURN
    ENDIF
    Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
    Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    PhotonPhotonFlux_pA_hardsphere_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
    RETURN
  END FUNCTION PhotonPhotonFlux_pA_hardsphere_fxn

  FUNCTION PhotonPhotonFlux_pA_WoodsSaxon(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_pA_WoodsSaxon
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER::MX_save=0,MY_save=0
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD,ZD_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    INTEGER::gstate=0,gsel
    SAVE gstate
    REAL(KIND(1d0)),PARAMETER::ZUNSET=-HUGE(1d0) ! (M15) a node of the grid not computed yet
    LOGICAL::lazy
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_pA_WoodsSaxon=0d0
       RETURN
    ENDIF
    IF(gstate.NE.photon_final_state)THEN
       ! a new grid for another final state
       init=0
       MX_save=0
       MY_save=0
       gstate=photon_final_state
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(ffe_supported())THEN
       ! kernel-factorised engine (the correction beyond the proton dipole form factor
       ! is only available in the legacy integration)
       PhotonPhotonFlux_pA_WoodsSaxon=ffe_api(2,x1,x2,1,force_pnohad1,GENERATE_PhotonPhotonFlux_GRID)
       RETURN
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin).AND..NOT.force_pnohad1&
         .AND.GENERATE_PhotonPhotonFlux_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=0 ! (M16) as PhotonPhotonFlux_pp: nsegd nodes per decade
       DO I=0,log10xmin+1,-1
          MX=MX+nsegd(I)
       ENDDO
       MY=MX
       ! try to deallocate first
       IF(ALLOCATED(XD_1D))THEN
          DEALLOCATE(XD_1D)
       ENDIF
       ALLOCATE(XD_1D(MX+1))
       IF(ALLOCATED(YD_1D))THEN
          DEALLOCATE(YD_1D)
       ENDIF
       ALLOCATE(YD_1D(MY+1))
       IF(ALLOCATED(ZD))THEN
          DEALLOCATE(ZD)
       ENDIF
       ALLOCATE(ZD(MX+1,MY+1))
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nsegd(I)
             log10x1=-1d0/DBLE(nsegd(I))*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x1
             YD_1D(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP 1
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       ! (M15) states 1-3 (pp_grid_lazy): the nodes are computed when a lookup first needs them (fill_block), as pp
       lazy=pp_grid_lazy.AND.photon_final_state.NE.0
       IF(.NOT.lazy.AND.(photon_final_state.EQ.1.OR.photon_final_state.EQ.2))CALL prepare_block(MX+1,10d0**XD_1D,&
            MY+1,10d0**YD_1D,0d0)
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1.AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_save(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             IF(lazy)THEN
                ZD(I,J)=ZUNSET
                CYCLE
             ENDIF
             ZD(I,J)=PhotonPhotonFlux_pA_WoodsSaxon_eval(xx1,xx2)
             IF(MX_save.EQ.0.AND.MY_save.EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   WRITE(*,*)"INFO: generate grid of photon-photon flux in pA or Ap (will take a few minutes)"
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF(ALLOCATED(ZD_save))THEN
          DEALLOCATE(ZD_save)
       ENDIF
       ALLOCATE(ZD_save(MX+1,MY+1))
       MX_save=MX
       MY_save=MY
       DO I=1,MX+1
          DO J=1,MY+1
             ZD_save(I,J)=ZD(I,J)
          ENDDO
       ENDDO
       init=1
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(GENERATE_PhotonPhotonFlux_GRID)THEN
          XI(1)=DLOG10(x1)
          YI(1)=DLOG10(x2)
          IF(XI(1).GE.0d0)THEN
             ilog10x1=-1
          ELSE
             ilog10x1=FLOOR(XI(1))
          ENDIF
          IF(YI(1).GE.0d0)THEN
             ilog10x2=-1
          ELSE
             ilog10x2=FLOOR(YI(1))
          ENDIF
          !K=nseg*(-ilog10x1-1)
          !DO I=1,nseg+1
          !   XD2_1D(I)=XD_1D(K+I)
          !ENDDO
          !L=nseg*(-ilog10x2-1)
          !DO I=1,nseg+1
          !   YD2_1D(I)=YD_1D(L+I)
          !ENDDO
          !DO I=1,nseg+1
          !   DO J=1,nseg+1
          !      ZD2(I,J)=ZD(K+I,L+J)
          !   ENDDO
          !ENDDO
          !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
          DO I=1,MX_save+1 ! (M16) the first node <= log10 x1 (nodes decreasing)
             IF(XD_1D(I).LE.XI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save+1)THEN
             K=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             K=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save+1)THEN
             K=MX_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #1 !"
             STOP 1
          ENDIF
          DO I=1,MY_save+1
             IF(YD_1D(I).LE.YI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save+1)THEN
             L=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             L=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save+1)THEN
             L=MY_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #2 !"
             STOP 1
          ENDIF
          IF(pp_grid_lazy.AND.photon_final_state.NE.0)CALL fill_block(K,L)
          DO I=1,n_interp
             XD2_1D(I)=XD_1D(K+I)
             YD2_1D(I)=YD_1D(L+I)
          ENDDO
          DO I=1,n_interp
             DO J=1,n_interp
                ZD2(I,J)=ZD(K+I,L+J)
             ENDDO
          ENDDO
          CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
               YD2_1D,ZD2,1,XI,YI,ZI)
       ELSE
          ZI(1)=PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    IF((photon_final_state.EQ.1.OR.photon_final_state.EQ.2).AND..NOT.force_pnohad1)THEN
       ! the reference only checks the value (zero, NaN, sign): the elastic one (as PhotonPhotonFlux_pp)
       gsel=photon_final_state
       photon_final_state=0
       IF(ffe_supported().AND.GENERATE_PhotonPhotonFlux_GRID)THEN
          ! (M11) from FFE's tables, as PhotonPhotonFlux_pp
          pnohadval=ffe_api(2,x1,x2,1,.TRUE.,.TRUE.)
       ELSE
          pnohadval=PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2,.TRUE.)
       ENDIF
       photon_final_state=gsel
    ELSE
       pnohadval=PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2,.TRUE.)
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pA_WoodsSaxon=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.&
            (DABS(ZI(1)/pnohadval).GT.1D2.AND.photon_final_state.NE.1.AND.photon_final_state.NE.2))THEN
          PhotonPhotonFlux_pA_WoodsSaxon=pnohadval
       ELSE
          PhotonPhotonFlux_pA_WoodsSaxon=ZI(1)
       ENDIF
    ELSE
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pA_WoodsSaxon=0d0
       ELSE
          PhotonPhotonFlux_pA_WoodsSaxon=pnohadval
       ENDIF
    ENDIF
    RETURN
  CONTAINS
    INTEGER FUNCTION nsegd(id) ! (M16) as in PhotonPhotonFlux_pp
      INTEGER,INTENT(IN)::id
      nsegd=nseg
      IF(pp_grid_fine.AND.photon_final_state.NE.0)THEN
         IF(id.EQ.0)nsegd=40
         IF(id.EQ.-1)nsegd=20
      ENDIF
    END FUNCTION nsegd
    SUBROUTINE fill_block(K0,L0)
      ! (M15) as PhotonPhotonFlux_pp: the nodes of the stencil K0+1..K0+n_interp, L0+1..L0+n_interp not computed yet,
      ! the transforms (or the b-space rows) of their x and scales first
      INTEGER,INTENT(IN)::K0,L0
      REAL(KIND(1d0)),DIMENSION(n_interp)::xa,xb
      INTEGER::ii,jj
      IF(ALL(ZD(K0+1:K0+n_interp,L0+1:L0+n_interp).NE.ZUNSET))RETURN
      xa=10d0**XD_1D(K0+1:K0+n_interp)
      xb=10d0**YD_1D(L0+1:L0+n_interp)
      IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)CALL prepare_block(n_interp,xa,n_interp,xb,&
           photon_mu2(xa(n_interp),xb(n_interp))*(1d0-1d-12))
      DO jj=1,n_interp
         DO ii=1,n_interp
            IF(ZD(K0+ii,L0+jj).NE.ZUNSET)CYCLE
            ZD(K0+ii,L0+jj)=PhotonPhotonFlux_pA_WoodsSaxon_eval(xa(ii),xb(jj))
            ZD_save(K0+ii,L0+jj)=ZD(K0+ii,L0+jj)
         ENDDO
      ENDDO
    END SUBROUTINE fill_block
    SUBROUTINE prepare_block(na,xa,nb,xb,mu2lo)
      ! the transforms (momentum space) or the b-space rows of the proton (legacy) for the beam-1 x of xa and the
      ! beam-2 x of xb: momentum space at the scales >= mu2lo of the grid rows and columns through them, legacy at
      ! the scales of the pairs xa x xb
      INTEGER,INTENT(IN)::na,nb
      REAL(KIND(1d0)),DIMENSION(na),INTENT(IN)::xa
      REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::xb
      REAL(KIND(1d0)),INTENT(IN)::mu2lo
      REAL(KIND(1d0))::d
      IF(qspace_ok())THEN
         ! the nucleus of _eval, set at its first call with 0 < x < 1 (a stencil may start at x = 1)
         d=PhotonPhotonFlux_pA_WoodsSaxon_eval(1d-3,1d-3,.TRUE.)
         CALL pA_qspace_setup_common()
         IF(nuclearA_beam1.NE.0)THEN
            CALL pA_qspace_prepare_grid(nb,xb,MX+1,10d0**XD_1D,na,xa,mu2lo)
         ELSE
            CALL pA_qspace_prepare_grid(na,xa,MY+1,10d0**YD_1D,nb,xb,mu2lo)
         ENDIF
      ELSEIF(nuclearA_beam1.NE.0)THEN
         CALL pp_inelastic_prepare(nb,xb,na,xa,.TRUE.)
      ELSE
         CALL pp_inelastic_prepare(na,xa,nb,xb,.TRUE.)
      ENDIF
    END SUBROUTINE prepare_block
  END FUNCTION PhotonPhotonFlux_pA_WoodsSaxon

  FUNCTION PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_pA_WoodsSaxon_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0))::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    REAL(KIND(1d0))::aaVal_common,wVal_common ! parameters in Woods-Saxon potential
    COMMON/PhotonPhoton_pA_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,aaVal_common,wVal_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact1,bfact2
    SAVE bfact1,bfact2
    REAL(KIND(1d0)),PARAMETER::bupper=2d0
    CHARACTER(len=7)::Aname
    REAL(KIND(1d0))::cmenergy
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    REAL(KIND(1d0))::elp,elA
    INTEGER,DIMENSION(3)::subn
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_pA_WoodsSaxon_eval=0d0
       RETURN
    ENDIF
    IF(photon_final_state.NE.0.AND.ANY(neutron_tagging.NE.-2))THEN
       WRITE(*,*)"ERROR: photon_final_state /= 0 is not available with neutron tagging"
       STOP
    ENDIF
    CALC_DELTABQ=.FALSE.
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_common=ebeam_PN(2)/mproton
          gamma2_common=ebeam_PN(1)/mN
       ELSE
          gamma1_common=ebeam_PN(1)/mproton
          gamma2_common=ebeam_PN(2)/mN
       ENDIF
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! proton Z is 1
       Z1=1d0
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       IF(nuclearA_beam1.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       ELSEIF(nuclearA_beam2.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
       ELSE
          WRITE(*,*)"ERROR: please set nuclearA_beam1/nuclearZ_beam1 or nuclearA_beam2/nuclearZ_beam2 nonzero first !"
          STOP 1
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in p+"//TRIM(Aname)//" collisions"
       CALL GetNuclearInfo(Aname,A_common,Z2,RA_common,aaVal_common,wVal_common)
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
       sigNN_inel_common=sigma_inelastic(cmenergy)
       sigNN_inel_common=sigNN_inel_common*0.1d0 ! from mb to fm^2
       ! 0.1973 is from fm to GeV-1
       bfact1=Rproton/GeVm12fm*mproton
       IF(IEPSILON_EDFF.EQ.0)THEN
          bfact2=RA_common/GeVm12fm*mN
       ELSE
          bfact2=(RA_common+aaVal_common)/GeVm12fm*mN
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          bfact1=bfact1*LOWER_BFactor_Limit
          bfact2=bfact2*LOWER_BFactor_Limit
       ENDIF
       bbx(1)=bupper
       bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we should increase the number of segments
          sub_num(1)=sub_num(1)*SUB_FACTOR
          sub_num(2)=sub_num(2)*SUB_FACTOR
       ENDIF
       sub_num(3)=10
       init=1
    ENDIF
    IF(nuclearA_beam1.NE.0)THEN
       ! swap two beams
       x1_common=x2
       x2_common=x1
       E1_common=ebeam_PN(2)*x2
       E2_common=ebeam_PN(1)*x1
    ELSE
       x1_common=x1
       x2_common=x2
       E1_common=ebeam_PN(1)*x1
       E2_common=ebeam_PN(2)*x2
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(ffe_supported())THEN
       ! kernel-factorised engine (the correction beyond the proton dipole form factor
       ! is only available in the legacy integration)
       PhotonPhotonFlux_pA_WoodsSaxon_eval=ffe_api(2,x1,x2,1,force_pnohad1,.FALSE.)
       RETURN
    ENDIF
    CALL set_photon_scale(x1,x2)
    IF(force_pnohad1.OR.photon_final_state.EQ.3)THEN
       ! we only use PNOHAD=1
       IF(photon_final_state.NE.0.AND.qspace_ok())THEN
          ! (M15) the final states with the momentum-space path: the exact collinear flux of the nucleus
          elp=PhotonFlux_proton_nob(x1_common,gamma1_common)
          CALL pA_qspace_setup(RA_common,wVal_common,aaVal_common,A_common,sigNN_inel_common)
          elA=Z2**2*alpha/x2_common*pA_nucleus_coll(x2_common)
       ELSEIF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          elp=PhotonFlux_proton_nob(x1_common,gamma1_common)
          elA=PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common,wVal_common,aaVal_common)
       ELSE
          elp=PhotonFlux_proton_EDFF_nob(x1_common,gamma1_common)
          elA=PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
       ENDIF
       ! (only the proton, beam 1 here, has an inelastic photon)
       PhotonPhotonFlux_pA_WoodsSaxon_eval=final_state_nob(elp*ael_photon,elA*ael_photon,x1_common,x2_common,1)
       PhotonPhotonFlux_pA_WoodsSaxon_eval=MAX(PhotonPhotonFlux_pA_WoodsSaxon_eval,0d0)
       RETURN
    ENDIF
    IF((photon_final_state.EQ.1.OR.photon_final_state.EQ.2).AND.qspace_ok())THEN
       ! (M15) momentum space, as pp (the proton is beam 1 here): 1: inel-el P = inel-el - inel-el (1-P), 2: (1-P)
       ! (el-el + inel-el)
       CALL pA_qspace_setup(RA_common,wVal_common,aaVal_common,A_common,sigNN_inel_common)
       CALL pA_qspace_L(x1_common,x2_common,elp,elA) ! elp = L(el_p,el_A), elA = L(in_p,el_A) (unit charges, alpha(0))
       IF(photon_final_state.EQ.1)THEN
          integral=-alpha**2*Z1**2*Z2**2/(x1*x2)*ael_photon*elA
          elA=Z2**2*alpha/x2_common*pA_nucleus_coll(x2_common) ! the exact collinear flux, as state 3
          PhotonPhotonFlux_pA_WoodsSaxon_eval=integral+final_state_nob(0d0,elA*ael_photon,x1_common,x2_common,1)
       ELSE
          PhotonPhotonFlux_pA_WoodsSaxon_eval=alpha**2*Z1**2*Z2**2/(x1*x2)*(ael_photon**2*elp+ael_photon*elA)
       ENDIF
       RETURN
    ENDIF
    IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
       CALL pp_inelastic_prepare(1,(/x1_common/),1,(/x2_common/),.TRUE.)
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact1*x1_common)
    IF(USE_CHARGEFORMFACTOR4PHOTON.AND.USE_EXACT_ELASTIC4PHOTON)aax(1)=aax(1)+DLOG(EXACT_BLOW) ! as pp
    ! b2*x2*mN = Exp(bA(2))
    aax(2)=DLOG(bfact2*x2_common)
    ! we should also choose the upper limit dynamically for the ion beam
    bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common/GeVm12fm*DSQRT(x2_common*mN)))
    subn=sub_num
    IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
       ! the proton b range down to the b of its inelastic photons at large k_T (as pp)
       CALL pp_inelastic_range(x1_common,aax(1),subn(1))
    ENDIF
    CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax,bbx,3,subn,1,1d-5,&
         integral,ind,eval_num)
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0.AND.photon_final_state.NE.1)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          !bbx(1)=bupper*2d0**(iter)
          !bbx(2)=bupper*2d0**(iter)
          bbx(1)=bbx(1)*2d0
          bbx(2)=bbx(2)*2d0
          CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax,bbx,3,subn,1,1d-5,&
               integral,ind,eval_num)
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       !bbx(2)=bupper
       bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common/GeVm12fm*DSQRT(x2_common*mN)))
    ENDIF
    IF(integral.LT.0d0.AND.photon_final_state.NE.1)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon flux at (x1,x2)=",x1_common,x2_common
          WRITE(*,*)"WARNING: use PNOHAD=1 approx. instead (most probably need to increase bupper)"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for approximation
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonFlux_proton_nob(x1_common,gamma1_common)
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval*&
            PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common,wVal_common,aaVal_common)
       ELSE
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonFlux_proton_EDFF_nob(x1_common,&
               gamma1_common)
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval*&
               PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
       ENDIF
       PhotonPhotonFlux_pA_WoodsSaxon_eval=MAX(PhotonPhotonFlux_pA_WoodsSaxon_eval,0d0)*ael_photon**2
    ELSE
       PhotonPhotonFlux_pA_WoodsSaxon_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
       IF(photon_final_state.EQ.1)THEN
          ! 1: inel-el P = inel-el - inel-el (1-P), the first from the collinear fluxes, the second integrated above
          ! (as pp); 2: no inel-inel in pA (the nucleus is elastic only)
          IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
             elA=PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common,wVal_common,aaVal_common)
          ELSE
             elA=PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
          ENDIF
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval+&
               final_state_nob(0d0,elA*ael_photon,x1_common,x2_common,1)
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_pA_WoodsSaxon_eval

  FUNCTION DeltaB_pA_at_x1x2(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::DeltaB_pA_at_x1x2
    REAL(KIND(1d0))::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    LOGICAL::force_pnohad1
    REAL(KIND(1d0))::num,den
    LOGICAL::use_grid_bu
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    ! numerator
    ! PhotonPhotonFlux*Delta b in unit of fm
    num=PhotonPhotonDeltaB_pA_WoodsSaxon(x1,x2,force_pnohad1)
    ! denominator
    ! PhotonPhotonFlux
    use_grid_bu=GENERATE_PhotonPhotonFlux_GRID
    ! Let us make it to be coherent with DELTAB grid case
    GENERATE_PhotonPhotonFlux_GRID=GENERATE_DELTAB_GRID
    den=PhotonPhotonFlux_pA_WoodsSaxon(x1,x2,force_pnohad1)
    GENERATE_PhotonPhotonFlux_GRID=use_grid_bu ! recover the original one
    IF(den.LE.0d0.OR.num.EQ.0d0)THEN
       DeltaB_pA_at_x1x2=0d0
    ELSE
       DeltaB_pA_at_x1x2=num/den
    ENDIF
    RETURN
  END FUNCTION DeltaB_pA_at_x1x2

  FUNCTION PhotonPhotonDeltaB_pA_WoodsSaxon(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonDeltaB_pA_WoodsSaxon
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER,DIMENSION(2)::MX_save=(/0,0/),MY_save=(/0,0/)
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::ZD
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD1_save,ZD2_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD1_save,ZD2_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    LOGICAL,DIMENSION(2)::gridready=(/.FALSE.,.FALSE./)
    SAVE gridready
    INTEGER::igrid
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonDeltaB_pA_WoodsSaxon=0d0
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(ffe_supported())THEN
       ! kernel-factorised engine (the correction beyond the proton dipole form factor
       ! is only available in the legacy integration)
       PhotonPhotonDeltaB_pA_WoodsSaxon=ffe_api(2,x1,x2,2,force_pnohad1,GENERATE_DELTAB_GRID)
       RETURN
    ENDIF
    IF(force_pnohad1)THEN
       igrid=2
    ELSE
       igrid=1
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin.OR..NOT.gridready(igrid)).AND.GENERATE_DELTAB_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=nseg*(-log10xmin)
       MY=MX
       ! try to deallocate first
       IF(init.EQ.0.OR.log10xmin_before.GT.log10xmin)THEN
          IF(ALLOCATED(XD_1D))THEN
             DEALLOCATE(XD_1D)
          ENDIF
          ALLOCATE(XD_1D(MX+1))
          IF(ALLOCATED(YD_1D))THEN
             DEALLOCATE(YD_1D)
          ENDIF
          ALLOCATE(YD_1D(MY+1))
          IF(ALLOCATED(ZD))THEN
             WRITE(*,*)"INFO: in PhotonPhotonDeltaB_pA_WoodsSaxon, the xmin of the grid has been updated to ",xmin
             ! save the values that have been calculated before
             IF(gridready(1))THEN
                ! the saved grid has the size before the xmin update
                MX_save(1)=nseg*(-log10xmin_before)
                MY_save(1)=MX_save(1)
                IF(ALLOCATED(ZD1_save))THEN
                   DEALLOCATE(ZD1_save)
                ENDIF
                ALLOCATE(ZD1_save(MX_save(1)+1,MY_save(1)+1))
                DO I=1,MX_save(1)+1
                   DO J=1,MY_save(1)+1
                      ZD1_save(I,J)=ZD(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(gridready(2))THEN
                ! the saved grid has the size before the xmin update
                MX_save(2)=nseg*(-log10xmin_before)
                MY_save(2)=MX_save(2)
                IF(ALLOCATED(ZD2_save))THEN
                   DEALLOCATE(ZD2_save)
                ENDIF
                ALLOCATE(ZD2_save(MX_save(2)+1,MY_save(2)+1))
                DO I=1,MX_save(2)+1
                   DO J=1,MY_save(2)+1
                      ZD2_save(I,J)=ZD(2,I,J)
                   ENDDO
                ENDDO
             ENDIF
             DEALLOCATE(ZD)
          ENDIF
          ALLOCATE(ZD(2,MX+1,MY+1))
          K=0
          DO I=0,log10xmin+1,-1
             DO J=1,nseg
                log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
                K=K+1
                XD_1D(K)=log10x1
                YD_1D(K)=log10x1
             ENDDO
          ENDDO
          IF(K.NE.MX)THEN
             WRITE(*,*)"ERROR: K != MX"
             STOP 1
          ENDIF
          XD_1D(MX+1)=DBLE(log10xmin)
          YD_1D(MY+1)=DBLE(log10xmin)
          gridready(1:2)=.FALSE.
       ENDIF
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! save the values that have been calculated before
             IF(MX_save(igrid).GT.0.AND.MY_save(igrid).GT.0.AND.I.LE.MX_save(igrid)+1&
                  .AND.J.LE.MY_save(igrid)+1)THEN
                IF(igrid.EQ.1)THEN
                   ZD(igrid,I,J)=ZD1_save(I,J)
                ELSEIF(igrid.EQ.2)THEN
                   ZD(igrid,I,J)=ZD2_save(I,J)
                ELSE
                   WRITE(*,*)"ERROR: do not know igrid=",igrid
                ENDIF
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             ZD(igrid,I,J)=PhotonPhotonDeltaB_pA_WoodsSaxon_eval(xx1,xx2,force_pnohad1)
             IF(MX_save(igrid).EQ.0.AND.MY_save(igrid).EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   IF(igrid.EQ.1)THEN
                      WRITE(*,*)"INFO: generate grid of photon-photon <delta b> in pA or Ap (will take a few minutes)"
                   ELSE
                      WRITE(*,*)"INFO: generate grid of photon-photon <delta b> in pA or Ap with PNOHAD=1 (will take a few minutes)"
                   ENDIF
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       gridready(igrid)=.TRUE.
       init=1
    ENDIF
    IF(GENERATE_DELTAB_GRID)THEN
       XI(1)=DLOG10(x1)
       YI(1)=DLOG10(x2)
       IF(XI(1).GE.0d0)THEN
          ilog10x1=-1
       ELSE
          ilog10x1=FLOOR(XI(1))
       ENDIF
       IF(YI(1).GE.0d0)THEN
          ilog10x2=-1
       ELSE
          ilog10x2=FLOOR(YI(1))
       ENDIF
       !K=nseg*(-ilog10x1-1)
       !DO I=1,nseg+1
       !   XD2_1D(I)=XD_1D(K+I)
       !ENDDO
       !L=nseg*(-ilog10x2-1)
       !DO I=1,nseg+1
       !   YD2_1D(I)=YD_1D(L+I)
       !ENDDO
       !DO I=1,nseg+1
       !   DO J=1,nseg+1
       !      ZD2(I,J)=ZD(igrid,K+I,L+J)
       !   ENDDO
       !ENDDO
       !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
       K=nseg*(-ilog10x1-1)
       DO I=K+1,K+nseg+1
          if(XD_1D(I).LE.XI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX+1)THEN
          K=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MX+1)THEN
          K=MX+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #1 !"
          STOP 1
       ENDIF
       L=nseg*(-ilog10x2-1)
       DO I=L+1,L+nseg+1
          if(YD_1D(I).LE.YI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY+1)THEN
          L=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          L=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MY+1)THEN
          L=MY+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #2 !"
          STOP 1
       ENDIF
       DO I=1,n_interp
          XD2_1D(I)=XD_1D(K+I)
          YD2_1D(I)=YD_1D(L+I)
       ENDDO
       DO I=1,n_interp
          DO J=1,n_interp
             ZD2(I,J)=ZD(igrid,K+I,L+J)
          ENDDO
       ENDDO
       CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
            YD2_1D,ZD2,1,XI,YI,ZI)
       IF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0)THEN
          PhotonPhotonDeltaB_pA_WoodsSaxon=0d0
       ELSE
          PhotonPhotonDeltaB_pA_WoodsSaxon=ZI(1)
       ENDIF
    ELSE
       PhotonPhotonDeltaB_pA_WoodsSaxon=PhotonPhotonDeltaB_pA_WoodsSaxon_eval(x1,x2,force_pnohad1)
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonDeltaB_pA_WoodsSaxon

  FUNCTION PhotonPhotonDeltaB_pA_WoodsSaxon_eval(x1,x2,FORCEPNOHAD1)
    ! this is PhotonPhotonFlux*Delta b, where Delta b is in unit of fm
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonDeltaB_pA_WoodsSaxon_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0))::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    REAL(KIND(1d0))::aaVal_common,wVal_common ! parameters in Woods-Saxon potential
    COMMON/PhotonPhoton_pA_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,aaVal_common,wVal_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact1,bfact2
    SAVE bfact1,bfact2
    REAL(KIND(1d0)),PARAMETER::bupper=2d0
    CHARACTER(len=7)::Aname
    REAL(KIND(1d0))::cmenergy
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonDeltaB_pA_WoodsSaxon_eval=0d0
       RETURN
    ENDIF
    CALC_DELTABQ=.TRUE.
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_common=ebeam_PN(2)/mproton
          gamma2_common=ebeam_PN(1)/mN
       ELSE
          gamma1_common=ebeam_PN(1)/mproton
          gamma2_common=ebeam_PN(2)/mN
       ENDIF
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! proton Z is 1
       Z1=1d0
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       IF(nuclearA_beam1.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       ELSEIF(nuclearA_beam2.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
       ELSE
          WRITE(*,*)"ERROR: please set nuclearA_beam1/nuclearZ_beam1 or nuclearA_beam2/nuclearZ_beam2 nonzero first !"
          STOP 1
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in p+"//TRIM(Aname)//" collisions"
       CALL GetNuclearInfo(Aname,A_common,Z2,RA_common,aaVal_common,wVal_common)
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
       sigNN_inel_common=sigma_inelastic(cmenergy)
       sigNN_inel_common=sigNN_inel_common*0.1d0 ! from mb to fm^2
       ! 0.1973 is from fm to GeV-1
       bfact1=Rproton/GeVm12fm*mproton
       IF(IEPSILON_EDFF.EQ.0)THEN
          bfact2=RA_common/GeVm12fm*mN
       ELSE
          bfact2=(RA_common+aaVal_common)/GeVm12fm*mN
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          bfact1=bfact1*LOWER_BFactor_Limit
          bfact2=bfact2*LOWER_BFactor_Limit
       ENDIF
       bbx(1)=bupper
       bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we should increase the number of segments
          sub_num(1)=sub_num(1)*SUB_FACTOR
          sub_num(2)=sub_num(2)*SUB_FACTOR
       ENDIF
       sub_num(3)=10
       init=1
    ENDIF
    IF(nuclearA_beam1.NE.0)THEN
       ! swap two beams
       x1_common=x2
       x2_common=x1
       E1_common=ebeam_PN(2)*x2
       E2_common=ebeam_PN(1)*x1
    ELSE
       x1_common=x1
       x2_common=x2
       E1_common=ebeam_PN(1)*x1
       E2_common=ebeam_PN(2)*x2
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(ffe_supported())THEN
       ! kernel-factorised engine (the correction beyond the proton dipole form factor
       ! is only available in the legacy integration)
       PhotonPhotonDeltaB_pA_WoodsSaxon_eval=ffe_api(2,x1,x2,2,force_pnohad1,.FALSE.)
       RETURN
    ENDIF
    CALL set_photon_scale(x1,x2)
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact1*x1_common)
    ! b2*x2*mN = Exp(bA(2))
    aax(2)=DLOG(bfact2*x2_common)
    ! we should also choose the upper limit dynamically for the ion beam
    bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common/GeVm12fm*DSQRT(x2_common*mN)))
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax(1:2),bbx(1:2),2,sub_num(1:2),1,1d-5,&
            integral,ind,eval_num)
    ELSE
       CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
            integral,ind,eval_num)
    ENDIF
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          !bbx(1)=bupper*2d0**(iter)
          !bbx(2)=bupper*2d0**(iter)
          bbx(1)=bbx(1)*2d0
          bbx(2)=bbx(2)*2d0
          IF(force_pnohad1)THEN
             ! we only use PNOHAD=1
             CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax(1:2),bbx(1:2),2,sub_num(1:2),1,1d-5,&
                  integral,ind,eval_num)
          ELSE
             CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
                  integral,ind,eval_num)
          ENDIF
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       !bbx(2)=bupper
       bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common/GeVm12fm*DSQRT(x2_common*mN)))
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon <delta b> at (x1,x2)=",x1_common,x2_common
          WRITE(*,*)"WARNING: most probably need to increase bupper"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       PhotonPhotonDeltaB_pA_WoodsSaxon_eval=0d0
    ELSE
       PhotonPhotonDeltaB_pA_WoodsSaxon_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonDeltaB_pA_WoodsSaxon_eval

  FUNCTION PhotonPhotonFlux_pA_WoodsSaxon_fxn(dim_num,bA)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_pA_WoodsSaxon_fxn
    INTEGER,INTENT(IN)::dim_num ! should be 3 or 2
    ! 1/0.1973d0 from fm to GeV-1 for b
    ! x1*b1*mproton=Exp(bA(1))
    ! x2*b2*mproton=Exp(bA(2))
    ! bA(3) = theta_{12}
    REAL(KIND(1d0)),DIMENSION(dim_num),INTENT(IN)::bA
    REAL(KIND(1d0))::b1,b2,b12,costh,pnohad
    REAL(KIND(1d0))::b12GeVm1,pbreak
    REAL(KIND(1d0))::Ngamma1,Ngamma2,Ninel1,blo_el
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! as PhotonPhotonFlux_pA_WoodsSaxon_eval
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0))::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    REAL(KIND(1d0))::aaVal_common,wVal_common ! parameters in Woods-Saxon potential
    COMMON/PhotonPhoton_pA_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,aaVal_common,wVal_common
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::RR,aaa
    REAL(KIND(1d0)),EXTERNAL::elliptic_em
    REAL(KIND(1d0))::arg1,arg2
    IF(dim_num.NE.3.AND.dim_num.NE.2)THEN
       WRITE(*,*)"ERROR: PhotonPhotonFlux_pA_WoodsSaxon_fxn is not a three or two dimensional function"
       STOP 1
    ENDIF
    ! in unit of GeV-1
    ! x1*b1*mproton=Exp(bA(1))
    b1=DEXP(bA(1))/x1_common/mproton
    ! x2*b2*mN=Exp(bA(2))
    b2=DEXP(bA(2))/x2_common/mN
    IF(dim_num.EQ.3)THEN
       costh=DCOS(bA(3))
       b12=DSQRT(b1**2+b2**2-2d0*b1*b2*costh)
       b12GeVm1=b12 ! in unit of GeV-1
       b12=b12*GeVm12fm ! from GeV-1 to fm
       pnohad=PNOHAD_pA_WoodsSaxon(b12,0d0,RA_common,wVal_common,aaVal_common,&
            A_common,sigNN_inel_common)
       IF(nuclearA_beam1.NE.0)THEN
          ! A+p
          IF(neutron_tagging(1).NE.-2)THEN
             ! multiply the forward neutron tagging probability
             pbreak=Pbreak_pA_WoodsSaxon(neutron_tagging(1),b12GeVm1)
             pnohad=pnohad*pbreak
          ENDIF
       ELSE
          ! p+A
          IF(neutron_tagging(2).NE.-2)THEN
             ! multiply the forward neutron tagging probability
             pbreak=Pbreak_pA_WoodsSaxon(neutron_tagging(2),b12GeVm1)
             pnohad=pnohad*pbreak
          ENDIF
       ENDIF
       ! hadronic interaction (1 and 2: the parts with 1-P, see PhotonPhotonFlux_pA_WoodsSaxon_eval)
       IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)pnohad=1d0-pnohad
       IF(CALC_DELTABQ)THEN
          pnohad=pnohad*b12
       ENDIF
    ELSE
       ! dim_num = 2
       ! it must be a case with pnohad=1
       ! and CALC_DELTABQ=.TRUE
       IF(.NOT.CALC_DELTABQ)THEN
          WRITE(*,*)"ERROR: cannot reach here with CALC_DELTABQ=.FALSE."
       ENDIF
       ! pnohad is Integrate[Sqrt[b1^2 + b2^2 - 2 b1*b2*Cos[phi]], {phi, 0, 2 Pi}]
       arg2=4d0*b1*b2/(b1+b2)**2
       IF(b1.NE.b2.AND.arg2.LT.1d0)THEN
          arg1=-4d0*b1*b2/(b1-b2)**2
          ! elliptic_em(m) is EllipticE[m] with m<1
          pnohad=DABS(b1-b2)*elliptic_em(arg1)
          pnohad=pnohad+(b1+b2)*elliptic_em(arg2)
       ELSE
          pnohad=4d0*b1
       ENDIF
       pnohad=pnohad*2d0*GeVm12fm ! convert from GeV-1 to fm
    ENDIF
    IF(pnohad.LE.0d0)THEN
       PhotonPhotonFlux_pA_WoodsSaxon_fxn=0d0
       RETURN
    ENDIF
    blo_el=0d0
    IF(photon_final_state.EQ.1.OR.photon_final_state.EQ.2)THEN
       ! the elastic proton photon on the b range of photon_final_state=0 (as pp); the range is extended for the
       ! inelastic photon only (pp_inelastic_range)
       blo_el=Rproton/GeVm12fm
       IF(USE_CHARGEFORMFACTOR4PHOTON)blo_el=blo_el*LOWER_BFactor_Limit
       IF(USE_CHARGEFORMFACTOR4PHOTON.AND.USE_EXACT_ELASTIC4PHOTON)blo_el=blo_el*EXACT_BLOW
    ENDIF
    IF(.NOT.USE_CHARGEFORMFACTOR4PHOTON)THEN
       Ngamma1=0d0
       IF(b1.GE.blo_el)Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
       Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    ELSE
       ! converting from fm to GeV-1
       RR=RA_common/GeVm12fm
       aaa=aaVal_common/GeVm12fm
       ! the proton is always 1 and the nucleus 2 here, whatever the beam order: the callers (_eval) already
       ! swapped x, E and gamma when the nucleus is beam 1 (v2.0, fix of the A+p branch). The nucleus uses the
       ! grid 2, as PhotonFlux_nucleus_ChFF_nob(2,...) in _eval.
       Ngamma1=0d0
       IF(b1.GE.blo_el)Ngamma1=PhotonNumberDensity_proton(b1,E1_common,gamma1_common,1)
       Ngamma2=PhotonNumberDensity_ChargeFormFactor_WS(b2,E2_common,gamma2_common,&
            RR,wVal_common,aaa,3d0,0.7d0,2)
    ENDIF
    ! the elastic photons normalised to alpha(mu2) with running alpha (set_photon_scale; the nucleus photon has no
    ! running alpha(Q2), Q2 <~ 1/R_A**2, only this normalisation)
    Ngamma1=Ngamma1*ael_photon
    Ngamma2=Ngamma2*ael_photon
    IF(photon_final_state.EQ.0)THEN
       PhotonPhotonFlux_pA_WoodsSaxon_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
    ELSE
       Ninel1=PhotonNumberDensity_inelastic_beam(b1,E1_common,gamma1_common,1)
       IF(photon_final_state.EQ.1)THEN
          ! -(1-P) inel-el; inel-el is added in PhotonPhotonFlux_pA_WoodsSaxon_eval
          PhotonPhotonFlux_pA_WoodsSaxon_fxn=-b1**2*b2**2*pnohad*Ninel1*Ngamma2
       ELSE
          ! 2: (1-P) (el-el + inel-el)
          PhotonPhotonFlux_pA_WoodsSaxon_fxn=b1**2*b2**2*pnohad*(Ngamma1+Ninel1)*Ngamma2
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_pA_WoodsSaxon_fxn

  FUNCTION PhotonPhotonFlux_AB_hardsphere(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_AB_hardsphere
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER::MX_save=0,MY_save=0
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD,ZD_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_AB_hardsphere=0d0
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin).AND..NOT.force_pnohad1&
         .AND.GENERATE_PhotonPhotonFlux_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=nseg*(-log10xmin)
       MY=MX
       ! try to deallocate first
       IF(ALLOCATED(XD_1D))THEN
          DEALLOCATE(XD_1D)
       ENDIF
       ALLOCATE(XD_1D(MX+1))
       IF(ALLOCATED(YD_1D))THEN
          DEALLOCATE(YD_1D)
       ENDIF
       ALLOCATE(YD_1D(MY+1))
       IF(ALLOCATED(ZD))THEN
          WRITE(*,*)"INFO: in PhotonPhotonFlux_AB_hardsphere, the xmin of the grid has been updated to ",xmin
          DEALLOCATE(ZD)
       ENDIF
       ALLOCATE(ZD(MX+1,MY+1))
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nseg
             log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x1
             YD_1D(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP 1
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! save the values that have been calculated before
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1.AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_SAVE(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             ZD(I,J)=PhotonPhotonFlux_AB_hardsphere_eval(xx1,xx2)
             IF(MX_SAVE.EQ.0.AND.MY_SAVE.EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   WRITE(*,*)"INFO: generate grid of photon-photon flux in AB (will take tens of seconds)"
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF(ALLOCATED(ZD_save))THEN
          DEALLOCATE(ZD_save)
       ENDIF
       ALLOCATE(ZD_save(MX+1,MY+1))
       ! save the values that have been calculated before
       MX_save=MX
       MY_save=MY
       DO I=1,MX+1
          DO J=1,MY+1
             ZD_save(I,J)=ZD(I,J)
          ENDDO
       ENDDO
       init=1
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(GENERATE_PhotonPhotonFlux_GRID)THEN
          XI(1)=DLOG10(x1)
          YI(1)=DLOG10(x2)
          IF(XI(1).GE.0d0)THEN
             ilog10x1=-1
          ELSE
             ilog10x1=FLOOR(XI(1))
          ENDIF
          IF(YI(1).GE.0d0)THEN
             ilog10x2=-1
          ELSE
             ilog10x2=FLOOR(YI(1))
          ENDIF
          !K=nseg*(-ilog10x1-1)
          !DO I=1,nseg+1
          !   XD2_1D(I)=XD_1D(K+I)
          !ENDDO
          !L=nseg*(-ilog10x2-1)
          !DO I=1,nseg+1
          !   YD2_1D(I)=YD_1D(L+I)
          !ENDDO
          !DO I=1,nseg+1
          !   DO J=1,nseg+1
          !      ZD2(I,J)=ZD(K+I,L+J)
          !   ENDDO
          !ENDDO
          !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
          K=nseg*(-ilog10x1-1)
          DO I=K+1,K+nseg+1
             if(XD_1D(I).LE.XI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save+1)THEN
             K=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             K=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save+1)THEN
             K=MX_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #1 !"
             STOP 1
          ENDIF
          L=nseg*(-ilog10x2-1)
          DO I=L+1,L+nseg+1
             if(YD_1D(I).LE.YI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save+1)THEN
             L=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             L=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save+1)THEN
             L=MY_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #2 !"
             STOP 1
          ENDIF
          DO I=1,n_interp
             XD2_1D(I)=XD_1D(K+I)
             YD2_1D(I)=YD_1D(L+I)
          ENDDO
          DO I=1,n_interp
             DO J=1,n_interp
                ZD2(I,J)=ZD(K+I,L+J)
             ENDDO
          ENDDO
          CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
               YD2_1D,ZD2,1,XI,YI,ZI)
       ELSE
          ZI(1)=PhotonPhotonFlux_AB_hardsphere_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    pnohadval=PhotonPhotonFlux_AB_hardsphere_eval(x1,x2,.TRUE.)
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_AB_hardsphere=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.DABS(ZI(1)/pnohadval).GT.1D2)THEN
          PhotonPhotonFlux_AB_hardsphere=pnohadval
       ELSE
          PhotonPhotonFlux_AB_hardsphere=ZI(1)
       ENDIF
    ELSE
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_AB_hardsphere=0d0
       ELSE
          PhotonPhotonFlux_AB_hardsphere=pnohadval
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_AB_hardsphere

  FUNCTION PhotonPhotonFlux_AB_hardsphere_eval(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_AB_hardsphere_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, only evaluate with PNOHAD=1
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0)),DIMENSION(2)::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    COMMON/PhotonPhoton_AB_HS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact1,bfact2
    SAVE bfact1,bfact2
    REAL(KIND(1d0)),PARAMETER::bupper=3d0
    REAL(KIND(1d0))::aaVal,wVal
    CHARACTER(len=7)::Aname1,Aname2
    REAL(KIND(1d0))::cmenergy
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_AB_hardsphere_eval=0d0
       RETURN
    ENDIF
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP 1
       ENDIf
       gamma1_common=ebeam_PN(1)/mN
       gamma2_common=ebeam_PN(2)/mN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       Aname1=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       CALL GetNuclearInfo(Aname1,A_common(1),Z1,RA_common(1),aaval,wval)
       IF(nuclearA_beam2.NE.nuclearA_beam1.OR.nuclearZ_beam1.NE.nuclearZ_beam2)THEN
          Aname2=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
          CALL GetNuclearInfo(Aname2,A_common(2),Z2,RA_common(2),aaval,wval)
       ELSE
          Aname2=Aname1
          A_common(2)=A_common(1)
          Z2=Z1
          RA_common(2)=RA_common(1)
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in "//TRIM(Aname1)//"+"//TRIM(Aname2)//" collisions"
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
       sigNN_inel_common=sigma_inelastic(cmenergy)
       sigNN_inel_common=sigNN_inel_common*0.1d0 ! from mb to fm^2
       ! 0.1973 is from fm to GeV-1
       bfact1=RA_common(1)/GeVm12fm*mN
       bfact2=RA_common(2)/GeVm12fm*mN
       bbx(1)=bupper
       bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       sub_num(3)=10
       init=1
    ENDIF
    x1_common=x1
    x2_common=x2
    E1_common=ebeam_PN(1)*x1
    E2_common=ebeam_PN(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       PhotonPhotonFlux_AB_hardsphere_eval=PhotonFlux_nucleus_EDFF_nob(x1_common,&
            gamma1_common,Z1,RA_common(1))
       PhotonPhotonFlux_AB_hardsphere_eval=PhotonPhotonFlux_AB_hardsphere_eval*&
            PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common(2))
       PhotonPhotonFlux_AB_hardsphere_eval=MAX(PhotonPhotonFlux_AB_hardsphere_eval,0d0)
       RETURN
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mN = Exp(bA(1))
    aax(1)=DLOG(bfact1*x1_common)
    ! b2*x2*mN = Exp(bA(2))
    aax(2)=DLOG(bfact2*x2_common)
    CALL ROMBERG_ND(PhotonPhotonFlux_AB_hardsphere_fxn,aax,bbx,3,sub_num,1,1d-5,&
         integral,ind,eval_num)
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          bbx(1)=bupper*2d0**(iter)
          bbx(2)=bupper*2d0**(iter)
          CALL ROMBERG_ND(PhotonPhotonFlux_AB_hardsphere_fxn,aax,bbx,3,sub_num,1,1d-5,&
               integral,ind,eval_num)
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       bbx(2)=bupper
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon flux at (x1,x2)=",x1,x2
          WRITE(*,*)"WARNING: use PNOHAD=1 approx. instead (most probably need to increase bupper)"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       PhotonPhotonFlux_AB_hardsphere_eval=PhotonFlux_nucleus_EDFF_nob(x1_common,&
            gamma1_common,Z1,RA_common(1))
       PhotonPhotonFlux_AB_hardsphere_eval=PhotonPhotonFlux_AB_hardsphere_eval*&
            PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common(2))
       PhotonPhotonFlux_AB_hardsphere_eval=MAX(PhotonPhotonFlux_AB_hardsphere_eval,0d0)
    ELSE
       PhotonPhotonFlux_AB_hardsphere_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_AB_hardsphere_eval

  FUNCTION PhotonPhotonFlux_AB_hardsphere_fxn(dim_num,bA)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_AB_hardsphere_fxn
    INTEGER,INTENT(IN)::dim_num ! should be 3
    ! 1/0.1973d0 from fm to GeV-1 for b
    ! x1*b1*mN=Exp(bA(1))
    ! x2*b2*mN=Exp(bA(2))
    ! bA(3) = theta_{12}
    REAL(KIND(1d0)),DIMENSION(dim_num),INTENT(IN)::bA
    REAL(KIND(1d0))::b1,b2,b12,costh,pnohad
    REAL(KIND(1d0))::Ngamma1,Ngamma2
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0)),DIMENSION(2)::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    COMMON/PhotonPhoton_AB_HS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    IF(dim_num.NE.3)THEN
       WRITE(*,*)"ERROR: PhotonPhotonFlux_AB_hardsphere_fxn is not a three dimensional function"
       STOP 1
    ENDIF
    costh=DCOS(bA(3))
    ! in unit of GeV-1
    ! x1*b1*mN=Exp(bA(1))                                                                                                         
    b1=DEXP(bA(1))/x1_common/mN
    ! x2*b2*mN=Exp(bA(2))
    b2=DEXP(bA(2))/x2_common/mN
    b12=DSQRT(b1**2+b2**2-2d0*b1*b2*costh)
    b12=b12*GeVm12fm ! from GeV-1 to fm
    pnohad=PNOHAD_AB_hardsphere(b12,0d0,A_common(1)*A_common(2),RA_common,&
         sigNN_inel_common)
    IF(pnohad.LE.0d0)THEN
       PhotonPhotonFlux_AB_hardsphere_fxn=0d0
       RETURN
    ENDIF
    Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
    Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    PhotonPhotonFlux_AB_hardsphere_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
    RETURN
  END FUNCTION PhotonPhotonFlux_AB_hardsphere_fxn

  FUNCTION PhotonPhotonFlux_AB_WoodsSaxon(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_AB_WoodsSaxon
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER::MX_save=0,MY_save=0
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD,ZD_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_AB_WoodsSaxon=0d0
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(PhotonPhotonFlux_Algorithm.EQ.2)THEN
       ! kernel-factorised engine (see the FFE section at the end of this module)
       PhotonPhotonFlux_AB_WoodsSaxon=ffe_api(3,x1,x2,1,force_pnohad1,GENERATE_PhotonPhotonFlux_GRID)
       RETURN
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin).AND..NOT.force_pnohad1&
         .AND.GENERATE_PhotonPhotonFlux_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=nseg*(-log10xmin)
       MY=MX
       ! try to deallocate first
       IF(ALLOCATED(XD_1D))THEN
          DEALLOCATE(XD_1D)
       ENDIF
       ALLOCATE(XD_1D(MX+1))
       IF(ALLOCATED(YD_1D))THEN
          DEALLOCATE(YD_1D)
       ENDIF
       ALLOCATE(YD_1D(MY+1))
       IF(ALLOCATED(ZD))THEN
          WRITE(*,*)"INFO: in PhotonPhotonFlux_AB_WoodsSaxon, the xmin of the grid has been updated to ",xmin
          DEALLOCATE(ZD)
       ENDIF
       ALLOCATE(ZD(MX+1,MY+1))
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nseg
             log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x1
             YD_1D(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP 1
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! save the values that have been calculated before
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1&
                  .AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_save(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             ZD(I,J)=PhotonPhotonFlux_AB_WoodsSaxon_eval(xx1,xx2)
             IF(MX_save.EQ.0.AND.MY_save.EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   WRITE(*,*)"INFO: generate grid of photon-photon flux in AB (will take a few minutes)"
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF(ALLOCATED(ZD_save))THEN
          DEALLOCATE(ZD_save)
       ENDIF
       ALLOCATE(ZD_save(MX+1,MY+1))
       ! save the values that have been calculated before
       MX_save=MX
       MY_save=MY
       DO I=1,MX+1
          DO J=1,MY+1
             ZD_save(I,J)=ZD(I,J)
          ENDDO
       ENDDO
       init=1
    ENDIF
    IF(.NOT.force_pnohad1)THEN
       IF(GENERATE_PhotonPhotonFlux_GRID)THEN
          XI(1)=DLOG10(x1)
          YI(1)=DLOG10(x2)
          IF(XI(1).GE.0d0)THEN
             ilog10x1=-1
          ELSE
             ilog10x1=FLOOR(XI(1))
          ENDIF
          IF(YI(1).GE.0d0)THEN
             ilog10x2=-1
          ELSE
             ilog10x2=FLOOR(YI(1))
          ENDIF
          !K=nseg*(-ilog10x1-1)
          !DO I=1,nseg+1
          !   XD2_1D(I)=XD_1D(K+I)
          !ENDDO
          !L=nseg*(-ilog10x2-1)
          !DO I=1,nseg+1
          !   YD2_1D(I)=YD_1D(L+I)
          !ENDDO
          !DO I=1,nseg+1
          !   DO J=1,nseg+1
          !      ZD2(I,J)=ZD(K+I,L+J)
          !   ENDDO
          !ENDDO
          !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
          K=nseg*(-ilog10x1-1)
          DO I=K+1,K+nseg+1
             if(XD_1D(I).LE.XI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save+1)THEN
             K=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             K=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save+1)THEN
             K=MX_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #1 !"
             STOP 1
          ENDIF
          L=nseg*(-ilog10x2-1)
          DO I=L+1,L+nseg+1
             if(YD_1D(I).LE.YI(1))EXIT
          ENDDO
          IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save+1)THEN
             L=I-n_interp/2-1
          ELSEIF(I-n_interp/2.LT.1)THEN
             L=0
          ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save+1)THEN
             L=MY_save+1-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here #2 !"
             STOP 1
          ENDIF
          DO I=1,n_interp
             XD2_1D(I)=XD_1D(K+I)
             YD2_1D(I)=YD_1D(L+I)
          ENDDO
          DO I=1,n_interp
             DO J=1,n_interp
                ZD2(I,J)=ZD(K+I,L+J)
             ENDDO
          ENDDO
          CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
               YD2_1D,ZD2,1,XI,YI,ZI)
       ELSE
          ZI(1)=PhotonPhotonFlux_AB_WoodsSaxon_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    pnohadval=PhotonPhotonFlux_AB_WoodsSaxon_eval(x1,x2,.TRUE.)
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_AB_WoodsSaxon=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.DABS(ZI(1)/pnohadval).GT.1D2)THEN
          PhotonPhotonFlux_AB_WoodsSaxon=pnohadval
       ELSE
          !IF(DABS(ZI(1)/pnohadval).GT.1D2.OR.(DABS(ZI(1)/pnohadval).LT.1D-2))THEN
          !   PRINT *, "WARNING:",x1,x2, ZI(1), pnohadval
          !   !STOP
          !ENDIF
          PhotonPhotonFlux_AB_WoodsSaxon=ZI(1)
       ENDIF
    ELSE
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_AB_WoodsSaxon=0d0
       ELSE
          PhotonPhotonFlux_AB_WoodsSaxon=pnohadval
       ENDIF
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_AB_WoodsSaxon

  FUNCTION PhotonPhotonFlux_AB_WoodsSaxon_eval(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonFlux_AB_WoodsSaxon_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0)),DIMENSION(2)::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    REAL(KIND(1d0)),DIMENSION(2)::Z_common
    REAL(KIND(1d0)),DIMENSION(2)::aaVal_common,wVal_common ! parameters in Woods-Saxon potential
    COMMON/PhotonPhoton_AB_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,Z_common,aaVal_common,wVal_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact1,bfact2
    SAVE bfact1,bfact2
    REAL(KIND(1d0)),PARAMETER::bupper=2d0
    CHARACTER(len=7)::Aname1,Aname2
    REAL(KIND(1d0))::cmenergy
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_AB_WoodsSaxon_eval=0d0
       RETURN
    ENDIF
    CALC_DELTABQ=.FALSE.
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP 1
       ENDIf
       gamma1_common=ebeam_PN(1)/mN
       gamma2_common=ebeam_PN(2)/mN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       Aname1=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       CALL GetNuclearInfo(Aname1,A_common(1),Z1,RA_common(1),aaVal_common(1),wVal_common(1))
       Z_common(1)=Z1
       IF(nuclearA_beam2.NE.nuclearA_beam1.OR.nuclearZ_beam1.NE.nuclearZ_beam2)THEN
          Aname2=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
          CALL GetNuclearInfo(Aname2,A_common(2),Z2,RA_common(2),aaVal_common(2),wVal_common(2))
          Z_common(2)=Z2
       ELSE
          Aname2=Aname1
          A_common(2)=A_common(1)
          Z_common(2)=Z_common(1)
          Z2=Z1
          RA_common(2)=RA_common(1)
          aaVal_common(2)=aaVal_common(1)
          wVal_common(2)=wVal_common(1)
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in "//TRIM(Aname1)//"+"//TRIM(Aname2)//" collisions"
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
       sigNN_inel_common=sigma_inelastic(cmenergy)
       sigNN_inel_common=sigNN_inel_common*0.1d0 ! from mb to fm^2
       ! 0.1973 is from fm to GeV-1
       IF(IEPSILON_EDFF.EQ.0)THEN
          bfact1=RA_common(1)/GeVm12fm*mN
          bfact2=RA_common(2)/GeVm12fm*mN
       ELSE
          bfact1=(RA_common(1)+aaVal_common(1))/GeVm12fm*mN
          bfact2=(RA_common(2)+aaVal_common(2))/GeVm12fm*mN
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we can integrate b down to zero
          bfact1=bfact1*LOWER_BFactor_Limit
          bfact2=bfact2*LOWER_BFactor_Limit
       ENDIF
       !bbx(1)=bupper
       !bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we should increase the number of segments
          sub_num(1)=sub_num(1)*SUB_FACTOR
          sub_num(2)=sub_num(2)*SUB_FACTOR
       ENDIF
       sub_num(3)=10
       init=1
    ENDIF
    x1_common=x1
    x2_common=x2
    E1_common=ebeam_PN(1)*x1
    E2_common=ebeam_PN(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(PhotonPhotonFlux_Algorithm.EQ.2)THEN
       ! kernel-factorised engine (see the FFE section at the end of this module)
       PhotonPhotonFlux_AB_WoodsSaxon_eval=ffe_api(3,x1,x2,1,force_pnohad1,.FALSE.)
       RETURN
    ENDIF
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonFlux_nucleus_ChFF_nob(1,x1_common,&
               gamma1_common,Z1,RA_common(1),wVal_common(1),aaVal_common(1))
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonPhotonFlux_AB_WoodsSaxon_eval*&
               PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common(2),&
               wVal_common(2),aaVal_common(2))
       ELSE
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonFlux_nucleus_EDFF_nob(x1_common,&
               gamma1_common,Z1,RA_common(1))
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonPhotonFlux_AB_WoodsSaxon_eval*&
               PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common(2))
       ENDIF
       PhotonPhotonFlux_AB_WoodsSaxon_eval=MAX(PhotonPhotonFlux_AB_WoodsSaxon_eval,0d0)
       RETURN
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mN = Exp(bA(1)) = b1*E_gamma1/gamma1 = b1tilde
    aax(1)=DLOG(bfact1*x1_common)
    ! b2*x2*mN = Exp(bA(2)) = b2*E_gamma2/gamma2 = b2tilde
    aax(2)=DLOG(bfact2*x2_common)
    ! we should also choose the upper limit dynamically
    bbx(1)=MAX(1.5d0,DLOG(0.8d0*RA_common(1)/GeVm12fm*DSQRT(x1_common*mN)))
    bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common(2)/GeVm12fm*DSQRT(x2_common*mN)))
    CALL ROMBERG_ND(PhotonPhotonFlux_AB_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
         integral,ind,eval_num)
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          bbx(1)=bbx(1)*2d0
          bbx(2)=bbx(2)*2d0
          CALL ROMBERG_ND(PhotonPhotonFlux_AB_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
               integral,ind,eval_num)
          iter=iter+1
       ENDDO
       !bbx(1)=bupper
       !bbx(2)=bupper
       bbx(1)=MAX(1.5d0,DLOG(0.8d0*RA_common(1)/GeVm12fm*DSQRT(x1_common*mN)))
       bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common(2)/GeVm12fm*DSQRT(x2_common*mN)))
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon flux at (x1,x2)=",x1,x2
          WRITE(*,*)"WARNING: use PNOHAD=1 approx. instead (most probably need to increase bupper)"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonFlux_nucleus_ChFF_nob(1,x1_common,&
               gamma1_common,Z1,RA_common(1),wVal_common(1),aaVal_common(1))
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonPhotonFlux_AB_WoodsSaxon_eval*&
               PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common(2),&
               wVal_common(2),aaVal_common(2))
       ELSE
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonFlux_nucleus_EDFF_nob(x1_common,&
               gamma1_common,Z1,RA_common(1))
          PhotonPhotonFlux_AB_WoodsSaxon_eval=PhotonPhotonFlux_AB_WoodsSaxon_eval*&
               PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common(2))
       ENDIF
       PhotonPhotonFlux_AB_WoodsSaxon_eval=MAX(PhotonPhotonFlux_AB_WoodsSaxon_eval,0d0)
    ELSE
       PhotonPhotonFlux_AB_WoodsSaxon_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonFlux_AB_WoodsSaxon_eval

  FUNCTION DeltaB_AB_at_x1x2(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::DeltaB_AB_at_x1x2
    REAL(KIND(1d0))::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1 ! If true, it only evaluates with PNOHAD=1
    LOGICAL::force_pnohad1
    REAL(KIND(1d0))::num,den
    LOGICAL::use_grid_bu
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    ! numerator
    ! PhotonPhotonFlux*Delta b in unit of fm
    num=PhotonPhotonDeltaB_AB_WoodsSaxon(x1,x2,force_pnohad1)
    ! denominator
    ! PhotonPhotonFlux
    use_grid_bu=GENERATE_PhotonPhotonFlux_GRID
    ! Let us make it to be coherent with DELTAB grid case
    GENERATE_PhotonPhotonFlux_GRID=GENERATE_DELTAB_GRID
    den=PhotonPhotonFlux_AB_WoodsSaxon(x1,x2,force_pnohad1)
    GENERATE_PhotonPhotonFlux_GRID=use_grid_bu ! recover the original one
    IF(den.LE.0d0.OR.num.EQ.0d0)THEN
       DeltaB_AB_at_x1x2=0d0
    ELSE
       DeltaB_AB_at_x1x2=num/den
    ENDIF
    RETURN
  END FUNCTION DeltaB_AB_at_x1x2

  FUNCTION PhotonPhotonDeltaB_AB_WoodsSaxon(x1,x2,FORCEPNOHAD1)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonDeltaB_AB_WoodsSaxon
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x1,log10x2
    INTEGER::ilog10x1,ilog10x2
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER,PARAMETER::n_interp=6
    INTEGER::MX,MY,I,J,K,L
    INTEGER,DIMENSION(2)::MX_save=(/0,0/),MY_save=(/0,0/)
    SAVE MX_save,MY_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D,YD_1D
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::ZD
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZD1_save,ZD2_save
    SAVE MX,MY,XD_1D,YD_1D,ZD,ZD1_save,ZD2_save
    !REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D,YD2_1D
    !REAL(KIND(1d0)),DIMENSION(nseg+1,nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D,YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,n_interp)::ZD2
    REAL(KIND(1d0))::xx1,xx2
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0))::pnohadval
    LOGICAL::force_pnohad1
    LOGICAL,DIMENSION(2)::gridready=(/.FALSE.,.FALSE./)
    SAVE gridready
    INTEGER::igrid
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonDeltaB_AB_WoodsSaxon=0d0
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(PhotonPhotonFlux_Algorithm.EQ.2)THEN
       ! kernel-factorised engine (see the FFE section at the end of this module)
       PhotonPhotonDeltaB_AB_WoodsSaxon=ffe_api(3,x1,x2,2,force_pnohad1,GENERATE_DELTAB_GRID)
       RETURN
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(force_pnohad1)THEN
       igrid=2
    ELSE
       igrid=1
    ENDIF
    IF((init.EQ.0.OR.x1.LT.xmin.OR.x2.LT.xmin.OR..NOT.gridready(igrid)).AND.GENERATE_DELTAB_GRID)THEN
       ! initialisation
       log10xmin_before=INT(DLOG10(xmin))
       IF(x1.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x1))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x1.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       IF(x2.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x2))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x2.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a 2-dim grid [xmin,1]x[xmin,1] first
       MX=nseg*(-log10xmin)
       MY=MX
       ! try to deallocate first
       IF(init.EQ.0.OR.log10xmin_before.GT.log10xmin)THEN
          IF(ALLOCATED(XD_1D))THEN
             DEALLOCATE(XD_1D)
          ENDIF
          ALLOCATE(XD_1D(MX+1))
          IF(ALLOCATED(YD_1D))THEN
             DEALLOCATE(YD_1D)
          ENDIF
          ALLOCATE(YD_1D(MY+1))
          IF(ALLOCATED(ZD))THEN
             WRITE(*,*)"INFO: in PhotonPhotonDeltaB_AB_WoodsSaxon, the xmin of the grid has been updated to ",xmin
             ! save the values that have been calculated before
             IF(gridready(1))THEN
                ! the saved grid has the size before the xmin update
                MX_save(1)=nseg*(-log10xmin_before)
                MY_save(1)=MX_save(1)
                IF(ALLOCATED(ZD1_save))THEN
                   DEALLOCATE(ZD1_save)
                ENDIF
                ALLOCATE(ZD1_save(MX_save(1)+1,MY_save(1)+1))
                DO I=1,MX_save(1)+1
                   DO J=1,MY_save(1)+1
                      ZD1_save(I,J)=ZD(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(gridready(2))THEN
                ! the saved grid has the size before the xmin update
                MX_save(2)=nseg*(-log10xmin_before)
                MY_save(2)=MX_save(2)
                IF(ALLOCATED(ZD2_save))THEN
                   DEALLOCATE(ZD2_save)
                ENDIF
                ALLOCATE(ZD2_save(MX_save(2)+1,MY_save(2)+1))
                DO I=1,MX_save(2)+1
                   DO J=1,MY_save(2)+1
                      ZD2_save(I,J)=ZD(2,I,J)
                   ENDDO
                ENDDO
             ENDIF
             DEALLOCATE(ZD)
          ENDIF
          ALLOCATE(ZD(2,MX+1,MY+1))
          K=0
          DO I=0,log10xmin+1,-1
             DO J=1,nseg
                log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
                K=K+1
                XD_1D(K)=log10x1
                YD_1D(K)=log10x1
             ENDDO
          ENDDO
          IF(K.NE.MX)THEN
             WRITE(*,*)"ERROR: K != MX"
             STOP 1
          ENDIF
          XD_1D(MX+1)=DBLE(log10xmin)
          YD_1D(MY+1)=DBLE(log10xmin)
          gridready(1:2)=.FALSE.
       ENDIF
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! save the values that have been calculated before
             IF(MX_save(igrid).GT.0.AND.MY_save(igrid).GT.0.AND.I.LE.MX_save(igrid)+1&
                  .AND.J.LE.MY_save(igrid)+1)THEN
                IF(igrid.EQ.1)THEN
                   ZD(igrid,I,J)=ZD1_save(I,J)
                ELSEIF(igrid.EQ.2)THEN
                   ZD(igrid,I,J)=ZD2_save(I,J)
                ELSE
                   WRITE(*,*)"ERROR: do not know igrid=",igrid
                ENDIF
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
             ZD(igrid,I,J)=PhotonPhotonDeltaB_AB_WoodsSaxon_eval(xx1,xx2,force_pnohad1)
             IF(MX_save(igrid).EQ.0.AND.MY_save(igrid).EQ.0)THEN
                IF(I.EQ.2.AND.J.EQ.2)THEN
                   IF(igrid.EQ.1)THEN
                      WRITE(*,*)"INFO: generate grid of photon-photon <delta b> in AB (will take a few minutes)"
                   ELSE
                      WRITE(*,*)"INFO: generate grid of photon-photon <delta b> in AB with PNOHAD=1 (will take a few minutes)"
                   ENDIF
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50,.TRUE.)
                ELSEIF(.NOT.(I.EQ.1.OR.J.EQ.1))THEN
                   ! we skip I=1 or J=1 which just returns 0
                   CALL progress(INT(((I-1)*(MY+1)+J)*50d0/((MX+1)*(MY+1))),50)
                ENDIF
             ENDIF
           ENDDO
       ENDDO
       gridready(igrid)=.TRUE.
       init=1
    ENDIF
    IF(GENERATE_DELTAB_GRID)THEN
       XI(1)=DLOG10(x1)
       YI(1)=DLOG10(x2)
       IF(XI(1).GE.0d0)THEN
          ilog10x1=-1
       ELSE
          ilog10x1=FLOOR(XI(1))
       ENDIF
       IF(YI(1).GE.0d0)THEN
          ilog10x2=-1
       ELSE
          ilog10x2=FLOOR(YI(1))
       ENDIF
       !K=nseg*(-ilog10x1-1)
       !DO I=1,nseg+1
       !   XD2_1D(I)=XD_1D(K+I)
       !ENDDO
       !L=nseg*(-ilog10x2-1)
       !DO I=1,nseg+1
       !   YD2_1D(I)=YD_1D(L+I)
       !ENDDO
       !DO I=1,nseg+1
       !   DO J=1,nseg+1
       !      ZD2(I,J)=ZD(igrid,K+I,L+J)
       !   ENDDO
       !ENDDO
       !CALL lagrange_interp_2d(nseg,nseg,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
       K=nseg*(-ilog10x1-1)
       DO I=K+1,K+nseg+1
          if(XD_1D(I).LE.XI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX+1)THEN
          K=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MX+1)THEN
          K=MX+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #1 !"
          STOP 1
       ENDIF
       L=nseg*(-ilog10x2-1)
       DO I=L+1,L+nseg+1
          if(YD_1D(I).LE.YI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY+1)THEN
          L=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          L=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MY+1)THEN
          L=MY+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #2 !"
          STOP 1
       ENDIF
       DO I=1,n_interp
          XD2_1D(I)=XD_1D(K+I)
          YD2_1D(I)=YD_1D(L+I)
       ENDDO
       DO I=1,n_interp
          DO J=1,n_interp
             ZD2(I,J)=ZD(igrid,K+I,L+J)
          ENDDO
       ENDDO
       CALL lagrange_interp_2d(n_interp-1,n_interp-1,XD2_1D,&
            YD2_1D,ZD2,1,XI,YI,ZI)
       IF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0)THEN
          PhotonPhotonDeltaB_AB_WoodsSaxon=0d0
       ELSE
          PhotonPhotonDeltaB_AB_WoodsSaxon=ZI(1)
       ENDIF
    ELSE
       PhotonPhotonDeltaB_AB_WoodsSaxon=PhotonPhotonDeltaB_AB_WoodsSaxon_eval(x1,x2,force_pnohad1)
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonDeltaB_AB_WoodsSaxon

  FUNCTION PhotonPhotonDeltaB_AB_WoodsSaxon_eval(x1,x2,FORCEPNOHAD1)
    ! this is PhotonPhotonFlux*Delta b, where Delta b is in unit of fm
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::PhotonPhotonDeltaB_AB_WoodsSaxon_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN),OPTIONAL::FORCEPNOHAD1
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0)),DIMENSION(2)::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    REAL(KIND(1d0)),DIMENSION(2)::Z_common
    REAL(KIND(1d0)),DIMENSION(2)::aaVal_common,wVal_common ! parameters in Woods-Saxon potential
    COMMON/PhotonPhoton_AB_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,Z_common,aaVal_common,wVal_common
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::integral,Z1,Z2
    SAVE Z1,Z2
    REAL(KIND(1d0)),DIMENSION(3)::aax,bbx
    INTEGER,DIMENSION(3)::sub_num
    INTEGER::ind,eval_num
    SAVE aax,bbx,sub_num
    REAL(KIND(1d0))::bfact1,bfact2
    SAVE bfact1,bfact2
    REAL(KIND(1d0)),PARAMETER::bupper=2d0
    CHARACTER(len=7)::Aname1,Aname2
    REAL(KIND(1d0))::cmenergy
    INTEGER,PARAMETER::itermax=5
    INTEGER::printnum=0,iter
    SAVE printnum
    LOGICAL::force_pnohad1
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonDeltaB_AB_WoodsSaxon_eval=0d0
       RETURN
    ENDIF
    CALC_DELTABQ=.TRUE.
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP 1
       ENDIf
       gamma1_common=ebeam_PN(1)/mN
       gamma2_common=ebeam_PN(2)/mN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       Aname1=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       CALL GetNuclearInfo(Aname1,A_common(1),Z1,RA_common(1),aaVal_common(1),wVal_common(1))
       Z_common(1)=Z1
       IF(nuclearA_beam2.NE.nuclearA_beam1.OR.nuclearZ_beam1.NE.nuclearZ_beam2)THEN
          Aname2=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
          CALL GetNuclearInfo(Aname2,A_common(2),Z2,RA_common(2),aaVal_common(2),wVal_common(2))
          Z_common(2)=Z2
       ELSE
          Aname2=Aname1
          A_common(2)=A_common(1)
          Z_common(2)=Z_common(1)
          Z2=Z1
          RA_common(2)=RA_common(1)
          aaVal_common(2)=aaVal_common(1)
          wVal_common(2)=wVal_common(1)
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in "//TRIM(Aname1)//"+"//TRIM(Aname2)//" collisions"
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
       sigNN_inel_common=sigma_inelastic(cmenergy)
       sigNN_inel_common=sigNN_inel_common*0.1d0 ! from mb to fm^2
       ! 0.1973 is from fm to GeV-1
       IF(IEPSILON_EDFF.EQ.0)THEN
          bfact1=RA_common(1)/GeVm12fm*mN
          bfact2=RA_common(2)/GeVm12fm*mN
       ELSE
          bfact1=(RA_common(1)+aaVal_common(1))/GeVm12fm*mN
          bfact2=(RA_common(2)+aaVal_common(2))/GeVm12fm*mN
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we can integrate b down to zero
          bfact1=bfact1*LOWER_BFactor_Limit
          bfact2=bfact2*LOWER_BFactor_Limit
       ENDIF
       !bbx(1)=bupper
       !bbx(2)=bupper
       aax(3)=0d0
       bbx(3)=TWOPI
       sub_num(1)=30
       sub_num(2)=30
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for the charge form factor
          ! we should increase the number of segments
          sub_num(1)=sub_num(1)*SUB_FACTOR
          sub_num(2)=sub_num(2)*SUB_FACTOR
       ENDIF
       sub_num(3)=10
       init=1
    ENDIF
    x1_common=x1
    x2_common=x2
    E1_common=ebeam_PN(1)*x1
    E2_common=ebeam_PN(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(PhotonPhotonFlux_Algorithm.EQ.2)THEN
       ! kernel-factorised engine (see the FFE section at the end of this module)
       PhotonPhotonDeltaB_AB_WoodsSaxon_eval=ffe_api(3,x1,x2,2,force_pnohad1,.FALSE.)
       RETURN
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mN = Exp(bA(1)) = b1*E_gamma1/gamma1 = b1tilde
    aax(1)=DLOG(bfact1*x1_common)
    ! b2*x2*mN = Exp(bA(2)) = b2*E_gamma2/gamma2 = b2tilde
    aax(2)=DLOG(bfact2*x2_common)
    ! we should also choose the upper limit dynamically
    bbx(1)=MAX(1.5d0,DLOG(0.8d0*RA_common(1)/GeVm12fm*DSQRT(x1_common*mN)))
    bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common(2)/GeVm12fm*DSQRT(x2_common*mN)))
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       CALL ROMBERG_ND(PhotonPhotonFlux_AB_WoodsSaxon_fxn,aax(1:2),bbx(1:2),2,sub_num(1:2),1,1d-5,&
            integral,ind,eval_num)
    ELSE
       CALL ROMBERG_ND(PhotonPhotonFlux_AB_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
            integral,ind,eval_num)
    ENDIF
    IF(ind.EQ.-1)THEN
       WRITE(*,*)"WARNING: the precision 1e-5 is not achieved"
    ENDIF
    IF(integral.LT.0d0)THEN
       ! try to rescue it by increasing bupper
       iter=1
       DO WHILE(integral.LT.0d0.AND.iter.LE.itermax)
          bbx(1)=bbx(1)*2d0
          bbx(2)=bbx(2)*2d0
          IF(force_pnohad1)THEN
             ! we only use PNOHAD=1
             CALL ROMBERG_ND(PhotonPhotonFlux_AB_WoodsSaxon_fxn,aax(1:2),bbx(1:2),2,sub_num(1:2),1,1d-5,&
                  integral,ind,eval_num)
          ELSE
             CALL ROMBERG_ND(PhotonPhotonFlux_AB_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
                  integral,ind,eval_num)
          ENDIF
          iter=iter+1
       ENDDO
       !bbx(1)=bupper
       !bbx(2)=bupper
       bbx(1)=MAX(1.5d0,DLOG(0.8d0*RA_common(1)/GeVm12fm*DSQRT(x1_common*mN)))
       bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common(2)/GeVm12fm*DSQRT(x2_common*mN)))
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon <delta b> at (x1,x2)=",x1,x2
          WRITE(*,*)"WARNING: most probably need to increase bupper"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       PhotonPhotonDeltaB_AB_WoodsSaxon_eval=0d0
    ELSE
       PhotonPhotonDeltaB_AB_WoodsSaxon_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
    ENDIF
    RETURN
  END FUNCTION PhotonPhotonDeltaB_AB_WoodsSaxon_eval

  FUNCTION PhotonPhotonFlux_AB_WoodsSaxon_fxn(dim_num,bA)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonPhotonFlux_AB_WoodsSaxon_fxn
    INTEGER,INTENT(IN)::dim_num ! should be 3 or 2
    ! 1/0.1973d0 from fm to GeV-1 for b
    ! x1*b1*mN=Exp(bA(1))
    ! x2*b2*mN=Exp(bA(2))
    ! bA(3) = theta_{12}
    REAL(KIND(1d0)),DIMENSION(dim_num),INTENT(IN)::bA
    REAL(KIND(1d0))::b1,b2,b12,costh,pnohad
    REAL(KIND(1d0))::b12GeVm1,pbreak
    REAL(KIND(1d0))::Ngamma1,Ngamma2
    REAL(KIND(1d0))::gamma1_common,gamma2_common ! Lorentz boost factors
    REAL(KIND(1d0))::E1_common,E2_common      ! energies of two photons
    REAL(KIND(1d0))::x1_common,x2_common      ! x1 and x2 of two photons
    REAL(KIND(1d0))::sigNN_inel_common        ! NN inelastic cross section
    REAL(KIND(1d0)),DIMENSION(2)::RA_common, A_common      ! radius of nuclei and atom number of nuclei
    REAL(KIND(1d0)),DIMENSION(2)::Z_common
    REAL(KIND(1d0)),DIMENSION(2)::aaVal_common,wVal_common ! parameters in Woods-Saxon potential
    COMMON/PhotonPhoton_AB_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,Z_common,aaVal_common,wVal_common
    REAL(KIND(1d0)),DIMENSION(2)::RR,aaa
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),EXTERNAL::elliptic_em
    REAL(KIND(1d0))::arg1,arg2
    IF(dim_num.NE.3.AND.dim_num.NE.2)THEN
       WRITE(*,*)"ERROR: PhotonPhotonFlux_AB_WoodsSaxon_fxn is not a three or two dimensional function"
       STOP 1
    ENDIF
    ! in unit of GeV-1
    ! x1*b1*mN=Exp(bA(1))
    b1=DEXP(bA(1))/x1_common/mN
    ! x2*b2*mN=Exp(bA(2))
    b2=DEXP(bA(2))/x2_common/mN
    IF(dim_num.EQ.3)THEN
       costh=DCOS(bA(3))
       b12=DSQRT(b1**2+b2**2-2d0*b1*b2*costh)
       b12GeVm1=b12 ! in unit of GeV-1
       b12=b12*GeVm12fm ! from GeV-1 to fm

       pnohad=PNOHAD_AB_WoodsSaxon(b12,0d0,RA_common,wVal_common,aaVal_common,&
            A_common,Z_common,sigNN_inel_common)
       IF(neutron_tagging(1).NE.-2.OR.neutron_tagging(2).NE.-2)THEN
          ! multiply the forward neutron tagging probability
          pbreak=Pbreak_AB_WoodsSaxon(neutron_tagging(1),neutron_tagging(2),&
               b12GeVm1,neutron_conj_sum)
          pnohad=pnohad*pbreak
       ENDIF
       IF(CALC_DELTABQ)THEN
          pnohad=pnohad*b12
       ENDIF
    ELSE
       ! dim_num = 2
       ! it must be a case with pnohad=1
       ! and CALC_DELTABQ=.TRUE
       IF(.NOT.CALC_DELTABQ)THEN
          WRITE(*,*)"ERROR: cannot reach here with CALC_DELTABQ=.FALSE."
       ENDIF
       ! pnohad is Integrate[Sqrt[b1^2 + b2^2 - 2 b1*b2*Cos[phi]], {phi, 0, 2 Pi}]
       arg2=4d0*b1*b2/(b1+b2)**2
       IF(b1.NE.b2.AND.arg2.LT.1d0)THEN
          arg1=-4d0*b1*b2/(b1-b2)**2
          ! elliptic_em(m) is EllipticE[m] with m<1
          pnohad=DABS(b1-b2)*elliptic_em(arg1)
          pnohad=pnohad+(b1+b2)*elliptic_em(arg2)
       ELSE
          pnohad=4d0*b1
       ENDIF
       pnohad=pnohad*2d0*GeVm12fm ! convert from GeV-1 to fm
    ENDIF
    IF(pnohad.LE.0d0)THEN
       PhotonPhotonFlux_AB_WoodsSaxon_fxn=0d0
       RETURN
    ENDIF
    IF(.NOT.USE_CHARGEFORMFACTOR4PHOTON)THEN
       Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
       Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    ELSE
       ! converting from fm to GeV-1
       RR(1)=RA_common(1)/GeVm12fm
       RR(2)=RA_common(2)/GeVm12fm
       aaa(1)=aaVal_common(1)/GeVm12fm
       aaa(2)=aaVal_common(2)/GeVm12fm
       Ngamma1=PhotonNumberDensity_ChargeFormFactor_WS(b1,E1_common,gamma1_common,&
            RR(1),wVal_common(1),aaa(1),3d0,0.7d0,1)
       Ngamma2=PhotonNumberDensity_ChargeFormFactor_WS(b2,E2_common,gamma2_common,&
            RR(2),wVal_common(2),aaa(2),3d0,0.7d0,2)
    ENDIF
    PhotonPhotonFlux_AB_WoodsSaxon_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
    RETURN
  END FUNCTION PhotonPhotonFlux_AB_WoodsSaxon_fxn

  ! photon-photon Luminosity
  ! see Eq.(7.13) in my notes OpticalGlauber.pdf
  FUNCTION Lgammagamma_UPC(scale,icoll,iprofile)
    ! icoll: 1 - pp; 2 - pA; 3 - AB
    ! iprofile: 0: P_{NOHAD}=1; 1 - Woods-Saxon; 2 - hard-sphere
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::Lgammagamma_UPC
    INTEGER,INTENT(IN)::icoll,iprofile
    REAL(KIND(1d0)),INTENT(IN)::scale ! shat=scale**2
    INTEGER::npoints
    REAL(KIND(1d0))::tau_common,s,log1oxmax
    INTEGER::collision_type_common,profile_type_common
    COMMON/Lgammagamma_UPC_com/collision_type_common,profile_type_common,&
         tau_common
    IF(photon_final_state.NE.0.AND.(icoll.EQ.3.OR.(icoll.EQ.2.AND.iprofile.EQ.2)))THEN
       WRITE(*,*)"ERROR: photon_final_state /= 0 is only available for pp and pA (Woods-Saxon)"
       STOP 1
    ENDIF
    ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
    ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
    s=4d0*ebeam_PN(1)*ebeam_PN(2)
    tau_common=scale**2/s
    collision_type_common=icoll
    profile_type_common=iprofile
    log1oxmax=DLOG(1d0/tau_common)
    IF(GENERATE_PhotonPhotonFlux_GRID)THEN
       npoints=1000
    ELSE
       npoints=100
    ENDIF
    CALL trapezoid_integration(npoints,Lgammagamma_UPC_fxn,&
         log1oxmax,Lgammagamma_UPC)
    RETURN
  END FUNCTION Lgammagamma_UPC

  FUNCTION Lgammagamma_UPC_fxn(log1ox)
    IMPLICIT NONE
    REAL(KIND(1d0))::Lgammagamma_UPC_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1ox ! log(1/x)
    REAL(KIND(1d0))::x1,x2
    REAL(KIND(1d0))::tau_common
    INTEGER::collision_type_common,profile_type_common
    COMMON/Lgammagamma_UPC_com/collision_type_common,profile_type_common,&
         tau_common
    x1=DEXP(-log1ox)
    x2=tau_common/x1
    IF(collision_type_common.EQ.1)THEN
       ! pp
       IF(profile_type_common.EQ.0)THEN
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_pp(x1,x2,.TRUE.)
       ELSE
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_pp(x1,x2)
       ENDIF
    ELSEIF(collision_type_common.EQ.2)THEN
       ! pA or Ap
       IF(profile_type_common.EQ.1)THEN
          ! Woods-Saxon
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_pA_WoodsSaxon(x1,x2)
       ELSEIF(profile_type_common.EQ.2)THEN
          ! Hard-Sphere
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_pA_hardsphere(x1,x2)
       ELSEIF(profile_type_common.EQ.0)THEN
          ! P_{NOHAD}=1
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_pA_WoodsSaxon(x1,x2,.TRUE.)
       ELSE
          WRITE(*,*)"ERROR: do not know the profile type = ",profile_type_common
          STOP 1
       ENDIF
    ELSEIF(collision_type_common.EQ.3)THEN
       ! AB
       IF(profile_type_common.EQ.1)THEN
          ! Woods-Saxon
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_AB_WoodsSaxon(x1,x2)
       ELSEIF(profile_type_common.EQ.2)THEN
          ! Hard-Sphere
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_AB_hardsphere(x1,x2)
       ELSEIF(profile_type_common.EQ.0)THEN
          ! P_{NOHAD}=1
          Lgammagamma_UPC_fxn=PhotonPhotonFlux_AB_WoodsSaxon(x1,x2,.TRUE.)
       ELSE
          WRITE(*,*)"ERROR: do not know the profile type = ",profile_type_common
          STOP 1
       ENDIF
    ELSE
       WRITE(*,*)"ERROR: do not know the collision type = ",collision_type_common
       STOP 1
    ENDIF
    RETURN
  END FUNCTION Lgammagamma_UPC_fxn
  
  ! dL/dW at W=scale
  ! dL/dW=Lgammagamma*2W/s
  ! it is used in hep-ph/0112211
  FUNCTION dLgammagammadW_UPC(scale,icoll,iprofile)
    ! icoll: 1 - pp; 2 - pA; 3 - AB
    ! iprofile: 0; P_{NOHAD}=1; 1 - Woods-Saxon; 2 - hard-sphere
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::dLgammagammadW_UPC
    INTEGER,INTENT(IN)::icoll,iprofile
    REAL(KIND(1d0)),INTENT(IN)::scale ! scale=W
    REAL(KIND(1d0))::s
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
    ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
    s=4d0*ebeam_PN(1)*ebeam_PN(2)
    dLgammagammadW_UPC=2d0*scale/s
    dLgammagammadW_UPC=dLgammagammadW_UPC*&
         Lgammagamma_UPC(scale,icoll,iprofile)
    RETURN
  END FUNCTION dLgammagammadW_UPC

  ! <delta b> at W=scale
  FUNCTION DeltaB_UPC_at_W(scale,icoll,iprofile)
    ! icoll: 1 - pp; 2 - pA; 3 - AB
    ! iprofile: 0; P_{NOHAD}=1; 1 - Woods-Saxon; 2 - hard-sphere (does not supported)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::DeltaB_UPC_at_W
    INTEGER,INTENT(IN)::icoll,iprofile
    REAL(KIND(1d0)),INTENT(IN)::scale ! scale=W
    REAL(KIND(1d0))::num,den
    LOGICAL::use_grid_bu
    IF(.NOT.print_banner)THEN
       WRITE(*,*)"==============================================================="
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|                       __    __  _______    ______           |"
       WRITE(*,*)"|                      |  \  |  \|       \  /      \          |"
       WRITE(*,*)"|   __      __         | $$  | $$| $$$$$$$\|  $$$$$$\         |"
       WRITE(*,*)"|  |  \    /  \ ______ | $$  | $$| $$__/ $$| $$   \$$         |"
       WRITE(*,*)"|   \$$ \/  $$ |      \| $$  | $$| $$    $$| $$               |"
       WRITE(*,*)"|     \$$  $$   \$$$$$$| $$  | $$| $$$$$$$ | $$   __          |"
       WRITE(*,*)"|      \$$$$           | $$__/ $$| $$      | $$__/  \         |"
       WRITE(*,*)"|      | $$             \$$    $$| $$       \$$    $$         |"
       WRITE(*,*)"|       \$$              \$$$$$$  \$$        \$$$$$$          |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    A library for exclusive photon-photon processes in       |"
       WRITE(*,*)"|    ultraperipheral proton and nuclear collisions (v2.0)     |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    By Hua-Sheng Shao (LPTHE) and David d'Enterria (CERN)    |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"|    Please cite arXiv:2207.03012 [JHEP 09 (2022) 248]        |"
       WRITE(*,*)"|                                                             |"
       WRITE(*,*)"==============================================================="
       print_banner=.TRUE.
    ENDIF
    ! numerator
    ! PhotonPhotonFlux*Delta b in unit of fm
    num=LgammagammaDeltaB_UPC(scale,icoll,iprofile)
    ! denominator
    ! PhotonPhotonFlux
    use_grid_bu=GENERATE_PhotonPhotonFlux_GRID
    ! Let us make it to be coherent with DELTAB grid case
    GENERATE_PhotonPhotonFlux_GRID=GENERATE_DELTAB_GRID
    den=Lgammagamma_UPC(scale,icoll,iprofile)
    GENERATE_PhotonPhotonFlux_GRID=use_grid_bu ! recover the original one
    IF(den.LE.0d0.OR.num.EQ.0d0)THEN
       DeltaB_UPC_at_W=0d0
    ELSE
       DeltaB_UPC_at_W=num/den
    ENDIF
    RETURN
  END FUNCTION DeltaB_UPC_at_W

  ! photon-photon Luminosity*Delta B
  FUNCTION LgammagammaDeltaB_UPC(scale,icoll,iprofile)
    ! icoll: 1 - pp; 2 - pA; 3 - AB
    ! iprofile: 0: P_{NOHAD}=1; 1 - Woods-Saxon; 2 - hard-sphere (does not supported)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::LgammagammaDeltaB_UPC
    INTEGER,INTENT(IN)::icoll,iprofile
    REAL(KIND(1d0)),INTENT(IN)::scale ! shat=scale**2
    INTEGER::npoints
    REAL(KIND(1d0))::tau_common,s,log1oxmax
    INTEGER::collision_type_common,profile_type_common
    COMMON/Lgammagamma_UPC_com/collision_type_common,profile_type_common,&
         tau_common
    ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
    ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
    s=4d0*ebeam_PN(1)*ebeam_PN(2)
    tau_common=scale**2/s
    collision_type_common=icoll
    profile_type_common=iprofile
    log1oxmax=DLOG(1d0/tau_common)
    IF(GENERATE_DELTAB_GRID)THEN
       npoints=1000
    ELSE
       npoints=100
    ENDIF
    CALL trapezoid_integration(npoints,LgammagammaDeltaB_UPC_fxn,&
         log1oxmax,LgammagammaDeltaB_UPC)
    RETURN
  END FUNCTION LgammagammaDeltaB_UPC

  FUNCTION LgammagammaDeltaB_UPC_fxn(log1ox)
    IMPLICIT NONE
    REAL(KIND(1d0))::LgammagammaDeltaB_UPC_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1ox ! log(1/x)
    REAL(KIND(1d0))::x1,x2
    REAL(KIND(1d0))::tau_common
    INTEGER::collision_type_common,profile_type_common
    COMMON/Lgammagamma_UPC_com/collision_type_common,profile_type_common,&
         tau_common
    x1=DEXP(-log1ox)
    x2=tau_common/x1
    IF(collision_type_common.EQ.1)THEN
       ! pp
       IF(profile_type_common.EQ.0)THEN
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_pp(x1,x2,.TRUE.)
       ELSE
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_pp(x1,x2)
       ENDIF
    ELSEIF(collision_type_common.EQ.2)THEN
       ! pA or Ap
       IF(profile_type_common.EQ.1)THEN
          ! Woods-Saxon
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_pA_WoodsSaxon(x1,x2)
       ELSEIF(profile_type_common.EQ.2)THEN
          ! Hard-Sphere
          WRITE(*,*)"ERROR: LgammagammaDeltaB_UPC does not support hard-sphere"
          STOP 1
       ELSEIF(profile_type_common.EQ.0)THEN
          ! P_{NOHAD}=1
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_pA_WoodsSaxon(x1,x2,.TRUE.)
       ELSE
          WRITE(*,*)"ERROR: do not know the profile type = ",profile_type_common
          STOP 1
       ENDIF
    ELSEIF(collision_type_common.EQ.3)THEN
       ! AB
       IF(profile_type_common.EQ.1)THEN
          ! Woods-Saxon
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_AB_WoodsSaxon(x1,x2)
       ELSEIF(profile_type_common.EQ.2)THEN
          ! Hard-Sphere
          WRITE(*,*)"ERROR: LgammagammaDeltaB_UPC does not support hard-sphere"
          STOP 1
       ELSEIF(profile_type_common.EQ.0)THEN
          ! P_{NOHAD}=1
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_AB_WoodsSaxon(x1,x2,.TRUE.)
       ELSE
          WRITE(*,*)"ERROR: do not know the profile type = ",profile_type_common
          STOP 1
       ENDIF
    ELSE
       WRITE(*,*)"ERROR: do not know the collision type = ",collision_type_common
       STOP 1
    ENDIF
    RETURN
  END FUNCTION LgammagammaDeltaB_UPC_fxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                        !
!   Tagging forward neutrons                                             !
!                                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Pbreak_AB_WoodsSaxon(in1,in2,Db,CONJUGATESUM)
    ! It gives us the probability Pbreak for Coulomb breakup of nuclei in nuclear-nulcear collisions
    ! Db is in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    ! in1, in2 = -2(0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)
    ! CONJUGATESUM is an optional logical variable (default is .FALSE.)
    ! If it is true, for in1=!=in2, it will sum (in1,in2)+(in2,in1)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::Pbreak_AB_WoodsSaxon
    INTEGER,INTENT(IN)::in1,in2
    REAL(KIND(1d0)),INTENT(IN)::Db
    LOGICAL,INTENT(IN),OPTIONAL::CONJUGATESUM ! If it is true, for in1=!=in2, it will sum (in1,in2)+(in2,in1) [default .FALSE.]
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),SAVE::gamma1_save,gamma2_save ! Lorentz boost factors in the lab frame
    REAL(KIND(1d0)),DIMENSION(2),SAVE::RA_save ! radius of nuclei
    INTEGER,DIMENSION(2),SAVE::A_save,Z_save ! atomic mass and atomic number of nuclei
    REAL(KIND(1d0)),DIMENSION(2),SAVE::aaVal_save,wVal_save ! parameters in Woods-Saxon potential
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0))::Z1,Z2,A1,A2
    INTEGER::i
    CHARACTER(len=7)::Aname1,Aname2
    LOGICAL::conj_sum
    REAL(KIND(1d0))::Pbreak_conj
    REAL(KIND(1d0))::PXn_LO_proj1,P1n_LO_proj1,P2n_LO_proj1,P3n_LO_proj1,P4n_LO_proj1
    REAL(KIND(1d0))::PXn_LO_proj2,P1n_LO_proj2,P2n_LO_proj2,P3n_LO_proj2,P4n_LO_proj2
    REAL(KIND(1d0))::PXn_proj1,P0n_proj1,P1n_proj1,P2n_proj1,P3n_proj1,P4n_proj1
    REAL(KIND(1d0))::PXn_proj2,P0n_proj2,P1n_proj2,P2n_proj2,P3n_proj2,P4n_proj2
    
    IF(in1.LT.-2.OR.in1.GT.4)THEN
       WRITE(*,*)"ERROR: only the following options for in1: -2 (0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)"
       STOP 1
    ENDIF
    IF(in2.LT.-2.OR.in2.GT.4)THEN
       WRITE(*,*)"ERROR: only the following options for in2: -2 (0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)"
       STOP 1
    ENDIF

    IF(in1.EQ.-2.AND.in2.EQ.-2)THEN
       Pbreak_AB_WoodsSaxon=1d0
       RETURN
    ENDIF
    
    conj_sum=.FALSE.
    IF(PRESENT(CONJUGATESUM))THEN
       conj_sum=CONJUGATESUM
    ENDIF
    IF(in1.EQ.in2)conj_sum=.FALSE.
    
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP 1
       ENDIf
       gamma1_save=ebeam_PN(1)/mN
       gamma2_save=ebeam_PN(2)/mN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       Aname1=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       CALL GetNuclearInfo(Aname1,A1,Z1,RA_save(1),aaVal_save(1),wVal_save(1))
       Z_save(1)=INT(Z1)
       A_save(1)=INT(A1)
       IF(nuclearA_beam2.NE.nuclearA_beam1.OR.nuclearZ_beam1.NE.nuclearZ_beam2)THEN
          Aname2=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
          CALL GetNuclearInfo(Aname2,A2,Z2,RA_save(2),aaVal_save(2),wVal_save(2))
          Z_save(2)=INT(Z2)
          A_save(2)=INT(A2)
       ELSE
          Aname2=Aname1
          A_save(2)=A_save(1)
          Z_save(2)=Z_save(1)
          Z2=Z1
          RA_save(2)=RA_save(1)
          aaVal_save(2)=aaVal_save(1)
          wVal_save(2)=wVal_save(1)
       ENDIF
       ! 0.1973 is from fm to GeV-1
       DO i=1,2
          RA_save(i)=RA_save(i)/GeVm12fm
          aaVal_save(i)=aaVal_save(i)/GeVm12fm
       ENDDO
       
       init=1
    ENDIF

    Pbreak_AB_WoodsSaxon=1d0
    Pbreak_conj=1d0
    ! For in2, projectile is beam 1 and target is beam 2
    ! P^{Xn,(1)}
    IF(in2.NE.-2.OR.conj_sum)THEN
       PXn_LO_proj1=PbreakXn_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
            RA_save(1),wVal_save(1),aaVal_save(1),Emin_Xn,Emax_Xn,1)
       PXn_LO_proj1=PXn_LO_proj1*Z_save(1)**2*alpha
       P0n_proj1=DEXP(-PXn_LO_proj1)
       IF(in2.EQ.-1.OR.(conj_sum.AND.in1.EQ.-1))THEN
          ! Xn
          PXn_proj1=MAX(1d0-P0n_proj1,0d0)
          IF(in2.EQ.-1)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*PXn_proj1
          ENDIF
          IF(conj_sum.AND.in1.EQ.-1)THEN
             Pbreak_conj=Pbreak_conj*PXn_proj1
          ENDIF
       ENDIF
       IF(in2.EQ.0.OR.(conj_sum.AND.in1.EQ.0))THEN
          ! 0n
          IF(in2.EQ.0)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P0n_proj1
          ENDIF
          IF(conj_sum.AND.in1.EQ.0)THEN
             Pbreak_conj=Pbreak_conj*P0n_proj1
          ENDIF
       ENDIF
       IF(in2.EQ.1.OR.(conj_sum.AND.in1.EQ.1))THEN
          ! 1n
          P1n_LO_proj1=Pbreak1n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_1n,Emax_1n,1)
          P1n_LO_proj1=P1n_LO_proj1*Z_save(1)**2*alpha
          P1n_proj1=P1n_LO_proj1*P0n_proj1
          IF(in2.EQ.1)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P1n_proj1
          ENDIF
          IF(conj_sum.AND.in1.EQ.1)THEN
             Pbreak_conj=Pbreak_conj*P1n_proj1
          ENDIF
       ENDIF
       IF(in2.EQ.2.OR.(conj_sum.AND.in1.EQ.2))THEN
          ! 2n
          P1n_LO_proj1=Pbreak1n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_1n,Emax_1n,1)
          P1n_LO_proj1=P1n_LO_proj1*Z_save(1)**2*alpha
          P2n_LO_proj1=Pbreak2n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_2n,Emax_2n,1)
          P2n_LO_proj1=P2n_LO_proj1*Z_save(1)**2*alpha
          P2n_proj1=(P1n_LO_proj1**2/2d0+P2n_LO_proj1)*P0n_proj1
          IF(in2.EQ.2)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P2n_proj1
          ENDIF
          IF(conj_sum.AND.in1.EQ.2)THEN
             Pbreak_conj=Pbreak_conj*P2n_proj1
          ENDIF
       ENDIF
       IF(in2.EQ.3.OR.(conj_sum.AND.in1.EQ.3))THEN
          ! 3n
          P1n_LO_proj1=Pbreak1n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_1n,Emax_1n,1)
          P1n_LO_proj1=P1n_LO_proj1*Z_save(1)**2*alpha
          P2n_LO_proj1=Pbreak2n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_2n,Emax_2n,1)
          P2n_LO_proj1=P2n_LO_proj1*Z_save(1)**2*alpha
          P3n_LO_proj1=Pbreak3n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_3n,Emax_3n,1)
          P3n_LO_proj1=P3n_LO_proj1*Z_save(1)**2*alpha
          P3n_proj1=(P1n_LO_proj1**3/6d0+P1n_LO_proj1*P2n_LO_proj1+P3n_LO_proj1)*P0n_proj1
          IF(in2.EQ.3)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P3n_proj1
          ENDIF
          IF(conj_sum.AND.in1.EQ.3)THEN
             Pbreak_conj=Pbreak_conj*P3n_proj1
          ENDIF
       ENDIF
       IF(in2.EQ.4.OR.(conj_sum.AND.in1.EQ.4))THEN
          ! 4n
          P1n_LO_proj1=Pbreak1n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_1n,Emax_1n,1)
          P1n_LO_proj1=P1n_LO_proj1*Z_save(1)**2*alpha
          P2n_LO_proj1=Pbreak2n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_2n,Emax_2n,1)
          P2n_LO_proj1=P2n_LO_proj1*Z_save(1)**2*alpha
          P3n_LO_proj1=Pbreak3n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_3n,Emax_3n,1)
          P3n_LO_proj1=P3n_LO_proj1*Z_save(1)**2*alpha
          P4n_LO_proj1=Pbreak4n_LO_AB_WS(Db,A_save(2),Z_save(2),gamma2_save,gamma1_save,&
               RA_save(1),wVal_save(1),aaVal_save(1),Emin_4n,Emax_4n,1)
          P4n_LO_proj1=P4n_LO_proj1*Z_save(1)**2*alpha
          P4n_proj1=(P1n_LO_proj1**4/24d0+P1n_LO_proj1**2*P2n_LO_proj1/2d0&
               +P1n_LO_proj1*P3n_LO_proj1+P2n_LO_proj1**2/2d0+P4n_LO_proj1)*P0n_proj1
          IF(in2.EQ.4)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P4n_proj1
          ENDIF
          IF(conj_sum.AND.in1.EQ.4)THEN
             Pbreak_conj=Pbreak_conj*P4n_proj1
          ENDIF
       ENDIF
    ENDIF
    ! For in1, projectile is beam 2 and target is beam 1
    IF(in1.NE.-2.OR.conj_sum)THEN
       PXn_LO_proj2=PbreakXn_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
            RA_save(2),wVal_save(2),aaVal_save(2),Emin_Xn,Emax_Xn,2)
       PXn_LO_proj2=PXn_LO_proj2*Z_save(2)**2*alpha
       P0n_proj2=DEXP(-PXn_LO_proj2)
       IF(in1.EQ.-1.OR.(conj_sum.AND.in2.EQ.-1))THEN
          ! Xn
          PXn_proj2=MAX(1d0-P0n_proj2,0d0)
          IF(in1.EQ.-1)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*PXn_proj2
          ENDIF
          IF(conj_sum.AND.in2.EQ.-1)THEN
             Pbreak_conj=Pbreak_conj*PXn_proj2
          ENDIF
       ENDIF
       IF(in1.EQ.0.OR.(conj_sum.AND.in2.EQ.0))THEN
          ! 0n
          IF(in1.EQ.0)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P0n_proj2
          ENDIF
          IF(conj_sum.AND.in2.EQ.0)THEN
             Pbreak_conj=Pbreak_conj*P0n_proj2
          ENDIF
       ENDIF
       IF(in1.EQ.1.OR.(conj_sum.AND.in2.EQ.1))THEN
          ! 1n
          P1n_LO_proj2=Pbreak1n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_1n,Emax_1n,2)
          P1n_LO_proj2=P1n_LO_proj2*Z_save(2)**2*alpha
          P1n_proj2=P1n_LO_proj2*P0n_proj2
          IF(in1.EQ.1)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P1n_proj2
          ENDIF
          IF(conj_sum.AND.in2.EQ.1)THEN
             Pbreak_conj=Pbreak_conj*P1n_proj2
          ENDIF
       ENDIF
       IF(in1.EQ.2.OR.(conj_sum.AND.in2.EQ.2))THEN
          ! 2n
          P1n_LO_proj2=Pbreak1n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_1n,Emax_1n,2)
          P1n_LO_proj2=P1n_LO_proj2*Z_save(2)**2*alpha
          P2n_LO_proj2=Pbreak2n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_2n,Emax_2n,2)
          P2n_LO_proj2=P2n_LO_proj2*Z_save(2)**2*alpha
          P2n_proj2=(P1n_LO_proj2**2/2d0+P2n_LO_proj2)*P0n_proj2
          IF(in1.EQ.2)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P2n_proj2
          ENDIF
          IF(conj_sum.AND.in2.EQ.2)THEN
             Pbreak_conj=Pbreak_conj*P2n_proj2
          ENDIF
       ENDIF
       IF(in1.EQ.3.OR.(conj_sum.AND.in2.EQ.3))THEN
          ! 3n
          P1n_LO_proj2=Pbreak1n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_1n,Emax_1n,2)
          P1n_LO_proj2=P1n_LO_proj2*Z_save(2)**2*alpha
          P2n_LO_proj2=Pbreak2n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_2n,Emax_2n,2)
          P2n_LO_proj2=P2n_LO_proj2*Z_save(2)**2*alpha
          P3n_LO_proj2=Pbreak3n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_3n,Emax_3n,2)
          P3n_LO_proj2=P3n_LO_proj2*Z_save(2)**2*alpha
          P3n_proj2=(P1n_LO_proj2**3/6d0+P1n_LO_proj2*P2n_LO_proj2+P3n_LO_proj2)*P0n_proj2
          IF(in1.EQ.3)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P3n_proj2
          ENDIF
          IF(conj_sum.AND.in2.EQ.3)THEN
             Pbreak_conj=Pbreak_conj*P3n_proj2
          ENDIF
       ENDIF
       IF(in1.EQ.4.OR.(conj_sum.AND.in2.EQ.4))THEN
          ! 4n
          P1n_LO_proj2=Pbreak1n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_1n,Emax_1n,2)
          P1n_LO_proj2=P1n_LO_proj2*Z_save(2)**2*alpha
          P2n_LO_proj2=Pbreak2n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_2n,Emax_2n,2)
          P2n_LO_proj2=P2n_LO_proj2*Z_save(2)**2*alpha
          P3n_LO_proj2=Pbreak3n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_3n,Emax_3n,2)
          P3n_LO_proj2=P3n_LO_proj2*Z_save(2)**2*alpha
          P4n_LO_proj2=Pbreak4n_LO_AB_WS(Db,A_save(1),Z_save(1),gamma1_save,gamma2_save,&
               RA_save(2),wVal_save(2),aaVal_save(2),Emin_4n,Emax_4n,2)
          P4n_LO_proj2=P4n_LO_proj2*Z_save(2)**2*alpha
          P4n_proj2=(P1n_LO_proj2**4/24d0+P1n_LO_proj2**2*P2n_LO_proj2/2d0&
               +P1n_LO_proj2*P3n_LO_proj2+P2n_LO_proj2**2/2d0+P4n_LO_proj2)*P0n_proj2
          IF(in1.EQ.4)THEN
             Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon*P4n_proj2
          ENDIF
          IF(conj_sum.AND.in2.EQ.4)THEN
             Pbreak_conj=Pbreak_conj*P4n_proj2
          ENDIF
       ENDIF
    ENDIF
    ! to sum the conjugated one
    IF(conj_sum)THEN
       Pbreak_AB_WoodsSaxon=Pbreak_AB_WoodsSaxon+Pbreak_conj
    ENDIF
    
    RETURN
  END FUNCTION Pbreak_AB_WoodsSaxon

  FUNCTION Pbreak_pA_WoodsSaxon(in2,Db)
    ! It gives us the probability Pbreak for Coulomb breakup of nuclei in proton-ion collisions
    ! Db is in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    ! in2 = -2(0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::Pbreak_pA_WoodsSaxon
    INTEGER,INTENT(IN)::in2
    REAL(KIND(1d0)),INTENT(IN)::Db
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),SAVE::gamma1_save,gamma2_save ! Lorentz boost factors
    REAL(KIND(1d0)),SAVE::RA_save ! radius of nuclei
    INTEGER,SAVE::A_save,Z_save ! atomic mass and atomic number of nuclei
    REAL(KIND(1d0)),SAVE::aaVal_save,wVal_save ! parameters in Woods-Saxon potential
    REAL(KIND(1d0))::alpha
    INTEGER::init=0
    SAVE init,alpha
    REAL(KIND(1d0))::Z1,Z2,A1,A2
    CHARACTER(len=7)::Aname
    REAL(KIND(1d0))::PXn_LO_proj1,P1n_LO_proj1,P2n_LO_proj1,P3n_LO_proj1,P4n_LO_proj1
    REAL(KIND(1d0))::PXn_proj1,P0n_proj1,P1n_proj1,P2n_proj1,P3n_proj1,P4n_proj1

    IF(in2.LT.-2.OR.in2.GT.4)THEN
       WRITE(*,*)"ERROR: only the following options for in2: -2 (0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)"
       STOP 1
    ENDIF

    IF(in2.EQ.-2)THEN
       Pbreak_pA_WoodsSaxon=1d0
       RETURN
    ENDIF
    
    IF(init.EQ.0)THEN
       IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
          nuclearA_beam1=0
          nuclearZ_beam1=0
       ELSE
          nuclearA_beam1=nb_proton(1)+nb_neutron(1)
          nuclearZ_beam1=nb_proton(1)
       ENDIF
       IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
          nuclearA_beam2=0
          nuclearZ_beam2=0
       ELSE
          nuclearA_beam2=nb_proton(2)+nb_neutron(2)
          nuclearZ_beam2=nb_proton(2)
       ENDIF
       ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
       ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_save=ebeam_PN(2)/mproton
          gamma2_save=ebeam_PN(1)/mN
       ELSE
          gamma1_save=ebeam_PN(1)/mproton
          gamma2_save=ebeam_PN(2)/mN
       ENDIF
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! proton Z is 1
       Z1=1d0
       ! read the nuclei information
       !nuclear_dir="./nuclear/"
       IF(nuclearA_beam1.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       ELSEIF(nuclearA_beam2.NE.0)THEN
          Aname=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
       ELSE
          WRITE(*,*)"ERROR: please set nuclearA_beam1/nuclearZ_beam1 or nuclearA_beam2/nuclearZ_beam2 nonzero first !"
          STOP 1
       ENDIF
       CALL GetNuclearInfo(Aname,A2,Z2,RA_save,aaVal_save,wVal_save)
       A_save=INT(A2)
       Z_save=INT(Z2)
       ! 0.1973 is from fm to GeV-1
       RA_save=RA_save/GeVm12fm
       aaVal_save=aaVal_save/GeVm12fm
       init=1
    ENDIF

    Pbreak_pA_WoodsSaxon=1d0
    PXn_LO_proj1=PbreakXn_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_Xn,Emax_Xn)
    PXn_LO_proj1=PXn_LO_proj1*alpha
    P0n_proj1=DEXP(-PXn_LO_proj1)
    IF(in2.EQ.-1)THEN
       ! Xn
       PXn_proj1=MAX(1d0-P0n_proj1,0d0)
       Pbreak_pA_WoodsSaxon=Pbreak_pA_WoodsSaxon*PXn_proj1
    ELSEIF(in2.EQ.0)THEN
       ! 0n
       Pbreak_pA_WoodsSaxon=Pbreak_pA_WoodsSaxon*P0n_proj1
    ELSEIF(in2.EQ.1)THEN
       ! 1n
       P1n_LO_proj1=Pbreak1n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_1n,Emax_1n)
       P1n_LO_proj1=P1n_LO_proj1*alpha
       P1n_proj1=P1n_LO_proj1*P0n_proj1
       Pbreak_pA_WoodsSaxon=Pbreak_pA_WoodsSaxon*P1n_proj1
    ELSEIF(in2.EQ.2)THEN
       ! 2n
       P1n_LO_proj1=Pbreak1n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_1n,Emax_1n)
       P1n_LO_proj1=P1n_LO_proj1*alpha
       P2n_LO_proj1=Pbreak2n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_2n,Emax_2n)
       P2n_LO_proj1=P2n_LO_proj1*alpha
       P2n_proj1=(P1n_LO_proj1**2/2d0+P2n_LO_proj1)*P0n_proj1
       Pbreak_pA_WoodsSaxon=Pbreak_pA_WoodsSaxon*P2n_proj1
    ELSEIF(in2.EQ.3)THEN
       ! 3n
       P1n_LO_proj1=Pbreak1n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_1n,Emax_1n)
       P1n_LO_proj1=P1n_LO_proj1*alpha
       P2n_LO_proj1=Pbreak2n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_2n,Emax_2n)
       P2n_LO_proj1=P2n_LO_proj1*alpha
       P3n_LO_proj1=Pbreak3n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_3n,Emax_3n)
       P3n_LO_proj1=P3n_LO_proj1*alpha
       P3n_proj1=(P1n_LO_proj1**3/6d0+P1n_LO_proj1*P2n_LO_proj1+P3n_LO_proj1)*P0n_proj1
       Pbreak_pA_WoodsSaxon=Pbreak_pA_WoodsSaxon*P3n_proj1
    ELSEIF(in2.EQ.4)THEN
       ! 4n
       P1n_LO_proj1=Pbreak1n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_1n,Emax_1n)
       P1n_LO_proj1=P1n_LO_proj1*alpha
       P2n_LO_proj1=Pbreak2n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_2n,Emax_2n)
       P2n_LO_proj1=P2n_LO_proj1*alpha
       P3n_LO_proj1=Pbreak3n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_3n,Emax_3n)
       P3n_LO_proj1=P3n_LO_proj1*alpha
       P4n_LO_proj1=Pbreak4n_LO_pA(Db,A_save,Z_save,gamma2_save,gamma1_save,Emin_4n,Emax_4n)
       P4n_LO_proj1=P4n_LO_proj1*alpha
       P4n_proj1=(P1n_LO_proj1**4/24d0+P1n_LO_proj1**2*P2n_LO_proj1/2d0&
            +P1n_LO_proj1*P3n_LO_proj1+P2n_LO_proj1**2/2d0+P4n_LO_proj1)*P0n_proj1
       Pbreak_pA_WoodsSaxon=Pbreak_pA_WoodsSaxon*P4n_proj1
    ENDIF

    RETURN
  END FUNCTION Pbreak_pA_WoodsSaxon

  FUNCTION PbreakXn_LO_AB_WS(Db,A_target,Z_target,gamma_target,&
       gamma_proj,RR_proj,w_proj,aa_proj,Egamma_min_target,Egamma_max_target,&
       ibeam_proj,integ_method)
    ! It gives us the lowest-order probability P^{Xn,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b,RR,aa should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ibeam_proj
    REAL(KIND(1d0))::PbreakXn_LO_AB_WS
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::RR_proj,w_proj,aa_proj ! nucleu parameters for projectile
    
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 1 MeV (0.001 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/PbreakXn_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/PbreakXn_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=4
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),DIMENSION(2)::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0)),DIMENSION(2)::R_proj_save,w_proj_save,aa_proj_save
    REAL(KIND(1d0)),DIMENSION(2)::gamma_proj_save,gamma_target_save
    INTEGER,DIMENSION(2)::A_target_save,Z_target_save
    SAVE R_proj_save,w_proj_save,aa_proj_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save
    
    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D
    
    IF(ibeam_proj.GT.2.OR.ibeam_proj.LT.1)THEN
       WRITE(*,*)"Error: ibeam_proj=/=1,2 in PbreakXn_LO_AB_WS"
       STOP 1
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mN)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)
    
    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    ibeam_proj_common=ibeam_proj
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj
    R_proj_common=RR_proj
    w_proj_common=w_proj
    aa_proj_common=aa_proj
  
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init(ibeam_proj).EQ.0.OR.(integ_method6.LT.0))THEN
       R_proj_save(ibeam_proj)=RR_proj
       w_proj_save(ibeam_proj)=w_proj
       aa_proj_save(ibeam_proj)=aa_proj
       gamma_proj_save(ibeam_proj)=gamma_proj
       gamma_target_save(ibeam_proj)=gamma_target
       A_target_save(ibeam_proj)=A_target
       Z_target_save(ibeam_proj)=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             PbreakXn_LO_AB_WS=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,PbreakXn_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(PbreakXn_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints) 
             ENDIF
             PbreakXn_LO_AB_WS=integral
             RETURN
          ENDIF
       ENDIF
       
       rescaling_bmax_save(ibeam_proj)=MAX(db,bmaxoR*RR_proj)
       
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(2,NXA))
       ENDIF
       
       dbb=9d0/DBLE(NBSEG)
       IF(init(3-ibeam_proj).EQ.1.AND.&
            rescaling_bmax_save(ibeam_proj).EQ.rescaling_bmax_save(3-ibeam_proj).AND.&
            R_proj_save(ibeam_proj).EQ.R_proj_save(3-ibeam_proj).AND.&
            w_proj_save(ibeam_proj).EQ.w_proj_save(3-ibeam_proj).AND.&
            aa_proj_save(ibeam_proj).EQ.aa_proj_save(3-ibeam_proj).AND.&
            gamma_proj_save(ibeam_proj).EQ.gamma_proj_save(3-ibeam_proj).AND.&
            gamma_target_save(ibeam_proj).EQ.gamma_target_save(3-ibeam_proj).AND.&
            A_target_save(ibeam_proj).EQ.A_target_save(3-ibeam_proj).AND.&
            Z_target_save(ibeam_proj).EQ.Z_target_save(3-ibeam_proj))THEN
          DO k=1,NXA
             XA(ibeam_proj,k)=XA(3-ibeam_proj,k)
             YA(ibeam_proj,k)=YA(3-ibeam_proj,k)
          ENDDO
       ELSE
          WRITE(*,*)"INFO: generate grid of LO Pbreak(Xn) from charge form factor of beam_proj=",ibeam_proj
          WRITE(*,*)"INFO: it will take a few seconds !"
          k=0
          DO n=0,nbmax
             ! from 10**(-n-1)*bmax to 10**(-n)*bmax
             DO i=1,NBSEG
                k=NBSEG*n+i
                ! these are b in unit GeV-1
                XA(ibeam_proj,NXA-k+1)=(10d0**(-n-1))*(one&
                     +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save(ibeam_proj)
             ENDDO
          ENDDO
          IF(k+1.NE.NXA)THEN
             WRITE(*,*)"ERROR: mismatching k+1 and NXA in PbreakXn_LO_AB_WS"
             STOP 1
          ENDIF
          XA(ibeam_proj,1)=0d0
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam_proj,I).EQ.0d0)THEN
                YA(ibeam_proj,I)=0d0
                CYCLE
             ENDIF

             npoints=30000
             b_proj_common=XA(ibeam_proj,I)
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,PbreakXn_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(PbreakXn_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             YA(ibeam_proj,I)=integral
          ENDDO
       ENDIF
       init(ibeam_proj)=1
    ENDIF
  
    IF(R_proj_save(ibeam_proj).NE.RR_proj&
         .OR.w_proj_save(ibeam_proj).NE.w_proj&
         .OR.aa_proj_save(ibeam_proj).NE.aa_proj)THEN
       WRITE(*,*)"ERROR: RA,wA,aA are not consistent in PbreakXn_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (RA,wA,aA)=",R_proj_save(ibeam_proj),&
            w_proj_save(ibeam_proj),aa_proj_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (RA,wA,aA)=",RR_proj,w_proj,aa_proj
       STOP 1
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in PbreakXn_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in PbreakXn_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF
  
    IF(db.GT.rescaling_bmax_save(ibeam_proj).OR.db.LE.0d0)THEN
       PbreakXn_LO_AB_WS=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save(ibeam_proj),1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(ibeam_proj,k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(ibeam_proj,i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(ibeam_proj,K+I)
          YD2_1D(I)=YA(ibeam_proj,K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       PbreakXn_LO_AB_WS=integral
       
    ENDIF
    RETURN
  END FUNCTION PbreakXn_LO_AB_WS

  FUNCTION PbreakXn_LO_AB_WS_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::PbreakXn_LO_AB_WS_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/PbreakXn_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/PbreakXn_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,SAVE::init=0
    
    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif
    
    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       PbreakXn_LO_AB_WS_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mN*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_WS(b_proj_common,Ega_lab,&
            gamma_proj_common,R_proj_common,w_proj_common,aa_proj_common,&
            3d0,0.7d0,ibeam_proj_common)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    PbreakXn_LO_AB_WS_fxn=GAMMAABS(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    PbreakXn_LO_AB_WS_fxn=PbreakXn_LO_AB_WS_fxn*CFF
    
    RETURN
  END FUNCTION PbreakXn_LO_AB_WS_fxn

  FUNCTION PbreakXn_LO_pA(Db,A_target,Z_target,gamma_target,&
       gamma_proj,Egamma_min_target,Egamma_max_target,&
       integ_method)
    ! It gives us the lowest-order probability P^{Xn,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::PbreakXn_LO_pA
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 1 MeV (0.001 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    
    INTEGER::A_target_common,Z_target_common
    COMMON/PbreakXn_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/PbreakXn_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common

    INTEGER::init=0
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::RA=35d0 ! in unit of GeV-1
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=5
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0))::gamma_proj_save,gamma_target_save
    INTEGER::A_target_save,Z_target_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save

    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D

    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mproton)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)

    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj

    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init.EQ.0.OR.integ_method6.LT.0)THEN
       gamma_proj_save=gamma_proj
       gamma_target_save=gamma_target
       A_target_save=A_target
       Z_target_save=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             PbreakXn_LO_pA=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
		Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,PbreakXn_LO_pA_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(PbreakXn_LO_pA_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             PbreakXn_LO_pA=integral
             RETURN
          ENDIF
       ENDIF

       rescaling_bmax_save=MAX(db,bmaxoR*RA)

       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(NXA))
       ENDIF

       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(NXA))
       ENDIF

       dbb=9d0/DBLE(NBSEG)
       WRITE(*,*)"INFO: generate grid of LO Pbreak(Xn) from charge form factor of pA"
       WRITE(*,*)"INFO: it will take a few seconds !"
       k=0
       DO n=0,nbmax
          ! from 10**(-n-1)*bmax to 10**(-n)*bmax
          DO i=1,NBSEG
             k=NBSEG*n+i
             ! these are b in unit GeV-1
             XA(NXA-k+1)=(10d0**(-n-1))*(one&
                  +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save
          ENDDO
       ENDDO
       IF(k+1.NE.NXA)THEN
          WRITE(*,*)"ERROR: mismatching k+1 and NXA in PbreakXn_LO_pA"
          STOP 1
       ENDIF
       XA(1)=0d0
       DO I=1,NXA
          IF(I.EQ.1)THEN
             CALL progress(INT(I*50d0/NXA),50,.TRUE.)
          ELSE
             CALL progress(INT(I*50d0/NXA),50)
          ENDIF
          IF(XA(I).EQ.0d0)THEN
             YA(I)=0d0
             CYCLE
          ENDIF
          npoints=30000
          b_proj_common=XA(I)
          IF(Egamma_max_target.LT.0d0)THEN
             ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0       
             Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
             Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
             xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
             log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
          ENDIF
          IF(log1oxgamma_min.EQ.0d0)THEN
             CALL trapezoid_integration(npoints,PbreakXn_LO_pA_fxn,&
                  log1oxgamma_max,integral)
          ELSE
             CALL simpson(PbreakXn_LO_pA_fxn,log1oxgamma_min,&
                  log1oxgamma_max,integral,npoints)
          ENDIF
          YA(I)=integral
       ENDDO
       init=1
    ENDIF

    IF(gamma_proj_save.NE.gamma_proj&
         .OR.gamma_target_save.NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in PbreakXn_LO_pA"
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save,&
            gamma_target_save
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in PbreakXn_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF

    IF(db.GT.rescaling_bmax_save.OR.db.LE.0d0)THEN
       PbreakXn_LO_pA=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save,1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(K+I)
          YD2_1D(I)=YA(K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       PbreakXn_LO_pA=integral

    ENDIF
    RETURN
  END FUNCTION PbreakXn_LO_pA

  FUNCTION PbreakXn_LO_pA_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::PbreakXn_LO_pA_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::A_target_common,Z_target_common
    COMMON/PbreakXn_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/PbreakXn_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif

    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       PbreakXn_LO_pA_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mproton*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_proton(b_proj_common,Ega_lab,&
            gamma_proj_common,proton_FF_correction,.FALSE.)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    PbreakXn_LO_pA_fxn=GAMMAABS(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    PbreakXn_LO_pA_fxn=PbreakXn_LO_pA_fxn*CFF

    RETURN
  END FUNCTION PbreakXn_LO_pA_fxn

  FUNCTION Pbreak1n_LO_AB_WS(Db,A_target,Z_target,gamma_target,&
       gamma_proj,RR_proj,w_proj,aa_proj,Egamma_min_target,Egamma_max_target,&
       ibeam_proj,integ_method)
    ! It gives us the lowest-order probability P^{1n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b,RR,aa should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ibeam_proj
    REAL(KIND(1d0))::Pbreak1n_LO_AB_WS
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::RR_proj,w_proj,aa_proj ! nucleu parameters for projectile
    
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 1 MeV (0.001 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 40 MeV (0.04 GeV).
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak1n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak1n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=4
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),DIMENSION(2)::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0)),DIMENSION(2)::R_proj_save,w_proj_save,aa_proj_save
    REAL(KIND(1d0)),DIMENSION(2)::gamma_proj_save,gamma_target_save
    INTEGER,DIMENSION(2)::A_target_save,Z_target_save
    SAVE R_proj_save,w_proj_save,aa_proj_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save
    
    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D
    
    IF(ibeam_proj.GT.2.OR.ibeam_proj.LT.1)THEN
       WRITE(*,*)"Error: ibeam_proj=/=1,2 in Pbreak1n_LO_AB_WS"
       STOP 1
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mN)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)
    
    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    ibeam_proj_common=ibeam_proj
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj
    R_proj_common=RR_proj
    w_proj_common=w_proj
    aa_proj_common=aa_proj
  
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init(ibeam_proj).EQ.0.OR.(integ_method6.LT.0))THEN
       R_proj_save(ibeam_proj)=RR_proj
       w_proj_save(ibeam_proj)=w_proj
       aa_proj_save(ibeam_proj)=aa_proj
       gamma_proj_save(ibeam_proj)=gamma_proj
       gamma_target_save(ibeam_proj)=gamma_target
       A_target_save(ibeam_proj)=A_target
       Z_target_save(ibeam_proj)=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak1n_LO_AB_WS=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak1n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak1n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints) 
             ENDIF
             Pbreak1n_LO_AB_WS=integral
             RETURN
          ENDIF
       ENDIF
       
       rescaling_bmax_save(ibeam_proj)=MAX(db,bmaxoR*RR_proj)
       
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(2,NXA))
       ENDIF
       
       dbb=9d0/DBLE(NBSEG)
       IF(init(3-ibeam_proj).EQ.1.AND.&
            rescaling_bmax_save(ibeam_proj).EQ.rescaling_bmax_save(3-ibeam_proj).AND.&
            R_proj_save(ibeam_proj).EQ.R_proj_save(3-ibeam_proj).AND.&
            w_proj_save(ibeam_proj).EQ.w_proj_save(3-ibeam_proj).AND.&
            aa_proj_save(ibeam_proj).EQ.aa_proj_save(3-ibeam_proj).AND.&
            gamma_proj_save(ibeam_proj).EQ.gamma_proj_save(3-ibeam_proj).AND.&
            gamma_target_save(ibeam_proj).EQ.gamma_target_save(3-ibeam_proj).AND.&
            A_target_save(ibeam_proj).EQ.A_target_save(3-ibeam_proj).AND.&
            Z_target_save(ibeam_proj).EQ.Z_target_save(3-ibeam_proj))THEN
          DO k=1,NXA
             XA(ibeam_proj,k)=XA(3-ibeam_proj,k)
             YA(ibeam_proj,k)=YA(3-ibeam_proj,k)
          ENDDO
       ELSE
          WRITE(*,*)"INFO: generate grid of LO Pbreak(1n) from charge form factor of beam_proj=",ibeam_proj
          WRITE(*,*)"INFO: it will take a few seconds !"
          k=0
          DO n=0,nbmax
             ! from 10**(-n-1)*bmax to 10**(-n)*bmax
             DO i=1,NBSEG
                k=NBSEG*n+i
                ! these are b in unit GeV-1
                XA(ibeam_proj,NXA-k+1)=(10d0**(-n-1))*(one&
                     +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save(ibeam_proj)
             ENDDO
          ENDDO
          IF(k+1.NE.NXA)THEN
             WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak1n_LO_AB_WS"
             STOP 1
          ENDIF
          XA(ibeam_proj,1)=0d0
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam_proj,I).EQ.0d0)THEN
                YA(ibeam_proj,I)=0d0
                CYCLE
             ENDIF

             npoints=30000
             b_proj_common=XA(ibeam_proj,I)
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak1n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak1n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             YA(ibeam_proj,I)=integral
          ENDDO
       ENDIF
       init(ibeam_proj)=1
    ENDIF
  
    IF(R_proj_save(ibeam_proj).NE.RR_proj&
         .OR.w_proj_save(ibeam_proj).NE.w_proj&
         .OR.aa_proj_save(ibeam_proj).NE.aa_proj)THEN
       WRITE(*,*)"ERROR: RA,wA,aA are not consistent in Pbreak1n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (RA,wA,aA)=",R_proj_save(ibeam_proj),&
            w_proj_save(ibeam_proj),aa_proj_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (RA,wA,aA)=",RR_proj,w_proj,aa_proj
       STOP 1
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak1n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak1n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF
  
    IF(db.GT.rescaling_bmax_save(ibeam_proj).OR.db.LE.0d0)THEN
       Pbreak1n_LO_AB_WS=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save(ibeam_proj),1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(ibeam_proj,k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(ibeam_proj,i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(ibeam_proj,K+I)
          YD2_1D(I)=YA(ibeam_proj,K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak1n_LO_AB_WS=integral
       
    ENDIF
    RETURN
  END FUNCTION Pbreak1n_LO_AB_WS

  FUNCTION Pbreak1n_LO_AB_WS_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak1n_LO_AB_WS_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak1n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak1n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif
    
    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak1n_LO_AB_WS_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mN*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_WS(b_proj_common,Ega_lab,&
            gamma_proj_common,R_proj_common,w_proj_common,aa_proj_common,&
            3d0,0.7d0,ibeam_proj_common)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak1n_LO_AB_WS_fxn=GAMMATO1N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak1n_LO_AB_WS_fxn=Pbreak1n_LO_AB_WS_fxn*CFF
    
    RETURN
  END FUNCTION Pbreak1n_LO_AB_WS_fxn

    FUNCTION Pbreak1n_LO_pA(Db,A_target,Z_target,gamma_target,&
       gamma_proj,Egamma_min_target,Egamma_max_target,&
       integ_method)
    ! It gives us the lowest-order probability P^{1n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak1n_LO_pA
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 1 MeV (0.001 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 40 MeV (0.04 GeV)
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak1n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak1n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common

    INTEGER::init=0
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::RA=35d0 ! in unit of GeV-1
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=5
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0))::gamma_proj_save,gamma_target_save
    INTEGER::A_target_save,Z_target_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save

    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D

    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mproton)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)

    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj

    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init.EQ.0.OR.integ_method6.LT.0)THEN
       gamma_proj_save=gamma_proj
       gamma_target_save=gamma_target
       A_target_save=A_target
       Z_target_save=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak1n_LO_pA=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
		Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak1n_LO_pA_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak1n_LO_pA_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             Pbreak1n_LO_pA=integral
             RETURN
          ENDIF
       ENDIF

       rescaling_bmax_save=MAX(db,bmaxoR*RA)

       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(NXA))
       ENDIF

       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(NXA))
       ENDIF

       dbb=9d0/DBLE(NBSEG)
       WRITE(*,*)"INFO: generate grid of LO Pbreak(1n) from charge form factor of pA"
       WRITE(*,*)"INFO: it will take a few seconds !"
       k=0
       DO n=0,nbmax
          ! from 10**(-n-1)*bmax to 10**(-n)*bmax
          DO i=1,NBSEG
             k=NBSEG*n+i
             ! these are b in unit GeV-1
             XA(NXA-k+1)=(10d0**(-n-1))*(one&
                  +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save
          ENDDO
       ENDDO
       IF(k+1.NE.NXA)THEN
          WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak1n_LO_pA"
          STOP 1
       ENDIF
       XA(1)=0d0
       DO I=1,NXA
          IF(I.EQ.1)THEN
             CALL progress(INT(I*50d0/NXA),50,.TRUE.)
          ELSE
             CALL progress(INT(I*50d0/NXA),50)
          ENDIF
          IF(XA(I).EQ.0d0)THEN
             YA(I)=0d0
             CYCLE
          ENDIF
          npoints=30000
          b_proj_common=XA(I)
          IF(Egamma_max_target.LT.0d0)THEN
             ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0       
             Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
             Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
             xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
             log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
          ENDIF
          IF(log1oxgamma_min.EQ.0d0)THEN
             CALL trapezoid_integration(npoints,Pbreak1n_LO_pA_fxn,&
                  log1oxgamma_max,integral)
          ELSE
             CALL simpson(Pbreak1n_LO_pA_fxn,log1oxgamma_min,&
                  log1oxgamma_max,integral,npoints)
          ENDIF
          YA(I)=integral
       ENDDO
       init=1
    ENDIF

    IF(gamma_proj_save.NE.gamma_proj&
         .OR.gamma_target_save.NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak1n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save,&
            gamma_target_save
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak1n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF

    IF(db.GT.rescaling_bmax_save.OR.db.LE.0d0)THEN
       Pbreak1n_LO_pA=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save,1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(K+I)
          YD2_1D(I)=YA(K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak1n_LO_pA=integral

    ENDIF
    RETURN
  END FUNCTION Pbreak1n_LO_pA

  FUNCTION Pbreak1n_LO_pA_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak1n_LO_pA_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak1n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak1n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif

    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak1n_LO_pA_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mproton*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_proton(b_proj_common,Ega_lab,&
            gamma_proj_common,proton_FF_correction,.FALSE.)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak1n_LO_pA_fxn=GAMMATO1N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak1n_LO_pA_fxn=Pbreak1n_LO_pA_fxn*CFF

    RETURN
  END FUNCTION Pbreak1n_LO_pA_fxn

  FUNCTION Pbreak2n_LO_AB_WS(Db,A_target,Z_target,gamma_target,&
       gamma_proj,RR_proj,w_proj,aa_proj,Egamma_min_target,Egamma_max_target,&
       ibeam_proj,integ_method)
    ! It gives us the lowest-order probability P^{2n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b,RR,aa should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ibeam_proj
    REAL(KIND(1d0))::Pbreak2n_LO_AB_WS
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::RR_proj,w_proj,aa_proj ! nucleu parameters for projectile
    
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 13 MeV (0.013 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 140 MeV (0.14 GeV).
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak2n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak2n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=4
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),DIMENSION(2)::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0)),DIMENSION(2)::R_proj_save,w_proj_save,aa_proj_save
    REAL(KIND(1d0)),DIMENSION(2)::gamma_proj_save,gamma_target_save
    INTEGER,DIMENSION(2)::A_target_save,Z_target_save
    SAVE R_proj_save,w_proj_save,aa_proj_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save
    
    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D
    
    IF(ibeam_proj.GT.2.OR.ibeam_proj.LT.1)THEN
       WRITE(*,*)"Error: ibeam_proj=/=1,2 in Pbreak2n_LO_AB_WS"
       STOP 1
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mN)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)
    
    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    ibeam_proj_common=ibeam_proj
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj
    R_proj_common=RR_proj
    w_proj_common=w_proj
    aa_proj_common=aa_proj
  
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init(ibeam_proj).EQ.0.OR.(integ_method6.LT.0))THEN
       R_proj_save(ibeam_proj)=RR_proj
       w_proj_save(ibeam_proj)=w_proj
       aa_proj_save(ibeam_proj)=aa_proj
       gamma_proj_save(ibeam_proj)=gamma_proj
       gamma_target_save(ibeam_proj)=gamma_target
       A_target_save(ibeam_proj)=A_target
       Z_target_save(ibeam_proj)=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak2n_LO_AB_WS=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak2n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak2n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints) 
             ENDIF
             Pbreak2n_LO_AB_WS=integral
             RETURN
          ENDIF
       ENDIF
       
       rescaling_bmax_save(ibeam_proj)=MAX(db,bmaxoR*RR_proj)
       
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(2,NXA))
       ENDIF
       
       dbb=9d0/DBLE(NBSEG)
       IF(init(3-ibeam_proj).EQ.1.AND.&
            rescaling_bmax_save(ibeam_proj).EQ.rescaling_bmax_save(3-ibeam_proj).AND.&
            R_proj_save(ibeam_proj).EQ.R_proj_save(3-ibeam_proj).AND.&
            w_proj_save(ibeam_proj).EQ.w_proj_save(3-ibeam_proj).AND.&
            aa_proj_save(ibeam_proj).EQ.aa_proj_save(3-ibeam_proj).AND.&
            gamma_proj_save(ibeam_proj).EQ.gamma_proj_save(3-ibeam_proj).AND.&
            gamma_target_save(ibeam_proj).EQ.gamma_target_save(3-ibeam_proj).AND.&
            A_target_save(ibeam_proj).EQ.A_target_save(3-ibeam_proj).AND.&
            Z_target_save(ibeam_proj).EQ.Z_target_save(3-ibeam_proj))THEN
          DO k=1,NXA
             XA(ibeam_proj,k)=XA(3-ibeam_proj,k)
             YA(ibeam_proj,k)=YA(3-ibeam_proj,k)
          ENDDO
       ELSE
          WRITE(*,*)"INFO: generate grid of LO Pbreak(2n) from charge form factor of beam_proj=",ibeam_proj
          WRITE(*,*)"INFO: it will take a few seconds !"
          k=0
          DO n=0,nbmax
             ! from 10**(-n-1)*bmax to 10**(-n)*bmax
             DO i=1,NBSEG
                k=NBSEG*n+i
                ! these are b in unit GeV-1
                XA(ibeam_proj,NXA-k+1)=(10d0**(-n-1))*(one&
                     +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save(ibeam_proj)
             ENDDO
          ENDDO
          IF(k+1.NE.NXA)THEN
             WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak2n_LO_AB_WS"
             STOP 1
          ENDIF
          XA(ibeam_proj,1)=0d0
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam_proj,I).EQ.0d0)THEN
                YA(ibeam_proj,I)=0d0
                CYCLE
             ENDIF

             npoints=30000
             b_proj_common=XA(ibeam_proj,I)
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak2n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak2n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             YA(ibeam_proj,I)=integral
          ENDDO
       ENDIF
       init(ibeam_proj)=1
    ENDIF
  
    IF(R_proj_save(ibeam_proj).NE.RR_proj&
         .OR.w_proj_save(ibeam_proj).NE.w_proj&
         .OR.aa_proj_save(ibeam_proj).NE.aa_proj)THEN
       WRITE(*,*)"ERROR: RA,wA,aA are not consistent in Pbreak2n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (RA,wA,aA)=",R_proj_save(ibeam_proj),&
            w_proj_save(ibeam_proj),aa_proj_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (RA,wA,aA)=",RR_proj,w_proj,aa_proj
       STOP 1
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak2n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak2n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF
  
    IF(db.GT.rescaling_bmax_save(ibeam_proj).OR.db.LE.0d0)THEN
       Pbreak2n_LO_AB_WS=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save(ibeam_proj),1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(ibeam_proj,k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(ibeam_proj,i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(ibeam_proj,K+I)
          YD2_1D(I)=YA(ibeam_proj,K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak2n_LO_AB_WS=integral
       
    ENDIF
    RETURN
  END FUNCTION Pbreak2n_LO_AB_WS

  FUNCTION Pbreak2n_LO_AB_WS_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak2n_LO_AB_WS_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak2n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak2n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif
    
    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak2n_LO_AB_WS_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mN*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_WS(b_proj_common,Ega_lab,&
            gamma_proj_common,R_proj_common,w_proj_common,aa_proj_common,&
            3d0,0.7d0,ibeam_proj_common)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak2n_LO_AB_WS_fxn=GAMMATO2N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak2n_LO_AB_WS_fxn=Pbreak2n_LO_AB_WS_fxn*CFF
    
    RETURN
  END FUNCTION Pbreak2n_LO_AB_WS_fxn

  FUNCTION Pbreak2n_LO_pA(Db,A_target,Z_target,gamma_target,&
       gamma_proj,Egamma_min_target,Egamma_max_target,&
       integ_method)
    ! It gives us the lowest-order probability P^{2n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak2n_LO_pA
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 14 MeV (0.014 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 140 MeV (0.14 GeV)
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak2n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak2n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common

    INTEGER::init=0
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::RA=35d0 ! in unit of GeV-1
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=5
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0))::gamma_proj_save,gamma_target_save
    INTEGER::A_target_save,Z_target_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save

    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D

    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mproton)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)

    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj

    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init.EQ.0.OR.integ_method6.LT.0)THEN
       gamma_proj_save=gamma_proj
       gamma_target_save=gamma_target
       A_target_save=A_target
       Z_target_save=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak2n_LO_pA=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
		Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak2n_LO_pA_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak2n_LO_pA_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             Pbreak2n_LO_pA=integral
             RETURN
          ENDIF
       ENDIF

       rescaling_bmax_save=MAX(db,bmaxoR*RA)

       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(NXA))
       ENDIF

       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(NXA))
       ENDIF

       dbb=9d0/DBLE(NBSEG)
       WRITE(*,*)"INFO: generate grid of LO Pbreak(2n) from charge form factor of pA"
       WRITE(*,*)"INFO: it will take a few seconds !"
       k=0
       DO n=0,nbmax
          ! from 10**(-n-1)*bmax to 10**(-n)*bmax
          DO i=1,NBSEG
             k=NBSEG*n+i
             ! these are b in unit GeV-1
             XA(NXA-k+1)=(10d0**(-n-1))*(one&
                  +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save
          ENDDO
       ENDDO
       IF(k+1.NE.NXA)THEN
          WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak2n_LO_pA"
          STOP 1
       ENDIF
       XA(1)=0d0
       DO I=1,NXA
          IF(I.EQ.1)THEN
             CALL progress(INT(I*50d0/NXA),50,.TRUE.)
          ELSE
             CALL progress(INT(I*50d0/NXA),50)
          ENDIF
          IF(XA(I).EQ.0d0)THEN
             YA(I)=0d0
             CYCLE
          ENDIF
          npoints=30000
          b_proj_common=XA(I)
          IF(Egamma_max_target.LT.0d0)THEN
             ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0       
             Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
             Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
             xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
             log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
          ENDIF
          IF(log1oxgamma_min.EQ.0d0)THEN
             CALL trapezoid_integration(npoints,Pbreak2n_LO_pA_fxn,&
                  log1oxgamma_max,integral)
          ELSE
             CALL simpson(Pbreak2n_LO_pA_fxn,log1oxgamma_min,&
                  log1oxgamma_max,integral,npoints)
          ENDIF
          YA(I)=integral
       ENDDO
       init=1
    ENDIF

    IF(gamma_proj_save.NE.gamma_proj&
         .OR.gamma_target_save.NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak2n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save,&
            gamma_target_save
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak2n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF

    IF(db.GT.rescaling_bmax_save.OR.db.LE.0d0)THEN
       Pbreak2n_LO_pA=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save,1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(K+I)
          YD2_1D(I)=YA(K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak2n_LO_pA=integral

    ENDIF
    RETURN
  END FUNCTION Pbreak2n_LO_pA

  FUNCTION Pbreak2n_LO_pA_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak2n_LO_pA_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak2n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak2n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif

    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak2n_LO_pA_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mproton*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_proton(b_proj_common,Ega_lab,&
            gamma_proj_common,proton_FF_correction,.FALSE.)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak2n_LO_pA_fxn=GAMMATO2N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak2n_LO_pA_fxn=Pbreak2n_LO_pA_fxn*CFF

    RETURN
  END FUNCTION Pbreak2n_LO_pA_fxn

  FUNCTION Pbreak3n_LO_AB_WS(Db,A_target,Z_target,gamma_target,&
       gamma_proj,RR_proj,w_proj,aa_proj,Egamma_min_target,Egamma_max_target,&
       ibeam_proj,integ_method)
    ! It gives us the lowest-order probability P^{3n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b,RR,aa should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ibeam_proj
    REAL(KIND(1d0))::Pbreak3n_LO_AB_WS
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::RR_proj,w_proj,aa_proj ! nucleu parameters for projectile
    
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 22 MeV (0.022 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 140 MeV (0.14 GeV).
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak3n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak3n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=4
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),DIMENSION(2)::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0)),DIMENSION(2)::R_proj_save,w_proj_save,aa_proj_save
    REAL(KIND(1d0)),DIMENSION(2)::gamma_proj_save,gamma_target_save
    INTEGER,DIMENSION(2)::A_target_save,Z_target_save
    SAVE R_proj_save,w_proj_save,aa_proj_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save
    
    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D
    
    IF(ibeam_proj.GT.2.OR.ibeam_proj.LT.1)THEN
       WRITE(*,*)"Error: ibeam_proj=/=1,2 in Pbreak3n_LO_AB_WS"
       STOP 1
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mN)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)
    
    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    ibeam_proj_common=ibeam_proj
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj
    R_proj_common=RR_proj
    w_proj_common=w_proj
    aa_proj_common=aa_proj
  
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init(ibeam_proj).EQ.0.OR.(integ_method6.LT.0))THEN
       R_proj_save(ibeam_proj)=RR_proj
       w_proj_save(ibeam_proj)=w_proj
       aa_proj_save(ibeam_proj)=aa_proj
       gamma_proj_save(ibeam_proj)=gamma_proj
       gamma_target_save(ibeam_proj)=gamma_target
       A_target_save(ibeam_proj)=A_target
       Z_target_save(ibeam_proj)=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak3n_LO_AB_WS=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak3n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak3n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints) 
             ENDIF
             Pbreak3n_LO_AB_WS=integral
             RETURN
          ENDIF
       ENDIF
       
       rescaling_bmax_save(ibeam_proj)=MAX(db,bmaxoR*RR_proj)
       
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(2,NXA))
       ENDIF
       
       dbb=9d0/DBLE(NBSEG)
       IF(init(3-ibeam_proj).EQ.1.AND.&
            rescaling_bmax_save(ibeam_proj).EQ.rescaling_bmax_save(3-ibeam_proj).AND.&
            R_proj_save(ibeam_proj).EQ.R_proj_save(3-ibeam_proj).AND.&
            w_proj_save(ibeam_proj).EQ.w_proj_save(3-ibeam_proj).AND.&
            aa_proj_save(ibeam_proj).EQ.aa_proj_save(3-ibeam_proj).AND.&
            gamma_proj_save(ibeam_proj).EQ.gamma_proj_save(3-ibeam_proj).AND.&
            gamma_target_save(ibeam_proj).EQ.gamma_target_save(3-ibeam_proj).AND.&
            A_target_save(ibeam_proj).EQ.A_target_save(3-ibeam_proj).AND.&
            Z_target_save(ibeam_proj).EQ.Z_target_save(3-ibeam_proj))THEN
          DO k=1,NXA
             XA(ibeam_proj,k)=XA(3-ibeam_proj,k)
             YA(ibeam_proj,k)=YA(3-ibeam_proj,k)
          ENDDO
       ELSE
          WRITE(*,*)"INFO: generate grid of LO Pbreak(3n) from charge form factor of beam_proj=",ibeam_proj
          WRITE(*,*)"INFO: it will take a few seconds !"
          k=0
          DO n=0,nbmax
             ! from 10**(-n-1)*bmax to 10**(-n)*bmax
             DO i=1,NBSEG
                k=NBSEG*n+i
                ! these are b in unit GeV-1
                XA(ibeam_proj,NXA-k+1)=(10d0**(-n-1))*(one&
                     +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save(ibeam_proj)
             ENDDO
          ENDDO
          IF(k+1.NE.NXA)THEN
             WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak3n_LO_AB_WS"
             STOP 1
          ENDIF
          XA(ibeam_proj,1)=0d0
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam_proj,I).EQ.0d0)THEN
                YA(ibeam_proj,I)=0d0
                CYCLE
             ENDIF

             npoints=30000
             b_proj_common=XA(ibeam_proj,I)
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak3n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak3n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             YA(ibeam_proj,I)=integral
          ENDDO
       ENDIF
       init(ibeam_proj)=1
    ENDIF
  
    IF(R_proj_save(ibeam_proj).NE.RR_proj&
         .OR.w_proj_save(ibeam_proj).NE.w_proj&
         .OR.aa_proj_save(ibeam_proj).NE.aa_proj)THEN
       WRITE(*,*)"ERROR: RA,wA,aA are not consistent in Pbreak3n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (RA,wA,aA)=",R_proj_save(ibeam_proj),&
            w_proj_save(ibeam_proj),aa_proj_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (RA,wA,aA)=",RR_proj,w_proj,aa_proj
       STOP 1
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak3n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak3n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF
  
    IF(db.GT.rescaling_bmax_save(ibeam_proj).OR.db.LE.0d0)THEN
       Pbreak3n_LO_AB_WS=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save(ibeam_proj),1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(ibeam_proj,k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(ibeam_proj,i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(ibeam_proj,K+I)
          YD2_1D(I)=YA(ibeam_proj,K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak3n_LO_AB_WS=integral
       
    ENDIF
    RETURN
  END FUNCTION Pbreak3n_LO_AB_WS

  FUNCTION Pbreak3n_LO_AB_WS_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak3n_LO_AB_WS_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak3n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak3n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif
    
    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak3n_LO_AB_WS_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mN*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_WS(b_proj_common,Ega_lab,&
            gamma_proj_common,R_proj_common,w_proj_common,aa_proj_common,&
            3d0,0.7d0,ibeam_proj_common)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak3n_LO_AB_WS_fxn=GAMMATO3N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak3n_LO_AB_WS_fxn=Pbreak3n_LO_AB_WS_fxn*CFF
    
    RETURN
  END FUNCTION Pbreak3n_LO_AB_WS_fxn

  FUNCTION Pbreak3n_LO_pA(Db,A_target,Z_target,gamma_target,&
       gamma_proj,Egamma_min_target,Egamma_max_target,&
       integ_method)
    ! It gives us the lowest-order probability P^{3n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak3n_LO_pA
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 22 MeV (0.022 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 140 MeV (0.14 GeV)
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak3n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak3n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common

    INTEGER::init=0
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::RA=35d0 ! in unit of GeV-1
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=5
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0))::gamma_proj_save,gamma_target_save
    INTEGER::A_target_save,Z_target_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save

    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D

    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mproton)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)

    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj

    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init.EQ.0.OR.integ_method6.LT.0)THEN
       gamma_proj_save=gamma_proj
       gamma_target_save=gamma_target
       A_target_save=A_target
       Z_target_save=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak3n_LO_pA=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
		Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak3n_LO_pA_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak3n_LO_pA_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             Pbreak3n_LO_pA=integral
             RETURN
          ENDIF
       ENDIF

       rescaling_bmax_save=MAX(db,bmaxoR*RA)

       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(NXA))
       ENDIF

       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(NXA))
       ENDIF

       dbb=9d0/DBLE(NBSEG)
       WRITE(*,*)"INFO: generate grid of LO Pbreak(3n) from charge form factor of pA"
       WRITE(*,*)"INFO: it will take a few seconds !"
       k=0
       DO n=0,nbmax
          ! from 10**(-n-1)*bmax to 10**(-n)*bmax
          DO i=1,NBSEG
             k=NBSEG*n+i
             ! these are b in unit GeV-1
             XA(NXA-k+1)=(10d0**(-n-1))*(one&
                  +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save
          ENDDO
       ENDDO
       IF(k+1.NE.NXA)THEN
          WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak3n_LO_pA"
          STOP 1
       ENDIF
       XA(1)=0d0
       DO I=1,NXA
          IF(I.EQ.1)THEN
             CALL progress(INT(I*50d0/NXA),50,.TRUE.)
          ELSE
             CALL progress(INT(I*50d0/NXA),50)
          ENDIF
          IF(XA(I).EQ.0d0)THEN
             YA(I)=0d0
             CYCLE
          ENDIF
          npoints=30000
          b_proj_common=XA(I)
          IF(Egamma_max_target.LT.0d0)THEN
             ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0       
             Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
             Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
             xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
             log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
          ENDIF
          IF(log1oxgamma_min.EQ.0d0)THEN
             CALL trapezoid_integration(npoints,Pbreak3n_LO_pA_fxn,&
                  log1oxgamma_max,integral)
          ELSE
             CALL simpson(Pbreak3n_LO_pA_fxn,log1oxgamma_min,&
                  log1oxgamma_max,integral,npoints)
          ENDIF
          YA(I)=integral
       ENDDO
       init=1
    ENDIF

    IF(gamma_proj_save.NE.gamma_proj&
         .OR.gamma_target_save.NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak3n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save,&
            gamma_target_save
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak3n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF

    IF(db.GT.rescaling_bmax_save.OR.db.LE.0d0)THEN
       Pbreak3n_LO_pA=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save,1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(K+I)
          YD2_1D(I)=YA(K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak3n_LO_pA=integral

    ENDIF
    RETURN
  END FUNCTION Pbreak3n_LO_pA

  FUNCTION Pbreak3n_LO_pA_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak3n_LO_pA_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak3n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak3n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif

    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak3n_LO_pA_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mproton*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_proton(b_proj_common,Ega_lab,&
            gamma_proj_common,proton_FF_correction,.FALSE.)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak3n_LO_pA_fxn=GAMMATO3N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak3n_LO_pA_fxn=Pbreak3n_LO_pA_fxn*CFF

    RETURN
  END FUNCTION Pbreak3n_LO_pA_fxn

  FUNCTION Pbreak4n_LO_AB_WS(Db,A_target,Z_target,gamma_target,&
       gamma_proj,RR_proj,w_proj,aa_proj,Egamma_min_target,Egamma_max_target,&
       ibeam_proj,integ_method)
    ! It gives us the lowest-order probability P^{4n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b,RR,aa should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ibeam_proj
    REAL(KIND(1d0))::Pbreak4n_LO_AB_WS
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::RR_proj,w_proj,aa_proj ! nucleu parameters for projectile
    
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 32 MeV (0.032 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 140 MeV (0.14 GeV).
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak4n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak4n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=4
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),DIMENSION(2)::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0)),DIMENSION(2)::R_proj_save,w_proj_save,aa_proj_save
    REAL(KIND(1d0)),DIMENSION(2)::gamma_proj_save,gamma_target_save
    INTEGER,DIMENSION(2)::A_target_save,Z_target_save
    SAVE R_proj_save,w_proj_save,aa_proj_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save
    
    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D
    
    IF(ibeam_proj.GT.2.OR.ibeam_proj.LT.1)THEN
       WRITE(*,*)"Error: ibeam_proj=/=1,2 in Pbreak4n_LO_AB_WS"
       STOP 1
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mN)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)
    
    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    ibeam_proj_common=ibeam_proj
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj
    R_proj_common=RR_proj
    w_proj_common=w_proj
    aa_proj_common=aa_proj
  
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init(ibeam_proj).EQ.0.OR.(integ_method6.LT.0))THEN
       R_proj_save(ibeam_proj)=RR_proj
       w_proj_save(ibeam_proj)=w_proj
       aa_proj_save(ibeam_proj)=aa_proj
       gamma_proj_save(ibeam_proj)=gamma_proj
       gamma_target_save(ibeam_proj)=gamma_target
       A_target_save(ibeam_proj)=A_target
       Z_target_save(ibeam_proj)=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak4n_LO_AB_WS=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak4n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak4n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints) 
             ENDIF
             Pbreak4n_LO_AB_WS=integral
             RETURN
          ENDIF
       ENDIF
       
       rescaling_bmax_save(ibeam_proj)=MAX(db,bmaxoR*RR_proj)
       
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(2,NXA))
       ENDIF
       
       dbb=9d0/DBLE(NBSEG)
       IF(init(3-ibeam_proj).EQ.1.AND.&
            rescaling_bmax_save(ibeam_proj).EQ.rescaling_bmax_save(3-ibeam_proj).AND.&
            R_proj_save(ibeam_proj).EQ.R_proj_save(3-ibeam_proj).AND.&
            w_proj_save(ibeam_proj).EQ.w_proj_save(3-ibeam_proj).AND.&
            aa_proj_save(ibeam_proj).EQ.aa_proj_save(3-ibeam_proj).AND.&
            gamma_proj_save(ibeam_proj).EQ.gamma_proj_save(3-ibeam_proj).AND.&
            gamma_target_save(ibeam_proj).EQ.gamma_target_save(3-ibeam_proj).AND.&
            A_target_save(ibeam_proj).EQ.A_target_save(3-ibeam_proj).AND.&
            Z_target_save(ibeam_proj).EQ.Z_target_save(3-ibeam_proj))THEN
          DO k=1,NXA
             XA(ibeam_proj,k)=XA(3-ibeam_proj,k)
             YA(ibeam_proj,k)=YA(3-ibeam_proj,k)
          ENDDO
       ELSE
          WRITE(*,*)"INFO: generate grid of LO Pbreak(4n) from charge form factor of beam_proj=",ibeam_proj
          WRITE(*,*)"INFO: it will take a few seconds !"
          k=0
          DO n=0,nbmax
             ! from 10**(-n-1)*bmax to 10**(-n)*bmax
             DO i=1,NBSEG
                k=NBSEG*n+i
                ! these are b in unit GeV-1
                XA(ibeam_proj,NXA-k+1)=(10d0**(-n-1))*(one&
                     +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save(ibeam_proj)
             ENDDO
          ENDDO
          IF(k+1.NE.NXA)THEN
             WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak4n_LO_AB_WS"
             STOP 1
          ENDIF
          XA(ibeam_proj,1)=0d0
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam_proj,I).EQ.0d0)THEN
                YA(ibeam_proj,I)=0d0
                CYCLE
             ENDIF

             npoints=30000
             b_proj_common=XA(ibeam_proj,I)
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
                Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mN),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak4n_LO_AB_WS_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak4n_LO_AB_WS_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             YA(ibeam_proj,I)=integral
          ENDDO
       ENDIF
       init(ibeam_proj)=1
    ENDIF
  
    IF(R_proj_save(ibeam_proj).NE.RR_proj&
         .OR.w_proj_save(ibeam_proj).NE.w_proj&
         .OR.aa_proj_save(ibeam_proj).NE.aa_proj)THEN
       WRITE(*,*)"ERROR: RA,wA,aA are not consistent in Pbreak4n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (RA,wA,aA)=",R_proj_save(ibeam_proj),&
            w_proj_save(ibeam_proj),aa_proj_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (RA,wA,aA)=",RR_proj,w_proj,aa_proj
       STOP 1
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak4n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak4n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF
  
    IF(db.GT.rescaling_bmax_save(ibeam_proj).OR.db.LE.0d0)THEN
       Pbreak4n_LO_AB_WS=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save(ibeam_proj),1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(ibeam_proj,k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(ibeam_proj,i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(ibeam_proj,K+I)
          YD2_1D(I)=YA(ibeam_proj,K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak4n_LO_AB_WS=integral
       
    ENDIF
    RETURN
  END FUNCTION Pbreak4n_LO_AB_WS

  FUNCTION Pbreak4n_LO_AB_WS_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak4n_LO_AB_WS_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    INTEGER::A_target_common,Z_target_common,ibeam_proj_common
    COMMON/Pbreak4n_AB_CFF_WS_INT/A_target_common,Z_target_common,ibeam_proj_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common,R_proj_common,w_proj_common,aa_proj_common
    COMMON/Pbreak4n_AB_CFF_WS_REAL/gamma_target_common,gamma_proj_common,b_proj_common,&
         R_proj_common,w_proj_common,aa_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif
    
    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak4n_LO_AB_WS_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mN*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_WS(b_proj_common,Ega_lab,&
            gamma_proj_common,R_proj_common,w_proj_common,aa_proj_common,&
            3d0,0.7d0,ibeam_proj_common)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak4n_LO_AB_WS_fxn=GAMMATO4N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak4n_LO_AB_WS_fxn=Pbreak4n_LO_AB_WS_fxn*CFF
    
    RETURN
  END FUNCTION Pbreak4n_LO_AB_WS_fxn

  FUNCTION Pbreak4n_LO_pA(Db,A_target,Z_target,gamma_target,&
       gamma_proj,Egamma_min_target,Egamma_max_target,&
       integ_method)
    ! It gives us the lowest-order probability P^{4n,(1)}
    ! for Coulomb breakup of a nucleus with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak4n_LO_pA
    REAL(KIND(1d0)),INTENT(IN)::Db ! Delta b
    INTEGER,INTENT(IN)::A_target,Z_target ! atomic mass and atomic numbers for the target
    ! note Lorentz gamma should be >= 1
    REAL(KIND(1d0)),INTENT(IN)::gamma_proj ! Lorentz gamma for projectile in lab frame
    REAL(KIND(1d0)),INTENT(IN)::gamma_target ! Lorentz gamma for target in lab frame
    REAL(KIND(1d0)),INTENT(IN)::Egamma_min_target ! Minimal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 32 MeV (0.032 GeV).
    REAL(KIND(1d0)),INTENT(IN)::Egamma_max_target ! Maximal Egamma in the target frame in unit of GeV
                                                  ! A nominal value is 140 MeV (0.14 GeV)
    ! If it is negative, we take Abs(Egamma_max_target)*(gamma_proj_target/Db),
    ! where gamma_proj_target=gamma_proj*gamma_target+sqrt((gamma_proj**2-1)*(gamma_target**2-1))
    ! which is the Lorentz gamma of projectile in the target frame
    ! If it is zero, xgamma upper limit is 1
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::gamma_proj_target
    REAL(KIND(1d0))::Egamma_min_lab,Egamma_max_lab
    REAL(KIND(1d0))::xgamma_min_lab,xgamma_max_lab
    REAL(KIND(1d0))::log1oxgamma_max, log1oxgamma_min
    REAL(KIND(1d0))::integral
    
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak4n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak4n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common

    INTEGER::init=0
    INTEGER::NXA,i,j,n,k,l
    SAVE init,NXA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::YA
    SAVE XA,YA
    REAL(KIND(1d0)),PARAMETER::RA=35d0 ! in unit of GeV-1
    REAL(KIND(1d0)),PARAMETER::bmaxoR=1d3
    INTEGER,PARAMETER::NBMAX=5
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=50
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::rescaling_bmax_save
    REAL(KIND(1d0))::dbb,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0))::gamma_proj_save,gamma_target_save
    INTEGER::A_target_save,Z_target_save
    SAVE gamma_proj_save,gamma_target_save
    SAVE A_target_save, Z_target_save

    INTEGER,PARAMETER::n_interp=10
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::YD2_1D

    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP 1
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP 1
    ENDIF
    gamma_proj_target=gamma_proj*gamma_target&
         +dsqrt((gamma_proj**2-one)*(gamma_target**2-one))

    ! lowest gamma energy [GeV] in the lab frame
    Egamma_min_lab=MAX(Egamma_min_target,0d0)/(gamma_target+dsqrt(gamma_target**2-one))
    xgamma_min_lab=Egamma_min_lab/(gamma_proj*mproton)
    log1oxgamma_max=DLOG(1d0/xgamma_min_lab)

    IF(Egamma_max_target.GT.0d0)THEN
       IF(Egamma_max_target.LE.Egamma_min_target)THEN
          WRITE(*,*)"Error: Egamma_max_target <= Egamma_min_target"
          STOP 1
       ENDIF
       ! highest gamma energy [GeV] in the lab frame if Egamma_max_target > 0
       Egamma_max_lab=Egamma_max_target/(gamma_target+dsqrt(gamma_target**2-one))
       xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
       log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
    ELSE
       log1oxgamma_min=0d0
    ENDIF

    A_target_common=A_target
    Z_target_common=Z_target
    gamma_target_common=gamma_target
    gamma_proj_common=gamma_proj

    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF

    IF(init.EQ.0.OR.integ_method6.LT.0)THEN
       gamma_proj_save=gamma_proj
       gamma_target_save=gamma_target
       A_target_save=A_target
       Z_target_save=Z_target

       IF(integ_method6.LT.0)THEN
          IF(db.EQ.0d0)THEN
             Pbreak4n_LO_pA=0d0
             RETURN
          ELSE
             npoints=30000
             b_proj_common=db
             IF(Egamma_max_target.LT.0d0)THEN
                ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0
                Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/db
		Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
                xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
                log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
             ENDIF
             IF(log1oxgamma_min.EQ.0d0)THEN
                CALL trapezoid_integration(npoints,Pbreak4n_LO_pA_fxn,&
                     log1oxgamma_max,integral)
             ELSE
                CALL simpson(Pbreak4n_LO_pA_fxn,log1oxgamma_min,&
                     log1oxgamma_max,integral,npoints)
             ENDIF
             Pbreak4n_LO_pA=integral
             RETURN
          ENDIF
       ENDIF

       rescaling_bmax_save=MAX(db,bmaxoR*RA)

       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(NXA))
       ENDIF

       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(NXA))
       ENDIF

       dbb=9d0/DBLE(NBSEG)
       WRITE(*,*)"INFO: generate grid of LO Pbreak(4n) from charge form factor of pA"
       WRITE(*,*)"INFO: it will take a few seconds !"
       k=0
       DO n=0,nbmax
          ! from 10**(-n-1)*bmax to 10**(-n)*bmax
          DO i=1,NBSEG
             k=NBSEG*n+i
             ! these are b in unit GeV-1
             XA(NXA-k+1)=(10d0**(-n-1))*(one&
                  +DBLE(NBSEG+1-i)*dbb)*rescaling_bmax_save
          ENDDO
       ENDDO
       IF(k+1.NE.NXA)THEN
          WRITE(*,*)"ERROR: mismatching k+1 and NXA in Pbreak4n_LO_pA"
          STOP 1
       ENDIF
       XA(1)=0d0
       DO I=1,NXA
          IF(I.EQ.1)THEN
             CALL progress(INT(I*50d0/NXA),50,.TRUE.)
          ELSE
             CALL progress(INT(I*50d0/NXA),50)
          ENDIF
          IF(XA(I).EQ.0d0)THEN
             YA(I)=0d0
             CYCLE
          ENDIF
          npoints=30000
          b_proj_common=XA(I)
          IF(Egamma_max_target.LT.0d0)THEN
             ! highest gamma energy [GeV] in the lab frame if Egamma_max_target < 0       
             Egamma_max_lab=DABS(Egamma_max_target)*gamma_proj_target/b_proj_common
             Egamma_max_lab=Egamma_max_lab/(gamma_target+dsqrt(gamma_target**2-one))
             xgamma_max_lab=MIN(Egamma_max_lab/(gamma_proj*mproton),1d0)
             log1oxgamma_min=DLOG(1d0/xgamma_max_lab)
          ENDIF
          IF(log1oxgamma_min.EQ.0d0)THEN
             CALL trapezoid_integration(npoints,Pbreak4n_LO_pA_fxn,&
                  log1oxgamma_max,integral)
          ELSE
             CALL simpson(Pbreak4n_LO_pA_fxn,log1oxgamma_min,&
                  log1oxgamma_max,integral,npoints)
          ENDIF
          YA(I)=integral
       ENDDO
       init=1
    ENDIF

    IF(gamma_proj_save.NE.gamma_proj&
         .OR.gamma_target_save.NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak4n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save,&
            gamma_target_save
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP 1
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak4n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP 1
    ENDIF

    IF(db.GT.rescaling_bmax_save.OR.db.LE.0d0)THEN
       Pbreak4n_LO_pA=0d0
    ELSE
       
       dbb=MIN(db/rescaling_bmax_save,1d0)
       N=-FLOOR(DLOG10(dbb))-1 ! db is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(k-NBSEG).GT.db)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(i).GE.db)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
       ENDIF

       DO I=1,n_interp
          XD2_1D(I)=XA(K+I)
          YD2_1D(I)=YA(K+I)
       ENDDO

       CALL SPLINE_INTERPOLATE(XD2_1D,YD2_1D,n_interp,db,integral)
       Pbreak4n_LO_pA=integral

    ENDIF
    RETURN
  END FUNCTION Pbreak4n_LO_pA

  FUNCTION Pbreak4n_LO_pA_fxn(log1oxgamma)
    ! log1oxgamma = log(1/xgamma)
    ! where xgamma is the momentum fraction in the lab frame
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::Pbreak4n_LO_pA_fxn
    REAL(KIND(1d0)),INTENT(IN)::log1oxgamma
    REAL(KIND(1d0))::xgamma,Ega_lab,Ega_target,CFF
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::A_target_common,Z_target_common
    COMMON/Pbreak4n_pA_CFF_INT/A_target_common,Z_target_common
    REAL(KIND(1d0))::gamma_proj_common,gamma_target_common
    REAL(KIND(1d0))::b_proj_common
    COMMON/Pbreak4n_pA_CFF_REAL/gamma_target_common,gamma_proj_common,b_proj_common
    INTEGER,SAVE::init=0

    if(init.eq.0)then
       ! uncertainty for neutron tagging
       XSIGMA=NEUTRON_XSIGMA
       init=1
    endif

    xgamma=DEXP(-log1oxgamma)
    IF(xgamma.GE.1d0.OR.xgamma.LE.0d0)THEN
       Pbreak4n_LO_pA_fxn=0d0
       RETURN
    ENDIF

    ! gamma energy in lab [in unit of GeV]
    Ega_lab=xgamma*mproton*gamma_proj_common

    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       CFF=PhotonNumberDensity_ChargeFormFactor_proton(b_proj_common,Ega_lab,&
            gamma_proj_common,proton_FF_correction,.FALSE.)
    ELSE
       CFF=PhotonNumberDensity(b_proj_common,Ega_lab,gamma_proj_common)
    ENDIF

    ! gamma energy in target frame [in unit of MeV]
    Ega_target=Ega_lab*(gamma_target_common+DSQRT(gamma_target_common**2-1d0))
    Ega_target=Ega_target*1d3 ! from GeV to MeV

    ! in unit of GeV**(-2)
    Pbreak4n_LO_pA_fxn=GAMMATO4N(Ega_target,&
         A_target_common,Z_target_common)*0.1d0/GeVm12fm**2

    Pbreak4n_LO_pA_fxn=Pbreak4n_LO_pA_fxn*CFF

    RETURN
  END FUNCTION Pbreak4n_LO_pA_fxn
  

  ! erf^-1(x)
  ! from http://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
  ! there is a built-in Fortran function of erf (derf)
  function erfinv(x)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    REAL(KIND(1d0))::erfinv
    REAL(KIND(1d0))::w,p
    IF(DABS(x).GE.1d0)THEN
       WRITE(*,*)"ERROR: |x|>= 1 in erfinv"
       STOP 1
    ENDIF
    w=-DLOG(1d0-x**2)
    IF(w.LT.5d0)THEN
       w=w-2.5d0
       p=2.81022636d-8
       p=3.43273939d-7+p*w
       p=-3.5233877d-6+p*w
       p=-4.39150654d-6+p*w
       p=0.00021858087d0+p*w
       p=-0.00125372503d0+p*w
       p=-0.00417768164d0+p*w
       p=0.246640727d0+p*w
       p=1.50140941d0+p*w
    ELSE
       w=dsqrt(w)-3d0
       p=-0.000200214257d0
       p=0.000100950558d0+p*w
       p=0.00134934322d0+p*w
       p=-0.00367342844d0+p*w
       p=0.00573950773d0+p*w
       p=-0.0076224613d0+p*w
       p=0.00943887047d0+p*w
       p=1.00167406d0+p*w
       p=2.83297682d0+p*w
    ENDIF
    erfinv=p*x
    RETURN
  END function erfinv

  SUBROUTINE InitialMomentumReshuffle_PhotonPhoton(NP,PA,Ecm,x1,x2,Q1,Q2,&
       PO,reshuffled)
    ! work in the rest frame of the two colliding nucleons
    IMPLICIT NONE
    integer,intent(in)::NP
    ! Ecm**2=s_{NN}
    real(kind(1d0)),intent(in)::Ecm,x1,x2,Q1,Q2
    real(kind(1d0)),dimension(NP,0:3),intent(in)::PA
    real(kind(1d0)),dimension(NP,0:3),intent(out)::PO
    LOGICAL,INTENT(out)::reshuffled
    real(kind(1d0))::r1,r2,ph1,ph2,x1bar,x2bar
    real(kind(1d0)),parameter::pipi=3.14159265358979323846264338328d0
    real(kind(1d0))::kperp2,kperp2max,Q
    real(kind(1d0)),dimension(0:3)::PBOO1,PBOO2
    integer::j
    reshuffled=.FALSE. ! if reshuffled is .False., we need to regenerate Q1 and Q2
    r1=rand() ! a random float number between 0 and 1
    r2=rand() ! a random float number between 0 and 1
    ph1=2d0*pipi*r1
    ph2=2d0*pipi*r2
    kperp2=Q1**2+Q2**2+2d0*Q1*Q2*DCOS(ph1-ph2)
    kperp2max=Ecm**2*(MIN(1d0,x1/x2,x2/x1)-x1*x2)
    IF(kperp2.GE.kperp2max)RETURN
    x1bar=DSQRT(x1/x2*kperp2/Ecm**2+x1**2)
    x2bar=DSQRT(x2/x1*kperp2/Ecm**2+x2**2)
    IF(x1bar.GE.1d0.or.x2bar.GE.1d0)RETURN
    ! new initial state
    PO(1,0)=Ecm/2d0*x1bar
    PO(1,1)=Q1*DCOS(ph1)
    PO(1,2)=Q1*DSIN(ph1)
    PO(1,3)=Ecm/2d0*x1bar
    PO(2,0)=Ecm/2d0*x2bar
    PO(2,1)=Q2*DCOS(ph2)
    PO(2,2)=Q2*DSIN(ph2)
    PO(2,3)=-Ecm/2d0*x2bar
    ! new final state
    DO j=0,3
       PBOO1(j)=PA(1,j)+PA(2,j)
       PBOO2(j)=PO(1,j)+PO(2,j)
    ENDDO
    Q=DSQRT(x1*x2)*Ecm
    DO j=3,NP
       PO(j,0:3)=PA(j,0:3)
       CALL BOOSTL52(Q,PBOO1,PBOO2,PO(j,0:3))
    ENDDO
    reshuffled=.TRUE.
    RETURN
  END SUBROUTINE InitialMomentumReshuffle_PhotonPhoton

  ! When PBOO(0)=P(0) and PBOO(1:3)=-P(1:3), it boosts P to its rest frame (0,0,0,Q)
  SUBROUTINE BOOSTL5(Q,PBOO,P)
    ! momentums are in normal representation with fourth comp is the zero comp
    ! Boost P via PBOO(PBOO^2=Q^2) to PLB,and set to P
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::Q
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PBOO
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(INOUT)::P
    REAL(KIND(1d0)),DIMENSION(0:3)::PCM,PLB
    REAL(KIND(1d0))::FACT
    INTEGER::J
    PCM(0:3)=P(0:3)
    PLB(0)=(PBOO(0)*PCM(0)&
         +PBOO(3)*PCM(3)+PBOO(2)*PCM(2)+PBOO(1)*PCM(1))/Q
    FACT=(PLB(0)+PCM(0))/(Q+PBOO(0))
    DO J=1,3
       PLB(J)=PCM(J)+FACT*PBOO(J)
    ENDDO
    P(0:3)=PLB(0:3)
    RETURN
  END SUBROUTINE BOOSTL5

  SUBROUTINE BOOSTL52(Q,PBOO1,PBOO2,P)
    ! Boost P from PBOO1 (PBOO1^2=Q^2) to PBOO2 (PBOO2^2=Q^2) frame
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::Q
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PBOO1,PBOO2
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(INOUT)::P
    REAL(KIND(1d0)),DIMENSION(0:3)::PBOO10,PRES
    PBOO10(0)=PBOO1(0)
    PBOO10(1:3)=-PBOO1(1:3)
    PRES(0:3)=P(0:3)
    CALL BOOSTL5(Q,PBOO10,PRES) ! PRES is in (Q,0,0,0) frame
    CALL BOOSTL5(Q,PBOO2,PRES)  ! PRES is in PBOO2 frame
    P(0:3)=PRES(0:3)
    RETURN
  END SUBROUTINE BOOSTL52

  SUBROUTINE BOOSTTOECM(NP,E1,E2,PA)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NP
    REAL(KIND(1d0)),INTENT(IN)::E1,E2
    REAL(KIND(1d0)),DIMENSION(NP,0:3),INTENT(INOUT)::PA
    REAL(KIND(1d0))::Ecm
    REAL(KIND(1d0)),DIMENSION(0:3)::PBOO
    INTEGER::J
    Ecm=2d0*DSQRT(E1*E2)
    PBOO(0)=E1+E2
    PBOO(1:2)=0d0
    PBOO(3)=E2-E1
    DO J=1,NP
       CALL BOOSTL5(Ecm,PBOO,PA(J,0:3))
    ENDDO
    RETURN
  END SUBROUTINE BOOSTTOECM

  SUBROUTINE BOOSTFROMECM(NP,E1,E2,PA)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NP
    REAL(KIND(1d0)),INTENT(IN)::E1,E2
    REAL(KIND(1d0)),DIMENSION(NP,0:3),INTENT(INOUT)::PA
    REAL(KIND(1d0))::Ecm
    REAL(KIND(1d0)),DIMENSION(0:3)::PBOO
    INTEGER::J
    Ecm=2d0*DSQRT(E1*E2)
    PBOO(0)=E1+E2
    PBOO(1:2)=0d0
    PBOO(3)=E1-E2
    DO J=1,NP
       CALL BOOSTL5(Ecm,PBOO,PA(J,0:3))
    ENDDO
    RETURN
  END SUBROUTINE BOOSTFROMECM

  FUNCTION PhotonNumberDensity_proton(b,Ega,gamma,ibeam)
    ! the proton photon number density with ChFF (Z=1, alpha=1): small-x approximation
    ! (PhotonNumberDensity_ChargeFormFactor_proton) or not (USE_EXACT_ELASTIC4PHOTON),
    ! with the form factors of proton_FF_correction; the last value of each beam (ibeam) is kept
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_proton
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    INTEGER,INTENT(IN)::ibeam
    REAL(KIND(1d0)),DIMENSION(4,2)::last=-1d0 ! b,Ega,gamma,value
    INTEGER,DIMENSION(2)::lastset=-1 ! the switches of last
    SAVE last,lastset
    INTEGER::set
    set=proton_FF_correction+MERGE(10,0,USE_EXACT_ELASTIC4PHOTON)+MERGE(100,0,alpha_running4photon)
    IF(ALL(last(1:3,ibeam).EQ.(/b,Ega,gamma/)).AND.lastset(ibeam).EQ.set)THEN
       PhotonNumberDensity_proton=last(4,ibeam)
       RETURN
    ENDIF
    IF(USE_EXACT_ELASTIC4PHOTON)THEN
       PhotonNumberDensity_proton=PhotonNumberDensity_ChargeFormFactor_proton_exact(b,Ega,gamma,proton_FF_correction)
    ELSE
       PhotonNumberDensity_proton=PhotonNumberDensity_ChargeFormFactor_proton(b,Ega,gamma,proton_FF_correction)
    ENDIF
    last(:,ibeam)=(/b,Ega,gamma,PhotonNumberDensity_proton/)
    lastset(ibeam)=set
    RETURN
  END FUNCTION PhotonNumberDensity_proton

  FUNCTION PhotonFlux_proton_nob(x,gamma)
    ! the proton photon flux with ChFF and PNOHAD=1, as PhotonNumberDensity_proton (it depends on x only);
    ! with GENERATE_PROTON_FLUX_GRID, in the exact distribution or with running alpha (numerical integrals), for
    ! XGMIN<=x<=XGMAX from a cubic spline of log(x*dN/dx) in log(x/(1-x)) (NG nodes; log(x) at small x, log(1-x)
    ! near 1, where the flux vanishes), built at the first call with the current switches; for XGMAX<x<=XHMAX from a
    ! second spline (NH nodes, half the spacing, from NHO nodes below XGMAX to NHO nodes above XHMAX, so that its
    ! end conditions are outside the range where it is used; grid/direct-1 <= 2e-6, all switches), built when first
    ! needed; direct otherwise (x > XHMAX: only the end point x = 1 of the luminosity integrals)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_nob
    REAL(KIND(1d0)),INTENT(IN)::x,gamma
    INTEGER,PARAMETER::NG=361 ! 40 nodes per decade at small x
    REAL(KIND(1d0)),PARAMETER::XGMIN=1d-9,XGMAX=0.9d0,XHMAX=0.9999d0
    INTEGER,PARAMETER::NHO=10,NH=242 ! (logit(XHMAX)-logit(XGMAX))/(spacing/2) = 220.3, +2 NHO, +1
    REAL(KIND(1d0)),DIMENSION(NG)::lx,ly,cb,cc,cd
    REAL(KIND(1d0)),DIMENSION(NH)::hx,hy,hb,hc,hd
    INTEGER::gset=-1,hset=-1,set,i
    REAL(KIND(1d0))::dl
    SAVE lx,ly,cb,cc,cd,gset,hx,hy,hb,hc,hd,hset
    IF(x.LT.XGMIN.OR.x.GT.XHMAX.OR..NOT.GENERATE_PROTON_FLUX_GRID.OR.&
         .NOT.(USE_EXACT_ELASTIC4PHOTON.OR.alpha_running4photon))THEN
       PhotonFlux_proton_nob=direct(x)
       RETURN
    ENDIF
    set=proton_FF_correction+MERGE(10,0,USE_EXACT_ELASTIC4PHOTON)+MERGE(100,0,alpha_running4photon)
    dl=(logit(XGMAX)-logit(XGMIN))/(NG-1)
    IF(x.GT.XGMAX)THEN
       IF(hset.NE.set)THEN
          DO i=1,NH
             hx(i)=logit(XGMAX)+(i-1-NHO)*0.5d0*dl
             hy(i)=DLOG(xoft(hx(i))*direct(xoft(hx(i))))
          ENDDO
          CALL spline(hx,hy,hb,hc,hd,NH)
          hset=set
       ENDIF
       PhotonFlux_proton_nob=DEXP(ispline(logit(x),hx,hy,hb,hc,hd,NH))/x
       RETURN
    ENDIF
    IF(gset.NE.set)THEN
       DO i=1,NG
          lx(i)=logit(XGMIN)+(i-1)*(logit(XGMAX)-logit(XGMIN))/(NG-1)
          ly(i)=DLOG(xoft(lx(i))*direct(xoft(lx(i))))
       ENDDO
       CALL spline(lx,ly,cb,cc,cd,NG)
       gset=set
    ENDIF
    PhotonFlux_proton_nob=DEXP(ispline(logit(x),lx,ly,cb,cc,cd,NG))/x
    RETURN
  CONTAINS
    FUNCTION direct(xx)
      REAL(KIND(1d0))::direct
      REAL(KIND(1d0)),INTENT(IN)::xx
      IF(USE_EXACT_ELASTIC4PHOTON)THEN
         direct=PhotonFlux_proton_ChFF_exact_nob(xx,proton_FF_correction)
      ELSE
         direct=PhotonFlux_proton_ChFF_nob(xx,gamma,proton_FF_correction)
      ENDIF
    END FUNCTION direct
    FUNCTION logit(xx)
      REAL(KIND(1d0))::logit
      REAL(KIND(1d0)),INTENT(IN)::xx
      logit=DLOG(xx/(1d0-xx))
    END FUNCTION logit
    FUNCTION xoft(t) ! inverse of logit
      REAL(KIND(1d0))::xoft
      REAL(KIND(1d0)),INTENT(IN)::t
      xoft=1d0/(1d0+DEXP(-t))
    END FUNCTION xoft
  END FUNCTION PhotonFlux_proton_nob


  FUNCTION e_pn(ib)
    ! energy per nucleon (GeV) of beam ib from run90.inc (ebeamMG5: whole beam particle), as ebeam_PN of v2.0
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::e_pn
    INTEGER,INTENT(IN)::ib
    e_pn=ebeamMG5(ib)/(nb_proton(ib)+nb_neutron(ib))
  END FUNCTION e_pn

  FUNCTION photon_mu2(x1,x2)
    ! the scale of the inelastic photons in pp: mu = W_gammagamma, mu2 = x1 x2 s (GeV^2)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0))::photon_mu2
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    photon_mu2=x1*x2*(4d0*e_pn(1)*e_pn(2))
  END FUNCTION photon_mu2

  SUBROUTINE set_photon_scale(x1,x2)
    ! the scale mu2_photon of the photon-photon flux at (x1,x2) and, with alpha_running4photon, the normalisation
    ! alpha(0)/alpha(mu2) of the elastic photons (ael_photon), as the inelastic photon (LUXqed, MSbar at mu): every
    ! photon then carries alpha(Q2)**2/alpha(mu2) and the hard process uses alpha(mu2); fixed alpha: 1, alpha(0)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    mu2_photon=photon_mu2(x1,x2)
    ael_photon=alpha_norm_mu(mu2_photon)
  END SUBROUTINE set_photon_scale

  FUNCTION final_state_nob(el1,el2,x1,x2,np)
    ! d2N/dx1dx2 with P_NOHAD=1 for photon_final_state, from the elastic dN/dx el1, el2 and the inelastic ones at
    ! mu2_photon (for 2 only inel-inel: 1-P=0); np = 2: pp, np = 1: pA with the proton as beam 1 (no inelastic
    ! photon of the nucleus)
    IMPLICIT NONE
    REAL(KIND(1d0))::final_state_nob
    REAL(KIND(1d0)),INTENT(IN)::el1,el2,x1,x2
    INTEGER,INTENT(IN)::np
    REAL(KIND(1d0))::in1,in2
    IF(photon_final_state.EQ.0)THEN
       final_state_nob=el1*el2
       RETURN
    ENDIF
    in1=PhotonFlux_proton_inelastic_nob(x1,mu2_photon)
    in2=0d0
    IF(np.EQ.2)in2=PhotonFlux_proton_inelastic_nob(x2,mu2_photon)
    SELECT CASE(photon_final_state)
    CASE(1)
       final_state_nob=el1*in2+in1*el2
    CASE(2)
       final_state_nob=in1*in2
    CASE DEFAULT
       final_state_nob=(el1+in1)*(el2+in2)
    END SELECT
  END FUNCTION final_state_nob

  FUNCTION PhotonNumberDensity_inelastic_beam(b,Ega,gamma,ibeam)
    ! the inelastic proton photon density (Z=1, alpha=1) at mu2_photon; the last value of each beam (ibeam) is kept
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_inelastic_beam
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    INTEGER,INTENT(IN)::ibeam
    REAL(KIND(1d0)),DIMENSION(5,2)::last=-1d0 ! b,Ega,gamma,mu2,value
    SAVE last
    IF(ALL(last(1:4,ibeam).EQ.(/b,Ega,gamma,mu2_photon/)))THEN
       PhotonNumberDensity_inelastic_beam=last(5,ibeam)
       RETURN
    ENDIF
    PhotonNumberDensity_inelastic_beam=PhotonNumberDensity_inelastic_proton(b,Ega,gamma,mu2_photon)
    last(:,ibeam)=(/b,Ega,gamma,mu2_photon,PhotonNumberDensity_inelastic_beam/)
  END FUNCTION PhotonNumberDensity_inelastic_beam


  LOGICAL FUNCTION qspace_ok()
    ! (M11) the momentum-space (1-P) parts apply: pp, exact elastic photon with ChFF, pp_final_state_qspace
    IMPLICIT NONE
    qspace_ok=pp_final_state_qspace.AND.USE_CHARGEFORMFACTOR4PHOTON.AND.USE_EXACT_ELASTIC4PHOTON
  END FUNCTION qspace_ok

  SUBROUTINE qspace_setup()
    ! (M11) q nodes, weights with (1-P)^ (slope b0 of this energy), Chebyshev nodes and the interpolation matrix
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::t6=(/-0.932469514203152027812301554494d0,-0.661209386466264513661399595020d0,&
         -0.238619186083196908630501721681d0,0.238619186083196908630501721681d0,0.661209386466264513661399595020d0,&
         0.932469514203152027812301554494d0/)
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::w6=(/0.171324492379170345040296142173d0,0.360761573048138607569833513838d0,&
         0.467913934572691047389870343990d0,0.467913934572691047389870343990d0,0.360761573048138607569833513838d0,&
         0.171324492379170345040296142173d0/)
    REAL(KIND(1d0))::b0,s,h,t
    INTEGER::i,j,k,np
    s=4d0*e_pn(1)*e_pn(2)
    b0=9.81D0+0.211D0*DLOG(s)+0.0185D0*DLOG(s)**2 ! as PhotonPhotonFlux_pp_eval
    IF(b0.EQ.qs_b0)RETURN
    qs_b0=b0
    np=QS_NQ/6; h=(QS_HI-QS_LO)/np
    DO i=1,np
       DO j=1,6
          k=6*(i-1)+j
          t=QS_LO+(i-0.5d0)*h+0.5d0*h*t6(j)
          qs_q(k)=DEXP(t)
          qs_w(k)=0.5d0*h*w6(j)*qs_q(k)**2*(4d0*pi*b0*DEXP(-b0*qs_q(k)**2/2d0)-pi*b0*DEXP(-b0*qs_q(k)**2/4d0))/(2d0*pi)
       ENDDO
    ENDDO
    CALL qs_cheb_nodes(QS_NC1,QS_NC2,qs_qc)
    DO k=1,QS_NQ
       CALL qs_ci_row(qs_q(k),QS_NC1,QS_NC2,qs_ci(k,:))
    ENDDO
    ! new weights: the stored transforms stay valid (they do not depend on b0)
  END SUBROUTINE qspace_setup

  SUBROUTINE qs_cheb_nodes(n1,n2,qc)
    ! (M11) the Chebyshev nodes of the inelastic transforms: n1 in ln q on [1e-6, QS_QM], n2 in q on [QS_QM, 1.6] GeV
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n1,n2
    REAL(KIND(1d0)),DIMENSION(n1+n2),INTENT(OUT)::qc
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    INTEGER::i
    DO i=1,n1
       qc(i)=DEXP(0.5d0*(QS_LO+DLOG(QS_QM))-0.5d0*(DLOG(QS_QM)-QS_LO)*DCOS(pi*(i-1)/(n1-1)))
    ENDDO
    DO i=1,n2
       qc(n1+i)=0.5d0*(QS_QM+DEXP(QS_HI))-0.5d0*(DEXP(QS_HI)-QS_QM)*DCOS(pi*(i-1)/(n2-1))
    ENDDO
  END SUBROUTINE qs_cheb_nodes

  SUBROUTINE qs_ci_row(q,n1,n2,row)
    ! (M11) the row at q of the interpolation of the inelastic transforms from their Chebyshev nodes (n1 in ln q on
    ! [1e-6, QS_QM], n2 in q on [QS_QM, 1.6] GeV: qs_cheb_nodes)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::q
    INTEGER,INTENT(IN)::n1,n2
    REAL(KIND(1d0)),DIMENSION(n1+n2),INTENT(OUT)::row
    REAL(KIND(1d0))::u
    row=0d0
    IF(q.LE.QS_QM)THEN
       u=(2d0*DLOG(q)-QS_LO-DLOG(QS_QM))/(DLOG(QS_QM)-QS_LO)
       CALL cheb_row(n1,u,row(1:n1))
    ELSE
       u=(2d0*q-QS_QM-DEXP(QS_HI))/(DEXP(QS_HI)-QS_QM)
       CALL cheb_row(n2,u,row(n1+1:n1+n2))
    ENDIF
  CONTAINS
    SUBROUTINE cheb_row(n,uu,w)
      ! barycentric weights of the n Chebyshev-Lobatto points -cos(pi (i-1)/(n-1)) at uu in [-1,1]
      INTEGER,INTENT(IN)::n
      REAL(KIND(1d0)),INTENT(IN)::uu
      REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::w
      REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
      INTEGER::ii
      REAL(KIND(1d0))::xi
      DO ii=1,n
         xi=-DCOS(pi*(ii-1)/(n-1))
         IF(ABS(uu-xi).LT.1d-15)THEN
            w=0d0; w(ii)=1d0; RETURN
         ENDIF
         w(ii)=(-1d0)**(ii-1)/(uu-xi)
         IF(ii.EQ.1.OR.ii.EQ.n)w(ii)=0.5d0*w(ii)
      ENDDO
      w=w/SUM(w)
    END SUBROUTINE cheb_row
  END SUBROUTINE qs_ci_row

  SUBROUTINE qspace_keys()
    ! (M11) the stores are emptied when their settings change
    IMPLICIT NONE
    include 'run90.inc'
    INTEGER::k
    k=proton_FF_correction+MERGE(10,0,alpha_running4photon)
    IF(k.NE.qs_ekey)THEN
       qs_ekey=k; qs_ne=0
    ENDIF
    k=inel_rhohat_key()
    IF(k.NE.qs_ikey)THEN
       qs_ikey=k; qs_ni=0
    ENDIF
  END SUBROUTINE qspace_keys

  FUNCTION qs_near(a,c)
    IMPLICIT NONE
    LOGICAL::qs_near
    REAL(KIND(1d0)),INTENT(IN)::a,c
    qs_near=ABS(a-c).LE.1d-13*MAX(ABS(a),ABS(c))
  END FUNCTION qs_near

  INTEGER FUNCTION qs_el(x)
    ! (M11) the store index of the elastic transform at x (computed if missing)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::i
    DO i=qs_ne,1,-1
       IF(qs_near(qs_ex(i),x))THEN
          qs_el=i; RETURN
       ENDIF
    ENDDO
    IF(qs_ne.GE.QS_MAXE)qs_ne=0 ! full: start again
    qs_ne=qs_ne+1
    qs_ex(qs_ne)=x
    CALL elastic_rhohat_cached(x,proton_FF_correction,QS_NQ,qs_q,qs_er(:,qs_ne))
    qs_el=qs_ne
  END FUNCTION qs_el

  INTEGER FUNCTION qs_in_find(x,mu2)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,mu2
    INTEGER::i
    qs_in_find=0
    DO i=qs_ni,1,-1
       IF(qs_near(qs_ix(i),x).AND.qs_near(qs_imu2(i),mu2))THEN
          qs_in_find=i; RETURN
       ENDIF
    ENDDO
  END FUNCTION qs_in_find

  SUBROUTINE qs_in_add(x,nmu,mu2v,busy)
    ! (M11) the inelastic transforms at x and the scales mu2v that are not in the store
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nmu
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    LOGICAL,INTENT(OUT),OPTIONAL::busy ! (M16) the row is being computed by another job (disk cache lock)
    REAL(KIND(1d0)),DIMENSION(nmu)::mm
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::R
    INTEGER::i,j,n
    LOGICAL::dup
    IF(PRESENT(busy))busy=.FALSE.
    n=0
    DO i=1,nmu
       IF(qs_in_find(x,mu2v(i)).GT.0)CYCLE
       dup=.FALSE.
       DO j=1,n
          IF(qs_near(mm(j),mu2v(i)))dup=.TRUE.
       ENDDO
       IF(dup)CYCLE
       n=n+1; mm(n)=mu2v(i)
    ENDDO
    IF(n.EQ.0)RETURN
    ALLOCATE(R(QS_NC,n))
    CALL inel_rhohat_cached(x,QS_NC,qs_qc,n,mm(1:n),R,busy)
    IF(PRESENT(busy))THEN
       IF(busy)RETURN
    ENDIF
    IF(qs_ni+n.GT.QS_MAXI)qs_ni=0 ! full: start again
    DO i=1,n
       qs_ni=qs_ni+1
       qs_ix(qs_ni)=x; qs_imu2(qs_ni)=mm(i); qs_ir(:,qs_ni)=R(:,i)
    ENDDO
    DEALLOCATE(R)
  END SUBROUTINE qs_in_add

  SUBROUTINE pp_qspace_prepare(n1,x1s,n2,x2s)
    ! (M11) the elastic and inelastic transforms of all x of the grid and all scales of the pairs
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n1,n2
    REAL(KIND(1d0)),DIMENSION(n1),INTENT(IN)::x1s
    REAL(KIND(1d0)),DIMENSION(n2),INTENT(IN)::x2s
    REAL(KIND(1d0)),DIMENSION(MAX(n1,n2))::mu2v
    INTEGER::i,j,k,ib
    CALL qspace_setup()
    CALL qspace_keys()
    DO ib=1,2
       DO i=1,MERGE(n1,n2,ib.EQ.1)
          IF(ib.EQ.1)THEN
             IF(x1s(i).LE.0d0.OR.x1s(i).GE.1d0)CYCLE
             k=0
             DO j=1,n2
                IF(x2s(j).LE.0d0.OR.x2s(j).GE.1d0)CYCLE
                k=k+1; mu2v(k)=photon_mu2(x1s(i),x2s(j))
             ENDDO
             IF(k.GT.0)CALL qs_in_add(x1s(i),k,mu2v(1:k))
             k=qs_el(x1s(i))
          ELSE
             IF(x2s(i).LE.0d0.OR.x2s(i).GE.1d0)CYCLE
             k=0
             DO j=1,n1
                IF(x1s(j).LE.0d0.OR.x1s(j).GE.1d0)CYCLE
                k=k+1; mu2v(k)=photon_mu2(x1s(j),x2s(i))
             ENDDO
             IF(k.GT.0)CALL qs_in_add(x2s(i),k,mu2v(1:k))
             k=qs_el(x2s(i))
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE pp_qspace_prepare


  SUBROUTINE pp_qspace_prepare_grid(n1,x1s,m2,y2s,n2,x2s,m1,y1s,mu2lo)
    ! (M11) for each x1 of x1s the inelastic transforms at the scales with every y2 of y2s (its grid row), for each
    ! x2 of x2s those with every y1 of y1s, at scales mu2 >= mu2lo; the elastic transforms of all of them
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n1,m2,n2,m1
    REAL(KIND(1d0)),DIMENSION(n1),INTENT(IN)::x1s
    REAL(KIND(1d0)),DIMENSION(m2),INTENT(IN)::y2s
    REAL(KIND(1d0)),DIMENSION(n2),INTENT(IN)::x2s
    REAL(KIND(1d0)),DIMENSION(m1),INTENT(IN)::y1s
    REAL(KIND(1d0)),INTENT(IN)::mu2lo
    REAL(KIND(1d0)),DIMENSION(MAX(m1,m2))::mu2v
    REAL(KIND(1d0))::m2x
    INTEGER::i,j,k,ipass
    LOGICAL,DIMENSION(n1)::todo1
    LOGICAL,DIMENSION(n2)::todo2
    LOGICAL::busy
    CALL qspace_setup()
    CALL qspace_keys()
    ! (M16) pass 1: the rows that no other job is computing (disk cache locks); pass 2: the others (waiting)
    todo1=.TRUE.; todo2=.TRUE.
    DO ipass=1,2
    DO i=1,n1
       IF(.NOT.todo1(i))CYCLE
       IF(x1s(i).LE.0d0.OR.x1s(i).GE.1d0)CYCLE
       k=0
       DO j=1,m2
          IF(y2s(j).LE.0d0.OR.y2s(j).GE.1d0)CYCLE
          m2x=photon_mu2(x1s(i),y2s(j))
          IF(m2x.LT.mu2lo)CYCLE
          k=k+1; mu2v(k)=m2x
       ENDDO
       busy=.FALSE.
       IF(k.GT.0.AND.ipass.EQ.1)CALL qs_in_add(x1s(i),k,mu2v(1:k),busy)
       IF(k.GT.0.AND.ipass.EQ.2)CALL qs_in_add(x1s(i),k,mu2v(1:k))
       todo1(i)=busy
       k=qs_el(x1s(i))
    ENDDO
    DO i=1,n2
       IF(.NOT.todo2(i))CYCLE
       IF(x2s(i).LE.0d0.OR.x2s(i).GE.1d0)CYCLE
       k=0
       DO j=1,m1
          IF(y1s(j).LE.0d0.OR.y1s(j).GE.1d0)CYCLE
          m2x=photon_mu2(y1s(j),x2s(i))
          IF(m2x.LT.mu2lo)CYCLE
          k=k+1; mu2v(k)=m2x
       ENDDO
       busy=.FALSE.
       IF(k.GT.0.AND.ipass.EQ.1)CALL qs_in_add(x2s(i),k,mu2v(1:k),busy)
       IF(k.GT.0.AND.ipass.EQ.2)CALL qs_in_add(x2s(i),k,mu2v(1:k))
       todo2(i)=busy
       k=qs_el(x2s(i))
    ENDDO
    IF(.NOT.(ANY(todo1).OR.ANY(todo2)))EXIT
    ENDDO
  END SUBROUTINE pp_qspace_prepare_grid

  SUBROUTINE pp_qspace_L(x1,x2,Lee,Lei)
    ! (M11) at the scale mu2_photon: Lee = (1/2 pi) int q dq (1-P)^ el1^ el2^, Lei = the same with el1^ in2^ +
    ! in1^ el2^ (elastic at alpha(0), inelastic as InelasticFlux)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    REAL(KIND(1d0)),INTENT(OUT)::Lee,Lei
    REAL(KIND(1d0)),DIMENSION(QS_NQ)::e1,e2,i1,i2
    INTEGER::k
    CALL qspace_setup()
    CALL qspace_keys()
    e1=qs_er(:,qs_el(x1))
    e2=qs_er(:,qs_el(x2))
    k=qs_in_find(x1,mu2_photon)
    IF(k.EQ.0)THEN
       CALL qs_in_add(x1,1,(/mu2_photon/)); k=qs_in_find(x1,mu2_photon)
    ENDIF
    i1=MATMUL(qs_ci,qs_ir(:,k))
    k=qs_in_find(x2,mu2_photon)
    IF(k.EQ.0)THEN
       CALL qs_in_add(x2,1,(/mu2_photon/)); k=qs_in_find(x2,mu2_photon)
    ENDIF
    i2=MATMUL(qs_ci,qs_ir(:,k))
    Lee=SUM(qs_w*e1*e2)
    Lei=SUM(qs_w*(e1*i2+i1*e2))
  END SUBROUTINE pp_qspace_L

  SUBROUTINE pA_qspace_setup(RA,wA,aA,A,sigNN)
    ! (M15) the q nodes and weights of pA for the nucleus (Woods-Saxon R_A, w, a in fm, mass number A) and sigma_NN
    ! (fm^2): (1-P_pA)^(q) = 2 pi int b db J0(q b) (1 - P_pA(b)), P_pA = exp(-sigma_NN A T^_A(b)) as
    ! PNOHAD_pA_WoodsSaxon, with T^_A by Gauss-Legendre in z (normalised to 1) instead of its interpolation grid
    ! (TA_MC_Glauber with use_MC_Glauber), on Gauss-Legendre 8 panels in b of width <= min(a/2, pi/(2 q_max)) up to
    ! R_A + 40 a (1-P < 1e-16 beyond); the stores are emptied when the nucleus or sigma_NN change
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::RA,wA,aA,A,sigNN
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    REAL(KIND(1d0)),DIMENSION(6)::key
    REAL(KIND(1d0)),DIMENSION(1)::rc
    REAL(KIND(1d0)),DIMENSION(QA_NG)::tg,wq
    REAL(KIND(1d0)),DIMENSION(8)::xg,wg
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::bn,wb,omp
    REAL(KIND(1d0))::R,q1,h1,h2,t,bmax,hb,zmax,hz,z,tn,chi,Pih
    INTEGER::n1,n2,i,j,k,nb,nz,iz,g
    key=(/RA,wA,aA,A,sigNN,MERGE(1d0,0d0,use_MC_Glauber)/)
    IF(ALL(key.EQ.qa_key))RETURN
    qa_key=key
    R=RA/GeVm12fm
    q1=QA_Q1R/R
    n1=CEILING((DLOG(q1)-QS_LO)/1.1d0); h1=(DLOG(q1)-QS_LO)/n1
    n2=CEILING((DEXP(QS_HI)-q1)/(QA_DQR/R)); h2=(DEXP(QS_HI)-q1)/n2
    qa_nq=QA_NG*(n1+n2)
    IF(ALLOCATED(qa_q))DEALLOCATE(qa_q,qa_w,qa_ci,qa_er,qa_nr,qa_ex,qa_nx)
    ALLOCATE(qa_q(qa_nq+1),qa_w(qa_nq+1),qa_ci(qa_nq+1,QA_NC),qa_er(qa_nq+1,QA_MAXE),qa_nr(qa_nq+1,QA_MAXN),&
         qa_ex(QA_MAXE),qa_nx(QA_MAXN))
    IF(.NOT.ALLOCATED(qa_ir))THEN
       ALLOCATE(qa_ix(QA_MAXI),qa_imu2(QA_MAXI),qa_ir(QA_NC,QA_MAXI))
       CALL qs_cheb_nodes(QA_NC1,QA_NC2,qa_qc)
    ENDIF
    qa_q(qa_nq+1)=0d0; qa_w(qa_nq+1)=0d0; qa_ci(qa_nq+1,:)=0d0
    CALL ffe_gauleg(QA_NG,tg,wq)
    DO i=1,n1
       DO j=1,QA_NG
          k=QA_NG*(i-1)+j
          t=QS_LO+(i-0.5d0)*h1+0.5d0*h1*tg(j)
          qa_q(k)=DEXP(t); qa_w(k)=0.5d0*h1*wq(j)*qa_q(k)**2 ! q dq = q^2 dln q
       ENDDO
    ENDDO
    DO i=1,n2
       DO j=1,QA_NG
          k=QA_NG*(n1+i-1)+j
          qa_q(k)=q1+(i-0.5d0)*h2+0.5d0*h2*tg(j); qa_w(k)=0.5d0*h2*wq(j)*qa_q(k)
       ENDDO
    ENDDO
    ! 1 - P_pA on the b nodes (GeV^-1)
    CALL ffe_gauleg(8,xg,wg)
    bmax=(RA+40d0*aA)/GeVm12fm
    hb=MIN(0.5d0*aA/GeVm12fm,0.5d0*pi/DEXP(QS_HI))
    nb=CEILING(bmax/hb); hb=bmax/nb
    ALLOCATE(bn(8*nb),wb(8*nb),omp(8*nb))
    DO i=1,nb
       DO g=1,8
          k=8*(i-1)+g
          bn(k)=(i-0.5d0)*hb+0.5d0*hb*xg(g); wb(k)=0.5d0*hb*wg(g)
       ENDDO
    ENDDO
    IF(use_MC_Glauber)THEN
       DO k=1,8*nb
          chi=TA_MC_Glauber(bn(k)*GeVm12fm,0d0,A)*sigNN ! as PNOHAD_pA_WoodsSaxon
          omp(k)=MERGE(1d0,-expm1_(-chi),chi.GT.500d0)
       ENDDO
    ELSE
       ! T_A(b) = int dz rho(sqrt(b^2+z^2)) (in units of the normalisation of rho_WoodsSaxon), Gauss-Legendre
       ! 8 on panels of width a/2 in z; normalised to int d^2b T^_A = 1 on the b nodes
       DO k=1,8*nb
          zmax=DSQRT(MAX((RA+40d0*aA)**2-(bn(k)*GeVm12fm)**2,0d0))
          nz=MAX(1,CEILING(zmax/(0.5d0*aA))); hz=zmax/nz
          tn=0d0
          DO iz=1,nz
             DO g=1,8
                z=(iz-0.5d0)*hz+0.5d0*hz*xg(g)
                tn=tn+hz*wg(g)*rho_WoodsSaxon(DSQRT((bn(k)*GeVm12fm)**2+z**2),RA,wA,aA,A) ! both signs of z
             ENDDO
          ENDDO
          omp(k)=tn
       ENDDO
       tn=2d0*pi*SUM(wb*bn*omp)*GeVm12fm**2 ! int d^2b T_A, b in fm
       DO k=1,8*nb
          chi=sigNN*A*omp(k)/tn
          omp(k)=MERGE(1d0,-expm1_(-chi),chi.GT.500d0)
       ENDDO
    ENDIF
    DO k=1,qa_nq
       Pih=2d0*pi*SUM(wb*bn*BESSEL_J0(qa_q(k)*bn)*omp)
       qa_w(k)=qa_w(k)*Pih/(2d0*pi)
       CALL qs_ci_row(qa_q(k),QA_NC1,QA_NC2,qa_ci(k,:))
    ENDDO
    qa_ne=0; qa_nn=0
    ! the collinear table of the nucleus (Parseval: one q = 0 transform per node)
    qa_cn=NINT(QA_CPD*DLOG10(QA_CXHI/QA_CXLO))+1; qa_ch=DLOG(QA_CXHI/QA_CXLO)/(qa_cn-1)
    IF(ALLOCATED(qa_ct))DEALLOCATE(qa_ct)
    ALLOCATE(qa_ct(qa_cn))
    ! (M16) through the disk cache of inel_cache_dir (one file per nucleus and table)
    IF(.NOT.vec_cache_get('nuct',(/1d0,R,wA,aA,QA_CXLO,QA_CXHI,DBLE(QA_CPD)/),qa_ct))THEN
       DO i=1,qa_cn
          CALL nucleus_rhohat(QA_CXLO*DEXP((i-1)*qa_ch),R,wA,aA/GeVm12fm,1,(/0d0/),rc)
          qa_ct(i)=DLOG(rc(1))
       ENDDO
       CALL vec_cache_put('nuct',(/1d0,R,wA,aA,QA_CXLO,QA_CXHI,DBLE(QA_CPD)/),qa_ct)
    ENDIF
  CONTAINS
    FUNCTION expm1_(y) ! exp(y)-1 without cancellation at small |y|
      REAL(KIND(1d0))::expm1_
      REAL(KIND(1d0)),INTENT(IN)::y
      IF(ABS(y).LT.1d-5)THEN
         expm1_=y*(1d0+y*(0.5d0+y/6d0))
      ELSE
         expm1_=DEXP(y)-1d0
      ENDIF
    END FUNCTION expm1_
  END SUBROUTINE pA_qspace_setup

  SUBROUTINE pA_qspace_setup_common()
    ! (M15) pA_qspace_setup with the nucleus of PhotonPhotonFlux_pA_WoodsSaxon_eval (its COMMON, set at its first call)
    IMPLICIT NONE
    REAL(KIND(1d0))::gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,sigNN_inel_common,RA_common,&
         A_common,aaVal_common,wVal_common
    COMMON/PhotonPhoton_pA_WS/gamma1_common,gamma2_common,E1_common,E2_common,x1_common,x2_common,&
         sigNN_inel_common,RA_common,A_common,aaVal_common,wVal_common
    CALL pA_qspace_setup(RA_common,wVal_common,aaVal_common,A_common,sigNN_inel_common)
  END SUBROUTINE pA_qspace_setup_common

  SUBROUTINE pA_qspace_keys()
    ! (M15) the stores of the proton transforms are emptied when their settings change
    IMPLICIT NONE
    include 'run90.inc'
    INTEGER::k
    k=proton_FF_correction+MERGE(10,0,alpha_running4photon)
    IF(k.NE.qa_ekey)THEN
       qa_ekey=k; qa_ne=0
    ENDIF
    k=inel_rhohat_key()
    IF(k.NE.qa_ikey)THEN
       qa_ikey=k; qa_ni=0
    ENDIF
  END SUBROUTINE pA_qspace_keys

  INTEGER FUNCTION qa_in_find(x,mu2)
    ! (M15) the store index of the inelastic transform of pA at x, mu2 (0: missing)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,mu2
    INTEGER::i
    qa_in_find=0
    DO i=qa_ni,1,-1
       IF(qs_near(qa_ix(i),x).AND.qs_near(qa_imu2(i),mu2))THEN
          qa_in_find=i; RETURN
       ENDIF
    ENDDO
  END FUNCTION qa_in_find

  SUBROUTINE qa_in_add(x,nmu,mu2v,busy)
    ! (M15) the inelastic transforms of pA (Chebyshev nodes qa_qc) at x and the scales mu2v not in the store
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nmu
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    LOGICAL,INTENT(OUT),OPTIONAL::busy ! (M16) the row is being computed by another job (disk cache lock)
    REAL(KIND(1d0)),DIMENSION(nmu)::mm
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::R
    INTEGER::i,j,n
    LOGICAL::dup
    IF(PRESENT(busy))busy=.FALSE.
    n=0
    DO i=1,nmu
       IF(qa_in_find(x,mu2v(i)).GT.0)CYCLE
       dup=.FALSE.
       DO j=1,n
          IF(qs_near(mm(j),mu2v(i)))dup=.TRUE.
       ENDDO
       IF(dup)CYCLE
       n=n+1; mm(n)=mu2v(i)
    ENDDO
    IF(n.EQ.0)RETURN
    ALLOCATE(R(QA_NC,n))
    CALL inel_rhohat_cached(x,QA_NC,qa_qc,n,mm(1:n),R,busy)
    IF(PRESENT(busy))THEN
       IF(busy)RETURN
    ENDIF
    IF(qa_ni+n.GT.QA_MAXI)qa_ni=0 ! full: start again
    DO i=1,n
       qa_ni=qa_ni+1
       qa_ix(qa_ni)=x; qa_imu2(qa_ni)=mm(i); qa_ir(:,qa_ni)=R(:,i)
    ENDDO
    DEALLOCATE(R)
  END SUBROUTINE qa_in_add

  INTEGER FUNCTION qa_el(x)
    ! (M15) the store index of the proton elastic transform at x on the pA nodes (computed if missing)
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::i
    DO i=qa_ne,1,-1
       IF(qs_near(qa_ex(i),x))THEN
          qa_el=i; RETURN
       ENDIF
    ENDDO
    IF(qa_ne.GE.QA_MAXE)qa_ne=0 ! full: start again
    qa_ne=qa_ne+1
    qa_ex(qa_ne)=x
    CALL elastic_rhohat_cached(x,proton_FF_correction,qa_nq+1,qa_q,qa_er(:,qa_ne))
    qa_el=qa_ne
  END FUNCTION qa_el

  INTEGER FUNCTION qa_nuc(x)
    ! (M15) the store index of the nucleus transform at x on the pA nodes (computed if missing)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::i
    DO i=qa_nn,1,-1
       IF(qs_near(qa_nx(i),x))THEN
          qa_nuc=i; RETURN
       ENDIF
    ENDDO
    IF(qa_nn.GE.QA_MAXN)qa_nn=0 ! full: start again
    qa_nn=qa_nn+1
    qa_nx(qa_nn)=x
    CALL nucleus_rhohat_cached(x,qa_key(1)/GeVm12fm,qa_key(2),qa_key(3)/GeVm12fm,qa_nq+1,qa_q,qa_nr(:,qa_nn))
    qa_nuc=qa_nn
  END FUNCTION qa_nuc

  SUBROUTINE pA_qspace_L(xp,xA,Lee,Lie)
    ! (M15) at the scale mu2_photon, the proton at xp, the nucleus at xA: Lee = (1/2 pi) int q dq (1-P_pA)^ el_p^
    ! el_A^, Lie = the same with in_p^ (elastic photons at alpha(0) and unit charge, inelastic as InelasticFlux);
    ! after pA_qspace_setup
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::xp,xA
    REAL(KIND(1d0)),INTENT(OUT)::Lee,Lie
    REAL(KIND(1d0)),DIMENSION(qa_nq)::eA
    INTEGER::k
    CALL pA_qspace_keys()
    eA=qa_w(1:qa_nq)*qa_nr(1:qa_nq,qa_nuc(xA))
    Lee=SUM(eA*qa_er(1:qa_nq,qa_el(xp)))
    k=qa_in_find(xp,mu2_photon)
    IF(k.EQ.0)THEN
       CALL qa_in_add(xp,1,(/mu2_photon/)); k=qa_in_find(xp,mu2_photon)
    ENDIF
    Lie=SUM(eA*MATMUL(qa_ci(1:qa_nq,:),qa_ir(:,k)))
  END SUBROUTINE pA_qspace_L

  FUNCTION pA_nucleus_coll(xA)
    ! (M15) int d^2b N_A(b) of the nucleus at xA (unit charge, alpha = 1): its transform at q = 0, from the table of
    ! pA_qspace_setup; Z^2 alpha/xA times this is the collinear flux of PhotonFlux_nucleus_ChFF_nob, whose table is off
    ! by up to 6e-3 at xA < 1e-4 (drivers/m15_nuc_coll); after pA_qspace_setup
    IMPLICIT NONE
    REAL(KIND(1d0))::pA_nucleus_coll
    REAL(KIND(1d0)),INTENT(IN)::xA
    REAL(KIND(1d0)),DIMENSION(8)::wt
    REAL(KIND(1d0)),DIMENSION(1)::r
    INTEGER::i0
    IF(xA.GE.QA_CXLO.AND.xA.LE.QA_CXHI)THEN
       ! the table (the transforms at all q of a new x cost 0.75 s; the luminosity integral asks for many x)
       CALL ffe_lagu(DLOG(xA/QA_CXLO)/qa_ch,qa_cn,8,i0,wt)
       pA_nucleus_coll=DEXP(SUM(wt*qa_ct(i0+1:i0+8)))
    ELSE
       CALL nucleus_rhohat(xA,qa_key(1)/GeVm12fm,qa_key(2),qa_key(3)/GeVm12fm,1,(/0d0/),r)
       pA_nucleus_coll=r(1)
    ENDIF
  END FUNCTION pA_nucleus_coll

  SUBROUTINE pA_qspace_prepare_grid(np,xps,my,ys,nA,xAs,mu2lo)
    ! (M15) for each proton x of xps the inelastic transforms at the scales with every nucleus x of ys (its grid row)
    ! at scales mu2 >= mu2lo, and its elastic transform; the transforms of the nucleus at xAs; after pA_qspace_setup
    IMPLICIT NONE
    INTEGER,INTENT(IN)::np,my,nA
    REAL(KIND(1d0)),DIMENSION(np),INTENT(IN)::xps
    REAL(KIND(1d0)),DIMENSION(my),INTENT(IN)::ys
    REAL(KIND(1d0)),DIMENSION(nA),INTENT(IN)::xAs
    REAL(KIND(1d0)),INTENT(IN)::mu2lo
    REAL(KIND(1d0)),DIMENSION(my)::mu2v
    REAL(KIND(1d0))::m2x
    INTEGER::i,j,k,ipass
    LOGICAL,DIMENSION(np)::todo
    LOGICAL::busy
    CALL pA_qspace_keys()
    todo=.TRUE. ! (M16) pass 1: the rows no other job is computing (disk cache locks); pass 2: the others
    DO ipass=1,2
    DO i=1,np
       IF(.NOT.todo(i))CYCLE
       IF(xps(i).LE.0d0.OR.xps(i).GE.1d0)CYCLE
       k=0
       DO j=1,my
          IF(ys(j).LE.0d0.OR.ys(j).GE.1d0)CYCLE
          m2x=photon_mu2(xps(i),ys(j))
          IF(m2x.LT.mu2lo)CYCLE
          k=k+1; mu2v(k)=m2x
       ENDDO
       busy=.FALSE.
       IF(k.GT.0.AND.ipass.EQ.1)CALL qa_in_add(xps(i),k,mu2v(1:k),busy)
       IF(k.GT.0.AND.ipass.EQ.2)CALL qa_in_add(xps(i),k,mu2v(1:k))
       todo(i)=busy
       k=qa_el(xps(i))
    ENDDO
    IF(.NOT.ANY(todo))EXIT
    ENDDO
    DO i=1,nA
       IF(xAs(i).LE.0d0.OR.xAs(i).GE.1d0)CYCLE
       k=qa_nuc(xAs(i))
    ENDDO
  END SUBROUTINE pA_qspace_prepare_grid

  SUBROUTINE pp_inelastic_prepare(n1,x1s,n2,x2s,ONLYBEAM1)
    ! ONLYBEAM1 (pA, the proton as beam 1): only the rows of beam 1
    ! fills the rows of the inelastic densities of beam 1 at x1s and of beam 2 at x2s at the
    ! b nodes of PhotonPhotonFlux_pp_eval (midpoints of the panels in ln(x b mproton), pp_inelastic_range, up to bupper)
    ! and the scales photon_mu2 of all pairs, in one call (inel_fast_prepare_rows: parallel with OpenMP); the values are
    ! those computed on demand otherwise
    IMPLICIT NONE
    include 'run90.inc'
    INTEGER,INTENT(IN)::n1,n2
    REAL(KIND(1d0)),DIMENSION(n1),INTENT(IN)::x1s
    REAL(KIND(1d0)),DIMENSION(n2),INTENT(IN)::x2s
    LOGICAL,INTENT(IN),OPTIONAL::ONLYBEAM1
    ! as in PhotonPhotonFlux_pp_eval
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0,bupper=2d0
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::bv,mu2v
    REAL(KIND(1d0)),DIMENSION(n1+n2)::xr
    INTEGER,DIMENSION(n1+n2)::nbr,nmr
    REAL(KIND(1d0))::alo,a,gamma
    INTEGER::nsub,ib,i,j,k,nr,ldb,nbeam
    nbeam=2
    IF(PRESENT(ONLYBEAM1))THEN
       IF(ONLYBEAM1)nbeam=1
    ENDIF
    ldb=0
    DO i=1,MERGE(n1+n2,n1,nbeam.EQ.2)
       IF(i.LE.n1)THEN
          a=x1s(i)
       ELSE
          a=x2s(i-n1)
       ENDIF
       IF(a.LE.0d0.OR.a.GE.1d0)CYCLE
       CALL pp_inelastic_range(a,alo,nsub)
       ldb=MAX(ldb,nsub)
    ENDDO
    IF(ldb.EQ.0)RETURN
    ALLOCATE(bv(ldb,n1+n2),mu2v(MAX(n1,n2),n1+n2))
    nr=0
    DO ib=1,nbeam
       gamma=e_pn(ib)/mproton
       DO i=1,MERGE(n1,n2,ib.EQ.1)
          IF(ib.EQ.1)THEN
             a=x1s(i)
          ELSE
             a=x2s(i)
          ENDIF
          IF(a.LE.0d0.OR.a.GE.1d0)CYCLE
          nr=nr+1
          ! x as in PhotonNumberDensity_inelastic_proton: Ega/(gamma mp)
          xr(nr)=(e_pn(ib)*a)/(gamma*mproton)
          CALL pp_inelastic_range(a,alo,nsub)
          nbr(nr)=nsub
          DO k=1,nsub
             bv(k,nr)=DEXP((DBLE(2*nsub-2*k+1)*alo+DBLE(2*k-1)*bupper)/DBLE(2*nsub))/a/mproton
          ENDDO
          nmr(nr)=0
          DO j=1,MERGE(n2,n1,ib.EQ.1)
             IF(ib.EQ.1)THEN
                IF(x2s(j).LE.0d0.OR.x2s(j).GE.1d0)CYCLE
                nmr(nr)=nmr(nr)+1; mu2v(nmr(nr),nr)=photon_mu2(a,x2s(j))
             ELSE
                IF(x1s(j).LE.0d0.OR.x1s(j).GE.1d0)CYCLE
                nmr(nr)=nmr(nr)+1; mu2v(nmr(nr),nr)=photon_mu2(x1s(j),a)
             ENDIF
          ENDDO
          IF(nmr(nr).EQ.0)nr=nr-1
       ENDDO
    ENDDO
    IF(nr.GT.0)CALL inel_fast_prepare_rows(nr,xr(1:nr),ldb,nbr(1:nr),bv(:,1:nr),MAX(n1,n2),nmr(1:nr),mu2v(:,1:nr))
  END SUBROUTINE pp_inelastic_prepare

  SUBROUTINE pp_inelastic_range(x,alo,nsub)
    ! photon_final_state 1, 2: the lower limit alo of ln(x b mproton) in PhotonPhotonFlux_pp_eval, down to
    ! b = 1d-3/mu_max (mu_max = sqrt(x s), the largest W at this x), as the inelastic photon at k_T ~ mu lives at
    ! b ~ 1/mu, and the number of midpoint panels nsub, of width at most 0.25 (the midpoint rule in ln b converges as
    ! exp(-pi^2/width) for these integrands), with the old lower limit (where the elastic photon stops in
    ! PhotonPhotonFlux_pp_fxn) as a panel edge; it depends on x only
    IMPLICIT NONE
    include 'run90.inc'
    REAL(KIND(1d0)),INTENT(IN)::x
    REAL(KIND(1d0)),INTENT(OUT)::alo
    INTEGER,INTENT(OUT)::nsub
    ! as in PhotonPhotonFlux_pp_eval
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0,Rproton=0.877d0,bupper=2d0
    REAL(KIND(1d0))::bfact,a0,h
    INTEGER::m,k
    bfact=Rproton/GeVm12fm*mproton
    nsub=30
    IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
       bfact=bfact*LOWER_BFactor_Limit
       nsub=nsub*SUB_FACTOR
    ENDIF
    a0=DLOG(bfact*x)
    IF(USE_CHARGEFORMFACTOR4PHOTON.AND.USE_EXACT_ELASTIC4PHOTON)a0=a0+DLOG(EXACT_BLOW) ! as PhotonPhotonFlux_pp_eval
    m=CEILING((bupper-a0)/0.25d0-1d-9)
    h=(bupper-a0)/m
    alo=MIN(a0,DLOG(x*mproton*1d-3/DSQRT(x*4d0*e_pn(1)*e_pn(2))))
    k=CEILING((a0-alo)/h-1d-9)
    alo=a0-k*h
    nsub=m+k
  END SUBROUTINE pp_inelastic_range
  SUBROUTINE ffe_gauleg(n,x,w)
    ! Gauss-Legendre nodes and weights on [-1,1] (Newton iteration on P_n), cached for n <= 64
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::x,w
    REAL(KIND(1d0)),DIMENSION(64,64),SAVE::xc,wc
    LOGICAL,DIMENSION(64),SAVE::done=.FALSE.
    INTEGER::i,j,m,it
    REAL(KIND(1d0))::z,z1,p1,p2,p3,pp
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    IF(n.LT.1.OR.n.GT.64)THEN
       WRITE(*,*)"ERROR: ffe_gauleg: n out of range, n=",n
       STOP 1
    ENDIF
    IF(.NOT.done(n))THEN
       m=(n+1)/2
       DO i=1,m
          z=DCOS(pi*(DBLE(i)-0.25d0)/(DBLE(n)+0.5d0))
          DO it=1,100
             p1=1d0
             p2=0d0
             DO j=1,n
                p3=p2
                p2=p1
                p1=((2d0*j-1d0)*z*p2-(j-1d0)*p3)/DBLE(j)
             ENDDO
             pp=DBLE(n)*(z*p1-p2)/(z*z-1d0)
             z1=z
             z=z1-p1/pp
             IF(DABS(z-z1).LT.1d-16)EXIT
          ENDDO
          p1=1d0
          p2=0d0
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2d0*j-1d0)*z*p2-(j-1d0)*p3)/DBLE(j)
          ENDDO
          pp=DBLE(n)*(z*p1-p2)/(z*z-1d0)
          xc(i,n)=-z
          xc(n+1-i,n)=z
          wc(i,n)=2d0/((1d0-z*z)*pp*pp)
          wc(n+1-i,n)=wc(i,n)
       ENDDO
       IF(MOD(n,2).EQ.1)xc(m,n)=0d0
       done(n)=.TRUE.
    ENDIF
    x(1:n)=xc(1:n,n)
    w(1:n)=wc(1:n,n)
  END SUBROUTINE ffe_gauleg

  SUBROUTINE ffe_lagw(vg,v,p,i0,wt)
    ! Lagrange weights of the p-point stencil of the piecewise-uniform grid vg around v;
    ! the interpolant at v is sum_m wt(m) f(i0+m). Stencils never cross a segment boundary.
    IMPLICIT NONE
    TYPE(ffe_vgrid_t),INTENT(IN)::vg
    REAL(KIND(1d0)),INTENT(IN)::v
    INTEGER,INTENT(IN)::p
    INTEGER,INTENT(OUT)::i0
    REAL(KIND(1d0)),DIMENSION(p),INTENT(OUT)::wt
    INTEGER::k,lo,hi,mid,nloc,j0,m,l
    REAL(KIND(1d0))::u,num
    lo=1
    hi=vg%nseg
    DO WHILE(hi.GT.lo)
       mid=(lo+hi+1)/2
       IF(vg%sv0(mid).LE.v)THEN
          lo=mid
       ELSE
          hi=mid-1
       ENDIF
    ENDDO
    k=lo
    nloc=vg%shi(k)-vg%slo(k)+1
    u=(v-vg%sv0(k))/vg%sh(k)
    j0=FLOOR(u)-p/2+1
    j0=MAX(0,MIN(j0,nloc-p))
    CALL ffe_lagcore(u-DBLE(j0),p,wt)
    i0=vg%slo(k)-1+j0
  END SUBROUTINE ffe_lagw

  SUBROUTINE ffe_lagcore(u,p,wt)
    ! Lagrange weights on the nodes 0..p-1 at u: prefix/suffix products, O(p)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::u
    INTEGER,INTENT(IN)::p
    REAL(KIND(1d0)),DIMENSION(p),INTENT(OUT)::wt
    REAL(KIND(1d0)),DIMENSION(0:16)::pre,suf,den
    INTEGER,SAVE::pden=0
    REAL(KIND(1d0)),DIMENSION(16),SAVE::rden
    INTEGER::m,l
    IF(pden.NE.p)THEN
       DO m=0,p-1
          den(m)=1d0
          DO l=0,p-1
             IF(l.NE.m)den(m)=den(m)*DBLE(m-l)
          ENDDO
          rden(m+1)=1d0/den(m)
       ENDDO
       pden=p
    ENDIF
    pre(0)=1d0
    DO m=1,p-1
       pre(m)=pre(m-1)*(u-DBLE(m-1))
    ENDDO
    suf(p-1)=1d0
    DO m=p-2,0,-1
       suf(m)=suf(m+1)*(u-DBLE(m+1))
    ENDDO
    DO m=0,p-1
       wt(m+1)=pre(m)*suf(m)*rden(m+1)
    ENDDO
  END SUBROUTINE ffe_lagcore

  SUBROUTINE ffe_lagu(u,n,p,i0,wt)
    ! Lagrange weights on a uniform grid with nodes 0..n-1 (coordinate u in node units)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::u
    INTEGER,INTENT(IN)::n,p
    INTEGER,INTENT(OUT)::i0
    REAL(KIND(1d0)),DIMENSION(p),INTENT(OUT)::wt
    INTEGER::m,l
    REAL(KIND(1d0))::num
    i0=FLOOR(u)-p/2+1
    i0=MAX(0,MIN(i0,n-p))
    CALL ffe_lagcore(u-DBLE(i0),p,wt)
  END SUBROUTINE ffe_lagu

  SUBROUTINE ffe_sort(n,a)
    ! insertion sort (the lists are short)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(INOUT)::a
    INTEGER::i,j
    REAL(KIND(1d0))::t
    DO i=2,n
       t=a(i)
       j=i-1
       DO WHILE(j.GE.1)
          IF(a(j).LE.t)EXIT
          a(j+1)=a(j)
          j=j-1
       ENDDO
       a(j+1)=t
    ENDDO
  END SUBROUTINE ffe_sort

  SUBROUTINE ffe_breaks(n,a,lo,hi,eps,nout)
    ! sort a(1:n), keep the values strictly inside (lo,hi), add lo and hi, merge values closer
    ! than eps; a must have room for n+2 entries. On exit a(1)=lo < a(2) < ... < a(nout)=hi.
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(*),INTENT(INOUT)::a
    REAL(KIND(1d0)),INTENT(IN)::lo,hi,eps
    INTEGER,INTENT(OUT)::nout
    INTEGER::i,m
    m=0
    DO i=1,n
       IF(a(i).GT.lo+eps.AND.a(i).LT.hi-eps)THEN
          m=m+1
          a(m)=a(i)
       ENDIF
    ENDDO
    CALL ffe_sort(m,a(1:m))
    DO i=m,1,-1
       a(i+1)=a(i)
    ENDDO
    a(1)=lo
    nout=1
    DO i=2,m+1
       IF(a(i)-a(nout).GT.eps)THEN
          nout=nout+1
          a(nout)=a(i)
       ENDIF
    ENDDO
    nout=nout+1
    a(nout)=hi
  END SUBROUTINE ffe_breaks

  FUNCTION ffe_rho(r,ib)
    ! normalised Woods-Saxon (3pF) density of beam ib, r in GeV^-1, result in GeV^3
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_rho
    REAL(KIND(1d0)),INTENT(IN)::r
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0))::arg
    arg=(r-ffe_bm(ib)%R)/ffe_bm(ib)%a
    IF(arg.GT.700d0)THEN
       ffe_rho=0d0
    ELSE
       ffe_rho=ffe_bm(ib)%rho0*(1d0+ffe_bm(ib)%w*(r/ffe_bm(ib)%R)**2)/(1d0+DEXP(arg))
    ENDIF
  END FUNCTION ffe_rho

  SUBROUTINE ffe_rbreaks(ib,rmax,nr,rb)
    ! radial breakpoints of the Woods-Saxon density: coarse in the flat interior, a/2 beyond R-8a
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0)),INTENT(IN)::rmax
    INTEGER,INTENT(OUT)::nr
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,INTENT(OUT)::rb
    REAL(KIND(1d0))::r0,h
    INTEGER::k,nmax
    r0=MAX(0d0,ffe_bm(ib)%R-8d0*ffe_bm(ib)%a)
    h=0.5d0*ffe_bm(ib)%a
    nmax=10+INT((rmax-r0)/h)+2
    ALLOCATE(rb(nmax))
    nr=1
    rb(1)=0d0
    IF(r0.GT.0d0)THEN
       DO k=1,4
          nr=nr+1
          rb(nr)=r0*DBLE(k)/4d0
       ENDDO
    ENDIF
    DO k=1,nmax
       IF(r0+DBLE(k)*h.GE.rmax)EXIT
       nr=nr+1
       rb(nr)=r0+DBLE(k)*h
    ENDDO
    nr=nr+1
    rb(nr)=rmax
  END SUBROUTINE ffe_rbreaks

  SUBROUTINE ffe_normalise(ib)
    ! rho0 such that int d^3r rho = 1 (numerical, any w)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ib
    INTEGER::nr,i,g
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::rb
    REAL(KIND(1d0)),DIMENSION(12)::xg,wg
    REAL(KIND(1d0))::s,r,hh,cc
    REAL(KIND(1d0)),PARAMETER::fourpi=12.5663706143591729538505735331d0
    ffe_bm(ib)%rho0=1d0
    CALL ffe_rbreaks(ib,ffe_bm(ib)%R+60d0*ffe_bm(ib)%a,nr,rb)
    CALL ffe_gauleg(12,xg,wg)
    s=0d0
    DO i=1,nr-1
       hh=0.5d0*(rb(i+1)-rb(i))
       cc=0.5d0*(rb(i+1)+rb(i))
       DO g=1,12
          r=cc+hh*xg(g)
          s=s+hh*wg(g)*r*r*ffe_rho(r,ib)
       ENDDO
    ENDDO
    ffe_bm(ib)%rho0=1d0/(fourpi*s)
  END SUBROUTINE ffe_normalise

  SUBROUTINE ffe_build_vgrid(ib,vlo,vhi)
    ! piecewise-uniform grid in v = ln b: fine near the lower edge (where W of large x sits),
    ! fine around ln R for ChFF (the knee of the density), h0 = 0.05 elsewhere
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0)),INTENT(IN)::vlo,vhi
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(16)::bk
    INTEGER::nb,i,k,n,ntot,g,i0,m
    INTEGER,DIMENSION(16)::nint
    REAL(KIND(1d0))::mid,h,lnR,sc,a,b
    REAL(KIND(1d0)),DIMENSION(p)::wt
    REAL(KIND(1d0)),DIMENSION(8)::xg,wg
    sc=MAX(FFE_NodeScale,1d-3)
    nb=0
    nb=nb+1
    bk(nb)=vlo+0.1d0
    nb=nb+1
    bk(nb)=vlo+0.3d0
    nb=nb+1
    bk(nb)=vlo+1.0d0
    lnR=-100d0
    IF(ffe_bm(ib)%chff.AND..NOT.ffe_bm(ib)%proton)THEN
       lnR=DLOG(ffe_bm(ib)%R)
       nb=nb+1
       bk(nb)=lnR-0.7d0
       nb=nb+1
       bk(nb)=lnR+0.7d0
    ENDIF
    CALL ffe_breaks(nb,bk,vlo,vhi,1d-6,n)
    ntot=1
    DO i=1,n-1
       mid=0.5d0*(bk(i)+bk(i+1))
       h=0.05d0
       IF(mid.LT.vlo+0.1d0)THEN
          h=0.005d0
       ELSEIF(mid.LT.vlo+0.3d0)THEN
          h=0.01d0
       ELSEIF(mid.LT.vlo+1d0)THEN
          h=0.02d0
       ENDIF
       IF(DABS(mid-lnR).LT.0.7d0)h=MIN(h,0.02d0)
       h=h/sc
       nint(i)=MAX(p-1,CEILING((bk(i+1)-bk(i))/h-1d-9))
       ntot=ntot+nint(i)
    ENDDO
    IF(ALLOCATED(ffe_vg(ib)%v))THEN
       DEALLOCATE(ffe_vg(ib)%v,ffe_vg(ib)%b,ffe_vg(ib)%qw,ffe_vg(ib)%slo,ffe_vg(ib)%shi,&
            ffe_vg(ib)%sv0,ffe_vg(ib)%sh)
    ENDIF
    ffe_vg(ib)%n=ntot
    ffe_vg(ib)%nseg=n-1
    ALLOCATE(ffe_vg(ib)%v(ntot),ffe_vg(ib)%b(ntot),ffe_vg(ib)%qw(ntot))
    ALLOCATE(ffe_vg(ib)%slo(n-1),ffe_vg(ib)%shi(n-1),ffe_vg(ib)%sv0(n-1),ffe_vg(ib)%sh(n-1))
    k=1
    ffe_vg(ib)%v(1)=bk(1)
    DO i=1,n-1
       ffe_vg(ib)%slo(i)=k
       ffe_vg(ib)%sv0(i)=bk(i)
       ffe_vg(ib)%sh(i)=(bk(i+1)-bk(i))/DBLE(nint(i))
       DO m=1,nint(i)
          k=k+1
          ffe_vg(ib)%v(k)=bk(i)+DBLE(m)*ffe_vg(ib)%sh(i)
       ENDDO
       ffe_vg(ib)%v(k)=bk(i+1)
       ffe_vg(ib)%shi(i)=k
    ENDDO
    ffe_vg(ib)%b=DEXP(ffe_vg(ib)%v)
    ! quadrature weights of the interpolant: qw(j) = int L_j(v) dv
    ffe_vg(ib)%qw=0d0
    CALL ffe_gauleg(8,xg,wg)
    DO i=1,ntot-1
       a=ffe_vg(ib)%v(i)
       b=ffe_vg(ib)%v(i+1)
       DO g=1,8
          CALL ffe_lagw(ffe_vg(ib),0.5d0*(a+b)+0.5d0*(b-a)*xg(g),p,i0,wt)
          ffe_vg(ib)%qw(i0+1:i0+p)=ffe_vg(ib)%qw(i0+1:i0+p)+0.5d0*(b-a)*wg(g)*wt
       ENDDO
    ENDDO
  END SUBROUTINE ffe_build_vgrid
  SUBROUTINE ffe_chff_setup(ib)
    ! x-independent data of the position-space ChFF amplitude of nucleus ib:
    ! s-panels (Woods-Saxon breakpoints + every table node b_j < smax), s-nodes, z-nodes with
    ! weight*2*rho(sqrt(s^2+z^2)) and the z-moments of the cos(qz) expansion
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ib
    INTEGER,PARAMETER::nsg=4,nzg=6,nmom=12
    INTEGER::nr,np,npts,i,j,k,g,is,nz,cnt,lo,hi,mid
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::rb,pk,zb
    REAL(KIND(1d0))::smax,s,z0,z1,zz,hh,cc
    REAL(KIND(1d0)),DIMENSION(nsg)::xs,wsg
    REAL(KIND(1d0)),DIMENSION(nzg)::xz,wz
    TYPE(ffe_chff_t)::cf
    smax=ffe_bm(ib)%R+40d0*ffe_bm(ib)%a
    cf%smax=smax
    CALL ffe_rbreaks(ib,smax,nr,rb)
    ALLOCATE(pk(nr+ffe_vg(ib)%n+2))
    np=nr
    pk(1:nr)=rb(1:nr)
    DO j=1,ffe_vg(ib)%n
       IF(ffe_vg(ib)%b(j).LT.smax)THEN
          np=np+1
          pk(np)=ffe_vg(ib)%b(j)
       ENDIF
    ENDDO
    CALL ffe_breaks(np,pk,0d0,smax,1d-12*smax,npts)
    cf%npan=npts-1
    ALLOCATE(cf%pa(cf%npan),cf%pb(cf%npan))
    cf%pa(1:cf%npan)=pk(1:npts-1)
    cf%pb(1:cf%npan)=pk(2:npts)
    ! table node -> panel ending at it (nearest panel end)
    ALLOCATE(cf%jpan(ffe_vg(ib)%n))
    cf%jpan=0
    DO j=1,ffe_vg(ib)%n
       IF(ffe_vg(ib)%b(j).GE.smax)CYCLE
       lo=1
       hi=cf%npan
       DO WHILE(hi.GT.lo)
          mid=(lo+hi)/2
          IF(cf%pb(mid).LT.ffe_vg(ib)%b(j)*(1d0-1d-10))THEN
             lo=mid+1
          ELSE
             hi=mid
          ENDIF
       ENDDO
       cf%jpan(j)=lo
    ENDDO
    ! s-nodes
    CALL ffe_gauleg(nsg,xs,wsg)
    CALL ffe_gauleg(nzg,xz,wz)
    cf%ns=nsg*cf%npan
    ALLOCATE(cf%s(cf%ns),cf%ws(cf%ns),cf%zoff(cf%ns+1))
    DO i=1,cf%npan
       hh=0.5d0*(cf%pb(i)-cf%pa(i))
       cc=0.5d0*(cf%pb(i)+cf%pa(i))
       DO g=1,nsg
          cf%s((i-1)*nsg+g)=cc+hh*xs(g)
          cf%ws((i-1)*nsg+g)=hh*wsg(g)
       ENDDO
    ENDDO
    ! z-nodes: breakpoints z_k = sqrt(r_k^2-s^2) for the radial breakpoints r_k > s
    ALLOCATE(zb(nr+2))
    cnt=0
    DO is=1,cf%ns
       s=cf%s(is)
       nz=0
       DO k=1,nr
          IF(rb(k).GT.s)THEN
             nz=nz+1
             zb(nz)=DSQRT((rb(k)-s)*(rb(k)+s))
          ENDIF
       ENDDO
       cnt=cnt+nzg*nz
    ENDDO
    cf%nz=cnt
    cf%nmom=nmom
    ALLOCATE(cf%z(cnt),cf%wrho(cnt),cf%mom(0:nmom,cf%ns))
    cnt=0
    DO is=1,cf%ns
       s=cf%s(is)
       cf%zoff(is)=cnt
       z0=0d0
       DO k=1,nr
          IF(rb(k).LE.s)CYCLE
          z1=DSQRT((rb(k)-s)*(rb(k)+s))
          hh=0.5d0*(z1-z0)
          cc=0.5d0*(z1+z0)
          DO g=1,nzg
             zz=cc+hh*xz(g)
             cnt=cnt+1
             cf%z(cnt)=zz
             cf%wrho(cnt)=2d0*hh*wz(g)*ffe_rho(DSQRT(s*s+zz*zz),ib)
          ENDDO
          z0=z1
       ENDDO
       DO k=0,nmom
          cf%mom(k,is)=SUM(cf%wrho(cf%zoff(is)+1:cnt)*cf%z(cf%zoff(is)+1:cnt)**(2*k))
       ENDDO
    ENDDO
    cf%zoff(cf%ns+1)=cnt
    ffe_cf(ib)=cf
  END SUBROUTINE ffe_chff_setup

  SUBROUTINE ffe_W_row(ib,x,W,OWNQ)
    ! W_j = b_j^2 N(x, b_j) on the v-grid of beam ib (unit charge, alpha = 1); OWNQ: the closed-form dipole of the
    ! proton also with the options of the inelastic branch (the proxies of ffe_build_K)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0)),INTENT(IN)::x
    REAL(KIND(1d0)),DIMENSION(ffe_vg(ib)%n),INTENT(OUT)::W
    LOGICAL,INTENT(IN),OPTIONAL::OWNQ
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    REAL(KIND(1d0)),EXTERNAL::GUPC_BESSK0E,GUPC_BESSK1E,GUPC_BESSI0E,GUPC_BESSI1E
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::sig,Ch,Dh
    INTEGER::j,n,is,k,ip,nsg
    REAL(KIND(1d0))::q,b,z,ig2,Ap,Al,qb,e2,fac,cq,qq,ss,ap1,aa,qp,qpb
    LOGICAL::epa,own
    n=ffe_vg(ib)%n
    q=x*ffe_bm(ib)%mass
    ig2=1d0/ffe_bm(ib)%gamma**2
    W=0d0
    IF(q.LE.0d0)RETURN
    own=.FALSE.
    IF(PRESENT(OWNQ))own=OWNQ
    IF(ffe_bm(ib)%proton.AND.ffe_bm(ib)%chff.AND..NOT.own.AND.(USE_EXACT_ELASTIC4PHOTON.OR.proton_FF_correction.NE.0&
         .OR.alpha_running4photon))THEN
       ! proton with the options of the inelastic branch: the densities of ElasticFlux (exact distribution or
       ! small-x approximation, dipole/dde1/luxqed form factors, running alpha(Q2) inside), normalised to
       ! alpha(0); the factor alpha(0)/alpha(mu2) of each photon is applied per (x1,x2) (ffe_anorm2)
       IF(x.GE.1d0)RETURN
       IF(USE_EXACT_ELASTIC4PHOTON)THEN
          ! the whole row at once (proton_exact_density_row); W = 0 where q b > 340, as below
          k=COUNT(q*ffe_vg(ib)%b(1:n).LE.340d0)
          CALL proton_exact_density_row(x,proton_FF_correction,k,ffe_vg(ib)%b(1:k),W(1:k))
          W(1:k)=W(1:k)*ffe_vg(ib)%b(1:k)**2
          RETURN
       ENDIF
       k=COUNT(q*ffe_vg(ib)%b(1:n).LE.340d0)
       IF(alpha_running4photon)THEN
          ! the running-alpha part of the dipole amplitude for the whole row (its poles once per x)
          ALLOCATE(sig(k))
          CALL dipole_running_amplitude_row(x*ffe_bm(ib)%ebeam/ffe_bm(ib)%gamma,k,ffe_vg(ib)%b(1:k),sig)
       ENDIF
       DO j=1,k
          b=ffe_vg(ib)%b(j)
          IF(alpha_running4photon)THEN
             W(j)=b*b*PhotonNumberDensity_ChargeFormFactor_proton(b,x*ffe_bm(ib)%ebeam,ffe_bm(ib)%gamma,&
                  proton_FF_correction,RUNAMP=sig(j))
          ELSE
             W(j)=b*b*PhotonNumberDensity_ChargeFormFactor_proton(b,x*ffe_bm(ib)%ebeam,ffe_bm(ib)%gamma,&
                  proton_FF_correction)
          ENDIF
       ENDDO
       RETURN
    ENDIF
    IF(ffe_bm(ib)%proton.AND.ffe_bm(ib)%chff)THEN
       ! proton, dipole form factor (closed form): q' = sqrt(q^2+1/a^2), a = 1/sqrt(0.71 GeV^2)
       aa=1.1867816581938533d0
       qp=DSQRT(q*q+1d0/aa**2)
       DO j=1,n
          b=ffe_vg(ib)%b(j)
          qb=q*b
          qpb=qp*b
          IF(qb.GT.340d0)CYCLE
          Ap=q*GUPC_BESSK1E(qb)*DEXP(-qb)-qp*GUPC_BESSK1E(qpb)*DEXP(-qpb)&
               -b/(2d0*aa**2)*GUPC_BESSK0E(qpb)*DEXP(-qpb)
          Al=q*(GUPC_BESSK0E(qb)*DEXP(-qb)-GUPC_BESSK0E(qpb)*DEXP(-qpb)&
               -b/(2d0*aa**2*qp)*GUPC_BESSK1E(qpb)*DEXP(-qpb))
          IF(FFE_Legacy_ChFF_Definitions)Al=0d0
          W(j)=b*b*(Ap*Ap+Al*Al*ig2)/pi**2
       ENDDO
       RETURN
    ENDIF
    IF(.NOT.ffe_bm(ib)%chff)THEN
       ! EDFF: point-like, b >= bmin is the lower edge of the grid
       DO j=1,n
          z=q*ffe_vg(ib)%b(j)
          IF(z.GT.340d0)CYCLE
          W(j)=z*z*DEXP(-2d0*z)*(GUPC_BESSK1E(z)**2+ig2*GUPC_BESSK0E(z)**2)/pi**2
       ENDDO
       RETURN
    ENDIF
    ! ChFF nucleus, position space (see the header of this section)
    ALLOCATE(sig(ffe_cf(ib)%ns),Ch(0:ffe_cf(ib)%npan),Dh(1:ffe_cf(ib)%npan+1))
    IF(q*ffe_cf(ib)%smax.LE.0.4d0)THEN
       qq=-q*q
       DO is=1,ffe_cf(ib)%ns
          ss=0d0
          cq=1d0
          DO k=0,ffe_cf(ib)%nmom
             ss=ss+cq*ffe_cf(ib)%mom(k,is)
             cq=cq*qq/DBLE((2*k+1)*(2*k+2))
          ENDDO
          sig(is)=ss
       ENDDO
    ELSE
       DO is=1,ffe_cf(ib)%ns
          k=ffe_cf(ib)%zoff(is)
          sig(is)=SUM(ffe_cf(ib)%wrho(k+1:ffe_cf(ib)%zoff(is+1))&
               *DCOS(q*ffe_cf(ib)%z(k+1:ffe_cf(ib)%zoff(is+1))))
       ENDDO
    ENDIF
    nsg=ffe_cf(ib)%ns/ffe_cf(ib)%npan
    ! Chat(b) = exp(-q b) C(b), forward; Dhat(b) = exp(q b) D(b), backward
    Ch(0)=0d0
    DO ip=1,ffe_cf(ib)%npan
       ss=0d0
       DO k=1,nsg
          is=(ip-1)*nsg+k
          ss=ss+ffe_cf(ib)%ws(is)*ffe_cf(ib)%s(is)*sig(is)*GUPC_BESSI0E(q*ffe_cf(ib)%s(is))&
               *DEXP(-q*(ffe_cf(ib)%pb(ip)-ffe_cf(ib)%s(is)))
       ENDDO
       Ch(ip)=DEXP(-q*(ffe_cf(ib)%pb(ip)-ffe_cf(ib)%pa(ip)))*Ch(ip-1)+ss
    ENDDO
    Dh(ffe_cf(ib)%npan+1)=0d0
    DO ip=ffe_cf(ib)%npan,1,-1
       ss=0d0
       DO k=1,nsg
          is=(ip-1)*nsg+k
          ss=ss+ffe_cf(ib)%ws(is)*ffe_cf(ib)%s(is)*sig(is)*GUPC_BESSK0E(q*ffe_cf(ib)%s(is))&
               *DEXP(-q*(ffe_cf(ib)%s(is)-ffe_cf(ib)%pa(ip)))
       ENDDO
       Dh(ip)=DEXP(-q*(ffe_cf(ib)%pb(ip)-ffe_cf(ib)%pa(ip)))*Dh(ip+1)+ss
    ENDDO
    DO j=1,n
       b=ffe_vg(ib)%b(j)
       qb=q*b
       IF(qb.GT.340d0)CYCLE
       epa=b.GE.ffe_cf(ib)%smax
       IF(FFE_Legacy_ChFF_Definitions)THEN
          epa=epa.OR.b.GT.3d0*ffe_bm(ib)%R.OR.(x.LT.1d-4.AND.b.GT.ffe_bm(ib)%R)
       ENDIF
       IF(epa)THEN
          e2=DEXP(-2d0*qb)
          W(j)=b*b*q*q*e2*(GUPC_BESSK1E(qb)**2+ig2*GUPC_BESSK0E(qb)**2)/pi**2
       ELSE
          ip=ffe_cf(ib)%jpan(j)
          Ap=twopi*q*(GUPC_BESSK1E(qb)*Ch(ip)-GUPC_BESSI1E(qb)*Dh(ip+1))
          Al=twopi*q*(GUPC_BESSK0E(qb)*Ch(ip)+GUPC_BESSI0E(qb)*Dh(ip+1))
          IF(FFE_Legacy_ChFF_Definitions)Al=0d0
          W(j)=b*b*(Ap*Ap+Al*Al*ig2)/pi**2
       ENDIF
    ENDDO
    IF(FFE_Verbose.GE.3)THEN
       ap1=twopi*Ch(ffe_cf(ib)%npan)*DEXP(q*ffe_cf(ib)%smax)
       WRITE(*,*)"FFE: ChFF normalisation 2 pi C(smax) - 1 at x =",x,":",ap1-1d0
    ENDIF
    DEALLOCATE(sig,Ch,Dh)
  END SUBROUTINE ffe_W_row

  FUNCTION ffe_TA(b,ib)
    ! thickness function of nucleus ib (normalised to 1), b in GeV^-1, result in GeV^2
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_TA
    REAL(KIND(1d0)),INTENT(IN)::b
    INTEGER,INTENT(IN)::ib
    INTEGER::nr,k,g
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::rb
    REAL(KIND(1d0)),DIMENSION(8)::xg,wg
    REAL(KIND(1d0))::z0,z1,hh,cc,zz,smax
    smax=ffe_bm(ib)%R+60d0*ffe_bm(ib)%a
    ffe_TA=0d0
    IF(b.GE.smax)RETURN
    CALL ffe_rbreaks(ib,smax,nr,rb)
    CALL ffe_gauleg(8,xg,wg)
    z0=0d0
    DO k=1,nr
       IF(rb(k).LE.b)CYCLE
       z1=DSQRT((rb(k)-b)*(rb(k)+b))
       hh=0.5d0*(z1-z0)
       cc=0.5d0*(z1+z0)
       DO g=1,8
          zz=cc+hh*xg(g)
          ffe_TA=ffe_TA+2d0*hh*wg(g)*ffe_rho(DSQRT(b*b+zz*zz),ib)
       ENDDO
       z0=z1
    ENDDO
  END FUNCTION ffe_TA

  SUBROUTINE ffe_build_ptable()
    ! chi(b) = sigma_NN * A * B * That_AB(b) (AB) or sigma_NN * A * That_A(b) (pA) on a uniform
    ! b grid (GeV^-1); That_AB by the Hankel transform of the numerical Woods-Saxon form factors.
    ! P_had = exp(-chi) for b < ffe_bcut, 1 beyond (chi < 1e-14 there).
    IMPLICIT NONE
    INTEGER::nb,i,k,g,nk,nr,ia,nkp,icut,ib,inuc
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::bt,chi,kk,wk,FA,FB,rb
    REAL(KIND(1d0)),DIMENSION(8)::xg,wg
    REAL(KIND(1d0))::amin,amax,bend,kmax,dk,hh,cc,r,s,kr,TAB,b1,b2,bm
    REAL(KIND(1d0)),DIMENSION(2)::Aarr,Zarr
    REAL(KIND(1d0)),PARAMETER::fourpi=12.5663706143591729538505735331d0
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    LOGICAL::jump
    CALL ffe_gauleg(8,xg,wg)
    jump=.FALSE.
    IF(ffe_icoll.EQ.1)THEN
       ! pp: P = (1-exp(-u))^2, u = b^2/(2 b0), in closed form (ffe_Phad); 1-P <= 2 exp(-u) < 1e-12
       ! beyond bcut. The ln chi table (chi = -ln P) only serves the test hooks.
       ffe_pjump=.FALSE.
       ffe_bcut=DSQRT(2d0*ffe_b0pp*DLOG(2d12))
       ffe_bcut_fm=ffe_bcut*GeVm12fm
       ffe_hb=DSQRT(ffe_b0pp)/16d0
       ffe_nP=INT(ffe_bcut/ffe_hb)+1
       IF(ALLOCATED(ffe_lnchi))DEALLOCATE(ffe_lnchi)
       ALLOCATE(ffe_lnchi(ffe_nP))
       DO i=1,ffe_nP
          b1=ffe_Phad(DBLE(i-1)*ffe_hb)
          IF(b1.LE.0d0)THEN
             ffe_lnchi(i)=690d0
          ELSE
             ffe_lnchi(i)=DLOG(MAX(-DLOG(b1),1d-30))
          ENDIF
       ENDDO
       ! breakpoint spacing: P varies on the scale sqrt(2 b0)
       CALL ffe_build_beta(0.5d0*DSQRT(ffe_b0pp))
       RETURN
    ENDIF
    IF(ffe_icoll.EQ.3)THEN
       amin=MIN(ffe_bm(1)%a,ffe_bm(2)%a)
       amax=MAX(ffe_bm(1)%a,ffe_bm(2)%a)
       bend=ffe_bm(1)%R+ffe_bm(2)%R+60d0*amax
    ELSE
       amin=ffe_bm(2)%a
       amax=amin
       bend=ffe_bm(2)%R+60d0*amax
    ENDIF
    ffe_hb=amin/8d0
    nb=INT(bend/ffe_hb)+2
    ALLOCATE(bt(nb),chi(nb))
    DO i=1,nb
       bt(i)=DBLE(i-1)*ffe_hb
    ENDDO
    IF(ffe_icoll.EQ.3.AND.use_MC_Glauber)THEN
       Aarr=(/ffe_bm(1)%anum,ffe_bm(2)%anum/)
       Zarr=(/ffe_bm(1)%znum,ffe_bm(2)%znum/)
       DO i=1,nb
          ! TAB_MC_Glauber takes b in fm and returns T_AB in fm^-2
          chi(i)=ffe_sigNN*TAB_MC_Glauber(bt(i)*GeVm12fm,0d0,Aarr,Zarr)*GeVm12fm**2
       ENDDO
       jump=.TRUE.
    ELSEIF(ffe_icoll.EQ.2.AND.use_MC_Glauber)THEN
       DO i=1,nb
          chi(i)=ffe_sigNN*TA_MC_Glauber(bt(i)*GeVm12fm,0d0,ffe_bm(2)%anum)*GeVm12fm**2
       ENDDO
       jump=.TRUE.
    ELSEIF(ffe_icoll.EQ.2)THEN
       DO i=1,nb
          chi(i)=ffe_sigNN*ffe_bm(2)%anum*ffe_TA(bt(i),2)
       ENDDO
    ELSE
       ! numerical form factors F(k) = 4 pi int r^2 rho j0(kr) dr on the k-nodes
       kmax=40d0/(3.14159265358979323846d0*(ffe_bm(1)%a+ffe_bm(2)%a))
       dk=1.5707963267948966d0/bend
       nkp=INT(kmax/dk)+1
       nk=8*nkp
       ALLOCATE(kk(nk),wk(nk),FA(nk),FB(nk))
       DO i=1,nkp
          DO g=1,8
             kk((i-1)*8+g)=(DBLE(i)-0.5d0+0.5d0*xg(g))*dk
             wk((i-1)*8+g)=0.5d0*dk*wg(g)
          ENDDO
       ENDDO
       DO ib=1,2
          IF(ib.EQ.2.AND.ffe_bm(2)%anum.EQ.ffe_bm(1)%anum.AND.ffe_bm(2)%znum.EQ.ffe_bm(1)%znum)THEN
             FB=FA
             EXIT
          ENDIF
          CALL ffe_rbreaks(ib,ffe_bm(ib)%R+60d0*ffe_bm(ib)%a,nr,rb)
          DO k=1,nk
             s=0d0
             DO i=1,nr-1
                ! sub-panels resolve sin(kr) as well
                ia=MAX(1,CEILING((rb(i+1)-rb(i))*kk(k)/1.5d0))
                hh=0.5d0*(rb(i+1)-rb(i))/DBLE(ia)
                DO inuc=1,ia
                   cc=rb(i)+(2d0*DBLE(inuc)-1d0)*hh
                   DO g=1,8
                      r=cc+hh*xg(g)
                      kr=kk(k)*r
                      s=s+hh*wg(g)*r*DSIN(kr)/kk(k)*ffe_rho(r,ib)
                   ENDDO
                ENDDO
             ENDDO
             IF(ib.EQ.1)THEN
                FA(k)=fourpi*s
             ELSE
                FB(k)=fourpi*s
             ENDIF
          ENDDO
          DEALLOCATE(rb)
       ENDDO
       DO i=1,nb
          TAB=0d0
          DO k=1,nk
             TAB=TAB+wk(k)*kk(k)*BESSEL_J0(kk(k)*bt(i))*FA(k)*FB(k)
          ENDDO
          chi(i)=ffe_sigNN*ffe_bm(1)%anum*ffe_bm(2)%anum*TAB/twopi
       ENDDO
    ENDIF
    ! P = 1 (to 1e-12) from the first node beyond the maximum of chi where chi < 1e-12 (for w < 0
    ! nuclei chi turns negative further out: those nodes are not kept in the table)
    icut=nb
    k=MAXLOC(chi(1:nb),DIM=1)
    DO i=k,nb
       IF(chi(i).LT.1d-12)THEN
          icut=i-1
          EXIT
       ENDIF
    ENDDO
    ffe_pjump=jump
    IF(jump)THEN
       ! closed-form fit with a hard cut-off: P jumps to 1 at the first zero
       bm=bt(icut)
       b1=bt(icut)
       b2=bt(MIN(icut+1,nb))
       DO k=1,60
          bm=0.5d0*(b1+b2)
          IF(ffe_icoll.EQ.3)THEN
             TAB=TAB_MC_Glauber(bm*GeVm12fm,0d0,Aarr,Zarr)
          ELSE
             TAB=TA_MC_Glauber(bm*GeVm12fm,0d0,ffe_bm(2)%anum)
          ENDIF
          IF(TAB.GT.0d0)THEN
             b1=bm
          ELSE
             b2=bm
          ENDIF
       ENDDO
       ffe_bcut=b1
       ffe_nP=icut
    ELSE
       ffe_nP=icut
       ffe_bcut=bt(MIN(nb,icut+1))
    ENDIF
    ffe_bcut_fm=ffe_bcut*GeVm12fm
    IF(ALLOCATED(ffe_lnchi))DEALLOCATE(ffe_lnchi)
    ALLOCATE(ffe_lnchi(ffe_nP))
    DO i=1,ffe_nP
       ffe_lnchi(i)=DLOG(MAX(chi(i),1d-30))
    ENDDO
    CALL ffe_build_beta(amin)
  END SUBROUTINE ffe_build_ptable

  FUNCTION ffe_Phad(b)
    ! hadronic survival probability at impact parameter b (GeV^-1)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_Phad
    REAL(KIND(1d0)),INTENT(IN)::b
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(p)::wt
    INTEGER::i0
    REAL(KIND(1d0))::y
    IF(b.GE.ffe_bcut)THEN
       ffe_Phad=1d0
       RETURN
    ENDIF
    IF(ffe_icoll.EQ.1)THEN
       ! pp, closed form (1-exp(-u))^2 (series for small u)
       y=0.5d0*b*b/ffe_b0pp
       IF(y.LT.1d-3)THEN
          y=y*(1d0-y*(0.5d0-y*(1d0/6d0-y/24d0)))
       ELSE
          y=1d0-DEXP(-y)
       ENDIF
       ffe_Phad=y*y
       RETURN
    ENDIF
    CALL ffe_lagu(b/ffe_hb,ffe_nP,p,i0,wt)
    y=SUM(wt*ffe_lnchi(i0+1:i0+p))
    IF(y.GT.6.6d0)THEN
       ffe_Phad=0d0
    ELSE
       ffe_Phad=MIN(1d0,DEXP(-DEXP(y)))
    ENDIF
  END FUNCTION ffe_Phad

  FUNCTION ffe_Ptag8(b)
    ! neutron-tagging probability at impact parameter b (GeV^-1), 8-point Lagrange on its ln b table
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_Ptag8
    REAL(KIND(1d0)),INTENT(IN)::b
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(p)::wt
    INTEGER::i0
    REAL(KIND(1d0))::u
    u=(DLOG(b)-ffe_Tv0)/ffe_Th
    u=MAX(0d0,MIN(u,DBLE(ffe_nT-1)))
    CALL ffe_lagu(u,ffe_nT,p,i0,wt)
    ffe_Ptag8=SUM(wt*ffe_PtagT(i0+1:i0+p))
  END FUNCTION ffe_Ptag8

  SUBROUTINE ffe_build_tagfine()
    ! the 8-point interpolant of the P_tag table sampled 4x finer, for the 6-point lookups of
    ! ffe_Ptag inside the kernel (the general 8-point weights dominated the tagged K build)
    IMPLICIT NONE
    INTEGER::i
    ffe_ThF=0.25d0*ffe_Th
    ffe_nTF=4*(ffe_nT-1)+1
    IF(ALLOCATED(ffe_PtagF))DEALLOCATE(ffe_PtagF)
    ALLOCATE(ffe_PtagF(ffe_nTF))
    DO i=1,ffe_nTF
       ffe_PtagF(i)=ffe_Ptag8(DEXP(ffe_Tv0+DBLE(i-1)*ffe_ThF))
    ENDDO
  END SUBROUTINE ffe_build_tagfine

  FUNCTION ffe_Ptag(b)
    ! neutron-tagging probability at impact parameter b (GeV^-1): 6-point Lagrange on the fine
    ! ln b table (constant beyond its ends)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_Ptag
    REAL(KIND(1d0)),INTENT(IN)::b
    REAL(KIND(1d0))::u,x0,x1,x2,x3,x4,x5
    INTEGER::i0
    IF(.NOT.ffe_tag)THEN
       ffe_Ptag=1d0
       RETURN
    ENDIF
    u=(DLOG(b)-ffe_Tv0)/ffe_ThF
    u=MAX(0d0,MIN(u,DBLE(ffe_nTF-1)))
    i0=MAX(0,MIN(INT(u)-2,ffe_nTF-6))
    x0=u-DBLE(i0)
    x1=x0-1d0
    x2=x0-2d0
    x3=x0-3d0
    x4=x0-4d0
    x5=x0-5d0
    ffe_Ptag=-ffe_PtagF(i0+1)*x1*x2*x3*x4*x5/120d0+ffe_PtagF(i0+2)*x0*x2*x3*x4*x5/24d0&
         -ffe_PtagF(i0+3)*x0*x1*x3*x4*x5/12d0+ffe_PtagF(i0+4)*x0*x1*x2*x4*x5/12d0&
         -ffe_PtagF(i0+5)*x0*x1*x2*x3*x5/24d0+ffe_PtagF(i0+6)*x0*x1*x2*x3*x4/120d0
  END FUNCTION ffe_Ptag

  SUBROUTINE ffe_build_beta(amin)
    ! breakpoints of the kernel: impact parameters where P_had crosses fixed levels, a/2
    ! subdivisions up to bcut, and (tagging) a geometric continuation to the end of the domain
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::amin
    INTEGER,PARAMETER::nlev=18
    REAL(KIND(1d0)),DIMENSION(nlev),PARAMETER::lev=(/1d-14,1d-10,1d-6,1d-4,1d-3,1d-2,5d-2,0.15d0,&
         0.3d0,0.5d0,0.7d0,0.85d0,0.95d0,0.99d0,0.999d0,0.9999d0,0.999999d0,0.99999999d0/)
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::bb
    INTEGER::nmax,n,i,k,j,nsub
    REAL(KIND(1d0))::b1,b2,bm,dmax,bstart
    nmax=nlev+INT(ffe_bcut/(0.5d0*amin))+400
    ALLOCATE(bb(nmax))
    n=0
    DO k=1,nlev
       ! first table node with P >= level, then bisection
       b2=-1d0
       DO i=1,ffe_nP
          IF(DBLE(i-1)*ffe_hb.GE.ffe_bcut)EXIT
          IF(ffe_Phad(DBLE(i-1)*ffe_hb).GE.lev(k))THEN
             b2=DBLE(i-1)*ffe_hb
             EXIT
          ENDIF
       ENDDO
       IF(b2.LE.0d0)CYCLE
       b1=MAX(0d0,b2-ffe_hb)
       DO j=1,50
          bm=0.5d0*(b1+b2)
          IF(ffe_Phad(bm).GE.lev(k))THEN
             b2=bm
          ELSE
             b1=bm
          ENDIF
       ENDDO
       n=n+1
       bb(n)=0.5d0*(b1+b2)
    ENDDO
    n=n+1
    bb(n)=ffe_bcut
    CALL ffe_sort(n,bb(1:n))
    ! subdivide: consecutive breakpoints <= a/2 apart up to the last level, <= 3a in the tail
    bstart=bb(1)
    nsub=n
    DO i=1,nsub-1
       dmax=0.5d0*amin
       IF(i.GE.nsub-1)dmax=3d0*amin
       IF(bb(i+1)-bb(i).GT.dmax)THEN
          k=CEILING((bb(i+1)-bb(i))/dmax)
          DO j=1,k-1
             n=n+1
             bb(n)=bb(i)+(bb(i+1)-bb(i))*DBLE(j)/DBLE(k)
          ENDDO
       ENDIF
    ENDDO
    IF(ffe_tag)THEN
       dmax=DEXP(ffe_vg(1)%v(ffe_vg(1)%n))+DEXP(ffe_vg(2)%v(ffe_vg(2)%n))
       bm=ffe_bcut
       DO WHILE(bm.LT.dmax.AND.n.LT.nmax)
          bm=bm*2d0
          n=n+1
          bb(n)=bm
       ENDDO
    ENDIF
    CALL ffe_sort(n,bb(1:n))
    IF(ALLOCATED(ffe_beta))DEALLOCATE(ffe_beta)
    ALLOCATE(ffe_beta(n))
    ffe_beta(1:n)=bb(1:n)
    ffe_nbeta=n
  END SUBROUTINE ffe_build_beta
  SUBROUTINE ffe_build_pfast()
    ! fine table of P(b) itself (uniform in b, step a/32) for fast 6-point lookups in Phi, and the
    ! moment tables g_2n(d) = int_0^tauB tau^2n (1-P(sqrt(d^2+tau^2))) [* b12 in fm] dtau of the
    ! far-zone expansion Phi = 2pi - (2/sqrt(G)) sum_n c_n g_2n(|d|)/(4G)^n, G = b1 b2 >= bcut^2
    IMPLICIT NONE
    INTEGER::i,k,n,np,g,m
    REAL(KIND(1d0))::d,tb,t0,t1,hh,cc,tau,b12,om,pw,wv
    REAL(KIND(1d0)),DIMENSION(8)::xg,wg
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::tk
    ffe_hf=MIN(ffe_bm(1)%a,ffe_bm(2)%a)/32d0
    IF(ffe_icoll.EQ.2)ffe_hf=ffe_bm(2)%a/32d0
    IF(ffe_icoll.EQ.1)ffe_hf=DSQRT(ffe_b0pp)/64d0
    ffe_nf=INT(ffe_bcut/ffe_hf)+8
    ! a jump at bcut: keep the stencils on the nodes below it
    IF(ffe_pjump)ffe_nf=INT(ffe_bcut/ffe_hf*(1d0-1d-12))+1
    IF(ALLOCATED(ffe_Pf))DEALLOCATE(ffe_Pf)
    ALLOCATE(ffe_Pf(ffe_nf))
    DO i=1,ffe_nf
       ffe_Pf(i)=ffe_Phad(DBLE(i-1)*ffe_hf)
    ENDDO
    ! moment tables
    ffe_hg=ffe_hf*2d0
    ffe_ng=INT(ffe_bcut/ffe_hg)+8
    IF(ALLOCATED(ffe_gt))DEALLOCATE(ffe_gt)
    ALLOCATE(ffe_gt(0:FFE_NMOM,ffe_ng,2),tk(ffe_nbeta+3))
    CALL ffe_gauleg(8,xg,wg)
    ffe_gt=0d0
    DO i=1,ffe_ng
       d=DBLE(i-1)*ffe_hg
       IF(d.GE.ffe_bcut)CYCLE
       tb=DSQRT((ffe_bcut-d)*(ffe_bcut+d))
       np=0
       DO k=1,ffe_nbeta
          IF(ffe_beta(k).GT.d.AND.ffe_beta(k).LT.ffe_bcut)THEN
             np=np+1
             tk(np)=DSQRT((ffe_beta(k)-d)*(ffe_beta(k)+d))
          ENDIF
       ENDDO
       CALL ffe_breaks(np,tk,0d0,tb,1d-12*tb,n)
       DO k=1,n-1
          t0=tk(k)
          t1=tk(k+1)
          hh=0.5d0*(t1-t0)
          cc=0.5d0*(t1+t0)
          DO g=1,8
             tau=cc+hh*xg(g)
             b12=DSQRT(d*d+tau*tau)
             om=1d0-ffe_Pfast(b12)
             pw=hh*wg(g)*om
             DO m=0,FFE_NMOM
                ffe_gt(m,i,1)=ffe_gt(m,i,1)+pw
                ffe_gt(m,i,2)=ffe_gt(m,i,2)+pw*b12*GeVm12fm
                pw=pw*tau*tau
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! tagging: inner moments h_2n(d) = int_0^tauB tau^2n P_had P_tag (b12 < bcut) [* b12 in fm];
    ! h_2n = tauB^(2n+1) H_n(w), w = tauB^2 = bcut^2 - d^2, with H_n = int_0^1 sigma^2n P(b12) dsigma,
    ! b12^2 = d^2 + w sigma^2, smooth in w (the edge P(bcut) P_tag(bcut) /= 0 is in tauB^(2n+1))
    IF(ffe_tag)THEN
       ! w in [0, bcut^2] (the stencils are one-sided at the ends)
       ffe_nw=INT(ffe_bcut/ffe_hg)+1
       ffe_hw=ffe_bcut**2/DBLE(ffe_nw-1)
       IF(ALLOCATED(ffe_ht))DEALLOCATE(ffe_ht)
       ALLOCATE(ffe_ht(0:FFE_NMOM,ffe_nw,2))
       ffe_ht=0d0
       DO i=1,ffe_nw
          wv=MIN(DBLE(i-1)*ffe_hw,ffe_bcut**2)
          d=DSQRT(MAX(0d0,ffe_bcut**2-wv))
          tb=DSQRT(wv)
          np=0
          IF(tb.GT.0d0)THEN
             DO k=1,ffe_nbeta
                IF(ffe_beta(k).GT.d.AND.ffe_beta(k).LT.ffe_bcut)THEN
                   np=np+1
                   tk(np)=DSQRT((ffe_beta(k)-d)*(ffe_beta(k)+d))/tb
                ENDIF
             ENDDO
          ENDIF
          CALL ffe_breaks(np,tk,0d0,1d0,1d-12,n)
          DO k=1,n-1
             t0=tk(k)
             t1=tk(k+1)
             hh=0.5d0*(t1-t0)
             cc=0.5d0*(t1+t0)
             DO g=1,8
                tau=cc+hh*xg(g)
                b12=DSQRT(d*d+wv*tau*tau)
                om=ffe_Pfast(b12)
                IF(om.GT.0d0)om=om*ffe_Ptag(b12)
                pw=hh*wg(g)*om
                DO m=0,FFE_NMOM
                   ffe_ht(m,i,1)=ffe_ht(m,i,1)+pw
                   ffe_ht(m,i,2)=ffe_ht(m,i,2)+pw*b12*GeVm12fm
                   pw=pw*tau*tau
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    DEALLOCATE(tk)
  END SUBROUTINE ffe_build_pfast

  FUNCTION ffe_Pfast(b)
    ! P_had(b) (b in GeV^-1) from the fine table, 6-point Lagrange
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_Pfast
    REAL(KIND(1d0)),INTENT(IN)::b
    REAL(KIND(1d0))::u,x0,x1,x2,x3,x4,x5
    INTEGER::i0
    IF(b.GE.ffe_bcut)THEN
       ffe_Pfast=1d0
       RETURN
    ENDIF
    u=b/ffe_hf
    i0=MAX(0,MIN(INT(u)-2,ffe_nf-6))
    x0=u-DBLE(i0)
    x1=x0-1d0
    x2=x0-2d0
    x3=x0-3d0
    x4=x0-4d0
    x5=x0-5d0
    ffe_Pfast=-ffe_Pf(i0+1)*x1*x2*x3*x4*x5/120d0+ffe_Pf(i0+2)*x0*x2*x3*x4*x5/24d0&
         -ffe_Pf(i0+3)*x0*x1*x3*x4*x5/12d0+ffe_Pf(i0+4)*x0*x1*x2*x4*x5/12d0&
         -ffe_Pf(i0+5)*x0*x1*x2*x3*x5/24d0+ffe_Pf(i0+6)*x0*x1*x2*x3*x4/120d0
    ffe_Pfast=MAX(0d0,MIN(1d0,ffe_Pfast))
  END FUNCTION ffe_Pfast

  FUNCTION ffe_PhiD1(b1,b2)
    ! int_0^{2pi} dphi b12 = 4 (b1+b2) E(m), m = 4 b1 b2/(b1+b2)^2, in fm (b in GeV^-1)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_PhiD1
    REAL(KIND(1d0)),INTENT(IN)::b1,b2
    REAL(KIND(1d0)),EXTERNAL::elliptic_em
    REAL(KIND(1d0))::m
    m=4d0*b1*b2/(b1+b2)**2
    IF(m.GE.1d0-1d-15)THEN
       ffe_PhiD1=4d0*(b1+b2)*GeVm12fm
    ELSE
       ffe_PhiD1=4d0*(b1+b2)*elliptic_em(m)*GeVm12fm
    ENDIF
  END FUNCTION ffe_PhiD1

  SUBROUTINE ffe_phi(nk,b1,b2,d,phi)
    ! kernels Phi_k(b1,b2) = int_0^{2pi} dphi P(b12) w_k(b12): w_1 = 1, w_2 = b12 (fm) [<Delta b>],
    ! and Phi_3 = int dphi b12 (fm) [<Delta b> with PNOHAD=1]; b12^2 = d^2 + 4 b1 b2 sin^2(phi/2),
    ! d = b1 - b2 (GeV^-1). Without tagging P = 1 beyond bcut.
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nk
    REAL(KIND(1d0)),INTENT(IN)::b1,b2,d
    REAL(KIND(1d0)),DIMENSION(3),INTENT(OUT)::phi
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    INTEGER,PARAMETER::ng=4
    REAL(KIND(1d0)),DIMENSION(ng)::xg,wg
    REAL(KIND(1d0)),DIMENSION(8)::xt,wtt
    REAL(KIND(1d0)),DIMENSION(16)::xgl,wgl
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,SAVE::ph
    INTEGER::i,np,ig,m,i0
    REAL(KIND(1d0))::Gp,G4,ad,Dsum,hh,cc,b12,sn,pp,s1,s2,arg,pd1,y,cn,sum1,sum2,u,wv,php
    REAL(KIND(1d0)),DIMENSION(8)::wt
    LOGICAL::longtail
    ! phi breakpoints: at most one per kernel breakpoint plus the two ends
    IF(ALLOCATED(ph))THEN
       IF(SIZE(ph).LT.ffe_nbeta+16)DEALLOCATE(ph)
    ENDIF
    IF(.NOT.ALLOCATED(ph))ALLOCATE(ph(MAX(5000,ffe_nbeta+16)))
    phi=0d0
    Gp=b1*b2
    G4=4d0*Gp
    ad=DABS(d)
    Dsum=b1+b2
    pd1=0d0
    IF(nk.GE.2)pd1=ffe_PhiD1(b1,b2)
    IF(nk.GE.3)phi(3)=pd1
    longtail=ffe_tag
    IF(.NOT.longtail)THEN
       IF(ad.GE.ffe_bcut)THEN
          phi(1)=2d0*pi
          phi(2)=pd1
          RETURN
       ENDIF
       IF(Dsum.LE.ffe_beta(1))THEN
          ! P < 1e-14 everywhere on the circle
          RETURN
       ENDIF
       IF(Gp.GE.ffe_bcut**2)THEN
          ! far zone: expansion in tau^2/(4G) <= 1/4 with the moment tables (26 terms: < 1e-16)
          u=ad/ffe_hg
          CALL ffe_lagu(u,ffe_ng,8,i0,wt)
          y=1d0/G4
          cn=1d0
          sum1=0d0
          sum2=0d0
          DO m=0,FFE_NMOM
             sum1=sum1+cn*SUM(wt*ffe_gt(m,i0+1:i0+8,1))
             IF(nk.GE.2)sum2=sum2+cn*SUM(wt*ffe_gt(m,i0+1:i0+8,2))
             cn=cn*y*DBLE(2*m+1)/DBLE(2*m+2)
          ENDDO
          phi(1)=2d0*pi-2d0*sum1/DSQRT(Gp)
          phi(2)=pd1-2d0*sum2/DSQRT(Gp)
          RETURN
       ENDIF
    ENDIF
    CALL ffe_gauleg(8,xt,wtt)
    IF(longtail.AND.Gp.GE.ffe_bcut**2)THEN
       ! far zone with tagging: the arc b12 < bcut from the inner moments (expansion in
       ! tau^2/(4G) <= 1/4), the rest (P = P_tag) by 8-point panels between the tagging ladder
       sum1=0d0
       sum2=0d0
       php=0d0
       IF(ad.LT.ffe_bcut)THEN
          wv=(ffe_bcut-ad)*(ffe_bcut+ad)
          u=wv/ffe_hw
          CALL ffe_lagu(u,ffe_nw,8,i0,wt)
          y=wv/G4
          cn=1d0
          DO m=0,FFE_NMOM
             sum1=sum1+cn*SUM(wt*ffe_ht(m,i0+1:i0+8,1))
             IF(nk.GE.2)sum2=sum2+cn*SUM(wt*ffe_ht(m,i0+1:i0+8,2))
             cn=cn*y*DBLE(2*m+1)/DBLE(2*m+2)
          ENDDO
          sum1=sum1*DSQRT(wv)
          sum2=sum2*DSQRT(wv)
          php=2d0*DASIN(MIN(1d0,DSQRT(wv/G4)))
       ENDIF
       np=1
       ph(1)=php
       DO i=1,ffe_nbeta
          IF(ffe_beta(i).LE.ffe_bcut*(1d0+1d-12).OR.ffe_beta(i).LE.ad)CYCLE
          IF(ffe_beta(i).GE.Dsum)EXIT
          arg=(ffe_beta(i)-ad)*(ffe_beta(i)+ad)/G4
          IF(arg.GE.1d0)EXIT
          np=np+1
          ph(np)=2d0*DASIN(DSQRT(arg))
       ENDDO
       np=np+1
       ph(np)=pi
       s1=0d0
       s2=0d0
       DO i=1,np-1
          IF(ph(i+1).LE.ph(i))CYCLE
          hh=0.5d0*(ph(i+1)-ph(i))
          cc=0.5d0*(ph(i+1)+ph(i))
          DO ig=1,8
             sn=DSIN(0.5d0*(cc+hh*xt(ig)))
             b12=DSQRT(d*d+G4*sn*sn)
             pp=ffe_Ptag(b12)
             s1=s1+hh*wtt(ig)*pp
             IF(nk.GE.2)s2=s2+hh*wtt(ig)*pp*b12
          ENDDO
       ENDDO
       phi(1)=2d0*s1+2d0*sum1/DSQRT(Gp)
       IF(nk.GE.2)phi(2)=2d0*s2*GeVm12fm+2d0*sum2/DSQRT(Gp)
       RETURN
    ENDIF
    CALL ffe_gauleg(ng,xg,wg)
    CALL ffe_gauleg(16,xgl,wgl)
    np=1
    ph(1)=0d0
    DO i=1,ffe_nbeta
       IF(ffe_beta(i).LE.ad)CYCLE
       IF(ffe_beta(i).GE.Dsum)EXIT
       arg=(ffe_beta(i)-ad)*(ffe_beta(i)+ad)/G4
       IF(arg.GE.1d0)EXIT
       np=np+1
       ph(np)=2d0*DASIN(DSQRT(arg))
    ENDDO
    np=np+1
    ph(np)=pi
    s1=0d0
    s2=0d0
    DO i=1,np-1
       IF(ph(i+1).LE.ph(i))CYCLE
       hh=0.5d0*(ph(i+1)-ph(i))
       cc=0.5d0*(ph(i+1)+ph(i))
       IF(.NOT.longtail.AND.i.GT.1)THEN
          sn=DSIN(0.5d0*ph(i))
          IF(d*d+G4*sn*sn.GE.ffe_bcut**2*(1d0-1d-14))THEN
             ! beyond bcut: P = 1
             s1=s1+2d0*hh
             IF(nk.GE.2)THEN
                DO ig=1,16
                   sn=DSIN(0.5d0*(cc+hh*xgl(ig)))
                   s2=s2+hh*wgl(ig)*DSQRT(d*d+G4*sn*sn)
                ENDDO
             ENDIF
             CYCLE
          ENDIF
       ENDIF
       IF(longtail)THEN
          sn=DSIN(0.5d0*ph(i))
          IF(d*d+G4*sn*sn.GE.ffe_bcut**2*(1d0-1d-14))THEN
             ! tagging tail beyond bcut (panels between the geometric beta, ratio 2): 8 nodes,
             ! 4 leave -7e-7 in the XnXn flux
             DO ig=1,8
                sn=DSIN(0.5d0*(cc+hh*xt(ig)))
                b12=DSQRT(d*d+G4*sn*sn)
                pp=ffe_Ptag(b12)
                s1=s1+hh*wtt(ig)*pp
                IF(nk.GE.2)s2=s2+hh*wtt(ig)*pp*b12
             ENDDO
             CYCLE
          ENDIF
       ENDIF
       DO ig=1,ng
          sn=DSIN(0.5d0*(cc+hh*xg(ig)))
          b12=DSQRT(d*d+G4*sn*sn)
          pp=ffe_Pfast(b12)
          IF(ffe_tag.AND.pp.GT.0d0)pp=pp*ffe_Ptag(b12)
          s1=s1+hh*wg(ig)*pp
          IF(nk.GE.2)s2=s2+hh*wg(ig)*pp*b12
       ENDDO
    ENDDO
    phi(1)=2d0*s1
    IF(nk.GE.2)phi(2)=2d0*s2*GeVm12fm
  END SUBROUTINE ffe_phi
  SUBROUTINE ffe_tnodes(s,ndb,far,nt,tt,wt)
    ! composite Gauss-Legendre rule in t = (v1-v2)/2 at fixed s = (v1+v2)/2: breakpoints at the
    ! domain edges (plus graded lines near the lower edges), the v-grid segment boundaries, the
    ! kernel lines |b1-b2| = beta_k and b1+b2 = beta_k, and (ndb) graded lines around t = 0
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::s
    LOGICAL,INTENT(IN)::ndb,far
    INTEGER,INTENT(OUT)::nt
    REAL(KIND(1d0)),DIMENSION(:),INTENT(OUT)::tt,wt
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,SAVE::br
    REAL(KIND(1d0)),DIMENSION(10),PARAMETER::del=(/0.002d0,0.005d0,0.01d0,0.02d0,0.05d0,0.1d0,&
         0.2d0,0.3d0,0.5d0,1d0/)
    REAL(KIND(1d0)),DIMENSION(8)::x8,w8
    REAL(KIND(1d0)),DIMENSION(5)::x5,w5
    INTEGER::nb,k,i,np,ns,j,g,nj,nbmax
    REAL(KIND(1d0))::V1lo,V1hi,V2lo,V2hi,tlo,thi,es,x,len,h,a,ht,u,tj(4),ta,tb
    LOGICAL::lsing,rsing
    ! t breakpoints: 1 + 20 graded lines + segment boundaries + 4 per kernel breakpoint
    ! + 80 around t = 0 + 4 jump lines, and the two ends added by ffe_breaks
    nbmax=27+ffe_vg(1)%nseg+ffe_vg(2)%nseg+4*ffe_nbeta+80
    IF(ALLOCATED(br))THEN
       IF(SIZE(br).LT.nbmax)DEALLOCATE(br)
    ENDIF
    IF(.NOT.ALLOCATED(br))ALLOCATE(br(MAX(40000,nbmax)))
    V1lo=ffe_vg(1)%v(1)
    V1hi=ffe_vg(1)%v(ffe_vg(1)%n)
    V2lo=ffe_vg(2)%v(1)
    V2hi=ffe_vg(2)%v(ffe_vg(2)%n)
    tlo=MAX(V1lo-s,s-V2hi)
    thi=MIN(V1hi-s,s-V2lo)
    nt=0
    IF(thi-tlo.LE.1d-12)RETURN
    es=DEXP(s)
    nb=1
    br(1)=0d0
    DO k=1,10
       br(nb+1)=V1lo+del(k)-s
       br(nb+2)=s-(V2lo+del(k))
       nb=nb+2
    ENDDO
    DO k=2,ffe_vg(1)%nseg
       nb=nb+1
       br(nb)=ffe_vg(1)%sv0(k)-s
    ENDDO
    DO k=2,ffe_vg(2)%nseg
       nb=nb+1
       br(nb)=s-ffe_vg(2)%sv0(k)
    ENDDO
    DO i=1,ffe_nbeta
       IF(far.AND.MOD(i,4).NE.0.AND.i.NE.ffe_nbeta)CYCLE
       x=ffe_beta(i)/(2d0*es)
       a=ASINH(x)
       IF(a.GT.thi-tlo)CYCLE
       br(nb+1)=a
       br(nb+2)=-a
       nb=nb+2
       IF(x.GT.1d0)THEN
          a=ACOSH(x)
          br(nb+1)=a
          br(nb+2)=-a
          nb=nb+2
       ENDIF
    ENDDO
    IF(ndb)THEN
       x=0.5d0
       DO k=1,40
          x=0.5d0*x
          IF(x.LT.1d-9)EXIT
          br(nb+1)=x
          br(nb+2)=-x
          nb=nb+2
       ENDDO
    ENDIF
    ! lines where P jumps (closed-form fits with a hard cut-off): Phi ~ sqrt(distance) there
    nj=0
    IF(ffe_pjump)THEN
       x=ffe_bcut/(2d0*es)
       tj(1)=ASINH(x)
       tj(2)=-tj(1)
       nj=2
       IF(x.GT.1d0)THEN
          tj(3)=ACOSH(x)
          tj(4)=-tj(3)
          nj=4
       ENDIF
       DO k=1,nj
          nb=nb+1
          br(nb)=tj(k)
       ENDDO
    ENDIF
    CALL ffe_breaks(nb,br,tlo,thi,1d-13,np)
    CALL ffe_gauleg(8,x8,w8)
    CALL ffe_gauleg(5,x5,w5)
    ht=0.25d0/MAX(FFE_NodeScale,1d-3)
    ! the panels below add at most 24 nodes each plus 8 per length ht
    IF(24*(np-1)+8*CEILING((thi-tlo)/ht)+8.GT.SIZE(tt))THEN
       WRITE(*,*)"ERROR: ffe_tnodes: up to",24*(np-1)+8*CEILING((thi-tlo)/ht)+8,&
            " t nodes for a work array of",SIZE(tt)
       STOP 1
    ENDIF
    DO i=1,np-1
       len=br(i+1)-br(i)
       ns=MAX(1,CEILING(len/ht-1d-9))
       h=len/DBLE(ns)
       lsing=.FALSE.
       rsing=.FALSE.
       DO k=1,nj
          IF(DABS(br(i)-tj(k)).LT.1d-12)lsing=.TRUE.
          IF(DABS(br(i+1)-tj(k)).LT.1d-12)rsing=.TRUE.
       ENDDO
       IF((lsing.OR.rsing).AND.len.GT.2d0*ht)THEN
          ! a long panel at a jump line: u^2 substitution on the sub-panel of length ht next to
          ! each singular end (8 nodes), regular panels on the rest (the far side of a jump line
          ! is smooth, and the integrand varies with W1 W2 on the scale of the whole panel)
          ta=br(i)
          tb=br(i+1)
          IF(lsing)THEN
             DO g=1,8
                u=0.5d0*(1d0+x8(g))
                nt=nt+1
                tt(nt)=ta+ht*u*u
                wt(nt)=0.5d0*w8(g)*2d0*ht*u
             ENDDO
             ta=ta+ht
          ENDIF
          IF(rsing)THEN
             DO g=1,8
                u=0.5d0*(1d0+x8(g))
                nt=nt+1
                tt(nt)=tb-ht*u*u
                wt(nt)=0.5d0*w8(g)*2d0*ht*u
             ENDDO
             tb=tb-ht
          ENDIF
          ns=MAX(1,CEILING((tb-ta)/ht-1d-9))
          h=(tb-ta)/DBLE(ns)
          DO j=1,ns
             a=ta+(DBLE(j)-0.5d0)*h
             DO g=1,8
                nt=nt+1
                tt(nt)=a+0.5d0*h*x8(g)
                wt(nt)=0.5d0*h*w8(g)
             ENDDO
          ENDDO
          CYCLE
       ENDIF
       IF(lsing.OR.rsing)THEN
          ! u^2 substitution towards the singular end(s), 16 nodes
          CALL ffe_gauleg(8,x8,w8)
          ta=br(i)
          tb=br(i+1)
          IF(lsing.AND.rsing)THEN
             ! both ends: split in the middle
             DO g=1,8
                u=0.5d0*(1d0+x8(g))
                nt=nt+1
                tt(nt)=ta+0.5d0*len*u*u
                wt(nt)=0.5d0*w8(g)*len*u*0.5d0*2d0
                nt=nt+1
                tt(nt)=tb-0.5d0*len*u*u
                wt(nt)=0.5d0*w8(g)*len*u*0.5d0*2d0
             ENDDO
          ELSE
             DO g=1,16/2
                u=0.5d0*(1d0+x8(g))
                nt=nt+1
                IF(lsing)THEN
                   tt(nt)=ta+len*u*u
                ELSE
                   tt(nt)=tb-len*u*u
                ENDIF
                wt(nt)=0.5d0*w8(g)*2d0*len*u
             ENDDO
          ENDIF
          CYCLE
       ENDIF
       DO j=1,ns
          a=br(i)+(DBLE(j)-0.5d0)*h
          IF(h.GT.0.02d0.OR.FFE_NodeScale.GT.1.5d0)THEN
             DO g=1,8
                nt=nt+1
                tt(nt)=a+0.5d0*h*x8(g)
                wt(nt)=0.5d0*h*w8(g)
             ENDDO
          ELSE
             DO g=1,5
                nt=nt+1
                tt(nt)=a+0.5d0*h*x5(g)
                wt(nt)=0.5d0*h*w5(g)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE ffe_tnodes

  SUBROUTINE ffe_spoint(s,ws,nk,mode,npx,ipx,prox,off,cnt)
    ! all inner nodes at fixed s: kernel values and proxy integrands; mode 1 adds ws*c*L1(x)L2 to
    ! the K matrices directly (far zone), mode 0 stores (v1, v2, c) in the node pool (near zone)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::s,ws
    INTEGER,INTENT(IN)::nk,npx,mode
    INTEGER,DIMENSION(2,npx),INTENT(IN)::ipx
    REAL(KIND(1d0)),DIMENSION(npx),INTENT(OUT)::prox
    INTEGER,INTENT(OUT)::off,cnt
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,SAVE::tt,wt
    REAL(KIND(1d0)),DIMENSION(p)::w1,w2
    REAL(KIND(1d0)),DIMENSION(3)::phi,c
    INTEGER::nt,i,i1,i2,m
    REAL(KIND(1d0))::es,t,v1,v2,b1,b2,d
    IF(.NOT.ALLOCATED(tt))ALLOCATE(tt(200000),wt(200000))
    prox=0d0
    CALL ffe_tnodes(s,nk.GE.2,mode.EQ.1,nt,tt,wt)
    es=DEXP(s)
    off=ffe_npool
    cnt=0
    IF(mode.EQ.0)CALL ffe_pool_reserve(nt)
    DO i=1,nt
       t=tt(i)
       v1=s+t
       v2=s-t
       b1=DEXP(v1)
       b2=DEXP(v2)
       d=2d0*es*DSINH(t)
       CALL ffe_phi(nk,b1,b2,d,phi)
       c(1:nk)=2d0*wt(i)*phi(1:nk)
       IF(MAXVAL(DABS(c(1:nk))).EQ.0d0)CYCLE
       CALL ffe_lagw(ffe_vg(1),v1,p,i1,w1)
       CALL ffe_lagw(ffe_vg(2),v2,p,i2,w2)
       IF(mode.EQ.1)THEN
          CALL ffe_scatter(nk,ws*c,i1,w1,i2,w2)
       ELSE
          ffe_npool=ffe_npool+1
          cnt=cnt+1
          ffe_pool(1,ffe_npool)=v1
          ffe_pool(2,ffe_npool)=v2
          ffe_pool(3:2+nk,ffe_npool)=c(1:nk)
       ENDIF
       DO m=1,npx
          prox(m)=prox(m)+c(1)*SUM(w1*ffe_W(ipx(1,m),i1+1:i1+p,1))*SUM(w2*ffe_W(ipx(2,m),i2+1:i2+p,2))
       ENDDO
    ENDDO
  END SUBROUTINE ffe_spoint

  SUBROUTINE ffe_scatter(nk,c,i1,w1,i2,w2)
    ! K(i1+a, i2+b, k) += c(k) w1(a) w2(b)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nk,i1,i2
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::c
    REAL(KIND(1d0)),DIMENSION(8),INTENT(IN)::w1,w2
    INTEGER::k,b
    REAL(KIND(1d0))::f
    DO k=1,nk
       IF(c(k).EQ.0d0)CYCLE
       DO b=1,8
          f=c(k)*w2(b)
          ffe_K(i1+1:i1+8,i2+b,k)=ffe_K(i1+1:i1+8,i2+b,k)+f*w1
       ENDDO
    ENDDO
    ffe_nodes=ffe_nodes+1
  END SUBROUTINE ffe_scatter

  SUBROUTINE ffe_pool_reserve(n)
    ! room for n more pooled nodes (v1, v2, c1..c3)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::tmp
    INTEGER::nn
    IF(.NOT.ALLOCATED(ffe_pool))THEN
       ALLOCATE(ffe_pool(5,MAX(n,500000)))
       ffe_npool=0
    ENDIF
    IF(ffe_npool+n.GT.SIZE(ffe_pool,2))THEN
       nn=MAX(2*SIZE(ffe_pool,2),ffe_npool+n)
       ALLOCATE(tmp(5,nn))
       tmp(:,1:ffe_npool)=ffe_pool(:,1:ffe_npool)
       CALL MOVE_ALLOC(tmp,ffe_pool)
    ENDIF
  END SUBROUTINE ffe_pool_reserve

  SUBROUTINE ffe_build_K(nk)
    ! K matrices of kinds 1..nk: fixed Gauss-Legendre panels in s in the far zone, global adaptive
    ! Gauss-Kronrod 10/21 in s in the near zone, driven by x-independent proxy integrals; the
    ! near-zone nodes are pooled during the adaptation and scattered once at the end
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nk
    INTEGER,PARAMETER::npmax=12,nimax=4000
    REAL(KIND(1d0)),DIMENSION(11),PARAMETER::xgk=(/0.995657163025808080735527280689003d0,&
         0.973906528517171720077964012084452d0,0.930157491355708226001207180059508d0,&
         0.865063366688984510732096688423493d0,0.780817726586416897063717578345042d0,&
         0.679409568299024406234327365114874d0,0.562757134668604683339000099272694d0,&
         0.433395394129247190799265943165784d0,0.294392862701460198131126603103866d0,&
         0.148874338981631210884826001129720d0,0d0/)
    REAL(KIND(1d0)),DIMENSION(11),PARAMETER::wgk=(/0.011694638867371874278064396062192d0,&
         0.032558162307964727478818972459390d0,0.054755896574351996031381300244580d0,&
         0.075039674810919952767043140916190d0,0.093125454583697605535065465083366d0,&
         0.109387158802297641899210590325805d0,0.123491976262065851077958109831074d0,&
         0.134709217311473325928054001771707d0,0.142775938577060080797094273138717d0,&
         0.147739104901338491374841515972068d0,0.149445554002916905664936468389821d0/)
    REAL(KIND(1d0)),DIMENSION(5),PARAMETER::wg10=(/0.066671344308688137593568809893332d0,&
         0.149451349150580593145776339657697d0,0.219086362515982043995534934228163d0,&
         0.269266719309996355091226921569469d0,0.295524224714752870173892994651338d0/)
    INTEGER,DIMENSION(2,npmax)::ipx
    REAL(KIND(1d0)),DIMENSION(npmax)::prox,Tf,Tot
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::ia,ibb
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::res,err
    INTEGER,DIMENSION(:,:),ALLOCATABLE::poff,pcnt
    REAL(KIND(1d0)),DIMENSION(8)::x8,w8,w1,w2
    REAL(KIND(1d0)),DIMENSION(3)::cc3
    INTEGER::npx,n1,n2,i,g,nint,iw,it,m,jj,off,cnt,n,i1,i2
    REAL(KIND(1d0))::slo,shi,sad,hs,a,c,h,tol,errsum,worst,e,hl,ws,g1,g2,tk0,tk1
    REAL(KIND(1d0)),DIMENSION(6)::xl
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::Wsv
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Wr
    LOGICAL::ourpx
    INTEGER::ib
    CALL ffe_build_tables()
    n1=ffe_vg(1)%n
    n2=ffe_vg(2)%n
    IF(.NOT.ALLOCATED(ffe_K))THEN
       ALLOCATE(ffe_K(n1,n2,3))
       ffe_K=0d0
    ENDIF
    ffe_K(:,:,1:nk)=0d0
    ffe_nodes=0
    ffe_npool=0
    ! proxies: pairs of grid x values (indices into ffe_W) inside the physically relevant range
    ! x_i <= g_i = min(0.3, 3/(m_i R_i)) (photon energies up to 3 gamma/R)
    g1=MIN(0.3d0,3d0/(ffe_bm(1)%mass*ffe_bm(1)%R))
    g2=MIN(0.3d0,3d0/(ffe_bm(2)%mass*ffe_bm(2)%R))
    xl=(/1d-6,1d-4,1d-2,0d0,0d0,0d0/)
    npx=0
    DO i=1,5
       npx=npx+1
       IF(i.LE.3)THEN
          ipx(:,npx)=(/ffe_xindex(MIN(xl(i),g1)),ffe_xindex(MIN(xl(i),g2))/)
       ELSEIF(i.EQ.4)THEN
          ipx(:,npx)=(/ffe_xindex(g1/3d0),ffe_xindex(g2/3d0)/)
       ELSE
          ipx(:,npx)=(/ffe_xindex(g1),ffe_xindex(g2)/)
       ENDIF
    ENDDO
    npx=npx+1
    ipx(:,npx)=(/ffe_xindex(1d-6),ffe_xindex(MIN(1d-2,g2))/)
    npx=npx+1
    ipx(:,npx)=(/ffe_xindex(MIN(1d-2,g1)),ffe_xindex(1d-6)/)
    npx=npx+1
    ipx(:,npx)=(/ffe_xindex(1d-4),ffe_xindex(g2)/)
    npx=npx+1
    ipx(:,npx)=(/ffe_xindex(g1),ffe_xindex(1d-4)/)
    ! with the proton W rows of the inelastic branch (ffe_ourW: kinks of the interpolated correction grid, the cut of
    ! the correction at 10 R_p), the proxies use the closed-form dipole rows of the proton: K does not depend on W,
    ! and smooth proxies steer the adaptation as for the dipole (the rows of ffe_W are restored at the end)
    ourpx=ffe_ourW()
    IF(ourpx)THEN
       ALLOCATE(Wsv(npx,SIZE(ffe_W,2),2))
       DO ib=1,2
          DO m=1,npx
             Wsv(m,:,ib)=ffe_W(ipx(ib,m),:,ib)
          ENDDO
          IF(.NOT.ffe_bm(ib)%proton)CYCLE
          ALLOCATE(Wr(ffe_vg(ib)%n))
          DO m=1,npx
             CALL ffe_W_row(ib,10d0**(-ffe_xg%v(ipx(ib,m))),Wr,.TRUE.)
             ffe_W(ipx(ib,m),1:ffe_vg(ib)%n,ib)=Wr
          ENDDO
          DEALLOCATE(Wr)
       ENDDO
    ENDIF
    slo=0.5d0*(ffe_vg(1)%v(1)+ffe_vg(2)%v(1))
    shi=0.5d0*(ffe_vg(1)%v(n1)+ffe_vg(2)%v(n2))
    sad=MAX(slo+0.2d0,MIN(shi,DLOG(2d0*ffe_bcut)+0.5d0))
    ! far zone: fixed panels, scattered directly
    CALL CPU_TIME(tk0)
    CALL ffe_gauleg(8,x8,w8)
    hs=0.25d0/MAX(FFE_NodeScale,1d-3)
    Tf=0d0
    IF(shi.GT.sad)THEN
       nint=MAX(1,CEILING((shi-sad)/hs))
       h=(shi-sad)/DBLE(nint)
       DO i=1,nint
          a=sad+DBLE(i-1)*h
          DO g=1,8
             CALL ffe_spoint(a+0.5d0*h*(1d0+x8(g)),0.5d0*h*w8(g),nk,1,npx,ipx(:,1:npx),prox(1:npx),off,cnt)
             Tf(1:npx)=Tf(1:npx)+0.5d0*h*w8(g)*prox(1:npx)
          ENDDO
       ENDDO
    ENDIF
    CALL CPU_TIME(tk1)
    IF(FFE_Verbose.GE.2)WRITE(*,*)"FFE: far zone CPU [s] =",tk1-tk0,", nodes =",ffe_nodes
    tk0=tk1
    ! near zone: global adaptive Gauss-Kronrod 10/21 with pooled nodes
    ALLOCATE(ia(nimax),ibb(nimax),res(npmax,nimax),err(npmax,nimax),poff(21,nimax),pcnt(21,nimax))
    nint=MAX(4,CEILING((sad-slo)/0.2d0))
    h=(sad-slo)/DBLE(nint)
    DO i=1,nint
       ia(i)=slo+DBLE(i-1)*h
       ibb(i)=ia(i)+h
    ENDDO
    IF(ffe_pjump)THEN
       ! a jump of P at bcut: the line b1+b2 = bcut touches t = 0 at s = ln(bcut/2), where the
       ! inner integral has a square-root singularity -> make it an interval end
       c=DLOG(0.5d0*ffe_bcut)
       DO i=1,nint
          IF(c.GT.ia(i)+1d-9.AND.c.LT.ibb(i)-1d-9)THEN
             nint=nint+1
             ia(nint)=c
             ibb(nint)=ibb(i)
             ibb(i)=c
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,nint
       CALL ffe_gk21(i)
    ENDDO
    tol=1d-9/MAX(FFE_NodeScale,1d-3)**2
    ! kernels with a jump (closed-form MC-Glauber fits with a hard cut-off) converge only
    ! algebraically near the jump lines: 1e-7 is still two orders below the 1e-5 target
    IF(ffe_pjump)tol=1d-7/MAX(FFE_NodeScale,1d-3)**2
    DO it=1,nimax
       Tot(1:npx)=Tf(1:npx)+SUM(res(1:npx,1:nint),DIM=2)
       errsum=0d0
       worst=-1d0
       iw=1
       DO i=1,nint
          e=0d0
          DO m=1,npx
             IF(Tot(m).GT.0d0)e=MAX(e,err(m,i)/Tot(m))
          ENDDO
          errsum=errsum+e
          IF(e.GT.worst)THEN
             worst=e
             iw=i
          ENDIF
       ENDDO
       IF(errsum.LE.tol.OR.nint.GE.nimax-1)EXIT
       IF(FFE_Verbose.GE.2.AND.MOD(it,25).EQ.1)THEN
          WRITE(*,'(A,I6,A,ES10.3,A,2F10.5,A,ES10.3)')" FFE: adapt nint=",nint," errsum=",errsum,&
               " worst interval=",ia(iw),ibb(iw)," worst=",worst
          WRITE(*,'(A,12ES10.2)')"      per-proxy rel err:",(SUM(err(m,1:nint))/MAX(Tot(m),1d-300),m=1,npx)
       ENDIF
       c=0.5d0*(ia(iw)+ibb(iw))
       nint=nint+1
       ia(nint)=c
       ibb(nint)=ibb(iw)
       ibb(iw)=c
       CALL ffe_gk21(iw)
       CALL ffe_gk21(nint)
    ENDDO
    IF(FFE_Verbose.GE.1)WRITE(*,*)"FFE: near-zone s intervals =",nint,", estimated rel. error =",errsum,&
         ", pooled nodes =",ffe_npool
    CALL CPU_TIME(tk1)
    IF(FFE_Verbose.GE.2)WRITE(*,*)"FFE: near zone adaptation CPU [s] =",tk1-tk0
    tk0=tk1
    ! scatter the pooled nodes of the final near-zone intervals
    DO i=1,nint
       hl=0.5d0*(ibb(i)-ia(i))
       DO jj=1,21
          ws=hl*wgk(MIN(jj,22-jj))
          off=poff(jj,i)
          cnt=pcnt(jj,i)
          DO n=off+1,off+cnt
             CALL ffe_lagw(ffe_vg(1),ffe_pool(1,n),8,i1,w1)
             CALL ffe_lagw(ffe_vg(2),ffe_pool(2,n),8,i2,w2)
             cc3=0d0
             cc3(1:nk)=ws*ffe_pool(3:2+nk,n)
             CALL ffe_scatter(nk,cc3,i1,w1,i2,w2)
          ENDDO
       ENDDO
    ENDDO
    CALL CPU_TIME(tk1)
    IF(FFE_Verbose.GE.2)WRITE(*,*)"FFE: near zone scatter CPU [s] =",tk1-tk0,", nodes =",ffe_nodes
    ffe_Kready(1:nk)=.TRUE.
    IF(ourpx)THEN
       ! the rows of the branch back (in reverse order: an x index may serve several proxies)
       DO ib=1,2
          DO m=npx,1,-1
             ffe_W(ipx(ib,m),:,ib)=Wsv(m,:,ib)
          ENDDO
       ENDDO
       DEALLOCATE(Wsv)
    ENDIF
    DEALLOCATE(ia,ibb,res,err,poff,pcnt)
    IF(ALLOCATED(ffe_pool))DEALLOCATE(ffe_pool)
    ffe_npool=0
  CONTAINS
    SUBROUTINE ffe_gk21(ii)
      ! Gauss-Kronrod 10/21 on interval ii; the 21 points are numbered 1..11 (left, centre last)
      ! as cc - hl xgk(j) and 12..21 as cc + hl xgk(j-11)... stored in poff/pcnt(1:21,ii)
      INTEGER,INTENT(IN)::ii
      REAL(KIND(1d0)),DIMENSION(npmax)::rg,pp
      REAL(KIND(1d0))::hl0,cc0
      INTEGER::j,o,ct
      hl0=0.5d0*(ibb(ii)-ia(ii))
      cc0=0.5d0*(ibb(ii)+ia(ii))
      res(1:npx,ii)=0d0
      rg(1:npx)=0d0
      DO j=1,11
         CALL ffe_spoint(cc0-hl0*xgk(j),0d0,nk,0,npx,ipx(:,1:npx),pp(1:npx),o,ct)
         poff(j,ii)=o
         pcnt(j,ii)=ct
         res(1:npx,ii)=res(1:npx,ii)+hl0*wgk(j)*pp(1:npx)
         IF(MOD(j,2).EQ.0)rg(1:npx)=rg(1:npx)+hl0*wg10(j/2)*pp(1:npx)
         IF(j.LT.11)THEN
            CALL ffe_spoint(cc0+hl0*xgk(j),0d0,nk,0,npx,ipx(:,1:npx),pp(1:npx),o,ct)
            poff(22-j,ii)=o
            pcnt(22-j,ii)=ct
            res(1:npx,ii)=res(1:npx,ii)+hl0*wgk(j)*pp(1:npx)
            IF(MOD(j,2).EQ.0)rg(1:npx)=rg(1:npx)+hl0*wg10(j/2)*pp(1:npx)
         ENDIF
      ENDDO
      err(1:npx,ii)=DABS(res(1:npx,ii)-rg(1:npx))
    END SUBROUTINE ffe_gk21
  END SUBROUTINE ffe_build_K

  FUNCTION ffe_xindex(x)
    ! grid node (1 = x=1, MX+1 = xmin) closest to x
    IMPLICIT NONE
    INTEGER::ffe_xindex
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::lo,hi,mid
    REAL(KIND(1d0))::u
    u=-DLOG10(x)
    lo=1
    hi=ffe_MX+1
    DO WHILE(hi-lo.GT.1)
       mid=(lo+hi)/2
       IF(ffe_xg%v(mid).LE.u)THEN
          lo=mid
       ELSE
          hi=mid
       ENDIF
    ENDDO
    ffe_xindex=lo
    IF(DABS(ffe_xg%v(hi)-u).LT.DABS(ffe_xg%v(lo)-u))ffe_xindex=hi
  END FUNCTION ffe_xindex

  SUBROUTINE ffe_build_grid(k)
    ! ZL(i,j,k) = ln(x1 x2 Z) on the grid: Z = 2 pi alpha^2 Z1^2 Z2^2/(x1 x2) W1 K W2^T
    IMPLICIT NONE
    INTEGER,INTENT(IN)::k
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::T,S
    REAL(KIND(1d0))::pref,u
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::an,xv
    INTEGER::i,j,kk
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    IF(.NOT.ffe_Kready(k))CALL ffe_build_K(k)
    ALLOCATE(T(ffe_MX+1,ffe_vg(2)%n),S(ffe_MX+1,ffe_MX+1))
    T=MATMUL(ffe_W(:,1:ffe_vg(1)%n,1),ffe_K(:,:,k))
    S=MATMUL(T,TRANSPOSE(ffe_W(:,1:ffe_vg(2)%n,2)))
    pref=twopi*ffe_alpha**2*ffe_bm(1)%znum**2*ffe_bm(2)%znum**2
    IF(.NOT.ALLOCATED(ffe_ZL))ALLOCATE(ffe_ZL(ffe_MX+1,ffe_MX+1,3))
    IF(alpha_running4photon.AND.ffe_icoll.NE.3)THEN
       ! ffe_anorm2 depends on x1 x2 only: the grid has log10 x on multiples of 1/80, so u_i + u_j = kk/80 repeats;
       ! the value is kept per kk (pairs with the same kk differ in x1 x2 by rounding only, ~1e-16)
       ALLOCATE(an(0:NINT(160d0*ffe_xg%v(ffe_MX+1))+1),xv(ffe_MX+1))
       an=-1d0
       xv=10d0**(-ffe_xg%v(1:ffe_MX+1))
       DO j=1,ffe_MX+1
          DO i=1,ffe_MX+1
             u=80d0*(ffe_xg%v(i)+ffe_xg%v(j))
             kk=NINT(u)
             IF(ABS(u-kk).LT.1d-9.AND.kk.GE.0.AND.kk.LE.UBOUND(an,1))THEN
                IF(an(kk).LT.0d0)an(kk)=ffe_anorm2(xv(i),xv(j))
                S(i,j)=S(i,j)*an(kk)
             ELSE
                S(i,j)=S(i,j)*ffe_anorm2(xv(i),xv(j))
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(an,xv)
    ENDIF
    DO j=1,ffe_MX+1
       DO i=1,ffe_MX+1
          IF(S(i,j).GT.0d0)THEN
             ffe_ZL(i,j,k)=DLOG(pref*S(i,j))
          ELSE
             ffe_ZL(i,j,k)=-745d0
          ENDIF
       ENDDO
    ENDDO
    ffe_Zready(k)=.TRUE.
    DEALLOCATE(T,S)
  END SUBROUTINE ffe_build_grid

  FUNCTION ffe_grid_val(k,x1,x2,zero)
    ! grid lookup: 8x8 Lagrange in log10 x of ln(x1 x2 Z); zero (optional): a node of the stencil has Z = 0
    ! (ln = -745: the exact elastic photon vanishes at x = 1)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_grid_val
    INTEGER,INTENT(IN)::k
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(OUT),OPTIONAL::zero
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(p)::w1,w2
    INTEGER::i0,j0,a
    REAL(KIND(1d0))::L
    CALL ffe_lagw(ffe_xg,-DLOG10(x1),p,i0,w1)
    CALL ffe_lagw(ffe_xg,-DLOG10(x2),p,j0,w2)
    IF(PRESENT(zero))zero=ANY(ffe_ZL(i0+1:i0+p,j0+1:j0+p,k).LE.-700d0)
    L=0d0
    DO a=1,p
       L=L+w1(a)*SUM(w2*ffe_ZL(i0+a,j0+1:j0+p,k))
    ENDDO
    ffe_grid_val=DEXP(L)/(x1*x2)
  END FUNCTION ffe_grid_val

  FUNCTION ffe_direct(k,x1,x2)
    ! direct evaluation W1(x1)^T K W2(x2) (no x interpolation)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_direct
    INTEGER,INTENT(IN)::k
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::W1,W2
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    IF(.NOT.ffe_Kready(k))CALL ffe_build_K(k)
    ALLOCATE(W1(ffe_vg(1)%n),W2(ffe_vg(2)%n))
    CALL ffe_W_row(1,x1,W1)
    CALL ffe_W_row(2,x2,W2)
    ffe_direct=twopi*ffe_alpha**2*ffe_bm(1)%znum**2*ffe_bm(2)%znum**2*DOT_PRODUCT(W1,MATMUL(ffe_K(:,:,k),W2))/(x1*x2)
    IF(alpha_running4photon.AND.ffe_icoll.NE.3)ffe_direct=ffe_direct*ffe_anorm2(x1,x2)
    DEALLOCATE(W1,W2)
  END FUNCTION ffe_direct

  FUNCTION ffe_pnohad1(x1,x2,direct)
    ! (2 pi)^2 alpha^2 Z1^2 Z2^2/(x1 x2) Nt1 Nt2 with Nt = int W dv (PNOHAD = 1, no b truncation)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_pnohad1
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN)::direct
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(p)::wt
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::W
    REAL(KIND(1d0))::n1,n2
    INTEGER::i0
    LOGICAL::dir
    REAL(KIND(1d0)),PARAMETER::fourpi2=39.4784176043574344753379639995d0
    dir=direct
    IF(.NOT.dir.AND.ffe_ourW())THEN
       ! as ffe_flux: the proton photon of the inelastic branch vanishes at x = 1 (patch, else direct)
       n1=ffe_N_val(1,x1)
       n2=ffe_N_val(2,x2)
       ffe_pnohad1=fourpi2*ffe_alpha**2*ffe_bm(1)%znum**2*ffe_bm(2)%znum**2*n1*n2/(x1*x2)
       IF(alpha_running4photon.AND.ffe_icoll.NE.3)ffe_pnohad1=ffe_pnohad1*ffe_anorm2(x1,x2)
       RETURN
    ENDIF
    IF(dir)THEN
       ALLOCATE(W(ffe_vg(1)%n))
       CALL ffe_W_row(1,x1,W)
       n1=DOT_PRODUCT(ffe_vg(1)%qw,W)
       DEALLOCATE(W)
       ALLOCATE(W(ffe_vg(2)%n))
       CALL ffe_W_row(2,x2,W)
       n2=DOT_PRODUCT(ffe_vg(2)%qw,W)
       DEALLOCATE(W)
    ELSE
       CALL ffe_lagw(ffe_xg,-DLOG10(x1),p,i0,wt)
       n1=DEXP(SUM(wt*ffe_lnN(i0+1:i0+p,1)))
       CALL ffe_lagw(ffe_xg,-DLOG10(x2),p,i0,wt)
       n2=DEXP(SUM(wt*ffe_lnN(i0+1:i0+p,2)))
    ENDIF
    ffe_pnohad1=fourpi2*ffe_alpha**2*ffe_bm(1)%znum**2*ffe_bm(2)%znum**2*n1*n2/(x1*x2)
    IF(alpha_running4photon.AND.ffe_icoll.NE.3)ffe_pnohad1=ffe_pnohad1*ffe_anorm2(x1,x2)
  END FUNCTION ffe_pnohad1
  SUBROUTINE ffe_init(icoll,xreq)
    ! (re)build every x-independent table of the engine when the configuration changed
    IMPLICIT NONE
    include 'run90.inc'
    INTEGER,INTENT(IN)::icoll
    REAL(KIND(1d0)),INTENT(IN)::xreq
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0
    REAL(KIND(1d0)),DIMENSION(40)::sig
    REAL(KIND(1d0))::A,Z,R,aa,w,cmenergy,t0,t1,tt0,base
    CHARACTER(len=7)::aname
    INTEGER::ib,inuc,lx
    ! MadGraph: the beams come from the MadGraph commons of run90.inc (ebeamMG5 is the energy
    ! of the whole beam particle), converted as in the legacy functions
    IF(nb_proton(1).EQ.1.AND.nb_neutron(1).EQ.0)THEN
       nuclearA_beam1=0
       nuclearZ_beam1=0
    ELSE
       nuclearA_beam1=nb_proton(1)+nb_neutron(1)
       nuclearZ_beam1=nb_proton(1)
    ENDIF
    IF(nb_proton(2).EQ.1.AND.nb_neutron(2).EQ.0)THEN
       nuclearA_beam2=0
       nuclearZ_beam2=0
    ELSE
       nuclearA_beam2=nb_proton(2)+nb_neutron(2)
       nuclearZ_beam2=nb_proton(2)
    ENDIF
    ebeam_PN(1)=ebeamMG5(1)/(nb_proton(1)+nb_neutron(1))
    ebeam_PN(2)=ebeamMG5(2)/(nb_proton(2)+nb_neutron(2))
    IF(xreq.LT.ffe_xmin)THEN
       lx=FLOOR(DLOG10(xreq))
       DO WHILE(10d0**lx.GT.xreq)
          lx=lx-1
       ENDDO
       ffe_xmin=10d0**lx
    ENDIF
    sig=0d0
    sig(1)=DBLE(icoll)
    sig(2)=DBLE(nuclearA_beam1)
    sig(3)=DBLE(nuclearZ_beam1)
    sig(4)=DBLE(nuclearA_beam2)
    sig(5)=DBLE(nuclearZ_beam2)
    sig(6)=ebeam_PN(1)
    sig(7)=ebeam_PN(2)
    sig(8)=MERGE(1d0,0d0,USE_CHARGEFORMFACTOR4PHOTON)
    sig(9)=DBLE(neutron_tagging(1))
    sig(10)=DBLE(neutron_tagging(2))
    sig(11)=MERGE(1d0,0d0,neutron_conj_sum)
    sig(12)=DBLE(IEPSILON_EDFF)
    sig(13)=MERGE(1d0,0d0,use_MC_Glauber)
    sig(14)=alphaem_elasticphoton
    sig(15)=MERGE(1d0,0d0,FFE_Legacy_ChFF_Definitions)
    sig(16)=FFE_NodeScale
    sig(17)=NEUTRON_XSIGMA
    sig(18)=DBLE(sigmaNN_inelastic_eval)
    sig(19)=ffe_xmin
    sig(20)=DBLE(proton_FF_correction)
    sig(21)=MERGE(1d0,0d0,FFE_Tag_Legacy)
    sig(22)=MERGE(1d0,0d0,USE_EXACT_ELASTIC4PHOTON) ! (0 for the defaults: their cache files are unchanged)
    sig(23)=MERGE(1d0,0d0,alpha_running4photon)
    IF(ffe_ready.AND.ALL(sig.EQ.ffe_sig))RETURN
    CALL CPU_TIME(tt0)
    t0=tt0
    ffe_sig=sig
    ffe_icoll=icoll
    ffe_Kready=.FALSE.
    ffe_Zready=.FALSE.
    ffe_WHready=.FALSE.
    ffe_ZHready=.FALSE.
    IF(ALLOCATED(ffe_K))DEALLOCATE(ffe_K)
    IF(ALLOCATED(ffe_ZL))DEALLOCATE(ffe_ZL)
    IF(alphaem_elasticphoton.LT.0d0)THEN
       IF(aqedup.GT.0d0)THEN
          ffe_alpha=aqedup
       ELSE
          ffe_alpha=0.0072992701d0
       ENDIF
    ELSE
       ffe_alpha=alphaem_elasticphoton
    ENDIF
    ! ------------------------------------------------------------------ beams --
    ffe_swap=.FALSE.
    IF(icoll.EQ.2)THEN
       ! beam 1 of the engine is always the proton
       ffe_swap=nuclearA_beam1.NE.0
       ffe_bm(1)%proton=.TRUE.
       ffe_bm(1)%chff=USE_CHARGEFORMFACTOR4PHOTON
       ffe_bm(1)%anum=1d0
       ffe_bm(1)%znum=1d0
       ffe_bm(1)%mass=mproton
       IF(ffe_swap)THEN
          ffe_bm(1)%ebeam=ebeam_PN(2)
          ffe_bm(2)%ebeam=ebeam_PN(1)
          aname=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
       ELSE
          ffe_bm(1)%ebeam=ebeam_PN(1)
          ffe_bm(2)%ebeam=ebeam_PN(2)
          aname=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
       ENDIF
       ffe_bm(1)%gamma=ffe_bm(1)%ebeam/mproton
       ffe_bm(1)%Rfm=Rproton
       ffe_bm(1)%R=Rproton/GeVm12fm
       ffe_bm(1)%afm=0d0
       ffe_bm(1)%a=0d0
       ffe_bm(1)%w=0d0
       ffe_bm(1)%bmin=ffe_bm(1)%R
       IF(ffe_bm(1)%chff)ffe_bm(1)%bmin=ffe_bm(1)%R*MERGE(0.1d0,1d-3,FFE_Legacy_ChFF_Definitions)
       inuc=2
       CALL GetNuclearInfo(aname,A,Z,R,aa,w)
       CALL ffe_set_nucleus(2,A,Z,R,aa,w)
    ELSEIF(icoll.EQ.1)THEN
       ! pp: two protons (EDFF for b >= R_p, or the dipole charge form factor)
       DO ib=1,2
          ffe_bm(ib)%proton=.TRUE.
          ffe_bm(ib)%chff=USE_CHARGEFORMFACTOR4PHOTON
          ffe_bm(ib)%anum=1d0
          ffe_bm(ib)%znum=1d0
          ffe_bm(ib)%mass=mproton
          ffe_bm(ib)%ebeam=ebeam_PN(ib)
          ffe_bm(ib)%gamma=ebeam_PN(ib)/mproton
          ffe_bm(ib)%Rfm=Rproton
          ffe_bm(ib)%R=Rproton/GeVm12fm
          ffe_bm(ib)%afm=0d0
          ffe_bm(ib)%a=0d0
          ffe_bm(ib)%w=0d0
          ffe_bm(ib)%bmin=ffe_bm(ib)%R
          IF(ffe_bm(ib)%chff)ffe_bm(ib)%bmin=ffe_bm(ib)%R*MERGE(0.1d0,1d-3,FFE_Legacy_ChFF_Definitions)
       ENDDO
    ELSE
       DO ib=1,2
          IF(ib.EQ.1)THEN
             aname=GetASymbol(nuclearA_beam1,nuclearZ_beam1)
          ELSE
             aname=GetASymbol(nuclearA_beam2,nuclearZ_beam2)
          ENDIF
          CALL GetNuclearInfo(aname,A,Z,R,aa,w)
          ffe_bm(ib)%ebeam=ebeam_PN(ib)
          CALL ffe_set_nucleus(ib,A,Z,R,aa,w)
       ENDDO
    ENDIF
    cmenergy=2d0*DSQRT(ebeam_PN(1)*ebeam_PN(2))
    ffe_sigNN=sigma_inelastic(cmenergy)*0.1d0/GeVm12fm**2 ! mb -> fm^2 -> GeV^-2
    ! pp: elastic slope of the DdE fit (arXiv:2207.03012), as in PhotonPhotonFlux_pp_eval
    ffe_b0pp=9.81d0+0.211d0*DLOG(cmenergy**2)+0.0185d0*DLOG(cmenergy**2)**2
    CALL CPU_TIME(t1)
    ffe_time(1)=t1-t0
    ffe_time(2:5)=0d0
    IF(icoll.EQ.3)THEN
       ffe_tag=neutron_tagging(1).NE.-2.OR.neutron_tagging(2).NE.-2
    ELSEIF(icoll.EQ.1)THEN
       ffe_tag=.FALSE.
    ELSE
       IF(ffe_swap)THEN
          ffe_tag=neutron_tagging(1).NE.-2
       ELSE
          ffe_tag=neutron_tagging(2).NE.-2
       ENDIF
    ENDIF
    ! the x grid of the final grids
    CALL ffe_build_xgrid(DBLE(NINT(-DLOG10(ffe_xmin))))
    ffe_MX=ffe_xg%n-1
    ! the tables are built on demand (ffe_build_tables): not at all when the grid is in the cache
    ffe_tables=.FALSE.
    ffe_lnNready=.FALSE.
    ffe_ready=.TRUE.
    IF(FFE_Grid_Cache)THEN
       IF(ffe_cache_load(1).AND.FFE_Verbose.GE.1)WRITE(*,*)"FFE: flux grid read from ",TRIM(ffe_cache_file(1))
    ENDIF
  CONTAINS
    SUBROUTINE ffe_set_nucleus(ib0,anum0,znum0,r0,dif0,w0)
      INTEGER,INTENT(IN)::ib0
      REAL(KIND(1d0)),INTENT(IN)::anum0,znum0,r0,dif0,w0
      ffe_bm(ib0)%proton=.FALSE.
      ffe_bm(ib0)%chff=USE_CHARGEFORMFACTOR4PHOTON
      ffe_bm(ib0)%anum=anum0
      ffe_bm(ib0)%znum=znum0
      ffe_bm(ib0)%Rfm=r0
      ffe_bm(ib0)%afm=dif0
      ffe_bm(ib0)%w=w0
      ffe_bm(ib0)%R=r0/GeVm12fm
      ffe_bm(ib0)%a=dif0/GeVm12fm
      ffe_bm(ib0)%mass=mN
      ffe_bm(ib0)%gamma=ffe_bm(ib0)%ebeam/mN
      IF(IEPSILON_EDFF.EQ.0)THEN
         base=ffe_bm(ib0)%R
      ELSE
         base=ffe_bm(ib0)%R+ffe_bm(ib0)%a
      ENDIF
      IF(ffe_bm(ib0)%chff)base=base*MERGE(0.1d0,1d-3,FFE_Legacy_ChFF_Definitions)
      ffe_bm(ib0)%bmin=base
      CALL ffe_normalise(ib0)
    END SUBROUTINE ffe_set_nucleus
  END SUBROUTINE ffe_init

  SUBROUTINE ffe_build_tables()
    ! the x-independent tables of the current configuration, built on demand (not needed when the
    ! grid comes from the cache): v-grids and ChFF data, survival and tagging tables, the photon
    ! densities W on the grid and the 1D tables of the PNOHAD=1 flux
    IMPLICIT NONE
    INTEGER::ib,i
    REAL(KIND(1d0))::vlo,vhi,t0,t1
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Wr
    IF(ffe_tables)RETURN
    CALL CPU_TIME(t0)
    ! v-grids: from bmin to 45/(xmin mass) (exp(-90) of the photon density beyond)
    DO ib=1,2
       vlo=DLOG(ffe_bm(ib)%bmin)
       vhi=DLOG(45d0/(ffe_xmin*ffe_bm(ib)%mass))
       CALL ffe_build_vgrid(ib,vlo,vhi)
       IF(ffe_bm(ib)%chff.AND..NOT.ffe_bm(ib)%proton)CALL ffe_chff_setup(ib)
    ENDDO
    CALL CPU_TIME(t1)
    ffe_time(1)=ffe_time(1)+t1-t0
    t0=t1
    ! kernel ingredients
    CALL ffe_build_ptable()
    IF(ffe_tag)THEN
       CALL ffe_build_tagtable()
       CALL ffe_build_tagfine()
    ENDIF
    CALL ffe_build_pfast()
    CALL CPU_TIME(t1)
    ffe_time(2)=t1-t0
    t0=t1
    ! photon densities on the grid
    IF(ALLOCATED(ffe_W))DEALLOCATE(ffe_W)
    ALLOCATE(ffe_W(ffe_MX+1,MAX(ffe_vg(1)%n,ffe_vg(2)%n),2))
    IF(ALLOCATED(ffe_lnN))DEALLOCATE(ffe_lnN)
    ALLOCATE(ffe_lnN(ffe_MX+1,2))
    ffe_W=0d0
    DO ib=1,2
       IF(ib.EQ.2.AND.ffe_same_beams())THEN
          ! identical proton beams with the rows of the inelastic branch: the rows of beam 1
          ffe_W(:,:,2)=ffe_W(:,:,1)
          ffe_lnN(:,2)=ffe_lnN(:,1)
          CYCLE
       ENDIF
       ALLOCATE(Wr(ffe_vg(ib)%n))
       DO i=1,ffe_MX+1
          CALL ffe_W_row(ib,10d0**(-ffe_xg%v(i)),Wr)
          ffe_W(i,1:ffe_vg(ib)%n,ib)=Wr
          ffe_lnN(i,ib)=DLOG(MAX(DOT_PRODUCT(ffe_vg(ib)%qw,Wr),1d-300))
       ENDDO
       DEALLOCATE(Wr)
    ENDDO
    CALL CPU_TIME(t1)
    ffe_time(3)=t1-t0
    ffe_tables=.TRUE.
    ffe_lnNready=.TRUE.
    IF(FFE_Verbose.GE.1)THEN
       WRITE(*,*)"FFE: v-grid nodes =",ffe_vg(1)%n,ffe_vg(2)%n,", P table bcut [fm] =",ffe_bcut_fm,&
            ", kernel breakpoints =",ffe_nbeta
       WRITE(*,*)"FFE: CPU setup, P, W [s] =",ffe_time(1:3)
    ENDIF
  END SUBROUTINE ffe_build_tables

  SUBROUTINE ffe_grid_prepare(k)
    ! the grid of kind k: from the cache if possible, else built here (and written to the cache).
    ! With the cache, the first process takes a lock (an atomic mkdir) and builds; the others wait
    ! for its file. A lock older than 300 s is stale (a killed process) and is removed.
    IMPLICIT NONE
    INTEGER,INTENT(IN)::k
    INTEGER::it
    LOGICAL::mine
    REAL(KIND(1d0))::t0,t1
    mine=.FALSE.
    IF(FFE_Grid_Cache)THEN
       IF(ffe_cache_load(k))RETURN
       mine=ffe_cache_lock(k)
       IF(.NOT.mine)THEN
          DO it=1,1500
             CALL EXECUTE_COMMAND_LINE('sleep 0.2')
             IF(ffe_cache_load(k))RETURN
             IF(ffe_cache_lock(k))THEN
                ! the lock was released (or was stale) without a usable file: build here
                mine=.TRUE.
                EXIT
             ENDIF
             ! nobody holds the lock and it cannot be taken (e.g. the directory does not exist or
             ! is not writable): build here without the cache
             IF(.NOT.ffe_cache_locked(k))EXIT
          ENDDO
       ENDIF
    ENDIF
    CALL ffe_build_tables()
    CALL CPU_TIME(t0)
    IF(.NOT.ffe_Kready(k))THEN
       CALL ffe_build_K(k)
       CALL CPU_TIME(t1)
       ffe_time(4)=ffe_time(4)+t1-t0
       t0=t1
    ENDIF
    CALL ffe_build_grid(k)
    CALL CPU_TIME(t1)
    ffe_time(5)=ffe_time(5)+t1-t0
    IF(FFE_Verbose.GE.1)WRITE(*,*)"FFE: CPU K, grid [s] =",ffe_time(4:5),", nodes =",ffe_nodes
    IF(FFE_Grid_Cache)THEN
       CALL ffe_cache_save(k)
       IF(mine)CALL ffe_cache_unlock(k)
    ENDIF
  END SUBROUTINE ffe_grid_prepare

  FUNCTION ffe_cache_file(k)
    ! the cache file of grid kind k for the current configuration: a 62-bit hash of the format
    ! version, the kind and the configuration signature (stored in the file and checked on reading)
    IMPLICIT NONE
    CHARACTER(len=600)::ffe_cache_file
    INTEGER,INTENT(IN)::k
    INTEGER(KIND=1),DIMENSION(8*40)::by
    INTEGER(KIND=8)::h1,h2
    INTEGER::i,n
    CHARACTER(len=16)::hx
    by=TRANSFER(ffe_sig,by)
    h1=INT(ffe_cache_version,8)*31_8+INT(k,8)
    h2=INT(ffe_cache_version,8)*131_8+INT(k,8)
    DO i=1,SIZE(by)
       h1=MOD(h1*65599_8+INT(by(i),8)+128_8,2147483647_8)
       h2=MOD(h2*69069_8+INT(by(i),8)+128_8,2147483629_8)
    ENDDO
    WRITE(hx,'(Z8.8,Z8.8)')INT(h1,4),INT(h2,4)
    n=LEN_TRIM(FFE_Grid_Cache_Dir)
    IF(n.EQ.0)THEN
       ffe_cache_file='gammaUPC_grid_'//hx//'_k'//ACHAR(48+k)//'.bin'
    ELSEIF(FFE_Grid_Cache_Dir(n:n).EQ.'/')THEN
       ffe_cache_file=FFE_Grid_Cache_Dir(1:n)//'gammaUPC_grid_'//hx//'_k'//ACHAR(48+k)//'.bin'
    ELSE
       ffe_cache_file=FFE_Grid_Cache_Dir(1:n)//'/gammaUPC_grid_'//hx//'_k'//ACHAR(48+k)//'.bin'
    ENDIF
  END FUNCTION ffe_cache_file

  FUNCTION ffe_cache_load(k)
    ! read the grid of kind k (and the 1D tables) if the file matches the configuration exactly
    IMPLICIT NONE
    LOGICAL::ffe_cache_load
    INTEGER,INTENT(IN)::k
    CHARACTER(len=8)::tag
    INTEGER::u,ios,ver,kk,mx,n
    LOGICAL::ex
    REAL(KIND(1d0)),DIMENSION(40)::sig
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Z,L
    REAL(KIND(1d0))::chk
    ffe_cache_load=.FALSE.
    INQUIRE(FILE=TRIM(ffe_cache_file(k)),EXIST=ex)
    IF(.NOT.ex)RETURN
    OPEN(NEWUNIT=u,FILE=TRIM(ffe_cache_file(k)),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',&
         ACTION='READ',IOSTAT=ios)
    IF(ios.NE.0)RETURN
    READ(u,IOSTAT=ios)tag,ver,kk,mx,n,sig
    IF(ios.NE.0.OR.tag.NE.'GUPCGRID'.OR.ver.NE.ffe_cache_version.OR.kk.NE.k.OR.mx.NE.ffe_MX&
         .OR.n.NE.ffe_MX+1.OR.ANY(sig.NE.ffe_sig))THEN
       CLOSE(u)
       RETURN
    ENDIF
    ALLOCATE(Z(mx+1,mx+1),L(mx+1,2))
    READ(u,IOSTAT=ios)Z,L,chk,tag
    CLOSE(u)
    IF(ios.NE.0.OR.tag.NE.'GUPCEND!'.OR.chk.NE.SUM(Z)+3d0*SUM(L))RETURN
    IF(.NOT.ALLOCATED(ffe_ZL))ALLOCATE(ffe_ZL(ffe_MX+1,ffe_MX+1,3))
    ffe_ZL(:,:,k)=Z
    IF(.NOT.ffe_lnNready)THEN
       IF(ALLOCATED(ffe_lnN))DEALLOCATE(ffe_lnN)
       ALLOCATE(ffe_lnN(ffe_MX+1,2))
       ffe_lnN=L
       ffe_lnNready=.TRUE.
    ENDIF
    ffe_Zready(k)=.TRUE.
    ffe_cache_load=.TRUE.
  END FUNCTION ffe_cache_load

  SUBROUTINE ffe_cache_save(k)
    ! write the grid of kind k to a temporary file and rename it into place, so that no process
    ! reads a partial file (failures only mean that there is no cache)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::k
    CHARACTER(len=600)::f
    CHARACTER(len=700)::tmp
    CHARACTER(len=24)::sfx
    INTEGER::u,ios,es,cs
    INTEGER(KIND=8)::cnt
    REAL::r
    IF(.NOT.ffe_Zready(k).OR..NOT.ffe_lnNready)RETURN
    f=ffe_cache_file(k)
    CALL SYSTEM_CLOCK(cnt)
    CALL RANDOM_NUMBER(r)
    WRITE(sfx,'(Z16.16,I0)')cnt,INT(r*1d6)
    tmp=TRIM(f)//'.tmp'//TRIM(sfx)
    OPEN(NEWUNIT=u,FILE=TRIM(tmp),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE',&
         IOSTAT=ios)
    IF(ios.NE.0)RETURN
    WRITE(u,IOSTAT=ios)'GUPCGRID',ffe_cache_version,k,ffe_MX,ffe_MX+1,ffe_sig,ffe_ZL(:,:,k),ffe_lnN,&
         SUM(ffe_ZL(:,:,k))+3d0*SUM(ffe_lnN),'GUPCEND!'
    CLOSE(u)
    IF(ios.NE.0)THEN
       CALL EXECUTE_COMMAND_LINE('rm -f "'//TRIM(tmp)//'"',WAIT=.TRUE.,EXITSTAT=es,CMDSTAT=cs)
       RETURN
    ENDIF
    CALL EXECUTE_COMMAND_LINE('mv -f "'//TRIM(tmp)//'" "'//TRIM(f)//'"',WAIT=.TRUE.,EXITSTAT=es,CMDSTAT=cs)
    IF(cs.NE.0.OR.es.NE.0)CALL EXECUTE_COMMAND_LINE('rm -f "'//TRIM(tmp)//'"',WAIT=.TRUE.,&
         EXITSTAT=es,CMDSTAT=cs)
  END SUBROUTINE ffe_cache_save

  FUNCTION ffe_cache_lock(k)
    ! take the build lock of grid kind k (mkdir is atomic); a lock older than 300 s is removed first
    IMPLICIT NONE
    LOGICAL::ffe_cache_lock
    INTEGER,INTENT(IN)::k
    CHARACTER(len=610)::d
    INTEGER::es,cs,u,ios
    REAL(KIND(1d0))::t
    LOGICAL::ex
    d=TRIM(ffe_cache_file(k))//'.lock'
    INQUIRE(FILE=TRIM(d)//'/t',EXIST=ex)
    IF(ex)THEN
       OPEN(NEWUNIT=u,FILE=TRIM(d)//'/t',STATUS='OLD',ACTION='READ',IOSTAT=ios)
       IF(ios.EQ.0)THEN
          READ(u,*,IOSTAT=ios)t
          CLOSE(u)
          IF(ios.EQ.0.AND.ffe_now()-t.GT.300d0)CALL EXECUTE_COMMAND_LINE('rm -rf "'//TRIM(d)//'"',&
               WAIT=.TRUE.,EXITSTAT=es,CMDSTAT=cs)
       ENDIF
    ENDIF
    CALL EXECUTE_COMMAND_LINE('mkdir "'//TRIM(d)//'" 2>/dev/null',WAIT=.TRUE.,EXITSTAT=es,CMDSTAT=cs)
    ffe_cache_lock=cs.EQ.0.AND.es.EQ.0
    IF(ffe_cache_lock)THEN
       OPEN(NEWUNIT=u,FILE=TRIM(d)//'/t',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios)
       IF(ios.EQ.0)THEN
          WRITE(u,*)ffe_now()
          CLOSE(u)
       ENDIF
    ENDIF
  END FUNCTION ffe_cache_lock

  FUNCTION ffe_cache_locked(k)
    ! another process holds the build lock of grid kind k (its time stamp is written right after
    ! the mkdir)
    IMPLICIT NONE
    LOGICAL::ffe_cache_locked
    INTEGER,INTENT(IN)::k
    INQUIRE(FILE=TRIM(ffe_cache_file(k))//'.lock/t',EXIST=ffe_cache_locked)
  END FUNCTION ffe_cache_locked

  SUBROUTINE ffe_cache_unlock(k)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::k
    INTEGER::es,cs
    CALL EXECUTE_COMMAND_LINE('rm -rf "'//TRIM(ffe_cache_file(k))//'.lock"',WAIT=.TRUE.,EXITSTAT=es,&
         CMDSTAT=cs)
  END SUBROUTINE ffe_cache_unlock

  FUNCTION ffe_now()
    ! wall-clock time in seconds (days since 1970 from the civil date; for lock ages only)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_now
    INTEGER,DIMENSION(8)::v
    INTEGER::y,m,d,era,yoe,doy,doe
    CALL DATE_AND_TIME(VALUES=v)
    y=v(1)
    m=v(2)
    d=v(3)
    IF(m.LE.2)y=y-1
    era=y/400
    yoe=y-era*400
    doy=(153*(m+MERGE(-3,9,m.GT.2))+2)/5+d-1
    doe=yoe*365+yoe/4-yoe/100+doy
    ffe_now=86400d0*DBLE(era*146097+doe-719468)+3600d0*v(5)+60d0*v(6)+v(7)+1d-3*v(8)-60d0*v(4)
  END FUNCTION ffe_now

  SUBROUTINE ffe_build_tagtable()
    ! P_tag(b) on a uniform ln b grid (b in GeV^-1). The lowest-order emission probabilities
    ! m_k(b) = Z_p^2 alpha int dln(1/x) sigma_k(E_target) N_p(b,x) of each projectile -> target
    ! direction are integrated by Gauss-Legendre in ln E_target with breakpoints at the jumps and
    ! kinks of the photo-absorption cross sections (sigma evaluated once per node and reused for
    ! every b), then combined with the channel combinatorics of Pbreak_AB_WoodsSaxon /
    ! Pbreak_pA_WoodsSaxon (Xn: 1-exp(-m_Xn), 0n: exp(-m_Xn), 1n: m_1 exp(-m_Xn), ...).
    ! FFE_Tag_Legacy = .TRUE. fills the table from the legacy Pbreak functions instead.
    IMPLICIT NONE
    include 'run90.inc'
    INTEGER::i,k,in1,in2,tagn,kn1,kn2
    REAL(KIND(1d0))::vmax,b,f1,f2,fc1,fc2
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::m
    LOGICAL::conj
    ffe_Th=0.02d0
    ffe_Tv0=DLOG(0.5d0/GeVm12fm)
    vmax=DLOG(2d0*(DEXP(ffe_vg(1)%v(ffe_vg(1)%n))+DEXP(ffe_vg(2)%v(ffe_vg(2)%n))))
    ffe_nT=INT((vmax-ffe_Tv0)/ffe_Th)+2
    IF(ALLOCATED(ffe_PtagT))DEALLOCATE(ffe_PtagT)
    ALLOCATE(ffe_PtagT(ffe_nT))
    IF(FFE_Tag_Legacy)THEN
       DO i=1,ffe_nT
          b=DEXP(ffe_Tv0+DBLE(i-1)*ffe_Th)
          IF(ffe_icoll.EQ.3)THEN
             ffe_PtagT(i)=Pbreak_AB_WoodsSaxon(neutron_tagging(1),neutron_tagging(2),b,neutron_conj_sum)
          ELSEIF(ffe_swap)THEN
             ffe_PtagT(i)=Pbreak_pA_WoodsSaxon(neutron_tagging(1),b)
          ELSE
             ffe_PtagT(i)=Pbreak_pA_WoodsSaxon(neutron_tagging(2),b)
          ENDIF
       ENDDO
       RETURN
    ENDIF
    ALLOCATE(m(ffe_nT,0:4,2))
    m=0d0
    IF(ffe_icoll.EQ.3)THEN
       in1=neutron_tagging(1)
       in2=neutron_tagging(2)
       conj=neutron_conj_sum.AND.in1.NE.in2
       ! direction 1: beam 1 -> beam 2 (tag of beam 2); direction 2: beam 2 -> beam 1
       kn2=-2
       IF(in2.NE.-2)kn2=MAX(kn2,in2)
       IF(conj.AND.in1.NE.-2)kn2=MAX(kn2,in1)
       kn1=-2
       IF(in1.NE.-2)kn1=MAX(kn1,in1)
       IF(conj.AND.in2.NE.-2)kn1=MAX(kn1,in2)
       IF(kn2.NE.-2)CALL ffe_mtab(1,2,MAX(kn2,0),m(:,:,1))
       IF(kn1.NE.-2)CALL ffe_mtab(2,1,MAX(kn1,0),m(:,:,2))
       DO i=1,ffe_nT
          f2=ffe_ptag_f(in2,m(i,:,1))
          f1=ffe_ptag_f(in1,m(i,:,2))
          ffe_PtagT(i)=f1*f2
          IF(conj)THEN
             fc2=ffe_ptag_f(in1,m(i,:,1))
             fc1=ffe_ptag_f(in2,m(i,:,2))
             ffe_PtagT(i)=ffe_PtagT(i)+fc1*fc2
          ENDIF
       ENDDO
    ELSE
       ! the proton (engine beam 1) breaks up the nucleus (engine beam 2)
       IF(ffe_swap)THEN
          tagn=neutron_tagging(1)
       ELSE
          tagn=neutron_tagging(2)
       ENDIF
       CALL ffe_mtab(1,2,MAX(tagn,0),m(:,:,1))
       DO i=1,ffe_nT
          ffe_PtagT(i)=ffe_ptag_f(tagn,m(i,:,1))
       ENDDO
    ENDIF
    DEALLOCATE(m)
  END SUBROUTINE ffe_build_tagtable

  FUNCTION ffe_ptag_f(itag,mm)
    ! probability of tag itag of the target from the lowest-order m_k (k = 0: Xn, 1..4: 1n..4n)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_ptag_f
    INTEGER,INTENT(IN)::itag
    REAL(KIND(1d0)),DIMENSION(0:4),INTENT(IN)::mm
    REAL(KIND(1d0))::p0
    p0=DEXP(-mm(0))
    SELECT CASE(itag)
    CASE(-2)
       ffe_ptag_f=1d0
    CASE(-1)
       ! 1-exp(-m) without cancellation at small m
       IF(mm(0).LT.1d-3)THEN
          ffe_ptag_f=mm(0)*(1d0-mm(0)/2d0*(1d0-mm(0)/3d0*(1d0-mm(0)/4d0)))
       ELSE
          ffe_ptag_f=1d0-p0
       ENDIF
    CASE(0)
       ffe_ptag_f=p0
    CASE(1)
       ffe_ptag_f=mm(1)*p0
    CASE(2)
       ffe_ptag_f=(mm(1)**2/2d0+mm(2))*p0
    CASE(3)
       ffe_ptag_f=(mm(1)**3/6d0+mm(1)*mm(2)+mm(3))*p0
    CASE DEFAULT
       ffe_ptag_f=(mm(1)**4/24d0+mm(1)**2*mm(2)/2d0+mm(1)*mm(3)+mm(2)**2/2d0+mm(4))*p0
    END SELECT
  END FUNCTION ffe_ptag_f

  SUBROUTINE ffe_mtab(ip,it,kmax,mt)
    ! m_k(b_i), k = 0 (Xn), 1..kmax (1n..), for projectile ip -> target it on the ln b table
    USE photoabsorption
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ip,it,kmax
    REAL(KIND(1d0)),DIMENSION(ffe_nT,0:4),INTENT(OUT)::mt
    REAL(KIND(1d0)),DIMENSION(0:4)::emin,emax
    INTEGER,PARAMETER::nscan=200000
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::br,tn,wn,Wr,Nbt
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::sg
    REAL(KIND(1d0)),DIMENSION(8)::xg,wg
    REAL(KIND(1d0))::fac,gt,tlo,thi,dt,t,x,b,q,z,pref,hp,cc,smax,u,ig2
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::sc
    REAL(KIND(1d0))::smx,d1,tj
    INTEGER::k,i,j,nb,np,nn,g,ns,nsub,i0
    REAL(KIND(1d0)),DIMENSION(8)::wt
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    REAL(KIND(1d0)),EXTERNAL::GUPC_BESSK0E,GUPC_BESSK1E
    XSIGMA=NEUTRON_XSIGMA
    emin=(/Emin_Xn,Emin_1n,Emin_2n,Emin_3n,Emin_4n/)
    emax=(/Emax_Xn,Emax_1n,Emax_2n,Emax_3n,Emax_4n/)
    gt=ffe_bm(it)%gamma
    fac=gt+DSQRT(gt*gt-1d0)
    ! E_target = x * mass_p * gamma_p * fac (GeV); integrate in t = ln E_target
    pref=ffe_bm(ip)%mass*ffe_bm(ip)%gamma*fac
    tlo=DLOG(MINVAL(emin(0:kmax)))
    thi=DLOG(pref)
    ALLOCATE(br(4*nscan+40),sc(nscan))
    nb=0
    DO k=0,kmax
       nb=nb+1
       br(nb)=DLOG(emin(k))
       IF(emax(k).GT.0d0)THEN
          nb=nb+1
          br(nb)=DLOG(MIN(emax(k),pref))
       ENDIF
    ENDDO
    ! jumps and kinks of each cross section from a fine scan in t; a jump candidate is confirmed
    ! by bisection (the step must persist down to 1e-12 in t)
    dt=(thi-tlo)/DBLE(nscan-1)
    DO k=0,kmax
       DO i=1,nscan
          sc(i)=ffe_sigma(k,DEXP(tlo+DBLE(i-1)*dt),it)
       ENDDO
       smx=MAXVAL(DABS(sc))
       DO i=3,nscan-2
          d1=DABS(sc(i+1)-sc(i))
          IF(d1.GT.4d0*MAX(DABS(sc(i)-sc(i-1)),DABS(sc(i+2)-sc(i+1))).AND.d1.GT.1d-9*smx)THEN
             tj=ffe_locate_jump(k,it,tlo+DBLE(i-1)*dt,tlo+DBLE(i)*dt)
             IF(DABS(ffe_sigma(k,DEXP(tj+1d-12),it)-ffe_sigma(k,DEXP(tj-1d-12),it)).GT.0.5d0*d1)THEN
                nb=nb+1
                br(nb)=tj
                CYCLE
             ENDIF
          ENDIF
          IF(DABS(sc(i+1)-2d0*sc(i)+sc(i-1)).GT.&
               4d0*MAX(DABS(sc(i)-2d0*sc(i-1)+sc(i-2)),DABS(sc(i+2)-2d0*sc(i+1)+sc(i)))&
               .AND.DABS(sc(i+1)-2d0*sc(i)+sc(i-1)).GT.1d-9*smx)THEN
             nb=nb+1
             br(nb)=tlo+DBLE(i-1)*dt
          ENDIF
       ENDDO
    ENDDO
    CALL ffe_breaks(nb,br,tlo,thi,1d-9,np)
    ! nodes: panels <= 0.05 below E = 0.2 GeV (giant resonances), <= 0.25 above
    CALL ffe_gauleg(8,xg,wg)
    nn=0
    ALLOCATE(tn(8*(np+INT((thi-tlo)/0.05d0)+10)),wn(8*(np+INT((thi-tlo)/0.05d0)+10)))
    DO j=1,np-1
       hp=MERGE(0.05d0,0.25d0,br(j).LT.DLOG(0.2d0))
       nsub=MAX(1,CEILING((br(j+1)-br(j))/hp-1d-9))
       DO ns=1,nsub
          cc=br(j)+(DBLE(ns)-0.5d0)*(br(j+1)-br(j))/DBLE(nsub)
          DO g=1,8
             nn=nn+1
             tn(nn)=cc+0.5d0*(br(j+1)-br(j))/DBLE(nsub)*xg(g)
             wn(nn)=0.5d0*(br(j+1)-br(j))/DBLE(nsub)*wg(g)
          ENDDO
       ENDDO
    ENDDO
    ! cross sections at the nodes (GeV^-2), zero outside each channel's energy window
    ALLOCATE(sg(nn,0:4))
    sg=0d0
    DO i=1,nn
       DO k=0,kmax
          IF(DEXP(tn(i)).LT.emin(k)*(1d0-1d-12))CYCLE
          IF(emax(k).GT.0d0.AND.DEXP(tn(i)).GT.emax(k)*(1d0+1d-12))CYCLE
          sg(i,k)=ffe_sigma(k,DEXP(tn(i)),it)*0.1d0/GeVm12fm**2
       ENDDO
    ENDDO
    ! photon number density of the projectile at the table b's, node by node
    mt=0d0
    ALLOCATE(Nbt(ffe_nT),Wr(ffe_vg(ip)%n))
    ig2=1d0/ffe_bm(ip)%gamma**2
    smax=0d0
    IF(ffe_bm(ip)%chff.AND..NOT.ffe_bm(ip)%proton)smax=ffe_cf(ip)%smax
    DO i=1,nn
       x=DEXP(tn(i))/pref
       q=x*ffe_bm(ip)%mass
       IF(ffe_bm(ip)%chff)THEN
          CALL ffe_W_row(ip,x,Wr)
       ENDIF
       DO j=1,ffe_nT
          b=DEXP(ffe_Tv0+DBLE(j-1)*ffe_Th)
          z=q*b
          IF(z.GT.340d0)THEN
             Nbt(j)=0d0
          ELSEIF(ffe_bm(ip)%chff.AND.(ffe_bm(ip)%proton.OR.b.LT.smax).AND.&
               b.GE.ffe_vg(ip)%b(1).AND.b.LE.ffe_vg(ip)%b(ffe_vg(ip)%n))THEN
             CALL ffe_lagw(ffe_vg(ip),DLOG(b),8,i0,wt)
             Nbt(j)=SUM(wt*Wr(i0+1:i0+8))/(b*b)
          ELSE
             ! point-like (EDFF, or ChFF outside the charge)
             Nbt(j)=q*q*DEXP(-2d0*z)*(GUPC_BESSK1E(z)**2+ig2*GUPC_BESSK0E(z)**2)/pi**2
          ENDIF
       ENDDO
       DO k=0,kmax
          IF(sg(i,k).EQ.0d0)CYCLE
          mt(:,k)=mt(:,k)+wn(i)*sg(i,k)*Nbt
       ENDDO
    ENDDO
    mt=mt*ffe_bm(ip)%znum**2*ffe_alpha
    DEALLOCATE(br,tn,wn,sg,Nbt,Wr,sc)
  END SUBROUTINE ffe_mtab

  FUNCTION ffe_sigma(k,E,it)
    ! photo-absorption cross sections of target it (mb) at E (GeV, target frame)
    USE photoabsorption
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_sigma
    INTEGER,INTENT(IN)::k,it
    REAL(KIND(1d0)),INTENT(IN)::E
    INTEGER::A,Z
    A=NINT(ffe_bm(it)%anum)
    Z=NINT(ffe_bm(it)%znum)
    SELECT CASE(k)
    CASE(0)
       ffe_sigma=GAMMAABS(E*1d3,A,Z)
    CASE(1)
       ffe_sigma=GAMMATO1N(E*1d3,A,Z)
    CASE(2)
       ffe_sigma=GAMMATO2N(E*1d3,A,Z)
    CASE(3)
       ffe_sigma=GAMMATO3N(E*1d3,A,Z)
    CASE DEFAULT
       ffe_sigma=GAMMATO4N(E*1d3,A,Z)
    END SELECT
  END FUNCTION ffe_sigma

  FUNCTION ffe_locate_jump(k,it,t1,t2)
    ! bisection for the position of a jump of sigma_k between ln E = t1 and t2
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_locate_jump
    INTEGER,INTENT(IN)::k,it
    REAL(KIND(1d0)),INTENT(IN)::t1,t2
    REAL(KIND(1d0))::a,b,c,sa,sb,sc
    INTEGER::j
    a=t1
    b=t2
    sa=ffe_sigma(k,DEXP(a),it)
    sb=ffe_sigma(k,DEXP(b),it)
    DO j=1,50
       c=0.5d0*(a+b)
       sc=ffe_sigma(k,DEXP(c),it)
       IF(DABS(sc-sa).LT.DABS(sc-sb))THEN
          a=c
          sa=sc
       ELSE
          b=c
          sb=sc
       ENDIF
    ENDDO
    ffe_locate_jump=0.5d0*(a+b)
  END FUNCTION ffe_locate_jump

  FUNCTION ffe_flux(icoll,x1,x2,knd,mode)
    ! knd 1: flux; 2: <Delta b> numerator (fm); 3: <Delta b> numerator with PNOHAD=1 (fm)
    ! mode 0: grid; 1: direct; 2: PNOHAD=1 flux from the 1D tables; 3: PNOHAD=1 flux, direct
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_flux
    INTEGER,INTENT(IN)::icoll,knd,mode
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    REAL(KIND(1d0))::xa,xb,t0,t1
    LOGICAL::zero
    ffe_flux=0d0
    IF(x1.LE.0d0.OR.x2.LE.0d0.OR.x1.GE.1d0.OR.x2.GE.1d0)RETURN
    CALL ffe_init(icoll,MIN(x1,x2))
    IF(icoll.EQ.2.AND.ffe_swap)THEN
       xa=x2
       xb=x1
    ELSE
       xa=x1
       xb=x2
    ENDIF
    SELECT CASE(mode)
    CASE(0)
       IF(.NOT.ffe_Zready(knd))CALL ffe_grid_prepare(knd)
       IF(ffe_ourW())THEN
          ! the proton photon of the inelastic branch vanishes at x = 1: the patch in s = -log10(1-x) near x = 1
          ! (ffe_patch_val), else the grid; direct evaluation where neither applies (a zero node in the stencil)
          IF(ffe_patch_side(xa,xb).GT.0)THEN
             ffe_flux=ffe_patch_val(knd,xa,xb,zero)
          ELSE
             ffe_flux=ffe_grid_val(knd,xa,xb,zero)
          ENDIF
          IF(zero)THEN
             CALL ffe_build_tables()
             IF(.NOT.ffe_Kready(knd))CALL ffe_build_K(knd)
             ffe_flux=ffe_direct(knd,xa,xb)
          ENDIF
       ELSE
          ffe_flux=ffe_grid_val(knd,xa,xb)
       ENDIF
    CASE(1)
       CALL ffe_build_tables()
       IF(.NOT.ffe_Kready(knd))THEN
          CALL CPU_TIME(t0)
          CALL ffe_build_K(knd)
          CALL CPU_TIME(t1)
          ffe_time(4)=ffe_time(4)+t1-t0
       ENDIF
       ffe_flux=ffe_direct(knd,xa,xb)
    CASE(2)
       IF(.NOT.ffe_lnNready)CALL ffe_build_tables()
       ffe_flux=ffe_pnohad1(xa,xb,.FALSE.)
    CASE DEFAULT
       CALL ffe_build_tables()
       ffe_flux=ffe_pnohad1(xa,xb,.TRUE.)
    END SELECT
  END FUNCTION ffe_flux

  FUNCTION FFE_Flux_AB(x1,x2,mode)
    ! public entry for tests: AB flux (mode as in ffe_flux)
    IMPLICIT NONE
    REAL(KIND(1d0))::FFE_Flux_AB
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    INTEGER,INTENT(IN)::mode
    FFE_Flux_AB=ffe_flux(3,x1,x2,1,mode)
  END FUNCTION FFE_Flux_AB

  FUNCTION FFE_DeltaB_AB(x1,x2,mode,pnohad1)
    ! public entry for tests: <Delta b> numerator (fm), mode 0 grid / 1 direct
    IMPLICIT NONE
    REAL(KIND(1d0))::FFE_DeltaB_AB
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    INTEGER,INTENT(IN)::mode
    LOGICAL,INTENT(IN)::pnohad1
    FFE_DeltaB_AB=ffe_flux(3,x1,x2,MERGE(3,2,pnohad1),mode)
  END FUNCTION FFE_DeltaB_AB

  SUBROUTINE FFE_Reset()
    IMPLICIT NONE
    ffe_ready=.FALSE.
    ffe_Kready=.FALSE.
    ffe_Zready=.FALSE.
    ffe_xmin=1d-8
  END SUBROUTINE FFE_Reset

  SUBROUTINE FFE_Info(times,nodes,nv1,nv2,nguard)
    ! CPU times (setup, P, W, K, grid) and sizes of the last build; nguard: number of grid
    ! values replaced by the PNOHAD=1 flux by the last-resort guard since the start of the run
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(5),INTENT(OUT)::times
    INTEGER,INTENT(OUT)::nodes,nv1,nv2
    INTEGER,INTENT(OUT),OPTIONAL::nguard
    IF(PRESENT(nguard))nguard=ffe_nguard
    times=ffe_time(1:5)
    nodes=ffe_nodes
    nv1=ffe_vg(1)%n
    nv2=ffe_vg(2)%n
  END SUBROUTINE FFE_Info


  FUNCTION FFE_Test_Value(what,icoll,x,b)
    ! test hooks (benchmark/bench_blocks.F90): after initialising the configuration of icoll,
    !  what = 1: That_AB(b) (or That_A for pA) in fm^-2 from the chi table, b in fm
    !         2: P_had(b), b in fm
    !         3: N(x,b) of engine beam 1 (unit charge, alpha = 1), b in GeV^-1
    !         4: N(x,b) of engine beam 2
    !         5: P_tag(b), b in GeV^-1
    !         6: the node of the ln b grid of engine beam 1 nearest to b (GeV^-1)
    IMPLICIT NONE
    REAL(KIND(1d0))::FFE_Test_Value
    INTEGER,INTENT(IN)::what,icoll
    REAL(KIND(1d0)),INTENT(IN)::x,b
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Wr
    REAL(KIND(1d0)),DIMENSION(8)::wt
    REAL(KIND(1d0))::y,bg
    INTEGER::i0,ib
    CALL ffe_init(icoll,MIN(1d-8,MAX(x,1d-300)))
    CALL ffe_build_tables()
    SELECT CASE(what)
    CASE(1,2)
       bg=b/GeVm12fm
       IF(bg.GE.ffe_bcut)THEN
          y=-1d300
       ELSE
          CALL ffe_lagu(bg/ffe_hb,ffe_nP,8,i0,wt)
          y=SUM(wt*ffe_lnchi(i0+1:i0+8))
       ENDIF
       IF(what.EQ.2)THEN
          FFE_Test_Value=ffe_Phad(bg)
       ELSEIF(y.LT.-690d0)THEN
          FFE_Test_Value=0d0
       ELSEIF(ffe_icoll.EQ.3)THEN
          FFE_Test_Value=DEXP(y)/(ffe_sigNN*ffe_bm(1)%anum*ffe_bm(2)%anum)/GeVm12fm**2
       ELSE
          FFE_Test_Value=DEXP(y)/(ffe_sigNN*ffe_bm(2)%anum)/GeVm12fm**2
       ENDIF
    CASE(3,4)
       ib=what-2
       ALLOCATE(Wr(ffe_vg(ib)%n))
       CALL ffe_W_row(ib,x,Wr)
       CALL ffe_lagw(ffe_vg(ib),DLOG(b),8,i0,wt)
       FFE_Test_Value=SUM(wt*Wr(i0+1:i0+8))/(b*b)
       DEALLOCATE(Wr)
    CASE(5)
       FFE_Test_Value=ffe_Ptag(b)
    CASE(6)
       i0=MINLOC(DABS(DLOG(ffe_vg(1)%b(1:ffe_vg(1)%n)/b)),1)
       FFE_Test_Value=ffe_vg(1)%b(i0)
    CASE DEFAULT
       FFE_Test_Value=0d0
    END SELECT
  END FUNCTION FFE_Test_Value


  FUNCTION ffe_api(icoll,x1,x2,knd,force,grid)
    ! dispatch of the public flux and <Delta b> functions to the engine
    ! knd 1: flux (force: PNOHAD=1 flux); knd 2: <Delta b> numerator in fm (force: PNOHAD=1)
    ! grid: interpolate the x grid (else evaluate directly)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_api
    INTEGER,INTENT(IN)::icoll,knd
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN)::force,grid
    REAL(KIND(1d0))::v,p1
    ffe_api=0d0
    IF(x1.LE.0d0.OR.x2.LE.0d0.OR.x1.GE.1d0.OR.x2.GE.1d0)RETURN
    IF(knd.EQ.1)THEN
       IF(force)THEN
          ffe_api=ffe_flux(icoll,x1,x2,1,MERGE(2,3,grid))
       ELSEIF(grid)THEN
          v=ffe_flux(icoll,x1,x2,1,0)
          p1=ffe_flux(icoll,x1,x2,1,2)
          IF(v.NE.v.OR.v.LT.0d0.OR.v.GT.1d2*p1)THEN
             ! last-resort guard, as in the legacy code: warned (first 5 times) and counted
             ! (FFE_Info, argument nguard)
             ffe_nguard=ffe_nguard+1
             IF(ffe_nguard.LE.5)WRITE(*,*)"WARNING: photon-photon flux grid unusable at (x1,x2) =",&
                  x1,x2,"; the PNOHAD=1 flux is used instead"
             IF(ffe_nguard.EQ.5)WRITE(*,*)"WARNING: further warnings of this kind are suppressed"
             v=p1
          ENDIF
          ffe_api=v
       ELSE
          ffe_api=ffe_flux(icoll,x1,x2,1,1)
       ENDIF
    ELSE
       ffe_api=ffe_flux(icoll,x1,x2,MERGE(3,2,force),MERGE(0,1,grid))
    ENDIF
  END FUNCTION ffe_api

  LOGICAL FUNCTION ffe_patch_beam(ib)
    ! the patch near x = 1 is used for beam ib: a proton with the W rows of the inelastic branch
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ib
    ffe_patch_beam=ffe_ourW().AND.ffe_bm(ib)%proton
  END FUNCTION ffe_patch_beam

  INTEGER FUNCTION ffe_patch_side(xa,xb)
    ! 1: xa (engine beam 1) in the patch range, xb on the x grid (xb < FFE_XHLO); 2: the reverse; 0: neither
    ! (both below FFE_XHLO: the grid); -1: not covered (both near 1, or x > FFE_XHHI): direct evaluation
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::xa,xb
    LOGICAL::na,nb
    na=xa.GE.FFE_XHLO.AND.ffe_patch_beam(1)
    nb=xb.GE.FFE_XHLO.AND.ffe_patch_beam(2)
    ffe_patch_side=0
    IF(na.AND.nb)THEN
       ffe_patch_side=-1
    ELSEIF(na)THEN
       ffe_patch_side=1
       IF(xa.GT.FFE_XHHI.OR.xb.GE.FFE_XHLO)ffe_patch_side=-1
    ELSEIF(nb)THEN
       ffe_patch_side=2
       IF(xb.GT.FFE_XHHI.OR.xa.GE.FFE_XHLO)ffe_patch_side=-1
    ENDIF
  END FUNCTION ffe_patch_side

  SUBROUTINE ffe_build_patchW()
    ! W rows and ln(int W dv) at the patch x of the beams with the patch
    IMPLICIT NONE
    INTEGER::ib,k,nv
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Wr
    IF(ffe_WHready)RETURN
    CALL ffe_build_tables()
    nv=MAX(ffe_vg(1)%n,ffe_vg(2)%n)
    IF(ALLOCATED(ffe_WH))DEALLOCATE(ffe_WH)
    IF(ALLOCATED(ffe_lnNH))DEALLOCATE(ffe_lnNH)
    ALLOCATE(ffe_WH(FFE_NH,nv,2),ffe_lnNH(FFE_NH,2))
    ffe_WH=0d0
    ffe_lnNH=-690d0
    DO ib=1,2
       IF(.NOT.ffe_patch_beam(ib))CYCLE
       IF(ib.EQ.2.AND.ffe_same_beams())THEN
          ffe_WH(:,:,2)=ffe_WH(:,:,1)
          ffe_lnNH(:,2)=ffe_lnNH(:,1)
          CYCLE
       ENDIF
       ALLOCATE(Wr(ffe_vg(ib)%n))
       DO k=1,FFE_NH
          CALL ffe_W_row(ib,ffe_xh(k),Wr)
          ffe_WH(k,1:ffe_vg(ib)%n,ib)=Wr
          ffe_lnNH(k,ib)=DLOG(MAX(DOT_PRODUCT(ffe_vg(ib)%qw,Wr),1d-300))
       ENDDO
       DEALLOCATE(Wr)
    ENDDO
    ffe_WHready=.TRUE.
  END SUBROUTINE ffe_build_patchW

  FUNCTION ffe_xh(k) ! x of node k of the patch
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_xh
    INTEGER,INTENT(IN)::k
    ffe_xh=1d0-10d0**(-(FFE_SH0+(k-1)*FFE_HH))
  END FUNCTION ffe_xh

  SUBROUTINE ffe_build_patch(k)
    ! ln(x1 x2 Z) of kind k with x1 (side 1) or x2 (side 2) at the patch nodes and the other x on the x grid
    IMPLICIT NONE
    INTEGER,INTENT(IN)::k
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::T,S
    REAL(KIND(1d0))::pref,xa,xb
    INTEGER::i,j,n1,n2
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    CALL ffe_build_patchW()
    IF(.NOT.ffe_Kready(k))CALL ffe_build_K(k)
    n1=ffe_vg(1)%n
    n2=ffe_vg(2)%n
    pref=twopi*ffe_alpha**2*ffe_bm(1)%znum**2*ffe_bm(2)%znum**2
    IF(.NOT.ALLOCATED(ffe_ZH))ALLOCATE(ffe_ZH(FFE_NH,ffe_MX+1,3,2))
    ffe_ZH(:,:,k,:)=-745d0
    IF(ffe_patch_beam(1))THEN
       ALLOCATE(T(FFE_NH,n2),S(FFE_NH,ffe_MX+1))
       T=MATMUL(ffe_WH(:,1:n1,1),ffe_K(:,:,k))
       S=MATMUL(T,TRANSPOSE(ffe_W(:,1:n2,2)))
       DO j=1,ffe_MX+1
          xb=10d0**(-ffe_xg%v(j))
          DO i=1,FFE_NH
             xa=ffe_xh(i)
             IF(alpha_running4photon.AND.ffe_icoll.NE.3)S(i,j)=S(i,j)*ffe_anorm2(xa,xb)
             IF(S(i,j).GT.0d0)ffe_ZH(i,j,k,1)=DLOG(pref*S(i,j))
          ENDDO
       ENDDO
       DEALLOCATE(T,S)
    ENDIF
    IF(ffe_patch_beam(2))THEN
       ALLOCATE(T(ffe_MX+1,n2),S(ffe_MX+1,FFE_NH))
       T=MATMUL(ffe_W(:,1:n1,1),ffe_K(:,:,k))
       S=MATMUL(T,TRANSPOSE(ffe_WH(:,1:n2,2)))
       DO j=1,FFE_NH
          xb=ffe_xh(j)
          DO i=1,ffe_MX+1
             xa=10d0**(-ffe_xg%v(i))
             IF(alpha_running4photon.AND.ffe_icoll.NE.3)S(i,j)=S(i,j)*ffe_anorm2(xa,xb)
             IF(S(i,j).GT.0d0)ffe_ZH(j,i,k,2)=DLOG(pref*S(i,j))
          ENDDO
       ENDDO
       DEALLOCATE(T,S)
    ENDIF
    ffe_ZHready(k)=.TRUE.
  END SUBROUTINE ffe_build_patch

  FUNCTION ffe_patch_val(k,xa,xb,zero)
    ! lookup in the patch: 8-point Lagrange in s = -log10(1-x) for the x near 1, 8-point on the x grid for the
    ! other; zero: a node of the stencil has Z = 0, or the point is not covered (then: direct evaluation)
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_patch_val
    INTEGER,INTENT(IN)::k
    REAL(KIND(1d0)),INTENT(IN)::xa,xb
    LOGICAL,INTENT(OUT)::zero
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(p)::ws,wg
    REAL(KIND(1d0))::xs,xo,L
    INTEGER::side,i0,j0,a
    ffe_patch_val=0d0
    side=ffe_patch_side(xa,xb)
    zero=side.LE.0
    IF(zero)RETURN
    IF(.NOT.ffe_ZHready(k))CALL ffe_build_patch(k)
    IF(side.EQ.1)THEN
       xs=xa; xo=xb
    ELSE
       xs=xb; xo=xa
    ENDIF
    CALL ffe_lagu((-DLOG10(1d0-xs)-FFE_SH0)/FFE_HH,FFE_NH,p,i0,ws)
    CALL ffe_lagw(ffe_xg,-DLOG10(xo),p,j0,wg)
    zero=ANY(ffe_ZH(i0+1:i0+p,j0+1:j0+p,k,side).LE.-700d0)
    IF(zero)RETURN
    L=0d0
    DO a=1,p
       L=L+ws(a)*SUM(wg*ffe_ZH(i0+a,j0+1:j0+p,k,side))
    ENDDO
    ffe_patch_val=DEXP(L)/(xa*xb)
  END FUNCTION ffe_patch_val

  FUNCTION ffe_N_val(ib,x)
    ! int W dv of beam ib at x for the PNOHAD=1 flux: the 1D table, the patch near x = 1, or directly
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_N_val
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,PARAMETER::p=8
    REAL(KIND(1d0)),DIMENSION(p)::wt
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Wr
    INTEGER::i0
    LOGICAL::dir
    dir=.FALSE.
    IF(x.GE.FFE_XHLO.AND.ffe_patch_beam(ib))THEN
       IF(x.LE.FFE_XHHI)THEN
          CALL ffe_build_patchW()
          CALL ffe_lagu((-DLOG10(1d0-x)-FFE_SH0)/FFE_HH,FFE_NH,p,i0,wt)
          dir=ANY(ffe_lnNH(i0+1:i0+p,ib).LE.-600d0)
          IF(.NOT.dir)ffe_N_val=DEXP(SUM(wt*ffe_lnNH(i0+1:i0+p,ib)))
       ELSE
          dir=.TRUE.
       ENDIF
    ELSE
       CALL ffe_lagw(ffe_xg,-DLOG10(x),p,i0,wt)
       dir=ANY(ffe_lnN(i0+1:i0+p,ib).LE.-600d0)
       IF(.NOT.dir)ffe_N_val=DEXP(SUM(wt*ffe_lnN(i0+1:i0+p,ib)))
    ENDIF
    IF(dir)THEN
       CALL ffe_build_tables()
       ALLOCATE(Wr(ffe_vg(ib)%n))
       CALL ffe_W_row(ib,x,Wr)
       ffe_N_val=DOT_PRODUCT(ffe_vg(ib)%qw,Wr)
       DEALLOCATE(Wr)
    ENDIF
  END FUNCTION ffe_N_val

  FUNCTION ffe_anorm2(x1,x2)
    ! running alpha (pp, pA): [alpha(0)/alpha(mu2)]^2, mu2 = x1 x2 s, the normalisation of the two photons of
    ! set_photon_scale (each photon carries alpha(Q2)^2/alpha(mu2), the hard process alpha(mu2))
    IMPLICIT NONE
    REAL(KIND(1d0))::ffe_anorm2
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    ffe_anorm2=alpha_norm_mu(x1*x2*4d0*ffe_bm(1)%ebeam*ffe_bm(2)%ebeam)**2
  END FUNCTION ffe_anorm2

  LOGICAL FUNCTION ffe_ourW()
    ! the W rows of the proton come from the densities of ElasticFlux (ffe_W_row): the options of the inelastic
    ! branch (exact distribution, dde1/luxqed form factors, running alpha)
    IMPLICIT NONE
    ffe_ourW=ffe_icoll.NE.3.AND.USE_CHARGEFORMFACTOR4PHOTON.AND.(USE_EXACT_ELASTIC4PHOTON.OR.proton_FF_correction.NE.0&
         .OR.alpha_running4photon)
  END FUNCTION ffe_ourW

  LOGICAL FUNCTION ffe_same_beams()
    ! two identical proton beams (same energy, same v-grid) with the W rows of the inelastic branch (ffe_ourW): the
    ! rows of beam 2 are those of beam 1 (FFE's own rows are computed for each beam, as in v2.0)
    IMPLICIT NONE
    ffe_same_beams=.FALSE.
    IF(.NOT.ffe_ourW())RETURN
    IF(.NOT.(ffe_bm(1)%proton.AND.ffe_bm(2)%proton.AND.ffe_bm(1)%chff.AND.ffe_bm(2)%chff))RETURN
    IF(ffe_bm(1)%ebeam.NE.ffe_bm(2)%ebeam.OR.ffe_bm(1)%gamma.NE.ffe_bm(2)%gamma.OR.ffe_bm(1)%mass.NE.ffe_bm(2)%mass)RETURN
    IF(ffe_vg(1)%n.NE.ffe_vg(2)%n)RETURN
    ffe_same_beams=ALL(ffe_vg(1)%b.EQ.ffe_vg(2)%b)
  END FUNCTION ffe_same_beams

  LOGICAL FUNCTION ffe_supported()
    ! pp and pA: the kernel-factorised engine (PhotonPhotonFlux_Algorithm=2), with the elastic proton photon of
    ! any option (EDFF, dipole ChFF, exact distribution, dde1/luxqed form factors, running alpha: ffe_W_row) for
    ! photon_final_state 0; the final states with the inelastic photon (1-3) use the legacy integration
    IMPLICIT NONE
    ffe_supported=PhotonPhotonFlux_Algorithm.EQ.2.AND.photon_final_state.EQ.0
  END FUNCTION ffe_supported


  SUBROUTINE ffe_build_xgrid(umax)
    ! x grid in u = -log10 x on [0, umax]: h = 1/80 for u < 2, h = 1/20 beyond
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::umax
    INTEGER::n1,n2,k,m
    REAL(KIND(1d0))::sc
    sc=MAX(FFE_NodeScale,1d0)
    IF(ALLOCATED(ffe_xg%v))DEALLOCATE(ffe_xg%v,ffe_xg%b,ffe_xg%qw,ffe_xg%slo,ffe_xg%shi,ffe_xg%sv0,ffe_xg%sh)
    n1=NINT(160d0*sc)
    n2=NINT((umax-2d0)*20d0*sc)
    ffe_xg%nseg=2
    ffe_xg%n=1+n1+n2
    ALLOCATE(ffe_xg%v(ffe_xg%n),ffe_xg%b(ffe_xg%n),ffe_xg%qw(ffe_xg%n))
    ALLOCATE(ffe_xg%slo(2),ffe_xg%shi(2),ffe_xg%sv0(2),ffe_xg%sh(2))
    ffe_xg%sv0=(/0d0,2d0/)
    ffe_xg%sh=(/2d0/DBLE(n1),(umax-2d0)/DBLE(n2)/)
    ffe_xg%slo=(/1,1+n1/)
    ffe_xg%shi=(/1+n1,ffe_xg%n/)
    k=0
    DO m=0,n1
       k=k+1
       ffe_xg%v(k)=DBLE(m)*ffe_xg%sh(1)
    ENDDO
    DO m=1,n2
       k=k+1
       ffe_xg%v(k)=2d0+DBLE(m)*ffe_xg%sh(2)
    ENDDO
    ffe_xg%v(ffe_xg%n)=umax
    ffe_xg%b=10d0**(-ffe_xg%v)
    ffe_xg%qw=0d0
  END SUBROUTINE ffe_build_xgrid
END MODULE PhotonPhotonFlux
