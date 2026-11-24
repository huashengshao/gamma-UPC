MODULE InelasticFlux
  ! Inelastic (incoherent) photon distribution of the proton, in b space and collinear, as the Python code
  ! (photon_pdf_components/Inelastic.py, Rho.py; Draft/Magnetic_EPA/incoherent_short.tex):
  !   x f(x,b) = int dmx^2/(2 pi) [distribution at fixed mx],  mx >= mp + mpi the mass of the final state,
  !   D = kT^2 + A^2, A^2 = x (mx^2 - mp^2) + x^2 mp^2, Q2 = D/(1-x), b~ = b/(1-x),
  !   F1 = sqrt(2 pi/(Q2 + mx^2 - mp^2)) (GE_h + tau GM_h)/(1 + tau),
  !   F2 = sqrt(2 pi/(Q2 + mx^2 - mp^2)) (GM_h - GE_h)/(1 + tau),   tau = Q2/(mx + mp)^2,
  !   GE_h^2 = (Q2 + mx^2 - mp^2)^2/(Q2 (Q2 + (mx - mp)^2)) tau FL_h,
  !   GM_h^2 = (1 + tau)/tau F2_h - (Q2 + mx^2 - mp^2)^2/(Q2 (Q2 + (mx - mp)^2)) FL_h,
  ! F2_h, FL_h of structure_functions; GE_h, GM_h tabulated in ln Q2 on the mx nodes (monotone cubic, PCHIP).
  ! Running alpha (alpha_running4photon of ElasticFlux): alpha(Q2)/alpha(0) in F1, F2, alpha(0)/alpha(mu2) overall.
  ! rho (t dependence of the hadronic tensor, rho_model) acts in the continuum (W^2 step of the structure
  ! functions; Python rho_region = 'continuum'), with the source rule across the step.
  ! Normalisation as the elastic densities of ElasticFlux: x f/alpha, b in GeV^-1; the scale mu2 (GeV^2) cuts
  ! kT^2 + x^2 mp^2 < mu2 (required; HUGE(1d0): no cut).
  USE ElasticFlux,ONLY:proton_exact_kspace,proton_exact_edges,proton_exact_comb_P,&
       mproton,alphaem_elasticphoton,aqedup,alpha_running4photon,GENERATE_PROTON_FLUX_GRID
  USE alpha_running,ONLY:alpha_ratio_tab
  USE OpticalGlauber_Geometry,ONLY:ChargeFormFactor_WoodsSaxon
  USE structure_functions,ONLY:F2_FL_h,mpi,W2lo,W2hi,Q2_pdf,Q2_pdf_ramp,rho_W2,low_Q2_model,sf_check,sf_version,&
       sf_pdf_id,sf_xc,sf_Q2c
  USE bessel_integration,ONLY:bessel_integrate_finite,sort_unique
  IMPLICIT NONE
  PRIVATE
  PUBLIC::PhotonNumberDensity_inelastic_proton,PhotonFlux_proton_inelastic_nob
  PUBLIC::xfxQ2_PDFinel_mx,GE_GM_h,coll_row ! at fixed mx (and the form factors), all mu of one x: for checks
  ! rho_model (Rho.py): 0 'delta' (no t dependence), 1 'gaussian', 2 'regge_shrink', 3 'dfjk' (default)
  INTEGER,PUBLIC::rho_model=3
  ! number of collinear k_T integrals whose |K15-G7| exceeded rtol_inel*int|integrand| (the Python recomputes them
  ! with quad)
  INTEGER,PUBLIC::inel_quad_fallbacks=0
  ! The b-space density is built in rows (all b and mu of one x at once) with k_T panels independent of b (edges at
  ! the ln Q2 nodes of the G^h tables, the mu cuts) and the Bessel oscillations integrated per panel (Levin
  ! collocation, Clenshaw-Curtis on narrow panels). rho = delta (lev_row): 1e-8 of a converged reference (the Python's
  ! panels: 3e-6 from it). rho_model /= 0 (rho_row): the rho pairs from momentum-space overlaps. Rows filled on demand
  ! or, all at once, by inel_fast_prepare(_rows) (serial; in parallel only if built with the optional make OPENMP=1).
  PUBLIC::inel_fast_prepare,inel_fast_prepare_rows
  ! two-dimensional Fourier transforms of the b-space densities (momentum space, the (1-P) parts of the final states of
  ! PhotonPhotonFlux)
  PUBLIC::inel_rhohat,elastic_rhohat,inel_rhohat_key,nucleus_rhohat
  LOGICAL,PARAMETER::rhohat_fast=.TRUE. ! rho_phat_fast in inel_rhohat (.FALSE.: rho_phat_mu, the reference)
  REAL(KIND(1d0)),PARAMETER::rhohat_kb2=30d0 ! (M12) bands of cuts rhohat_kb2 D <= kc < rhohat_kb D to O(D^2) (<= 1e-6)
  REAL(KIND(1d0)),PARAMETER::rhohat_kb=300d0 ! bands of cuts kc >= rhohat_kb D at leading order (<= 1.4e-6)
  REAL(KIND(1d0)),PARAMETER::rhohat_kf=3d0 ! cuts kc < rhohat_kf D by rho_phat_mu (coarse outer panels)
  ! disk cache of the rows (fast_row_compute): directory (default '.', the working directory; '' disables). One file
  ! per computed block of a row, named by a hash of the settings (rho_model, low_Q2_model, alpha_running4photon,
  ! sf_pdf_id and its continuation thresholds, CACHE_VERSION), x and the b and mu2 lists; read back only if all of
  ! them agree exactly, written through a temporary file and a rename (safe for parallel jobs sharing the folder).
  ! inel_cache_hits, inel_cache_writes: blocks read from and written to the cache in this run
  CHARACTER(LEN=512),PUBLIC::inel_cache_dir='.'
  INTEGER,PUBLIC::inel_cache_hits=0,inel_cache_writes=0
  INTEGER,PARAMETER::CACHE_VERSION=2 ! 2 (M15): rho_row without the Delta^2 interpolant, 6 Delta panels per decade
  LOGICAL::cache_warned=.FALSE.
  REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0,twopi=2d0*pi
  REAL(KIND(1d0)),PARAMETER::mp=mproton,mp2=mp**2
  ! k_T range and mx nodes (Inelastic.py): upper limit Q2 = Q2_fac (mx^2 + Q2 at kT = 0); nodes dropped when
  ! A b~ > Ab_max (rho = delta rows) or A^2 > A2_mu mu2
  REAL(KIND(1d0)),PARAMETER::Q2_fac=1d3,rtol_inel=1d-2,Ab_max=50d0,A2_mu=1d3
  REAL(KIND(1d0)),PARAMETER::lnQ2_min=-25d0,lnQ2_max=25d0
  INTEGER,PARAMETER::NLQ=369 ! ln Q2 nodes of the tables (smooth models): 40 + 99 + 230
  INTEGER,PARAMETER::n_res=8,n_cont=6
  REAL(KIND(1d0)),PARAMETER::du_mx=1d0,u_max_mx=27.6310211159285482d0 ! ln(1e12)
  ! Gauss-Legendre nodes and weights: 6, 8 points
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::t6=(/-0.932469514203152027812301554494d0,-0.661209386466264513661399595020d0,&
       -0.238619186083196908630501721681d0,0.238619186083196908630501721681d0,0.661209386466264513661399595020d0,&
       0.932469514203152027812301554494d0/)
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::w6=(/0.171324492379170345040296142173d0,0.360761573048138607569833513838d0,&
       0.467913934572691047389870343990d0,0.467913934572691047389870343990d0,0.360761573048138607569833513838d0,&
       0.171324492379170345040296142173d0/)
  REAL(KIND(1d0)),DIMENSION(8),PARAMETER::t8=(/-0.960289856497536231683560868569d0,-0.796666477413626739591553936476d0,&
       -0.525532409916328985817739049189d0,-0.183434642495649804939476142360d0,0.183434642495649804939476142360d0,&
       0.525532409916328985817739049189d0,0.796666477413626739591553936476d0,0.960289856497536231683560868569d0/)
  REAL(KIND(1d0)),DIMENSION(8),PARAMETER::w8=(/0.101228536290376259152531354310d0,0.222381034453374470544355994426d0,&
       0.313706645877887287337962201987d0,0.362683783378361982965150449277d0,0.362683783378361982965150449277d0,&
       0.313706645877887287337962201987d0,0.222381034453374470544355994426d0,0.101228536290376259152531354310d0/)
  ! the six k_T integrals (I_rows of the Python): (form, nk, nj) = (F1,2,J1),(F1,1,J0),(F2,3,J0),(F2,3,J2),(F2,2,J1),(F2,1,J0)
  INTEGER,DIMENSION(6),PARAMETER::row_nj=(/1,0,0,2,1,0/) ! (form and power of k_T: rows_at)
  ! tables: ln Q2 nodes, mx nodes and weights of int dmx^2/(2 pi), GE_h, GM_h and their PCHIP slopes per node
  LOGICAL::init=.FALSE.
  INTEGER::NMX,tab_model=-1 ! low_Q2_model of the tables (rebuilt when it changes, as set_low_Q2_model of the Python)
  INTEGER::tab_version=-1 ! sf_version (PDF set) of the tables (rebuilt when it changes)
  REAL(KIND(1d0)),DIMENSION(NLQ)::lq
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::mx_nodes,mx_weights
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::tab ! (NLQ, 4: GE, GM, dGE, dGM, node)
  ! tables of mx values that are not nodes (last NEXTRA)
  INTEGER,PARAMETER::NEXTRA=16
  ! fast rows: NLV Chebyshev-Lobatto points per panel (differentiation matrices, Clenshaw-Curtis weights),
  ! Clenshaw-Curtis for panels with b~ (width) <= lev_theta, Levin beyond
  INTEGER,PARAMETER::NLV=8
  REAL(KIND(1d0)),PARAMETER::lev_theta=3d0
  ! fast rows with rho (rho_row): Delta panels rho_dens per decade from rho_dlo kappa to Dmax, the overlaps exact on
  ! every panel (M15: no interpolant in Delta^2; 6 panels per decade, 3 left gaussian and regge_shrink 3e-5 of the
  ! peak from refined rows, 6 give <= 4e-7); inner psi integral: one Gauss-Legendre 4 panel when
  ! Delta < rho_inner_short k1, else rho_inner_npan Gauss-Legendre 8 panels (each setting <= 1e-8 of the peak against
  ! tighter ones)
  REAL(KIND(1d0)),PARAMETER::rho_dens=6d0,rho_dlo=1d-4,rho_inner_short=0.3d0
  INTEGER,PARAMETER::rho_inner_npan=2
  ! inner psi integrals of rho_phat_mu (inner_phi_cc): NCL (one panel when Delta < rho_inner_short k1: NCS)
  ! Chebyshev-Lobatto points per panel, values -> Chebyshev coefficients (cml, cms), integrals of T_k on [-1,1] (wtl, wts)
  INTEGER,PARAMETER::NCL=17,NCS=9
  REAL(KIND(1d0)),DIMENSION(NCL)::xcl,wtl,wvl ! wvl, wvs: Clenshaw-Curtis weights of the values
  REAL(KIND(1d0)),DIMENSION(NCS)::xcs,wts,wvs
  REAL(KIND(1d0)),DIMENSION(NCL,NCL)::cml
  REAL(KIND(1d0)),DIMENSION(NCS,NCS)::cms
  REAL(KIND(1d0)),DIMENSION(NLV)::tlv,wcc
  REAL(KIND(1d0)),DIMENSION(NLV)::wbar ! barycentric weights of the Chebyshev-Lobatto points: (-1)^(j-1), halved at the ends
  REAL(KIND(1d0)),DIMENSION(NLV,NLV)::Dlv,D2lv
  LOGICAL::lv_init=.FALSE.
  ! store of the fast rows: per x, the b and mu2 computed and x f/alpha there (last NFROW rows kept)
  INTEGER,PARAMETER::NFROW=256
  TYPE frow_t
     REAL(KIND(1d0))::x=-1d0
     INTEGER::nb=0,nmu=0,key=-1
     REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::b,mu2
     REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::N
  END TYPE frow_t
  TYPE(frow_t),DIMENSION(NFROW),SAVE::frows
  INTEGER::next_frow=0
  REAL(KIND(1d0)),DIMENSION(NEXTRA)::mx_extra=-1d0
  REAL(KIND(1d0)),DIMENSION(NLQ,4,NEXTRA)::tab_extra
  INTEGER::next_extra=0
  ! state of the integrands (x, mx, the k_T weights)
  REAL(KIND(1d0))::xc=0d0,mxc=-1d0,sregc=1d0
  INTEGER::itab=0 ! table of mxc: > 0 node, < 0 extra
  INTEGER::lqhint=0 ! ln Q2 interval of the last GE_GM_core
  !$OMP THREADPRIVATE(xc,mxc,sregc,itab,lqhint)
  ! grid of the collinear flux (with GENERATE_PROTON_FLUX_GRID of ElasticFlux): ln(x f/alpha) of coll_row (cps) on
  ! v = ln(x/(1-x)) (CVH: 80 per decade at small x), CXMIN <= x <= CXMAX, and
  ! t = ln(mu^2 - x^2 mp^2) (CTH: 20 per decade of k_c = sqrt(mu^2 - x^2 mp^2)), CKMIN <= k_c <= CKMAX; a v row (all
  ! t) is computed when first needed (coll_row: all mu at the cost of one); 6-point Lagrange in v and t (with the
  ! steps of the mx-node cut taken out, below); coll_row directly outside the grid or where a value of the stencil is
  ! zero.
  REAL(KIND(1d0)),PARAMETER::CXMIN=1d-7,CXMAX=0.9999d0,CKMIN=1d0,CKMAX=1d5
  ! The flux is a staircase in mu: each mx node is cut at A^2 > A2_mu mu2 (as the Python), a step that grows as
  ! x -> 1 (up to 2% between the t nodes at x = 0.9, mu = 10 TeV). So each node of the grid keeps, besides the number
  ! of active mx nodes (cnact), the flux with k = -KC..KC nodes fewer or more than active (cps, ln; partial sums in
  ! node order, k = 0 the value itself), and the interpolation takes at every node of the stencil the value with the
  ! active nodes of the point (nactive): no steps inside the stencil (grid/direct-1 <= 2e-6, 1e-6 <= x <= 0.9999,
  ! 5 GeV <= mu <= 10 TeV). 12 MB, allocated at the first use.
  INTEGER,PARAMETER::KC=8
  INTEGER,DIMENSION(:,:),ALLOCATABLE::cnact
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::cps
  REAL(KIND(1d0)),PARAMETER::CVH=2.30258509299404568d0/80d0,CTH=2d0*2.30258509299404568d0/20d0
  INTEGER,PARAMETER::NCV=886,NCT=102 ! v(NCV-2) >= ln(CXMAX/(1-CXMAX)), t(NCT) >= 2 ln(CKMAX)
  LOGICAL,DIMENSION(NCV)::cdone=.FALSE.
  INTEGER::ckey=-1
CONTAINS
  ! ================================================================ public: b space and collinear
  FUNCTION PhotonNumberDensity_inelastic_proton(b,Ega,gamma,mu2)
    ! x f(x,b)/alpha of the inelastic photon, x = Ega/(gamma mp), b in GeV^-1 (as the elastic densities)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_inelastic_proton
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    REAL(KIND(1d0)),INTENT(IN)::mu2 ! GeV^2 (HUGE(1d0): no cut, the Python default; for checks only)
    REAL(KIND(1d0))::x,m2
    PhotonNumberDensity_inelastic_proton=0d0
    x=Ega/(gamma*mp)
    IF(x.GE.1d0)RETURN
    m2=mu2
    PhotonNumberDensity_inelastic_proton=fast_density(x,b,m2)
  END FUNCTION PhotonNumberDensity_inelastic_proton

  FUNCTION PhotonFlux_proton_inelastic_nob(x,mu2)
    ! dN/dx of the inelastic photon (int d^2b, P_NOHAD = 1) at the scale mu2 (GeV^2, needed: it grows as ln mu2)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_inelastic_nob
    REAL(KIND(1d0)),INTENT(IN)::x,mu2
    REAL(KIND(1d0))::alpha,F(1)
    alpha=alphaem_elasticphoton
    IF(alpha.LT.0d0)alpha=MERGE(aqedup,0.0072992701d0,aqedup.GT.0d0)
    PhotonFlux_proton_inelastic_nob=0d0
    IF(x.GE.1d0)RETURN
    IF(.NOT.coll_grid(x,mu2,F(1)))CALL coll_row(x,1,(/mu2/),F)
    PhotonFlux_proton_inelastic_nob=alpha/x*F(1)
  END FUNCTION PhotonFlux_proton_inelastic_nob

  FUNCTION coll_grid(x,mu2,F)
    ! x f/alpha (collinear) from the grid (cps); .FALSE. when not available (grids off, outside, zero values)
    IMPLICIT NONE
    LOGICAL::coll_grid
    REAL(KIND(1d0)),INTENT(IN)::x,mu2
    REAL(KIND(1d0)),INTENT(OUT)::F
    REAL(KIND(1d0)),PARAMETER::v0=-16.1180955509583157d0 ! ln(CXMIN/(1-CXMIN))
    REAL(KIND(1d0))::kc2,v,t,wv(6),wt(6),mu2v(NCT),row(NCT),xi,ps(-KC:KC,NCT),z(6,6)
    INTEGER::i0,j0,i,j,key,na(NCT),ns,k
    coll_grid=.FALSE.; F=0d0
    IF(.NOT.GENERATE_PROTON_FLUX_GRID.OR.mu2.GE.HUGE(1d0).OR.x.LT.CXMIN.OR.x.GT.CXMAX)RETURN
    kc2=mu2-x**2*mp2
    IF(kc2.LT.CKMIN**2.OR.kc2.GT.CKMAX**2)RETURN
    CALL inel_init()
    key=low_Q2_model+10*MERGE(1,0,alpha_running4photon)+100*tab_version
    IF(key.NE.ckey)THEN
       cdone=.FALSE.; ckey=key
    ENDIF
    v=DLOG(x/(1d0-x)); t=DLOG(kc2)
    i0=MIN(MAX(FLOOR((v-v0)/CVH)-1,1),NCV-5)
    j0=MIN(MAX(FLOOR((t-2d0*DLOG(CKMIN))/CTH)-1,1),NCT-5)
    DO i=i0,i0+5
       IF(cdone(i))CYCLE
       xi=1d0/(1d0+DEXP(-(v0+(i-1)*CVH)))
       DO j=1,NCT
          mu2v(j)=DEXP(2d0*DLOG(CKMIN)+(j-1)*CTH)+xi**2*mp2
       ENDDO
       IF(.NOT.ALLOCATED(cps))ALLOCATE(cnact(NCT,NCV),cps(-KC:KC,NCT,NCV))
       CALL coll_row(xi,NCT,mu2v,row,na,ps)
       cnact(:,i)=na
       cps(:,:,i)=-HUGE(1d0)
       WHERE(ps.GT.0d0)cps(:,:,i)=DLOG(ps)
       cdone(i)=.TRUE.
    ENDDO
    ! the values at the nodes with the active mx nodes of the point (ns): no steps inside the stencil
    ns=nactive(x,mu2)
    DO i=1,6
       DO j=1,6
          k=ns-cnact(j0+j-1,i0+i-1)
          IF(ABS(k).GT.KC)RETURN
          z(j,i)=cps(k,j0+j-1,i0+i-1)
       ENDDO
    ENDDO
    IF(ANY(z.EQ.-HUGE(1d0)))RETURN
    wv=lagrange6((v-v0)/CVH-(i0-1))
    wt=lagrange6((t-2d0*DLOG(CKMIN))/CTH-(j0-1))
    F=DEXP(DOT_PRODUCT(wv,MATMUL(wt,z)))
    coll_grid=.TRUE.
  CONTAINS
    FUNCTION lagrange6(p) ! weights of the 6 nodes 0..5 (unit spacing) at p
      REAL(KIND(1d0)),DIMENSION(6)::lagrange6
      REAL(KIND(1d0)),INTENT(IN)::p
      INTEGER::k,l
      lagrange6=1d0
      DO k=1,6
         DO l=1,6
            IF(l.NE.k)lagrange6(k)=lagrange6(k)*(p-(l-1))/DBLE(k-l)
         ENDDO
      ENDDO
    END FUNCTION lagrange6
  END FUNCTION coll_grid

  INTEGER FUNCTION nactive(x,mu2) ! number of mx nodes with A^2 <= A2_mu mu2 (the first ones: A^2 grows with mx)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,mu2
    INTEGER::lo,hi,mid
    lo=0; hi=NMX
    DO WHILE(hi.GT.lo)
       mid=(lo+hi+1)/2
       IF(DSQRT(A2(x,mx_nodes(mid)))**2.GT.A2_mu*mu2)THEN ! the test of coll_row (A**2, A = sqrt(A2))
          hi=mid-1
       ELSE
          lo=mid
       ENDIF
    ENDDO
    nactive=lo
  END FUNCTION nactive

  SUBROUTINE coll_row(x,nmu,mu2v,F,NACT,PS)
    ! collinear x f/alpha (int d^2b) at one x and all mu2v (GeV^2; HUGE(1d0): no cut): int dmx^2/(2 pi) on the mx
    ! nodes (mx nodes with A^2 > A2_mu mu2 dropped per mu, as mx_integral of the Python) of 2 pi int dkT/kT of the
    ! distribution with the integrals replaced by their integrands (xfxQ2_PDFinel_mx). Per mx node the k_T panels of
    ! lev_row (geometric 10 per decade from 1e-3 A, the k of the ln Q2 nodes of the G^h tables, the mu cuts, the
    ! upper limit), the integrand at the NLV Chebyshev-Lobatto points with the interpolatory weights wcc, summed up
    ! to each mu cut: all mu at the cost of one.
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nmu
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(OUT)::F
    ! optional (coll_grid): per mu the number of active mx nodes and the flux with k = -KC..KC nodes fewer or more
    ! (partial sums in node order; k = 0: F; nodes beyond the last one computed: -1)
    INTEGER,DIMENSION(nmu),INTENT(OUT),OPTIONAL::NACT
    REAL(KIND(1d0)),DIMENSION(-KC:KC,nmu),INTENT(OUT),OPTIONAL::PS
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::cnode
    INTEGER::last,n
    REAL(KIND(1d0))::cs
    REAL(KIND(1d0)),DIMENSION(nmu)::kmu,Im
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::e
    REAL(KIND(1d0)),DIMENSION(6)::v
    REAL(KIND(1d0))::A,hi0,klo,l0,l1,c0,h,kk,Q2,Icum,pk
    INTEGER::i,m,ngeo,ne,np,k,jj
    CALL inel_init()
    IF(.NOT.lv_init)CALL levin_setup()
    DO m=1,nmu
       kmu(m)=HUGE(1d0)
       IF(mu2v(m).LT.HUGE(1d0))kmu(m)=DSQRT(MAX(mu2v(m)-x**2*mp2,0d0))
    ENDDO
    F=0d0
    ALLOCATE(cnode(NMX,MERGE(nmu,1,PRESENT(PS)))); cnode=0d0
    last=0
    DO i=1,NMX
       A=DSQRT(A2(x,mx_nodes(i)))
       IF(ALL(mu2v.LT.HUGE(1d0).AND.A**2.GT.A2_mu*mu2v))CYCLE
       last=i
       hi0=DSQRT(Q2_fac*((1d0-x)*mx_nodes(i)**2+A**2))
       CALL set_state(x,mx_nodes(i),1d0)
       klo=1d-3*MIN(A,hi0)
       ngeo=INT(10*LOG10(hi0/klo))+2
       ALLOCATE(e(ngeo+NLQ+nmu+4))
       l0=LOG10(klo); l1=LOG10(hi0)
       DO k=0,ngeo-1
          e(k+1)=10d0**(l0+k*((l1-l0)/(ngeo-1)))
       ENDDO
       ne=ngeo+1; e(ne)=0d0
       DO k=1,NLQ
          Q2=DEXP(lq(k))*(1d0-x)-A**2
          IF(Q2.GT.0d0)THEN
             IF(DSQRT(Q2).LT.hi0)THEN
                ne=ne+1; e(ne)=DSQRT(Q2)
             ENDIF
          ENDIF
       ENDDO
       DO m=1,nmu
          IF(kmu(m).LT.hi0)THEN
             ne=ne+1; e(ne)=kmu(m)
          ENDIF
       ENDDO
       ne=ne+1; e(ne)=hi0
       CALL sort_unique(e,ne,np)
       Icum=0d0; Im=0d0
       DO k=1,np-1
          DO m=1,nmu
             IF(kmu(m).EQ.e(k))Im(m)=Icum
          ENDDO
          c0=0.5d0*(e(k)+e(k+1)); h=0.5d0*(e(k+1)-e(k))
          pk=0d0
          DO jj=1,NLV
             kk=c0+h*tlv(jj)
             IF(kk.LE.0d0)CYCLE ! the integrand vanishes at kT = 0
             CALL rows_at(kk,v)
             pk=pk+h*wcc(jj)*combination_v(v,x,mx_nodes(i))/kk
          ENDDO
          Icum=Icum+pk
       ENDDO
       DO m=1,nmu
          IF(kmu(m).GE.hi0)Im(m)=Icum
          IF(PRESENT(PS))cnode(i,m)=mx_weights(i)*Im(m) ! (only with PS: cnode has nmu columns then)
          IF(mu2v(m).LT.HUGE(1d0).AND.A**2.GT.A2_mu*mu2v(m))CYCLE
          F(m)=F(m)+mx_weights(i)*Im(m)
       ENDDO
       DEALLOCATE(e)
    ENDDO
    IF(PRESENT(PS))THEN
       DO m=1,nmu
          n=nactive(x,mu2v(m))
          NACT(m)=n
          PS(:,m)=-1d0
          cs=0d0
          DO i=1,MIN(n+KC,last)
             cs=cs+cnode(i,m)
             IF(i.GE.n-KC)PS(i-n,m)=cs*twopi*alpha_norm(mu2v(m))/pi**2
          ENDDO
          DO k=-KC,MIN(-1,-n) ! fewer than no node
             PS(k,m)=-1d0
          ENDDO
          IF(n.LE.KC.AND.n-KC.LE.0)PS(-n,m)=0d0
       ENDDO
    ENDIF
    DEALLOCATE(cnode)
    DO m=1,nmu
       F(m)=F(m)*twopi*alpha_norm(mu2v(m))/pi**2
    ENDDO
  END SUBROUTINE coll_row

  ! ================================================================ distribution at fixed mx
  FUNCTION combination_v(v,x,mx)
    ! combination of the pairs P(i,j) = v(i)*v(j) (the same arithmetic, without the 6x6 products)
    IMPLICIT NONE
    REAL(KIND(1d0))::combination_v
    REAL(KIND(1d0)),DIMENSION(6),INTENT(IN)::v
    REAL(KIND(1d0)),INTENT(IN)::x,mx
    REAL(KIND(1d0))::c
    c=mx-(1d0-x)*mp
    combination_v=(1d0+(1d0-x)**2)/2d0*(v(1)*v(1))+((v(3)*v(3))+(v(4)*v(4)))/2d0/(mx+mp)**2&
         +x*(c/(mx+mp)*((v(3)*v(6))+(v(2)*v(3))+x*(mx+mp)*c*(v(2)*v(6)))-((1d0-x)*mx-mp)/(mx+mp)*(v(1)*v(5)))&
         +x**2/2d0*(c**2*((v(2)*v(2))+(v(6)*v(6)))+(mx**2+mp2)/(mx+mp)**2*(v(5)*v(5)))
  END FUNCTION combination_v

  FUNCTION combination(P,x,mx)
    ! the quadratic form of the distribution in the pairs P(i,j) of the six integrals (b-space rows and collinear)
    IMPLICIT NONE
    REAL(KIND(1d0))::combination
    REAL(KIND(1d0)),DIMENSION(6,6),INTENT(IN)::P
    REAL(KIND(1d0)),INTENT(IN)::x,mx
    REAL(KIND(1d0))::c
    c=mx-(1d0-x)*mp
    combination=(1d0+(1d0-x)**2)/2d0*P(1,1)+(P(3,3)+P(4,4))/2d0/(mx+mp)**2&
         +x*(c/(mx+mp)*(P(3,6)+P(2,3)+x*(mx+mp)*c*P(2,6))-((1d0-x)*mx-mp)/(mx+mp)*P(1,5))&
         +x**2/2d0*(c**2*(P(2,2)+P(6,6))+(mx**2+mp2)/(mx+mp)**2*P(5,5))
  END FUNCTION combination

  SUBROUTINE rows_at(kT,v)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::kT
    REAL(KIND(1d0)),DIMENSION(6),INTENT(OUT)::v
    REAL(KIND(1d0))::D,Q2,F(2),k2,k3
    D=kT**2+A2(xc,mxc)
    Q2=D/(1d0-xc)
    CALL F1F2(mxc,Q2,F(1),F(2))
    ! v(i) = kT**row_nk(i)/D*F(row_form(i)), the powers as the integer power of the compiler (kT*(kT*kT) for 3)
    k2=kT*kT; k3=kT*k2
    v(1)=k2/D*F(1); v(2)=kT/D*F(1); v(3)=k3/D*F(2); v(4)=k3/D*F(2); v(5)=k2/D*F(2); v(6)=kT/D*F(2)
  END SUBROUTINE rows_at

  FUNCTION xfxQ2_PDFinel_mx(x,mx,mu2)
    ! collinear x f/alpha at fixed mx: 2 pi (1-x)^2 int dkT/kT of the distribution with the integrals replaced by
    ! their integrands (xfxQ2_PDFinel_mx of the Python; independent of rho_model, as there)
    IMPLICIT NONE
    REAL(KIND(1d0))::xfxQ2_PDFinel_mx
    REAL(KIND(1d0)),INTENT(IN)::x,mx,mu2
    REAL(KIND(1d0)),DIMENSION(1)::I,err,absres
    CALL inel_init()
    CALL set_state(x,mx,1d0)
    CALL bessel_integrate_finite(collinear_g,1,(/-1/),0d0,DSQRT(A2(x,mx)),kT_hi(x,mx,mu2),kT_edges(x,mx),&
         I,err,absres)
    IF(err(1).GT.rtol_inel*absres(1))inel_quad_fallbacks=inel_quad_fallbacks+1
    xfxQ2_PDFinel_mx=twopi*(1d0-x)**2*alpha_norm(mu2)/pi**2/(1d0-x)**2*I(1)
  END FUNCTION xfxQ2_PDFinel_mx

  SUBROUTINE collinear_g(kT,n,v)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::kT
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::v
    REAL(KIND(1d0)),DIMENSION(6)::h
    REAL(KIND(1d0)),DIMENSION(6,6)::P
    INTEGER::k
    CALL rows_at(kT,h)
    DO k=1,6
       P(:,k)=h*h(k)
    ENDDO
    v(1)=combination(P,xc,mxc)/kT
  END SUBROUTINE collinear_g

  ! ================================================================ form factors
  SUBROUTINE F1F2(mx,Q2,F1,F2)
    ! Dirac and Pauli form factors at the amplitude level (times alpha(Q2)/alpha(0) with running alpha)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::mx,Q2
    REAL(KIND(1d0)),INTENT(OUT)::F1,F2
    REAL(KIND(1d0))::GE,GM,tau,c
    CALL GE_GM_core(mx,Q2,GE,GM)
    tau=Q2/(mx+mp)**2
    c=DSQRT(twopi/(Q2+mx**2-mp2))/(1d0+tau)
    IF(alpha_running4photon)c=c*alpha_ratio_tab(Q2) ! (M14) the table of alpha_running (= alpha_ratio_ph to 1.8e-11)
    F1=c*(GE+tau*GM)
    F2=c*(GM-GE)
  END SUBROUTINE F1F2

  SUBROUTINE GE_GM_h(mx,Q2,GE,GM)
    ! GE_h, GM_h from the table of mx (built on first use), PCHIP in ln Q2 clipped to the table
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::mx,Q2
    REAL(KIND(1d0)),INTENT(OUT)::GE,GM
    CALL inel_init()
    CALL GE_GM_core(mx,Q2,GE,GM)
  END SUBROUTINE GE_GM_h

  SUBROUTINE GE_GM_core(mx,Q2,GE,GM)
    ! GE_GM_h after inel_init: the interval of the last call (lqhint) tried first, else the bisection (the same
    ! interval as the bisection in both cases)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::mx,Q2
    REAL(KIND(1d0)),INTENT(OUT)::GE,GM
    REAL(KIND(1d0))::v,t,h,m
    INTEGER::i,lo,hi,mid,k
    IF(mx.NE.mxc.OR.itab.EQ.0)CALL find_table(mx)
    v=MIN(MAX(DLOG(Q2),lnQ2_min),lnQ2_max)
    lo=0
    IF(lqhint.GE.1.AND.lqhint.LT.NLQ)THEN ! lq(lo) <= v < lq(lo+1), as the bisection finds it
       IF(v.GE.lq(lqhint).AND.v.LT.lq(lqhint+1))lo=lqhint
    ENDIF
    IF(lo.EQ.0)THEN
       lo=1; hi=NLQ
       DO WHILE(hi-lo.GT.1) ! lq(lo) <= v < lq(hi)
          mid=(lo+hi)/2
          IF(v.GE.lq(mid))THEN
             lo=mid
          ELSE
             hi=mid
          ENDIF
       ENDDO
       lqhint=lo
    ENDIF
    i=MIN(lo,NLQ-1)
    t=v-lq(i); h=lq(i+1)-lq(i)
    DO k=1,2
       IF(itab.GT.0)THEN
          CALL hermite(tab(i,k,itab),tab(i+1,k,itab),tab(i,k+2,itab),tab(i+1,k+2,itab),h,t,m)
       ELSE
          CALL hermite(tab_extra(i,k,-itab),tab_extra(i+1,k,-itab),tab_extra(i,k+2,-itab),tab_extra(i+1,k+2,-itab),h,t,m)
       ENDIF
       IF(k.EQ.1)THEN
          GE=m
       ELSE
          GM=m
       ENDIF
    ENDDO
  END SUBROUTINE GE_GM_core

  SUBROUTINE hermite(y0,y1,d0,d1,h,t,y) ! cubic Hermite on [0,h] at t (scipy PPoly coefficients)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::y0,y1,d0,d1,h,t
    REAL(KIND(1d0)),INTENT(OUT)::y
    REAL(KIND(1d0))::s,c2,c3
    s=(y1-y0)/h
    c2=(3d0*s-2d0*d0-d1)/h
    c3=(d0+d1-2d0*s)/h**2
    y=((c3*t+c2)*t+d0)*t+y0
  END SUBROUTINE hermite

  SUBROUTINE find_table(mx)
    ! sets mxc, itab for mx: a node (itab > 0) or an extra table (itab < 0, built if missing)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::mx
    INTEGER::i
    mxc=mx
    DO i=1,NMX
       IF(mx_nodes(i).EQ.mx)THEN
          itab=i
          RETURN
       ENDIF
    ENDDO
    DO i=1,NEXTRA
       IF(mx_extra(i).EQ.mx)THEN
          itab=-i
          RETURN
       ENDIF
    ENDDO
    next_extra=MOD(next_extra,NEXTRA)+1
    mx_extra(next_extra)=mx
    CALL build_table(mx,tab_extra(:,:,next_extra))
    itab=-next_extra
  END SUBROUTINE find_table

  SUBROUTINE build_table(mx,t)
    ! GE_h, GM_h at the ln Q2 nodes (sqrt of the squares clipped at 0; NaN -> 0) and the PCHIP slopes
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::mx
    REAL(KIND(1d0)),DIMENSION(NLQ,4),INTENT(OUT)::t
    REAL(KIND(1d0))::Q2,F2,FL,r,tau,GE2,GM2
    INTEGER::i
    DO i=1,NLQ
       Q2=DEXP(lq(i))
       CALL F2_FL_h(Q2/(Q2+mx**2-mp2),Q2,F2,FL)
       r=(Q2+mx**2-mp2)**2/(Q2*(Q2+(mx-mp)**2))
       tau=Q2/(mx+mp)**2
       GE2=r*tau*FL
       GM2=(1d0+tau)/tau*F2-r*FL
       IF(GE2.NE.GE2)GE2=0d0
       IF(GM2.NE.GM2)GM2=0d0
       t(i,1)=DSQRT(MAX(GE2,0d0))
       t(i,2)=DSQRT(MAX(GM2,0d0))
    ENDDO
    CALL pchip_slopes(lq,t(:,1),t(:,3))
    CALL pchip_slopes(lq,t(:,2),t(:,4))
  END SUBROUTINE build_table

  SUBROUTINE pchip_slopes(x,y,d) ! derivatives of scipy PchipInterpolator (Fritsch-Butland, scipy edge rule)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::x,y
    REAL(KIND(1d0)),DIMENSION(:),INTENT(OUT)::d
    REAL(KIND(1d0)),DIMENSION(SIZE(x)-1)::h,m
    REAL(KIND(1d0))::w1,w2
    INTEGER::n,k
    n=SIZE(x)
    h=x(2:n)-x(1:n-1)
    m=(y(2:n)-y(1:n-1))/h
    d=0d0
    DO k=2,n-1
       IF(SIGN(1d0,m(k-1)).NE.SIGN(1d0,m(k)).OR.m(k).EQ.0d0.OR.m(k-1).EQ.0d0)CYCLE
       w1=2d0*h(k)+h(k-1)
       w2=h(k)+2d0*h(k-1)
       d(k)=1d0/((w1/m(k-1)+w2/m(k))/(w1+w2))
    ENDDO
    d(1)=edge(h(1),h(2),m(1),m(2))
    d(n)=edge(h(n-1),h(n-2),m(n-1),m(n-2))
  CONTAINS
    FUNCTION edge(h0,h1,m0,m1)
      REAL(KIND(1d0))::edge
      REAL(KIND(1d0)),INTENT(IN)::h0,h1,m0,m1
      edge=((2d0*h0+h1)*m0-h0*m1)/(h0+h1)
      IF(sgn(edge).NE.sgn(m0))THEN
         edge=0d0
      ELSEIF(sgn(m0).NE.sgn(m1).AND.ABS(edge).GT.3d0*ABS(m0))THEN
         edge=3d0*m0
      ENDIF
    END FUNCTION edge
    FUNCTION sgn(v) ! np.sign
      REAL(KIND(1d0))::sgn
      REAL(KIND(1d0)),INTENT(IN)::v
      sgn=0d0
      IF(v.GT.0d0)sgn=1d0
      IF(v.LT.0d0)sgn=-1d0
    END FUNCTION sgn
  END SUBROUTINE pchip_slopes

  ! ================================================================ kinematics
  FUNCTION A2(x,mx)
    REAL(KIND(1d0))::A2
    REAL(KIND(1d0)),INTENT(IN)::x,mx
    A2=x*(mx**2-mp2)+x**2*mp2
  END FUNCTION A2

  FUNCTION kT_hi(x,mx,mu2) ! min(sqrt(Q2_fac ((1-x) mx^2 + A^2)), the mu cut)
    REAL(KIND(1d0))::kT_hi
    REAL(KIND(1d0)),INTENT(IN)::x,mx,mu2
    kT_hi=DSQRT(Q2_fac*((1d0-x)*mx**2+A2(x,mx)))
    IF(mu2.LT.HUGE(1d0))kT_hi=MIN(kT_hi,DSQRT(MAX(mu2-x**2*mp2,0d0)))
  END FUNCTION kT_hi

  FUNCTION kT_edges(x,mx) ! kT at the ends of the Q2_pdf window, for mx^2 > W2lo
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::kT_edges
    REAL(KIND(1d0)),INTENT(IN)::x,mx
    REAL(KIND(1d0))::k2(2)
    k2=-1d0
    IF(mx**2.GT.W2lo)k2=(1d0-x)*(/Q2_pdf/Q2_pdf_ramp,Q2_pdf*Q2_pdf_ramp/)-A2(x,mx)
    kT_edges=DSQRT(PACK(k2,k2.GT.0d0))
  END FUNCTION kT_edges

  FUNCTION alpha_norm(mu2) ! alpha(0)/alpha(mu2) with running alpha and finite mu2, else 1
    REAL(KIND(1d0))::alpha_norm
    REAL(KIND(1d0)),INTENT(IN)::mu2
    alpha_norm=1d0
    IF(alpha_running4photon.AND.mu2.LT.HUGE(1d0))alpha_norm=1d0/alpha_ratio_tab(mu2)
  END FUNCTION alpha_norm

  SUBROUTINE set_state(x,mx,sreg)
    REAL(KIND(1d0)),INTENT(IN)::x,mx,sreg
    xc=x; sregc=sreg
    IF(mx.NE.mxc.OR.itab.EQ.0)CALL find_table(mx)
  END SUBROUTINE set_state

  ! ================================================================ initialisation: nodes and tables
  SUBROUTINE inel_init()
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(8),PARAMETER::res_edges=(/0d0,1.15d0,1.30d0,1.45d0,1.60d0,0d0,1.87d0,0d0/)
    REAL(KIND(1d0))::e(8),step,u0,u
    INTEGER::i,j,k,npan
    CALL sf_check()
    IF(init)THEN
       IF(low_Q2_model.EQ.tab_model.AND.sf_version.EQ.tab_version)RETURN
       ! new low-Q2 model or PDF set: the tables and the extra tables are rebuilt
       tab_model=low_Q2_model; tab_version=sf_version
       DO k=1,NMX
          CALL build_table(mx_nodes(k),tab(:,:,k))
       ENDDO
       mx_extra=-1d0; next_extra=0
       mxc=-1d0; itab=0
       RETURN
    ENDIF
    init=.TRUE.
    tab_model=low_Q2_model; tab_version=sf_version
    ! ln Q2: arange(-25,-5,0.5), linspace(-5,ln Q2_pdf,100) without its last point, linspace(ln Q2_pdf,25,230)
    DO i=0,39
       lq(i+1)=-25d0+0.5d0*i
    ENDDO
    step=(DLOG(Q2_pdf)+5d0)/99d0
    DO i=0,98
       lq(41+i)=i*step-5d0
    ENDDO
    step=(lnQ2_max-DLOG(Q2_pdf))/229d0
    DO i=0,228
       lq(140+i)=i*step+DLOG(Q2_pdf)
    ENDDO
    lq(NLQ)=lnQ2_max
    ! mx nodes: resonances, Gauss-Legendre 8 on the panels of res_edges; continuum, Gauss-Legendre 6 in
    ! u = ln(mx^2 - mp^2), panels of width du_mx from ln(W2hi - mp^2) beyond ln(1e12)
    e=res_edges; e(1)=mp+mpi; e(6)=DSQRT(W2lo); e(8)=DSQRT(W2hi)
    u0=DLOG(W2hi-mp2)
    npan=CEILING((u_max_mx-u0)/du_mx)
    NMX=7*n_res+npan*n_cont
    ALLOCATE(mx_nodes(NMX),mx_weights(NMX),tab(NLQ,4,NMX))
    k=0
    DO i=1,7
       DO j=1,n_res
          k=k+1
          mx_nodes(k)=0.5d0*(e(i)+e(i+1))+0.5d0*(e(i+1)-e(i))*t8(j)
          mx_weights(k)=0.5d0*(e(i+1)-e(i))*w8(j)*2d0*mx_nodes(k)/twopi
       ENDDO
    ENDDO
    DO i=0,npan-1
       DO j=1,n_cont
          k=k+1
          u=(u0+du_mx*i)+0.5d0*du_mx*(1d0+t6(j))
          mx_nodes(k)=DSQRT(DEXP(u)+mp2)
          mx_weights(k)=0.5d0*du_mx*w6(j)*DEXP(u)/twopi
       ENDDO
    ENDDO
    DO k=1,NMX
       CALL build_table(mx_nodes(k),tab(:,:,k))
    ENDDO
    mxc=-1d0; itab=0
  END SUBROUTINE inel_init

  ! ================================================================ rho
  FUNCTION B_slope(xbj) ! rho = exp(-B Delta^2) (Rho.py: gaussian, regge_shrink, dfjk with flavour u)
    REAL(KIND(1d0))::B_slope
    REAL(KIND(1d0)),INTENT(IN)::xbj
    SELECT CASE(rho_model)
    CASE(1)
       B_slope=2d0
    CASE(2)
       B_slope=(2d0+0.25d0*DLOG(1d0/xbj))*(1d0-xbj)
    CASE DEFAULT
       B_slope=0.9d0*(1d0-xbj)**3*DLOG(1d0/xbj)+0.59d0*(1d0-xbj)**3+1.22d0*xbj*(1d0-xbj)**2
    END SELECT
  END FUNCTION B_slope

  FUNCTION W_step(mx) ! W^2 step of the structure functions (0 resonances, 1 continuum)
    REAL(KIND(1d0))::W_step
    REAL(KIND(1d0)),INTENT(IN)::mx
    IF(mx**2.LE.W2lo)THEN
       W_step=0d0
    ELSEIF(mx**2.GE.W2hi)THEN
       W_step=1d0
    ELSE
       W_step=rho_W2(MIN(MAX(mx**2,W2lo),W2hi))
    ENDIF
  END FUNCTION W_step

  FUNCTION xB(mx,Q2)
    REAL(KIND(1d0))::xB
    REAL(KIND(1d0)),INTENT(IN)::mx,Q2
    xB=Q2/(Q2+mx**2-mp2)
  END FUNCTION xB

  ! ================================================================ fast rows (the b-space density)
  FUNCTION fast_density(x,b,mu2)
    ! x f/alpha at (x, b, mu2) from the store, computed if missing
    IMPLICIT NONE
    REAL(KIND(1d0))::fast_density
    REAL(KIND(1d0)),INTENT(IN)::x,b,mu2
    INTEGER::ir,ib,im
    INTEGER::last=0
    SAVE last
    ! the row of the last call first
    ir=0
    IF(last.GT.0)THEN
       IF(near(frows(last)%x,x).AND.frows(last)%key.EQ.fast_key())ir=last
    ENDIF
    IF(ir.EQ.0)ir=find_frow(x)
    ib=0; im=0
    IF(ir.GT.0)THEN
       ib=find_near(frows(ir)%b,frows(ir)%nb,b)
       im=find_near(frows(ir)%mu2,frows(ir)%nmu,mu2)
    ENDIF
    IF(ib.EQ.0.OR.im.EQ.0)THEN
       CALL inel_fast_prepare(x,1,(/b/),1,(/mu2/))
       ir=find_frow(x)
       ib=find_near(frows(ir)%b,frows(ir)%nb,b)
       im=find_near(frows(ir)%mu2,frows(ir)%nmu,mu2)
    ENDIF
    last=ir
    fast_density=frows(ir)%N(ib,im)
  END FUNCTION fast_density

  LOGICAL FUNCTION near(a,c) ! the same x, b or mu2 of the store up to rounding (1e-13 relative)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::a,c
    near=ABS(a-c).LE.1d-13*MAX(ABS(a),ABS(c))
  END FUNCTION near

  INTEGER FUNCTION find_near(v,n,a) ! first i <= n with near(v(i),a), else 0
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::v
    REAL(KIND(1d0)),INTENT(IN)::a
    INTEGER::i
    find_near=0
    DO i=1,n
       IF(near(v(i),a))THEN
          find_near=i
          RETURN
       ENDIF
    ENDDO
  END FUNCTION find_near

  FUNCTION find_frow(x) ! the row of x (0 if none, or if made with other settings)
    IMPLICIT NONE
    INTEGER::find_frow
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::i
    find_frow=0
    DO i=1,NFROW
       IF(near(frows(i)%x,x).AND.frows(i)%key.EQ.fast_key())THEN
          find_frow=i
          RETURN
       ENDIF
    ENDDO
  END FUNCTION find_frow

  INTEGER FUNCTION inel_rhohat_key() ! (M11) the settings inel_rhohat depends on (as the rows)
    IMPLICIT NONE
    inel_rhohat_key=fast_key()
  END FUNCTION inel_rhohat_key

  FUNCTION fast_key() ! the settings the rows depend on
    IMPLICIT NONE
    INTEGER::fast_key
    CALL inel_init()
    fast_key=low_Q2_model+10*MERGE(1,0,alpha_running4photon)+1000*rho_model+10000*tab_version
  END FUNCTION fast_key

  SUBROUTINE inel_fast_prepare(x,nb,bv,nmu,mu2v)
    ! makes sure the row of x holds x f/alpha at all b in bv (GeV^-1) and mu2 in mu2v (GeV^2, HUGE: no cut)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nb,nmu
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    CALL inel_fast_prepare_rows(1,(/x/),nb,(/nb/),RESHAPE(bv,(/nb,1/)),nmu,(/nmu/),RESHAPE(mu2v,(/nmu,1/)))
  END SUBROUTINE inel_fast_prepare

  SUBROUTINE inel_fast_prepare_rows(nx,xv,ldb,nbv,bv,ldm,nmv,mu2v)
    ! as inel_fast_prepare for nx values of x: row of xv(i) at the b bv(1:nbv(i),i) and the mu2 mu2v(1:nmv(i),i) (equal
    ! x merged). Per row, missing mu2 are computed at the b already there and missing b at all mu2 (fast_row_compute);
    ! the rows are computed in parallel with OpenMP (make OPENMP=1), the smallest x (the slowest) first, at most
    ! NFROW/2 rows at a time
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nx,ldb,ldm
    REAL(KIND(1d0)),DIMENSION(nx),INTENT(IN)::xv
    INTEGER,DIMENSION(nx),INTENT(IN)::nbv,nmv
    REAL(KIND(1d0)),DIMENSION(ldb,nx),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(ldm,nx),INTENT(IN)::mu2v
    TYPE task_t
       INTEGER::ir=0,ob=0,om=0,nnb=0,nnm=0
       REAL(KIND(1d0))::x=0d0
       REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Nm,Nb
    END TYPE task_t
    TYPE(task_t),DIMENSION(:),ALLOCATABLE::tk
    INTEGER,DIMENSION(nx)::ord
    INTEGER::i,j,k,c0,c1,nt,t,ir,ob,om,nnb,nnm
    CALL inel_init()
    IF(.NOT.lv_init)CALL levin_setup()
    DO i=1,nx
       ord(i)=i
    ENDDO
    DO i=2,nx ! x ascending
       k=ord(i); j=i-1
       DO WHILE(j.GE.1)
          IF(xv(ord(j)).LE.xv(k))EXIT
          ord(j+1)=ord(j); j=j-1
       ENDDO
       ord(j+1)=k
    ENDDO
    DO c0=1,nx,NFROW/2
       c1=MIN(nx,c0+NFROW/2-1)
       ALLOCATE(tk(c1-c0+1)); nt=0
       ! the new b and mu2 are appended to the rows now (counts included), their values filled below
       DO i=c0,c1
          k=ord(i)
          ir=find_frow(xv(k))
          IF(ir.EQ.0)THEN
             next_frow=MOD(next_frow,NFROW)+1
             ir=next_frow
             frows(ir)%x=xv(k); frows(ir)%key=fast_key(); frows(ir)%nb=0; frows(ir)%nmu=0
             IF(ALLOCATED(frows(ir)%b))DEALLOCATE(frows(ir)%b,frows(ir)%mu2,frows(ir)%N)
             ALLOCATE(frows(ir)%b(64),frows(ir)%mu2(16),frows(ir)%N(64,16))
          ENDIF
          t=0
          DO j=1,nt
             IF(tk(j)%ir.EQ.ir)t=j
          ENDDO
          IF(t.EQ.0)THEN
             nt=nt+1; t=nt
             tk(t)%ir=ir; tk(t)%x=frows(ir)%x; tk(t)%ob=frows(ir)%nb; tk(t)%om=frows(ir)%nmu
          ENDIF
          CALL grow_frow(ir,frows(ir)%nb+nbv(k),frows(ir)%nmu+nmv(k))
          DO j=1,nmv(k)
             IF(find_near(frows(ir)%mu2,frows(ir)%nmu,mu2v(j,k)).GT.0)CYCLE
             frows(ir)%nmu=frows(ir)%nmu+1; frows(ir)%mu2(frows(ir)%nmu)=mu2v(j,k)
             tk(t)%nnm=tk(t)%nnm+1
          ENDDO
          DO j=1,nbv(k)
             IF(find_near(frows(ir)%b,frows(ir)%nb,bv(j,k)).GT.0)CYCLE
             frows(ir)%nb=frows(ir)%nb+1; frows(ir)%b(frows(ir)%nb)=bv(j,k)
             tk(t)%nnb=tk(t)%nnb+1
          ENDDO
       ENDDO
       !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(t,ir,ob,om,nnb,nnm)
       DO t=1,nt
          ir=tk(t)%ir; ob=tk(t)%ob; om=tk(t)%om; nnb=tk(t)%nnb; nnm=tk(t)%nnm
          IF(nnm.GT.0.AND.ob.GT.0)THEN ! new mu2 at the old b
             ALLOCATE(tk(t)%Nm(ob,nnm))
             CALL fast_row_compute(tk(t)%x,ob,frows(ir)%b(1:ob),nnm,frows(ir)%mu2(om+1:om+nnm),tk(t)%Nm)
          ENDIF
          IF(nnb.GT.0)THEN ! new b at all mu2
             ALLOCATE(tk(t)%Nb(nnb,om+nnm))
             CALL fast_row_compute(tk(t)%x,nnb,frows(ir)%b(ob+1:ob+nnb),om+nnm,frows(ir)%mu2(1:om+nnm),tk(t)%Nb)
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       DO t=1,nt
          ir=tk(t)%ir; ob=tk(t)%ob; om=tk(t)%om; nnb=tk(t)%nnb; nnm=tk(t)%nnm
          IF(ALLOCATED(tk(t)%Nm))frows(ir)%N(1:ob,om+1:om+nnm)=tk(t)%Nm
          IF(ALLOCATED(tk(t)%Nb))frows(ir)%N(ob+1:ob+nnb,1:om+nnm)=tk(t)%Nb
       ENDDO
       DEALLOCATE(tk)
    ENDDO
  END SUBROUTINE inel_fast_prepare_rows

  SUBROUTINE grow_frow(ir,nb,nmu) ! capacity of row ir at least nb x nmu
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ir,nb,nmu
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::tb,tm
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::tN
    INTEGER::cb,cm
    cb=SIZE(frows(ir)%b); cm=SIZE(frows(ir)%mu2)
    IF(nb.LE.cb.AND.nmu.LE.cm)RETURN
    cb=MAX(cb,2*nb); cm=MAX(cm,2*nmu)
    ALLOCATE(tb(cb),tm(cm),tN(cb,cm))
    tb(1:frows(ir)%nb)=frows(ir)%b(1:frows(ir)%nb)
    tm(1:frows(ir)%nmu)=frows(ir)%mu2(1:frows(ir)%nmu)
    tN(1:frows(ir)%nb,1:frows(ir)%nmu)=frows(ir)%N(1:frows(ir)%nb,1:frows(ir)%nmu)
    CALL MOVE_ALLOC(tb,frows(ir)%b); CALL MOVE_ALLOC(tm,frows(ir)%mu2); CALL MOVE_ALLOC(tN,frows(ir)%N)
  END SUBROUTINE grow_frow

  SUBROUTINE levin_setup()
    ! Chebyshev-Lobatto points t_i (ascending), barycentric differentiation matrix, weights int_{-1}^{1} l_j
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(NLV)::wb
    REAL(KIND(1d0))::l
    INTEGER::i,j,q,k
    DO i=1,NLV
       tlv(i)=-DCOS(pi*(i-1)/(NLV-1))
    ENDDO
    DO i=1,NLV
       wbar(i)=MERGE(0.5d0,1d0,i.EQ.1.OR.i.EQ.NLV)*(-1d0)**(i-1)
    ENDDO
    DO j=1,NLV
       wb(j)=1d0
       DO k=1,NLV
          IF(k.NE.j)wb(j)=wb(j)/(tlv(j)-tlv(k))
       ENDDO
    ENDDO
    Dlv=0d0
    DO i=1,NLV
       DO j=1,NLV
          IF(i.NE.j)Dlv(i,j)=(wb(j)/wb(i))/(tlv(i)-tlv(j))
       ENDDO
       Dlv(i,i)=-SUM(Dlv(i,:))
    ENDDO
    D2lv=MATMUL(Dlv,Dlv)
    wcc=0d0
    DO j=1,NLV
       DO q=1,8
          l=1d0
          DO k=1,NLV
             IF(k.NE.j)l=l*(t8(q)-tlv(k))/(tlv(j)-tlv(k))
          ENDDO
          wcc(j)=wcc(j)+w8(q)*l
       ENDDO
    ENDDO
    CALL cheb_setup(NCL,xcl,wtl,cml)
    CALL cheb_setup(NCS,xcs,wts,cms)
    wvl=MATMUL(wtl,cml); wvs=MATMUL(wts,cms)
    lv_init=.TRUE.
  END SUBROUTINE levin_setup

  SUBROUTINE solve_lin(M,R,n,nr) ! Gaussian elimination with partial pivoting: M X = R, X returned in R(:,1:nr)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n,nr
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(INOUT)::M,R
    REAL(KIND(1d0))::f,t
    INTEGER::i,j,p,l
    DO i=1,n
       p=i
       DO j=i+1,n
          IF(ABS(M(j,i)).GT.ABS(M(p,i)))p=j
       ENDDO
       IF(p.NE.i)THEN
          DO l=i,n
             t=M(i,l); M(i,l)=M(p,l); M(p,l)=t
          ENDDO
          DO l=1,nr
             t=R(i,l); R(i,l)=R(p,l); R(p,l)=t
          ENDDO
       ENDIF
       DO j=i+1,n
          f=M(j,i)/M(i,i)
          IF(f.NE.0d0)THEN
             DO l=i+1,n
                M(j,l)=M(j,l)-f*M(i,l)
             ENDDO
             DO l=1,nr
                R(j,l)=R(j,l)-f*R(i,l)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    DO i=n,1,-1
       DO l=1,nr
          t=R(i,l)
          DO j=i+1,n
             t=t-M(i,j)*R(j,l)
          ENDDO
          R(i,l)=t/M(i,i)
       ENDDO
    ENDDO
  END SUBROUTINE solve_lin

  SUBROUTINE lev_row(x,nb,bv,nmu,mu2v,Nrow)
    ! x f/alpha (rho = delta) at all b (GeV^-1) and mu2 for one x. Panels independent of b: geometric (10 per decade
    ! from 1e-3 A), the k of the ln Q2 nodes of the G^h tables (the integrand is analytic between them), the mu cuts,
    ! the upper limit. Per panel the six integrands at NLV Chebyshev-Lobatto points (once); per b: Clenshaw-Curtis if
    ! b~ (panel width) <= lev_theta, else Levin collocation (p' + A^T p = (g,0), integral = p.(J_n,J_n+1) at the ends).
    ! The Levin system of order n is (M_n/b~ + b~ I) p2 = g with M_n independent of b: M_n = Q H Q^T (Hessenberg, per
    ! panel and n, once), so that per b only (H + b~^2 I) y = b~ Q^T g is solved (O(NLV^2)); the Bessel functions at the
    ! panel ends are shared by neighbouring panels, J_2 and J_3 by recurrence for z >= 2
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nb,nmu
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nb,nmu),INTENT(OUT)::Nrow
    REAL(KIND(1d0)),DIMENSION(nmu)::kmu
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::e
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::G ! (NLV,6,panel)
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::LH,LR,LC ! Levin: H, end rows (Q(1,:),Q(NLV,:),D-rows), Q^T g
    REAL(KIND(1d0)),DIMENSION(6,nmu)::Im
    REAL(KIND(1d0)),DIMENSION(6)::Icum,pk
    REAL(KIND(1d0))::A,hi0,bt,klo,l0,l1,c0,h,z,Jz(0:3),pre,Q2,bt2
    REAL(KIND(1d0))::Ja(0:3),Jc(0:3),kn(NLV)
    REAL(KIND(1d0)),DIMENSION(NLV,NLV)::Sm,Qm
    REAL(KIND(1d0)),DIMENSION(NLV,3)::Rs
    REAL(KIND(1d0))::pa1,pc1,p2a,p2c
    INTEGER::i,j,m,ngeo,ne,np,k,jj,q,nord,kl
    INTEGER,DIMENSION(:),ALLOCATABLE::mlo,mhi
    INTEGER,DIMENSION(nmu)::mord
    ! tree of panels per mu segment (leaves 1..np are the panels): range, children, k range, widest panel, moments of
    ! the Clenshaw-Curtis sums on the NCT Chebyshev-Lobatto points of the node
    INTEGER,PARAMETER::NCT=24
    REAL(KIND(1d0)),PARAMETER::PCT=12d0 ! largest phase b~ (k range) of a node used as a whole (J_n: ~1e-12)
    REAL(KIND(1d0)),DIMENSION(NCT),SAVE::tct,wct
    REAL(KIND(1d0)),DIMENSION(NLV,NCT),SAVE::LLT ! Lagrange basis of the tct at the Lobatto points tlv
    LOGICAL,SAVE::tinit=.FALSE.
    INTEGER,DIMENSION(:),ALLOCATABLE::nc1,nc2,nlo,nhi,roots,lst
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::na,nbd,nwm
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::MM
    REAL(KIND(1d0)),DIMENSION(NCT)::lw
    INTEGER,DIMENSION(64)::stk
    INTEGER::nn,nseg,nl,nl2,ip,ns,nd,ic
    REAL(KIND(1d0))::u,ca,cb
    INTEGER,DIMENSION(0:2),PARAMETER::nrq=(/3,2,1/)
    INTEGER,DIMENSION(3,0:2),PARAMETER::iqn=RESHAPE((/2,3,6,1,5,0,4,0,0/),(/3,3/)) ! the rows of each order
    CALL inel_init()
    IF(.NOT.lv_init)CALL levin_setup()
    IF(.NOT.tinit)THEN
       DO m=1,NCT
          tct(m)=-DCOS(pi*(m-1)/(NCT-1))
       ENDDO
       DO m=1,NCT
          wct(m)=1d0
          DO q=1,NCT
             IF(q.NE.m)wct(m)=wct(m)/(tct(m)-tct(q))
          ENDDO
       ENDDO
       DO jj=1,NLV
          CALL bary_weights(tlv(jj),lw)
          LLT(jj,:)=lw
       ENDDO
       tinit=.TRUE.
    ENDIF
    DO m=1,nmu
       kmu(m)=HUGE(1d0)
       IF(mu2v(m).LT.HUGE(1d0))kmu(m)=DSQRT(MAX(mu2v(m)-x**2*mp2,0d0))
    ENDDO
    Nrow=0d0
    DO i=1,NMX
       A=DSQRT(A2(x,mx_nodes(i)))
       hi0=DSQRT(Q2_fac*((1d0-x)*mx_nodes(i)**2+A**2))
       CALL set_state(x,mx_nodes(i),1d0)
       ! b-independent panels
       klo=1d-3*MIN(A,hi0)
       ngeo=INT(10*LOG10(hi0/klo))+2
       ALLOCATE(e(ngeo+NLQ+nmu+4))
       l0=LOG10(klo); l1=LOG10(hi0)
       DO k=0,ngeo-1
          e(k+1)=10d0**(l0+k*((l1-l0)/(ngeo-1)))
       ENDDO
       ne=ngeo+1; e(ne)=0d0
       DO k=1,NLQ
          Q2=DEXP(lq(k))*(1d0-x)-A**2
          IF(Q2.GT.0d0)THEN
             IF(DSQRT(Q2).LT.hi0)THEN
                ne=ne+1; e(ne)=DSQRT(Q2)
             ENDIF
          ENDIF
       ENDDO
       DO m=1,nmu
          IF(kmu(m).LT.hi0)THEN
             ne=ne+1; e(ne)=kmu(m)
          ENDIF
       ENDDO
       ne=ne+1; e(ne)=hi0
       CALL sort_unique(e,ne,np)
       np=np-1 ! panels e(k),e(k+1), k=1..np
       ALLOCATE(G(NLV,6,np),LH(NLV,NLV,0:2,np),LR(NLV,4,0:2,np),LC(NLV,3,0:2,np),mlo(np),mhi(np))
       ! the mu cuts at each panel start e(k): mord(mlo(k):mhi(k))
       q=0
       DO k=1,np
          mlo(k)=q+1
          DO m=1,nmu
             IF(kmu(m).EQ.e(k))THEN
                q=q+1; mord(q)=m
             ENDIF
          ENDDO
          mhi(k)=q
       ENDDO
       DO k=1,np
          c0=0.5d0*(e(k)+e(k+1)); h=0.5d0*(e(k+1)-e(k))
          DO jj=1,NLV
             CALL rows_at(c0+h*tlv(jj),G(jj,:,k))
          ENDDO
          IF(e(k).LE.0d0)CYCLE ! Clenshaw-Curtis only
          kn=c0+h*tlv
          DO nord=0,2
             DO jj=1,NLV
                DO q=1,NLV
                   Sm(jj,q)=D2lv(jj,q)/h**2-Dlv(jj,q)/h*(nord+1)/kn(q)+nord/kn(jj)*Dlv(jj,q)/h
                ENDDO
                Sm(jj,jj)=Sm(jj,jj)-nord*(nord+1)/kn(jj)**2
             ENDDO
             CALL hessenberg(NLV,Sm,Qm)
             LH(:,:,nord,k)=Sm
             LR(:,1,nord,k)=Qm(1,:); LR(:,2,nord,k)=Qm(NLV,:)
             LR(:,3,nord,k)=MATMUL(Dlv(1,:),Qm)/h; LR(:,4,nord,k)=MATMUL(Dlv(NLV,:),Qm)/h
             DO q=1,nrq(nord)
                LC(:,q,nord,k)=MATMUL(G(:,iqn(q,nord),k),Qm) ! Q^T g
             ENDDO
          ENDDO
       ENDDO
       ! the trees: segments start at k = 1 and at the panels with mu cuts; leaves' moments from the CC sums, a parent's
       ! from its children (the parent basis at the children's points: exact, a polynomial of degree NCT-1)
       ALLOCATE(nc1(2*np),nc2(2*np),nlo(2*np),nhi(2*np),na(2*np),nbd(2*np),nwm(2*np),MM(6,NCT,2*np),roots(np),lst(np))
       DO k=1,np
          nc1(k)=0; nc2(k)=0; nlo(k)=k; nhi(k)=k; na(k)=e(k); nbd(k)=e(k+1); nwm(k)=e(k+1)-e(k)
          h=0.5d0*(e(k+1)-e(k))
          DO q=1,6
             MM(q,:,k)=h*MATMUL(wcc*G(:,q,k),LLT)
          ENDDO
       ENDDO
       nn=np; nseg=0; k=1
       DO WHILE(k.LE.np)
          ! segment k .. ke
          ip=k+1
          DO WHILE(ip.LE.np)
             IF(mhi(ip).GE.mlo(ip))EXIT
             ip=ip+1
          ENDDO
          nl=0
          DO q=k,ip-1
             nl=nl+1; lst(nl)=q
          ENDDO
          DO WHILE(nl.GT.1)
             nl2=0
             DO q=1,nl-1,2
                nn=nn+1
                nc1(nn)=lst(q); nc2(nn)=lst(q+1)
                nlo(nn)=nlo(lst(q)); nhi(nn)=nhi(lst(q+1)); na(nn)=na(lst(q)); nbd(nn)=nbd(lst(q+1))
                nwm(nn)=MAX(nwm(lst(q)),nwm(lst(q+1)))
                MM(:,:,nn)=0d0
                DO ic=1,2
                   nd=MERGE(lst(q),lst(q+1),ic.EQ.1)
                   DO m=1,NCT
                      ! child point m in the parent's coordinate
                      u=(2d0*(0.5d0*(na(nd)+nbd(nd))+0.5d0*(nbd(nd)-na(nd))*tct(m))-na(nn)-nbd(nn))/(nbd(nn)-na(nn))
                      CALL bary_weights(u,lw)
                      DO jj=1,NCT
                         MM(:,jj,nn)=MM(:,jj,nn)+MM(:,m,nd)*lw(jj)
                      ENDDO
                   ENDDO
                ENDDO
                nl2=nl2+1; lst(nl2)=nn
             ENDDO
             IF(MOD(nl,2).EQ.1)THEN
                nl2=nl2+1; lst(nl2)=lst(nl)
             ENDIF
             nl=nl2
          ENDDO
          nseg=nseg+1; roots(nseg)=lst(1)
          k=ip
       ENDDO
       DO j=1,nb
          bt=bv(j)/(1d0-x)
          IF(A*bt.GT.Ab_max)CYCLE ! (the b of a row need not be ascending)
          bt2=bt*bt
          Icum=0d0; Im=0d0
          kl=-1 ! the panel end where Jc was last computed
          DO ns=1,nseg
           k=nlo(roots(ns))
           DO m=mlo(k),mhi(k) ! the mu cuts at the start of the segment
              Im(:,mord(m))=Icum
           ENDDO
           ip=1; stk(1)=roots(ns)
           DO WHILE(ip.GT.0)
             nd=stk(ip); ip=ip-1
             pk=0d0
             IF(nc1(nd).GT.0)THEN
                IF(bt*(nbd(nd)-na(nd)).LE.PCT.AND.bt*nwm(nd).LE.lev_theta)THEN
                   ! the whole node: Clenshaw-Curtis sums of its panels with J_n interpolated on its Chebyshev points
                   ca=0.5d0*(na(nd)+nbd(nd)); cb=0.5d0*(nbd(nd)-na(nd))
                   DO m=1,NCT
                      CALL bessel_j0123(bt*(ca+cb*tct(m)),Jz,2)
                      DO q=1,6
                         pk(q)=pk(q)+MM(q,m,nd)*Jz(row_nj(q))
                      ENDDO
                   ENDDO
                   Icum=Icum+pk
                ELSE
                   ip=ip+1; stk(ip)=nc2(nd); ip=ip+1; stk(ip)=nc1(nd)
                ENDIF
                CYCLE
             ENDIF
             k=nd
             c0=0.5d0*(e(k)+e(k+1)); h=0.5d0*(e(k+1)-e(k))
             IF(bt*2d0*h.LE.lev_theta.OR.e(k).LE.0d0)THEN ! Clenshaw-Curtis on the Lobatto points
                DO jj=1,NLV
                   z=bt*(c0+h*tlv(jj))
                   CALL bessel_j0123(z,Jz,2)
                   DO q=1,6
                      pk(q)=pk(q)+h*wcc(jj)*G(jj,q,k)*Jz(row_nj(q))
                   ENDDO
                ENDDO
             ELSE ! Levin
                IF(kl.EQ.k)THEN
                   Ja=Jc
                ELSE
                   CALL bessel_j0123(bt*e(k),Ja,3)
                ENDIF
                CALL bessel_j0123(bt*e(k+1),Jc,3)
                kl=k+1
                DO nord=0,2
                   Sm=LH(:,:,nord,k)
                   DO jj=1,NLV
                      Sm(jj,jj)=Sm(jj,jj)+bt2
                   ENDDO
                   Rs(:,1:nrq(nord))=bt*LC(:,1:nrq(nord),nord,k)
                   CALL solve_hess(NLV,Sm,Rs,nrq(nord))
                   DO q=1,nrq(nord)
                      ! p2 and p1 = (D/h - (n+1)/k) p2 / b~ at k = a (row 1) and k = c (row NLV)
                      p2a=DOT_PRODUCT(LR(:,1,nord,k),Rs(:,q)); p2c=DOT_PRODUCT(LR(:,2,nord,k),Rs(:,q))
                      pa1=(DOT_PRODUCT(LR(:,3,nord,k),Rs(:,q))-(nord+1)/e(k)*p2a)/bt
                      pc1=(DOT_PRODUCT(LR(:,4,nord,k),Rs(:,q))-(nord+1)/e(k+1)*p2c)/bt
                      pk(iqn(q,nord))=pc1*Jc(nord)+p2c*Jc(nord+1)-pa1*Ja(nord)-p2a*Ja(nord+1)
                   ENDDO
                ENDDO
             ENDIF
             Icum=Icum+pk
           ENDDO
          ENDDO
          DO m=1,nmu
             IF(kmu(m).GE.hi0)Im(:,m)=Icum
             IF(mu2v(m).LT.HUGE(1d0).AND.A**2.GT.A2_mu*mu2v(m))CYCLE
             pre=alpha_norm(mu2v(m))/pi**2/(1d0-x)**2
             Nrow(j,m)=Nrow(j,m)+mx_weights(i)*pre*combination_v(Im(:,m),x,mx_nodes(i))
          ENDDO
       ENDDO
       DEALLOCATE(e,G,LH,LR,LC,mlo,mhi,nc1,nc2,nlo,nhi,na,nbd,nwm,MM,roots,lst)
    ENDDO
  CONTAINS
    SUBROUTINE bary_weights(uu,ww) ! Lagrange basis of the tct at uu
      REAL(KIND(1d0)),INTENT(IN)::uu
      REAL(KIND(1d0)),DIMENSION(NCT),INTENT(OUT)::ww
      INTEGER::mm
      mm=MINLOC(ABS(uu-tct),1)
      IF(ABS(uu-tct(mm)).LT.1d-14)THEN
         ww=0d0; ww(mm)=1d0
      ELSE
         ww=wct/(uu-tct); ww=ww/SUM(ww)
      ENDIF
    END SUBROUTINE bary_weights
  END SUBROUTINE lev_row

  SUBROUTINE bessel_j0123(z,J,nmax)
    ! J_0 .. J_nmax (nmax <= 3) at z >= 0: J_2 = 2 J_1/z - J_0 for z >= 1e-3 (absolute error ~1e-16, as J_0, J_1),
    ! its series below; J_3 = 4 J_2/z - J_1 for z >= 2, BESSEL_JN below
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::z
    INTEGER,INTENT(IN)::nmax
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::J
    J=0d0
    J(0)=BESSEL_J0(z); J(1)=BESSEL_J1(z)
    IF(z.GE.1d-3)THEN
       J(2)=2d0*J(1)/z-J(0)
    ELSE
       J(2)=z*z/8d0*(1d0-z*z/12d0)
    ENDIF
    IF(nmax.GE.3)THEN
       IF(z.GE.2d0)THEN
          J(3)=4d0*J(2)/z-J(1)
       ELSE
          J(3)=BESSEL_JN(3,z)
       ENDIF
    ENDIF
  END SUBROUTINE bessel_j0123

  SUBROUTINE hessenberg(n,A,Q) ! A -> H = Q^T A Q upper Hessenberg (Householder), Q orthogonal
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n,n),INTENT(INOUT)::A
    REAL(KIND(1d0)),DIMENSION(n,n),INTENT(OUT)::Q
    REAL(KIND(1d0)),DIMENSION(n)::v,w
    REAL(KIND(1d0))::al,vn
    INTEGER::k,i
    Q=0d0
    DO i=1,n
       Q(i,i)=1d0
    ENDDO
    DO k=1,n-2
       al=DSQRT(SUM(A(k+1:n,k)**2))
       IF(al.EQ.0d0)CYCLE
       IF(A(k+1,k).GT.0d0)al=-al
       v=0d0
       v(k+1:n)=A(k+1:n,k)
       v(k+1)=v(k+1)-al
       vn=SUM(v(k+1:n)**2)
       IF(vn.EQ.0d0)CYCLE
       ! A = P A P, Q = Q P, P = I - 2 v v^T/vn
       w=MATMUL(v(k+1:n),A(k+1:n,:))
       DO i=k+1,n
          A(i,:)=A(i,:)-(2d0/vn)*v(i)*w
       ENDDO
       w=MATMUL(A(:,k+1:n),v(k+1:n))
       DO i=k+1,n
          A(:,i)=A(:,i)-(2d0/vn)*w*v(i)
       ENDDO
       w=MATMUL(Q(:,k+1:n),v(k+1:n))
       DO i=k+1,n
          Q(:,i)=Q(:,i)-(2d0/vn)*w*v(i)
       ENDDO
       A(k+2:n,k)=0d0
    ENDDO
  END SUBROUTINE hessenberg

  SUBROUTINE solve_hess(n,H,R,nr) ! H X = R for upper Hessenberg H (pivoting between adjacent rows), X in R(:,1:nr)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n,nr
    REAL(KIND(1d0)),DIMENSION(n,n),INTENT(INOUT)::H
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(INOUT)::R
    REAL(KIND(1d0))::f,t
    INTEGER::i,j,l
    DO i=1,n-1
       IF(ABS(H(i+1,i)).GT.ABS(H(i,i)))THEN
          DO l=i,n
             t=H(i,l); H(i,l)=H(i+1,l); H(i+1,l)=t
          ENDDO
          DO l=1,nr
             t=R(i,l); R(i,l)=R(i+1,l); R(i+1,l)=t
          ENDDO
       ENDIF
       f=H(i+1,i)/H(i,i)
       IF(f.NE.0d0)THEN
          DO l=i+1,n
             H(i+1,l)=H(i+1,l)-f*H(i,l)
          ENDDO
          DO l=1,nr
             R(i+1,l)=R(i+1,l)-f*R(i,l)
          ENDDO
       ENDIF
    ENDDO
    DO i=n,1,-1
       DO l=1,nr
          t=R(i,l)
          DO j=i+1,n
             t=t-H(i,j)*R(j,l)
          ENDDO
          R(i,l)=t/H(i,i)
       ENDDO
    ENDDO
  END SUBROUTINE solve_hess

  ! ================================================================ fast rows with rho (rho_model /= 0)
  SUBROUTINE fast_row_compute(x,nb,bv,nmu,mu2v,Nrow)
    ! x f/alpha at all b (GeV^-1) and mu2 of one x: lev_row (rho = delta) or rho_row
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nb,nmu
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nb,nmu),INTENT(OUT)::Nrow
    INTEGER,DIMENSION(7)::hi
    REAL(KIND(1d0)),DIMENSION(3)::hr
    CHARACTER(LEN=600)::fname
    LOGICAL::cached
    cached=LEN_TRIM(inel_cache_dir).GT.0
    IF(cached)THEN
       hi=(/CACHE_VERSION,rho_model,low_Q2_model,MERGE(1,0,alpha_running4photon),sf_pdf_id,nb,nmu/)
       hr=(/sf_xc,MIN(sf_Q2c,1d300),x/)
       fname=TRIM(inel_cache_dir)//'/gammaUPC_inel_'//cache_hash(hi,hr,bv,mu2v)//'.bin'
       IF(cache_read(fname,hi,hr,nb,bv,nmu,mu2v,Nrow))THEN
          !$OMP ATOMIC
          inel_cache_hits=inel_cache_hits+1
          RETURN
       ENDIF
    ENDIF
    IF(rho_model.EQ.0)THEN
       CALL lev_row(x,nb,bv,nmu,mu2v,Nrow)
    ELSE
       CALL rho_row(x,nb,bv,nmu,mu2v,Nrow)
    ENDIF
    IF(cached)CALL cache_write(fname,hi,hr,nb,bv,nmu,mu2v,Nrow)
  END SUBROUTINE fast_row_compute

  FUNCTION cache_hash(hi,hr,bv,mu2v) ! two 32-bit FNV-1a hashes of the bytes of the key, as 16 hex digits
    IMPLICIT NONE
    CHARACTER(LEN=16)::cache_hash
    INTEGER,DIMENSION(:),INTENT(IN)::hi
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::hr,bv,mu2v
    INTEGER(KIND=8)::h1,h2
    h1=2166136261_8; h2=84696351_8
    CALL fnv(TRANSFER(hi,(/0_1/))); CALL fnv(TRANSFER(hr,(/0_1/)))
    CALL fnv(TRANSFER(bv,(/0_1/))); CALL fnv(TRANSFER(mu2v,(/0_1/)))
    WRITE(cache_hash,'(2Z8.8)')h1,h2
  CONTAINS
    SUBROUTINE fnv(by)
      INTEGER(KIND=1),DIMENSION(:),INTENT(IN)::by
      INTEGER::i
      DO i=1,SIZE(by)
         h1=IAND(IEOR(h1,IAND(INT(by(i),8),255_8))*16777619_8,4294967295_8)
         h2=IAND(IEOR(h2,IAND(INT(by(i),8),255_8))*16777619_8+1_8,4294967295_8)
      ENDDO
    END SUBROUTINE fnv
  END FUNCTION cache_hash

  LOGICAL FUNCTION cache_read(fname,hi,hr,nb,bv,nmu,mu2v,Nrow)
    ! the block of fname if its key (hi, hr, bv, mu2v) agrees exactly
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::fname
    INTEGER,INTENT(IN)::nb,nmu
    INTEGER,DIMENSION(7),INTENT(IN)::hi
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::hr
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nb,nmu),INTENT(OUT)::Nrow
    CHARACTER(LEN=8)::magic
    INTEGER,DIMENSION(7)::fi
    REAL(KIND(1d0)),DIMENSION(3)::fr
    REAL(KIND(1d0)),DIMENSION(nb)::fb
    REAL(KIND(1d0)),DIMENSION(nmu)::fm
    INTEGER::u,ios,endmark
    cache_read=.FALSE.
    OPEN(NEWUNIT=u,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=ios)
    IF(ios.NE.0)RETURN
    READ(u,IOSTAT=ios)magic,fi,fr
    IF(ios.EQ.0.AND.magic.EQ.'GUPCINEL'.AND.ALL(fi.EQ.hi).AND.ALL(fr.EQ.hr))THEN
       READ(u,IOSTAT=ios)fb,fm
       IF(ios.EQ.0.AND.ALL(fb.EQ.bv).AND.ALL(fm.EQ.mu2v))THEN
          READ(u,IOSTAT=ios)Nrow,endmark
          cache_read=ios.EQ.0.AND.endmark.EQ.20261121
       ENDIF
    ENDIF
    CLOSE(u)
  END FUNCTION cache_read

  SUBROUTINE cache_write(fname,hi,hr,nb,bv,nmu,mu2v,Nrow)
    ! fname through a temporary file (process id) and a rename; a warning (once) if the folder is not writable
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::fname
    INTEGER,INTENT(IN)::nb,nmu
    INTEGER,DIMENSION(7),INTENT(IN)::hi
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::hr
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nb,nmu),INTENT(IN)::Nrow
    CHARACTER(LEN=640)::tmp
    CHARACTER(LEN=12)::pid
    INTEGER::u,ios
    WRITE(pid,'(I12)')GETPID()
    tmp=TRIM(fname)//'.'//TRIM(ADJUSTL(pid))//'.tmp'
    OPEN(NEWUNIT=u,FILE=TRIM(tmp),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios)
    IF(ios.EQ.0)THEN
       WRITE(u,IOSTAT=ios)'GUPCINEL',hi,hr,bv,mu2v,Nrow,20261121
       CLOSE(u)
    ENDIF
    IF(ios.EQ.0)CALL RENAME(TRIM(tmp),TRIM(fname),ios)
    IF(ios.EQ.0)THEN
       !$OMP ATOMIC
       inel_cache_writes=inel_cache_writes+1
    ELSEIF(.NOT.cache_warned)THEN
       !$OMP CRITICAL(inel_cache_warning)
       IF(.NOT.cache_warned)WRITE(*,*)"WARNING: InelasticFlux: cannot write the row cache in ",TRIM(inel_cache_dir),&
            " (inel_cache_dir; '' disables it)"
       cache_warned=.TRUE.
       !$OMP END CRITICAL(inel_cache_warning)
    ENDIF
  END SUBROUTINE cache_write

  SUBROUTINE node_panels_mu(x,mx,hi,nmu,kcv,Bexp,e,np,G)
    ! b-independent panels of one mx node up to hi (geometric 10/decade from 1e-3 A, the k of the ln Q2 table nodes,
    ! the cuts kcv) and the six integrands at the NLV Chebyshev-Lobatto points of each panel, times exp(-Bexp k^2)
    ! if Bexp > 0
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,mx,hi,Bexp
    INTEGER,INTENT(IN)::nmu
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::kcv
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,INTENT(OUT)::e
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)::G
    INTEGER,INTENT(OUT)::np
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::t
    REAL(KIND(1d0))::A,klo,l0,l1,Q2,c0,h,kk
    INTEGER::k,ngeo,ne,jj
    CALL set_state(x,mx,1d0)
    A=DSQRT(A2(x,mx))
    klo=1d-3*MIN(A,hi)
    ngeo=INT(10*LOG10(hi/klo))+2
    ALLOCATE(t(ngeo+NLQ+nmu+3))
    l0=LOG10(klo); l1=LOG10(hi)
    DO k=0,ngeo-1
       t(k+1)=10d0**(l0+k*((l1-l0)/(ngeo-1)))
    ENDDO
    ne=ngeo+1; t(ne)=0d0
    DO k=1,NLQ
       Q2=DEXP(lq(k))*(1d0-x)-A**2
       IF(Q2.GT.0d0)THEN
          IF(DSQRT(Q2).LT.hi)THEN
             ne=ne+1; t(ne)=DSQRT(Q2)
          ENDIF
       ENDIF
    ENDDO
    DO k=1,nmu
       IF(kcv(k).LT.hi)THEN
          ne=ne+1; t(ne)=kcv(k)
       ENDIF
    ENDDO
    ne=ne+1; t(ne)=hi
    CALL sort_unique(t,ne,np)
    ALLOCATE(e(np)); e=t(1:np)
    np=np-1
    ALLOCATE(G(NLV,6,np))
    DO k=1,np
       c0=0.5d0*(e(k)+e(k+1)); h=0.5d0*(e(k+1)-e(k))
       DO jj=1,NLV
          kk=c0+h*tlv(jj)
          CALL rows_at(kk,G(jj,:,k))
          IF(Bexp.GT.0d0)G(jj,:,k)=G(jj,:,k)*DEXP(-Bexp*kk**2)
       ENDDO
    ENDDO
  END SUBROUTINE node_panels_mu

  SUBROUTINE cheb_eval(k,e,np,G,f,hint)
    ! the six integrands at k from their values G(:,q,panel) at the NLV Chebyshev-Lobatto points of each panel
    ! (barycentric interpolation)
    ! hint (optional): the panel of the previous call, tried first (k increasing along the inner integrals)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::k
    INTEGER,INTENT(IN)::np
    REAL(KIND(1d0)),DIMENSION(np+1),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(NLV,6,np),INTENT(IN)::G
    REAL(KIND(1d0)),DIMENSION(6),INTENT(OUT)::f
    INTEGER,INTENT(INOUT),OPTIONAL::hint
    REAL(KIND(1d0))::t,w,sw
    INTEGER::lo,hi,mid,j
    lo=1; hi=np+1
    IF(PRESENT(hint))THEN
       IF(hint.GE.1.AND.hint.LE.np)THEN
          IF(k.GE.e(hint))THEN
             lo=hint
             IF(hint.LT.np)THEN
                IF(k.LT.e(hint+1))THEN
                   hi=hint+1
                ELSEIF(hint+1.LT.np)THEN
                   IF(k.LT.e(hint+2))THEN
                      lo=hint+1; hi=hint+2
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    DO WHILE(hi-lo.GT.1)
       mid=(lo+hi)/2
       IF(k.GE.e(mid))THEN
          lo=mid
       ELSE
          hi=mid
       ENDIF
    ENDDO
    IF(PRESENT(hint))hint=lo
    t=(2d0*k-e(lo)-e(lo+1))/(e(lo+1)-e(lo))
    f=0d0; sw=0d0
    DO j=1,NLV
       IF(t.EQ.tlv(j))THEN
          f=G(j,:,lo); RETURN
       ENDIF
       w=wbar(j)/(t-tlv(j))
       f(1)=f(1)+w*G(j,1,lo); f(2)=f(2)+w*G(j,2,lo); f(3)=f(3)+w*G(j,3,lo)
       f(4)=f(4)+w*G(j,4,lo); f(5)=f(5)+w*G(j,5,lo); f(6)=f(6)+w*G(j,6,lo)
       sw=sw+w
    ENDDO
    f=f/sw
  END SUBROUTINE cheb_eval

  FUNCTION bary(t,y) ! barycentric interpolation on the NLV Chebyshev-Lobatto points
    REAL(KIND(1d0))::bary
    REAL(KIND(1d0)),INTENT(IN)::t
    REAL(KIND(1d0)),DIMENSION(NLV),INTENT(IN)::y
    REAL(KIND(1d0))::w,sw
    INTEGER::j
    bary=0d0; sw=0d0
    DO j=1,NLV
       IF(t.EQ.tlv(j))THEN
          bary=y(j); RETURN
       ENDIF
       w=wbar(j)/(t-tlv(j))
       bary=bary+w*y(j); sw=sw+w
    ENDDO
    bary=bary/sw
  END FUNCTION bary

  SUBROUTINE lev_transform(bt,e,np,Gp,nr,nj,I)
    ! I(r) = int dk g_r(k) J_nj(r)(bt k) over the panels e (g_r at the NLV Chebyshev-Lobatto points in Gp(:,r,panel)):
    ! Clenshaw-Curtis on panels with bt (width) <= lev_theta (or starting at 0), Levin collocation elsewhere (as lev_row)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::bt
    INTEGER,INTENT(IN)::np,nr
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(NLV,nr,np),INTENT(IN)::Gp
    INTEGER,DIMENSION(nr),INTENT(IN)::nj
    REAL(KIND(1d0)),DIMENSION(nr),INTENT(OUT)::I
    REAL(KIND(1d0))::c0,h,z,Jz(0:3),Ja(0:3),Jc(0:3),kn(NLV),Sm(NLV,NLV),pa1,pc1
    REAL(KIND(1d0)),DIMENSION(NLV,nr)::Rs
    INTEGER,DIMENSION(nr)::iq
    INTEGER::k,jj,q,nord,nn
    I=0d0
    DO k=1,np
       c0=0.5d0*(e(k)+e(k+1)); h=0.5d0*(e(k+1)-e(k))
       IF(bt*2d0*h.LE.lev_theta.OR.e(k).LE.0d0)THEN
          DO jj=1,NLV
             z=bt*(c0+h*tlv(jj))
             Jz(0)=BESSEL_J0(z); Jz(1)=BESSEL_J1(z); Jz(2)=BESSEL_JN(2,z)
             DO q=1,nr
                I(q)=I(q)+h*wcc(jj)*Gp(jj,q,k)*Jz(nj(q))
             ENDDO
          ENDDO
       ELSE
          kn=c0+h*tlv
          Ja(0)=BESSEL_J0(bt*e(k)); Ja(1)=BESSEL_J1(bt*e(k)); Ja(2)=BESSEL_JN(2,bt*e(k)); Ja(3)=BESSEL_JN(3,bt*e(k))
          Jc(0)=BESSEL_J0(bt*e(k+1)); Jc(1)=BESSEL_J1(bt*e(k+1)); Jc(2)=BESSEL_JN(2,bt*e(k+1)); Jc(3)=BESSEL_JN(3,bt*e(k+1))
          DO nord=0,2
             nn=0
             DO q=1,nr
                IF(nj(q).EQ.nord)THEN
                   nn=nn+1; iq(nn)=q
                ENDIF
             ENDDO
             IF(nn.EQ.0)CYCLE
             DO jj=1,NLV
                DO q=1,NLV
                   Sm(jj,q)=(D2lv(jj,q)/h**2-Dlv(jj,q)/h*(nord+1)/kn(q)+nord/kn(jj)*Dlv(jj,q)/h)/bt
                ENDDO
                Sm(jj,jj)=Sm(jj,jj)-nord*(nord+1)/(kn(jj)**2*bt)+bt
                DO q=1,nn
                   Rs(jj,q)=Gp(jj,iq(q),k)
                ENDDO
             ENDDO
             CALL solve_lin(Sm,Rs,NLV,nn)
             DO q=1,nn
                pa1=(DOT_PRODUCT(Dlv(1,:),Rs(:,q))/h-(nord+1)/kn(1)*Rs(1,q))/bt
                pc1=(DOT_PRODUCT(Dlv(NLV,:),Rs(:,q))/h-(nord+1)/kn(NLV)*Rs(NLV,q))/bt
                I(iq(q))=I(iq(q))+pc1*Jc(nord)+Rs(NLV,q)*Jc(nord+1)-pa1*Ja(nord)-Rs(1,q)*Ja(nord+1)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE lev_transform

  SUBROUTINE amps_mu(bt,nmu,kcv,e,np,G,I)
    ! the six I_p(bt) for every cut in kcv (partial sums at the cut edges; cuts beyond the last edge: the full sum)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::bt
    INTEGER,INTENT(IN)::nmu,np
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::kcv
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(:,:,:),INTENT(IN)::G
    REAL(KIND(1d0)),DIMENSION(6,nmu),INTENT(OUT)::I
    REAL(KIND(1d0)),DIMENSION(6)::Ic,pk
    INTEGER::k,m
    Ic=0d0; I=0d0
    DO k=1,np
       DO m=1,nmu
          IF(kcv(m).EQ.e(k))I(:,m)=Ic
       ENDDO
       CALL lev_transform(bt,e(k:k+1),1,G(:,:,k:k),6,row_nj,pk)
       Ic=Ic+pk
    ENDDO
    DO m=1,nmu
       IF(kcv(m).GE.e(np+1))I(:,m)=Ic
    ENDDO
  END SUBROUTINE amps_mu

  SUBROUTINE rho_phat_mu(nmu,kcv,hi0,nD,Dv,e,np,G,Ph,COARSE)
    ! P^_pq(Delta) = R(p,q) + R(q,p), R(p,q) = int dk1 f_p(k1) int dpsi 2 f_q(k2)/k2 cos(n phi) over the half plane
    ! k2 >= k1 (from int J_n J_n J_0 t dt = cos(n phi)/(2 pi area)), for all cuts kcv at once: for k1 < kc - Delta the
    ! inner integral is independent of kc, so a prefix sum over k1 serves all cuts and only k1 in [kc - Delta, kc) is
    ! integrated again per cut (inner_phi_cc: all cuts from one set of integrand values); Delta = 0: Parseval,
    ! 2 pi int dk f_p f_q/k. Outer k1 panels: Gauss-Legendre 8 between
    ! the node edges, 0.5 Delta, the cuts and the cuts - Delta, Delta - the cuts. At 0.5 Delta (where the half plane
    ! starts to cut the psi range) and at Delta - kc (where the cut kc starts to cut it) the outer integrand rises as
    ! the square root of the distance: on the panel that starts there k1 = a + (b-a) s^2, Gauss-Legendre in s
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nmu,nD,np
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::kcv
    REAL(KIND(1d0)),INTENT(IN)::hi0
    REAL(KIND(1d0)),DIMENSION(nD),INTENT(IN)::Dv
    REAL(KIND(1d0)),DIMENSION(np+1),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(NLV,6,np),INTENT(IN)::G
    REAL(KIND(1d0)),DIMENSION(6,6,nD,nmu),INTENT(OUT)::Ph
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::COARSE ! (M11) > 0: outer panels geometric 6 per decade from COARSE to hi0
                                                ! (and the special points) instead of the node panels
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::o,sq
    REAL(KIND(1d0)),DIMENSION(6)::f1,Phi
    REAL(KIND(1d0)),DIMENSION(6,6)::Rs
    REAL(KIND(1d0)),DIMENSION(6,6,nmu)::Rb,Rpart
    REAL(KIND(1d0))::D,k1,w1,kc(nmu)
    INTEGER::id,no,nn,k,j,m,p,q
    LOGICAL,DIMENSION(nmu)::based
    REAL(KIND(1d0)),DIMENSION(6,nmu)::Phic
    REAL(KIND(1d0)),DIMENSION(nmu)::kcc
    INTEGER,DIMENSION(nmu)::ic
    INTEGER::nc,c,npr,ipr,nsq,isq
    LOGICAL::sqe
    REAL(KIND(1d0))::sj
    INTEGER,DIMENSION(14)::pp,qq
    ! (M13, COARSE only) k1 < ks = 0.1 D with no cut inside the inner range and hi0 > D + ks: the inner integral
    ! interpolated from NS Chebyshev points of [0, ks] (k2 within D (1 +- 0.1): smooth in k1)
    INTEGER,PARAMETER::NS=6
    REAL(KIND(1d0)),DIMENSION(NS)::tsn,lw
    REAL(KIND(1d0)),DIMENSION(6,NS)::PhiS
    REAL(KIND(1d0)),DIMENSION(6,1)::Pd
    REAL(KIND(1d0))::ks
    INTEGER::js,jl
    kc=MIN(kcv,hi0)
    npr=0 ! the pairs with the same Bessel order
    DO p=1,6
       DO q=1,6
          IF(row_nj(p).EQ.row_nj(q))THEN
             npr=npr+1; pp(npr)=p; qq(npr)=q
          ENDIF
       ENDDO
    ENDDO
    DO id=1,nD
       D=Dv(id)
       IF(D.LE.0d0)THEN
          Rs=0d0
          DO m=1,nmu
             Rb(:,:,m)=-1d0
          ENDDO
          DO k=1,np
             DO m=1,nmu
                IF(kc(m).LE.e(k).AND.Rb(1,1,m).LT.0d0)Rb(:,:,m)=Rs
             ENDDO
             DO j=1,8
                k1=0.5d0*(e(k)+e(k+1))+0.5d0*(e(k+1)-e(k))*t8(j)
                w1=0.5d0*(e(k+1)-e(k))*w8(j)
                CALL cheb_eval(k1,e,np,G,f1)
                DO p=1,6
                   DO q=1,6
                      IF(row_nj(p).EQ.row_nj(q))Rs(p,q)=Rs(p,q)+twopi*w1*f1(p)*f1(q)/k1
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          DO m=1,nmu
             IF(Rb(1,1,m).LT.0d0)Rb(:,:,m)=Rs
             Ph(:,:,id,m)=Rb(:,:,m)
          ENDDO
          CYCLE
       ENDIF
       nsq=0
       no=0
       IF(PRESENT(COARSE))THEN
          nn=MAX(1,CEILING(6d0*LOG10(MAX(hi0/COARSE,10d0))))
          ALLOCATE(o(nn+3*nmu+5),sq(nmu+1))
          no=1; o(1)=0d0
          DO k=0,nn
             no=no+1; o(no)=COARSE*(hi0/COARSE)**(DBLE(k)/nn)
          ENDDO
          no=no-1 ! hi0 itself is added below
       ELSE
          ALLOCATE(o(np+3*nmu+3),sq(nmu+1))
          DO k=1,np+1
             IF(e(k).LT.hi0)THEN
                no=no+1; o(no)=e(k)
             ENDIF
          ENDDO
       ENDIF
       no=no+1; o(no)=hi0
       IF(0.5d0*D.LT.hi0)THEN
          no=no+1; o(no)=0.5d0*D
          nsq=nsq+1; sq(nsq)=0.5d0*D
       ENDIF
       DO m=1,nmu
          no=no+1; o(no)=kc(m)
          IF(kc(m)-D.GT.0d0)THEN
             no=no+1; o(no)=kc(m)-D
          ENDIF
          IF(D-kc(m).GT.0d0.AND.D-kc(m).LT.kc(m))THEN
             no=no+1; o(no)=D-kc(m)
             nsq=nsq+1; sq(nsq)=D-kc(m)
          ENDIF
       ENDDO
       ks=0d0
       IF(PRESENT(COARSE))THEN
          IF(hi0.GT.1.1d0*D.AND.0.1d0*D.GT.COARSE)THEN
             ks=0.1d0*D
             no=no+1; o(no)=ks
             DO js=1,NS
                tsn(js)=0.5d0*ks*(1d0-DCOS((2*js-1)*pi/(2*NS)))
                CALL inner_phi_cc(tsn(js),D,hi0,0,kc(1:0),e,np,G,PhiS(:,js),Pd(:,1:0))
             ENDDO
          ENDIF
       ENDIF
       CALL sort_unique(o,no,nn)
       Rs=0d0; Rpart=0d0; based=.FALSE.; Rb=0d0
       DO k=1,nn-1
          DO m=1,nmu ! the prefix up to kc-Delta serves cut m
             IF(.NOT.based(m).AND.o(k).GE.MAX(kc(m)-D,0d0))THEN
                Rb(:,:,m)=Rs; based(m)=.TRUE.
             ENDIF
          ENDDO
          sqe=.FALSE.
          DO isq=1,nsq
             IF(ABS(o(k)-sq(isq)).LE.1d-13*o(k))sqe=.TRUE.
          ENDDO
          DO j=1,8
             IF(sqe)THEN ! square-root onset at o(k)
                sj=0.5d0*(1d0+t8(j))
                k1=o(k)+(o(k+1)-o(k))*sj**2
                w1=(o(k+1)-o(k))*sj*w8(j)
             ELSE
                k1=0.5d0*(o(k)+o(k+1))+0.5d0*(o(k+1)-o(k))*t8(j)
                w1=0.5d0*(o(k+1)-o(k))*w8(j)
             ENDIF
             CALL cheb_eval(k1,e,np,G,f1)
             nc=0
             DO m=1,nmu ! k1 in [kc-Delta, kc): inner cut at kc
                IF(kc(m).GE.hi0)CYCLE
                IF(k1.LT.kc(m)-D.OR.k1.GE.kc(m))CYCLE
                nc=nc+1; ic(nc)=m; kcc(nc)=kc(m)
             ENDDO
             IF(k1.LT.ks.AND.nc.EQ.0)THEN
                DO js=1,NS
                   lw(js)=1d0
                   DO jl=1,NS
                      IF(jl.NE.js)lw(js)=lw(js)*(k1-tsn(jl))/(tsn(js)-tsn(jl))
                   ENDDO
                ENDDO
                Phi=MATMUL(PhiS,lw)
             ELSE
                CALL inner_phi_cc(k1,D,hi0,nc,kcc(1:nc),e,np,G,Phi,Phic(:,1:nc))
             ENDIF
             DO ipr=1,npr
                Rs(pp(ipr),qq(ipr))=Rs(pp(ipr),qq(ipr))+w1*f1(pp(ipr))*Phi(qq(ipr))
             ENDDO
             DO c=1,nc
                m=ic(c)
                DO ipr=1,npr
                   Rpart(pp(ipr),qq(ipr),m)=Rpart(pp(ipr),qq(ipr),m)+w1*f1(pp(ipr))*Phic(qq(ipr),c)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO m=1,nmu
          IF(kc(m).GE.hi0)THEN
             Ph(:,:,id,m)=Rs+TRANSPOSE(Rs)
          ELSE
             IF(.NOT.based(m))Rb(:,:,m)=Rs
             Ph(:,:,id,m)=Rb(:,:,m)+Rpart(:,:,m)+TRANSPOSE(Rb(:,:,m)+Rpart(:,:,m))
          ENDIF
       ENDDO
       DEALLOCATE(o,sq)
    ENDDO
  END SUBROUTINE rho_phat_mu

  SUBROUTINE rho_row(x,nb,bv,nmu,mu2v,Nrow)
    ! x f/alpha with rho (rho_model) at all b and mu2 of one x, per mx node with s = W_step(mx): resonances (s = 0)
    ! products of the amplitudes (amps_mu); continuum (s = 1) the rho pairs (1/2 pi) int Delta dDelta J0(b~ Delta)
    ! exp(-B Delta^2) P^(Delta) (rho_phat_mu; Levin/Clenshaw-Curtis on Delta panels rho_dens per decade from rho_dlo
    ! kappa to Dmax = sqrt(16 ln 10/B), kappa = A; P^ exact on every panel: (M15) the former interpolant of P^ in
    ! Delta^2 on [0, Dmax^2] (narrow nodes) or [0, (rho_d2frac kappa)^2] was off by up to 5e-5 of the peak with a
    ! finite mu cut, now <= 6e-7 against refined rows (drivers/rep_inel; a row costs 3.3 times more)); inside
    ! the W^2 step the source rule s^2 [pair] + s (1-s) [rho(k) F] G + (1-s) s F [rho(k) G] + (1-s)^2 F G. No A b cut
    ! (as the Python with rho). Per node only the mu2 it contributes to (A^2 <= A2_mu mu2) are computed
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nb,nmu
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nb,nmu),INTENT(OUT)::Nrow
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::e,e2,De
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::G,G2
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::Ph,Gpm
    REAL(KIND(1d0)),DIMENSION(nmu)::kcv
    INTEGER,DIMENSION(nmu)::ima
    INTEGER::nma,ma
    REAL(KIND(1d0)),DIMENSION(6,nmu)::Ia,Ir
    REAL(KIND(1d0)),DIMENSION(6,6)::Pm
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Ivm
    INTEGER,DIMENSION(14)::pp,qq
    REAL(KIND(1d0))::A,hi0,sw,Bs,kap,Dmax,c0,h,Dn(NLV),pre,bt,Dlo
    INTEGER::i,j,m,k,np,np2,ndp,npr,p,q,jj
    npr=0
    DO p=1,6
       DO q=1,6
          IF(row_nj(p).EQ.row_nj(q))THEN
             npr=npr+1; pp(npr)=p; qq(npr)=q
          ENDIF
       ENDDO
    ENDDO
    Nrow=0d0
    DO i=1,NMX
       A=DSQRT(A2(x,mx_nodes(i)))
       hi0=kT_hi(x,mx_nodes(i),HUGE(1d0))
       nma=0
       DO m=1,nmu
          IF(mu2v(m).LT.HUGE(1d0).AND.A**2.GT.A2_mu*mu2v(m))CYCLE
          nma=nma+1; ima(nma)=m; kcv(nma)=kT_hi(x,mx_nodes(i),mu2v(m))
       ENDDO
       IF(nma.EQ.0)CYCLE
       sw=W_step(mx_nodes(i))
       Bs=B_slope(xB(mx_nodes(i),A**2/(1d0-x))); kap=A
       CALL node_panels_mu(x,mx_nodes(i),hi0,nma,kcv(1:nma),0d0,e,np,G)
       IF(sw.GT.0d0)THEN
          Dmax=DSQRT(16d0*DLOG(10d0)/Bs)
          Dlo=rho_dlo*MIN(kap,Dmax)
          ndp=INT(rho_dens*LOG10(Dmax/Dlo))+2
          ALLOCATE(De(ndp+2))
          De(1)=0d0
          DO k=0,ndp-1
             De(k+2)=Dlo*(Dmax/Dlo)**(DBLE(k)/(ndp-1))
          ENDDO
          ALLOCATE(Ph(6,6,NLV,nma),Gpm(NLV,npr,ndp,nma),Ivm(npr,nma))
          DO k=1,ndp
             c0=0.5d0*(De(k)+De(k+1)); h=0.5d0*(De(k+1)-De(k))
             Dn=c0+h*tlv
             IF(Dn(1).LE.0d0)Dn(1)=0d0
             CALL rho_phat_mu(nma,kcv(1:nma),hi0,NLV,Dn,e,np,G,Ph)
             DO m=1,nma
                DO jj=1,NLV
                   DO q=1,npr
                      Gpm(jj,q,k,m)=Dn(jj)*DEXP(-Bs*Dn(jj)**2)*Ph(pp(q),qq(q),jj,m)/twopi
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       IF(sw.GT.0d0.AND.sw.LT.1d0)CALL node_panels_mu(x,mx_nodes(i),hi0,nma,kcv(1:nma),Bs,e2,np2,G2)
       DO j=1,nb
          bt=bv(j)/(1d0-x)
          IF(sw.LT.1d0)CALL amps_mu(bt,nma,kcv(1:nma),e,np,G,Ia)
          IF(sw.GT.0d0.AND.sw.LT.1d0)CALL amps_mu(bt,nma,kcv(1:nma),e2,np2,G2,Ir)
          IF(ALLOCATED(Gpm))CALL lev_transform_mu(bt,De,ndp,Gpm,npr,nma,Ivm) ! the rho pairs of all mu
          DO ma=1,nma
             m=ima(ma)
             IF(.NOT.ALLOCATED(Gpm))THEN ! sw = 0
                DO q=1,6
                   Pm(:,q)=Ia(:,ma)*Ia(q,ma)
                ENDDO
             ELSE
                Pm=0d0
                DO q=1,npr
                   Pm(pp(q),qq(q))=Ivm(q,ma)
                ENDDO
                IF(sw.LT.1d0)THEN
                   Pm=sw**2*Pm
                   DO q=1,6
                      Pm(:,q)=Pm(:,q)+sw*Ir(:,ma)*(1d0-sw)*Ia(q,ma)+(1d0-sw)*Ia(:,ma)*sw*Ir(q,ma)+(1d0-sw)**2*Ia(:,ma)*Ia(q,ma)
                   ENDDO
                ENDIF
             ENDIF
             pre=alpha_norm(mu2v(m))/pi**2/(1d0-x)**2
             Nrow(j,m)=Nrow(j,m)+mx_weights(i)*pre*combination(Pm,x,mx_nodes(i))
          ENDDO
       ENDDO
       IF(ALLOCATED(De))DEALLOCATE(De)
       IF(ALLOCATED(Ph))DEALLOCATE(Ph,Gpm,Ivm)
       IF(ALLOCATED(e2))DEALLOCATE(e2,G2)
       DEALLOCATE(e,G)
    ENDDO
  END SUBROUTINE rho_row


  SUBROUTINE inel_rhohat(x,nq,qv,nmu,mu2v,R)
    ! (M11) the two-dimensional Fourier transform of the b-space density at one x: R(iq,m) = int d^2b exp(-i q.b)
    ! x f/alpha(b) at q = qv(iq) (GeV) and mu2v(m), every rho_model: per mx node w alpha_norm/pi^2 combination of
    ! the momentum-space overlaps P^(Delta) at Delta = q (1-x) (rho_phat_mu), times exp(-B Delta^2) in the
    ! continuum, the source rule inside the W^2 step (rho /= delta); rho = delta: P^ for every node
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::nq,nmu
    REAL(KIND(1d0)),DIMENSION(nq),INTENT(IN)::qv
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::mu2v
    REAL(KIND(1d0)),DIMENSION(nq,nmu),INTENT(OUT)::R
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::e,e2
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::G,G2
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::Ph,Ph2,Phs
    REAL(KIND(1d0)),DIMENSION(nq)::Dv
    REAL(KIND(1d0)),DIMENSION(nmu)::kcv
    INTEGER,DIMENSION(nmu)::ima
    REAL(KIND(1d0)),DIMENSION(6,6)::Pm
    REAL(KIND(1d0))::A,hi0,sw,Bs,pre
    INTEGER::i,m,ma,nma,np,np2,iq
    CALL inel_init()
    IF(.NOT.lv_init)CALL levin_setup()
    Dv=qv*(1d0-x)
    R=0d0
    DO i=1,NMX
       A=DSQRT(A2(x,mx_nodes(i)))
       hi0=kT_hi(x,mx_nodes(i),HUGE(1d0))
       nma=0
       DO m=1,nmu
          IF(mu2v(m).LT.HUGE(1d0).AND.A**2.GT.A2_mu*mu2v(m))CYCLE
          nma=nma+1; ima(nma)=m; kcv(nma)=kT_hi(x,mx_nodes(i),mu2v(m))
       ENDDO
       IF(nma.EQ.0)CYCLE
       sw=0d0; Bs=0d0
       IF(rho_model.NE.0)THEN
          sw=W_step(mx_nodes(i))
          Bs=B_slope(xB(mx_nodes(i),A**2/(1d0-x)))
       ENDIF
       CALL node_panels_mu(x,mx_nodes(i),hi0,nma,kcv(1:nma),0d0,e,np,G)
       ALLOCATE(Ph(6,6,nq,nma),Ph2(6,6,nq,nma),Phs(6,6,nq,nma))
       Ph2=0d0; Phs=0d0
       IF(rhohat_fast)THEN
          CALL rho_phat_fast(nma,kcv(1:nma),hi0,A,nq,Dv,e,np,G,Ph)
       ELSE
          CALL rho_phat_mu(nma,kcv(1:nma),hi0,nq,Dv,e,np,G,Ph)
       ENDIF
       IF(sw.GT.0d0.AND.sw.LT.1d0)THEN
          ! the mixed products [rho(k) F] G + F [rho(k) G]: P^[G+G2] - P^[G] - P^[G2] (bilinear), G2 = G exp(-B k^2)
          CALL node_panels_mu(x,mx_nodes(i),hi0,nma,kcv(1:nma),Bs,e2,np2,G2)
          IF(rhohat_fast)THEN
             CALL rho_phat_fast(nma,kcv(1:nma),hi0,A,nq,Dv,e2,np2,G2,Ph2)
             G2=G2+G
             CALL rho_phat_fast(nma,kcv(1:nma),hi0,A,nq,Dv,e2,np2,G2,Phs)
          ELSE
             CALL rho_phat_mu(nma,kcv(1:nma),hi0,nq,Dv,e2,np2,G2,Ph2)
             G2=G2+G
             CALL rho_phat_mu(nma,kcv(1:nma),hi0,nq,Dv,e2,np2,G2,Phs)
          ENDIF
       ENDIF
       DO ma=1,nma
          m=ima(ma)
          pre=mx_weights(i)*alpha_norm(mu2v(m))/pi**2
          DO iq=1,nq
             IF(sw.LE.0d0)THEN
                Pm=Ph(:,:,iq,ma)
             ELSEIF(sw.GE.1d0)THEN
                Pm=DEXP(-Bs*Dv(iq)**2)*Ph(:,:,iq,ma)
             ELSE
                Pm=sw**2*DEXP(-Bs*Dv(iq)**2)*Ph(:,:,iq,ma)+sw*(1d0-sw)*(Phs(:,:,iq,ma)-Ph(:,:,iq,ma)-Ph2(:,:,iq,ma))&
                     +(1d0-sw)**2*Ph(:,:,iq,ma)
             ENDIF
             R(iq,m)=R(iq,m)+pre*combination(Pm,x,mx_nodes(i))
          ENDDO
       ENDDO
       DEALLOCATE(Ph,Ph2,Phs,e,G)
       IF(ALLOCATED(G2))DEALLOCATE(e2,G2)
    ENDDO
  END SUBROUTINE inel_rhohat


  SUBROUTINE rho_phat_fast(nmu,kcv,hi0,Anode,nD,Dv,e,np,G,Ph)
    ! (M11) rho_phat_mu with the outer k1 integral of R(p,q) = int dk1 f_p(k1) Phi_q(k1,D) split as
    !   pi int_0^kc f_p f_q/k                                  (Parseval half, cumulative per node)
    ! + int_0^kL f_p (Phi_q - pi f_q/k1)                       (exact: inner_phi_cc, kL = min(KF D, kc_min - D))
    ! + D int_kL^(kc-D) f_p (2 f_q' - f_q/k)/k
    ! + D^2 pi/4 int_kL^(kc-D) f_p (k f_q'' - f_q' + (1-n^2) f_q/k)/k^2   (the inner integral to O((D/k1)^2))
    ! + int_(kc-D)^kc f_p (Phi_q^(kc) - pi f_q/k1)             (exact band of each cut, k1 = a + (b-a) s^2)
    ! Outer panels of the exact part: geometric NGO per decade from 1e-3 min(D, A), 0.5 D (square-root onset); below
    ! 0.1 D the inner integral interpolated from 6 Chebyshev points (M12); D < A/KA: the expansion throughout.
    ! Bands: kc >= 300 D leading order -2 D f_p f_q/kc, 30 D <= kc < 300 D with the D^2 term from one exact band
    ! per cut (M12), else exact (Gauss-Legendre 5). The cuts with kc < 3 D (and D <= 0): rho_phat_mu
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nmu,nD,np
    REAL(KIND(1d0)),DIMENSION(nmu),INTENT(IN)::kcv
    REAL(KIND(1d0)),INTENT(IN)::hi0,Anode
    REAL(KIND(1d0)),DIMENSION(nD),INTENT(IN)::Dv
    REAL(KIND(1d0)),DIMENSION(np+1),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(NLV,6,np),INTENT(IN)::G
    REAL(KIND(1d0)),DIMENSION(6,6,nD,nmu),INTENT(OUT)::Ph
    REAL(KIND(1d0)),PARAMETER::KF=20d0,KA=50d0
    INTEGER,PARAMETER::NGO=6
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::Pcum ! (6,6,3,np+1) prefix sums of the panel integrals
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::Gd1,Gd2
    REAL(KIND(1d0)),DIMENSION(6,6,3)::c1,c2
    REAL(KIND(1d0)),DIMENSION(6,6)::Rs,Rm
    REAL(KIND(1d0)),DIMENSION(6)::f1,Phi
    REAL(KIND(1d0)),DIMENSION(6,1)::Phic
    REAL(KIND(1d0)),DIMENSION(nmu)::kc
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::Phs
    INTEGER,DIMENSION(nmu)::isl,ifa
    INTEGER::nsl,nfa,ms
    REAL(KIND(1d0))::kcf
    REAL(KIND(1d0)),DIMENSION(nmu)::kLm
    REAL(KIND(1d0)),DIMENSION(6,6,nmu)::Rsn
    REAL(KIND(1d0))::D,kL,klo,l0,l1,a,b,k1,w1,sj,h
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::oe
    INTEGER::id,k,j,m,no,nn,ipr,npr,hint,ip
    INTEGER,DIMENSION(14)::pp,qq
    LOGICAL::sqp
    ! (M12) cumulative kernels at each cut and the integrands there (independent of D); the region k1 < KSM D of
    ! the exact part from the inner integral at NS Chebyshev points of [0, KSM D] (smooth in k1 there: k2 within
    ! D (1 +- KSM)), interpolated
    REAL(KIND(1d0)),PARAMETER::KSM=0.1d0
    INTEGER,PARAMETER::NS=6
    REAL(KIND(1d0)),DIMENSION(6,6,3,nmu)::Ckc
    REAL(KIND(1d0)),DIMENSION(6,nmu)::fkc
    REAL(KIND(1d0)),DIMENSION(6,6,3)::cKF
    REAL(KIND(1d0)),DIMENSION(NS)::tsn,lw
    REAL(KIND(1d0)),DIMENSION(6,NS)::PhiS
    REAL(KIND(1d0))::ks
    INTEGER::js,jl
    ! (M12) bands of the cuts with rhohat_kb2 D <= kc < rhohat_kb D: -2 D f_p f_q/kc + c2 D^2, c2 from one exact band
    ! at D0 = kc/100 per cut
    REAL(KIND(1d0)),DIMENSION(6,6,nmu)::c2b
    LOGICAL,DIMENSION(nmu)::c2set
    REAL(KIND(1d0)),DIMENSION(5),PARAMETER::t5=(/-0.906179845938663992797626878299d0,-0.538469310105683091036314420700d0,&
         0d0,0.538469310105683091036314420700d0,0.906179845938663992797626878299d0/)
    REAL(KIND(1d0)),DIMENSION(5),PARAMETER::w5=(/0.236926885056189087514264040720d0,0.478628670499366468041291514836d0,&
         0.568888888888888888888888888889d0,0.478628670499366468041291514836d0,0.236926885056189087514264040720d0/)
    kc=MIN(kcv,hi0)
    npr=0
    DO ip=1,6
       DO k=1,6
          IF(row_nj(ip).EQ.row_nj(k))THEN
             npr=npr+1; pp(npr)=ip; qq(npr)=k
          ENDIF
       ENDDO
    ENDDO
    ! derivatives at the Lobatto points and the prefix sums of the three kernels
    ALLOCATE(Gd1(NLV,6,np),Gd2(NLV,6,np),Pcum(6,6,3,np+1))
    DO k=1,np
       h=0.5d0*(e(k+1)-e(k))
       Gd1(:,:,k)=MATMUL(Dlv,G(:,:,k))/h
       Gd2(:,:,k)=MATMUL(D2lv,G(:,:,k))/h**2
    ENDDO
    Pcum(:,:,:,1)=0d0
    DO k=1,np
       CALL panel_ints(k,e(k),e(k+1),c1)
       Pcum(:,:,:,k+1)=Pcum(:,:,:,k)+c1
    ENDDO
    DO m=1,nmu
       CALL cum_at(kc(m),Ckc(:,:,:,m))
       CALL cheb_eval(kc(m),e,np,G,fkc(:,m))
    ENDDO
    c2set=.FALSE.
    DO id=1,nD
       D=Dv(id)
       nsl=0; nfa=0
       DO m=1,nmu
          IF(D.LE.0d0.OR.kc(m).LT.rhohat_kf*D)THEN
             nsl=nsl+1; isl(nsl)=m
          ELSE
             nfa=nfa+1; ifa(nfa)=m
          ENDIF
       ENDDO
       IF(nsl.GT.0)THEN
          ALLOCATE(Phs(6,6,1,nsl))
          ! the largest of these cuts as the upper limit (the same values: k2 <= kc for every cut passed)
          CALL rho_phat_mu(nsl,kc(isl(1:nsl)),MAXVAL(kc(isl(1:nsl))),1,Dv(id:id),e,np,G,Phs,&
               COARSE=1d-3*MIN(D,Anode,MINVAL(kc(isl(1:nsl)))))
          DO ms=1,nsl
             Ph(:,:,id,isl(ms))=Phs(:,:,1,ms)
          ENDDO
          DEALLOCATE(Phs)
       ENDIF
       IF(nfa.EQ.0)CYCLE
       ! per fast cut the end of its exact part kL(m) = min(KF D, kc - D): one pass over [0, max kL] with the
       ! exact integral kept at each kL(m)
       kcf=0d0
       DO ms=1,nfa
          kLm(ms)=MIN(KF*D,kc(ifa(ms))-D); kcf=MAX(kcf,kLm(ms))
       ENDDO
       kL=kcf
       klo=1d-3*MIN(D,Anode)
       IF(ALLOCATED(oe))DEALLOCATE(oe)
       ALLOCATE(oe(NGO*CEILING(MAX(LOG10(MAX(kL,klo)/klo),1d0))+nfa+8))
       no=1; oe(1)=0d0
       IF(klo.LT.kL)THEN
          l0=LOG10(klo); l1=LOG10(kL)
          nn=MAX(1,CEILING(NGO*(l1-l0)))
          DO k=0,nn
             no=no+1; oe(no)=10d0**(l0+(l1-l0)*k/nn)
          ENDDO
       ELSE
          no=no+1; oe(no)=kL
       ENDIF
       IF(0.5d0*D.LT.kL)THEN
          no=no+1; oe(no)=0.5d0*D
       ENDIF
       ks=KSM*D
       IF(ks.LT.kL.AND.ks.GT.klo)THEN
          no=no+1; oe(no)=ks
          DO js=1,NS ! Chebyshev-Gauss points of [0,ks] and the inner integrals there
             tsn(js)=0.5d0*ks*(1d0-DCOS((2*js-1)*pi/(2*NS)))
             CALL inner_phi_cc(tsn(js),D,hi0,0,kc(1:0),e,np,G,PhiS(:,js),Phic(:,1:0),.TRUE.)
          ENDDO
       ELSE
          ks=0d0
       ENDIF
       DO ms=1,nfa
          no=no+1; oe(no)=kLm(ms)
       ENDDO
       CALL sort_unique(oe,no,nn)
       Rs=0d0; hint=0
       ! D << A: the expansion in D throughout (no exact part)
       IF(KA*D.LT.Anode)THEN
          nn=1; kLm(1:nfa)=0d0
       ENDIF
       DO ms=1,nfa
          IF(kLm(ms).LE.0d0)Rsn(:,:,ms)=0d0
       ENDDO
       DO k=1,nn-1
          a=oe(k); b=oe(k+1)
          sqp=ABS(a-0.5d0*D).LE.1d-13*D
          DO j=1,8
             IF(sqp)THEN
                sj=0.5d0*(1d0+t8(j)); k1=a+(b-a)*sj**2; w1=(b-a)*sj*w8(j)
             ELSE
                k1=0.5d0*(a+b)+0.5d0*(b-a)*t8(j); w1=0.5d0*(b-a)*w8(j)
             ENDIF
             CALL cheb_eval(k1,e,np,G,f1,hint)
             IF(k1.LT.ks)THEN
                DO js=1,NS
                   lw(js)=1d0
                   DO jl=1,NS
                      IF(jl.NE.js)lw(js)=lw(js)*(k1-tsn(jl))/(tsn(js)-tsn(jl))
                   ENDDO
                ENDDO
                Phi=MATMUL(PhiS,lw)
             ELSE
                CALL inner_phi_cc(k1,D,hi0,0,kc(1:0),e,np,G,Phi,Phic(:,1:0),.TRUE.)
             ENDIF
             DO ipr=1,npr
                Rs(pp(ipr),qq(ipr))=Rs(pp(ipr),qq(ipr))+w1*f1(pp(ipr))*(Phi(qq(ipr))-pi*f1(qq(ipr))/k1)
             ENDDO
          ENDDO
          DO ms=1,nfa
             IF(ABS(b-kLm(ms)).LE.1d-13*b)Rsn(:,:,ms)=Rs
          ENDDO
       ENDDO
       IF(KA*D.LT.Anode)THEN
          cKF=0d0
       ELSE
          CALL cum_at(KF*D,cKF)
       ENDIF
       DO ms=1,nfa
          m=ifa(ms)
          IF(kLm(ms).LE.0d0)THEN
             c1=0d0
          ELSEIF(ABS(kLm(ms)-KF*D).LE.1d-13*kLm(ms))THEN
             c1=cKF
          ELSE
             CALL cum_at(kLm(ms),c1)
          ENDIF
          Rm=Rsn(:,:,ms)+pi*Ckc(:,:,1,m)
          IF(kc(m).GE.rhohat_kb*D)THEN ! the tail kernels at kc - D as at kc: O((D/kc)^2)
             c2=Ckc(:,:,:,m)
          ELSE
             CALL cum_at(kc(m)-D,c2)
          ENDIF
          Rm=Rm+D*(c2(:,:,2)-c1(:,:,2))+D**2*0.25d0*pi*(c2(:,:,3)-c1(:,:,3))
          ! band [kc-D, kc): kc >= rhohat_kb D: -2 D f_p f_q/kc (leading order in D/kc), else exact
          IF(kc(m).GE.rhohat_kb*D)THEN
             DO ipr=1,npr
                Rm(pp(ipr),qq(ipr))=Rm(pp(ipr),qq(ipr))-2d0*D*fkc(pp(ipr),m)*fkc(qq(ipr),m)/kc(m)
             ENDDO
             Ph(:,:,id,m)=Rm+TRANSPOSE(Rm)
             CYCLE
          ENDIF
          IF(kc(m).GE.rhohat_kb2*D)THEN
             IF(.NOT.c2set(m))THEN
                CALL band_exact(m,0.01d0*kc(m),Rs)
                DO ipr=1,npr
                   c2b(pp(ipr),qq(ipr),m)=(Rs(pp(ipr),qq(ipr))+2d0*0.01d0*kc(m)*fkc(pp(ipr),m)*fkc(qq(ipr),m)/kc(m))&
                        /(0.01d0*kc(m))**2
                ENDDO
                c2set(m)=.TRUE.
             ENDIF
             DO ipr=1,npr
                Rm(pp(ipr),qq(ipr))=Rm(pp(ipr),qq(ipr))-2d0*D*fkc(pp(ipr),m)*fkc(qq(ipr),m)/kc(m)&
                     +c2b(pp(ipr),qq(ipr),m)*D**2
             ENDDO
          ELSE
             CALL band_exact(m,D,Rs)
             Rm=Rm+Rs
          ENDIF
          Ph(:,:,id,m)=Rm+TRANSPOSE(Rm)
       ENDDO
    ENDDO
    DEALLOCATE(Gd1,Gd2,Pcum)
  CONTAINS
    SUBROUTINE band_exact(mm,DD,Rb)
      ! int_(kc-DD)^kc f_p (Phi_q^(kc) - pi f_q/k1) dk1 of cut mm (Gauss-Legendre 5 in s, k1 = a + (b-a) s^2)
      INTEGER,INTENT(IN)::mm
      REAL(KIND(1d0)),INTENT(IN)::DD
      REAL(KIND(1d0)),DIMENSION(6,6),INTENT(OUT)::Rb
      REAL(KIND(1d0)),DIMENSION(6)::fb,Pb
      REAL(KIND(1d0)),DIMENSION(6,1)::Pc
      REAL(KIND(1d0))::aa,bb,kk,ww,ss
      INTEGER::jj,ii
      Rb=0d0
      aa=kc(mm)-DD; bb=kc(mm)
      DO jj=1,5
         ss=0.5d0*(1d0+t5(jj)); kk=aa+(bb-aa)*ss**2; ww=(bb-aa)*ss*w5(jj)
         CALL cheb_eval(kk,e,np,G,fb)
         IF(kc(mm).LT.hi0)THEN
            CALL inner_phi_cc(kk,DD,hi0,1,kc(mm:mm),e,np,G,Pb,Pc)
            Pb=Pc(:,1)
         ELSE
            CALL inner_phi_cc(kk,DD,hi0,0,kc(1:0),e,np,G,Pb,Pc(:,1:0))
         ENDIF
         DO ii=1,npr
            Rb(pp(ii),qq(ii))=Rb(pp(ii),qq(ii))+ww*fb(pp(ii))*(Pb(qq(ii))-pi*fb(qq(ii))/kk)
         ENDDO
      ENDDO
    END SUBROUTINE band_exact
    SUBROUTINE panel_ints(kp,ua,ub,cc)
      ! the three kernels integrated over [ua,ub] inside panel kp (Gauss-Legendre 8)
      INTEGER,INTENT(IN)::kp
      REAL(KIND(1d0)),INTENT(IN)::ua,ub
      REAL(KIND(1d0)),DIMENSION(6,6,3),INTENT(OUT)::cc
      REAL(KIND(1d0)),DIMENSION(6)::f,fd,fdd
      REAL(KIND(1d0))::kk,ww,tt
      INTEGER::jj,ii,p,q,nq
      cc=0d0
      IF(ub.LE.ua)RETURN
      DO jj=1,8
         kk=0.5d0*(ua+ub)+0.5d0*(ub-ua)*t8(jj); ww=0.5d0*(ub-ua)*w8(jj)
         tt=(2d0*kk-e(kp)-e(kp+1))/(e(kp+1)-e(kp))
         DO ii=1,6
            f(ii)=bary(tt,G(:,ii,kp)); fd(ii)=bary(tt,Gd1(:,ii,kp)); fdd(ii)=bary(tt,Gd2(:,ii,kp))
         ENDDO
         DO ii=1,npr
            p=pp(ii); q=qq(ii); nq=row_nj(q)
            cc(p,q,1)=cc(p,q,1)+ww*f(p)*f(q)/kk
            cc(p,q,2)=cc(p,q,2)+ww*f(p)*(2d0*fd(q)-f(q)/kk)/kk
            cc(p,q,3)=cc(p,q,3)+ww*f(p)*(kk*fdd(q)-fd(q)+(1-nq**2)*f(q)/kk)/kk**2
         ENDDO
      ENDDO
    END SUBROUTINE panel_ints
    SUBROUTINE cum_at(ku,cc)
      ! the three kernels integrated over [0,ku]
      REAL(KIND(1d0)),INTENT(IN)::ku
      REAL(KIND(1d0)),DIMENSION(6,6,3),INTENT(OUT)::cc
      REAL(KIND(1d0)),DIMENSION(6,6,3)::cp
      INTEGER::lo,hi,mid
      IF(ku.LE.e(1))THEN
         cc=0d0; RETURN
      ENDIF
      IF(ku.GE.e(np+1))THEN
         cc=Pcum(:,:,:,np+1); RETURN
      ENDIF
      lo=1; hi=np+1
      DO WHILE(hi-lo.GT.1)
         mid=(lo+hi)/2
         IF(ku.GE.e(mid))THEN
            lo=mid
         ELSE
            hi=mid
         ENDIF
      ENDDO
      CALL panel_ints(lo,e(lo),ku,cp)
      cc=Pcum(:,:,:,lo)+cp
    END SUBROUTINE cum_at
  END SUBROUTINE rho_phat_fast


  SUBROUTINE elastic_rhohat(x,iFF,nq,qv,R)
    ! (M11) the two-dimensional Fourier transform of the exact elastic density of the proton (alpha(0) normalisation,
    ! as PhotonNumberDensity_ChargeFormFactor_proton_exact) at q = qv (GeV): exact_combination of the pair overlaps
    ! P^(q (1-x)) of its six k_T integrands (rho_phat_fast), divided by pi^2. Panels: geometric 10/decade from
    ! 1e-3 x mp to k_T at Q2 = 1e4 GeV^2, the joins of the form factors
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::iFF,nq
    REAL(KIND(1d0)),DIMENSION(nq),INTENT(IN)::qv
    REAL(KIND(1d0)),DIMENSION(nq),INTENT(OUT)::R
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::t,ej
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::G
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::Ph
    REAL(KIND(1d0)),DIMENSION(6,NLV)::V
    REAL(KIND(1d0))::A,hi,klo,l0,l1,c0,h
    INTEGER::k,ngeo,ne,np,iq
    IF(.NOT.lv_init)CALL levin_setup()
    R=0d0
    IF(x.GE.1d0)RETURN
    A=x*mp
    hi=DSQRT(MAX(1d4*(1d0-x)-A**2,1d-6))
    klo=1d-3*MIN(A,hi)
    ngeo=INT(10*LOG10(hi/klo))+2
    ej=proton_exact_edges(x,iFF)
    ALLOCATE(t(ngeo+SIZE(ej)+2))
    l0=LOG10(klo); l1=LOG10(hi)
    DO k=0,ngeo-1
       t(k+1)=10d0**(l0+k*((l1-l0)/(ngeo-1)))
    ENDDO
    ne=ngeo+1; t(ne)=0d0
    DO k=1,SIZE(ej)
       IF(ej(k).LT.hi)THEN
          ne=ne+1; t(ne)=ej(k)
       ENDIF
    ENDDO
    CALL sort_unique(t,ne,np)
    np=np-1
    ALLOCATE(G(NLV,6,np),Ph(6,6,nq,1))
    DO k=1,np
       c0=0.5d0*(t(k)+t(k+1)); h=0.5d0*(t(k+1)-t(k))
       CALL proton_exact_kspace(x,iFF,NLV,c0+h*tlv,V)
       G(:,:,k)=TRANSPOSE(V)
    ENDDO
    CALL rho_phat_fast(1,(/hi/),hi,A,nq,qv*(1d0-x),t(1:np+1),np,G,Ph)
    DO iq=1,nq
       R(iq)=proton_exact_comb_P(x,Ph(:,:,iq,1))/pi**2
    ENDDO
  END SUBROUTINE elastic_rhohat

  SUBROUTINE nucleus_rhohat(x,RR,w,aa,nq,qv,R)
    ! (M15) the two-dimensional Fourier transform of the elastic photon density of a nucleus with the Woods-Saxon
    ! charge form factor F (unit charge, alpha = 1, as PhotonNumberDensity_ChargeFormFactor_WS) at q = qv (GeV): the
    ! overlap P^(q) of its k_T integrand k^2 F(sqrt(k^2+A^2))/(k^2+A^2), A = x mN, divided by pi^2. RR, aa (the
    ! Woods-Saxon radius and diffuseness) in GeV^-1. Panels: geometric 10/decade from 1e-3 min(A, 1/RR) to the cut
    ! 1e3 GeV, and of width pi/(4 RR) (F oscillates with period ~2 pi/RR) from 1/RR to 8/aa (F < e^-25 beyond: the
    ! integrand is 0 there). The overlaps by rho_phat_mu (about 2 ms per q): rho_phat_fast expands the inner integral
    ! in q/k beyond k = 20 q, which assumes that f varies on the scale k (the proton); F varies on the scale 1/RR and
    ! the expansion is off by up to 5e-6 at q ~ 1e-3 GeV (rho_phat_mu agrees with a direct two-dimensional quadrature
    ! to 1e-10: drivers/m15_nuc_rho, scripts/m15_nuc_rho_ref.py)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,RR,w,aa
    INTEGER,INTENT(IN)::nq
    REAL(KIND(1d0)),DIMENSION(nq),INTENT(IN)::qv
    REAL(KIND(1d0)),DIMENSION(nq),INTENT(OUT)::R
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! as PhotonPhotonFlux_pA_WoodsSaxon_eval
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::t
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::G
    REAL(KIND(1d0)),DIMENSION(:,:,:,:),ALLOCATABLE::Ph
    REAL(KIND(1d0))::A,hi,klo,l0,l1,c0,h,k1,k2,kk
    INTEGER::k,i,ngeo,nlin,ne,np
    IF(.NOT.lv_init)CALL levin_setup()
    R=0d0
    IF(x.GE.1d0.OR.x.LE.0d0)RETURN
    A=x*mN
    hi=1d3
    klo=1d-3*MIN(A,1d0/RR)
    ngeo=INT(10*LOG10(hi/klo))+2
    k1=1d0/RR; k2=MAX(8d0/aa,k1)
    nlin=CEILING((k2-k1)/(pi/(4d0*RR)))
    ALLOCATE(t(ngeo+nlin+2))
    l0=LOG10(klo); l1=LOG10(hi)
    DO k=0,ngeo-1
       t(k+1)=10d0**(l0+k*((l1-l0)/(ngeo-1)))
    ENDDO
    ne=ngeo+1; t(ne)=0d0
    DO k=0,nlin
       ne=ne+1; t(ne)=k1+(k2-k1)*k/nlin
    ENDDO
    CALL sort_unique(t,ne,np)
    np=np-1
    ALLOCATE(G(NLV,6,np),Ph(6,6,nq,1))
    G=0d0
    DO k=1,np
       IF(t(k).GE.k2)EXIT
       c0=0.5d0*(t(k)+t(k+1)); h=0.5d0*(t(k+1)-t(k))
       DO i=1,NLV
          kk=c0+h*tlv(i)
          G(i,1,k)=kk**2/(kk**2+A**2)*ChargeFormFactor_WoodsSaxon(DSQRT(kk**2+A**2),RR,w,aa,10)
       ENDDO
    ENDDO
    CALL rho_phat_mu(1,(/hi/),hi,nq,qv,t(1:np+1),np,G,Ph)
    R=Ph(1,1,:,1)/pi**2
  END SUBROUTINE nucleus_rhohat

  SUBROUTINE cheb_setup(n,xc,wt,cm)
    ! n Chebyshev-Lobatto points xc(i) = cos(pi (i-1)/(n-1)) on [-1,1], the matrix cm of values -> coefficients of
    ! sum_k c_k T_(k-1), and wt(k) = int_{-1}^{1} T_(k-1)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::xc,wt
    REAL(KIND(1d0)),DIMENSION(n,n),INTENT(OUT)::cm
    INTEGER::i,k,m
    m=n-1
    DO i=1,n
       xc(i)=DCOS(pi*(i-1)/m)
    ENDDO
    DO k=1,n
       DO i=1,n
          cm(k,i)=2d0/m*DCOS(pi*(k-1)*(i-1)/m)
          IF(i.EQ.1.OR.i.EQ.n)cm(k,i)=0.5d0*cm(k,i)
          IF(k.EQ.1.OR.k.EQ.n)cm(k,i)=0.5d0*cm(k,i)
       ENDDO
       wt(k)=0d0
       IF(MOD(k-1,2).EQ.0)wt(k)=2d0/(1d0-DBLE(k-1)**2)
    ENDDO
  END SUBROUTINE cheb_setup

  SUBROUTINE cheb_antider(n,u,a)
    ! a(k) = int_{-1}^{u} T_(k-1), k = 1..n
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),INTENT(IN)::u
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::a
    REAL(KIND(1d0)),DIMENSION(0:n)::T,Tm
    INTEGER::k
    T(0)=1d0; T(1)=u; Tm(0)=1d0; Tm(1)=-1d0
    DO k=2,n
       T(k)=2d0*u*T(k-1)-T(k-2)
       Tm(k)=-Tm(k-1)
    ENDDO
    a(1)=T(1)-Tm(1)
    IF(n.GE.2)a(2)=0.25d0*(T(2)-Tm(2))
    DO k=2,n-1
       a(k+1)=0.5d0*((T(k+1)-Tm(k+1))/(k+1)-(T(k-1)-Tm(k-1))/(k-1))
    ENDDO
  END SUBROUTINE cheb_antider

  SUBROUTINE inner_phi_cc(k1,D,kf,nc,kcs,e,np,G,Phi,Phic,SHORT2)
    ! Phi_q(k1,D) = int dpsi 2 f_q(k2)/k2 cos(n phi) over k2 >= k1, k2 <= kf (k2^2 = k1^2 + D^2 - 2 k1 D cos psi), and
    ! the same up to the cuts kcs(c) < kf in Phic(:,c), from one set of integrand values: rho_inner_npan uniform panels
    ! in psi of NCL Chebyshev-Lobatto points (one panel of NCS points when D < rho_inner_short k1); a panel by
    ! Clenshaw-Curtis, the panel of a cut by the integral of the Chebyshev interpolant up to it
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::k1,D,kf
    INTEGER,INTENT(IN)::nc,np
    REAL(KIND(1d0)),DIMENSION(nc),INTENT(IN)::kcs
    REAL(KIND(1d0)),DIMENSION(np+1),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(NLV,6,np),INTENT(IN)::G
    REAL(KIND(1d0)),DIMENSION(6),INTENT(OUT)::Phi
    REAL(KIND(1d0)),DIMENSION(6,nc),INTENT(OUT)::Phic
    LOGICAL,INTENT(IN),OPTIONAL::SHORT2 ! (M11) the short rule also for k1 < rho_inner_short D; (M12) NCS points
                                        ! per panel also elsewhere
    REAL(KIND(1d0)),DIMENSION(6,NCL,rho_inner_npan)::gv ! integrand values per panel
    REAL(KIND(1d0)),DIMENSION(6,rho_inner_npan)::Pp ! panel integrals
    REAL(KIND(1d0)),DIMENSION(NCL)::av,wq
    REAL(KIND(1d0))::pa,pb,h,psi,k2,cphi,cn(0:2),f2(6),u,cu,xi,wi
    INTEGER::npan,nn,ip,i,p,jc,iu,hint
    Phi=0d0; Phic=0d0
    pa=DACOS(MIN(1d0,D/(2d0*k1)))
    pb=DACOS(MAX(-1d0,(k1**2+D**2-kf**2)/(2d0*k1*D)))
    IF(pb.LE.pa)RETURN
    IF(D.LT.rho_inner_short*k1)THEN
       npan=1; nn=NCS
    ELSEIF(PRESENT(SHORT2).AND.k1.LT.rho_inner_short*D)THEN
       npan=1; nn=NCS
    ELSEIF(PRESENT(SHORT2))THEN ! (M12) rho_inner_npan panels of NCS points (as NCL to 3e-6 in the transforms)
       npan=rho_inner_npan; nn=NCS
    ELSE
       npan=rho_inner_npan; nn=NCL
    ENDIF
    h=(pb-pa)/npan
    hint=0
    DO ip=1,npan
       Pp(:,ip)=0d0
       DO i=1,nn
          IF(nn.EQ.NCL)THEN
             xi=xcl(i); wi=wvl(i)
          ELSE
             xi=xcs(i); wi=wvs(i)
          ENDIF
          psi=pa+(ip-0.5d0)*h+0.5d0*h*xi
          gv(:,i,ip)=0d0
          k2=DSQRT(MAX(k1**2+D**2-2d0*k1*D*DCOS(psi),0d0))
          IF(k2.LE.0d0)CYCLE
          CALL cheb_eval(k2,e,np,G,f2,hint)
          cphi=(k1**2+k2**2-D**2)/(2d0*k1*k2)
          cn(0)=1d0; cn(1)=cphi; cn(2)=2d0*cphi**2-1d0
          DO p=1,6
             gv(p,i,ip)=2d0*f2(p)/k2*cn(row_nj(p))
             Pp(p,ip)=Pp(p,ip)+0.5d0*h*wi*gv(p,i,ip)
          ENDDO
       ENDDO
       Phi=Phi+Pp(:,ip)
    ENDDO
    DO jc=1,nc
       cu=(k1**2+D**2-kcs(jc)**2)/(2d0*k1*D)
       IF(cu.LE.-1d0)THEN
          Phic(:,jc)=Phi
          CYCLE
       ENDIF
       u=DACOS(MIN(cu,1d0))
       IF(u.LE.pa)CYCLE
       iu=MIN(npan,INT((u-pa)/h)+1)
       DO ip=1,iu-1
          Phic(:,jc)=Phic(:,jc)+Pp(:,ip)
       ENDDO
       ! weights of the values for int_{-1}^{t} of the Chebyshev interpolant: av . cm(:,i)
       CALL cheb_antider(nn,MIN(1d0,2d0*(u-pa-(iu-1)*h)/h-1d0),av(1:nn))
       IF(nn.EQ.NCL)THEN
          DO i=1,nn
             wq(i)=0.5d0*h*DOT_PRODUCT(av(1:nn),cml(:,i))
          ENDDO
       ELSE
          DO i=1,nn
             wq(i)=0.5d0*h*DOT_PRODUCT(av(1:nn),cms(:,i))
          ENDDO
       ENDIF
       DO i=1,nn
          Phic(1,jc)=Phic(1,jc)+wq(i)*gv(1,i,iu); Phic(2,jc)=Phic(2,jc)+wq(i)*gv(2,i,iu)
          Phic(3,jc)=Phic(3,jc)+wq(i)*gv(3,i,iu); Phic(4,jc)=Phic(4,jc)+wq(i)*gv(4,i,iu)
          Phic(5,jc)=Phic(5,jc)+wq(i)*gv(5,i,iu); Phic(6,jc)=Phic(6,jc)+wq(i)*gv(6,i,iu)
       ENDDO
    ENDDO
  END SUBROUTINE inner_phi_cc
  SUBROUTINE lev_weights(bt,a,c,nord,y)
    ! int_a^c g(k) J_nord(bt k) dk = y . g at the NLV Chebyshev-Lobatto points of [a,c]: Clenshaw-Curtis when
    ! bt (c-a) <= lev_theta or a <= 0; else Levin collocation, the integral being linear in g: with S p2 = g (see
    ! lev_row) and I = w . p2, y = S^-T w (one solve for all integrands sharing bt, [a,c] and nord)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::bt,a,c
    INTEGER,INTENT(IN)::nord
    REAL(KIND(1d0)),DIMENSION(NLV),INTENT(OUT)::y
    REAL(KIND(1d0))::c0,h,z,Ja(2),Jc(2),kn(NLV),St(NLV,NLV),Y2(NLV,1)
    INTEGER::jj,q
    c0=0.5d0*(a+c); h=0.5d0*(c-a)
    kn=c0+h*tlv
    IF(bt*2d0*h.LE.lev_theta.OR.a.LE.0d0)THEN
       DO jj=1,NLV
          z=bt*kn(jj)
          y(jj)=h*wcc(jj)*BESSEL_JN(nord,z)
       ENDDO
       RETURN
    ENDIF
    Ja(1)=BESSEL_JN(nord,bt*a); Ja(2)=BESSEL_JN(nord+1,bt*a)
    Jc(1)=BESSEL_JN(nord,bt*c); Jc(2)=BESSEL_JN(nord+1,bt*c)
    ! w: I = Jc(n) pc1 + Jc(n+1) p2(NLV) - Ja(n) pa1 - Ja(n+1) p2(1), p1 = (D/h - (n+1)/k) p2 / bt at the ends
    y=Jc(1)*Dlv(NLV,:)/(h*bt)-Ja(1)*Dlv(1,:)/(h*bt)
    y(NLV)=y(NLV)-Jc(1)*(nord+1)/(kn(NLV)*bt)+Jc(2)
    y(1)=y(1)+Ja(1)*(nord+1)/(kn(1)*bt)-Ja(2)
    ! S^T (S as in lev_row)
    DO jj=1,NLV
       DO q=1,NLV
          St(q,jj)=(D2lv(jj,q)/h**2-Dlv(jj,q)/h*(nord+1)/kn(q)+nord/kn(jj)*Dlv(jj,q)/h)/bt
       ENDDO
       St(jj,jj)=St(jj,jj)-nord*(nord+1)/(kn(jj)**2*bt)+bt
    ENDDO
    Y2(:,1)=y
    CALL solve_lin(St,Y2,NLV,1)
    y=Y2(:,1)
  END SUBROUTINE lev_weights

  SUBROUTINE lev_transform_mu(bt,e,np,Gp,nr,nm,I)
    ! I(r,m) = int dk g_(r,m)(k) J_0(bt k) over the panels e (g at the NLV Chebyshev-Lobatto points in Gp(:,r,panel,m)),
    ! as lev_transform with one set of weights (lev_weights) per panel for all r and m
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::bt
    INTEGER,INTENT(IN)::np,nr,nm
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::e
    REAL(KIND(1d0)),DIMENSION(NLV,nr,np,nm),INTENT(IN)::Gp
    REAL(KIND(1d0)),DIMENSION(nr,nm),INTENT(OUT)::I
    REAL(KIND(1d0)),DIMENSION(NLV)::y
    INTEGER::k,q,m,jj
    REAL(KIND(1d0))::t
    I=0d0
    DO k=1,np
       CALL lev_weights(bt,e(k),e(k+1),0,y)
       DO m=1,nm
          DO q=1,nr
             t=0d0
             DO jj=1,NLV
                t=t+y(jj)*Gp(jj,q,k,m)
             ENDDO
             I(q,m)=I(q,m)+t
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE lev_transform_mu
END MODULE InelasticFlux
