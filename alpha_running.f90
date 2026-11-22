MODULE alpha_running
  ! alpha(Q2)/alpha(0) of a spacelike photon (q2=-Q2), with the on-shell running alpha(Q) of LbLatNLO
  ! (arXiv:2604.00099; evaluate_couplings.f90, ALPHAEW with alpha_scheme=2, alpha_nloop=2, inputs of
  ! run_default.inp: 1/alpha(0)=137.036, 1/alpha(mZ)=128.94) without its two-loop QCD term (switched on there
  ! abruptly at Q=1 GeV, at most 0.22% on the photon densities):
  !   alpha(Q)=alpha(0)/(K-Delta alpha_lep(Q)-Delta alpha_top(Q)-Delta alpha_had(Q)),
  ! one- and two-loop QED vacuum polarisation of leptons and quarks, evaluated at q2=-Q2 (the photon of the flux is
  ! spacelike, as in LUXqed; ALPHAEW is the timelike Re Pi at q2=Q2)
  USE nielsen_generalized_polylog
  IMPLICIT NONE
  PRIVATE
  PUBLIC::alpha_ratio_ph
  REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
  REAL(KIND(1d0)),PARAMETER::zeta3=1.20205690315959428539973816151d0
  REAL(KIND(1d0)),PARAMETER::alphaemm1=137.036d0
  ! K=1-Delta alpha_had(mZ) [from 1/alpha(mZ), Evaluate_DalphahadMZ]+Delta alpha_had^pert(mZ) [u,d masses 83 MeV]
  ! (timelike at mZ, from LbLatNLO: Claude/build/lbl_ref/ref_consts); alpha(0)/alpha_ratio_ph(0)=K
  REAL(KIND(1d0)),PARAMETER::Kconst=1.00002168859830620d0
  ! masses (GeV), charges**2 (quarks), as LbL_Global (umass_PDG=dmass_PDG=83 MeV set by ALPHAEW for alpha_nloop=2)
  REAL(KIND(1d0)),DIMENSION(3),PARAMETER::mlep=(/0.51099895d-3,0.1056583755d0,1.77686d0/)
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::mq=(/83d-3,83d-3,0.490d0,1.5d0,4.5d0,172.69d0/)
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::eq2=(/4d0,1d0,1d0,4d0,1d0,4d0/)/9d0
CONTAINS
  FUNCTION alpha_ratio_ph(Q2)
    ! alpha(Q2)/alpha(0) at q2=-Q2 (Q2 in GeV**2)
    IMPLICIT NONE
    REAL(KIND(1d0))::alpha_ratio_ph
    REAL(KIND(1d0)),INTENT(IN)::Q2
    REAL(KIND(1d0))::D
    INTEGER::i
    D=Kconst
    IF(Q2.GT.0d0)THEN
       DO i=1,3 ! leptons
          D=D+sigma1(mlep(i)**2/Q2)/(3d0*pi*alphaemm1)+sigma2(mlep(i)**2/Q2)/(pi*alphaemm1)**2
       ENDDO
       DO i=1,6 ! quarks (u,d,s,c,b: Delta alpha_had; t: Delta alpha_top)
          D=D+eq2(i)*(sigma1(mq(i)**2/Q2)/(pi*alphaemm1)+3d0*eq2(i)*sigma2(mq(i)**2/Q2)/(pi*alphaemm1)**2)
       ENDDO
    ENDIF
    alpha_ratio_ph=1d0/D
    RETURN
  END FUNCTION alpha_ratio_ph

  FUNCTION sigma1(r)
    ! ReSIGMAAA of LbLatNLO at q2=-Q2: one-loop Re Sigma_AA(q2)-Sigma_AA(0) per fermion over alpha/(3pi)Nc Qf**2,
    ! r=M**2/Q2 (MQ2Q2=-r); closed form, or the expansions for r<=1e-8 and r>=1e5 as in ReSIGMAAA
    IMPLICIT NONE
    REAL(KIND(1d0))::sigma1
    REAL(KIND(1d0)),INTENT(IN)::r
    REAL(KIND(1d0))::v,L,ri
    IF(r.LE.1d-8)THEN
       L=DLOG(r)
       sigma1=L+5d0/3d0-6d0*r+r**2*(3d0-6d0*L)+r**3*(16d0/3d0+8d0*L)-r**4*(33d0/2d0+18d0*L)
    ELSEIF(r.GE.1d5)THEN
       ri=-1d0/r
       sigma1=ri*(1d0/5d0+ri*(3d0/140d0+ri*(1d0/315d0+ri*(1d0/1848d0+ri*(1d0/10010d0+ri*(1d0/51480d0&
            +ri/255255d0))))))
    ELSE
       v=DSQRT(1d0+4d0*r)
       sigma1=5d0/3d0-4d0*r+v*(1d0-2d0*r)*DLOG((v-1d0)/(v+1d0))
    ENDIF
    RETURN
  END FUNCTION sigma1

  FUNCTION sigma2(r)
    ! ReSIGMAAA2L of LbLatNLO at q2=-Q2 (two loops), r=M**2/Q2 (MQ2Q2=-r); phi_n(x)=Li_n(x)+2Li_n(-x) real for
    ! 0<x<1; the expansion for r<=1e-8 is that of ReSIGMAAA2L without the pi**2 (zeta2) of Re log(-r)**2
    IMPLICIT NONE
    REAL(KIND(1d0))::sigma2
    REAL(KIND(1d0)),INTENT(IN)::r
    REAL(KIND(1d0))::v,x,L,phi1,phi2,phi3,ri
    IF(r.LE.1d-8)THEN
       L=DLOG(r)
       sigma2=5d0/24d0+0.25d0*L-zeta3-3d0*r*L+r**2*(1d0/6d0+4d0*zeta3+2.5d0*L-3d0*L**2)&
            -r**3*(248d0/27d0-188d0/27d0*L-58d0/9d0*L**2)
    ELSEIF(r.GE.1d3)THEN
       ri=-1d0/r
       sigma2=ri*(41d0/162d0+ri*(449d0/10800d0+ri*(62479d0/7938000d0+ri*(25993d0/16329600d0&
            +ri*(6756019d0/20170458000d0+ri*338452951d0/4674935865600d0)))))
    ELSE
       v=DSQRT(1d0+4d0*r)
       x=(v-1d0)/(v+1d0)
       L=DLOG(x)
       phi1=-DLOG(1d0-x)-2d0*DLOG(1d0+x)
       phi2=DBLE(Nielsen_PolyLog(1,1,x)+2d0*Nielsen_PolyLog(1,1,-x))
       phi3=DBLE(Nielsen_PolyLog(2,1,x)+2d0*Nielsen_PolyLog(2,1,-x))
       sigma2=(v**4-2d0*v**2-3d0)/12d0*(phi1*L**2-4d0*phi2*L+6d0*phi3+3d0*zeta3)&
            +(3d0*v-v**3)/12d0*(4d0*phi1*L-4d0*phi2+3d0*L**2)&
            +(5d0*v-3d0*v**3)/8d0*L+(7d0*v**4-22d0*v**2-33d0)/96d0*L**2&
            -13d0*(v**2-1d0)/24d0+5d0/24d0
    ENDIF
    RETURN
  END FUNCTION sigma2
END MODULE alpha_running
