! Independent brute-force reference for the AB (and pA) photon-photon flux of gamma-UPC.
!
!   ref_ab <namelist file>
!
! Same physics as the new engine (ChFF with the longitudinal 1/gamma^2 term, ChFF photon b >= 1e-3 R,
! EDFF b >= R or R+a, no truncation in b) but independent numerics:
!  * P_had: T_A(s) by adaptive integration in z of the normalised Woods-Saxon density, T_AB(b) by a
!    direct 2D convolution (adaptive in s and psi), tabulated in ln chi on a fine b grid
!    (the engine uses the Hankel transform of the form factors);
!  * ChFF photon amplitude: the k-space Fourier-Bessel integral of F(sqrt(k^2+q^2))/(k^2+q^2) with
!    F from the Woods-Saxon density (Taylor series at small Q, adaptive radial integral otherwise),
!    tabulated per x on a fine ln b grid (the engine uses the position-space Gauss-law form);
!  * the flux: nested globally adaptive Gauss-Kronrod 10/21 (gk21_adapt) in (ln b1, ln b2, phi),
!    with breakpoints where |b1-b2| or b1+b2 crosses the edge of P (the engine uses a fixed
!    two-dimensional rule in rotated logarithms and an x-independent kernel matrix).
! Neutron tagging uses the engine's P_tag table (validated separately).
MODULE ref_state
  IMPLICIT NONE
  REAL(KIND(1d0)),PARAMETER::GeVm12fm=0.1973d0,mN=0.9315d0,mproton=0.938272081d0
  REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
  ! beams (lengths in GeV^-1)
  LOGICAL,DIMENSION(2)::isp
  LOGICAL::chff=.FALSE.
  REAL(KIND(1d0)),DIMENSION(2)::RA,aa,ww,ZZ,AA_,gam,mass,bmin,rho0
  INTEGER::icoll_r=3
  ! P_had table: ln chi on b = (i-1) hb, P = 1 beyond bcut
  INTEGER::np_t
  REAL(KIND(1d0))::hb,bcut,sigNN
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::lnchi
  REAL(KIND(1d0)),DIMENSION(200)::beta
  INTEGER::nbeta
  ! tagging table on ln b
  LOGICAL::tag=.FALSE.
  INTEGER::nT
  REAL(KIND(1d0))::Tv0,Th
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Ptab
  ! T_A tables (fine grid in s) for the convolution
  INTEGER,DIMENSION(2)::nta
  REAL(KIND(1d0)),DIMENSION(2)::hta
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::tat
  ! ChFF amplitude tables per beam for the current x: A_perp, A_par on v = ln b
  INTEGER,DIMENSION(2)::nva
  REAL(KIND(1d0)),DIMENSION(2)::va0,vah,xcur
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::atab
  ! F(Q) Taylor moments and table
  REAL(KIND(1d0)),DIMENSION(0:12,2)::rmom
  INTEGER,DIMENSION(2)::nq
  REAL(KIND(1d0)),DIMENSION(2)::hq
  REAL(KIND(1d0)),DIMENSION(40001,2)::fq
  ! integration state
  INTEGER::ib_c,ib_cur
  REAL(KIND(1d0))::r_c,b_c,q_c,s_c,bb_c,d_c,G4_c,x_c1,x_c2,v1_c,b1_c
  INTEGER::neval_tot=0,nmid=0,nmid_ier=0,nphi=0,nphi_ier=0,nu1=0
  REAL(KIND(1d0))::reltol_outer=1d-8,wn2=0d0,wn2b=0d0
  LOGICAL::debug=.FALSE.
  REAL(KIND(1d0))::b0pp=0d0     ! pp: elastic slope (GeV^-2) of P = (1-exp(-b^2/(2 b0)))^2
  INTEGER::njump=0
  REAL(KIND(1d0)),DIMENSION(10)::bjump
  INTEGER::nf2=0
  LOGICAL::wdb=.FALSE.   ! <Delta b> numerator: kernel P(b12)*b12 (fm)
END MODULE ref_state

MODULE ref_funcs
  USE ref_state
  USE NINTLIB
  IMPLICIT NONE
CONTAINS
  FUNCTION rho(r,ib)
    REAL(KIND(1d0))::rho
    REAL(KIND(1d0)),INTENT(IN)::r
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0))::arg
    arg=(r-RA(ib))/aa(ib)
    IF(arg.GT.700d0)THEN
       rho=0d0
    ELSE
       rho=rho0(ib)*(1d0+ww(ib)*(r/RA(ib))**2)/(1d0+DEXP(arg))
    ENDIF
  END FUNCTION rho

  FUNCTION f_norm(r)
    REAL(KIND(1d0))::f_norm
    REAL(KIND(1d0)),INTENT(IN)::r
    f_norm=4d0*pi*r*r*rho(r,ib_c)
  END FUNCTION f_norm

  FUNCTION f_mom(r)
    ! 4 pi r^2 rho r^(2k), k = NINT(r_c)
    REAL(KIND(1d0))::f_mom
    REAL(KIND(1d0)),INTENT(IN)::r
    f_mom=4d0*pi*r*r*rho(r,ib_c)*r**(2*NINT(r_c))
  END FUNCTION f_mom

  FUNCTION f_tz(z)
    REAL(KIND(1d0))::f_tz
    REAL(KIND(1d0)),INTENT(IN)::z
    f_tz=2d0*rho(DSQRT(s_c*s_c+z*z),ib_c)
  END FUNCTION f_tz

  FUNCTION TA_direct(s,ib)
    ! thickness function by adaptive integration in z (normalised to 1), GeV^2
    REAL(KIND(1d0))::TA_direct
    REAL(KIND(1d0)),INTENT(IN)::s
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0))::res,err,zmax,pts(1)
    INTEGER::ier,ibs
    ibs=ib_c
    ib_c=ib
    s_c=s
    zmax=RA(ib)+70d0*aa(ib)
    pts(1)=MAX(0d0,DSQRT(MAX(RA(ib)**2-s*s,0d0)))
    CALL gk21_adapt(f_tz,0d0,zmax,0d0,1d-13,res,err,ier,1,pts)
    TA_direct=res
    ib_c=ibs
  END FUNCTION TA_direct

  FUNCTION TA_tab(s,ib)
    ! T_A from the fine table (8-point Lagrange)
    REAL(KIND(1d0))::TA_tab
    REAL(KIND(1d0)),INTENT(IN)::s
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0))::u
    IF(s.GE.hta(ib)*DBLE(nta(ib)-1))THEN
       TA_tab=0d0
       RETURN
    ENDIF
    u=s/hta(ib)
    TA_tab=lag8(u,nta(ib),tat(:,ib))
  END FUNCTION TA_tab

  FUNCTION lag8(u,n,y)
    ! 8-point Lagrange interpolation on nodes 0..n-1 (O(p) weights)
    REAL(KIND(1d0))::lag8
    REAL(KIND(1d0)),INTENT(IN)::u
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(IN)::y
    REAL(KIND(1d0)),DIMENSION(0:7),PARAMETER::rden=(/-1d0/5040d0,1d0/720d0,-1d0/240d0,1d0/144d0,&
         -1d0/144d0,1d0/240d0,-1d0/720d0,1d0/5040d0/)
    REAL(KIND(1d0)),DIMENSION(0:7)::pre,suf
    INTEGER::i0,m
    REAL(KIND(1d0))::v
    i0=MAX(0,MIN(FLOOR(u)-3,n-8))
    v=u-DBLE(i0)
    pre(0)=1d0
    DO m=1,7
       pre(m)=pre(m-1)*(v-DBLE(m-1))
    ENDDO
    suf(7)=1d0
    DO m=6,0,-1
       suf(m)=suf(m+1)*(v-DBLE(m+1))
    ENDDO
    lag8=0d0
    DO m=0,7
       lag8=lag8+pre(m)*suf(m)*rden(m)*y(i0+m+1)
    ENDDO
  END FUNCTION lag8

  FUNCTION f_conv_psi(psi)
    REAL(KIND(1d0))::f_conv_psi
    REAL(KIND(1d0)),INTENT(IN)::psi
    REAL(KIND(1d0))::r2
    r2=DSQRT(MAX(s_c*s_c+bb_c*bb_c-2d0*s_c*bb_c*DCOS(psi),0d0))
    f_conv_psi=TA_tab(r2,2)
  END FUNCTION f_conv_psi

  FUNCTION f_conv_s(s)
    REAL(KIND(1d0))::f_conv_s
    REAL(KIND(1d0)),INTENT(IN)::s
    REAL(KIND(1d0))::res,err,ta
    INTEGER::ier
    ta=TA_tab(s,1)
    IF(ta.LE.1d-16*tat(1,1))THEN
       f_conv_s=0d0
       RETURN
    ENDIF
    s_c=s
    CALL gk21_adapt(f_conv_psi,0d0,pi,1d-14*tat(1,2),1d-11,res,err,ier)
    f_conv_s=2d0*s*ta*res
  END FUNCTION f_conv_s

  SUBROUTINE build_P()
    ! bcut (P = 1 beyond, 1-P < 1e-12) and the breakpoints where P crosses fixed levels
    INTEGER::i,k,j
    REAL(KIND(1d0))::b,b1,b2,bm,amin
    REAL(KIND(1d0)),DIMENSION(14),PARAMETER::lev=(/1d-14,1d-10,1d-6,1d-4,1d-2,0.1d0,0.3d0,0.5d0,&
         0.7d0,0.9d0,0.99d0,0.9999d0,0.999999d0,0.99999999d0/)
    amin=MIN(aa(1),aa(2))
    hb=amin/16d0
    np_t=0
    bcut=1d300
    DO i=INT((RA(1)+RA(2)+80d0*amin)/hb),1,-1
       b=DBLE(i)*hb
       IF(1d0-Phad(b).GT.1d-12)THEN
          bcut=b+hb
          EXIT
       ENDIF
    ENDDO
    nbeta=0
    DO k=1,14
       b2=-1d0
       DO i=0,INT(bcut/hb)
          IF(Phad(DBLE(i)*hb).GE.lev(k))THEN
             b2=DBLE(i)*hb
             EXIT
          ENDIF
       ENDDO
       IF(b2.LE.0d0)CYCLE
       b1=MAX(0d0,b2-hb)
       DO j=1,60
          bm=0.5d0*(b1+b2)
          IF(Phad(bm).GE.lev(k))THEN
             b2=bm
          ELSE
             b1=bm
          ENDIF
       ENDDO
       nbeta=nbeta+1
       beta(nbeta)=0.5d0*(b1+b2)
    ENDDO
    nbeta=nbeta+1
    beta(nbeta)=bcut
    IF(icoll_r.EQ.1)RETURN   ! pp: P in closed form, no table
    hb=amin/64d0
    np_t=INT(bcut/hb)+9
    IF(ALLOCATED(lnchi))DEALLOCATE(lnchi)
    ALLOCATE(lnchi(np_t))
    DO i=1,np_t
       lnchi(i)=Phad_engine(DBLE(i-1)*hb)
    ENDDO
    ! discontinuities of P (MC-Glauber fits): a step much larger than both neighbouring steps;
    ! located by bisection, used as breakpoints, and P is taken from the engine next to them
    njump=0
    DO i=2,np_t-2
       b1=DABS(lnchi(i+1)-lnchi(i))
       IF(b1.GT.1d-6.AND.b1.GT.20d0*MAX(DABS(lnchi(i)-lnchi(i-1)),DABS(lnchi(i+2)-lnchi(i+1))))THEN
          b1=DBLE(i-1)*hb
          b2=DBLE(i)*hb
          bm=Phad_engine(b1)
          DO j=1,80
             IF(DABS(Phad_engine(0.5d0*(b1+b2))-bm).LT.0.5d0*DABS(lnchi(i+1)-lnchi(i)))THEN
                b1=0.5d0*(b1+b2)
             ELSE
                b2=0.5d0*(b1+b2)
             ENDIF
          ENDDO
          IF(njump.LT.10)THEN
             njump=njump+1
             bjump(njump)=b2
             nbeta=nbeta+1
             beta(nbeta)=b2
             WRITE(*,*)"REF: P jump at b [fm] =",b2*GeVm12fm,", from",lnchi(i)," to",lnchi(i+1)
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE build_P

  FUNCTION Phad_engine(b)
    ! hadronic survival probability: the engine's table (validated separately, row A_tab)
    USE ElasticPhotonPhotonFlux,ONLY:FFE_Test_Value
    REAL(KIND(1d0))::Phad_engine
    REAL(KIND(1d0)),INTENT(IN)::b
    Phad_engine=FFE_Test_Value(2,icoll_r,1d-3,b*GeVm12fm)
  END FUNCTION Phad_engine

  FUNCTION Phad(b)
    ! fine local table of Phad_engine (spacing a/64, 8-point Lagrange), P = 1 beyond bcut
    REAL(KIND(1d0))::Phad
    REAL(KIND(1d0)),INTENT(IN)::b
    INTEGER::i
    REAL(KIND(1d0))::u
    IF(icoll_r.EQ.1)THEN
       ! pp: (1-exp(-u))^2, u = b^2/(2 b0), exact (series for small u)
       u=0.5d0*b*b/b0pp
       IF(u.LT.1d-3)THEN
          u=u*(1d0-u*(0.5d0-u*(1d0/6d0-u/24d0)))
       ELSE
          u=1d0-DEXP(-u)
       ENDIF
       Phad=u*u
       RETURN
    ENDIF
    IF(b.GE.bcut)THEN
       Phad=1d0
       RETURN
    ENDIF
    IF(np_t.LE.0)THEN
       Phad=Phad_engine(b)
       RETURN
    ENDIF
    DO i=1,njump
       IF(DABS(b-bjump(i)).LT.6d0*hb)THEN
          Phad=Phad_engine(b)
          RETURN
       ENDIF
    ENDDO
    Phad=MIN(1d0,MAX(0d0,lag8(b/hb,np_t,lnchi)))
  END FUNCTION Phad

  FUNCTION Ptot(b)
    REAL(KIND(1d0))::Ptot
    REAL(KIND(1d0)),INTENT(IN)::b
    REAL(KIND(1d0))::u
    Ptot=Phad(b)
    IF(tag)THEN
       u=MAX(0d0,MIN((DLOG(b)-Tv0)/Th,DBLE(nT-1)))
       Ptot=Ptot*lag8(u,nT,Ptab)
    ENDIF
  END FUNCTION Ptot

  ! ------------------------------------------------------------ ChFF amplitude --
  FUNCTION FF(Q,ib)
    ! Woods-Saxon form factor, F(0)=1: Taylor series for Q*Rmax < 1, radial integral otherwise
    REAL(KIND(1d0))::FF
    REAL(KIND(1d0)),INTENT(IN)::Q
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0))::res,err,c,pts(1)
    INTEGER::k,ier
    IF(Q*(RA(ib)+10d0*aa(ib)).LT.1d0)THEN
       FF=0d0
       c=1d0
       DO k=0,12
          FF=FF+c*rmom(k,ib)
          c=-c*Q*Q
       ENDDO
    ELSEIF(Q.GE.hq(ib)*DBLE(nq(ib)-1))THEN
       FF=0d0
    ELSE
       FF=lag8(Q/hq(ib),nq(ib),fq(:,ib))
    ENDIF
  END FUNCTION FF

  SUBROUTINE build_FF(ib)
    ! F(Q) on a uniform Q grid by adaptive radial integration (interpolated in FF)
    INTEGER,INTENT(IN)::ib
    INTEGER::i,ier
    REAL(KIND(1d0))::res,err,pts(1),Q
    hq(ib)=0.05d0/(RA(ib)+20d0*aa(ib))
    nq(ib)=MIN(40001,INT(55d0/(pi*aa(ib))/hq(ib))+2)
    DO i=1,nq(ib)
       Q=DBLE(i-1)*hq(ib)
       IF(i.EQ.1)THEN
          fq(i,ib)=1d0
          CYCLE
       ENDIF
       ib_c=ib
       q_c=Q
       pts(1)=RA(ib)
       CALL gk21_adapt(f_ffr,0d0,RA(ib)+70d0*aa(ib),1d-17,1d-13,res,err,ier,1,pts)
       fq(i,ib)=res
    ENDDO
  END SUBROUTINE build_FF

  FUNCTION f_ffr(r)
    REAL(KIND(1d0))::f_ffr
    REAL(KIND(1d0)),INTENT(IN)::r
    f_ffr=4d0*pi*r*rho(r,ib_c)*DSIN(q_c*r)/q_c
  END FUNCTION f_ffr

  FUNCTION f_kperp(k)
    ! k^2 J1(kb) F(sqrt(k^2+q^2))/(k^2+q^2)
    REAL(KIND(1d0))::f_kperp
    REAL(KIND(1d0)),INTENT(IN)::k
    REAL(KIND(1d0))::Q2
    Q2=k*k+x_c1*x_c1
    f_kperp=k*k*BESSEL_J1(k*b_c)*FF(DSQRT(Q2),ib_cur)/Q2
  END FUNCTION f_kperp

  FUNCTION f_kpar(k)
    ! q k J0(kb) F/(k^2+q^2): longitudinal amplitude (= q K0(qb) for a point charge)
    REAL(KIND(1d0))::f_kpar
    REAL(KIND(1d0)),INTENT(IN)::k
    REAL(KIND(1d0))::Q2
    Q2=k*k+x_c1*x_c1
    f_kpar=x_c1*k*BESSEL_J0(k*b_c)*FF(DSQRT(Q2),ib_cur)/Q2
  END FUNCTION f_kpar

  SUBROUTINE build_amp(ib,x)
    ! A_perp, A_par of beam ib at x on a fine ln b grid, from the k-space integrals
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::i,ier,nk,j
    REAL(KIND(1d0))::b,res,err,kmax,pts(200),q,smax
    REAL(KIND(1d0)),EXTERNAL::BESSK0,BESSK1
    q=x*mass(ib)
    xcur(ib)=x
    smax=RA(ib)+40d0*aa(ib)
    va0(ib)=DLOG(bmin(ib))-0.2d0
    vah(ib)=1d0/80d0
    nva(ib)=INT((DLOG(smax)+0.2d0-va0(ib))/vah(ib))+2
    kmax=50d0/(pi*aa(ib))
    DO i=1,nva(ib)
       b=DEXP(va0(ib)+DBLE(i-1)*vah(ib))
       ib_cur=ib
       x_c1=q
       b_c=b
       ! breakpoints at the zeros scale of J1(kb) and at k ~ q
       nk=0
       DO j=1,6
          IF(q*10d0**(j-3).LT.kmax)THEN
             nk=nk+1
             pts(nk)=q*10d0**(j-3)
          ENDIF
       ENDDO
       DO j=1,190
          IF(DBLE(j)*pi/b.GE.kmax.OR.nk.GE.200)EXIT
          nk=nk+1
          pts(nk)=DBLE(j)*pi/b
       ENDDO
       CALL gk21_adapt(f_kperp,0d0,kmax,1d-16*q,1d-12,res,err,ier,nk,pts)
       atab(i,1,ib)=res
       CALL gk21_adapt(f_kpar,0d0,kmax,1d-16*q,1d-12,res,err,ier,nk,pts)
       atab(i,2,ib)=res
    ENDDO
    IF(debug)WRITE(*,*)"DBG amp ib,nva,nan count,first,last",ib,nva(ib),COUNT(atab(1:nva(ib),:,ib).NE.atab(1:nva(ib),:,ib)),&
         atab(1,1,ib),atab(nva(ib),1,ib)
  END SUBROUTINE build_amp

  FUNCTION Wfun(ib,x,b)
    ! W = b^2 N(x,b) for unit charge and alpha = 1
    REAL(KIND(1d0))::Wfun
    INTEGER,INTENT(IN)::ib
    REAL(KIND(1d0)),INTENT(IN)::x,b
    REAL(KIND(1d0)),EXTERNAL::BESSK0,BESSK1
    REAL(KIND(1d0))::q,z,ap,al,u,qp,a2
    q=x*mass(ib)
    z=q*b
    IF(z.GT.600d0)THEN
       Wfun=0d0
       RETURN
    ENDIF
    IF(.NOT.chff)THEN
       Wfun=z*z*(BESSK1(z)**2+BESSK0(z)**2/gam(ib)**2)/pi**2
    ELSEIF(isp(ib))THEN
       a2=1d0/0.71d0
       qp=DSQRT(q*q+1d0/a2)
       ap=q*BESSK1(z)-qp*BESSK1(qp*b)-b/(2d0*a2)*BESSK0(qp*b)
       al=q*(BESSK0(z)-BESSK0(qp*b)-b/(2d0*a2*qp)*BESSK1(qp*b))
       Wfun=b*b*(ap*ap+al*al/gam(ib)**2)/pi**2
    ELSEIF(b.GE.RA(ib)+40d0*aa(ib))THEN
       Wfun=z*z*(BESSK1(z)**2+BESSK0(z)**2/gam(ib)**2)/pi**2
    ELSE
       u=(DLOG(b)-va0(ib))/vah(ib)
       ap=lag8(u,nva(ib),atab(1:nva(ib),1,ib))
       al=lag8(u,nva(ib),atab(1:nva(ib),2,ib))
       Wfun=b*b*(ap*ap+al*al/gam(ib)**2)/pi**2
    ENDIF
  END FUNCTION Wfun

  ! ------------------------------------------------------------- the flux --
  FUNCTION f_phi(phi)
    REAL(KIND(1d0))::f_phi
    REAL(KIND(1d0)),INTENT(IN)::phi
    REAL(KIND(1d0))::sn
    REAL(KIND(1d0))::b12
    sn=DSIN(0.5d0*phi)
    b12=DSQRT(d_c*d_c+G4_c*sn*sn)
    f_phi=Ptot(b12)
    IF(wdb)f_phi=f_phi*b12*GeVm12fm
    neval_tot=neval_tot+1
  END FUNCTION f_phi

  FUNCTION Phi_k(b1,b2)
    REAL(KIND(1d0))::Phi_k
    REAL(KIND(1d0)),INTENT(IN)::b1,b2
    REAL(KIND(1d0))::res,err,pts(200),arg,ad
    INTEGER::ier,np,i
    d_c=b1-b2
    ad=DABS(d_c)
    G4_c=4d0*b1*b2
    IF(.NOT.tag.AND..NOT.wdb.AND.ad.GE.bcut)THEN
       Phi_k=2d0*pi
       RETURN
    ENDIF
    np=0
    DO i=1,nbeta
       arg=(beta(i)-ad)*(beta(i)+ad)/G4_c
       IF(beta(i).LE.ad.OR.arg.GE.1d0)CYCLE
       np=np+1
       pts(np)=2d0*DASIN(DSQRT(arg))
    ENDDO
    CALL gk21_adapt(f_phi,0d0,pi,1d-15,1d-12,res,err,ier,np,pts)
    nphi=nphi+1
    IF(ier.NE.0)nphi_ier=nphi_ier+1
    Phi_k=2d0*res
  END FUNCTION Phi_k

  FUNCTION f_u2(v2)
    REAL(KIND(1d0))::f_u2
    REAL(KIND(1d0)),INTENT(IN)::v2
    REAL(KIND(1d0))::b2,w2
    b2=DEXP(v2)
    nf2=nf2+1
    w2=Wfun(2,x_c2,b2)
    IF(w2.EQ.0d0)THEN
       f_u2=0d0
       RETURN
    ENDIF
    f_u2=w2*Phi_k(b1_c,b2)
  END FUNCTION f_u2

  FUNCTION f_w2n(v2)
    REAL(KIND(1d0))::f_w2n
    REAL(KIND(1d0)),INTENT(IN)::v2
    f_w2n=Wfun(2,x_c2,DEXP(v2))
  END FUNCTION f_w2n

  FUNCTION f_w2b(v2)
    REAL(KIND(1d0))::f_w2b
    REAL(KIND(1d0)),INTENT(IN)::v2
    f_w2b=Wfun(2,x_c2,DEXP(v2))*DEXP(v2)
  END FUNCTION f_w2b

  FUNCTION f_u1(v1)
    REAL(KIND(1d0))::f_u1
    REAL(KIND(1d0)),INTENT(IN)::v1
    REAL(KIND(1d0))::w1,res,err,pts(500),lo,hi,bb,eabs
    INTEGER::ier,np,i
    b1_c=DEXP(v1)
    nu1=nu1+1
    w1=Wfun(1,x_c1,b1_c)
    IF(w1.EQ.0d0)THEN
       f_u1=0d0
       RETURN
    ENDIF
    lo=DLOG(bmin(2))
    hi=DLOG(45d0/(x_c2*mass(2)))
    np=0
    DO i=1,nbeta
       bb=b1_c+beta(i)
       np=np+1
       pts(np)=DLOG(bb)
       bb=b1_c-beta(i)
       IF(bb.GT.0d0)THEN
          np=np+1
          pts(np)=DLOG(bb)
       ENDIF
       bb=beta(i)-b1_c
       IF(bb.GT.0d0)THEN
          np=np+1
          pts(np)=DLOG(bb)
       ENDIF
    ENDDO
    np=np+1
    pts(np)=v1
    IF(chff)THEN
       np=np+1
       pts(np)=DLOG(RA(2))
    ENDIF
    ! absolute tolerance: 1e-12 of the unsuppressed scale (Phi <= 2 pi, or 2 pi b12 for <Delta b>);
    ! without it the integral stalls at the noise floor where Phi ~ 0 (b1, b2 inside the overlap)
    IF(wdb)THEN
       eabs=1d-12*2d0*pi*GeVm12fm*(b1_c*wn2+wn2b)
    ELSE
       eabs=1d-12*2d0*pi*wn2
    ENDIF
    CALL gk21_adapt(f_u2,lo,hi,eabs,1d-10,res,err,ier,np,pts)
    nmid=nmid+1
    IF(ier.NE.0)nmid_ier=nmid_ier+1
    IF(debug)WRITE(*,'(A,ES12.4,A,I9,A,I9,A,I9,A,I2,A,ES12.4,A,ES12.4)')" DBG b1=",b1_c," nf2=",nf2," nphi=",nphi,&
         " nP=",neval_tot," ier=",ier," res=",res," err=",err
    f_u1=w1*res
  END FUNCTION f_u1
END MODULE ref_funcs

PROGRAM ref_ab
  USE ref_state
  USE ref_funcs
  USE OpticalGlauber_Geometry
  USE ElasticPhotonPhotonFlux
  IMPLICIT NONE
  include 'run90.inc'
  CHARACTER(len=200)::label='ref',points_file='',out_file='ref.json',code_label='ref'
  INTEGER::icoll=3,A1=208,Z1=82,A2=208,Z2=82,tag1=-2,tag2=-2,iepsilon=0
  REAL(KIND(1d0))::E1=2510d0,E2=2510d0
  LOGICAL::conj_sum=.FALSE.,mc_glauber=.FALSE.
  INTEGER::eval_kinds=8,nW=0,nlookup=0,nWdb=0,algorithm=-1
  LOGICAL::do_grid=.FALSE.,do_eval=.TRUE.,do_pnohad1=.FALSE.,deltab_grid=.FALSE.
  REAL(KIND(1d0))::Wmin=1d0,Wmax=1d0,Wdbmin=1d0,Wdbmax=1d0
  NAMELIST /config/ label,code_label,points_file,out_file,icoll,A1,Z1,A2,Z2,E1,E2,chff,&
       conj_sum,mc_glauber,tag1,tag2,iepsilon,do_grid,do_eval,do_pnohad1,eval_kinds,&
       nlookup,nW,Wmin,Wmax,nWdb,Wdbmin,Wdbmax,algorithm,deltab_grid,wdb,reltol_outer,debug
  CHARACTER(len=400)::arg
  CHARACTER(len=7)::aname
  INTEGER::iu,ios,np,i,kind_i,ib,ier,j,k
  REAL(KIND(1d0)),DIMENSION(2000)::px1,px2,val,er
  INTEGER,DIMENSION(2000)::pk
  REAL(KIND(1d0))::xr1,xr2,anum,znum,rr,adif,wpar,cme,res,err,lo,hi,t0,t1,alpha,pts(10)
  REAL(KIND(1d0)),DIMENSION(2)::xlast
  LOGICAL::swap
  CALL GET_COMMAND_ARGUMENT(1,arg)
  OPEN(NEWUNIT=iu,FILE=TRIM(arg),STATUS='OLD')
  READ(iu,NML=config)
  CLOSE(iu)
  icoll_r=icoll
  ebeam(1)=E1
  ebeam(2)=E2
  nuclearA_beam1=MERGE(0,A1,A1.EQ.1.AND.Z1.EQ.1)
  nuclearZ_beam1=MERGE(0,Z1,A1.EQ.1.AND.Z1.EQ.1)
  nuclearA_beam2=MERGE(0,A2,A2.EQ.1.AND.Z2.EQ.1)
  nuclearZ_beam2=MERGE(0,Z2,A2.EQ.1.AND.Z2.EQ.1)
  neutron_tagging(1)=tag1
  neutron_tagging(2)=tag2
  neutron_conj_sum=conj_sum
  USE_CHARGEFORMFACTOR4PHOTON=chff
  alpha=0.0072992700729927005d0
  use_MC_Glauber=mc_glauber
  IEPSILON_EDFF=iepsilon
  ! beams (engine convention for pA: beam 1 is the proton)
  swap=icoll.EQ.2.AND.nuclearA_beam1.NE.0
  DO ib=1,2
     IF(icoll.EQ.1)THEN
        ! pp: both beams are protons
        isp(ib)=.TRUE.
        AA_(ib)=1d0
        ZZ(ib)=1d0
        RA(ib)=0.877d0/GeVm12fm
        mass(ib)=mproton
        gam(ib)=ebeam(ib)/mproton
        bmin(ib)=RA(ib)
        IF(chff)bmin(ib)=1d-3*bmin(ib)
        CYCLE
     ENDIF
     IF(icoll.EQ.2.AND.ib.EQ.1)THEN
        isp(1)=.TRUE.
        AA_(1)=1d0
        ZZ(1)=1d0
        RA(1)=0.877d0/GeVm12fm
        mass(1)=mproton
        gam(1)=MERGE(E2,E1,swap)/mproton
        bmin(1)=RA(1)
        IF(chff)bmin(1)=1d-3*bmin(1)
        CYCLE
     ENDIF
     IF(icoll.EQ.2)THEN
        IF(swap)THEN
           aname=GetASymbol(A1,Z1)
        ELSE
           aname=GetASymbol(A2,Z2)
        ENDIF
     ELSEIF(ib.EQ.1)THEN
        aname=GetASymbol(A1,Z1)
     ELSE
        aname=GetASymbol(A2,Z2)
     ENDIF
     CALL GetNuclearInfo(aname,anum,znum,rr,adif,wpar)
     isp(ib)=.FALSE.
     AA_(ib)=anum
     ZZ(ib)=znum
     RA(ib)=rr/GeVm12fm
     aa(ib)=adif/GeVm12fm
     ww(ib)=wpar
     mass(ib)=mN
     IF(icoll.EQ.2)THEN
        gam(ib)=MERGE(E1,E2,swap)/mN
     ELSE
        gam(ib)=ebeam(ib)/mN
     ENDIF
     bmin(ib)=RA(ib)
     IF(iepsilon.EQ.1)bmin(ib)=RA(ib)+aa(ib)
     IF(chff)bmin(ib)=1d-3*bmin(ib)
     ! normalisation and radial moments (adaptive)
     rho0(ib)=1d0
     ib_c=ib
     pts(1)=RA(ib)
     CALL gk21_adapt(f_norm,0d0,RA(ib)+70d0*aa(ib),0d0,1d-13,res,err,ier,1,pts)
     IF(ier.NE.0.OR.res.LE.0d0)STOP 'REF: density normalisation failed'
     rho0(ib)=1d0/res
     DO k=0,12
        r_c=DBLE(k)
        CALL gk21_adapt(f_mom,0d0,RA(ib)+70d0*aa(ib),0d0,1d-13,res,err,ier,1,pts)
        IF(ier.EQ.6)STOP 'REF: invalid tolerance'
        rmom(k,ib)=res
        DO j=1,2*k+1
           rmom(k,ib)=rmom(k,ib)/DBLE(j)
        ENDDO
     ENDDO
  ENDDO
  IF(icoll.EQ.2)aa(1)=aa(2)
  IF(icoll.EQ.1)THEN
     ! pp: closed-form survival probability; its length scale sets the breakpoint search step
     b0pp=9.81d0+0.211d0*DLOG(4d0*E1*E2)+0.0185d0*DLOG(4d0*E1*E2)**2
     aa(1)=0.5d0*DSQRT(b0pp)
     aa(2)=aa(1)
  ENDIF
  IF(chff)THEN
     DO ib=1,2
        IF(.NOT.isp(ib))CALL build_FF(ib)
     ENDDO
  ENDIF
  cme=2d0*DSQRT(E1*E2)
  sigNN=sigma_inelastic(cme)*0.1d0/GeVm12fm**2
  CALL CPU_TIME(t0)
  CALL build_P()
  CALL CPU_TIME(t1)
  WRITE(*,*)"REF: P table built in [s]",t1-t0,", bcut [fm] =",bcut*GeVm12fm
  ! tagging: the engine's P_tag (validated separately: m_k against adaptive integration, row A_pbreak)
  tag=tag1.NE.-2.OR.tag2.NE.-2
  IF(tag)THEN
     Th=0.01d0
     Tv0=DLOG(0.5d0/GeVm12fm)
     nT=INT((DLOG(1d12)-Tv0)/Th)+2
     ALLOCATE(Ptab(nT))
     DO i=1,nT
        Ptab(i)=FFE_Test_Value(5,icoll,1d-3,DEXP(Tv0+DBLE(i-1)*Th))
     ENDDO
     ! the tagging factor is long-range: geometric breakpoints bcut 2^k beyond bcut (u2 and phi)
     DO k=1,60
        IF(nbeta.GE.200.OR.bcut*2d0**k.GT.1d12)EXIT
        nbeta=nbeta+1
        beta(nbeta)=bcut*2d0**k
     ENDDO
  ENDIF
  ! points
  np=0
  OPEN(NEWUNIT=iu,FILE=TRIM(points_file),STATUS='OLD')
  DO
     READ(iu,*,IOSTAT=ios)xr1,xr2,kind_i
     IF(ios.NE.0)EXIT
     IF(IAND(eval_kinds,ISHFT(1,kind_i)).EQ.0)CYCLE
     np=np+1
     px1(np)=xr1
     px2(np)=xr2
     pk(np)=kind_i
  ENDDO
  CLOSE(iu)
  ALLOCATE(atab(2000,2,2))
  xlast=-1d0
  DO i=1,np
     CALL CPU_TIME(t0)
     ! engine beam order: (x of beam 1, x of beam 2); for A+p beam 1 is the proton = user beam 2
     IF(swap)THEN
        xr1=px2(i)
        xr2=px1(i)
     ELSE
        xr1=px1(i)
        xr2=px2(i)
     ENDIF
     IF(chff)THEN
        IF(.NOT.isp(1).AND.xr1.NE.xlast(1))CALL build_amp(1,xr1)
        IF(.NOT.isp(2).AND.xr2.NE.xlast(2))CALL build_amp(2,xr2)
        xlast=(/xr1,xr2/)
     ENDIF
     x_c1=xr1
     x_c2=xr2
     pts(1)=DLOG(RA(2))
     CALL gk21_adapt(f_w2n,DLOG(bmin(2)),DLOG(45d0/(xr2*mass(2))),0d0,1d-10,wn2,err,ier,1,pts)
     CALL gk21_adapt(f_w2b,DLOG(bmin(2)),DLOG(45d0/(xr2*mass(2))),0d0,1d-10,wn2b,err,ier,1,pts)
     IF(debug)WRITE(*,*)"DBG wn2,wn2b",wn2,wn2b
     lo=DLOG(bmin(1))
     hi=DLOG(45d0/(xr1*mass(1)))
     pts(1)=DLOG(RA(1))
     CALL gk21_adapt(f_u1,lo,hi,0d0,reltol_outer,res,err,ier,1,pts)
     val(i)=2d0*pi*alpha**2*ZZ(1)**2*ZZ(2)**2*res/(px1(i)*px2(i))
     er(i)=2d0*pi*alpha**2*ZZ(1)**2*ZZ(2)**2*err/(px1(i)*px2(i))
     CALL CPU_TIME(t1)
     WRITE(*,'(A,I8,A,I6,A,I10,A,I8,A,I10,A,I3)')" REF: diag nmid=",nmid," mid_ier=",nmid_ier," nphi=",nphi,&
          " phi_ier=",nphi_ier," nP=",neval_tot," outer ier=",ier
     nmid=0
     nmid_ier=0
     nphi=0
     nphi_ier=0
     neval_tot=0
     WRITE(*,'(A,2ES12.4,A,ES24.16,A,ES9.2,A,I2,A,F8.2)')" REF: x1,x2=",px1(i),px2(i)," flux=",val(i),&
          " relerr=",er(i)/val(i)," ier=",ier," t[s]=",t1-t0
  ENDDO
  OPEN(NEWUNIT=iu,FILE=TRIM(out_file),STATUS='REPLACE')
  WRITE(iu,'(A)')'{"schema": "gupc-ref-1", "label": "'//TRIM(label)//'", "points": ['
  DO i=1,np
     WRITE(iu,'(A,ES26.17E3,A,ES26.17E3,A,I0,A,ES26.17E3,A,ES26.17E3,A)')'  {"x1": ',px1(i),', "x2": ',px2(i),&
          ', "kind": ',pk(i),', "value": ',val(i),', "error": ',er(i),'}'//TRIM(MERGE(',',' ',i.LT.np))
  ENDDO
  WRITE(iu,'(A)')']}'
  CLOSE(iu)
END PROGRAM ref_ab
