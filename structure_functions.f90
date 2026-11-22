MODULE structure_functions
  ! Proton structure functions F2(x,Q2), FL(x,Q2) (x = x_Bjorken) of the inelastic photon distribution, as the Python
  ! code of the inelastic photon (photon_pdf_components: Low_Q2.py, High_Q2.py, F2_parameterization, F2_FL_h of
  ! Inelastic.py), in the conventions of LUXqed17 (arXiv:1708.01256, sec. 8):
  !  low Q2 / resonances: CLAS data (spline) joined to the fits Christy-Bosted, SLAC-MIT, Breidenbach-Kuti;
  !  continuum: GD11-P (HERMES), and DGLAP at high Q2 (PDF set of LHAPDF, NLO/NNLO coefficient functions);
  !  FL = (1+4 mp^2 x^2/Q2) R/(1+R) F2 with R1998 for the fits.
  ! Smooth joins everywhere (W2 step, fits -> DGLAP ramp in ln Q2, R below Q2=0.34 with value and slope matched).
  ! low_Q2_model selects only the fits outside the CLAS data: 1 Christy-Bosted, SLAC-MIT, Breidenbach-Kuti with the
  ! lowest chi2 (Python 'spline best smooth', default); 2 SLAC-MIT -> Breidenbach-Kuti, the fits of LUXqed17
  ! (Python 'spline smooth', the 'luxqed' fits with the smooth joins).
  USE interpolation,ONLY:lagrange_interp_2d
  USE ElasticFlux,ONLY:mproton,GeVm12fm ! proton mass and GeV-1 -> fm of gammaUPC
  IMPLICIT NONE
  PRIVATE
  PUBLIC::F2_FL_h,F2_FL_DGLAP,sf_init,sf_check
  PUBLIC::mpi,W2lo,W2hi,Q2_pdf,Q2_pdf_ramp,rho_W2 ! also used by InelasticFlux
  INTEGER,PUBLIC::low_Q2_model=1
  ! PDF set of the DGLAP part by its LHAPDF ID (set ID + member, as in pdfsets.index); sf_init derives the set name
  ! (sf_pdf_set) and the member
  INTEGER,PUBLIC::sf_pdf_id=91700 ! PDF4LHC15_nnlo_100, member 0
  CHARACTER(LEN=128),PUBLIC::sf_pdf_set=''
  ! continuation of the PDFs at small x: below x_c = XMin of the set (all Q2), and below x_c = sf_x_safe for
  ! Q2 < sf_Q2_safe, x f(x,Q2) = x f(x_c,Q2) (x/x_c)^(-lambda), lambda = ln(x f(x_c)/x f(2 x_c))/ln 2 (x f(x_c) if either
  ! is not positive); a warning is printed once. sf_x_safe < 0 (default): automatic, from a scan of the set at its
  ! initialisation (pdf_scan: non-positive values or kinks of ln x f in ln x of g, u+ubar, d+dbar, s+sbar for
  ! XMin <= x <= 1e-3, 6 <= Q2 <= 50 GeV^2, the region of the DGLAP part at small x); e.g. PDF4LHC15_nnlo_100 has
  ! corrupted gluon knots at x < 1e-6, Q = 1.39-2.25 GeV, spread by the LHAPDF interpolation up to Q2 ~ 8 GeV^2.
  ! sf_Q2_safe < 0: all Q2 for a manual sf_x_safe
  REAL(KIND(1d0)),PUBLIC::sf_x_safe=-1d0,sf_Q2_safe=-1d0
  INTEGER,PUBLIC::sf_version=0 ! incremented by each (re)initialisation of the PDF set (for the caches of InelasticFlux)
  ! diagnostics of the continuation (reset by sf_init): number of continued PDF lookups and the smallest x requested;
  ! the continuation affects photons with x below about x_c (sf_xc, sf_Q2c: the thresholds in use)
  INTEGER(KIND=8),PUBLIC::sf_continued=0
  REAL(KIND(1d0)),PUBLIC::sf_continued_xmin=1d0,sf_xc=0d0,sf_Q2c=0d0
  ! constants: proton mass and fm -> GeV-1 of gammaUPC (ElasticFlux); pion mass, alpha and unit conversion of the fits
  REAL(KIND(1d0)),PARAMETER::mp=mproton,mp2=mp**2,mpi=0.135d0,alpha_QED=1d0/137.036d0
  REAL(KIND(1d0)),PARAMETER::fm2GeV=1d0/GeVm12fm,barn2GeV=2.566d-3,mdip2=0.71d0,mup=2.793d0
  REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
  ! transitions (Python global_variables): resonances -> continuum W2lo -> W2hi; fits -> DGLAP across Q2_pdf (ramp of
  ! ln Q2 by a factor Q2_pdf_ramp each side)
  REAL(KIND(1d0)),PARAMETER::W2lo=3d0,W2hi=4d0,Q2_pdf=9d0,Q2_pdf_ramp=1.5d0
  INCLUDE 'clas_spline.inc'
  ! LHAPDF set (DGLAP part)
  LOGICAL::pdf_init=.FALSE.,pdf_warned=.FALSE.
  INTEGER::sf_order=1,pdf_member=0,init_id=-1
  REAL(KIND(1d0))::mcharm=1.5d0,mbottom=4.75d0,pdf_xmin=0d0,x_c=0d0,Q2_c=0d0,init_xs=-2d0,init_qs=-2d0
CONTAINS
  SUBROUTINE F2_FL_h(x,Q2,F2,FL)
    ! F2_FL_h of Inelastic.py (smooth models): F2 = (1-r) F2lo + (r-s) F2_GD11-P + s F2_DGLAP, r = rho(W2)
    ! (W2lo -> W2hi), s = r ramp(ln Q2) across Q2_pdf
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),INTENT(OUT)::F2,FL
    REAL(KIND(1d0))::W2,r,s,F,R1,F2d,FLd
    W2=Q2*(1d0-x)/x+mp2
    r=MERGE(0d0,MERGE(1d0,rho_W2(MIN(MAX(W2,W2lo),W2hi)),W2.GE.W2hi),W2.LE.W2lo)
    s=r*ramp(DLOG(Q2),DLOG(Q2_pdf/Q2_pdf_ramp),DLOG(Q2_pdf*Q2_pdf_ramp))
    F2=0d0; FL=0d0
    IF(s.GT.0d0)THEN
       CALL F2_FL_DGLAP(x,Q2,F2d,FLd)
       F2=s*F2d; FL=s*FLd
    ENDIF
    IF(s.LT.1d0)THEN
       F=0d0
       IF(r.LT.1d0)F=(1d0-r)*F2lo(x,Q2)
       IF(r.GT.s)F=F+(r-s)*F2_gd(x,Q2)
       R1=Rext(x,Q2)
       F2=F2+F
       FL=FL+(1d0+4d0*mp2*x**2/Q2)*R1/(1d0+R1)*F
    ENDIF
    RETURN
  END SUBROUTINE F2_FL_h

  ! ---------------------------------------------------------------- steps, W2 transition
  FUNCTION step(w) ! 2w^2-w^4, w clipped to [0,1]
    REAL(KIND(1d0))::step
    REAL(KIND(1d0)),INTENT(IN)::w
    REAL(KIND(1d0))::v
    v=MIN(MAX(w,0d0),1d0)
    step=2d0*v**2-v**4
  END FUNCTION step
  FUNCTION ramp(v,a,b)
    REAL(KIND(1d0))::ramp
    REAL(KIND(1d0)),INTENT(IN)::v,a,b
    ramp=step((v-a)/(b-a))
  END FUNCTION ramp
  FUNCTION rho_W2(W2) ! rho of Low_Q2.py (W2lo -> W2hi), not clipped
    REAL(KIND(1d0))::rho_W2
    REAL(KIND(1d0)),INTENT(IN)::W2
    REAL(KIND(1d0))::o
    o=(W2-W2lo)/(W2hi-W2lo)
    rho_W2=2d0*o**2-o**4
  END FUNCTION rho_W2

  ! ---------------------------------------------------------------- low-Q2 models
  FUNCTION F2lo(x,Q2) ! smooth(spline_weight, F2_spline, outside): F2_best_smooth (1) or F2_fits(0.30,0.45) (2)
    REAL(KIND(1d0))::F2lo
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0))::w
    w=spline_weight(x,Q2)
    IF(w.EQ.1d0)THEN
       F2lo=F2_spline(x,Q2)
       RETURN
    ENDIF
    IF(low_Q2_model.EQ.1)THEN
       F2lo=F2_best_smooth(x,Q2)
    ELSE
       F2lo=F2_fits(x,Q2,0.30d0,0.45d0)
    ENDIF
    IF(w.GT.0d0)F2lo=w*F2_spline(x,Q2)+(1d0-w)*F2lo
  END FUNCTION F2lo
  FUNCTION spline_weight(x,Q2) ! ramps inside the low-x edge of the rows and at the Q2 edges of the data
    REAL(KIND(1d0))::spline_weight
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),PARAMETER::x_margin=0.5d0,x_ramp=2d0,Q2_ramp=1.15d0
    REAL(KIND(1d0))::lq,q_lo,q_hi,wQ,lo
    lq=DLOG(Q2); q_lo=DLOG(Q2rowclas(1)); q_hi=DLOG(Q2rowclas(NROWCLAS))
    wQ=ramp(lq,q_lo-DLOG(Q2_ramp),q_lo)*ramp(-lq,-q_hi-DLOG(Q2_ramp),-q_hi)
    spline_weight=0d0
    IF(wQ.EQ.0d0)RETURN
    lo=interp1(Q2,Q2rowclas,xminrowclas)
    spline_weight=wQ*ramp(x,lo+x_margin*dxclas,lo+(x_margin+x_ramp)*dxclas)
  END FUNCTION spline_weight
  FUNCTION F2_spline(x,Q2) ! the spline in (W,Q2); below the high-x edge W0 of the rows, cubic to 0 at the pion threshold
    REAL(KIND(1d0))::F2_spline
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0))::W,xe,W0,Wth,f0,d0,w1
    W=DSQRT(Q2*(1d0-x)/x+mp2)
    xe=interp1(Q2,Q2rowclas,xmaxrowclas)
    W0=DSQRT(Q2*(1d0-xe)/xe+mp2)
    IF(W.GE.W0)THEN
       F2_spline=MAX(clas_bispline(W,Q2,0),0d0)
       RETURN
    ENDIF
    Wth=mp+mpi
    f0=MAX(clas_bispline(W0,Q2,0),0d0)
    d0=MIN(clas_bispline(W0,Q2,1)*(W0-Wth),3d0*f0)
    w1=MAX((W-Wth)/(W0-Wth),0d0)
    F2_spline=w1**2*(f0*(3d0-2d0*w1)-d0*(1d0-w1))
  END FUNCTION F2_spline
  FUNCTION interp1(v,xs,ys) ! np.interp: linear, constant outside
    REAL(KIND(1d0))::interp1
    REAL(KIND(1d0)),INTENT(IN)::v
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::xs,ys
    INTEGER::i,n
    n=SIZE(xs)
    IF(v.LE.xs(1))THEN
       interp1=ys(1)
    ELSEIF(v.GE.xs(n))THEN
       interp1=ys(n)
    ELSE
       DO i=1,n-1
          IF(v.LT.xs(i+1))EXIT
       ENDDO
       interp1=ys(i)+(ys(i+1)-ys(i))*(v-xs(i))/(xs(i+1)-xs(i))
    ENDIF
  END FUNCTION interp1
  FUNCTION clas_bispline(W,Q2,dW) ! the CLAS spline (dW=1: d/dW), arguments clamped to the knot range as FITPACK:
    ! cubic x cubic Lagrange interpolation (interpolation.f90) on the 4x4 nodes of the knot cell, exact for the spline
    REAL(KIND(1d0))::clas_bispline
    REAL(KIND(1d0)),INTENT(IN)::W,Q2
    INTEGER,INTENT(IN)::dW
    REAL(KIND(1d0)),DIMENSION(1)::Wc,Qc,z
    REAL(KIND(1d0)),DIMENSION(4)::wn,qn
    REAL(KIND(1d0)),DIMENSION(4,4)::zn
    INTEGER::i,j
    Wc(1)=MIN(MAX(W,kWclas(1)),kWclas(NKWCLAS))
    Qc(1)=MIN(MAX(Q2,kQclas(1)),kQclas(NKQCLAS))
    i=1
    DO WHILE(i.LT.NKWCLAS-1.AND.Wc(1).GE.kWclas(i+1))
       i=i+1
    ENDDO
    j=1
    DO WHILE(j.LT.NKQCLAS-1.AND.Qc(1).GE.kQclas(j+1))
       j=j+1
    ENDDO
    i=3*i-2; j=3*j-2
    wn=WGclas(i:i+3); qn=QGclas(j:j+3)
    IF(dW.EQ.1)THEN
       zn=dSclas(i:i+3,j:j+3)
    ELSE
       zn=Sclas(i:i+3,j:j+3)
    ENDIF
    CALL lagrange_interp_2d(3,3,wn,qn,zn,1,Wc,Qc,z)
    clas_bispline=z(1)
  END FUNCTION clas_bispline

  FUNCTION F2_best_smooth(x,Q2) ! Christy-Bosted, SLAC-MIT rows 1, 2, Christy-Bosted, Breidenbach-Kuti, smooth steps
    REAL(KIND(1d0))::F2_best_smooth
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),PARAMETER::Q2_cb=8d0,Q2_bk_cb=6d0
    REAL(KIND(1d0)),DIMENSION(2,3),PARAMETER::win=RESHAPE((/0.04d0,0.05d0,0.10d0,0.14d0,0.20d0,0.225d0/),(/2,3/))
    REAL(KIND(1d0))::W2,s
    INTEGER::i
    W2=Q2*(1d0-x)/x+mp2
    IF(Q2.GE.Q2_cb)THEN
       F2_best_smooth=F2_bk_pos(W2,Q2)
       RETURN
    ENDIF
    IF(Q2.GE.Q2_bk_cb)THEN
       s=ramp(DLOG(Q2),DLOG(Q2_bk_cb),DLOG(Q2_cb))
       F2_best_smooth=(1d0-s)*F2_cb(W2,Q2)+s*F2_bk_pos(W2,Q2)
       RETURN
    ENDIF
    DO i=1,3
       IF(Q2.LT.win(2,i))THEN
          s=ramp(DLOG(Q2),DLOG(win(1,i)),DLOG(win(2,i)))
          F2_best_smooth=(1d0-s)*piece(i)
          IF(s.GT.0d0)F2_best_smooth=F2_best_smooth+s*piece(i+1)
          RETURN
       ENDIF
    ENDDO
    F2_best_smooth=F2_cb(W2,Q2)
  CONTAINS
    FUNCTION piece(k)
      REAL(KIND(1d0))::piece
      INTEGER,INTENT(IN)::k
      SELECT CASE(k)
      CASE(1,4); piece=F2_cb(W2,Q2)
      CASE(2); piece=F2_sa(W2,Q2,1)
      CASE DEFAULT; piece=F2_sa(W2,Q2,2)
      END SELECT
    END FUNCTION piece
  END FUNCTION F2_best_smooth
  FUNCTION F2_bk_pos(W2,Q2) ! Breidenbach-Kuti, scaled to zero at the pion threshold below W_bk_lo
    REAL(KIND(1d0))::F2_bk_pos
    REAL(KIND(1d0)),INTENT(IN)::W2,Q2
    REAL(KIND(1d0)),PARAMETER::W_bk_lo=1.25d0
    REAL(KIND(1d0))::W
    W=DSQRT(W2)
    IF(W.GE.W_bk_lo)THEN
       F2_bk_pos=F2_bk(W2,Q2)
    ELSE
       F2_bk_pos=F2_bk(W_bk_lo**2,Q2)*step((W-(mp+mpi))/(W_bk_lo-(mp+mpi)))
    ENDIF
  END FUNCTION F2_bk_pos
  FUNCTION F2_fits(x,Q2,sa_lo,sa_hi) ! SLAC-MIT row 2 -> row 4 (sa window) -> Breidenbach-Kuti (Q2 3.5-4.5, W 1.25-1.35)
    REAL(KIND(1d0))::F2_fits
    REAL(KIND(1d0)),INTENT(IN)::x,Q2,sa_lo,sa_hi
    REAL(KIND(1d0))::W2,s1,s2
    W2=Q2*(1d0-x)/x+mp2
    s1=ramp(DLOG(Q2),DLOG(sa_lo),DLOG(sa_hi))
    s2=ramp(DLOG(Q2),DLOG(3.5d0),DLOG(4.5d0))*ramp(DSQRT(W2),1.25d0,1.35d0)
    F2_fits=0d0
    IF(s1.LT.1d0)F2_fits=F2_fits+(1d0-s1)*F2_sa(W2,Q2,2)
    IF(s1.GT.0d0.AND.s2.LT.1d0)F2_fits=F2_fits+s1*(1d0-s2)*F2_sa(W2,Q2,4)
    IF(s1.GT.0d0.AND.s2.GT.0d0)F2_fits=F2_fits+s1*s2*F2_bk(W2,Q2)
  END FUNCTION F2_fits

  ! ---------------------------------------------------------------- R1998 (Rext of Low_Q2.py, x_R at fixed W2)
  FUNCTION Rext(x,Q2)
    REAL(KIND(1d0))::Rext
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),PARAMETER::Qr2=0.34d0
    REAL(KIND(1d0))::u,q1,q2b,R0,d
    IF(Q2.GE.Qr2)THEN
       Rext=Rav(x,Q2)
       RETURN
    ENDIF
    u=Q2/Qr2 ! smooth models: a u + b u^3 with the value and slope of Rav at Qr2
    q1=Qr2*(1d0+1d-6); q2b=Qr2*(1d0-1d-6)
    R0=Rav(x_R(x,Q2,Qr2),Qr2)
    d=(Rav(x_R(x,Q2,q1),q1)-Rav(x_R(x,Q2,q2b),q2b))/2d-6
    Rext=(3d0*R0-d)/2d0*u+(d-R0)/2d0*u**3
  END FUNCTION Rext
  FUNCTION x_R(x,Q2,Q2new) ! x at Q2new with the same W2
    REAL(KIND(1d0))::x_R
    REAL(KIND(1d0)),INTENT(IN)::x,Q2,Q2new
    x_R=Q2new/(Q2new+Q2*(1d0-x)/x)
  END FUNCTION x_R
  FUNCTION Rav(x,Q2)
    REAL(KIND(1d0))::Rav
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::a=(/0.0485d0,0.5470d0,2.0621d0,-0.3804d0,0.5090d0,-0.0285d0/)
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::b=(/0.0481d0,0.6114d0,-0.3509d0,-0.4611d0,0.7172d0,-0.0317d0/)
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::c=(/0.0577d0,0.4644d0,1.8288d0,12.3708d0,-43.1043d0,41.7415d0/)
    REAL(KIND(1d0))::th,Ra,Rb,Rc,L
    th=1d0+12d0*(Q2/(1d0+Q2))*(0.015625d0/(0.015625d0+x**2))
    L=DLOG(25d0*Q2)
    Ra=(a(1)/L)*th+(a(2)/(Q2**4+a(3)**4)**0.25d0)*(1d0+a(4)*x+a(5)*x**2)*(x**a(6))
    Rb=(b(1)/L)*th+((b(2)/Q2)+(b(3)/(Q2**2+0.09d0)))*(1d0+b(4)*x+b(5)*x**2)*(x**b(6))
    Rc=(c(1)/L)*th+c(2)*((Q2-(c(4)*x+c(5)*x**2+c(6)*x**3))**2+c(3)**2)**(-0.5d0)
    Rav=(Ra+Rb+Rc)/3d0
  END FUNCTION Rav

  ! ---------------------------------------------------------------- GD11-P (HERMES, arXiv:1103.5704; ALLM form)
  FUNCTION F2_gd(x,Q2)
    REAL(KIND(1d0))::F2_gd
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),PARAMETER::lQCD2=0.06527d0,Q20=1.374d0,m2=0.5063d0,mreg2=0.03190d0,mpom2=34.75d0
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER::ar=(/0.3425d0,1.0603d0,0.5164d0/),br=(/-10.408d0,14.857d0,0.07739d0/)
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER::cr=(/1.3633d0,2.256d0,2.209d0/)
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER::ap=(/-0.11895d0,-0.4783d0,1.353d0/),bp=(/1.0833d0,2.656d0,1.771d0/)
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER::cp=(/0.3638d0,0.1211d0,1.166d0/)
    REAL(KIND(1d0))::W2,t,xx,xreg,xpom
    W2=Q2*(1d0-x)/x+mp2
    t=DLOG(DLOG((Q20+Q2)/lQCD2)/DLOG(Q20/lQCD2))
    xreg=(Q2+mreg2)/(Q2+W2-(mp2-mreg2))
    xpom=(Q2+mpom2)/(Q2+W2-(mp2-mpom2))
    xx=Q2/(Q2+W2-mp2)
    F2_gd=(Q2/(Q2+m2))*(t1(t,cr)*xreg**t1(t,ar)*(1d0-xx)**t1(t,br)+t2(t,cp)*xpom**t2(t,ap)*(1d0-xx)**t1(t,bp))
  CONTAINS
    FUNCTION t1(tt,p)
      REAL(KIND(1d0))::t1
      REAL(KIND(1d0)),INTENT(IN)::tt,p(3)
      t1=p(1)+p(2)*tt**p(3)
    END FUNCTION t1
    FUNCTION t2(tt,p)
      REAL(KIND(1d0))::t2
      REAL(KIND(1d0)),INTENT(IN)::tt,p(3)
      t2=p(1)+(p(1)-p(2))*((1d0+tt**p(3))**(-1)-1d0)
    END FUNCTION t2
  END FUNCTION F2_gd

  ! ---------------------------------------------------------------- Breidenbach-Kuti (Phys. Lett. 41B (1972))
  FUNCTION F2_bk(W2,Q2)
    REAL(KIND(1d0))::F2_bk
    REAL(KIND(1d0)),INTENT(IN)::W2,Q2
    REAL(KIND(1d0)),DIMENSION(4),PARAMETER::mass=(/1.226d0,1.508d0,1.705d0,1.920d0/)
    REAL(KIND(1d0)),DIMENSION(4),PARAMETER::dw=(/0.115d0,0.08d0,0.085d0,0.22d0/)
    REAL(KIND(1d0)),DIMENSION(4),PARAMETER::a=(/0.7178d0,0.3999d0,0.3268d0,0.096d0/)
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER::c=(/0.8285d0,1.383d0,-2.003d0/),bb=(/-0.5291d0,1.26d0,0.22d0/)
    REAL(KIND(1d0)),PARAMETER::Wt=1.2172d0,R=0.8d0*fm2GeV
    REAL(KIND(1d0))::f2,res,bkg,m2,g
    INTEGER::i
    F2_bk=0d0
    IF(W2.LE.(mp+mpi)**2)RETURN
    f2=0d0
    DO i=1,3
       f2=f2+c(i)*(W2/(W2+Q2))**(i+2)
    ENDDO
    res=0d0
    DO i=1,4
       m2=mass(i)**2
       IF(i.EQ.1.AND.W2.LE.m2)THEN
          g=dw(1)*(V(W2)/V(m2))*DSQRT(k2(W2)/k2(m2))
          res=res+a(1)*(mp/DSQRT(W2))*((W2-mp2)/(m2-mp2))**2*(k2(m2)/k2(W2))**1.5d0&
               *(g**2*m2)/((W2-m2)**2+g**2*m2)
       ELSE
          res=res+a(i)*dw(i)**2*m2/((W2-m2)**2+dw(i)**2*m2)
       ENDIF
    ENDDO
    bkg=0d0
    DO i=1,3
       bkg=bkg+bb(i)/(1d0+DSQRT(W2)-Wt)**i
    ENDDO
    F2_bk=f2*(res+1d0-bkg)
  CONTAINS
    FUNCTION k2(w)
      REAL(KIND(1d0))::k2
      REAL(KIND(1d0)),INTENT(IN)::w
      k2=((w+mp2-mpi**2)**2)/(4d0*w)-mp2
    END FUNCTION k2
    FUNCTION V(w)
      REAL(KIND(1d0))::V
      REAL(KIND(1d0)),INTENT(IN)::w
      V=(k2(w)*R**2)/(1d0+k2(w)*R**2)
    END FUNCTION V
  END FUNCTION F2_bk

  ! ---------------------------------------------------------------- SLAC-MIT (Stein et al., PRD 12 (1975) 1884)
  FUNCTION F2_sa(W2,Q2,row)
    REAL(KIND(1d0))::F2_sa
    REAL(KIND(1d0)),INTENT(IN)::W2,Q2
    INTEGER,INTENT(IN)::row
    REAL(KIND(1d0)),DIMENSION(4),PARAMETER::X=(/0.160d0,0.350d0,0.350d0,0.350d0/)
    INTEGER,DIMENSION(4),PARAMETER::oam=(/1,2,3,3/),tam=(/1,1,2,2/)
    REAL(KIND(1d0)),DIMENSION(4,7),PARAMETER::mass=RESHAPE((/1.227d0,1.519d0,1.675d0,1d-15,1.231d0,1.516d0,1.681d0,1d-15,&
         1.225d0,1.512d0,1.692d0,1.950d0,1.224d0,1.509d0,1.696d0,1.950d0,1.226d0,1.507d0,1.7d0,1.950d0,&
         1.233d0,1.511d0,1.705d0,1.950d0,1.236d0,1.513d0,1.713d0,1.950d0/),(/4,7/))
    REAL(KIND(1d0)),DIMENSION(4,7),PARAMETER::dw=RESHAPE((/0.116d0,0.144d0,0.113d0,0d0,0.123d0,0.13d0,0.102d0,0d0,&
         0.119d0,0.114d0,0.115d0,0.2d0,0.119d0,0.097d0,0.105d0,0.2d0,0.12d0,0.088d0,0.094d0,0.2d0,&
         0.124d0,0.077d0,0.092d0,0.2d0,0.141d0,0.078d0,0.102d0,0.2d0/),(/4,7/))
    REAL(KIND(1d0)),DIMENSION(4,7),PARAMETER::amp=RESHAPE((/289.7d0,118.2d0,64.92d0,0d0,221.3d0,100.8d0,61.56d0,0d0,&
         167.1d0,82.94d0,61.06d0,21.87d0,137.7d0,72.52d0,57.32d0,12.57d0,118.8d0,70.62d0,51.18d0,10.26d0,&
         104.5d0,62.06d0,50.79d0,13.71d0,93.79d0,61.15d0,47d0,11.78d0/),(/4,7/))
    REAL(KIND(1d0)),DIMENSION(5),PARAMETER::a=(/1.0621d0,-2.2594d0,10.54d0,-15.8277d0,6.7931d0/)
    REAL(KIND(1d0)),DIMENSION(4,7),PARAMETER::c=RESHAPE((/351.2d0,-926.5d0,1857.6d0,-1191d0,167.2d0,-364.3d0,1270.2d0,&
         -1017.1d0,186.2d0,-481.4d0,1323.5d0,-1011.9d0,164.6d0,-316.4d0,912.1d0,-681.2d0,185.3d0,-333.4d0,879.5d0,&
         -651.6d0,54d0,584.2d0,-940.3d0,446.1d0,94.6d0,320.3d0,-482.9d0,204.4d0/),(/4,7/))
    REAL(KIND(1d0)),PARAMETER::w2lim=3.7014d0,f2lim=0.3081d0
    REAL(KIND(1d0))::tau,GE,GM,w2el,f2,BW,BWNR,m2,gt,gg,norm
    INTEGER::i
    F2_sa=0d0
    IF(W2.LE.(mp+mpi)**2)RETURN
    tau=Q2/(4d0*mp2); GE=1d0/(1d0+Q2/mdip2)**2; GM=mup*GE
    w2el=(GE**2+tau*GM**2)/(1d0+tau)
    f2=0d0
    DO i=1,5
       f2=f2+a(i)*(W2/(W2+Q2))**(i+2)
    ENDDO
    BW=0d0
    DO i=1,4
       m2=mass(i,row)**2
       gt=dw(i,row)*(k2(W2)/k2(m2))**(oam(i)+0.5d0)*((k2(m2)+X(i)**2)/(k2(W2)+X(i)**2))**oam(i)
       gg=dw(i,row)*(q2s(W2)/q2s(m2))**tam(i)*((q2s(m2)+X(i)**2)/(q2s(W2)+X(i)**2))**tam(i)
       BW=BW+(amp(i,row)*barn2GeV)*(q2s(m2)/q2s(W2))*m2*((gt*gg)/((m2-W2)**2+m2*gt**2))
    ENDDO
    BWNR=0d0
    DO i=1,4
       BWNR=BWNR+(c(i,row)*barn2GeV)*k2(W2)**(i/2d0)
    ENDDO
    norm=(1d0-w2el)*f2/(4d0*pi**2*alpha_QED*f2lim*w2lim)
    F2_sa=norm*(BW+BWNR)
  CONTAINS
    FUNCTION k2(w)
      REAL(KIND(1d0))::k2
      REAL(KIND(1d0)),INTENT(IN)::w
      k2=((w+mp2-mpi**2)**2)/(4d0*w)-mp2
    END FUNCTION k2
    FUNCTION q2s(w)
      REAL(KIND(1d0))::q2s
      REAL(KIND(1d0)),INTENT(IN)::w
      q2s=((w+mp2+Q2)**2)/(4d0*w)-mp2
    END FUNCTION q2s
  END FUNCTION F2_sa

  ! ---------------------------------------------------------------- Christy-Bosted (PRC 81 (2010) 055213)
  FUNCTION F2_cb(W2,Q2)
    REAL(KIND(1d0))::F2_cb
    REAL(KIND(1d0)),INTENT(IN)::W2,Q2
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::M=(/1.230d0,1.530d0,1.506d0,1.698d0,1.665d0,1.433d0,1.934d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::G=(/0.136d0,0.220d0,0.083d0,0.096d0,0.109d0,0.379d0,0.380d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::lr=(/1d0,0d0,2d0,3d0,0d0,1d0,3d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::Xr=(/0.1446d0,0.215d0,0.215d0,0.215d0,0.215d0,0.215d0,0.215d0/)
    REAL(KIND(1d0)),DIMENSION(7,3),PARAMETER::beta=RESHAPE((/1.0d0,0.45d0,0.65d0,0.65d0,0.4d0,0.65d0,0.5d0,&
         0d0,0.10d0,0.35d0,0.35d0,0.5d0,0.35d0,0.5d0,0d0,0.45d0,0d0,0d0,0.1d0,0d0,0d0/),(/7,3/))
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::AT0=(/7.780d0,6.335d0,0.603d0,2.330d0,1.979d0,0.0225d0,3.419d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::ac=(/4.229d0,6823.2d0,21.24d0,-0.288d0,-0.562d0,462.13d0,0d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::bc=(/1.260d0,33521.0d0,0.056d0,0.186d0,0.390d0,0.192d0,0d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::cc=(/2.124d0,2.569d0,2.489d0,0.064d0,0.549d0,1.914d0,1.0d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::AL0=(/29.4140d0,0d0,157.92d0,4.2160d0,13.764d0,5.5124d0,11.0d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::dc=(/19.910d0,0d0,97.046d0,0.038d0,0.314d0,0.054d0,1.895d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::ec=(/0.226d0,0d0,0.310d0,1.218d0,3.000d0,1.309d0,0.514d0/)
    REAL(KIND(1d0)),PARAMETER::meta=0.547d0,ub_GeV=389.4d0,Q20_T=0.05d0,Q20_L=0.125d0,m0_L=4.2802d0
    REAL(KIND(1d0)),DIMENSION(5,2),PARAMETER::nrT=RESHAPE((/246.1d0,0.0675d0,1.3501d0,0.1205d0,-0.0038d0,&
         -89.4d0,0.2098d0,1.5715d0,0.0907d0,0.0104d0/),(/5,2/))
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::nrL=(/86.7d0,0d0,4.0294d0,3.1285d0,0.3340d0,4.9623d0/)
    REAL(KIND(1d0))::W,BW(7),G1,G2,G3,Gtot,Kcm,Kcmr,Ggam,AT(7),AL(7),sT,sL,dW,xp,x,t,F1,FL,r,d
    INTEGER::i
    F2_cb=0d0
    IF(W2.LE.(mp+mpi)**2)RETURN
    W=DSQRT(W2)
    Kcm=(W2-mp2)/(2d0*W)
    DO i=1,7
       CALL ratio(pcm(W,mpi),pcm(M(i),mpi),r,d); G1=G(i)*r**(2d0*lr(i)+1d0)*d**lr(i)
       CALL ratio(pcm(W,2d0*mpi),pcm(M(i),2d0*mpi),r,d); G2=W*G(i)/M(i)*r**(2d0*lr(i)+4d0)*d**(lr(i)+2d0)
       CALL ratio(pcm(W,meta),pcm(M(i),meta),r,d); G3=G(i)*r**(2d0*lr(i)+1d0)*d**lr(i)
       Gtot=beta(i,1)*G1+beta(i,2)*G2+beta(i,3)*G3
       Kcmr=(M(i)**2-mp2)/(2d0*M(i))
       Ggam=G(i)*(Kcm/Kcmr)**2*(Kcmr**2+Xr(i)**2)/(Kcm**2+Xr(i)**2)
       BW(i)=((M(i)**2-mp2)/(W2-mp2))**2*W/M(i)*Gtot*Ggam/(G(i)*((W2-M(i)**2)**2+(M(i)*Gtot)**2))
    ENDDO
    AT=AT0/(1d0+Q2/0.91d0)**cc*(1d0+ac*Q2/(1d0+bc*Q2))
    AL=AL0*Q2/(1d0+dc*Q2)*DEXP(-ec*Q2)
    dW=W-(mp+mpi)
    sT=0d0
    DO i=1,2
       sT=sT+nrT(1,i)/(Q2+nrT(2,i))**(nrT(3,i)+nrT(4,i)*Q2+nrT(5,i)*Q2**2)*dW**(i+0.5d0)
    ENDDO
    sT=sT/(1d0+(W2-(mp+mpi)**2)/(Q2+Q20_T))
    x=Q2/(Q2+W2-mp2); xp=1d0/(1d0+(W2-(mp+mpi)**2)/(Q2+Q20_L))
    t=DLOG(DLOG((Q2+m0_L)/0.33d0**2)/DLOG(m0_L/0.33d0**2))
    sL=nrL(1)*(1d0-xp)**(nrL(2)*t+nrL(3))/(1d0-x)*Q2**nrL(4)/(Q2+Q20_L)**(1d0+nrL(4))*xp**(nrL(5)+nrL(6)*t)
    sT=W*SUM(BW*AT**2)+sT
    sL=W*SUM(BW*AL**2)+sL
    F1=(W2-mp2)*sT/(8d0*pi**2*alpha_QED*ub_GeV)
    FL=2d0*x*F1*sL/sT
    F2_cb=(2d0*x*F1+FL)/(1d0+4d0*mp2*x**2/Q2)
  CONTAINS
    FUNCTION pcm(ww,m) ! meson momentum in the centre of mass
      REAL(KIND(1d0))::pcm
      REAL(KIND(1d0)),INTENT(IN)::ww,m
      pcm=DSQRT(MAX((ww**2-mp2+m**2)**2/(4d0*ww**2)-m**2,0d0))
    END FUNCTION pcm
    SUBROUTINE ratio(p,pr,rr,dd)
      REAL(KIND(1d0)),INTENT(IN)::p,pr
      REAL(KIND(1d0)),INTENT(OUT)::rr,dd
      IF(pr.GT.0d0)THEN
         rr=p/pr; dd=(pr**2+Xr(i)**2)/(p**2+Xr(i)**2)
      ELSE
         rr=0d0; dd=(pr**2+Xr(i)**2)/(p**2+Xr(i)**2)
      ENDIF
    END SUBROUTINE ratio
  END FUNCTION F2_cb

  ! ---------------------------------------------------------------- DGLAP (High_Q2.py)
  SUBROUTINE sf_check()
    ! (re)initialises the PDF set when sf_pdf_id, sf_x_safe or sf_Q2_safe changed (or at the first call)
    IMPLICIT NONE
    IF(.NOT.pdf_init.OR.sf_pdf_id.NE.init_id.OR.sf_x_safe.NE.init_xs.OR.sf_Q2_safe.NE.init_qs)CALL sf_init()
  END SUBROUTINE sf_check

  SUBROUTINE sf_init()
    ! the set of the LHAPDF ID sf_pdf_id (pdfsets.index: the largest set ID <= sf_pdf_id; member = the difference),
    ! opened in LHAPDF slot 1; OrderQCD, MCharm, MBottom, XMin, NumMembers from its .info file; the continuation
    ! threshold x_c, Q2_c (sf_x_safe, sf_Q2_safe)
    IMPLICIT NONE
    CHARACTER(LEN=4096)::paths
    CHARACTER(LEN=512)::line,fname,dir
    CHARACTER(LEN=128)::name
    INTEGER::i0,i1,ios,u,id,best,nmem
    LOGICAL::found
    CALL lhapdf_getdatapath(paths)
    ! the set name from pdfsets.index
    best=-1; name=''
    i0=1
    DO WHILE(i0.LE.LEN_TRIM(paths))
       CALL next_dir()
       OPEN(NEWUNIT=u,FILE=TRIM(dir)//'/pdfsets.index',STATUS='OLD',ACTION='READ',IOSTAT=ios)
       IF(ios.EQ.0)THEN
          DO
             READ(u,'(A)',IOSTAT=ios)line
             IF(ios.NE.0)EXIT
             READ(line,*,IOSTAT=ios)id,fname
             IF(ios.NE.0)CYCLE
             IF(id.LE.sf_pdf_id.AND.id.GT.best)THEN
                best=id; name=fname(1:128)
             ENDIF
          ENDDO
          CLOSE(u)
       ENDIF
    ENDDO
    IF(best.LT.0)THEN
       WRITE(*,*)"ERROR: structure_functions: LHAPDF ID ",sf_pdf_id," not found in pdfsets.index"
       STOP
    ENDIF
    sf_pdf_set=name; pdf_member=sf_pdf_id-best
    ! its .info file
    found=.FALSE.; nmem=-1; pdf_xmin=0d0
    i0=1
    DO WHILE(.NOT.found.AND.i0.LE.LEN_TRIM(paths))
       CALL next_dir()
       fname=TRIM(dir)//'/'//TRIM(sf_pdf_set)//'/'//TRIM(sf_pdf_set)//'.info'
       OPEN(NEWUNIT=u,FILE=TRIM(fname),STATUS='OLD',ACTION='READ',IOSTAT=ios)
       IF(ios.EQ.0)THEN
          found=.TRUE.
          DO
             READ(u,'(A)',IOSTAT=ios)line
             IF(ios.NE.0)EXIT
             IF(INDEX(line,'OrderQCD:').EQ.1)READ(line(10:),*)sf_order
             IF(INDEX(line,'MCharm:').EQ.1)READ(line(8:),*)mcharm
             IF(INDEX(line,'MBottom:').EQ.1)READ(line(9:),*)mbottom
             IF(INDEX(line,'XMin:').EQ.1)READ(line(6:),*)pdf_xmin
             IF(INDEX(line,'NumMembers:').EQ.1)READ(line(12:),*)nmem
          ENDDO
          CLOSE(u)
       ENDIF
    ENDDO
    IF(.NOT.found)THEN
       WRITE(*,*)"ERROR: structure_functions: cannot read the .info file of ",TRIM(sf_pdf_set)
       STOP
    ENDIF
    IF(nmem.GT.0.AND.pdf_member.GE.nmem)THEN
       WRITE(*,*)"ERROR: structure_functions: LHAPDF ID ",sf_pdf_id," is not a member of ",TRIM(sf_pdf_set)
       STOP
    ENDIF
    CALL lhapdf_initpdfset_byname(1,TRIM(sf_pdf_set))
    sf_order=MIN(MAX(sf_order,1),2)
    ! continuation thresholds: XMin (all Q2) and x_c for Q2 < Q2_c (scan or manual)
    WRITE(*,'(A,I8,3A,I4,A,ES9.2)')" INFO: structure functions: LHAPDF ID",sf_pdf_id," (",TRIM(sf_pdf_set),&
         " member",pdf_member,"); PDFs continued below its XMin =",pdf_xmin
    IF(sf_x_safe.LT.0d0)THEN
       CALL pdf_scan(x_c,Q2_c)
       IF(x_c.GT.pdf_xmin)THEN
          WRITE(*,'(A,ES9.2,A,ES9.2,A)')" WARNING: structure functions: the PDF set is not smooth or not positive at"//&
               " small x (scan): continued below x =",x_c," for Q2 <",Q2_c," GeV^2"
       ELSE
          WRITE(*,'(A)')" INFO: structure functions: PDF scan at small x (6 <= Q2 <= 50 GeV^2): smooth and positive"
       ENDIF
    ELSE
       x_c=MAX(pdf_xmin,sf_x_safe); Q2_c=HUGE(1d0)
       IF(sf_Q2_safe.GE.0d0)Q2_c=sf_Q2_safe
       WRITE(*,'(A,ES9.2,A,ES9.2,A)')" INFO: structure functions: PDFs continued below x =",x_c," for Q2 <",&
            MIN(Q2_c,1d99)," GeV^2 (sf_x_safe, sf_Q2_safe)"
    ENDIF
    init_id=sf_pdf_id; init_xs=sf_x_safe; init_qs=sf_Q2_safe
    sf_xc=x_c; sf_Q2c=Q2_c; sf_continued=0; sf_continued_xmin=1d0
    sf_version=sf_version+1
    pdf_warned=.FALSE.
    pdf_init=.TRUE.
  CONTAINS
    SUBROUTINE next_dir() ! the next directory of the ':'-separated paths from i0
      i1=INDEX(paths(i0:),':')
      IF(i1.EQ.0)THEN
         i1=LEN_TRIM(paths)+1
      ELSE
         i1=i0+i1-1
      ENDIF
      dir=paths(i0:i1-1)
      i0=i1+1
    END SUBROUTINE next_dir
  END SUBROUTINE sf_init

  FUNCTION pdf_xf(pid,x,Q2)
    ! x f(x,Q2) of the set (slot 1, pdf_member), 0 for x >= 1; for x < x_c and Q2 < Q2_c the power-law continuation
    ! from x_c (XMin, or the scan/manual x_c below Q2_c: sf_init), with a warning at the first use
    IMPLICIT NONE
    REAL(KIND(1d0))::pdf_xf
    INTEGER,INTENT(IN)::pid
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0))::f0,f1,xc
    pdf_xf=0d0
    IF(x.GE.1d0)RETURN
    xc=pdf_xmin
    IF(Q2.LT.Q2_c)xc=MAX(xc,x_c)
    IF(x.LT.xc)THEN
       !$OMP ATOMIC
       sf_continued=sf_continued+1
       IF(x.LT.sf_continued_xmin)THEN
          !$OMP CRITICAL(sf_xmin)
          sf_continued_xmin=MIN(sf_continued_xmin,x)
          !$OMP END CRITICAL(sf_xmin)
       ENDIF
       IF(.NOT.pdf_warned)THEN
          !$OMP CRITICAL(sf_warning)
          IF(.NOT.pdf_warned)WRITE(*,'(A,ES9.2,A,ES9.2,3A)')" WARNING: structure functions: PDF requested at x =",x,&
               " below x_c =",xc," (",TRIM(sf_pdf_set),"): power-law continuation from x_c used (warned once)"
          pdf_warned=.TRUE.
          !$OMP END CRITICAL(sf_warning)
       ENDIF
       CALL lhapdf_xfxq2(1,pdf_member,pid,xc,Q2,f0)
       CALL lhapdf_xfxq2(1,pdf_member,pid,2d0*xc,Q2,f1)
       IF(f0.GT.0d0.AND.f1.GT.0d0)THEN
          pdf_xf=f0*(x/xc)**(-DLOG(f0/f1)/DLOG(2d0))
       ELSE
          pdf_xf=f0
       ENDIF
       RETURN
    ENDIF
    CALL lhapdf_xfxq2(1,pdf_member,pid,x,Q2,pdf_xf)
  END FUNCTION pdf_xf

  SUBROUTINE pdf_scan(xs,Q2s)
    ! scan of x f for g, u+ubar, d+dbar, s+sbar on x = XMin ... 1e-3 (20 per decade) and Q2 = 6 ... 50 GeV^2 (factor
    ! 1.1): a node is bad when x f <= 0 or |second difference of ln x f in ln x| > 0.2 (smooth small-x PDFs: ~0.01);
    ! xs = 2 x (largest bad x), at least XMin, and Q2s = 1.1 Q2 (largest Q2 with a bad node); xs = XMin, Q2s = 0 if none
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::xs,Q2s
    INTEGER,PARAMETER::NQ=23
    INTEGER,DIMENSION(4),PARAMETER::pids=(/21,2,1,3/)
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::xv,lf
    REAL(KIND(1d0))::Q2,a,b,xbad,qbad
    INTEGER::nx,i,j,p
    xbad=0d0; qbad=0d0
    nx=MAX(INT(20d0*DLOG10(1d-3/MAX(pdf_xmin,1d-12))),2)+1
    ALLOCATE(xv(nx),lf(nx))
    DO i=1,nx
       xv(i)=MAX(pdf_xmin,1d-12)*10d0**(DBLE(i-1)/20d0)
    ENDDO
    DO j=0,NQ-1
       Q2=6d0*1.1d0**j
       DO p=1,4
          DO i=1,nx
             CALL lhapdf_xfxq2(1,pdf_member,pids(p),xv(i),Q2,a)
             IF(pids(p).NE.21)THEN
                CALL lhapdf_xfxq2(1,pdf_member,-pids(p),xv(i),Q2,b)
                a=a+b
             ENDIF
             IF(a.LE.0d0)THEN
                lf(i)=-HUGE(1d0)
                xbad=MAX(xbad,xv(i)); qbad=MAX(qbad,Q2)
             ELSE
                lf(i)=DLOG(a)
             ENDIF
          ENDDO
          DO i=2,nx-1
             IF(lf(i-1).EQ.-HUGE(1d0).OR.lf(i).EQ.-HUGE(1d0).OR.lf(i+1).EQ.-HUGE(1d0))CYCLE
             IF(ABS(lf(i-1)-2d0*lf(i)+lf(i+1)).GT.0.2d0)THEN
                xbad=MAX(xbad,xv(i+1)); qbad=MAX(qbad,Q2)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    IF(xbad.GT.0d0)THEN
       xs=MAX(2d0*xbad,pdf_xmin); Q2s=1.1d0*qbad
    ELSE
       xs=pdf_xmin; Q2s=0d0
    ENDIF
  END SUBROUTINE pdf_scan

  SUBROUTINE F2_FL_DGLAP(x,Q2,F2,FL)
    ! F2, FL from the PDFs at Q (massless coefficient functions, NLO and NNLO for NNLO sets: F2_FL_array of
    ! High_Q2.py, same z nodes: 48-point Gauss-Legendre in ln z on [x, max(x,1/2)] and in ln(1-z) on [max(x,1/2), 1-1e-12])
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x,Q2
    REAL(KIND(1d0)),INTENT(OUT)::F2,FL
    INTEGER,PARAMETER::NGL=48
    REAL(KIND(1d0)),PARAMETER::CF=4d0/3d0,z_eps=1d-12
    REAL(KIND(1d0)),DIMENSION(5),PARAMETER::e2q=(/1d0,4d0,1d0,4d0,1d0/)/9d0
    REAL(KIND(1d0)),DIMENSION(NGL)::tgl,wgl
    REAL(KIND(1d0)),DIMENSION(2*NGL)::z,w
    REAL(KIND(1d0))::zmid,a,a2,Qs,Q1,G,G1,S,E2,nf,L1,Lx,quark,Cq_reg,plus,Cg,zz,ww,qp(5),q1p(5),Fg,FLs
    REAL(KIND(1d0))::sumF2NS,sumF2S,sumFL
    LOGICAL::active(5)
    INTEGER::nz,i,q
    SAVE tgl,wgl
    LOGICAL,SAVE::glinit=.FALSE.
    CALL sf_check()
    IF(.NOT.glinit)THEN
       CALL gauleg(NGL,tgl,wgl)
       glinit=.TRUE.
    ENDIF
    ! z nodes
    zmid=MAX(x,0.5d0)
    nz=0
    IF(zmid.GT.x)THEN
       DO i=1,NGL ! ln z on [ln x, ln zmid]
          zz=0.5d0*(DLOG(x)+DLOG(zmid))+0.5d0*(DLOG(zmid)-DLOG(x))*tgl(i)
          nz=nz+1; z(nz)=DEXP(zz); w(nz)=0.5d0*(DLOG(zmid)-DLOG(x))*wgl(i)*DEXP(zz)
       ENDDO
    ENDIF
    DO i=1,NGL ! u = ln(1-z) on [ln(1-zmid), ln z_eps]
       zz=0.5d0*(DLOG(1d0-zmid)+DLOG(z_eps))+0.5d0*(DLOG(z_eps)-DLOG(1d0-zmid))*tgl(i)
       nz=nz+1; z(nz)=1d0-DEXP(zz); w(nz)=-0.5d0*(DLOG(z_eps)-DLOG(1d0-zmid))*wgl(i)*DEXP(zz)
    ENDDO
    DO q=1,5
       active(q)=q.LT.4.OR.(q.EQ.4.AND.Q2.GT.mcharm**2).OR.(q.EQ.5.AND.Q2.GT.mbottom**2)
    ENDDO
    nf=COUNT(active)
    E2=SUM(e2q,MASK=active)
    CALL lhapdf_alphasq2(1,pdf_member,Q2,a)
    a=a/(2d0*pi)
    DO q=1,5
       q1p(q)=pdfx(q,x)+pdfx(-q,x)
    ENDDO
    Q1=SUM(e2q*q1p,MASK=active)
    G1=pdfx(21,x)
    Lx=DLOG(1d0-x)
    quark=0d0; Fg=0d0; FLs=0d0; sumF2NS=0d0; sumF2S=0d0; sumFL=0d0
    DO i=1,nz
       zz=z(i); ww=w(i)
       DO q=1,5
          qp(q)=pdfx(q,x/zz)+pdfx(-q,x/zz)
       ENDDO
       Qs=SUM(e2q*qp,MASK=active)
       S=SUM(qp,MASK=active)
       G=pdfx(21,x/zz)
       L1=DLOG(1d0-zz)
       Cq_reg=CF*(-(1d0+zz)*L1-(1d0+zz**2)/(1d0-zz)*DLOG(zz)+3d0+2d0*zz)
       plus=CF*(2d0*L1-1.5d0)/(1d0-zz)*(Qs-Q1)
       quark=quark+ww*(Cq_reg*Qs+plus)
       Cg=(zz**2+(1d0-zz)**2)*(L1-DLOG(zz))-1d0+8d0*zz*(1d0-zz)
       Fg=Fg+ww*Cg*G
       FLs=FLs+ww*(CF*2d0*zz*Qs+4d0*zz*(1d0-zz)*E2*G)
       IF(sf_order.GE.2)THEN
          sumF2NS=sumF2NS+ww*(C2NSP2A(zz,nf)*Qs+C2NS2B(zz,nf)*(Qs-Q1))
          sumF2S=sumF2S+ww*(C2PS2A(zz)*S+C2G2A(zz)*G)
          sumFL=sumFL+ww*(CLNSP2A(zz,nf)*Qs+E2*(CLPS2A(zz)*S+CLG2A(zz)*G))
       ENDIF
    ENDDO
    quark=quark+CF*Q1*(Lx**2-1.5d0*Lx-(4.5d0+pi**2/3d0))
    F2=Q1+a*(quark+E2*Fg)
    FL=a*FLs
    IF(sf_order.GE.2)THEN
       a2=(a/2d0)**2
       F2=F2+a2*(sumF2NS+C2NSP2C(x,nf)*Q1+E2*(sumF2S-0.28d0*G1))
       FL=FL+a2*(sumFL-0.164d0*Q1)
    ENDIF
  CONTAINS
    FUNCTION pdfx(pid,xx) ! x f(x,Q2), 0 for x >= 1 (continued at small x: pdf_xf)
      REAL(KIND(1d0))::pdfx
      INTEGER,INTENT(IN)::pid
      REAL(KIND(1d0)),INTENT(IN)::xx
      pdfx=pdf_xf(pid,xx,Q2)
    END FUNCTION pdfx
  END SUBROUTINE F2_FL_DGLAP

  ! NNLO coefficient functions (van Neerven-Vogt; as APFEL ZeroMassCoefficientFunctions; High_Q2.py)
  FUNCTION C2NSP2A(z,nf)
    REAL(KIND(1d0))::C2NSP2A
    REAL(KIND(1d0)),INTENT(IN)::z,nf
    REAL(KIND(1d0))::L,L1
    L=DLOG(z); L1=DLOG(1d0-z)
    C2NSP2A=-69.59d0-1008d0*z-2.835d0*L**3-17.08d0*L**2+5.986d0*L-17.19d0*L1**3+71.08d0*L1**2-660.7d0*L1&
         -174.8d0*L*L1**2+95.09d0*L**2*L1&
         +nf*(-5.691d0-37.91d0*z+2.244d0*L**2+5.770d0*L-1.707d0*L1**2+22.95d0*L1+3.036d0*L**2*L1+17.97d0*L*L1)
  END FUNCTION C2NSP2A
  FUNCTION C2NS2B(z,nf)
    REAL(KIND(1d0))::C2NS2B
    REAL(KIND(1d0)),INTENT(IN)::z,nf
    REAL(KIND(1d0))::L1
    L1=DLOG(1d0-z)
    C2NS2B=(14.2222d0*L1**3-61.3333d0*L1**2-31.105d0*L1+188.64d0+nf*(1.77778d0*L1**2-8.5926d0*L1+6.3489d0))/(1d0-z)
  END FUNCTION C2NS2B
  FUNCTION C2NSP2C(x,nf)
    REAL(KIND(1d0))::C2NSP2C
    REAL(KIND(1d0)),INTENT(IN)::x,nf
    REAL(KIND(1d0))::L1
    L1=DLOG(1d0-x)
    C2NSP2C=3.55555d0*L1**4-20.4444d0*L1**3-15.5525d0*L1**2+188.64d0*L1-338.531d0+0.485d0&
         +nf*(0.592593d0*L1**3-4.2963d0*L1**2+6.3489d0*L1+46.844d0-0.0035d0)
  END FUNCTION C2NSP2C
  FUNCTION C2PS2A(z)
    REAL(KIND(1d0))::C2PS2A
    REAL(KIND(1d0)),INTENT(IN)::z
    REAL(KIND(1d0))::L,L1
    L=DLOG(z); L1=DLOG(1d0-z)
    C2PS2A=5.290d0*(1d0/z-1d0)+4.310d0*L**3-2.086d0*L**2+39.78d0*L-0.101d0*(1d0-z)*L1**3&
         -(24.75d0-13.80d0*z)*L**2*L1+30.23d0*L*L1
  END FUNCTION C2PS2A
  FUNCTION C2G2A(z)
    REAL(KIND(1d0))::C2G2A
    REAL(KIND(1d0)),INTENT(IN)::z
    REAL(KIND(1d0))::L,L1
    L=DLOG(z); L1=DLOG(1d0-z)
    C2G2A=1d0/z*(11.90d0+1494d0*L1)+5.319d0*L**3-59.48d0*L**2-284.8d0*L+392.4d0-1483d0*L1&
         +(6.445d0+209.4d0*(1d0-z))*L1**3-24.00d0*L1**2-724.1d0*L**2*L1-871.8d0*L*L1**2
  END FUNCTION C2G2A
  FUNCTION CLNSP2A(z,nf)
    REAL(KIND(1d0))::CLNSP2A
    REAL(KIND(1d0)),INTENT(IN)::z,nf
    REAL(KIND(1d0))::L,L1
    L=DLOG(z); L1=DLOG(1d0-z)
    CLNSP2A=-40.41d0+97.48d0*z+(26.56d0*z-0.031d0)*L**2-14.85d0*L+13.62d0*L1**2-55.79d0*L1-150.5d0*L*L1&
         +nf*16d0/27d0*(6d0*z*L1-12d0*z*L-25d0*z+6d0)
  END FUNCTION CLNSP2A
  FUNCTION CLPS2A(z)
    REAL(KIND(1d0))::CLPS2A
    REAL(KIND(1d0)),INTENT(IN)::z
    REAL(KIND(1d0))::L,L1
    L=DLOG(z); L1=DLOG(1d0-z)
    CLPS2A=(15.94d0-5.212d0*z)*(1d0-z)**2*L1+(0.421d0+1.520d0*z)*L**2+28.09d0*(1d0-z)*L-(2.370d0/z-19.27d0)*(1d0-z)**3
  END FUNCTION CLPS2A
  FUNCTION CLG2A(z)
    REAL(KIND(1d0))::CLG2A
    REAL(KIND(1d0)),INTENT(IN)::z
    REAL(KIND(1d0))::L,L1
    L=DLOG(z); L1=DLOG(1d0-z)
    CLG2A=(94.74d0-49.20d0*z)*(1d0-z)*L1**2+864.8d0*(1d0-z)*L1+1161d0*z*L*L1+60.06d0*z*L**2&
         +39.66d0*(1d0-z)*L-5.333d0*(1d0/z-1d0)
  END FUNCTION CLG2A

  SUBROUTINE gauleg(n,t,w) ! Gauss-Legendre nodes and weights on [-1,1], increasing nodes (numpy leggauss)
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::t,w
    REAL(KIND(1d0))::z,z1,p1,p2,p3,pp
    INTEGER::i,j,it
    DO i=1,(n+1)/2
       z=DCOS(pi*(i-0.25d0)/(n+0.5d0))
       DO it=1,100
          p1=1d0; p2=0d0
          DO j=1,n
             p3=p2; p2=p1
             p1=((2d0*j-1d0)*z*p2-(j-1d0)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1d0)
          z1=z; z=z1-p1/pp
          IF(ABS(z-z1).LT.1d-16)EXIT
       ENDDO
       t(i)=-z; t(n+1-i)=z
       w(i)=2d0/((1d0-z*z)*pp*pp); w(n+1-i)=w(i)
    ENDDO
  END SUBROUTINE gauleg
END MODULE structure_functions
