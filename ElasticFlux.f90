MODULE ElasticFlux
  ! photon flux of a single elastic (coherent) source: proton or nucleus
  USE OpticalGlauber_Geometry
  USE NINTLIB ! for multiple dimensional integrations
  USE interpolation
  USE bessel_integration
  USE alpha_running
  IMPLICIT NONE
  PRIVATE
  REAL(KIND(1d0)),PUBLIC::alphaem_elasticphoton=0.0072992700729927005d0
  REAL(KIND(1d0)),PUBLIC::aqedup=0.0072992700729927005d0
  REAL(KIND(1d0)),PARAMETER,PUBLIC::LOWER_BFactor_Limit=1D-1
  REAL(KIND(1d0)),PARAMETER,PUBLIC::GeVm12fm=0.1973d0 ! GeV-1 to fm
  REAL(KIND(1d0)),PARAMETER,PUBLIC::mproton=0.938272081d0 ! the mass of proton (GeV), module-wide (new code, structure_functions)
  ! running alpha(Q2)/alpha(0) of the photon (as LUXqed) in each amplitude of the proton photon density with ChFF,
  ! small-x approximation and exact; alpha(Q2) of LbLatNLO at q2=-Q2 without its two-loop QCD term (alpha_running.f90)
  ! default .FALSE.: with it the photon-photon grid build is x22 slower (small-x approximation, x9 exact)
  LOGICAL,PUBLIC::alpha_running4photon=.FALSE.
  ! grids built on demand where the direct evaluation is numerical: the correction beyond the dipole of the exact
  ! distribution (exact_correction_grid) and, in PhotonPhotonFlux, the collinear proton flux; .FALSE.: direct
  LOGICAL,PUBLIC::GENERATE_PROTON_FLUX_GRID=.TRUE.
  ! alpha(Q2)/alpha(0)-1 = sum_j crun(j)*Q2/(Q2+srun(j)), j=0,nrun: alpha_ratio_ph as a sum of poles, for the closed forms
  INCLUDE 'alpha_running_poles.inc'
  ! x and form factors of the exact proton distribution, for its integrands
  REAL(KIND(1d0))::xexact
  INTEGER::iFFexact ! < 0: only the part beyond the dipole
  ! form factor of the correction beyond the dipole in the small-x approximation (1 dde1, 2 luxqed)
  INTEGER::iFFcorr=1
  ! Q2 (GeV**2) above which the integrand of the correction vanishes: G_E=G_dipole above 3.5 for dde1
  ! (see ratio_dde1); for luxqed G_E=rE*G_dipole above the A1 table (9.995), and the correction
  ! is (rE-1)*dipole plus the integral of G_E-rE*G_dipole below 9.995
  REAL(KIND(1d0)),DIMENSION(2),PARAMETER::Q2max_corr=(/3.5d0,9.995d0/)
  REAL(KIND(1d0)),PARAMETER::mupomuN=2.793d0 ! magnetic moment of proton
  PUBLIC::PhotonNumberDensity ! EDFF
  ! ChFF
  PUBLIC::PhotonNumberDensity_ChargeFormFactor_WS,PhotonNumberDensity_ChargeFormFactor_proton
  PUBLIC::PhotonFlux_proton_EDFF_nob,PhotonFlux_proton_ChFF_nob
  PUBLIC::PhotonFlux_nucleus_EDFF_nob,PhotonFlux_nucleus_ChFF_nob
  ! proton without the small-x approximation
  PUBLIC::PhotonNumberDensity_ChargeFormFactor_proton_exact,PhotonFlux_proton_ChFF_exact_nob
  PUBLIC::dipole_running_amplitude_row ! the running-alpha part of the small-x dipole amplitude at all b of one x
  PUBLIC::proton_exact_density_row ! the exact density at all b of one x at once (FFE rows)
  PUBLIC::proton_exact_kspace,proton_exact_edges,proton_exact_comb_P ! (M11) the exact density in momentum space
  PUBLIC::alpha_norm_mu ! alpha(0)/alpha(mu2) of an elastic photon at the scale mu with alpha_running4photon
  ! For generating virtualities of initial photons
  PUBLIC::generate_Q2_epa_proton_iWW,generate_Q2_epa_ion_ChFF
  PUBLIC::progress
CONTAINS
  FUNCTION PhotonNumberDensity(b,Ega,gamma)
    ! It gives us the photon number density with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    REAL(KIND(1d0)),PARAMETER::pi2=9.86960440108935861883449099988d0
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::xx,Egaoga
    REAL(KIND(1d0)),EXTERNAL::BESSK1,BESSK0
    Egaoga=Ega/gamma ! =x_gamma * mN
    xx=Egaoga*b
    PhotonNumberDensity=one/pi2*Egaoga**2*(BESSK1(xx)**2+one/gamma**2*BESSK0(xx)**2)
    RETURN
  END FUNCTION PhotonNumberDensity

  FUNCTION PhotonNumberDensity_ChargeFormFactor_proton(b,Ega,gamma,CORRECTIONQ,MINUSDIPOLESQUAREQ,MU2,RUNAMP)
    ! It gives us the photon number density with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_ChargeFormFactor_proton
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    INTEGER,INTENT(IN),OPTIONAL::CORRECTIONQ ! the correction beyond the dipole: 0 none (default),
                                             ! 1 dde1 (ratio_dde1), 2 luxqed (formfactor_A1)
    LOGICAL,INTENT(IN),OPTIONAL::MINUSDIPOLESQUAREQ ! whether we want to minus the dipole square
                                                    ! the default is no (.FALSE.)
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::MU2 ! scale (GeV**2) of the photon normalisation, see alpha_norm_mu
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::RUNAMP ! dipole_running_amplitude(b,Ega/gamma) when already known (rows)
    REAL(KIND(1d0)),PARAMETER::aa=1.1867816581938533d0 ! in unit of GeV-1 = 1/DSQRT(0.71 GeV^2)
    REAL(KIND(1d0))::btilde,atilde
    REAL(KIND(1d0))::Egaoga
    REAL(KIND(1d0)),EXTERNAL::BESSK1,BESSK0
    REAL(KIND(1d0))::integral,dipole,corr,sqrtterm
    REAL(KIND(1d0)),PARAMETER::pi2=9.86960440108935861883449099988d0
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::logb
    REAL(KIND(1d0)),PARAMETER::eulergamma=0.577215664901532860606512090082d0
    REAL(KIND(1d0)),PARAMETER::logtwo=0.693147180559945309417232121458d0
    IF(b.EQ.0d0)THEN
       PhotonNumberDensity_ChargeFormFactor_proton=0d0
       RETURN
    ENDIF
    Egaoga=Ega/gamma ! = x_gamma*mN
    ! we use the dipole form of proton form factor
    ! ChargeFormFactor_dipole_proton(Q)
    ! the analytic form can be fully integrated out
    ! see eq.(7.18) in /Users/erdissshaw/Works/Manuscript/OpticalGlauber
    ! first, let us define tilde functions
    btilde=b*Egaoga
    atilde=aa*Egaoga
    sqrtterm=DSQRT(1d0+atilde**(-2))
    IF(btilde.GE.1d-3.OR.btilde/atilde.GE.1d-3)THEN
       integral=Egaoga*(BESSK1(btilde)-sqrtterm*BESSK1(btilde*sqrtterm)&
            -btilde/(2d0*atilde**2)*BESSK0(btilde*sqrtterm))
    ELSE
       logb=DLOG(btilde)
       ! it is better to use the Taylor expansion
       ! for log(btilde) terms we sum up to higher order
       integral=(btilde**3/(16d0*atilde**4)+btilde**5*(3D0*atilde**2+2D0)/(384d0*atilde**6)&
            +btilde**7*(6d0*atilde**4+8d0*atilde**2+3d0)/(18432d0*atilde**8))*logb
       integral=integral+btilde*(one-atilde**2*sqrtterm)/(4d0*atilde**2)
       integral=integral-btilde**3/(64d0*atilde**4)*(2D0*(atilde**4-one)*sqrtterm&
            -2D0*atilde**2+(3d0-4d0*eulergamma+4d0*logtwo))
       integral=integral*Egaoga
    ENDIF
    ! the running of alpha (alpha_running4photon)
    IF(PRESENT(RUNAMP))THEN
       integral=integral+RUNAMP
    ELSE
       integral=integral+dipole_running_amplitude(b,Egaoga)
    ENDIF
    dipole=0d0
    IF(PRESENT(MINUSDIPOLESQUAREQ))THEN
       IF(MINUSDIPOLESQUAREQ)THEN
          dipole=integral
       ENDIF
    ENDIF
    corr=0d0
    IF(PRESENT(CORRECTIONQ))THEN
       IF(CORRECTIONQ.GT.0)THEN
          ! include the correction beyond the dipole
          iFFcorr=CORRECTIONQ
          corr=PhotonNumberDensity_ChargeFormFactor_proton_correction(b,Ega,gamma,1)
          integral=integral+corr
       ENDIF
    ENDIF
    IF(PRESENT(MINUSDIPOLESQUAREQ))THEN
       IF(MINUSDIPOLESQUAREQ)THEN
          integral=2d0*dipole*corr+corr**2
       ELSE
          integral=integral**2
       ENDIF
    ELSE
       integral=integral**2
    ENDIF
    PhotonNumberDensity_ChargeFormFactor_proton=one/pi2*integral*alpha_norm_mu(MU2)
    RETURN
  END FUNCTION PhotonNumberDensity_ChargeFormFactor_proton

  FUNCTION alpha_norm_mu(MU2)
    ! the normalisation of an elastic photon at the scale mu (MU2 in GeV**2) with alpha_running4photon: alpha(0)/alpha(MU2),
    ! as the inelastic photon (LUXqed, MSbar at mu): the photon then carries alpha(Q2)**2/alpha(mu2) and the hard process
    ! must use alpha(mu2); the scale only moves alpha between the flux and the hard vertex (the product keeps alpha(Q2)**2
    ! per photon). 1 with fixed alpha or without MU2 (mu = infinity: the photon normalised to alpha(0), hard vertex alpha(0))
    IMPLICIT NONE
    REAL(KIND(1d0))::alpha_norm_mu
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::MU2
    alpha_norm_mu=1d0
    IF(.NOT.alpha_running4photon.OR..NOT.PRESENT(MU2))RETURN
    IF(MU2.GT.0d0.AND.MU2.LT.HUGE(1d0))alpha_norm_mu=1d0/alpha_ratio_tab(MU2) ! (alpha_ratio_ph to 2e-11)
  END FUNCTION alpha_norm_mu

  FUNCTION PhotonNumberDensity_ChargeFormFactor_proton_correction(b,Ega,gamma,integ_method)
    ! It gives us the photon number density with Z=1 and alpha=1 (only the correction part beyond the dipole FF)
    ! b should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_ChargeFormFactor_proton_correction
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule
                                              ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::Egaoga,integral,xga,xga2,Egaoga2,Ega2
    REAL(KIND(1d0))::Egaoga_common,b_common
    COMMON/PND_CFF_proton_corr/Egaoga_common,b_common
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0))::btilde,atilde,db
    INTEGER::init=0
    INTEGER::NXA,NYA,i,j,n,k,l
    SAVE init,NXA,NYA
    INTEGER::iFF_grid=0 ! the form factor and running of the grid (corr_grid_key)
    SAVE iFF_grid
    INTEGER::NYA_save=0
    SAVE NYA_save
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XA,YA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZA,ZA_save
    SAVE XA,YA,ZA,ZA_save
    REAL(KIND(1d0)),PARAMETER::bmaxoR=10d0
    INTEGER,PARAMETER::NBMAX=2
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=12
    ! NXSEG for x_gamma from 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::NXSEG=8
    INTEGER::log10xmin,log10xmin_before,ilog10x
    REAL(KIND(1d0))::XMIN=1D-8
    SAVE XMIN
    REAL(KIND(1d0))::log10x1
    REAL(KIND(1d0))::rescaling_bmax_save
    SAVE rescaling_bmax_save
    INTEGER,PARAMETER::n_interp=6
    REAL(KIND(1d0)),DIMENSION(n_interp)::XA2,YA2
    REAL(KIND(1d0)),PARAMETER::PIPI=3.14159265358979323846264338328d0
    INTEGER::npoints
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::colv
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(nxseg+1)::YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,nxseg+1)::ZD2
    Egaoga=Ega/gamma ! gamma=Ebeam/mproton, xga=Ega/Ebeam, Egaoga=xga*mproton
    xga=Egaoga/mproton
    IF(xga.LE.1d-99)THEN
       PhotonNumberDensity_ChargeFormFactor_proton_correction=0d0
       RETURN
    ENDIF
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF
    IF(init.EQ.1.AND.corr_grid_key().NE.iFF_grid.AND.integ_method6.GE.0)THEN
       ! another form factor or running of alpha: generate the grid again
       init=0
       NYA_save=0
    ENDIF
    IF(init.EQ.0.OR.xga.LT.xmin.OR.(integ_method6.LT.0))THEN
       ! first do a rescaling
       Egaoga_common=Egaoga
       b_common=b
       IF(integ_method6.LT.0)THEN
          IF(b.EQ.0d0)THEN
             PhotonNumberDensity_ChargeFormFactor_proton_correction=0d0
             RETURN
          ELSE
             npoints=50000
             CALL trapezoid_integration(npoints,PND_ChargeFormFactor_proton_corr_fxn,&
                  DSQRT(Q2max_corr(iFFcorr)),integral)
             PhotonNumberDensity_ChargeFormFactor_proton_correction=integral+corr_luxqed_dipole(b,Ega,gamma)
             RETURN
          ENDIF
       ENDIF
       rescaling_bmax_save=MAX(b,bmaxoR*Rproton/GeVm12fm)
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(NXA))
       ENDIF
       log10xmin_before=INT(DLOG10(xmin))
       IF(xga.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(xga))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(xga.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
          IF(ALLOCATED(YA))THEN
             WRITE(*,*)" "
          ENDIF
          WRITE(*,*)"INFO: in ChFF photon number density of proton, the xmin of the grid has been updated to ",xmin
          IF(ALLOCATED(YA))THEN
             DEALLOCATE(YA)
          ENDIF
          IF(ALLOCATED(ZA))THEN
             IF(init.EQ.1)THEN
                IF(ALLOCATED(ZA_save))THEN
                   DEALLOCATE(ZA_save)
                ENDIF
                ALLOCATE(ZA_save(NXA,NYA))
                ! save the values that have been calculated before
                NYA_save=NYA
                DO I=1,NXA
                   DO J=1,NYA
                      ZA_save(I,J)=ZA(I,J)
                   ENDDO
                ENDDO
             ENDIF
             ! then deallocate it
             DEALLOCATE(ZA)
          ENDIF
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       NYA=NXSEG*(-log10xmin)+1
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(NYA))
       ENDIF
       IF(.NOT.ALLOCATED(ZA))THEN
          ALLOCATE(ZA(NXA,NYA))
       ENDIF
       db=9d0/DBLE(NBSEG)
       WRITE(*,*)"INFO: generate grid of photon number density correction from charge form factor"
       WRITE(*,*)"INFO: it will take a few seconds !"
       k=0
       DO n=0,nbmax
          ! from 10**(-n-1)*bmax to 10**(-n)*bmax
          DO i=1,NBSEG
             k=NBSEG*n+i
             ! these are b in unit GeV-1 (not multiplied Egaoga yet !)
             XA(NXA-k+1)=(10d0**(-n-1))*(1d0+DBLE(NBSEG+1-i)*db)*rescaling_bmax_save
          ENDDO
       ENDDO
       IF(k+1.NE.NXA)THEN
          WRITE(*,*)"ERROR: mismatching k+1 and NXA in PhotonNumberDensity_ChargeFormFactor_proton_correction"
          STOP 1
       ENDIF
       XA(1)=0d0
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nxseg
             log10x1=-1d0/DBLE(nxseg)*DBLE(J-1)+DBLE(I)
             K=K+1
             YA(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.NYA-1)THEN
          WRITE(*,*)"ERROR: K != NYA-1"
          STOP 1
       ENDIF
       YA(NYA)=DBLE(log10xmin)
       ! columns of fixed x: the k_T integrals for all b at once (corr_grid_column; was a trapezoid rule of 10000 to
       ! 20000+ points per node)
       ALLOCATE(colv(NXA-1))
       DO J=1,NYA
          IF(NYA_save.EQ.0)THEN
             IF(J.EQ.1)THEN
                CALL progress(INT(J*50d0/NYA),50,.TRUE.)
             ELSE
                CALL progress(INT(J*50d0/NYA),50)
             ENDIF
          ENDIF
          ! use the ones that have been calculated before
          IF(NYA_save.GT.0.AND.J.LE.NYA_save)THEN
             ZA(1:NXA,J)=ZA_save(1:NXA,J)
             CYCLE
          ENDIF
          xga2=10d0**(YA(J)) ! x_gamma
          Egaoga2=xga2*mproton
          Ega2=Egaoga2*gamma
          ! b=0 (XA(1)): the whole row vanishes
          ZA(1,J)=0d0
          CALL corr_grid_column(xga2,NXA-1,XA(2:NXA),colv)
          DO I=2,NXA
             IF(ISNAN(colv(I-1)))THEN
                WRITE(*,*)"ERROR: the integral is not stable (b,Ega,gamma)=",XA(I),Ega2,gamma
                STOP 1
             ENDIF
             ZA(I,J)=colv(I-1)+corr_luxqed_dipole(XA(I),Ega2,gamma)
          ENDDO
       ENDDO
       DEALLOCATE(colv)
       init=1
       iFF_grid=corr_grid_key()
    ENDIF
    IF(b.GT.rescaling_bmax_save.OR.b.LE.0d0)THEN
       PhotonNumberDensity_ChargeFormFactor_proton_correction=0d0
    ELSE
       XI(1)=b
       YI(1)=DLOG10(xga)
       
       db=MIN(b/rescaling_bmax_save,1d0)
       N=-FLOOR(DLOG10(db))-1 ! b is from 10**(-n-1)*bmax to 10**(-n)*bmax
       IF(N.LT.0)THEN
          K=NXA-n_interp
       ELSE
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(k-NBSEG).GT.b)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(k).LT.b)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(i).GE.b)EXIT
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
       IF(YI(1).GE.0d0)THEN
          ilog10x=-1
       ELSE
          ilog10x=FLOOR(YI(1))
       ENDIF
       L=NXSEG*(-ilog10x-1)

       DO I=1,n_interp
          XD2_1D(I)=XA(K+I)
       ENDDO
       DO I=1,NXSEG+1
          YD2_1D(I)=YA(L+I)
       ENDDO
       DO I=1,n_interp
          DO J=1,NXSEG+1
             ZD2(I,J)=ZA(K+I,L+J)
          ENDDO
       ENDDO
       CALL lagrange_interp_2d(n_interp-1,NXSEG,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
       integral=ZI(1)
       PhotonNumberDensity_ChargeFormFactor_proton_correction=integral
    ENDIF
    RETURN
  END FUNCTION PhotonNumberDensity_ChargeFormFactor_proton_correction

  FUNCTION PND_ChargeFormFactor_proton_corr_fxn(kT)
    ! x = ArcTan[kT*gamma/Ega]*2/Pi
    IMPLICIT NONE
    REAL(KIND(1d0))::PND_ChargeFormFactor_proton_corr_fxn
    REAL(KIND(1d0)),INTENT(IN)::kT ! in unit of GeV
    REAL(KIND(1d0))::CFF,pref,bkT,dipole
    REAL(KIND(1d0)),EXTERNAL::BESSJ1
    REAL(KIND(1d0))::Egaoga_common,b_common
    COMMON/PND_CFF_proton_corr/Egaoga_common,b_common
    REAL(KIND(1d0))::Q2 ! ,expo
    REAL(KIND(1d0)),PARAMETER::one=1d0
    ! the parameter in the dipole form factor
    REAL(KIND(1d0)),PARAMETER::aa=1.1867816581938533d0 ! in unit of GeV-1 = 1/DSQRT(0.71 GeV^2)
    !IF(x.GE.1d0.OR.x.LE.0d0)THEN
    IF(kT.LE.0d0)THEN
       PND_ChargeFormFactor_proton_corr_fxn=0d0
       RETURN
    ENDIF
    !pref=DTAN(PIo2*x)
    !kT=pref*Egaoga_common
    bkT=kT*b_common
    pref=kT**2/(kT**2+Egaoga_common**2)
    !pref=pref**2*PIo2*Egaoga_common
    Q2=kT**2+Egaoga_common**2
    IF(Q2.GE.Q2max_corr(iFFcorr))THEN
       PND_ChargeFormFactor_proton_corr_fxn=0d0
       RETURN
    ENDIF
    dipole=one/(one+Q2*aa**2)**2
    IF(iFFcorr.EQ.2)THEN
       CFF=formfactor_A1(Q2,1)-A1_end_ratio(1)*dipole ! G_E-rE*G_dipole
    ELSE
       CFF=dipole*(ratio_dde1(Q2)-one)
    ENDIF
    IF(ISNAN(CFF))THEN
       PRINT *, "ChargeFormFactor is NaN with ",DSQRT(Q2)
       STOP 1
       CFF=0d0
    ENDIF
    IF(alpha_running4photon)CFF=CFF*alpha_ratio_ph(Q2)
    PND_ChargeFormFactor_proton_corr_fxn=pref*BESSJ1(bkT)*CFF
    RETURN
  END FUNCTION PND_ChargeFormFactor_proton_corr_fxn

  FUNCTION corr_grid_key()
    ! the settings of the correction grids: form factor (iFFcorr) and running of alpha
    IMPLICIT NONE
    INTEGER::corr_grid_key
    corr_grid_key=iFFcorr+MERGE(10,0,alpha_running4photon)
    RETURN
  END FUNCTION corr_grid_key

  FUNCTION corr_luxqed_dipole(b,Ega,gamma)
    ! for luxqed (iFFcorr=2), the part (rE-1)*dipole of the correction (in the amplitude), with
    ! rE=G_E/G_dipole above the A1 table
    IMPLICIT NONE
    REAL(KIND(1d0))::corr_luxqed_dipole
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    REAL(KIND(1d0)),PARAMETER::PIPI=3.14159265358979323846264338328d0
    corr_luxqed_dipole=0d0
    IF(iFFcorr.NE.2)RETURN
    corr_luxqed_dipole=(A1_end_ratio(1)-1d0)*&
         PIPI*DSQRT(PhotonNumberDensity_ChargeFormFactor_proton(b,Ega,gamma))
    RETURN
  END FUNCTION corr_luxqed_dipole

  FUNCTION ratio_dde1(Q2)
    ! G_E/G_dipole of the proton used with proton_FF_correction (Q2 in GeV**2):
    ! polynomial a below Q2lim, polynomial b up to Q2fit, dipole (ratio 1) above.
    ! The pieces are joined by cubics matching value and slope across
    ! [Q2lim-dQ2lim,Q2lim+dQ2lim] and [Q2fit-dQ2fit,Q2fit+dQ2fit], so that there are no jumps
    IMPLICIT NONE
    REAL(KIND(1d0))::ratio_dde1
    REAL(KIND(1d0)),INTENT(IN)::Q2
    REAL(KIND(1d0)),PARAMETER::Q2lim=1.12d0,dQ2lim=0.1d0,Q2fit=3d0,dQ2fit=0.5d0
    ! fitted parameters (eq.(7.22) of the gamma-UPC paper had instead
    ! a = 0.999781, -0.347088, 1.96245, -5.3141, 4.91061 and a Gaussian with b = 1.07751, -2.04371, 5.35245)
    REAL(KIND(1d0)),DIMENSION(0:4),PARAMETER::a=(/0.999871d0,-0.215829d0,0.509109d0,-0.621597d0,0.246705d0/)
    REAL(KIND(1d0)),DIMENSION(0:4),PARAMETER::b=(/1.01809d0,-0.778974d-1,-0.184511d-1,0.289159d-2,-0.121585d-3/)
    IF(Q2.LE.Q2lim-dQ2lim)THEN
       ratio_dde1=poly(a,Q2,0)
    ELSEIF(Q2.LT.Q2lim+dQ2lim)THEN
       ratio_dde1=join(Q2,Q2lim-dQ2lim,Q2lim+dQ2lim,poly(a,Q2lim-dQ2lim,0),poly(a,Q2lim-dQ2lim,1),&
            poly(b,Q2lim+dQ2lim,0),poly(b,Q2lim+dQ2lim,1))
    ELSEIF(Q2.LE.Q2fit-dQ2fit)THEN
       ratio_dde1=poly(b,Q2,0)
    ELSEIF(Q2.LT.Q2fit+dQ2fit)THEN
       ratio_dde1=join(Q2,Q2fit-dQ2fit,Q2fit+dQ2fit,poly(b,Q2fit-dQ2fit,0),poly(b,Q2fit-dQ2fit,1),1d0,0d0)
    ELSE
       ratio_dde1=1d0
    ENDIF
    RETURN
  CONTAINS
    FUNCTION poly(c,q,id) ! sum_i c(i) q**i (id=0) or its derivative (id=1)
      REAL(KIND(1d0))::poly
      REAL(KIND(1d0)),DIMENSION(0:4),INTENT(IN)::c
      REAL(KIND(1d0)),INTENT(IN)::q
      INTEGER,INTENT(IN)::id
      IF(id.EQ.0)THEN
         poly=c(0)+q*(c(1)+q*(c(2)+q*(c(3)+q*c(4))))
      ELSE
         poly=c(1)+q*(2d0*c(2)+q*(3d0*c(3)+q*4d0*c(4)))
      ENDIF
    END FUNCTION poly
    FUNCTION join(q,lo,hi,f0,d0,f1,d1) ! cubic with values f0,f1 and slopes d0,d1 at lo,hi
      REAL(KIND(1d0))::join
      REAL(KIND(1d0)),INTENT(IN)::q,lo,hi,f0,d0,f1,d1
      REAL(KIND(1d0))::t,h
      h=hi-lo
      t=(q-lo)/h
      join=((2d0*t-3d0)*t*t+1d0)*f0+((t-2d0)*t+1d0)*t*h*d0+(3d0-2d0*t)*t*t*f1+(t-1d0)*t*t*h*d1
    END FUNCTION join
  END FUNCTION ratio_dde1

  FUNCTION PhotonNumberDensity_ChargeFormFactor_proton_exact(b,Ega,gamma,iFF,integ_method,MU2)
    ! It gives us the photon number density with Z=1 and alpha=1 of the proton without the small-x approximation:
    ! Dirac (F1) and Pauli (F2) form factors at Q2=(kT**2+x**2*mproton**2)/(1-x), b~=b/(1-x)
    ! (eqs.(44)-(46) of Magnetic_EPA/coherent_short.tex; xfxbQ2_GPDel in the Python code)
    ! b should be written in unit of GeV-1
    ! iFF: proton form factors, 0: dipole, 1: dde1 (as proton_FF_correction), 2: luxqed (A1 fit)
    ! the dipole part is analytic (exact_dipole_amplitudes), the correction beyond it is integrated
    ! numerically up to b = bmaxoR*Rproton (as in PhotonNumberDensity_ChargeFormFactor_proton_correction)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_ChargeFormFactor_proton_exact
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    INTEGER,INTENT(IN)::iFF
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! direct numerical integration of all when integ_method < 0
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::MU2 ! scale (GeV**2) of the photon normalisation, see alpha_norm_mu
    REAL(KIND(1d0)),PARAMETER::pi2=9.86960440108935861883449099988d0
    REAL(KIND(1d0)),PARAMETER::bmaxoR=10d0
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    INTEGER,DIMENSION(6),PARAMETER::nj=(/1,0,0,2,1,0/) ! the orders of the Bessel functions
    REAL(KIND(1d0)),DIMENSION(6)::integrals,corr
    LOGICAL::direct
    PhotonNumberDensity_ChargeFormFactor_proton_exact=0d0
    xexact=Ega/(gamma*mproton)
    IF(xexact.GE.1d0)RETURN
    direct=b.LE.0d0
    IF(PRESENT(integ_method))direct=direct.OR.integ_method.LT.0
    ! the six integrals int dkT kT**nk/D F_i J_n(b~ kT) of exact_integrands
    IF(direct)THEN
       iFFexact=iFF
       CALL bessel_integrate(exact_integrands,6,nj,b/(1d0-xexact),xexact*mproton,HUGE(1d0),&
            exact_edges(),integrals)
    ELSE
       CALL exact_dipole_amplitudes(b/(1d0-xexact),1d0,mupomuN,integrals)
       IF(iFF.NE.0.AND.b.LE.bmaxoR*Rproton/GeVm12fm)THEN
          ! F_i - F_i(dipole), zero above Q2max_corr for dde1; for luxqed F_i - F_i(rescaled dipole) is zero
          ! above Q2max_corr and the rescaled dipole (G_E=rE*G_D, G_M=rM*mupomuN*G_D) minus the dipole is analytic
          IF(GENERATE_PROTON_FLUX_GRID)THEN
             corr=exact_correction_grid(b/(1d0-xexact),xexact,iFF)
          ELSE
             corr=exact_correction(b/(1d0-xexact),xexact,iFF)
          ENDIF
          integrals=integrals+corr
          IF(iFF.EQ.2)THEN
             CALL exact_dipole_amplitudes(b/(1d0-xexact),A1_end_ratio(1)-1d0,mupomuN*(A1_end_ratio(2)-1d0),corr)
             integrals=integrals+corr
          ENDIF
       ENDIF
    ENDIF
    PhotonNumberDensity_ChargeFormFactor_proton_exact=exact_combination(integrals)/pi2/(1d0-xexact)**2*alpha_norm_mu(MU2)
    RETURN
  END FUNCTION PhotonNumberDensity_ChargeFormFactor_proton_exact

  FUNCTION exact_correction(bt,x,iFF)
    ! the six integrals of exact_integrands with F_i - F_i(dipole) (rescaled dipole for luxqed), zero above
    ! Q2max_corr, at b~=bt and x (sets xexact and iFFexact)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(6)::exact_correction
    REAL(KIND(1d0)),INTENT(IN)::bt,x
    INTEGER,INTENT(IN)::iFF
    INTEGER,DIMENSION(6),PARAMETER::nj=(/1,0,0,2,1,0/) ! the orders of the Bessel functions
    xexact=x
    iFFexact=-iFF
    CALL bessel_integrate(exact_integrands,6,nj,bt,x*mproton,DSQRT(MAX(Q2max_corr(iFF)*(1d0-x)-(x*mproton)**2,0d0)),&
         exact_edges(),exact_correction)
    RETURN
  END FUNCTION exact_correction

  FUNCTION exact_correction_grid(bt,x,iFF)
    ! exact_correction from a grid in v=log(x/(1-x)) and t=log(bt)+bt/BETA (log at small b~, uniform steps BETA*HT
    ! at large b~, where the correction oscillates with the joins of the form factors), XCMIN<=x<=XCMAX,
    ! BTMIN<=bt<=BTMAX; 6-point Lagrange interpolation in v and t; a node is computed when first needed (all six
    ! integrals, without the cut at bmaxoR*Rproton, applied by the caller). Direct evaluation outside.
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(6)::exact_correction_grid
    REAL(KIND(1d0)),INTENT(IN)::bt,x
    INTEGER,INTENT(IN)::iFF
    REAL(KIND(1d0)),PARAMETER::XCMIN=1d-9,XCMAX=0.3d0,BTMIN=1d-2,BTMAX=1d2,BETA=1d0 ! x<=0.3: luxqed 2e-5 at x>0.35
    REAL(KIND(1d0)),PARAMETER::HX=0.1151292546497023d0,HT=0.2d0 ! 20 per decade in x; 0.2 GeV-1 at large b~
    INTEGER,PARAMETER::NP=6,NXC=185,NTC=550 ! v(NXC)>log(XCMAX/(1-XCMAX)); t(NTC)>log(BTMAX)+BTMAX/BETA
    REAL(KIND(1d0)),DIMENSION(6,NXC,NTC)::C
    LOGICAL,DIMENSION(NXC,NTC)::done
    INTEGER::gkey=-1
    SAVE C,done,gkey
    REAL(KIND(1d0))::v,t,v0,t0,wx(NP),wt(NP),xsave,xn,bn
    INTEGER::i0,j0,i,j,k,isave
    IF(x.LT.XCMIN.OR.x.GT.XCMAX.OR.bt.LT.BTMIN.OR.bt.GT.BTMAX)THEN
       exact_correction_grid=exact_correction(bt,x,iFF)
       RETURN
    ENDIF
    IF(gkey.NE.iFF+MERGE(10,0,alpha_running4photon))THEN
       done=.FALSE.
       gkey=iFF+MERGE(10,0,alpha_running4photon)
    ENDIF
    v0=DLOG(XCMIN/(1d0-XCMIN)); t0=DLOG(BTMIN)+BTMIN/BETA
    v=DLOG(x/(1d0-x)); t=DLOG(bt)+bt/BETA
    i0=MIN(MAX(FLOOR((v-v0)/HX)+1-(NP/2-1),1),NXC-NP+1)
    j0=MIN(MAX(FLOOR((t-t0)/HT)+1-(NP/2-1),1),NTC-NP+1)
    xsave=xexact; isave=iFFexact
    DO i=i0,i0+NP-1
       DO j=j0,j0+NP-1
          IF(.NOT.done(i,j))THEN
             xn=1d0/(1d0+DEXP(-(v0+(i-1)*HX))) ! x of the node
             bn=bt_of_t(t0+(j-1)*HT)
             C(:,i,j)=exact_correction(bn,xn,iFF)
             done(i,j)=.TRUE.
          ENDIF
       ENDDO
    ENDDO
    xexact=xsave; iFFexact=isave
    wx=lagrange_weights(v,v0+(i0-1)*HX,HX)
    wt=lagrange_weights(t,t0+(j0-1)*HT,HT)
    DO k=1,6
       exact_correction_grid(k)=DOT_PRODUCT(wx,MATMUL(C(k,i0:i0+NP-1,j0:j0+NP-1),wt))
    ENDDO
    RETURN
  CONTAINS
    FUNCTION bt_of_t(tt) ! inverse of t=log(b~)+b~/BETA (Newton)
      REAL(KIND(1d0))::bt_of_t
      REAL(KIND(1d0)),INTENT(IN)::tt
      INTEGER::n
      bt_of_t=MIN(DEXP(tt),BETA*MAX(tt,1d0))
      DO n=1,100
         bt_of_t=bt_of_t-(DLOG(bt_of_t)+bt_of_t/BETA-tt)/(1d0/bt_of_t+1d0/BETA)
      ENDDO
    END FUNCTION bt_of_t
    FUNCTION lagrange_weights(y,y1,h) ! weights of the NP nodes y1+(k-1)*h at y
      REAL(KIND(1d0)),DIMENSION(NP)::lagrange_weights
      REAL(KIND(1d0)),INTENT(IN)::y,y1,h
      INTEGER::kk,l
      lagrange_weights=1d0
      DO kk=1,NP
         DO l=1,NP
            IF(l.NE.kk)lagrange_weights(kk)=lagrange_weights(kk)*(y-y1-(l-1)*h)/((kk-l)*h)
         ENDDO
      ENDDO
    END FUNCTION lagrange_weights
  END FUNCTION exact_correction_grid

  FUNCTION PhotonFlux_proton_ChFF_exact_nob(x,iFF,MU2)
    ! set PNOHARD=1
    ! for proton with the exact distribution: int d^2b of PhotonNumberDensity_ChargeFormFactor_proton_exact
    ! by Parseval, alpha/x*2/pi*int dkT/kT of exact_combination with the integrals replaced by their integrands
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_ChFF_exact_nob
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::iFF
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::MU2 ! scale (GeV**2) of the photon normalisation, see alpha_norm_mu
    REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
    REAL(KIND(1d0))::alpha,integral(1)
    alpha=alphaem_elasticphoton
    IF(alpha.LT.0d0)alpha=MERGE(aqedup,0.0072992701d0,aqedup.GT.0d0)
    PhotonFlux_proton_ChFF_exact_nob=0d0
    IF(x.GE.1d0)RETURN
    xexact=x
    iFFexact=iFF
    CALL bessel_integrate(exact_collinear_integrand,1,(/-1/),0d0,x*mproton,HUGE(1d0),exact_edges(),integral)
    PhotonFlux_proton_ChFF_exact_nob=alpha/x*2d0/pi*integral(1)*alpha_norm_mu(MU2)
    RETURN
  END FUNCTION PhotonFlux_proton_ChFF_exact_nob


  SUBROUTINE proton_exact_kspace(x,iFF,nk,kv,V)
    ! (M11) the six k_T integrands of the exact distribution (exact_integrands) at kv(1:nk): V(6,nk)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::iFF,nk
    REAL(KIND(1d0)),DIMENSION(nk),INTENT(IN)::kv
    REAL(KIND(1d0)),DIMENSION(6,nk),INTENT(OUT)::V
    INTEGER::i
    xexact=x
    iFFexact=iFF
    DO i=1,nk
       CALL exact_integrands(kv(i),6,V(:,i))
    ENDDO
  END SUBROUTINE proton_exact_kspace

  FUNCTION proton_exact_edges(x,iFF)
    ! (M11) the k_T of the joins of the form factors (exact_edges)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::proton_exact_edges
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::iFF
    xexact=x
    iFFexact=iFF
    proton_exact_edges=exact_edges()
  END FUNCTION proton_exact_edges

  FUNCTION proton_exact_comb_P(x,P)
    ! (M11) exact_combination with the products I(a) I(b) replaced by the pair overlaps P(a,b) (symmetric)
    IMPLICIT NONE
    REAL(KIND(1d0))::proton_exact_comb_P
    REAL(KIND(1d0)),INTENT(IN)::x
    REAL(KIND(1d0)),DIMENSION(6,6),INTENT(IN)::P
    REAL(KIND(1d0))::mp2
    mp2=mproton**2
    proton_exact_comb_P=(1d0+(1d0-x)**2)/2d0*P(1,1)+(P(3,3)+P(4,4))/(8d0*mp2)&
         +x**2/2d0*(P(3,6)+P(5,5)/2d0+P(1,5)+P(2,3)+2d0*x**2*mp2*P(2,6))&
         +x**4*mp2/2d0*(P(2,2)+P(6,6))
  END FUNCTION proton_exact_comb_P

  SUBROUTINE exact_integrands(kT,n,v)
    ! kT**nk/D*F_i(Q2) (times alpha(Q2)/alpha(0) with alpha_running4photon) for
    ! (F_i,nk) = (F1,2),(F1,1),(F2,3),(F2,3),(F2,2),(F2,1), D=kT**2+x**2*mproton**2, Q2=D/(1-x)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::kT
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::v
    REAL(KIND(1d0))::D,Q2,F1,F2,F1D,F2D,ar
    D=kT**2+(xexact*mproton)**2
    Q2=D/(1d0-xexact)
    CALL proton_F1F2(Q2,ABS(iFFexact),F1,F2)
    IF(iFFexact.LT.0)THEN ! only the part beyond the dipole (rescaled for luxqed)
       CALL proton_F1F2(Q2,MERGE(3,0,iFFexact.EQ.-2),F1D,F2D)
       F1=F1-F1D
       F2=F2-F2D
    ENDIF
    IF(alpha_running4photon)THEN
       ar=alpha_ratio_tab(Q2)
       F1=F1*ar
       F2=F2*ar
    ENDIF
    v=(/kT**2*F1,kT*F1,kT**3*F2,kT**3*F2,kT**2*F2,kT*F2/)/D
    RETURN
  END SUBROUTINE exact_integrands

  SUBROUTINE dipole_running_collinear_integrand(kT,n,v)
    ! (kT**2/Q2 G_D(Q2) alpha(Q2)/alpha(0))**2/kT, Q2=kT**2+x**2*mproton**2 (small-x approximation)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::kT
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::v
    REAL(KIND(1d0))::Q2
    Q2=kT**2+(xexact*mproton)**2
    v(1)=(kT**2/Q2/(1d0+Q2/0.71d0)**2*alpha_ratio_ph(Q2))**2/kT
    RETURN
  END SUBROUTINE dipole_running_collinear_integrand

  SUBROUTINE exact_collinear_integrand(kT,n,v)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::kT
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::v
    REAL(KIND(1d0)),DIMENSION(6)::h
    CALL exact_integrands(kT,6,h)
    v(1)=exact_combination(h)/kT
    RETURN
  END SUBROUTINE exact_collinear_integrand

  FUNCTION exact_combination(I)
    ! the bracket of eqs.(44)-(46), in the integrals I of exact_integrands
    IMPLICIT NONE
    REAL(KIND(1d0))::exact_combination
    REAL(KIND(1d0)),DIMENSION(6),INTENT(IN)::I ! F1J1,F1J0,F2J0_k3,F2J2,F2J1,F2J0_k
    REAL(KIND(1d0))::x,mp2
    x=xexact
    mp2=mproton**2
    exact_combination=(1d0+(1d0-x)**2)/2d0*I(1)**2+(I(3)**2+I(4)**2)/(8d0*mp2)&
         +x**2/2d0*(I(3)*I(6)+I(5)**2/2d0+I(1)*I(5)+I(2)*(I(3)+2d0*x**2*mp2*I(6)))&
         +x**4*mp2/2d0*(I(2)**2+I(6)**2)
    RETURN
  END FUNCTION exact_combination

  FUNCTION exact_edges()
    ! kT where the form factors at Q2=(kT**2+x**2*mproton**2)/(1-x) are not smooth (panel edges)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::exact_edges
    REAL(KIND(1d0)),DIMENSION(4),PARAMETER::Q2dde1=(/1.02d0,1.22d0,2.5d0,3.5d0/) ! joins of ratio_dde1
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::Q2A1=(/0.25d0,0.5d0,0.75d0,1d0,1.5d0,3d0,9.995d0/)
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::k2
    SELECT CASE(ABS(iFFexact))
    CASE(1)
       k2=(1d0-xexact)*Q2dde1-(xexact*mproton)**2
    CASE(2)
       k2=(1d0-xexact)*Q2A1-(xexact*mproton)**2
    CASE DEFAULT
       ALLOCATE(k2(0))
    END SELECT
    exact_edges=DSQRT(PACK(k2,k2.GT.0d0))
    RETURN
  END FUNCTION exact_edges

  SUBROUTINE proton_F1F2(Q2,iFF,F1,F2)
    ! Dirac and Pauli form factors of the proton from G_E=F1-tau*F2, G_M=F1+F2 (Q2 in GeV**2)
    ! iFF=0: dipole, G_E=G_D, G_M=mup*G_D; 1: G_E=G_D*ratio_dde1, G_M=mup*G_D; 2: A1 fit (luxqed);
    ! 3: G_E=rE*G_D, G_M=rM*mup*G_D (A1_end_ratio)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::Q2
    INTEGER,INTENT(IN)::iFF
    REAL(KIND(1d0)),INTENT(OUT)::F1,F2
    REAL(KIND(1d0))::tau,GD,GE,GM
    tau=Q2/(4d0*mproton**2)
    GD=1d0/(1d0+Q2/0.71d0)**2
    SELECT CASE(iFF)
    CASE(1)
       GE=GD*ratio_dde1(Q2)
       GM=mupomuN*GD
    CASE(2)
       GE=formfactor_A1(Q2,1)
       GM=mupomuN*formfactor_A1(Q2,2)
    CASE(3) ! the dipole rescaled to the A1 fit above its table
       GE=A1_end_ratio(1)*GD
       GM=mupomuN*A1_end_ratio(2)*GD
    CASE DEFAULT
       GE=GD
       GM=mupomuN*GD
    END SELECT
    F1=(GE+tau*GM)/(1d0+tau)
    F2=(GM-GE)/(1d0+tau)
    RETURN
  END SUBROUTINE proton_F1F2

  FUNCTION formfactor_A1(Q2,i)
    ! G_E (i=1) and G_M/mup (i=2) of the A1 world fit with polarised data (as LUXqed),
    ! cubic spline of the tabulated values; above the table the dipole shape normalised at its end
    IMPLICIT NONE
    REAL(KIND(1d0))::formfactor_A1
    REAL(KIND(1d0)),INTENT(IN)::Q2
    INTEGER,INTENT(IN)::i
    REAL(KIND(1d0)),DIMENSION(0:1000)::a1GE,a1GM
    INCLUDE 'proton_A1_formfactors.inc'
    REAL(KIND(1d0)),DIMENSION(0:1000)::q
    REAL(KIND(1d0)),DIMENSION(0:1000,2)::bs,cs,ds
    INTEGER::init=0,j
    SAVE init,q,bs,cs,ds
    IF(init.EQ.0)THEN
       q=(/0d0,(0.005d0+0.01d0*(j-1),j=1,1000)/)
       CALL spline(q,a1GE,bs(:,1),cs(:,1),ds(:,1),1001)
       CALL spline(q,a1GM,bs(:,2),cs(:,2),ds(:,2),1001)
       init=1
    ENDIF
    IF(Q2.GE.q(1000))THEN
       formfactor_A1=MERGE(a1GE(1000),a1GM(1000),i.EQ.1)*((1d0+q(1000)/0.71d0)/(1d0+Q2/0.71d0))**2
    ELSEIF(i.EQ.1)THEN
       formfactor_A1=ispline(Q2,q,a1GE,bs(:,1),cs(:,1),ds(:,1),1001)
    ELSE
       formfactor_A1=ispline(Q2,q,a1GM,bs(:,2),cs(:,2),ds(:,2),1001)
    ENDIF
    RETURN
  END FUNCTION formfactor_A1

  FUNCTION A1_end_ratio(i)
    ! G_E/G_D (i=1) and G_M/(mup*G_D) (i=2) of formfactor_A1 above the A1 table
    IMPLICIT NONE
    REAL(KIND(1d0))::A1_end_ratio
    INTEGER,INTENT(IN)::i
    A1_end_ratio=formfactor_A1(1d2,i)*(1d0+1d2/0.71d0)**2
    RETURN
  END FUNCTION A1_end_ratio

  SUBROUTINE partial_fractions(p,a2,n0,n1,C)
    ! (n0+n1*y)/((y+a2)**2*prod_i(y+p(i))) = sum_i C(i)/(y+p(i)) + C(n+1)/(y+a2) + C(n+2)/(y+a2)**2, n=SIZE(p)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::p
    REAL(KIND(1d0)),INTENT(IN)::a2,n0,n1
    REAL(KIND(1d0)),DIMENSION(SIZE(p)+2),INTENT(OUT)::C
    INTEGER::i,l,n
    n=SIZE(p)
    DO i=1,n
       C(i)=(n0-n1*p(i))/((a2-p(i))**2*PRODUCT(p-p(i),MASK=(/(l.NE.i,l=1,n)/)))
    ENDDO
    C(n+2)=(n0-n1*a2)/PRODUCT(p-a2)
    C(n+1)=(n1-C(n+2)*PRODUCT(p-a2)*SUM(1d0/(p-a2)))/PRODUCT(p-a2)
    RETURN
  END SUBROUTINE partial_fractions

  SUBROUTINE exact_dipole_amplitudes(bt,cE,cM,I)
    ! the integrals of exact_integrands at b~=bt in closed form, for G_E=cE*G_D and G_M=cM*G_D:
    ! F1 and F2 are rational in y=kT**2, with poles at y=-a2 (double, G_D) and -c2 (1+tau), and 1/D=1/(y+w2);
    ! by partial fractions the integrals are sums of int dkT kT**p J_n(bt kT)/(kT**2+m2)**(1 or 2),
    ! which are modified Bessel functions K_n (exact_dipole_basis).
    ! With alpha_running4photon, alpha(Q2)/alpha(0)=1+sum_j crun(j) Q2/(Q2+srun(j)) (alpha_running_poles.inc):
    ! each term gives the same integrals with 1/D replaced by 1/(y+w2+(1-x)*srun(j))
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::bt,cE,cM
    REAL(KIND(1d0)),DIMENSION(6),INTENT(OUT)::I
    REAL(KIND(1d0)),PARAMETER::Lambda2=0.71d0 ! the parameter in the dipole form factor (GeV**2)
    REAL(KIND(1d0))::s,w2,a2,c2,P,Q4,e2
    REAL(KIND(1d0)),DIMENSION(0:2)::Ka,Kc
    INTEGER::j
    s=1d0-xexact
    w2=(xexact*mproton)**2
    a2=w2+s*Lambda2
    c2=w2+4d0*s*mproton**2
    P=(s*Lambda2)**2
    Q4=4d0*s*mproton**2
    CALL besselK012(DSQRT(a2)*bt,Ka)
    CALL besselK012(DSQRT(c2)*bt,Kc)
    I=amplitudes(w2)
    IF(alpha_running4photon)THEN
       DO j=0,nrun
          e2=w2+s*srun(j)
          I=I+crun(j)*amplitudes(e2)
       ENDDO
    ENDIF
    RETURN
  CONTAINS
    FUNCTION amplitudes(q) ! the six integrals with 1/(y+q) in place of 1/D
      REAL(KIND(1d0)),DIMENSION(6)::amplitudes
      REAL(KIND(1d0)),INTENT(IN)::q
      REAL(KIND(1d0)),DIMENSION(0:2)::Kq
      REAL(KIND(1d0)),DIMENSION(4)::CF1,CF2 ! coefficients of 1/(y+c2),1/(y+q),1/(y+a2),1/(y+a2)**2 for F1, F2
      REAL(KIND(1d0)),DIMENSION(4,4)::L ! L(t,j): basis integral t for the j-th term
      INTEGER::t
      CALL besselK012(DSQRT(q)*bt,Kq)
      CALL partial_fractions((/c2,q/),a2,P*(cM*w2+Q4*cE),P*cM,CF1)
      CALL partial_fractions((/c2,q/),a2,P*(cM-cE)*Q4,0d0,CF2)
      DO t=1,4
         L(t,:)=(/exact_dipole_basis(t,c2,bt,Kc,.FALSE.),exact_dipole_basis(t,q,bt,Kq,.FALSE.),&
              exact_dipole_basis(t,a2,bt,Ka,.FALSE.),exact_dipole_basis(t,a2,bt,Ka,.TRUE.)/)
      ENDDO
      amplitudes=(/DOT_PRODUCT(CF1,L(1,:)),DOT_PRODUCT(CF1,L(2,:)),DOT_PRODUCT(CF2,L(3,:)),&
           DOT_PRODUCT(CF2,L(4,:)),DOT_PRODUCT(CF2,L(1,:)),DOT_PRODUCT(CF2,L(2,:))/)
    END FUNCTION amplitudes
  END SUBROUTINE exact_dipole_amplitudes

  FUNCTION dipole_running_amplitude(b,Egaoga)
    ! in the small-x approximation with alpha_running4photon, the part of the dipole amplitude
    ! int dkT kT**2/(kT**2+Egaoga**2) G_D(Q2) J1(b kT) [alpha(Q2)/alpha(0)-1], Q2=kT**2+Egaoga**2, in closed form
    ! (as exact_dipole_amplitudes): dipole_running_amplitude_row at one b
    IMPLICIT NONE
    REAL(KIND(1d0))::dipole_running_amplitude
    REAL(KIND(1d0)),INTENT(IN)::b,Egaoga
    REAL(KIND(1d0)),DIMENSION(1)::A
    dipole_running_amplitude=0d0
    IF(.NOT.alpha_running4photon.OR.b.LE.0d0)RETURN
    CALL dipole_running_amplitude_row(Egaoga,1,(/b/),A)
    dipole_running_amplitude=A(1)
    RETURN
  END FUNCTION dipole_running_amplitude

  SUBROUTINE dipole_running_amplitude_row(Egaoga,n,b,A)
    ! dipole_running_amplitude at all b(1:n) > 0 of one Egaoga = x mp. With the poles of alpha_running_poles.inc,
    ! sum_j crun(j) [C1(j) m_j K_1(m_j b) + C2(j) a K_1(a b) + C3(j) b K_0(a b)/2], m_j**2 = Egaoga**2+srun(j),
    ! a**2 = Egaoga**2+0.71 (partial fractions of G_D/(y+m_j**2)): the a terms summed over j once per x; a pole term
    ! is left out where it is below exp(-50) of the j=0 term (m_j-m_0) b > 50 (K_1(z) ~ exp(-z): < 2e-22 of it)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::Egaoga
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(IN)::b
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::A
    REAL(KIND(1d0)),PARAMETER::Lambda2=0.71d0 ! the parameter in the dipole form factor (GeV**2)
    REAL(KIND(1d0)),PARAMETER::ZSKIP=50d0
    REAL(KIND(1d0))::a2,S2,S3
    REAL(KIND(1d0)),DIMENSION(0:nrun)::e2,m,C1
    REAL(KIND(1d0)),DIMENSION(0:2)::Ka,Kq
    REAL(KIND(1d0)),DIMENSION(3)::C
    INTEGER::j,k
    A=0d0
    IF(.NOT.alpha_running4photon)RETURN
    a2=Egaoga**2+Lambda2
    S2=0d0; S3=0d0
    DO j=0,nrun
       e2(j)=Egaoga**2+srun(j)
       m(j)=DSQRT(e2(j))
       CALL partial_fractions((/e2(j)/),a2,Lambda2**2,0d0,C)
       C1(j)=crun(j)*C(1)
       S2=S2+crun(j)*C(2)
       S3=S3+crun(j)*C(3)
    ENDDO
    DO k=1,n
       IF(b(k).LE.0d0)CYCLE
       CALL besselK012(DSQRT(a2)*b(k),Ka)
       A(k)=S2*exact_dipole_basis(1,a2,b(k),Ka,.FALSE.)+S3*exact_dipole_basis(1,a2,b(k),Ka,.TRUE.)
       DO j=0,nrun
          IF((m(j)-m(0))*b(k).GT.ZSKIP)EXIT ! m(j) increases with j
          CALL besselK012(m(j)*b(k),Kq)
          A(k)=A(k)+C1(j)*m(j)*Kq(1)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE dipole_running_amplitude_row

  FUNCTION exact_dipole_basis(t,m2,bt,K,double)
    ! int_0^infinity dkT kT**p J_n(bt kT)/(kT**2+m2) for t=1: (n,p)=(1,2), 2: (0,1), 3: (0,3), 4: (2,3)
    ! (t=3,4 up to terms that cancel in exact_dipole_amplitudes), or with (kT**2+m2)**2 when double
    IMPLICIT NONE
    REAL(KIND(1d0))::exact_dipole_basis
    INTEGER,INTENT(IN)::t
    REAL(KIND(1d0)),INTENT(IN)::m2,bt
    REAL(KIND(1d0)),DIMENSION(0:2),INTENT(IN)::K ! K_n(sqrt(m2)*bt) (besselK012)
    LOGICAL,INTENT(IN)::double
    REAL(KIND(1d0))::m,z
    m=DSQRT(m2)
    z=m*bt
    SELECT CASE(MERGE(t+4,t,double))
    CASE(1); exact_dipole_basis=m*K(1)
    CASE(2); exact_dipole_basis=K(0)
    CASE(3); exact_dipole_basis=-m2*K(0)
    CASE(4); exact_dipole_basis=m2*K(2)
    CASE(5); exact_dipole_basis=bt*K(0)/2d0
    CASE(6); exact_dipole_basis=bt*K(1)/(2d0*m)
    CASE(7); exact_dipole_basis=K(0)-z*K(1)/2d0
    CASE DEFAULT; exact_dipole_basis=z*K(1)/2d0
    END SELECT
    RETURN
  END FUNCTION exact_dipole_basis

  SUBROUTINE proton_exact_density_row(x,iFF,n,b,dens)
    ! PhotonNumberDensity_ChargeFormFactor_proton_exact (alpha(0) normalisation) at all b(1:n) (ascending) of one x:
    ! the closed-form dipole amplitudes (exact_dipole_amplitudes) with the partial fractions and the running-alpha
    ! poles computed once for the row (linear in G_E, G_M: the rescaled dipole of luxqed at no cost), and the
    ! correction beyond the dipole (iFF > 0, b <= bmaxoR R_p, k_T up to the end of the fit, exact_correction) for all b
    ! at once on shared Gauss-Legendre nodes in k_T (panels between the joins of the form factors, width
    ! <= min(PH/b~, 0.25 GeV) for the b of each octave of b~)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER,INTENT(IN)::iFF,n
    REAL(KIND(1d0)),DIMENSION(n),INTENT(IN)::b
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::dens
    REAL(KIND(1d0)),PARAMETER::pi2=9.86960440108935861883449099988d0,Lambda2=0.71d0
    REAL(KIND(1d0)),PARAMETER::bmaxoR=10d0,Rproton=0.877d0
    REAL(KIND(1d0)),PARAMETER::PH=10d0 ! largest phase b~ k_T across a fine panel of NG Gauss-Legendre nodes
    REAL(KIND(1d0)),PARAMETER::PC=12d0 ! largest phase across a coarse panel of NC Chebyshev points (J_n: ~1e-12)
    REAL(KIND(1d0)),PARAMETER::HF1=0.25d0,HF2=0.05d0 ! fine panel width (GeV), dde1 / luxqed (A1 spline)
    INTEGER,PARAMETER::NG=20,NC=24
    REAL(KIND(1d0)),DIMENSION(4,0:nrun+1)::CU1E,CU1M,CU2E,CU2M ! partial fractions for (G_E,G_M)/G_D=(1,0),(0,1)
    REAL(KIND(1d0)),DIMENSION(0:nrun+1)::qv,wq
    REAL(KIND(1d0)),DIMENSION(6,n)::I,Ic
    REAL(KIND(1d0)),DIMENSION(6)::AE,AM,v
    REAL(KIND(1d0)),DIMENSION(0:2)::Ka,Kc,Kq
    REAL(KIND(1d0)),DIMENSION(4,4)::L,LS
    REAL(KIND(1d0)),DIMENSION(4,4)::SU ! (set 1E 1M 2E 2M, column 1..4) of the partial fractions summed over the poles
    REAL(KIND(1d0)),DIMENSION(NG),SAVE::tg,wg
    REAL(KIND(1d0)),DIMENSION(NC),SAVE::tc,wc
    REAL(KIND(1d0)),DIMENSION(NC,NC),SAVE::TL,TR
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::GM
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::ghc
    INTEGER,DIMENSION(:),ALLOCATABLE::glo,ghi,lmx,loff
    INTEGER,DIMENSION(:,:),ALLOCATABLE::lv
    INTEGER::ngr,lev
    REAL(KIND(1d0)),DIMENSION(2)::ce
    REAL(KIND(1d0)),DIMENSION(NC)::lw
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::G2
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::cn
    INTEGER,DIMENSION(:),ALLOCATABLE::pan
    REAL(KIND(1d0))::hf,hc,ca,cb,u
    INTEGER::ncp
    LOGICAL,SAVE::glinit=.FALSE.
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::e,kn,wk
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::g
    REAL(KIND(1d0))::s,w2,a2,c2,P,Q4,bt,rE,rM,kmax,bcut,bthi,h,J0,J1,J2,xs
    INTEGER::nq,j,iq,t,jlo,jhi,ip,m,ne,nk,ig,k,isave
    dens=0d0
    ! 1-x < 1e-14: the density underflows to 0 (exactly 0 in double precision for 1-x <= 1e-12, fixed and running
    ! alpha); below 3.4e-16 the running-alpha pole sums gave NaN (the end point x = 1 of the luminosity integrals)
    IF(x.GE.1d0-1d-14.OR.x.LE.0d0)RETURN
    xs=xexact; isave=iFFexact
    xexact=x
    s=1d0-x
    w2=(x*mproton)**2
    a2=w2+s*Lambda2
    c2=w2+4d0*s*mproton**2
    P=(s*Lambda2)**2
    Q4=4d0*s*mproton**2
    ! the 1/(y+q) terms: q = w2 (weight 1) and, with running alpha, w2+(1-x) srun(j) (weights crun(j))
    nq=0; qv(0)=w2; wq(0)=1d0
    IF(alpha_running4photon)THEN
       DO j=0,nrun
          nq=nq+1; qv(nq)=w2+s*srun(j); wq(nq)=crun(j)
       ENDDO
    ENDIF
    DO iq=0,nq
       CALL partial_fractions((/c2,qv(iq)/),a2,P*Q4,0d0,CU1E(:,iq))
       CALL partial_fractions((/c2,qv(iq)/),a2,P*w2,P,CU1M(:,iq))
       CALL partial_fractions((/c2,qv(iq)/),a2,-P*Q4,0d0,CU2E(:,iq))
       CALL partial_fractions((/c2,qv(iq)/),a2,P*Q4,0d0,CU2M(:,iq))
    ENDDO
    rE=0d0; rM=0d0
    IF(iFF.EQ.2)THEN
       rE=A1_end_ratio(1)-1d0; rM=mupomuN*(A1_end_ratio(2)-1d0)
    ENDIF
    bcut=bmaxoR*Rproton/GeVm12fm
    ! the terms in 1/(y+c2), 1/(y+a2), 1/(y+a2)**2 summed over the poles q once per row (SU), so that per b only the
    ! 1/(y+q) terms are summed; a pole term is left out where it is below exp(-50) of the q = w2 term
    ! ((sqrt(q)-sqrt(w2)) b~ > 50, as dipole_running_amplitude_row)
    SU(1,:)=MATMUL(CU1E(:,0:nq),wq(0:nq)); SU(2,:)=MATMUL(CU1M(:,0:nq),wq(0:nq))
    SU(3,:)=MATMUL(CU2E(:,0:nq),wq(0:nq)); SU(4,:)=MATMUL(CU2M(:,0:nq),wq(0:nq))
    I=0d0
    DO j=1,n
       bt=b(j)/s
       IF(bt.LE.0d0)CYCLE
       CALL besselK012(DSQRT(a2)*bt,Ka)
       CALL besselK012(DSQRT(c2)*bt,Kc)
       ! L(t,1..4): the integrals of type t with 1/(y+c2), 1/(y+q), 1/(y+a2), 1/(y+a2)**2; column 2 per pole
       DO t=1,4
          L(t,1)=exact_dipole_basis(t,c2,bt,Kc,.FALSE.)
          L(t,3)=exact_dipole_basis(t,a2,bt,Ka,.FALSE.)
          L(t,4)=exact_dipole_basis(t,a2,bt,Ka,.TRUE.)
       ENDDO
       ! (type t, set u) of the six amplitudes: E (u=1 for 1-2, 3 for 3-6), M (u=2, 4)
       DO t=1,4
          DO k=1,4
             LS(t,k)=SU(k,1)*L(t,1)+SU(k,3)*L(t,3)+SU(k,4)*L(t,4)
          ENDDO
       ENDDO
       DO iq=0,nq
          IF((DSQRT(qv(iq))-DSQRT(qv(0)))*bt.GT.50d0)CYCLE
          CALL besselK012(DSQRT(qv(iq))*bt,Kq)
          DO t=1,4
             u=wq(iq)*exact_dipole_basis(t,qv(iq),bt,Kq,.FALSE.)
             LS(t,1)=LS(t,1)+CU1E(2,iq)*u; LS(t,2)=LS(t,2)+CU1M(2,iq)*u
             LS(t,3)=LS(t,3)+CU2E(2,iq)*u; LS(t,4)=LS(t,4)+CU2M(2,iq)*u
          ENDDO
       ENDDO
       AE=(/LS(1,1),LS(2,1),LS(3,3),LS(4,3),LS(1,3),LS(2,3)/)
       AM=(/LS(1,2),LS(2,2),LS(3,4),LS(4,4),LS(1,4),LS(2,4)/)
       I(:,j)=AE+mupomuN*AM
       IF(iFF.EQ.2.AND.b(j).LE.bcut)I(:,j)=I(:,j)+rE*AE+rM*AM
    ENDDO
    ! the correction beyond the dipole: the integrand on fine Gauss-Legendre nodes of the row (panels between the
    ! joins of the form factors, width <= HF: the structure of the fits, and <= PH/b~max), then for each octave of b~
    ! projected on the NC Chebyshev points of coarse panels (phase b~ k_T <= PC across a panel: J_n interpolated to
    ! ~1e-11), so that each b costs NC Bessel evaluations per coarse panel; directly on the fine nodes where the coarse
    ! panels would be narrower than the fine ones
    Ic=0d0
    kmax=0d0
    IF(iFF.NE.0)kmax=DSQRT(MAX(Q2max_corr(iFF)*s-w2,0d0))
    IF(kmax.GT.0d0.AND.b(1).LE.bcut)THEN
       IF(.NOT.glinit)THEN
          CALL gauss_legendre(NG,tg,wg)
          DO m=1,NC
             tc(m)=-DCOS(3.14159265358979323846264338328d0*(m-1)/(NC-1))
          ENDDO
          DO m=1,NC
             wc(m)=1d0
             DO t=1,NC
                IF(t.NE.m)wc(m)=wc(m)/(tc(m)-tc(t))
             ENDDO
          ENDDO
          ! TL(m',m), TR(m',m): parent Lagrange basis m at the point m' of the left / right child
          DO t=1,NC
             DO k=1,2
                u=(tc(t)+MERGE(-1d0,1d0,k.EQ.1))/2d0
                IF(ANY(ABS(u-tc).LT.1d-14))THEN
                   lw=0d0; lw(MINLOC(ABS(u-tc),1))=1d0
                ELSE
                   lw=wc/(u-tc); lw=lw/SUM(lw)
                ENDIF
                IF(k.EQ.1)THEN
                   TL(t,:)=lw
                ELSE
                   TR(t,:)=lw
                ENDIF
             ENDDO
          ENDDO
          glinit=.TRUE.
       ENDIF
       iFFexact=-iFF
       e=exact_edges()
       e=(/0d0,PACK(e,e.LT.kmax),kmax/)
       ne=SIZE(e)
       jhi=COUNT(b.LE.bcut)
       hf=MIN(MERGE(HF2,HF1,iFF.EQ.2),PH/MAX(b(jhi)/s,1d-300))
       nk=0
       DO ip=1,ne-1
          nk=nk+NG*MAX(1,CEILING((e(ip+1)-e(ip))/hf))
       ENDDO
       ALLOCATE(kn(nk),wk(nk),g(6,nk),pan(nk))
       k=0
       DO ip=1,ne-1
          m=MAX(1,CEILING((e(ip+1)-e(ip))/hf))
          h=(e(ip+1)-e(ip))/m
          DO t=1,m
             DO ig=1,NG
                k=k+1
                kn(k)=e(ip)+h*(t-0.5d0+0.5d0*tg(ig))
                wk(k)=0.5d0*h*wg(ig)
                pan(k)=ip
             ENDDO
          ENDDO
       ENDDO
       DO k=1,nk
          CALL exact_integrands(kn(k),6,v)
          g(:,k)=wk(k)*v
       ENDDO
       ! groups of b: octaves of b~ (b <= bcut), each with coarse panels of width <= hc = PC/b~max
       ALLOCATE(glo(n),ghi(n),ghc(n))
       ngr=0
       jhi=0
       DO WHILE(jhi.LT.n)
          jlo=jhi+1
          IF(b(jlo).GT.bcut)EXIT
          bthi=2d0*b(jlo)/s
          jhi=jlo
          DO WHILE(jhi.LT.n)
             IF(b(jhi+1)/s.GT.bthi.OR.b(jhi+1).GT.bcut)EXIT
             jhi=jhi+1
          ENDDO
          ngr=ngr+1
          glo(ngr)=jlo; ghi(ngr)=jhi; ghc(ngr)=PC/(b(jhi)/s)
       ENDDO
       ! nested coarse panels: 2**l equal panels of [0,kmax] at level l (they interpolate J_n(b~ k), smooth across the
       ! joins of the form factors, which only the fine nodes need); a group uses the coarsest level with width <= hc.
       ! The moments of the fine integrand on the Chebyshev points are computed at the finest level (barycentric
       ! weights at the fine nodes) and passed to the coarser levels with the values of the parent Lagrange basis at
       ! the children's points (exact: the parent basis is a polynomial of degree NC-1 on each child)
       ce=(/0d0,kmax/)
       ALLOCATE(lv(1,ngr),lmx(1),loff(1))
       lmx=-1
       DO ig=1,ngr
          DO ip=1,1
             lv(ip,ig)=-1
             IF(ghc(ig).LE.2d0*hf)CYCLE
             lv(ip,ig)=MAX(0,CEILING(DLOG((ce(ip+1)-ce(ip))/ghc(ig))/DLOG(2d0)-1d-12))
             lmx(ip)=MAX(lmx(ip),lv(ip,ig))
          ENDDO
       ENDDO
       ncp=0
       DO ip=1,1
          loff(ip)=ncp
          IF(lmx(ip).GE.0)ncp=ncp+2**(lmx(ip)+1)-1
       ENDDO
       ! panel q (0-based) of level lev of interval ip: column loff(ip)+2**lev-1+q+1 of GM
       ALLOCATE(GM(6,NC,MAX(ncp,1)))
       GM=0d0
       DO k=1,nk
          ip=1
          IF(lmx(ip).LT.0)CYCLE
          m=2**lmx(ip)
          h=(ce(ip+1)-ce(ip))/m
          t=MIN(MAX(INT((kn(k)-ce(ip))/h),0),m-1)
          ca=ce(ip)+t*h
          u=2d0*(kn(k)-ca)/h-1d0
          lw=wc/(u-tc)
          IF(ANY(ABS(u-tc).LT.1d-14))THEN
             lw=0d0; lw(MINLOC(ABS(u-tc),1))=1d0
          ELSE
             lw=lw/SUM(lw)
          ENDIF
          iq=loff(ip)+m-1+t+1
          DO ig=1,NC
             GM(:,ig,iq)=GM(:,ig,iq)+g(:,k)*lw(ig)
          ENDDO
       ENDDO
       DO ip=1,1
          DO lev=lmx(ip)-1,0,-1
             DO t=0,2**lev-1
                iq=loff(ip)+2**lev-1+t+1
                k=loff(ip)+2**(lev+1)-1+2*t+1
                GM(:,:,iq)=MATMUL(GM(:,:,k),TL)+MATMUL(GM(:,:,k+1),TR)
             ENDDO
          ENDDO
       ENDDO
       DO ig=1,ngr
          jlo=glo(ig); jhi=ghi(ig)
          IF(ghc(ig).LE.2d0*hf)THEN
             ! direct sums on the fine nodes
             DO j=jlo,jhi
                bt=b(j)/s
                v=0d0
                DO k=1,nk
                   CALL bj012(bt*kn(k),J0,J1,J2)
                   v=v+g(:,k)*(/J1,J0,J0,J2,J1,J0/)
                ENDDO
                Ic(:,j)=v
             ENDDO
          ELSE
             ! the Chebyshev points of the group's level and their moments
             m=0
             DO ip=1,1
                m=m+NC*2**lv(ip,ig)
             ENDDO
             IF(ALLOCATED(G2))DEALLOCATE(G2,cn)
             ALLOCATE(G2(6,m),cn(m))
             m=0
             DO ip=1,1
                lev=lv(ip,ig)
                h=(ce(ip+1)-ce(ip))/2**lev
                DO t=0,2**lev-1
                   iq=loff(ip)+2**lev-1+t+1
                   DO k=1,NC
                      m=m+1
                      cn(m)=ce(ip)+(t+0.5d0)*h+0.5d0*h*tc(k)
                      G2(:,m)=GM(:,k,iq)
                   ENDDO
                ENDDO
             ENDDO
             DO j=jlo,jhi
                bt=b(j)/s
                v=0d0
                DO k=1,m
                   CALL bj012(bt*cn(k),J0,J1,J2)
                   v=v+G2(:,k)*(/J1,J0,J0,J2,J1,J0/)
                ENDDO
                Ic(:,j)=v
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    DO j=1,n
       IF(b(j).LE.0d0)CYCLE
       dens(j)=exact_combination(I(:,j)+Ic(:,j))/pi2/s**2
    ENDDO
    xexact=xs; iFFexact=isave
  CONTAINS
    SUBROUTINE bj012(z,j0,j1,j2) ! J_0, J_1, J_2 at z >= 0
      REAL(KIND(1d0)),INTENT(IN)::z
      REAL(KIND(1d0)),INTENT(OUT)::j0,j1,j2
      j0=BESSEL_J0(z); j1=BESSEL_J1(z)
      IF(z.GT.1d-4)THEN
         j2=2d0*j1/z-j0
      ELSE
         j2=z*z/8d0*(1d0-z*z/12d0)
      ENDIF
    END SUBROUTINE bj012
  END SUBROUTINE proton_exact_density_row

  SUBROUTINE gauss_legendre(nn,tt,ww)
    ! Gauss-Legendre nodes and weights on [-1,1] (Newton on P_nn)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nn
    REAL(KIND(1d0)),DIMENSION(nn),INTENT(OUT)::tt,ww
    REAL(KIND(1d0))::zz,p1,p2,p3,pp
    INTEGER::ii,jj,it
    DO ii=1,nn
       zz=DCOS(3.14159265358979323846264338328d0*(ii-0.25d0)/(nn+0.5d0))
       DO it=1,100
          p1=1d0; p2=0d0
          DO jj=1,nn
             p3=p2; p2=p1; p1=((2d0*jj-1d0)*zz*p2-(jj-1d0)*p3)/jj
          ENDDO
          pp=nn*(zz*p1-p2)/(zz*zz-1d0)
          zz=zz-p1/pp
          IF(ABS(p1/pp).LT.1d-16)EXIT
       ENDDO
       p1=1d0; p2=0d0
       DO jj=1,nn
          p3=p2; p2=p1; p1=((2d0*jj-1d0)*zz*p2-(jj-1d0)*p3)/jj
       ENDDO
       pp=nn*(zz*p1-p2)/(zz*zz-1d0)
       tt(ii)=-zz; ww(ii)=2d0/((1d0-zz*zz)*pp*pp)
    ENDDO
  END SUBROUTINE gauss_legendre

  SUBROUTINE corr_grid_column(xga,nb,bv,col)
    ! the k_T integral of PND_ChargeFormFactor_proton_corr_fxn (form factor iFFcorr, running alpha as set) at all
    ! b of bv(1:nb) for one x_gamma: composite Gauss-Legendre rule (NG nodes per panel) on panels between the joins of
    ! the form factor (dde1: ratio_dde1; luxqed: knots of the A1 table) of width <= HF (the structure of the fits) and
    ! <= PH/max(b) (phase of J_1 across a panel); the Bessel functions on the shared nodes for all b at once
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::xga
    INTEGER,INTENT(IN)::nb
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(IN)::bv
    REAL(KIND(1d0)),DIMENSION(nb),INTENT(OUT)::col
    REAL(KIND(1d0)),PARAMETER::PH=10d0,HF1=0.25d0,HF2=0.05d0
    INTEGER,PARAMETER::NG=20
    REAL(KIND(1d0)),DIMENSION(4),PARAMETER::Q2dde1=(/1.02d0,1.22d0,2.5d0,3.5d0/)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER::Q2A1=(/0.25d0,0.5d0,0.75d0,1d0,1.5d0,3d0,9.995d0/)
    REAL(KIND(1d0)),DIMENSION(NG),SAVE::tg,wg
    LOGICAL,SAVE::glinit=.FALSE.
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::q2e,e
    REAL(KIND(1d0)),DIMENSION(31)::eg
    REAL(KIND(1d0))::e2,kmax,hf,h,k,Q2,f,dipole,CFF,ek
    REAL(KIND(1d0)),PARAMETER::aa=1.1867816581938533d0 ! 1/DSQRT(0.71 GeV^2), as PND_ChargeFormFactor_proton_corr_fxn
    REAL(KIND(1d0))::rE
    INTEGER::ip,m,t,ig,ne
    col=0d0
    e2=(xga*mproton)**2
    IF(e2.GE.Q2max_corr(iFFcorr).OR.nb.LE.0)RETURN
    kmax=DSQRT(Q2max_corr(iFFcorr)-e2)
    IF(.NOT.glinit)THEN
       CALL gauss_legendre(NG,tg,wg)
       glinit=.TRUE.
    ENDIF
    IF(iFFcorr.EQ.2)THEN
       q2e=Q2A1
    ELSE
       q2e=Q2dde1
    ENDIF
    q2e=q2e-e2
    rE=A1_end_ratio(1)
    e=DSQRT(PACK(q2e,q2e.GT.0d0))
    ! k_T**2/(k_T**2+x**2 mp**2) varies on the scale x mp: edges at x mp 4**m below the first join (and the first A1 knot)
    ek=MIN(MINVAL(e),MERGE(DSQRT(MAX(0.005d0-e2,0d0)),1d300,iFFcorr.EQ.2))
    eg=(/(DSQRT(e2)*4d0**m,m=0,30)/)
    e=(/0d0,PACK(eg,eg.LT.MIN(ek,kmax)),PACK(e,e.LT.kmax),kmax/)
    ne=SIZE(e)
    hf=MIN(MERGE(HF2,HF1,iFFcorr.EQ.2),PH/MAX(MAXVAL(bv),1d-300))
    DO ip=1,ne-1
       m=MAX(1,CEILING((e(ip+1)-e(ip))/hf))
       h=(e(ip+1)-e(ip))/m
       DO t=1,m
          DO ig=1,NG
             k=e(ip)+h*(t-0.5d0+0.5d0*tg(ig))
             Q2=k*k+e2
             dipole=1d0/(1d0+Q2*aa**2)**2
             IF(iFFcorr.EQ.2)THEN
                CFF=formfactor_A1(Q2,1)-rE*dipole
             ELSE
                CFF=dipole*(ratio_dde1(Q2)-1d0)
             ENDIF
             IF(alpha_running4photon)CFF=CFF*alpha_ratio_tab(Q2)
             f=0.5d0*h*wg(ig)*k*k/Q2*CFF
             col=col+f*BESSEL_J1(bv*k)
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE corr_grid_column

  SUBROUTINE besselK012(z,K)
    ! modified Bessel functions K_0, K_1, K_2 at z > 0 to double precision: the Chebyshev expansions of tbessk.f90
    ! (exp-scaled, ~1e-15), K_2 = K_0 + 2 K_1/z; zero above z = 700 (as before)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::z
    REAL(KIND(1d0)),DIMENSION(0:2),INTENT(OUT)::K
    REAL(KIND(1d0)),EXTERNAL::GUPC_BESSK0E,GUPC_BESSK1E
    REAL(KIND(1d0))::ez
    K=0d0
    IF(z.GT.7d2)RETURN
    ez=DEXP(-z)
    K(0)=GUPC_BESSK0E(z)*ez
    K(1)=GUPC_BESSK1E(z)*ez
    K(2)=K(0)+2d0*K(1)/z
    RETURN
  END SUBROUTINE besselK012

  FUNCTION PhotonNumberDensity_ChargeFormFactor_WS(b,Ega,gamma,RR,w,aa,bcut,btilcut,ibeam,integ_method)
    ! It gives us the photon number density with Z=1 and alpha=1 
    ! b,RR,aa should be written in unit of GeV-1
    ! Ega should be in unit of GeV
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    INTEGER::ibeam
    REAL(KIND(1d0))::PhotonNumberDensity_ChargeFormFactor_WS
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma,RR,w,aa
    REAL(KIND(1d0)),INTENT(IN)::bcut ! if bcut > 0, when b > bcut*RR, it will simply use PhotonNumberDensity (not from form factor). 
                                     ! This might be necessary in order to improve the numerical efficiency
                                     ! A nominal value is 2-3.
    REAL(KIND(1d0)),INTENT(IN)::btilcut ! if btilcut > 0, when b*Ega/gamma > btilcut*RR, it will simply use  PhotonNumberDensity (necessary for numerical stability)
                                        ! A nominal value is 0.7d0.
    INTEGER,INTENT(IN),OPTIONAL::integ_method ! 1: direct trapezoid rule; 2: modified W transform with simpson (a bit slow)
    ! do not generate the grid but direct integration when integ_method < 0
    INTEGER::integ_method6
    REAL(KIND(1d0)),PARAMETER::one=1d0
    REAL(KIND(1d0))::Egaoga,integral,xga,Egaoga2,xga2,Ega2
    REAL(KIND(1d0))::Egaoga_common,b_common,R_common,w_common,aa_common
    COMMON/PND_CFF_WS/Egaoga_common,b_common,R_common,w_common,aa_common
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA,NYA,i,j,n,k,l
    SAVE init,NXA,NYA
    INTEGER,DIMENSION(2)::NYA_save=(/0,0/)
    SAVE NYA_save
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA,YA,ZA1_save,ZA2_save
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::ZA
    SAVE XA,YA,ZA,ZA1_save,ZA2_save
    REAL(KIND(1d0)),PARAMETER::bmaxoR=10d0
    INTEGER,PARAMETER::NBMAX=2
    ! From Nicolas Crepet:
    ! For x>0.1, it turns out that we only use NBSEG=28 and NXSEG=16 to be
    ! smooth enough for the curve
    ! 0 to 10**(-nbmax)*bmax
    ! 10**(-n-1)*bmax to 10**(-n)*bmax
    INTEGER,PARAMETER::NBSEG=12
    ! NXSEG for x_gamma from 10**(-n-1) to 10**(-n) 
    INTEGER,PARAMETER::NXSEG=8
    INTEGER::log10xmin,log10xmin_before,ilog10x
    REAL(KIND(1d0))::XMIN=1D-8
    SAVE XMIN
    REAL(KIND(1d0))::log10x1
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0)),DIMENSION(2)::rescaling_bmax_save
    REAL(KIND(1d0))::db,bb
    SAVE rescaling_bmax_save
    REAL(KIND(1d0)),DIMENSION(2)::R_save,w_save,aa_save
    SAVE R_save,w_save,aa_save
    INTEGER,PARAMETER::n_interp=6
    REAL(KIND(1d0)),DIMENSION(n_interp)::XA2,YA2
    REAL(KIND(1d0)),PARAMETER::PIPI=3.14159265358979323846264338328d0
    INTEGER::iter,npoints
    REAL(KIND(1d0))::integ1,integ2,integ3,integ4
    INTEGER,PARAMETER::itermax=12
    REAL(KIND(1d0)),DIMENSION(1)::XI,YI,ZI
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(nxseg+1)::YD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp,nxseg+1)::ZD2
    IF(ibeam.GT.2.OR.ibeam.LT.1)THEN
       WRITE(*,*)"Error: ibeam=/=1,2 in PhotonNumberDensity_ChargeFormFactor_WS"
       STOP 1
    ENDIF
    Egaoga=Ega/gamma ! gamma=Ebeam/mN, xga=Ega/Ebeam, Egaoga=xga*mN
    xga=Egaoga/mN
    IF(xga.LE.1d-99)THEN
       PhotonNumberDensity_ChargeFormFactor_WS=0d0
       RETURN
    ENDIF
    integ_method6=1
    IF(PRESENT(integ_method))THEN
       integ_method6=integ_method
    ENDIF
    IF(init(ibeam).EQ.0.OR.xga.LT.xmin.OR.(integ_method6.LT.0))THEN
       ! first do a rescaling (this can be explicitly verified from the analytic expressions)
       R_common=RR*Egaoga
       w_common=w ! for w I do not need the scaling
       aa_common=aa*Egaoga
       R_save(ibeam)=RR
       w_save(ibeam)=w
       aa_save(ibeam)=aa
       b_common=b*Egaoga
       Egaoga_common=one
       IF(integ_method6.LT.0)THEN
          IF(b.EQ.0d0)THEN
             PhotonNumberDensity_ChargeFormFactor_WS=0d0
             RETURN
          ELSE
             IF(btilcut.GT.0d0.AND.b_common.GT.btilcut*RR)THEN
                PhotonNumberDensity_ChargeFormFactor_WS=PhotonNumberDensity(b,Ega,gamma)
             ELSE
                npoints=30000
                CALL trapezoid_integration(npoints,PND_ChargeFormFactor_WS_fxn,&
                     one,integral)
                integral=integral*Egaoga/2d0
                integ2=PhotonNumberDensity_AnalyticInt4Series_WS(b,Ega,gamma,RR,w,aa,-1,10)/PIPI
                PhotonNumberDensity_ChargeFormFactor_WS=(integral+integ2)**2
             ENDIF
             RETURN
          ENDIF
       ENDIF
       rescaling_bmax_save(ibeam)=MAX(b,bmaxoR*RR)
       IF(bcut.GT.0d0)THEN
          rescaling_bmax_save(ibeam)=MIN(rescaling_bmax_save(ibeam),bcut*RR)
       ENDIF
       NXA=NBSEG*(nbmax+1)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       log10xmin_before=INT(DLOG10(xmin))
       IF(xga.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(xga))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(xga.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
          IF(ALLOCATED(YA))THEN
             WRITE(*,*)" "
          ENDIF
          WRITE(*,*)"INFO: in ChFF photon number density of nucleus, the xmin of the grid has been updated to ",xmin
          IF(ALLOCATED(YA))THEN
             DEALLOCATE(YA)
          ENDIF
          IF(ALLOCATED(ZA))THEN
             ! save the values that have been calculated before
             IF(init(1).EQ.1)THEN
                NYA_save(1)=NYA
                IF(ALLOCATED(ZA1_save))THEN
                   DEALLOCATE(ZA1_save)
                ENDIF
                ALLOCATE(ZA1_save(NXA,NYA))
                DO I=1,NXA
                   DO J=1,NYA
                      ZA1_save(I,J)=ZA(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(init(2).EQ.1)THEN
                NYA_save(2)=NYA
                IF(ALLOCATED(ZA2_save))THEN
                   DEALLOCATE(ZA2_save)
                ENDIF
                ALLOCATE(ZA2_save(NXA,NYA))
                DO I=1,NXA
                   DO J=1,NYA
                      ZA2_save(I,J)=ZA(2,I,J)
                   ENDDO
                ENDDO
             ENDIF
             ! then deallocate it
             DEALLOCATE(ZA)
          ENDIF
          init(1:2)=0
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       NYA=NXSEG*(-log10xmin)+1
       IF(.NOT.ALLOCATED(YA))THEN
          ALLOCATE(YA(2,NYA))
       ENDIF
       IF(.NOT.ALLOCATED(ZA))THEN
          ALLOCATE(ZA(2,NXA,NYA))
       ENDIF
       db=9d0/DBLE(NBSEG)
       IF(init(3-ibeam).EQ.1.AND.rescaling_bmax_save(ibeam).EQ.rescaling_bmax_save(3-ibeam).AND.&
            R_save(ibeam).EQ.R_save(3-ibeam).AND.w_save(ibeam).EQ.w_save(3-ibeam).AND.&
            aa_save(ibeam).EQ.aa_save(3-ibeam))THEN
          DO k=1,NXA
             XA(ibeam,k)=XA(3-ibeam,k)
          ENDDO
          DO k=1,NYA
             YA(ibeam,k)=YA(3-ibeam,k)
          ENDDO
          DO i=1,NXA
             DO j=1,NYA
                ZA(ibeam,i,j)=ZA(3-ibeam,i,j)
             ENDDO
          ENDDO
       ELSE
          WRITE(*,*)"INFO: generate grid of photon number density from charge form factor of beam=",ibeam
          WRITE(*,*)"INFO: it will take a few minutes !"
          k=0
          DO n=0,nbmax
             ! from 10**(-n-1)*bmax to 10**(-n)*bmax
             DO i=1,NBSEG
                k=NBSEG*n+i
                ! these are b in unit GeV-1 (not multiplied Egaoga yet !)
                XA(ibeam,NXA-k+1)=(10d0**(-n-1))*(1d0+DBLE(NBSEG+1-i)*db)*rescaling_bmax_save(ibeam)
             ENDDO
          ENDDO
          IF(k+1.NE.NXA)THEN
             WRITE(*,*)"ERROR: mismatching k+1 and NXA in PhotonNumberDensity_ChargeFormFactor_WS"
             STOP 1
          ENDIF
          XA(ibeam,1)=0d0
          K=0
          DO I=0,log10xmin+1,-1
             DO J=1,nxseg
                log10x1=-1d0/DBLE(nxseg)*DBLE(J-1)+DBLE(I)
                K=K+1
                YA(ibeam,K)=log10x1
             ENDDO
          ENDDO
          IF(K.NE.NYA-1)THEN
             WRITE(*,*)"ERROR: K != NYA-1"
             STOP 1
          ENDIF
          YA(ibeam,NYA)=DBLE(log10xmin)
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam,I).EQ.0d0)THEN
                ! b=0: the whole row vanishes (J is not the x index here)
                ZA(ibeam,I,1:NYA)=0d0
                CYCLE
             ENDIF
             DO J=1,NYA
                ! use the ones that have been calculated before
                IF(NYA_save(ibeam).GT.0.AND.J.LE.NYA_save(ibeam))THEN
                   IF(ibeam.EQ.1)THEN
                      ZA(ibeam,I,J)=ZA1_save(I,J)
                   ELSEIF(ibeam.EQ.2)THEN
                      ZA(ibeam,I,J)=ZA2_save(I,J)
                   ELSE
                      WRITE(*,*)"ERROR: do not know ibeam=",ibeam
                   ENDIF
                   CYCLE
                ENDIF
                xga2=10d0**(YA(ibeam,J)) ! x_gamma
                Egaoga2=xga2*mN
                Ega2=Egaoga2*gamma
                b_common=XA(ibeam,I)*Egaoga2 ! = b*x*mN
                R_common=RR*Egaoga2 ! =  R*x*mN
                aa_common=aa*Egaoga2 ! = aa*x*mN
                IF((btilcut.GT.0d0.AND.b_common.GT.btilcut*RR).OR.(xga2.LT.1D-4.AND.XA(ibeam,I).GT.RR))THEN
                   integ2=PhotonNumberDensity_AnalyticInt4Series_WS(XA(ibeam,I),Ega2,gamma,RR,w,aa,-1,10)/PIPI
                   integ3=PhotonNumberDensity(b_common/Egaoga2,Ega2,gamma)
                   integral=DSQRT(integ3)-integ2
                   ZA(ibeam,I,J)=integral
                ELSE
                   IF(integ_method6.EQ.2)THEN
                      ! use modified W transform (by Sidi) to calculate the integral
                      CALL mWT_integrate_PND_ChargeFormFactor_WS(integral)
                      integral=integral*Egaoga2/PIPI
                   ELSE
                      IF(xga2.LT.1D-7)THEN
                         npoints=500000
                      ELSE
                         npoints=10000
                      ENDIF
                      CALL trapezoid_integration(npoints,PND_ChargeFormFactor_WS_fxn,&
                           one,integral)
                      integral=integral*Egaoga2/2d0
                      integ4=integral
                      integ2=PhotonNumberDensity_AnalyticInt4Series_WS(XA(ibeam,I),Ega2,gamma,RR,w,aa,-1,10)/PIPI
                      integ1=(integral+integ2)**2
                      integ3=PhotonNumberDensity(b_common/Egaoga2,Ega2,gamma)
                      IF(b_common/Egaoga2.GE.2d0*RR.AND.(DABS(integ3/integ1).GT.1.5d0.OR.DABS(integ3/integ1).LT.0.67d0))THEN
                         ! when b_common/Egaoga > 2*RA, the EPA is expected to be good
                         iter=1
                         DO WHILE((DABS(integ3/integ1).GT.1.5d0.OR.DABS(integ3/integ1).LT.0.67d0.OR.&
                              DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)&
                              .AND.iter.LE.itermax)
                            integ4=integral
                            ! increase the points by a factor of 2
                            npoints=npoints*2
                            CALL trapezoid_integration(npoints,PND_ChargeFormFactor_WS_fxn,&
                                 one,integral)
                            integral=integral*Egaoga2/2d0
                            integ1=(integral+integ2)**2
                            iter=iter+1
                         END DO
                         IF(DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)THEN
                            WRITE(*,*)"WARNING: the integral is not stable (b,Ega,gamma,RA,wA,aA)=",&
                                 b_common/Egaoga2,Ega2,gamma,RR,w,aa
                            WRITE(*,*)"WARNING: integral in two iterations #1:",integ4,integral
                         ENDIF
                         IF(DABS(integ3/integ1).GT.1.5d0.OR.DABS(integ3/integ1).LT.0.67d0)THEN
                            WRITE(*,*)"WARNING: the EPA is not good at (b,Ega,gamma,RA,wA,aA)=",&
                                 b_common/Egaoga2,Ega2,gamma,RR,w,aa
                            WRITE(*,*)"WARNING: EPA, non-EPA #1:",integ3,integ1
                         ENDIF
                      ELSEIF(DABS(integ3/integ1).LT.0.67d0.AND.xga2.LT.0.2d0.AND.xga2.GE.1D-7)THEN
                         ! in general, we expect the charge form factor is smaller than EPA when xga is not too close to 1
                         iter=1
                         DO WHILE((DABS(integ3/integ1).LT.0.67d0.OR.&
                              DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)&
                              .AND.iter.LE.itermax)
                            integ4=integral
                            ! increase the points by a factor of 2
                            npoints=npoints*2
                            CALL trapezoid_integration(npoints,PND_ChargeFormFactor_WS_fxn,&
                                 one,integral)
                            integral=integral*Egaoga2/2d0
                            integ1=(integral+integ2)**2
                            iter=iter+1
                         END DO
                         IF(DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)THEN
                            WRITE(*,*)"WARNING: the integral is not stable (b,Ega,gamma,RA,wA,aA)=",&
                                 b_common/Egaoga2,Ega2,gamma,RR,w,aa
                            WRITE(*,*)"WARNING: integral in two iterations #2:",integ4,integral
                         ENDIF
                         IF(DABS(integ3/integ1).LT.0.67d0)THEN
                            WRITE(*,*)"WARNING: the EPA is not good at (b,Ega,gamma,RA,wA,aA)=",&
                                 b_common/Egaoga2,Ega2,gamma,RR,w,aa
                            WRITE(*,*)"WARNING: EPA, non-EPA #2:",integ3,integ1
                         ENDIF
                      ELSEIF(b_common.GT.1D-7.AND.xga2.GE.1D-7)THEN
                         ! we try to do some numerical improvement
                         iter=1
                         DO WHILE((iter.EQ.1.OR.&
                              DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)&
                              .AND.iter.LE.itermax)
                            integ4=integral
                            ! increase the points by a factor of 2
                            npoints=npoints*2
                            CALL trapezoid_integration(npoints,PND_ChargeFormFactor_WS_fxn,&
                                 one,integral)
                            integral=integral*Egaoga2/2d0
                            integ1=(integral+integ2)**2
                            iter=iter+1
                         END DO
                         IF(DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)THEN
                            WRITE(*,*)"WARNING: the integral is not stable (b,Ega,gamma,RA,wA,aA)=",&
                                 b_common/Egaoga2,Ega2,gamma,RR,w,aa
                            WRITE(*,*)"WARNING: integral in two iterations #3:",integ4,integral
                         ENDIF
                      ENDIF
                   ENDIF
                   IF(ISNAN(integral))THEN
                      WRITE(*,*)"ERROR: the integral is not stable (b,Ega,gamma,RA,wA,aA)=",&
                           b_common/Egaoga2,Ega2,gamma,RR,w,aa
                      STOP 1
                   ENDIF
                   ZA(ibeam,I,J)=integral
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       init(ibeam)=1
    ENDIF
    IF(R_save(ibeam).NE.RR.OR.w_save(ibeam).NE.w.OR.aa_save(ibeam).NE.aa)THEN
       WRITE(*,*)"ERROR: RA,wA,aA are not consistent in PhotonNumberDensity_ChargeFormFactor_WS"
       WRITE(*,*)"INFO: ibeam=",ibeam
       WRITE(*,*)"INFO: Saved ones (RA,wA,aA)=",R_save(ibeam),w_save(ibeam),aa_save(ibeam)
       WRITE(*,*)"INFO: New ones (RA,wA,aA)=",RR,w,aa
       STOP 1
    ENDIF
    IF(b.GT.rescaling_bmax_save(ibeam).OR.b.LE.0d0)THEN
       IF(bcut.LE.0d0.OR.b.LE.0d0)THEN
          PhotonNumberDensity_ChargeFormFactor_WS=0d0
       ELSE
          ! we simply use PhotonNumberDensity (EPA)
          PhotonNumberDensity_ChargeFormFactor_WS=PhotonNumberDensity(b,Ega,gamma)
       ENDIF
    ELSEIF((bcut.GT.0d0.AND.b.GT.bcut*RR).OR.(btilcut.GT.0d0.AND.b*Egaoga.GT.btilcut*RR).OR.&
         (xga.LT.1D-4.AND.b.GT.RR))THEN
       ! we simply use PhotonNumberDensity (EPA)
       PhotonNumberDensity_ChargeFormFactor_WS=PhotonNumberDensity(b,Ega,gamma)
    ELSE
       XI(1)=b
       YI(1)=DLOG10(xga)

       db=MIN(b/rescaling_bmax_save(ibeam),1d0)
       N=-FLOOR(DLOG10(db))-1 ! b is in 10**(-n-1) to 10**(-n)
       IF(N.LT.0)THEN
          ! b=rescaling_bmax_save(ibeam)
          K=NXA-n_interp
          !integral=YA(ibeam,NXA)
       ELSE
          ! NXA=NBSEG*(nbmax+1)+1
          IF(N.LT.NBMAX)THEN
             k=NXA-NBSEG*n
          ELSE
             k=NBSEG+1
          ENDIF
          IF(XA(ibeam,k-NBSEG).GT.b)THEN
             WRITE(*,*)"Error: k is not proper #1"
             STOP 1
          ENDIF
          IF(XA(ibeam,k).LT.b)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP 1
          ENDIF
          DO i=k-NBSEG,k
             IF(XA(ibeam,i).GE.b)EXIT
          ENDDO
          IF(i-n_interp/2+2.GE.1.AND.i-n_interp/2+1+n_interp.LE.NXA)THEN
             K=i-n_interp/2+1
             !DO j=1,n_interp
             !   XA2(j)=XA(ibeam,i-n_interp/2+1+j)
             !   YA2(j)=YA(ibeam,i-n_interp/2+1+j)
             !ENDDO
          ELSEIF(i-n_interp/2+2.LT.1)THEN
             K=0
             !DO j=1,n_interp
             !   XA2(j)=XA(ibeam,j)
             !   YA2(j)=YA(ibeam,j)
             !ENDDO
          ELSEIF(i-n_interp/2+1+n_interp.GT.NXA)THEN
             K=NXA-n_interp
             !DO j=1,n_interp
             !   XA2(n_interp-j+1)=XA(ibeam,NA+1-j)
             !   YA2(n_interp-j+1)=YA(ibeam,NA+1-j)
             !ENDDO
          ELSE
             WRITE(*,*)"Error: you cannot reach here !"
             STOP 1
          ENDIF
          !CALL SPLINE_INTERPOLATE(XA2,YA2,n_interp,bb,integral)
       ENDIF

       IF(YI(1).GE.0d0)THEN
          ilog10x=-1
       ELSE
          ilog10x=FLOOR(YI(1))
       ENDIF
       L=NXSEG*(-ilog10x-1)
       
       DO I=1,n_interp
          XD2_1D(I)=XA(ibeam,K+I)
       ENDDO
       DO I=1,NXSEG+1
          YD2_1D(I)=YA(ibeam,L+I)
       ENDDO
       DO I=1,n_interp
          DO J=1,NXSEG+1
             ZD2(I,J)=ZA(ibeam,K+I,L+J)
          ENDDO
       ENDDO
       CALL lagrange_interp_2d(n_interp-1,NXSEG,XD2_1D,YD2_1D,ZD2,1,XI,YI,ZI)
       integral=ZI(1)
       ! integral=integral*Egaoga/2d0 (the factor has been included in the grid)
       ! the series ones are fully known analytically
       ! Let us keep the first 10 terms
       ! If we want to keep all order terms for K1(btil), set NMIN > 0
       ! otherwise, set NMIN < 0
       integral=integral+PhotonNumberDensity_AnalyticInt4Series_WS(b,Ega,gamma,RR,w,aa,-1,10)/PIPI
       PhotonNumberDensity_ChargeFormFactor_WS=integral**2
    ENDIF
    !PhotonNumberDensity_ChargeFormFactor_WS=Egaoga**2/4d0*integral**2
    RETURN
  END FUNCTION PhotonNumberDensity_ChargeFormFactor_WS

  FUNCTION PND_ChargeFormFactor_WS_fxn(x)
    ! x = ArcTan[kT*gamma/Ega]*2/Pi
    IMPLICIT NONE
    REAL(KIND(1d0))::PND_ChargeFormFactor_WS_fxn
    REAL(KIND(1d0)),INTENT(IN)::x
    REAL(KIND(1d0))::Q,CFF,pref,bkT,kT
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    !REAL(KIND(1d0)),PARAMETER::PIo4=0.785398163397448309615660845820d0
    !REAL(KIND(1d0)),PARAMETER::sqrt2Pi=2.50662827463100050241576528481d0
    REAL(KIND(1d0)),EXTERNAL::BESSJ1
    REAL(KIND(1d0))::Egaoga_common,b_common,R_common,w_common,aa_common
    COMMON/PND_CFF_WS/Egaoga_common,b_common,R_common,w_common,aa_common
    IF(x.GE.1d0.OR.x.LE.0d0)THEN
       PND_ChargeFormFactor_WS_fxn=0d0
       RETURN
    ENDIF
    pref=DTAN(PIo2*x)
    kT=pref*Egaoga_common
    bkT=kT*b_common
    pref=pref**2
    Q=DSQRT(kT**2+Egaoga_common**2)
    ! 10 means including the series 10 terms
    !CFF=ChargeFormFactor_WoodsSaxon(Q,R_common,w_common,aa_common,10)
    ! Let us exclude the series terms, which can be integrated fully analytically
    CFF=ChargeFormFactor_WoodsSaxon(Q,R_common,w_common,aa_common,0)
    IF(ISNAN(CFF))THEN
       PRINT *, "ChargeFormFactor is NaN with ",Q, R_common, w_common, aa_common
       STOP 1
       CFF=0d0
    ENDIF
    !IF(bkT.LE.5d2)THEN
    PND_ChargeFormFactor_WS_fxn=pref*BESSJ1(bkT)*CFF
    !ELSE
    !   PND_ChargeFormFactor_WS_fxn=pref*CFF*(3d0/4d0/DSQRT(bkT)**3/sqrt2Pi*DSIN(PIo4+bkT)&
    !        -2d0/sqrt2Pi/DSQRT(bkT)*DCOS(PIo4+bkT))
    !ENDIF
    RETURN
  END FUNCTION PND_ChargeFormFactor_WS_fxn
  
  SUBROUTINE mWT_integrate_PND_ChargeFormFactor_WS(integral)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::integral
    REAL(KIND(1d0)),EXTERNAL::ZEROJP ! zeros of the derivative of J_n
    REAL(KIND(1d0))::Egaoga_common,b_common,R_common,w_common,aa_common
    COMMON/PND_CFF_WS/Egaoga_common,b_common,R_common,w_common,aa_common
    REAL(KIND(1d0))::btil,Rtil,atil,kTtil,integ,kTtilmax
    REAL(KIND(1d0)),PARAMETER::PIPI=3.14159265358979323846264338328d0
    INTEGER::kmin,kmax,nmax,kk,nn,i,ninterval
    INTEGER::pmax,pmax_real
    INTEGER::pmax_save=-2
    SAVE pmax_save
    INTEGER,PARAMETER::PMAXMAX=20,PMINMIN=15
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XS,PSIS,FS
    !REAL(KIND(1d0)),DIMENSION(-1:PMAXMAX+2)::XS
    !REAL(KIND(1d0)),DIMENSION(-1:PMAXMAX+1)::PSIS,FS
    SAVE XS, PSIS, FS
    REAL(KIND(1d0))::xk,xn,Mp0,Np0
    btil=b_common*Egaoga_common
    Rtil=R_common*Egaoga_common
    atil=aa_common*Egaoga_common
    ! the integrand (charge form factor) is exponentially suppressed via e^(-pi*q*aA)
    ! we stop at e^(-15)
    kTtilmax=DSQRT((15d0/(PIPI*atil))**2-1d0)
    ! the zeros of sin(q*RA)
    ! they are DSQRT(((k*Pi)/Rtil)**2-1d0)
    kmin=MAX(CEILING(Rtil/PIPI),1)
    kmax=FLOOR(DSQRT(kTtilmax**2+1d0)*Rtil/PIPI)
    ! the zeros of J1(kT*b)
    ! the ith zeros of J1(x) is ZEROJP(0,i) with ZEROJP(0,1)=0
    ! The fact is that ZEROJP(0,i) is close to (i-1)*Pi
    ! (i-1)*Pi <= ZEROJP(0,i) < i*Pi
    nmax=FLOOR(kTtilmax/(PIPI*btil))+1
    IF(kmax.GE.kmin)THEN
       pmax=nmax+(kmax-kmin+1)-4
    ELSE
       pmax=nmax-4
    ENDIF
    pmax=MIN(pmax,PMAXMAX)
    pmax=MAX(pmax,PMINMIN)
    IF(pmax.LT.-1)THEN
       ninterval=1000
       ! integrate over f from [a,b] with n intervals
       CALL simpson(PND_ChargeFormFactor_WS_fxn2,0d0,kTtilmax,integral,ninterval)
    ELSE
       IF(pmax_save.LT.pmax)THEN
          ! first let us allocate the arrays
          IF(ALLOCATED(XS))THEN
             DEALLOCATE(XS)
          ENDIF
          IF(ALLOCATED(PSIS))THEN
             DEALLOCATE(PSIS)
          ENDIF
          IF(ALLOCATED(FS))THEN
             DEALLOCATE(FS)
          ENDIF
          ALLOCATE(XS(-1:pmax+2))
          ALLOCATE(PSIS(-1:pmax+1))
          ALLOCATE(FS(-1:pmax+1))
          pmax_save=pmax
       ENDIF
       kk=kmin
       nn=1
       xk=DSQRT(((kk*PIPI)/Rtil)**2-1d0)
       xn=ZEROJP(0,nn)
       pmax_real=pmax
       IF(kmin.LE.kmax)THEN
          DO i=-1,pmax+2
             IF(xn.LT.xk)THEN
                XS(i)=xn
                nn=nn+1
                xn=ZEROJP(0,nn)
             ELSEIF(xk.LT.xn)THEN
                XS(i)=xk
                kk=kk+1
                xk=DSQRT(((kk*PIPI)/Rtil)**2-1d0)
             ELSE
                XS(i)=xn
                nn=nn+1
                xn=ZEROJP(0,nn)
                kk=kk+1
                xk=DSQRT(((kk*PIPI)/Rtil)**2-1d0)
                pmax_real=pmax_real-1
             ENDIF
          ENDDO
       ELSE
          DO i=-1,pmax+2
             XS(i)=ZEROJP(0,i+2)
          ENDDO
       ENDIF
       DO i=-1,pmax_real+1
          ninterval=200
          CALL simpson(PND_ChargeFormFactor_WS_fxn2,XS(i),XS(i+1),PSIS(i),ninterval)
          IF(i.EQ.-1)THEN
             FS(i)=0d0
          ELSE
             FS(i)=FS(i-1)+PSIS(i-1)
          ENDIF
       ENDDO
       Mp0=mWT_Mfun(pmax_real,0,pmax_real+2,XS(-1:pmax_real+2),PSIS(-1:pmax_real+1),&
            FS(-1:pmax_real+1))
       Np0=mWT_Nfun(pmax_real,0,pmax_real+2,XS(-1:pmax_real+2),PSIS(-1:pmax_real+1),&
            FS(-1:pmax_real+1))
       integral=Mp0/Np0
    ENDIF
    RETURN
  END SUBROUTINE mWT_integrate_PND_ChargeFormFactor_WS

  FUNCTION PND_ChargeFormFactor_WS_fxn2(x)
    ! x = kTtil
    IMPLICIT NONE
    REAL(KIND(1d0))::PND_ChargeFormFactor_WS_fxn2
    REAL(KIND(1d0)),INTENT(IN)::x
    REAL(KIND(1d0))::kTtil,Qtil,CFF,bkT,Rtil,atil,btil
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    REAL(KIND(1d0)),EXTERNAL::BESSJ1
    REAL(KIND(1d0))::Egaoga_common,b_common,R_common,w_common,aa_common
    COMMON/PND_CFF_WS/Egaoga_common,b_common,R_common,w_common,aa_common
    REAL(KIND(1d0)),PARAMETER::one=1d0
    IF(x.LE.0d0)THEN
       PND_ChargeFormFactor_WS_fxn2=0d0
       RETURN
    ENDIF
    kTtil=x ! kTtil=kT/Egaoga
    btil=b_common*Egaoga_common
    Rtil=R_common*Egaoga_common
    atil=aa_common*Egaoga_common
    bkT=kTtil*btil
    Qtil=DSQRT(kTtil**2+one)
    ! 10 means including the series 10 terms
    !CFF=ChargeFormFactor_WoodsSaxon(Q,R_common,w_common,aa_common,10)
    ! Let us exclude the series terms, which can be integrated fully analytically
    ! This is rescaling invariant by Egamma/gamma=x_gamma*mN 
    CFF=ChargeFormFactor_WoodsSaxon(Qtil,Rtil,w_common,atil,0)
    IF(ISNAN(CFF))THEN
       PRINT *, "ChargeFormFactor is NaN with ",Qtil, Rtil, w_common, atil
       STOP 1
       CFF=0d0
    ENDIF
    PND_ChargeFormFactor_WS_fxn2=BESSJ1(bkT)*CFF*kTtil**2/Qtil**2
    RETURN
  END FUNCTION PND_ChargeFormFactor_WS_fxn2

  ! Eq.(6) in 1607.06083 with Z=1 and alpha=1
  ! also see my notes OpticalGlauber.pdf
  FUNCTION NGAMMA(xi,gamma)
    IMPLICIT NONE
    REAL(KIND(1d0))::NGAMMA
    REAL(KIND(1d0)),INTENT(IN)::xi,gamma
    REAL(KIND(1d0)),EXTERNAL::BESSK1,BESSK0
    REAL(KIND(1d0)),PARAMETER::PIo2=1.57079632679489661923132169164d0
    NGAMMA=1d0/PIo2*(xi*BESSK0(xi)*BESSK1(xi)&
         -(1d0-1d0/gamma**2)*xi**2/2d0*(BESSK1(xi)**2-BESSK0(xi)**2))
    RETURN
  END FUNCTION NGAMMA

  ! alpha=1
  ! see my notes OpticalGlauber.pdf
  FUNCTION NGAMMA_proton_dipole(aptilde)
    IMPLICIT NONE
    REAL(KIND(1d0))::NGAMMA_proton_dipole
    REAL(KIND(1d0)),INTENT(IN)::aptilde ! aptilde=ap*x*mproton, ap=1/sqrt(0.71) GeV-1
    REAL(KIND(1d0)),PARAMETER::PIPI=3.14159265358979323846264338328d0
    NGAMMA_proton_dipole=1d0/PIPI*((1d0+4d0*aptilde**2)*DLOG(1d0+aptilde**(-2))&
         -(24d0*aptilde**4+42d0*aptilde**2+17d0)/(6d0*(1d0+aptilde**2)**2))
    RETURN
  END FUNCTION NGAMMA_PROTON_DIPOLE

  FUNCTION PhotonFlux_proton_EDFF_nob(x,gamma)
    ! set PNOHARD=1
    ! for proton with EDFF
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_EDFF_nob
    REAL(KIND(1d0)),INTENT(IN)::x,gamma
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    REAL(KIND(1d0))::Rp,alpha,Z
    SAVE Rp,alpha,Z
    REAL(KIND(1d0))::xi
    IF(init.EQ.0)THEN
       Rp=Rproton/GeVm12fm ! from fm to GeV-1
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! for proton
       Z=1d0
       init=1
    ENDIF
    xi=mproton*x*Rp
    PhotonFlux_proton_EDFF_nob=alpha/x*Z**2*NGAMMA(xi,gamma)
    RETURN
  END FUNCTION PhotonFlux_proton_EDFF_nob

  FUNCTION PhotonFlux_proton_dipole_nob(x)
    ! set PNOHARD=1 but with dipole ChFF
    ! for proton
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_dipole_nob
    REAL(KIND(1d0)),INTENT(IN)::x
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0)),PARAMETER::aa=1.1867816581938533d0 ! in unit of GeV-1 = 1/DSQRT(0.71 GeV^2)
    REAL(KIND(1d0))::alpha
    SAVE alpha
    REAL(KIND(1d0))::aptilde
    IF(init.EQ.0)THEN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       ! for proton
       init=1
    ENDIF
    aptilde=mproton*x*aa
    PhotonFlux_proton_dipole_nob=alpha/x*NGAMMA_proton_dipole(aptilde)
    RETURN
  END FUNCTION PhotonFlux_proton_dipole_nob

  FUNCTION PhotonFlux_proton_ChFF_nob(x,gamma,CORRECTIONQ,MU2)
    ! set PNOHARD=1
    ! for proton with ChFF
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_ChFF_nob
    REAL(KIND(1d0)),INTENT(IN)::x,gamma
    INTEGER,INTENT(IN),OPTIONAL::CORRECTIONQ ! the correction beyond the dipole: 0 none (default),
                                             ! 1 dde1, 2 luxqed (as in PhotonNumberDensity_ChargeFormFactor_proton)
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::MU2 ! scale (GeV**2) of the photon normalisation, see alpha_norm_mu
    LOGICAL::corrQ
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::alpha
    SAVE alpha
    REAL(KIND(1d0))::aptilde,Egaoga,runcol(1)
    IF(init.EQ.0)THEN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       init=1
    ENDIF
    IF(PRESENT(CORRECTIONQ))THEN
       corrQ=CORRECTIONQ.GT.0
       IF(corrQ)iFFcorr=CORRECTIONQ
    ELSE
       corrQ=.FALSE.
    ENDIF
    IF(corrQ)THEN
       ! include the correction beyond the dipole
       PhotonFlux_proton_ChFF_nob=PhotonFlux_proton_ChFF_corr_nob_eval(x,gamma)
       PhotonFlux_proton_ChFF_nob=alpha/x*PhotonFlux_proton_ChFF_nob
       PhotonFlux_proton_ChFF_nob=PhotonFlux_proton_ChFF_nob+PhotonFlux_proton_dipole_nob(x)
    ELSE
       PhotonFlux_proton_ChFF_nob=PhotonFlux_proton_dipole_nob(x)
    ENDIF
    IF(alpha_running4photon)THEN
       ! the dipole with the running of alpha, by Parseval: alpha/x*2/pi*int dkT/kT (kT**2/Q2 G_D alpha(Q2)/alpha(0))**2
       xexact=x
       CALL bessel_integrate(dipole_running_collinear_integrand,1,(/-1/),0d0,x*mproton,HUGE(1d0),(/0d0/),runcol)
       PhotonFlux_proton_ChFF_nob=PhotonFlux_proton_ChFF_nob-PhotonFlux_proton_dipole_nob(x)+&
            alpha/x*2d0/3.14159265358979323846264338328d0*runcol(1)
    ENDIF
    PhotonFlux_proton_ChFF_nob=PhotonFlux_proton_ChFF_nob*alpha_norm_mu(MU2)
    RETURN
  END FUNCTION PhotonFlux_proton_ChFF_nob

  FUNCTION PhotonFlux_proton_ChFF_corr_nob_eval(x,gamma)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_ChFF_corr_nob_eval
    REAL(KIND(1d0)),INTENT(IN)::x,gamma
    REAL(KIND(1d0))::Ega
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::npoints
    REAL(KIND(1d0))::bupper,Egaoga,xx1
    REAL(KIND(1d0))::Ega_common,gamma_common
    COMMON /proton_ChFF_corr_nob_fxn/ Ega_common,gamma_common
    INTEGER::init=0
    SAVE init
    INTEGER::MX_save=0
    SAVE MX_save
    INTEGER::iFF_grid=0 ! the form factor and running of the grid (corr_grid_key)
    SAVE iFF_grid
    REAL(KIND(1d0))::xmin=1D-8
    SAVE xmin
    INTEGER::log10xmin,log10xmin_before
    REAL(KIND(1d0))::log10x
    INTEGER::ilog10x
    ! nseg for 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::nseg=10
    INTEGER::MX,I,J,K,L
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::XD_1D
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::ZD,ZD_save
    SAVE MX,XD_1D,ZD,ZD_save
    REAL(KIND(1d0)),DIMENSION(nseg+1)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(nseg+1)::ZD2
    REAL(KIND(1d0)),DIMENSION(1)::XI,ZI
    REAL(KIND(1d0)),PARAMETER::bmaxoR=10d0
    REAL(KIND(1d0)),PARAMETER::Rproton=0.877d0 ! the charge radius of proton (in fm)
    IF(x.LE.1d-99.OR.gamma.LE.0d0)THEN
       PhotonFlux_proton_ChFF_corr_nob_eval=0d0
       RETURN
    ENDIF
    IF(init.EQ.1.AND.corr_grid_key().NE.iFF_grid)THEN
       ! another form factor or running of alpha: generate the grid again
       init=0
       MX_save=0
    ENDIF
    IF(init.EQ.0.OR.x.LT.xmin)THEN
       log10xmin_before=INT(DLOG10(xmin))
       IF(x.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       ! let us generate a grid first
       MX=nseg*(-log10xmin)
       ! try to deallocate first
       IF(ALLOCATED(XD_1D))THEN
          DEALLOCATE(XD_1D)
       ENDIF
       ALLOCATE(XD_1D(MX+1))
       IF(ALLOCATED(ZD))THEN
          DEALLOCATE(ZD)
       ENDIF
       ALLOCATE(ZD(MX+1))
       K=0
       DO I=0,log10xmin+1,-1
          DO J=1,nseg
             log10x=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP 1
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       bupper=bmaxoR*Rproton/GeVm12fm ! the hard cut in the grid of PhotonNumberDensity_ChargeFormFactor_proton_correction 
       DO I=1,MX+1
          IF(MX_save.GT.0.AND.I.LE.MX_save+1)THEN
             ZD(I)=ZD_save(I)
             CYCLE
          ENDIF
          xx1=10d0**(XD_1D(I))
          Ega=xx1*gamma*mproton
          npoints=5000
          Ega_common=Ega
          gamma_common=gamma
          !IF(Ega/gamma.LT.1d-6.AND.MOD(I,3).NE.0)THEN
          !IF(Ega/gamma.LT.1d-8)THEN
          !   ZD(I)=ZD(I-1)
          !ELSE
          CALL trapezoid_integration(npoints,PhotonFlux_proton_ChFF_corr_nob_fxn,&
               bupper,ZD(I))
          IF(MX_save.EQ.0)THEN
             IF(I.EQ.1)THEN
                WRITE(*,*)"INFO: generate grid of proton photon flux beyond dipole correction (will take a few seconds)"
                CALL progress(INT(I*50d0/(MX+1)),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/(MX+1)),50)
             ENDIF
          ENDIF
       ENDDO
       IF(ALLOCATED(ZD_save))THEN
          DEALLOCATE(ZD_save)
       ENDIF
       ALLOCATE(ZD_save(MX+1))
       ! save the values that have been calculated before
       MX_save=MX
       DO I=1,MX+1
          ZD_save(I)=ZD(I)
       ENDDO
       init=1
       iFF_grid=corr_grid_key()
    ENDIF
    XI(1)=DLOG10(x)
    IF(XI(1).GE.0d0)THEN
       ilog10x=-1
    ELSE
       ilog10x=FLOOR(XI(1))
    ENDIF
    K=nseg*(-ilog10x-1)
    ! it turns out to be important with spline
    ! to reverse the order for a fast dropping curve
    ! x is now always from smaller to larger values
    DO I=1,nseg+1
       XD2_1D(nseg+2-I)=XD_1D(K+I)
       ZD2(nseg+2-I)=ZD(K+I)
    ENDDO
    CALL SPLINE_INTERPOLATE(XD2_1D,ZD2,nseg+1,XI(1),PhotonFlux_proton_ChFF_corr_nob_eval)
    RETURN
  END FUNCTION PhotonFlux_proton_ChFF_corr_nob_eval

  FUNCTION PhotonFlux_proton_ChFF_corr_nob_fxn(b)
    ! b is in unit of GeV-1
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_ChFF_corr_nob_fxn
    REAL(KIND(1d0)),INTENT(IN)::b
    REAL(KIND(1d0))::Ega, gamma
    COMMON /proton_ChFF_corr_nob_fxn/ Ega, gamma
    REAL(KIND(1d0)),PARAMETER::twopi=6.2831853071795864769252867d0
    REAL(KIND(1d0))::Egaoga
    IF(b.LE.0d0)THEN
       PhotonFlux_proton_ChFF_corr_nob_fxn=0d0
       RETURN
    ENDIF
    PhotonFlux_proton_ChFF_corr_nob_fxn=PhotonNumberDensity_ChargeFormFactor_proton(b,Ega,gamma,iFFcorr,.TRUE.)
    PhotonFlux_proton_ChFF_corr_nob_fxn=PhotonFlux_proton_ChFF_corr_nob_fxn*twopi*b
    RETURN
  END FUNCTION PhotonFlux_proton_ChFF_corr_nob_fxn

  FUNCTION PhotonFlux_nucleus_EDFF_nob(x,gamma,Z,RA)
    ! set PNOHARD=1
    ! for nucleus with EDFF
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_nucleus_EDFF_nob
    REAL(KIND(1d0)),INTENT(IN)::x,gamma,Z,RA ! RA is the radius of nucleus (in unit of fm)
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    REAL(KIND(1d0))::alpha
    SAVE alpha
    REAL(KIND(1d0))::xi
    IF(init.EQ.0)THEN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
       init=1
    ENDIF
    ! 0.1973 is from fm to GeV-1
    xi=mN*x*RA/GeVm12fm
    PhotonFlux_nucleus_EDFF_nob=alpha/x*Z**2*NGAMMA(xi,gamma)
    RETURN
  END FUNCTION PhotonFlux_nucleus_EDFF_nob

  FUNCTION PhotonFlux_nucleus_ChFF_nob(ibeam,x,gamma,Z,RA,wA,aA)
    ! set PNOHARD=1
    ! for nucleus with ChFF
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_nucleus_ChFF_nob
    INTEGER,INTENT(IN)::ibeam
    REAL(KIND(1d0)),INTENT(IN)::x,gamma,Z
    REAL(KIND(1d0)),INTENT(IN)::RA,wA,aA ! RA is the radius of nucleus (in unit of fm)
                                         ! aA is the diffusivity of nucleus (in unit of fm)
    INTEGER,DIMENSION(2)::init=(/0,0/)
    INTEGER::NXA
    SAVE init,NXA
    INTEGER,DIMENSION(2)::NXA_save=(/0,0/)
    SAVE NXA_save
    REAL(KIND(1d0)),PARAMETER::mN=0.9315d0 ! average nucleaon mass in nuclei (GeV)
    ! NXSEG for x_gamma from 10**(-n-1) to 10**(-n)
    INTEGER,PARAMETER::NXSEG=20
    INTEGER::log10xmin,log10xmin_before,ilog10x,i,j,k,l,n
    REAL(KIND(1d0))::XMIN=1D-8
    SAVE XMIN
    REAL(KIND(1d0))::log10x1
    REAL(KIND(1d0)),DIMENSION(2)::R_save,w_save,aa_save,gamma_save
    SAVE R_save,w_save,aa_save,gamma_save
    REAL(KIND(1d0)),DIMENSION(1)::XI,ZI
    INTEGER,PARAMETER::n_interp=10
    REAL(KIND(1d0)),DIMENSION(n_interp)::XD2_1D
    REAL(KIND(1d0)),DIMENSION(n_interp)::ZD2
    REAL(KIND(1d0))::alpha
    SAVE alpha
    REAL(KIND(1d0))::Ega,xga2,Ega2,Egaoga2
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::XA
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ZA
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::ZA1_save,ZA2_save
    SAVE XA,ZA,ZA1_save,ZA2_save
    IF(ibeam.GT.2.OR.ibeam.LT.1)THEN
       WRITE(*,*)"Error: ibeam=/=1,2 in PhotonFlux_nucleus_ChFF_nob"
       STOP 1
    ENDIF
    IF(x.LT.1d-99)THEN
       PhotonFlux_nucleus_ChFF_nob=0d0
       RETURN
    ENDIF
    IF(init(1).EQ.0.AND.init(2).EQ.0)THEN
       IF(alphaem_elasticphoton.LT.0d0)THEN
          IF(aqedup.GT.0d0)THEN
             alpha=aqedup
          ELSE
             alpha = 0.0072992701d0
          ENDIF
       ELSE
          alpha=alphaem_elasticphoton
       ENDIF
    ENDIF
    IF(init(ibeam).EQ.0.OR.x.LT.xmin)THEN
       R_save(ibeam)=RA
       w_save(ibeam)=wA
       aa_save(ibeam)=aA
       gamma_save(ibeam)=gamma
       log10xmin_before=INT(DLOG10(xmin))
       IF(x.LT.xmin)THEN
          log10xmin=FLOOR(DLOG10(x))
          xmin=10d0**(log10xmin)
          ! to avoid some numerical artificial
          DO WHILE(x.LT.xmin)
             log10xmin=log10xmin-1
             xmin=10d0**(log10xmin)
          ENDDO
          !IF(ALLOCATED(XA))THEN
          !   WRITE(*,*)" "
          !ENDIF
          WRITE(*,*)"INFO: in PhotonFlux_nucleus_ChFF_nob, the xmin of the grid has been updated to ",xmin
          IF(ALLOCATED(XA))THEN
             DEALLOCATE(XA)
          ENDIF
          IF(ALLOCATED(ZA))THEN
             ! save the values that have been calculated before
             IF(init(1).EQ.1)THEN
                NXA_save(1)=NXA
                IF(ALLOCATED(ZA1_save))THEN
                   DEALLOCATE(ZA1_save)
                ENDIF
                ALLOCATE(ZA1_save(NXA))
                DO I=1,NXA
                   ZA1_save(I)=ZA(1,I)
                ENDDO
             ENDIF
             IF(init(2).EQ.1)THEN
                NXA_save(2)=NXA
                IF(ALLOCATED(ZA2_save))THEN
                   DEALLOCATE(ZA2_save)
                ENDIF
                ALLOCATE(ZA2_save(NXA))
                DO I=1,NXA
                   ZA2_save(I)=ZA(2,I)
                ENDDO
             ENDIF
             ! then deallocate it
             DEALLOCATE(ZA)
          ENDIF
          init(1:2)=0
       ENDIF
       log10xmin=INT(DLOG10(xmin))
       NXA=NXSEG*(-log10xmin)+1
       IF(.NOT.ALLOCATED(XA))THEN
          ALLOCATE(XA(2,NXA))
       ENDIF
       IF(.NOT.ALLOCATED(ZA))THEN
          ALLOCATE(ZA(2,NXA))
       ENDIF
       IF(init(3-ibeam).EQ.1.AND.&
            R_save(ibeam).EQ.R_save(3-ibeam).AND.&
            w_save(ibeam).EQ.w_save(3-ibeam).AND.&
            aa_save(ibeam).EQ.aa_save(3-ibeam).AND.&
            gamma_save(ibeam).EQ.gamma_save(3-ibeam))THEN
          DO k=1,NXA
             XA(ibeam,k)=XA(3-ibeam,k)
             ZA(ibeam,k)=ZA(3-ibeam,k)
          ENDDO
       ELSE
          K=0
          DO I=0,log10xmin+1,-1
             DO J=1,nxseg
                log10x1=-1d0/DBLE(nxseg)*DBLE(J-1)+DBLE(I)
                K=K+1
                XA(ibeam,K)=log10x1
             ENDDO
          ENDDO
          IF(K.NE.NXA-1)THEN
             WRITE(*,*)"ERROR: K != NXA-1"
             STOP 1
          ENDIF
          XA(ibeam,NXA)=DBLE(log10xmin)
          DO I=1,NXA
             ! use the ones that have been calculated before
             IF(NXA_save(ibeam).GT.0.AND.I.LE.NXA_save(ibeam))THEN
                IF(ibeam.EQ.1)THEN
                   ZA(ibeam,I)=ZA1_save(I)
                ELSEIF(ibeam.EQ.2)THEN
                   ZA(ibeam,I)=ZA2_save(I)
                ELSE
                   WRITE(*,*)"ERROR: do not known ibeam=",ibeam
                ENDIF
                CYCLE
             ENDIF
             xga2=10d0**(XA(ibeam,I)) ! x_gamma
             Egaoga2=xga2*mN
             Ega2=Egaoga2*gamma
             ZA(ibeam,I)=alpha*Z**2*&
                  PhotonFlux_nucleus_ChFF_nob_eval(ibeam,Ega2,gamma,RA,wA,aA)
             IF(NXA_save(ibeam).EQ.0)THEN
                IF(I.EQ.1)THEN
                   WRITE(*,*)"INFO: generate grid of b-integrated photon flux from charge form factor of beam=",ibeam
                   WRITE(*,*)"INFO: it will take tens of seconds !"
                   CALL progress(INT(I*50d0/NXA),50,.TRUE.)
                ELSE
                   CALL progress(INT(I*50d0/NXA),50)
                ENDIF
             ENDIF
          ENDDO
       ENDIF
       init(ibeam)=1
    ENDIF
    IF(R_save(ibeam).NE.RA.OR.w_save(ibeam).NE.wA.OR.aa_save(ibeam).NE.aA.OR.gamma_save(ibeam).NE.gamma)THEN
       WRITE(*,*)"ERROR: gamma,RA,wA,aA are not consistent in PhotonFlux_nucleus_ChFF_nob"
       WRITE(*,*)"INFO: ibeam=",ibeam
       WRITE(*,*)"INFO: Saved ones (gamma,RA,wA,aA)=",gamma_save(ibeam),R_save(ibeam),w_save(ibeam),aa_save(ibeam)
       WRITE(*,*)"INFO: New ones (gamma,RA,wA,aA)=",gamma,RA,wA,aA
       STOP 1
    ENDIF
    XI(1)=DLOG10(x)
    N=-FLOOR(XI(1))-1 ! x is from 10**(-n-1) to 10**(-n)
    IF(N.LT.0)THEN
       K=0
    ELSE
       IF(N.GE.-INT(DLOG10(xmin)))THEN
          k=NXA-NXSEG
       ELSE
          k=NXSEG*n+1
       ENDIF
       IF(XA(ibeam,k+NXSEG).GT.XI(1))THEN
          WRITE(*,*)"Error: k is not proper #1"
          STOP 1
       ENDIF
       IF(XA(ibeam,k).LT.XI(1))THEN
          WRITE(*,*)"Error: k is not proper #2"
          STOP 1
       ENDIF
       DO i=k,k+NXSEG
          IF(XA(ibeam,i).LE.XI(1))EXIT
       ENDDO
       IF(i-n_interp/2.GE.1.AND.i-n_interp/2-1+n_interp.LE.NXA)THEN
          K=i-n_interp/2-1
       ELSEIF(i-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(i-n_interp/2-1+n_interp.GT.NXA)THEN
          K=NXA-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here !"
          STOP 1
       ENDIF
    ENDIF
    ! it turns out to be important with spline
    ! to reverse the order for a fast dropping curve
    ! x is now always from smaller to larger values
    DO I=1,n_interp
       XD2_1D(n_interp+1-I)=XA(ibeam,K+I)
       ZD2(n_interp+1-I)=ZA(ibeam,K+I)
    ENDDO
    CALL SPLINE_INTERPOLATE(XD2_1D,ZD2,n_interp,XI(1),PhotonFlux_nucleus_ChFF_nob)
    PhotonFlux_nucleus_ChFF_nob=PhotonFlux_nucleus_ChFF_nob/x
    RETURN
  END FUNCTION PhotonFlux_nucleus_ChFF_nob

  FUNCTION PhotonFlux_nucleus_ChFF_nob_eval(ibeam,Ega,gamma,Rval,wVal,aaVal)
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_nucleus_ChFF_nob_eval
    INTEGER,INTENT(IN)::ibeam
    ! Ega in unit of GeV
    ! Rval and aaVal in unit of fm
    ! gamma and wVal are dimensionless
    REAL(KIND(1d0)),INTENT(IN)::Ega,gamma,Rval,wVal,aaVal
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    INTEGER::npoints
    REAL(KIND(1d0))::bupper,blower,Egaoga
    INTEGER::ibeam_common
    COMMON /nucleus_ChFF_nob_fxn_beam/ ibeam_common
    REAL(KIND(1d0))::Ega_common,gamma_common,Rval_common,wVal_common,aaVal_common
    COMMON /nucleus_ChFF_nob_fxn/ Ega_common, gamma_common, &
         Rval_common, wVal_common, aaVal_common
    npoints=50000
    ibeam_common=ibeam
    Ega_common=Ega ! in unit of GeV
    gamma_common=gamma
    ! converting from fm to GeV-1
    Rval_common=Rval/GeVm12fm
    wVal_common=wVal
    ! converting from fm to GeV-1
    aaVal_common=aaVal/GeVm12fm
    ! Exp(bA)=x*mN*b
    bupper=MAX(1.5d0,DLOG(0.8d0*Rval_common*DSQRT(Ega/gamma)))
    blower=DLOG(LOWER_BFactor_Limit*Rval_common*Ega/gamma)
    CALL simpson(PhotonFlux_nucleus_ChFF_nob_fxn,blower,bupper,&
         PhotonFlux_nucleus_ChFF_nob_eval,npoints)
    RETURN
  END FUNCTION PhotonFlux_nucleus_ChFF_nob_eval

  FUNCTION PhotonFlux_nucleus_ChFF_nob_fxn(bA)
    ! b is in unit of GeV-1
    ! Exp(bA)=x*mN*b
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_nucleus_ChFF_nob_fxn
    REAL(KIND(1d0)),INTENT(IN)::bA
    REAL(KIND(1d0))::b
    REAL(KIND(1d0))::btilde
    INTEGER::ibeam
    COMMON /nucleus_ChFF_nob_fxn_beam/ ibeam
    REAL(KIND(1d0))::Ega, gamma, Rval, wVal, aaVal
    COMMON /nucleus_ChFF_nob_fxn/ Ega, gamma, Rval, wVal, aaVal
    REAL(KIND(1d0)),PARAMETER::twopi=6.2831853071795864769252867d0
    b=DEXP(bA)/(Ega/gamma)
    IF(b.LE.0d0)THEN
       PhotonFlux_nucleus_ChFF_nob_fxn=0d0
       RETURN
    ENDIF
    PhotonFlux_nucleus_ChFF_nob_fxn=PhotonNumberDensity_ChargeFormFactor_WS(b,Ega,gamma,&
               Rval,wVal,aaVal,3d0,0.7d0,ibeam)
!    PhotonFlux_nucleus_ChFF_nob_fxn=PhotonNumberDensity(b,Ega,gamma)
    PhotonFlux_nucleus_ChFF_nob_fxn=PhotonFlux_nucleus_ChFF_nob_fxn*twopi*b**2
    RETURN
  END FUNCTION PhotonFlux_nucleus_ChFF_nob_fxn

  SUBROUTINE generate_Q2_epa_proton_iWW(x,Q2max,Q2)
    use interpolation
    implicit none
    real(kind(1d0)),intent(in)::x,Q2max
    real(kind(1d0)),intent(out)::Q2
    real(kind(1d0)),parameter::mp=0.938272081d0 ! proton mass in unit of GeV
    real(kind(1d0)),parameter::mupomuN=2.793d0
    real(kind(1d0)),parameter::Q02=0.71d0 ! in unit of GeV**2
    real(kind(1d0))::mp2,Q2min,xx
    real(kind(1d0))::xmax,xmaxQ2max
    real(kind(1d0))::logQ2oQ02max,logQ2oQ02min
    real(kind(1d0))::log10xmaxQ2maxm1,tlog10xm1
    real(kind(1d0))::log10xm1
    real(kind(1d0))::max_Q2,max_fun,logQ2oQ02now
    real(kind(1d0))::Q2max_save
    data Q2max_save /-1d0/
    save Q2max_save
    real(kind(1d0)),dimension(100)::p_x_array
    real(kind(1d0)),dimension(100)::p_xmax_array
    real(kind(1d0)),dimension(100)::p_fmax_array
    save p_x_array,p_xmax_array,p_fmax_array
    real(kind(1d0))::r1,r2,w
    integer::i,j
    if(x.ge.1d0.or.x.le.0d0)then
       write(*,*)"ERROR: x>=1 or x<=0"
       stop 1
    endif
    mp2=mp**2
    Q2min=mp2*x**2/(1d0-x)

    xmax=xmaxvalue(mp2,Q2max)

    if(Q2max.le.Q2min.or.x.ge.xmax)then
       Q2=Q2max
       return
    endif

    logQ2oQ02max=dlog(Q2max/Q02)
    logQ2oQ02min=dlog(Q2min/Q02)
    
    if(Q2max_save.ne.Q2max)then
       ! we need to generate the grid first
       Q2max_save=Q2max
       xmaxQ2max=xmaxvalue(mp2,Q2max)
       log10xmaxQ2maxm1=dlog10(1d0/xmaxQ2max)
       do i=0,9
          do j=0,9
             tlog10xm1=log10xmaxQ2maxm1+0.1D0*dble(j)+dble(i)
             p_x_array(1+j+i*10)=tlog10xm1
             xx=10d0**(-tlog10xm1)
             if(i.eq.0.and.j.eq.0)then
                max_Q2=logQ2oQ02max
                max_fun=proton_iWW_distfun(xx,max_Q2)
                p_xmax_array(1+j+i*10)=max_Q2
                p_fmax_array(1+j+i*10)=max_fun
             else
                call calc_proton_iWW_distfun_max(xx,logQ2oQ02min,&
                     logQ2oQ02max,max_Q2)
                max_fun=proton_iWW_distfun(xx,max_Q2)
                p_xmax_array(1+j+i*10)=max_Q2
                p_fmax_array(1+j+i*10)=max_fun
             endif
          enddo
       enddo
    endif
    log10xm1=DLOG10(1d0/x)
    call spline_interpolate(p_x_array,p_fmax_array,100,&
         log10xm1,max_fun)
    logQ2oQ02now=logQ2oQ02min
    DO WHILE(.TRUE.)
       r1=rand() ! a random float number between 0 and 1
       logQ2oQ02now=(logQ2oQ02max-logQ2oQ02min)*r1+logQ2oQ02min
       w=proton_iWW_distfun(x,logQ2oQ02now)/max_fun
       r2=rand() ! a random float number between 0 and 1
       if(r2.le.w)exit
    ENDDO
    Q2=dexp(logQ2oQ02now)*Q02
    return
  END SUBROUTINE generate_Q2_epa_proton_iWW

  SUBROUTINE generate_Q2_epa_ion_ChFF(ion_Form,ibeam,x,Q2max,RA,aA,wA,Q2)
    use interpolation
    implicit none
    integer,intent(in)::ion_Form ! Form1: Q**2=kT**2+(mn*x)**2, Qmin**2=(mn*x)**2;
                                 ! Form2: Q**2=kT**2/(1-x)+(mn*x)**2/(1-x), Qmin**2=(mn*x)**2/(1-x) [more exact]
    integer,intent(in)::ibeam ! first beam: ibeam=1; second beam: ibeam=2
    real(kind(1d0)),intent(in)::x,Q2max
    real(kind(1d0)),intent(in)::RA,aA,wA ! RA and aA are in fm
    real(kind(1d0)),intent(out)::Q2
    real(kind(1d0)),parameter::mn=0.9315d0 ! averaged nucleon mass in unit of GeV
    real(kind(1d0)),parameter::Q02=0.71d0 ! in unit of GeV**2
    real(kind(1d0))::mn2,Q2min,xx,RAA,aAA
    real(kind(1d0))::xmax,xmaxQ2max
    real(kind(1d0))::logQ2oQ02max,logQ2oQ02min
    real(kind(1d0))::log10xmaxQ2maxm1,tlog10xm1
    real(kind(1d0))::log10xm1
    real(kind(1d0))::max_Q2,max_fun,logQ2oQ02now
    real(kind(1d0)),dimension(2)::Q2max_save
    data Q2max_save /-1d0,-1d0/
    save Q2max_save
    real(kind(1d0)),dimension(2,100)::A_x_array
    real(kind(1d0)),dimension(2,100)::A_xmax_array
    real(kind(1d0)),dimension(2,100)::A_fmax_array
    save A_x_array,A_xmax_array,A_fmax_array
    real(kind(1d0))::r1,r2,w
    integer::i,j
    if(ion_Form.ne.1.and.ion_Form.ne.2)then
       write(*,*)"ERROR: ion_Form =/= 1 or 2"
       stop 1
    endif
    if(ibeam.ne.1.and.ibeam.ne.2)then
       write(*,*)"ERROR: ibeam =/= 1 or 2"
       stop 1
    endif
    if(x.ge.1d0.or.x.le.0d0)then
       write(*,*)"ERROR: x>=1 or x<=0"
       stop 1
    endif
    mn2=mn**2
    if(ion_Form.eq.2)then
       Q2min=mn2*x**2/(1d0-x)
    else
       Q2min=mn2*x**2
    endif

    RAA=RA/GeVm12fm ! from fm to GeV-1
    aAA=aA/GeVm12fm ! from fm to GeV-1
      

    xmax=xmaxvalue(mn2,Q2max)

    if(Q2max.le.Q2min.or.x.ge.xmax)then
       Q2=Q2max
       return
    endif

    logQ2oQ02max=dlog(Q2max/Q02)
    logQ2oQ02min=dlog(Q2min/Q02)

    if(Q2max_save(ibeam).ne.Q2max)then
       ! we need to generate the grid first
       Q2max_save(ibeam)=Q2max
       xmaxQ2max=xmaxvalue(mn2,Q2max)
       log10xmaxQ2maxm1=dlog10(1d0/xmaxQ2max)
       do i=0,9
          do j=0,9
             tlog10xm1=log10xmaxQ2maxm1+0.1D0*dble(j)+dble(i)
             A_x_array(ibeam,1+j+i*10)=tlog10xm1
             xx=10d0**(-tlog10xm1)
             if(i.eq.0.and.j.eq.0)then
                max_Q2=logQ2oQ02max
                max_fun=ion_ChFF_distfun(ion_Form,aAA,RAA,wA,&
                     Q2min,xx,max_Q2)
                A_xmax_array(ibeam,1+j+i*10)=max_Q2
                A_fmax_array(ibeam,1+j+i*10)=max_fun
             else
                call calc_ion_ChFF_distfun_max(ion_Form,aAA,RAA,wA,&
                     Q2min,xx,logQ2oQ02min,logQ2oQ02max,max_Q2)
                max_fun=ion_ChFF_distfun(ion_Form,aAA,RAA,wA,&
                     Q2min,xx,max_Q2)
                A_xmax_array(ibeam,1+j+i*10)=max_Q2
                A_fmax_array(ibeam,1+j+i*10)=max_fun
             endif
          enddo
       enddo
    endif
    log10xm1=DLOG10(1d0/x)
    call spline_interpolate(A_x_array(ibeam,1:100),&
         A_fmax_array(ibeam,1:100),100,log10xm1,max_fun)
    logQ2oQ02now=logQ2oQ02min
    DO WHILE(.TRUE.)
       r1=rand() ! a random float number between 0 and 1
       logQ2oQ02now=(logQ2oQ02max-logQ2oQ02min)*r1+logQ2oQ02min
       w=ion_ChFF_distfun(ion_Form,aAA,RAA,wA,Q2min,x,logQ2oQ02now)/max_fun
       r2=rand() ! a random float number between 0 and 1
       if(r2.le.w)exit
    ENDDO
    Q2=dexp(logQ2oQ02now)*Q02
    return
  END SUBROUTINE generate_Q2_epa_ion_ChFF

  function xmaxvalue(mp2,q2max)
    implicit none
    real(kind(1d0))::xmaxvalue
    real(kind(1d0)),intent(in)::mp2,q2max
    xmaxvalue=(DSQRT(Q2max*(4d0*mp2+Q2max))-Q2max)/(2d0*mp2)
    return
  end function xmaxvalue

  ! (1 - x)*FE + x^2/2*FM - (1 - x)*Q2min/Q2*FE
  function proton_iWW_distfun(xx,logQ2oQ02)
    implicit none
    real(kind(1d0))::proton_iWW_distfun
    real(kind(1d0)),intent(in)::xx,logQ2oQ02
    real(kind(1d0)),parameter::mp=0.938272081d0 ! proton mass in unit of GeV
    real(kind(1d0)),parameter::mupomuN=2.793d0
    real(kind(1d0)),parameter::Q02=0.71d0 ! in unit of GeV**2
    real(kind(1d0))::expp,mp2
    mp2=mp**2
    expp=DEXP(logQ2oQ02)
    proton_iWW_distfun=(-8D0*mp2**2*xx**2+expp**2*mupomuN**2*Q02**2*&
         (2D0-2D0*xx+xx**2)+2D0*expp*mp2*Q02*(4D0-4D0*xx&
         +mupomuN**2*xx**2))/(2D0*expp*(1D0+expp)**4*Q02&
         *(4D0*mp2+expp*Q02))
    return
  end function proton_iWW_distfun

  function proton_iWW_mdistfun_min(logQ2oQ02)
    implicit none
    real(kind(1d0))::proton_iWW_mdistfun_min
    real(kind(1d0)),intent(in)::logQ2oQ02
    real(kind(1d0))::x
    common/proton_iWW_dist/x
      
    proton_iWW_mdistfun_min=proton_iWW_distfun(x,logQ2oQ02)
    proton_iWW_mdistfun_min=-proton_iWW_mdistfun_min
    return
  end function proton_iWW_mdistfun_min

  subroutine calc_proton_iWW_distfun_max(xx,logQ2oQ02min,&
       logQ2oQ02max,logQ2oQ02)
    use fmin_module
    implicit none
    real(kind(1d0)),intent(in)::xx
    real(kind(1d0)),intent(out)::logQ2oQ02
    real(kind(1d0)),intent(in)::logQ2oQ02min,logQ2oQ02max
    real(kind(1d0)),parameter::tol=1d-8
    real(kind(1d0))::x
    common/proton_iWW_dist/x
    x=xx
    logQ2oQ02=fmin(proton_iWW_mdistfun_min,logQ2oQ02min,&
         logQ2oQ02max,tol)
    return
  end subroutine calc_proton_iWW_distfun_max

  function ion_ChFF_distfun(ion_Form,aAA,RAA,wA,Q2min,xx,logQ2oQ02)
    implicit none
    integer,intent(in)::ion_Form ! Form1: Q**2=kT**2+(mn*x)**2, Qmin**2=(mn*x)**2;
                                 ! Form2: Q**2=kT**2/(1-x)+(mn*x)**2/(1-x), Qmin**2=(mn*x)**2/(1-x) [more exact]
    real(kind(1d0)),intent(in)::Q2min,RAA,aAA,wA
    real(kind(1d0)),parameter::Q02=0.71d0 ! in unit of GeV**2
    real(kind(1d0))::ion_ChFF_distfun
    real(kind(1d0)),intent(in)::xx,logQ2oQ02
    real(kind(1d0))::qq,FchA,expp
    expp=DEXP(logQ2oQ02)*Q02
    if(ion_Form.eq.2)then
       qq=dsqrt((1d0-xx)*expp)
    else
       qq=dsqrt(expp)
    endif
    ! keep the 10 terms
    FchA=ChargeFormFactor_WoodsSaxon(qq,RAA,wA,aAA,10)
    ion_ChFF_distfun=(1d0-Q2min/expp)*FchA**2
    return
  end function ion_ChFF_distfun

  function ion_ChFF_mdistfun_min(logQ2oQ02)
    implicit none
    real(kind(1d0)),intent(in)::logQ2oQ02
    real(kind(1d0))::ion_ChFF_mdistfun_min
    integer::ion_Form
    common/ion_form_par/ion_Form
    real(kind(1d0))::aAA,RAA,wA,x,Q2min
    common/ion_ChFF_dist/aAA,RAA,wA,x,Q2min

    ion_ChFF_mdistfun_min=ion_ChFF_distfun(ion_Form,aAA,RAA,wA,&
         Q2min,x,logQ2oQ02)
    ion_ChFF_mdistfun_min=-ion_ChFF_mdistfun_min
    return
  end function ion_ChFF_mdistfun_min

  subroutine calc_ion_ChFF_distfun_max(ion_Form,aAA,RAA,wA,&
       Q2min,xx,logQ2oQ02min,logQ2oQ02max,logQ2oQ02)
    use fmin_module
    implicit none
    integer,intent(in)::ion_Form
    real(kind(1d0)),intent(in)::aAA,RAA,wA,Q2min
    real(kind(1d0)),intent(in)::xx
    real(kind(1d0)),intent(out)::logQ2oQ02
    real(kind(1d0)),intent(in)::logQ2oQ02min,logQ2oQ02max
    real(kind(1d0)),parameter::tol=1d-8
    integer::ion_Form_com
    common/ion_form_par/ion_Form_com
    real(kind(1d0))::aAA_com,RAA_com,wA_com,x,Q2min_com
    common/ion_ChFF_dist/aAA_com,RAA_com,wA_com,x,Q2min_com
    ion_Form_com=ion_Form
    Q2min_com=Q2min
    x=xx
    aAA_com=aAA
    RAA_com=RAA
    wA_com=wA
    logQ2oQ02=fmin(ion_ChFF_mdistfun_min,logQ2oQ02min,&
         logQ2oQ02max,tol)
    return
  end subroutine calc_ion_ChFF_distfun_max

  subroutine progress(j,nmax,forceinit)
    implicit none
    integer,intent(in)::j,nmax
    logical,intent(in),optional::forceinit
    integer::k
    character(:), allocatable :: bar, bar0
    character(5)::nmax_str
    !character(len=)::bar="???% |                                     |"
    integer::init=0
    save init,bar,bar0,nmax_str
    IF(present(forceinit))THEN
       IF(forceinit)THEN
          init=0
          IF(ALLOCATED(bar))THEN
             DEALLOCATE(bar)
          ENDIF
          IF(ALLOCATED(bar0))THEN
             DEALLOCATE(bar0)
          ENDIF
       ENDIF
    ENDIF
    IF(init.EQ.0)THEN
       allocate(character(nmax+7) :: bar)
       allocate(character(nmax+7) :: bar0)
       bar(1:6)="???% |"
       do k=1,nmax
          bar(6+k:6+k)=" "
       enddo
       bar(nmax+7:nmax+7)="|"
       bar0=bar
       !bar="???% |"//repeat(' ',nmax)//"|"
       write(unit=nmax_str,fmt="(i5)") nmax+7
       nmax_str=adjustl(nmax_str)
       init=1
    ENDIF
    bar=bar0
    write(unit=bar(1:3),fmt="(i3)") INT(100*DBLE(j)/DBLE(nmax))
    do k=1, j
       bar(6+k:6+k)="*"
    enddo
    ! print the progress bar.
    write(unit=6,fmt="(a1,a"//trim(nmax_str)//")",advance="no") char(13), bar
    if (j.NE.nmax) then
       flush(unit=6)
    else
       write(unit=6,fmt=*)
    endif
    return
  end subroutine progress
END MODULE ElasticFlux
