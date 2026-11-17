MODULE PhotonPhotonFlux
  USE ElasticFlux
  USE OpticalGlauber_Geometry
  USE NINTLIB ! for multiple dimensional integrations
  USE interpolation
  IMPLICIT NONE
  PRIVATE
  PUBLIC::alphaem_elasticphoton
  ! For UPCs photon-photon collisions
  LOGICAL,PUBLIC::USE_CHARGEFORMFACTOR4PHOTON=.FALSE.
  PUBLIC::PhotonPhotonFlux_pp,PhotonPhotonFlux_pp_eval
  PUBLIC::PhotonFlux_proton_EDFF_nob,PhotonFlux_proton_ChFF_nob
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
  LOGICAL,PUBLIC::proton_FF_correction=.FALSE. ! when USE_CHARGEFORMFACTOR4PHOTON=.TRUE., 
                                               ! it will include the correction beyond proton dipole form factor
                                               ! when proton_FF_correction=.TRUE.
                                               ! it will impact both pp and pA UPCs
  INTEGER,PUBLIC::IEPSILON_EDFF=0  ! 0: bmin = RA; 1: bmin = RA+a
  LOGICAL,PRIVATE::CALC_DELTABQ=.FALSE. ! whether we want to calcuate DeltaB (.TRUE.) or not (.FALSE.)
  LOGICAL,PUBLIC::GENERATE_DELTAB_GRID=.FALSE. ! if true, it will always first generate the grid
  LOGICAL,PUBLIC::GENERATE_PhotonPhotonFlux_GRID=.TRUE. ! if true, it will alwyas first generate the grid for PhotonPhotonFlux
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
    IF(.NOT.print_banner)THEN
       include "banner.inc"
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_pp=0d0
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
          DO J=1,nseg
             log10x1=-1d0/DBLE(nseg)*DBLE(J-1)+DBLE(I)
             K=K+1
             XD_1D(K)=log10x1
             YD_1D(K)=log10x1
          ENDDO
       ENDDO
       IF(K.NE.MX)THEN
          WRITE(*,*)"ERROR: K != MX"
          STOP
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             ! use the ones that have been calculated before
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1.AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_save(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
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
             STOP
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
             STOP
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
          ZI(1)=PhotonPhotonFlux_pp_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    pnohadval=PhotonPhotonFlux_pp_eval(x1,x2,.TRUE.)
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pp=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.DABS(ZI(1)/pnohadval).GT.1D2)THEN
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
    REAL(KIND(1d0))::cmenergy
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_pp_eval=0d0
       RETURN
    ENDIF
    CALC_DELTABQ=.FALSE.
    IF(init.EQ.0)THEN
       gamma1_common=ebeam(1)/mproton
       gamma2_common=ebeam(2)/mproton
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
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2)) ! in unit of GeV
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
    E1_common=ebeam(1)*x1
    E2_common=ebeam(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          PhotonPhotonFlux_pp_eval=PhotonFlux_proton_ChFF_nob(x1_common,gamma1_common,proton_FF_correction)
          PhotonPhotonFlux_pp_eval=PhotonPhotonFlux_pp_eval*&
               PhotonFlux_proton_ChFF_nob(x2_common,gamma2_common,proton_FF_correction)
       ELSE
          PhotonPhotonFlux_pp_eval=PhotonFlux_proton_EDFF_nob(x1_common,gamma1_common)
          PhotonPhotonFlux_pp_eval=PhotonPhotonFlux_pp_eval*&
               PhotonFlux_proton_EDFF_nob(x2_common,gamma2_common)
       ENDIF
       PhotonPhotonFlux_pp_eval=MAX(PhotonPhotonFlux_pp_eval,0d0)
       RETURN
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact*x1)
    ! b2*x2*mproton = Exp(bA(2))
    aax(2)=DLOG(bfact*x2)
    CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax,bbx,3,sub_num,1,1d-5,&
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
          CALL ROMBERG_ND(PhotonPhotonFlux_pp_fxn,aax,bbx,3,sub_num,1,1d-5,&
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
       PhotonPhotonFlux_pp_eval=MAX(PhotonPhotonFlux_pp_eval,0d0)
    ELSE
       PhotonPhotonFlux_pp_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
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
       include "banner.inc"
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
                MX_save(1)=MX
                MY_save(1)=MY
                IF(ALLOCATED(ZD1_save))THEN
                   DEALLOCATE(ZD1_save)
                ENDIF
                ALLOCATE(ZD1_save(MX+1,MY+1))
                DO I=1,MX+1
                   DO J=1,MY+1
                      ZD1_save(I,J)=ZD(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(gridready(2))THEN
                MX_save(2)=MX
                MY_save(2)=MY
                IF(ALLOCATED(ZD2_save))THEN
                   DEALLOCATE(ZD2_save)
                ENDIF
                ALLOCATE(ZD2_save(MX+1,MY+1))
                DO I=1,MX+1
                   DO J=1,MY+1
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
             STOP
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
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save(igrid)+1)THEN
          K=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save(igrid)+1)THEN
          K=MX_save(igrid)+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #1 !"
          STOP
       ENDIF
       L=nseg*(-ilog10x2-1)
       DO I=L+1,L+nseg+1
          if(YD_1D(I).LE.YI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save(igrid)+1)THEN
          L=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          L=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save(igrid)+1)THEN
          L=MY_save(igrid)+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #2 !"
          STOP
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
    CALC_DELTABQ=.TRUE.
    IF(init.EQ.0)THEN
       gamma1_common=ebeam(1)/mproton
       gamma2_common=ebeam(2)/mproton
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
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2)) ! in unit of GeV
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
    E1_common=ebeam(1)*x1
    E2_common=ebeam(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
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
    REAL(KIND(1d0))::Ngamma1,Ngamma2
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
       STOP
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
    IF(.NOT.USE_CHARGEFORMFACTOR4PHOTON)THEN
       ! EDFF
       Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
       Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    ELSE
       ! ChFF
       Ngamma1=PhotonNumberDensity_ChargeFormFactor_proton(b1,E1_common,gamma1_common,&
            proton_FF_correction)
       Ngamma2=PhotonNumberDensity_ChargeFormFactor_proton(b2,E2_common,gamma2_common,&
            proton_FF_correction)
    ENDIF
    PhotonPhotonFlux_pp_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
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
       include "banner.inc"
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
          STOP
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
             STOP
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
             STOP
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
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_common=ebeam(2)/mproton
          gamma2_common=ebeam(1)/mN
       ELSE
          gamma1_common=ebeam(1)/mproton
          gamma2_common=ebeam(2)/mN
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
          STOP
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in p+"//TRIM(Aname)//" collisions"
       CALL GetNuclearInfo(Aname,A_common,Z2,RA_common,aaval,wval)
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2))
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
       E1_common=ebeam(2)*x2
       E2_common=ebeam(1)*x1
    ELSE
       x1_common=x1
       x2_common=x2
       E1_common=ebeam(1)*x1
       E2_common=ebeam(2)*x2
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
       STOP
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
    IF(.NOT.print_banner)THEN
       include "banner.inc"
       print_banner=.TRUE.
    ENDIF
    IF(x1.LE.1d-99.OR.x2.LE.1d-99.OR.x1.GT.1d0.OR.x2.GT.1d0)THEN
       PhotonPhotonFlux_pA_WoodsSaxon=0d0
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
          STOP
       ENDIF
       XD_1D(MX+1)=DBLE(log10xmin)
       YD_1D(MY+1)=DBLE(log10xmin)
       DO I=1,MX+1
          xx1=10d0**(XD_1D(I))
          DO J=1,MY+1
             IF(MX_save.GT.0.AND.MY_save.GT.0.AND.I.LE.MX_save+1.AND.J.LE.MY_save+1)THEN
                ZD(I,J)=ZD_save(I,J)
                CYCLE
             ENDIF
             xx2=10d0**(YD_1D(J))
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
             STOP
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
             STOP
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
          ZI(1)=PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2)
       ENDIF
    ENDIF
    ! Let us always evaluate PNOHAD=1 as a reference to compare
    pnohadval=PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2,.TRUE.)
    IF(.NOT.force_pnohad1)THEN
       IF(ISNAN(pnohadval).OR.pnohadval.EQ.0d0)THEN
          PhotonPhotonFlux_pA_WoodsSaxon=0d0
       ELSEIF(ISNAN(ZI(1)).OR.ZI(1).LT.0d0.OR.DABS(ZI(1)/pnohadval).GT.1D2)THEN
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
    IF(x1.LE.0d0.OR.x1.GE.1d0.OR.x2.LE.0d0.OR.x2.GE.1d0)THEN
       PhotonPhotonFlux_pA_WoodsSaxon_eval=0d0
       RETURN
    ENDIF
    CALC_DELTABQ=.FALSE.
    IF(init.EQ.0)THEN
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_common=ebeam(2)/mproton
          gamma2_common=ebeam(1)/mN
       ELSE
          gamma1_common=ebeam(1)/mproton
          gamma2_common=ebeam(2)/mN
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
          STOP
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in p+"//TRIM(Aname)//" collisions"
       CALL GetNuclearInfo(Aname,A_common,Z2,RA_common,aaVal_common,wVal_common)
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2))
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
       E1_common=ebeam(2)*x2
       E2_common=ebeam(1)*x1
    ELSE
       x1_common=x1
       x2_common=x2
       E1_common=ebeam(1)*x1
       E2_common=ebeam(2)*x2
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
    IF(force_pnohad1)THEN
       ! we only use PNOHAD=1
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonFlux_proton_ChFF_nob(x1_common,&
               gamma1_common,proton_FF_correction)
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval*&
            PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common,wVal_common,aaVal_common)
       ELSE
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonFlux_proton_EDFF_nob(x1_common,&
               gamma1_common)
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval*&
               PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
       ENDIF
       PhotonPhotonFlux_pA_WoodsSaxon_eval=MAX(PhotonPhotonFlux_pA_WoodsSaxon_eval,0d0)
       RETURN
    ENDIF
    ! we should choose the lower limit dynamically
    ! b1*x1*mproton = Exp(bA(1))
    aax(1)=DLOG(bfact1*x1_common)
    ! b2*x2*mN = Exp(bA(2))
    aax(2)=DLOG(bfact2*x2_common)
    ! we should also choose the upper limit dynamically for the ion beam
    bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common/GeVm12fm*DSQRT(x2_common*mN)))
    CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
         integral,ind,eval_num)
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
          CALL ROMBERG_ND(PhotonPhotonFlux_pA_WoodsSaxon_fxn,aax,bbx,3,sub_num,1,1d-5,&
               integral,ind,eval_num)
          iter=iter+1
       ENDDO
       bbx(1)=bupper
       !bbx(2)=bupper
       bbx(2)=MAX(1.5d0,DLOG(0.8d0*RA_common/GeVm12fm*DSQRT(x2_common*mN)))
    ENDIF
    IF(integral.LT.0d0)THEN
       printnum=printnum+1
       IF(printnum.LE.5)THEN
          WRITE(*,*)"WARNING: negative photon flux at (x1,x2)=",x1_common,x2_common
          WRITE(*,*)"WARNING: use PNOHAD=1 approx. instead (most probably need to increase bupper)"
          IF(printnum.EQ.5)WRITE(*,*)"WARNING: Further warning will be suppressed"
       ENDIF
       IF(USE_CHARGEFORMFACTOR4PHOTON)THEN
          ! for approximation
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonFlux_proton_ChFF_nob(x1_common,&
               gamma1_common,proton_FF_correction)
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval*&
            PhotonFlux_nucleus_ChFF_nob(2,x2_common,gamma2_common,Z2,RA_common,wVal_common,aaVal_common)
       ELSE
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonFlux_proton_EDFF_nob(x1_common,&
               gamma1_common)
          PhotonPhotonFlux_pA_WoodsSaxon_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval*&
               PhotonFlux_nucleus_EDFF_nob(x2_common,gamma2_common,Z2,RA_common)
       ENDIF
       PhotonPhotonFlux_pA_WoodsSaxon_eval=MAX(PhotonPhotonFlux_pA_WoodsSaxon_eval,0d0)
    ELSE
       PhotonPhotonFlux_pA_WoodsSaxon_eval=TWOPI/(x1*x2)*alpha**2*Z1**2*Z2**2*integral
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
       include "banner.inc"
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
                MX_save(1)=MX
                MY_save(1)=MY
                IF(ALLOCATED(ZD1_save))THEN
                   DEALLOCATE(ZD1_save)
                ENDIF
                ALLOCATE(ZD1_save(MX+1,MY+1))
                DO I=1,MX+1
                   DO J=1,MY+1
                      ZD1_save(I,J)=ZD(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(gridready(2))THEN
                MX_save(2)=MX
                MY_save(2)=MY
                IF(ALLOCATED(ZD2_save))THEN
                   DEALLOCATE(ZD2_save)
                ENDIF
                ALLOCATE(ZD2_save(MX+1,MY+1))
                DO I=1,MX+1
                   DO J=1,MY+1
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
             STOP
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
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save(igrid)+1)THEN
          K=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save(igrid)+1)THEN
          K=MX_save(igrid)+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #1 !"
          STOP
       ENDIF
       L=nseg*(-ilog10x2-1)
       DO I=L+1,L+nseg+1
          if(YD_1D(I).LE.YI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save(igrid)+1)THEN
          L=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          L=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save(igrid)+1)THEN
          L=MY_save(igrid)+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #2 !"
          STOP
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
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_common=ebeam(2)/mproton
          gamma2_common=ebeam(1)/mN
       ELSE
          gamma1_common=ebeam(1)/mproton
          gamma2_common=ebeam(2)/mN
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
          STOP
       ENDIF
       WRITE(*,*)"INFO: Two photon UPCs in p+"//TRIM(Aname)//" collisions"
       CALL GetNuclearInfo(Aname,A_common,Z2,RA_common,aaVal_common,wVal_common)
       ! read the inelastic NN cross section
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2))
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
       E1_common=ebeam(2)*x2
       E2_common=ebeam(1)*x1
    ELSE
       x1_common=x1
       x2_common=x2
       E1_common=ebeam(1)*x1
       E2_common=ebeam(2)*x2
    ENDIF
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
    ENDIF
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
    REAL(KIND(1d0))::Ngamma1,Ngamma2
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
       STOP
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
    IF(.NOT.USE_CHARGEFORMFACTOR4PHOTON)THEN
       Ngamma1=PhotonNumberDensity(b1,E1_common,gamma1_common)
       Ngamma2=PhotonNumberDensity(b2,E2_common,gamma2_common)
    ELSE
       ! converting from fm to GeV-1
       RR=RA_common/GeVm12fm
       aaa=aaVal_common/GeVm12fm
       IF(nuclearA_beam1.NE.0)THEN
          Ngamma1=PhotonNumberDensity_ChargeFormFactor_WS(b1,E1_common,gamma1_common,&
               RR,wVal_common,aaa,3d0,0.7d0,1)
          Ngamma2=PhotonNumberDensity_ChargeFormFactor_proton(b2,E2_common,gamma2_common,&
               proton_FF_correction)
       ELSE
          Ngamma1=PhotonNumberDensity_ChargeFormFactor_proton(b1,E1_common,gamma1_common,&
               proton_FF_correction)
          Ngamma2=PhotonNumberDensity_ChargeFormFactor_WS(b2,E2_common,gamma2_common,&
               RR,wVal_common,aaa,3d0,0.7d0,2)
       ENDIF
    ENDIF
    PhotonPhotonFlux_pA_WoodsSaxon_fxn=b1**2*b2**2*pnohad*Ngamma1*Ngamma2
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
       include "banner.inc"
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
          STOP
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
             STOP
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
             STOP
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
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP
       ENDIf
       gamma1_common=ebeam(1)/mN
       gamma2_common=ebeam(2)/mN
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
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2))
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
    E1_common=ebeam(1)*x1
    E2_common=ebeam(2)*x2
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
       STOP
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
       include "banner.inc"
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
          STOP
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
             STOP
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
             STOP
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
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP
       ENDIF
       gamma1_common=ebeam(1)/mN
       gamma2_common=ebeam(2)/mN
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
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2))
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
    E1_common=ebeam(1)*x1
    E2_common=ebeam(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
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
       include "banner.inc"
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
                MX_save(1)=MX
                MY_save(1)=MY
                IF(ALLOCATED(ZD1_save))THEN
                   DEALLOCATE(ZD1_save)
                ENDIF
                ALLOCATE(ZD1_save(MX+1,MY+1))
                DO I=1,MX+1
                   DO J=1,MY+1
                      ZD1_save(I,J)=ZD(1,I,J)
                   ENDDO
                ENDDO
             ENDIF
             IF(gridready(2))THEN
                MX_save(2)=MX
                MY_save(2)=MY
                IF(ALLOCATED(ZD2_save))THEN
                   DEALLOCATE(ZD2_save)
                ENDIF
                ALLOCATE(ZD2_save(MX+1,MY+1))
                DO I=1,MX+1
                   DO J=1,MY+1
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
             STOP
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
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MX_save(igrid)+1)THEN
          K=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          K=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MX_save(igrid)+1)THEN
          K=MX_save(igrid)+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #1 !"
          STOP
       ENDIF
       L=nseg*(-ilog10x2-1)
       DO I=L+1,L+nseg+1
          if(YD_1D(I).LE.YI(1))EXIT
       ENDDO
       IF(I-n_interp/2.GE.1.AND.I-n_interp/2-1+n_interp.LE.MY_save(igrid)+1)THEN
          L=I-n_interp/2-1
       ELSEIF(I-n_interp/2.LT.1)THEN
          L=0
       ELSEIF(I-n_interp/2-1+n_interp.GT.MY_save(igrid)+1)THEN
          L=MY_save(igrid)+1-n_interp
       ELSE
          WRITE(*,*)"Error: you cannot reach here #2 !"
          STOP
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
          PhotonPhotonDeltaB_AB_WoodsSaxon=pnohadval
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
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP
       ENDIF
       gamma1_common=ebeam(1)/mN
       gamma2_common=ebeam(2)/mN
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
       cmenergy=2d0*DSQRT(ebeam(1)*ebeam(2))
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
    E1_common=ebeam(1)*x1
    E2_common=ebeam(2)*x2
    IF(.NOT.PRESENT(FORCEPNOHAD1))THEN
       force_pnohad1=.FALSE.
    ELSE
       force_pnohad1=FORCEPNOHAD1
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
       STOP
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
    s=4d0*ebeam(1)*ebeam(2)
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
          STOP
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
          STOP
       ENDIF
    ELSE
       WRITE(*,*)"ERROR: do not know the collision type = ",collision_type_common
       STOP
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
       include "banner.inc"
       print_banner=.TRUE.
    ENDIF
    s=4d0*ebeam(1)*ebeam(2)
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
       include "banner.inc"
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
    s=4d0*ebeam(1)*ebeam(2)
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
          STOP
       ELSEIF(profile_type_common.EQ.0)THEN
          ! P_{NOHAD}=1
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_pA_WoodsSaxon(x1,x2,.TRUE.)
       ELSE
          WRITE(*,*)"ERROR: do not know the profile type = ",profile_type_common
          STOP
       ENDIF
    ELSEIF(collision_type_common.EQ.3)THEN
       ! AB
       IF(profile_type_common.EQ.1)THEN
          ! Woods-Saxon
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_AB_WoodsSaxon(x1,x2)
       ELSEIF(profile_type_common.EQ.2)THEN
          ! Hard-Sphere
          WRITE(*,*)"ERROR: LgammagammaDeltaB_UPC does not support hard-sphere"
          STOP
       ELSEIF(profile_type_common.EQ.0)THEN
          ! P_{NOHAD}=1
          LgammagammaDeltaB_UPC_fxn=PhotonPhotonDeltaB_AB_WoodsSaxon(x1,x2,.TRUE.)
       ELSE
          WRITE(*,*)"ERROR: do not know the profile type = ",profile_type_common
          STOP
       ENDIF
    ELSE
       WRITE(*,*)"ERROR: do not know the collision type = ",collision_type_common
       STOP
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
       STOP
    ENDIF
    IF(in2.LT.-2.OR.in2.GT.4)THEN
       WRITE(*,*)"ERROR: only the following options for in2: -2 (0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)"
       STOP
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
       IF(nuclearA_beam1.EQ.0.OR.nuclearA_beam2.EQ.0)THEN
          WRITE(*,*)"ERROR: Please set two beams as heavy ions first"
          STOP
       ENDIF
       gamma1_save=ebeam(1)/mN
       gamma2_save=ebeam(2)/mN
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

    IF(in2.LT.-1.OR.in2.GT.4)THEN
       WRITE(*,*)"ERROR: only the following options for in2: -2 (0n+Xn), -1 (Xn), 0 (0n), 1 (1n), 2 (2n), 3 (3n), 4 (4n)"
       STOP
    ENDIF

    IF(in2.EQ.-2)THEN
       Pbreak_pA_WoodsSaxon=1d0
       RETURN
    ENDIF
    
    IF(init.EQ.0)THEN
       IF(nuclearA_beam1.NE.0)THEN
          gamma1_save=ebeam(2)/mproton
          gamma2_save=ebeam(1)/mN
       ELSE
          gamma1_save=ebeam(1)/mproton
          gamma2_save=ebeam(2)/mN
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
          STOP
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
       STOP
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in PbreakXn_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in PbreakXn_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
          STOP
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
       STOP
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in PbreakXn_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak1n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak1n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
          STOP
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
       STOP
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak1n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak2n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak2n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
          STOP
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
       STOP
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak2n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak3n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak3n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
          STOP
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
       STOP
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak3n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_target.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_target < 1"
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj_save(ibeam_proj).NE.gamma_proj&
         .OR.gamma_target_save(ibeam_proj).NE.gamma_target)THEN
       WRITE(*,*)"ERROR: gamma_proj,gamma_target are not consistent in Pbreak4n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (gamma_proj,gamma_target)=",gamma_proj_save(ibeam_proj),&
            gamma_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (gamma_proj,gamma_target)=",gamma_proj,gamma_target
       STOP
    ENDIF
    IF(A_target_save(ibeam_proj).NE.A_target&
         .OR.Z_target_save(ibeam_proj).NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak4n_LO_AB_WS"
       WRITE(*,*)"INFO: ibeam_proj=",ibeam_proj
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save(ibeam_proj),&
            Z_target_save(ibeam_proj)
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(ibeam_proj,k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
    ENDIF
    IF(gamma_proj.LT.1d0)THEN
       WRITE(*,*)"Error: gamma_proj < 1"
       STOP
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
          STOP
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
          STOP
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
       STOP
    ENDIF
    IF(A_target_save.NE.A_target&
         .OR.Z_target_save.NE.Z_target)THEN
       WRITE(*,*)"ERROR: A_target,Z_target are not consistent in Pbreak4n_LO_pA"
       WRITE(*,*)"INFO: Saved ones (A_target,Z_target)=",A_target_save,&
            Z_target_save
       WRITE(*,*)"INFO: New ones (A_target,Z_target)=",A_target,Z_target
       STOP
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
             STOP
          ENDIF
          IF(XA(k).LT.db)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
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

END MODULE PhotonPhotonFlux
