MODULE ElasticFlux
  ! photon flux of a single elastic (coherent) source: proton or nucleus
  USE OpticalGlauber_Geometry
  USE NINTLIB ! for multiple dimensional integrations
  USE interpolation
  IMPLICIT NONE
  PRIVATE
  REAL(KIND(1d0)),PUBLIC::alphaem_elasticphoton=0.0072992700729927005d0
  REAL(KIND(1d0)),PUBLIC::aqedup=0.0072992700729927005d0
  REAL(KIND(1d0)),PARAMETER,PUBLIC::LOWER_BFactor_Limit=1D-1
  REAL(KIND(1d0)),PARAMETER,PUBLIC::GeVm12fm=0.1973d0 ! GeV-1 to fm
  PUBLIC::PhotonNumberDensity ! EDFF
  ! ChFF
  PUBLIC::PhotonNumberDensity_ChargeFormFactor_WS,PhotonNumberDensity_ChargeFormFactor_proton
  PUBLIC::PhotonFlux_proton_EDFF_nob,PhotonFlux_proton_ChFF_nob
  PUBLIC::PhotonFlux_nucleus_EDFF_nob,PhotonFlux_nucleus_ChFF_nob
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

  FUNCTION PhotonNumberDensity_ChargeFormFactor_proton(b,Ega,gamma,CORRECTIONQ,MINUSDIPOLESQUAREQ)
    ! It gives us the photon number density with Z=1 and alpha=1
    ! b should be written in unit of GeV-1
    ! 1 GeV^-1 = 0.1973e-15 m = 0.1973 fm
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonNumberDensity_ChargeFormFactor_proton
    REAL(KIND(1d0)),INTENT(IN)::b,Ega,gamma
    LOGICAL,INTENT(IN),OPTIONAL::CORRECTIONQ ! whether we need to include the correction beyond the dipole
                                             ! the default is no (.FALSE.)
    LOGICAL,INTENT(IN),OPTIONAL::MINUSDIPOLESQUAREQ ! whether we want to minus the dipole square
                                                    ! the default is no (.FALSE.)
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
    dipole=0d0
    IF(PRESENT(MINUSDIPOLESQUAREQ))THEN
       IF(MINUSDIPOLESQUAREQ)THEN
          dipole=integral
       ENDIF
    ENDIF
    corr=0d0
    IF(PRESENT(CORRECTIONQ))THEN
       IF(CORRECTIONQ)THEN
          ! include the correction beyond the dipole
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
    PhotonNumberDensity_ChargeFormFactor_proton=one/pi2*integral
    RETURN
  END FUNCTION PhotonNumberDensity_ChargeFormFactor_proton

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
    REAL(KIND(1d0))::log10x1,kTmax
    REAL(KIND(1d0))::rescaling_bmax_save
    SAVE rescaling_bmax_save
    INTEGER,PARAMETER::n_interp=6
    REAL(KIND(1d0)),DIMENSION(n_interp)::XA2,YA2
    REAL(KIND(1d0)),PARAMETER::PIPI=3.14159265358979323846264338328d0
    INTEGER::iter,npoints
    REAL(KIND(1d0))::integ1,integ2,integ3,integ4
    INTEGER,PARAMETER::itermax=8
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
                  DSQRT(3.5d0),integral)
             PhotonNumberDensity_ChargeFormFactor_proton_correction=integral
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
          STOP
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
          STOP
       ENDIF
       kTmax=DSQRT(3.5d0) ! match with Q2max in PND_ChargeFormFactor_proton_corr_fxn
       YA(NYA)=DBLE(log10xmin)
       DO I=1,NXA
          IF(NYA_save.EQ.0)THEN
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
          ENDIF
          IF(XA(I).EQ.0d0)THEN
             ZA(I,J)=0d0
             CYCLE
          ENDIF
          DO J=1,NYA
             ! use the ones that have been calculated before
             IF(NYA_save.GT.0.AND.J.LE.NYA_save)THEN
                ZA(I,J)=ZA_save(I,J)
                CYCLE
             ENDIF
             xga2=10d0**(YA(J)) ! x_gamma
             Egaoga2=xga2*mproton
             Ega2=Egaoga2*gamma
             b_common=XA(I) ! = b
             Egaoga_common=Egaoga2
             npoints=10000
             CALL trapezoid_integration(npoints,PND_ChargeFormFactor_proton_corr_fxn,&
                  kTmax,integral)
             ! we try to do some numerical improvement
             integ4=0d0
             iter=1
             DO WHILE((iter.EQ.1.OR.&
                  DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)&
                  .AND.iter.LE.itermax)
                integ4=integral
                ! increase the points by a factor of 2
                npoints=npoints*2
                CALL trapezoid_integration(npoints,PND_ChargeFormFactor_proton_corr_fxn,&
                     kTmax,integral)
                iter=iter+1
             END DO
             IF(DABS(integ4/integral).GT.1.5d0.OR.DABS(integ4/integral).LT.0.67d0)THEN
                WRITE(*,*)"WARNING: the integral is not stable (b,Ega,gamma,npoints)=",&
                     b_common,Ega2,gamma,npoints
                WRITE(*,*)"WARNING: integral in two iterations #1:",integ4,integral
             ENDIF
             IF(ISNAN(integral))THEN
                WRITE(*,*)"ERROR: the integral is not stable (b,Ega,gamma,RA,wA,aA)=",&
                     b_common,Ega2,gamma
                STOP
             ENDIF
             ZA(I,J)=integral
          ENDDO
       ENDDO
       init=1
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
             STOP
          ENDIF
          IF(XA(k).LT.b)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
    REAL(KIND(1d0)),PARAMETER::Q2max=3.5d0 ! in GeV**2: the correction vanishes above (see ratio_dde1)
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
    IF(Q2.GE.Q2max)THEN
       PND_ChargeFormFactor_proton_corr_fxn=0d0
       RETURN
    ENDIF
    dipole=one/(one+Q2*aa**2)**2
    CFF=dipole*(ratio_dde1(Q2)-one)
    IF(ISNAN(CFF))THEN
       PRINT *, "ChargeFormFactor is NaN with ",DSQRT(Q2)
       STOP
       CFF=0d0
    ENDIF
    PND_ChargeFormFactor_proton_corr_fxn=pref*BESSJ1(bkT)*CFF
    RETURN
  END FUNCTION PND_ChargeFormFactor_proton_corr_fxn

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
       STOP
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
             STOP
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
             STOP
          ENDIF
          YA(ibeam,NYA)=DBLE(log10xmin)
          DO I=1,NXA
             IF(I.EQ.1)THEN
                CALL progress(INT(I*50d0/NXA),50,.TRUE.)
             ELSE
                CALL progress(INT(I*50d0/NXA),50)
             ENDIF
             IF(XA(ibeam,I).EQ.0d0)THEN
                ZA(ibeam,I,J)=0d0
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
                      STOP
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
       STOP
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
             STOP
          ENDIF
          IF(XA(ibeam,k).LT.b)THEN
             WRITE(*,*)"Error: k is not proper #2"
             STOP
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
             STOP
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
       STOP
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
       STOP
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

  FUNCTION PhotonFlux_proton_ChFF_nob(x,gamma,CORRECTIONQ)
    ! set PNOHARD=1
    ! for proton with ChFF
    IMPLICIT NONE
    REAL(KIND(1d0))::PhotonFlux_proton_ChFF_nob
    REAL(KIND(1d0)),INTENT(IN)::x,gamma
    LOGICAL,INTENT(IN),OPTIONAL::CORRECTIONQ ! whether we need to include the correction beyond the dipole
                                             ! the default is no (.FALSE.)
    LOGICAL::corrQ
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0)),PARAMETER::mproton=0.938272081d0 ! the mass of proton (GeV)
    REAL(KIND(1d0))::alpha
    SAVE alpha
    REAL(KIND(1d0))::aptilde,Egaoga
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
       corrQ=CORRECTIONQ
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
          STOP
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
    PhotonFlux_proton_ChFF_corr_nob_fxn=PhotonNumberDensity_ChargeFormFactor_proton(b,Ega,gamma,.TRUE.,.TRUE.)
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
       STOP
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
             STOP
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
       STOP
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
          STOP
       ENDIF
       IF(XA(ibeam,k).LT.XI(1))THEN
          WRITE(*,*)"Error: k is not proper #2"
          STOP
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
          STOP
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
       stop
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
       stop
    endif
    if(ibeam.ne.1.and.ibeam.ne.2)then
       write(*,*)"ERROR: ibeam =/= 1 or 2"
       stop
    endif
    if(x.ge.1d0.or.x.le.0d0)then
       write(*,*)"ERROR: x>=1 or x<=0"
       stop
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
