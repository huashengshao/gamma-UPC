! Dump of the new engine's building blocks for the independent checks of the benchmark report
! (family A): That_AB(b) and P_had(b), the photon number density N(x,b), and the neutron-tagging
! probability with the legacy value and, for EDFF, an independent adaptive integration of the
! lowest-order Xn probability.   bench_blocks <namelist file>   (same /config/ as bench_ab)
MODULE blocks_ref
  USE NINTLIB
  IMPLICIT NONE
  REAL(KIND(1d0))::bb_c,pref_c,gp_c,mp_c
  INTEGER::at_c,zt_c
CONTAINS
  FUNCTION f_mxn(t)
    ! sigma_abs(E) N_EDFF(b, x(E)) in t = ln E (E in GeV, target frame)
    USE photoabsorption
    REAL(KIND(1d0))::f_mxn
    REAL(KIND(1d0)),INTENT(IN)::t
    REAL(KIND(1d0)),EXTERNAL::BESSK0,BESSK1
    REAL(KIND(1d0))::E,x,q,z
    E=DEXP(t)
    x=E/pref_c
    q=x*mp_c
    z=q*bb_c
    IF(z.GT.600d0)THEN
       f_mxn=0d0
       RETURN
    ENDIF
    f_mxn=GAMMAABS(E*1d3,at_c,zt_c)*0.1d0/0.1973d0**2*q*q*(BESSK1(z)**2+BESSK0(z)**2/gp_c**2)&
         /9.86960440108935861883449099988d0
  END FUNCTION f_mxn
END MODULE blocks_ref

PROGRAM bench_blocks
  USE OpticalGlauber_Geometry
  USE ElasticPhotonPhotonFlux
  USE blocks_ref
  IMPLICIT NONE
  include 'run90.inc'
  CHARACTER(len=200)::label='blocks',code_label='new',points_file='',out_file='blocks.json'
  INTEGER::icoll=3,A1=208,Z1=82,A2=208,Z2=82,tag1=-2,tag2=-2,iepsilon=0
  REAL(KIND(1d0))::E1=2510d0,E2=2510d0
  LOGICAL::chff=.FALSE.,conj_sum=.FALSE.,mc_glauber=.FALSE.
  INTEGER::eval_kinds=8,nW=0,nlookup=0,nWdb=0,algorithm=-1
  LOGICAL::do_grid=.FALSE.,do_eval=.TRUE.,do_pnohad1=.FALSE.,deltab_grid=.FALSE.
  REAL(KIND(1d0))::Wmin=1d0,Wmax=1d0,Wdbmin=1d0,Wdbmax=1d0
  NAMELIST /config/ label,code_label,points_file,out_file,icoll,A1,Z1,A2,Z2,E1,E2,chff,&
       conj_sum,mc_glauber,tag1,tag2,iepsilon,do_grid,do_eval,do_pnohad1,eval_kinds,&
       nlookup,nW,Wmin,Wmax,nWdb,Wdbmin,Wdbmax,algorithm,deltab_grid
  CHARACTER(len=400)::arg
  CHARACTER(len=7)::aname
  INTEGER::iu,i,j,ier,nb
  REAL(KIND(1d0)),DIMENSION(9)::xs=(/1d-8,1d-6,1d-4,1d-3,3d-3,1d-2,3d-2,1d-1,0.3d0/)
  REAL(KIND(1d0)),DIMENSION(15)::bR=(/1d-3,1d-2,0.1d0,0.3d0,0.5d0,0.8d0,1d0,1.1d0,1.2d0,1.5d0,&
       2d0,3d0,5d0,10d0,40d0/)
  REAL(KIND(1d0))::anum,znum,rr,adif,wpar,b,v,pl,pe,res,err,mx,pts(60),gt,fac,alpha
  CHARACTER(len=1)::sep
  CALL GET_COMMAND_ARGUMENT(1,arg)
  OPEN(NEWUNIT=iu,FILE=TRIM(arg),STATUS='OLD')
  READ(iu,NML=config)
  CLOSE(iu)
  ebeam(1)=E1
  ebeam(2)=E2
  nuclearA_beam1=MERGE(0,A1,A1.EQ.1.AND.Z1.EQ.1)
  nuclearZ_beam1=MERGE(0,Z1,A1.EQ.1.AND.Z1.EQ.1)
  nuclearA_beam2=MERGE(0,A2,A2.EQ.1.AND.Z2.EQ.1)
  nuclearZ_beam2=MERGE(0,Z2,A2.EQ.1.AND.Z2.EQ.1)
  USE_CHARGEFORMFACTOR4PHOTON=chff
  neutron_tagging(1)=tag1
  neutron_tagging(2)=tag2
  neutron_conj_sum=conj_sum
  use_MC_Glauber=mc_glauber
  IEPSILON_EDFF=iepsilon
  alpha=0.0072992700729927005d0
  IF(icoll.EQ.3)THEN
     aname=GetASymbol(A1,Z1)
  ELSE
     aname=GetASymbol(MAX(A1,A2),MERGE(Z1,Z2,A1.GT.A2))
  ENDIF
  CALL GetNuclearInfo(aname,anum,znum,rr,adif,wpar)
  OPEN(NEWUNIT=iu,FILE=TRIM(out_file),STATUS='REPLACE')
  WRITE(iu,'(A)')'{"schema": "gupc-blocks-1", "label": "'//TRIM(label)//'",'
  WRITE(iu,'(A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A)')' "R_fm": ',rr,', "a_fm": ',adif,', "w": ',wpar,','
  ! That_AB / That_A and P_had on b = 0.25 fm steps
  WRITE(iu,'(A)')' "tab": ['
  nb=INT((2d0*rr+40d0*adif)/0.25d0)
  DO i=0,nb
     b=0.25d0*DBLE(i)
     sep=MERGE(',',' ',i.LT.nb)
     WRITE(iu,'(A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A)')'  [',b,', ',FFE_Test_Value(1,icoll,1d-3,b),', ',&
          FFE_Test_Value(2,icoll,1d-3,b),']'//sep
  ENDDO
  WRITE(iu,'(A)')' ],'
  ! photon number density of engine beam 1 (b in GeV^-1, R in GeV^-1 = rr/0.1973)
  WRITE(iu,'(A)')' "N": ['
  DO i=1,9
     DO j=1,15
        b=bR(j)*rr/0.1973d0
        IF(icoll.EQ.2)b=bR(j)*0.877d0/0.1973d0
        v=FFE_Test_Value(3,icoll,xs(i),b)
        sep=MERGE(',',' ',.NOT.(i.EQ.9.AND.j.EQ.15))
        WRITE(iu,'(A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A)')'  [',xs(i),', ',b,', ',v,']'//sep
     ENDDO
  ENDDO
  WRITE(iu,'(A)')' ],'
  ! neutron tagging
  WRITE(iu,'(A)')' "ptag": ['
  IF(tag1.NE.-2.OR.tag2.NE.-2)THEN
     gt=E2/0.9315d0
     fac=gt+DSQRT(gt*gt-1d0)
     pref_c=0.9315d0*(E1/0.9315d0)*fac
     mp_c=0.9315d0
     gp_c=E1/0.9315d0
     at_c=NINT(anum)
     zt_c=NINT(znum)
     DO i=0,40
        b=(1d1/0.1973d0)*10d0**(DBLE(i)/5d0)
        pe=FFE_Test_Value(5,icoll,1d-3,b)
        pl=Pbreak_AB_WoodsSaxon(tag1,tag2,b,conj_sum)
        mx=-1d0
        IF(.NOT.chff)THEN
           ! independent m_Xn(b) of beam 1 -> beam 2 (Emin_Xn = 5 MeV, up to x = 1)
           bb_c=b
           DO j=1,60
              pts(j)=DLOG(5d-3)+DBLE(j)*0.06d0
           ENDDO
           CALL gk21_adapt(f_mxn,DLOG(5d-3),DLOG(pref_c),0d0,1d-11,res,err,ier,60,pts)
           mx=res*znum**2*alpha
        ENDIF
        sep=MERGE(',',' ',i.LT.40)
        WRITE(iu,'(A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A)')'  [',b,', ',pe,', ',pl,', ',mx,']'//sep
     ENDDO
  ENDIF
  WRITE(iu,'(A)')' ]}'
  CLOSE(iu)
END PROGRAM bench_blocks
