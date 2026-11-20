! Benchmark driver for the gamma-UPC photon-photon flux grids (AB and pA).
! It uses only the public API, so the same source builds against the legacy
! library (commit 6f68064) and the new one (-DGUPC_NEW enables new-only switches).
!
! Usage: bench_ab <namelist file>
!   The namelist /config/ (see below) sets the beams, flags and what to measure;
!   points_file lists "x1 x2 kind" lines (kind: 0 node, 1 cell midpoint, 2 random,
!   3 reference point). The results are written as JSON to out_file.
PROGRAM bench_ab
  USE OpticalGlauber_Geometry
  USE ElasticPhotonPhotonFlux
  IMPLICIT NONE
  include 'run90.inc'
  ! ---------------------------------------------------------------- config --
  CHARACTER(len=200)::label='run',code_label='legacy',points_file='',out_file='bench.json'
  INTEGER::icoll=3            ! 1 pp, 2 pA/Ap, 3 AB
  INTEGER::A1=208,Z1=82,A2=208,Z2=82
  REAL(KIND(1d0))::E1=2510d0,E2=2510d0 ! per nucleon (or proton) energies in GeV
  LOGICAL::chff=.FALSE.,conj_sum=.FALSE.,mc_glauber=.FALSE.
  INTEGER::tag1=-2,tag2=-2,iepsilon=0
  LOGICAL::do_grid=.TRUE.     ! time the first grid call (must be the first flux call)
  LOGICAL::do_eval=.TRUE.     ! direct (non-grid) evaluation at the points
  LOGICAL::do_pnohad1=.TRUE.  ! PNOHAD=1 reference values at the points
  INTEGER::eval_kinds=8       ! bit mask of point kinds evaluated directly (1:node 2:mid 4:random 8:ref)
  INTEGER::nlookup=20000      ! lookups for the mean lookup time
  INTEGER::nW=40              ! dL/dW points (log-spaced in [Wmin,Wmax])
  REAL(KIND(1d0))::Wmin=1d0,Wmax=2000d0
  INTEGER::nWdb=0             ! DeltaB_UPC_at_W points (log-spaced in [Wdbmin,Wdbmax])
  REAL(KIND(1d0))::Wdbmin=2d0,Wdbmax=200d0
  INTEGER::algorithm=-1       ! new code only: PhotonPhotonFlux_Algorithm (-1 = library default)
  LOGICAL::deltab_grid=.FALSE. ! GENERATE_DELTAB_GRID
  LOGICAL::dbpoints=.FALSE.    ! evaluate the <Delta b> numerator (direct) instead of the flux at the points
  LOGICAL::legacy_defs=.FALSE. ! new code only: FFE_Legacy_ChFF_Definitions
  REAL(KIND(1d0))::node_scale=1d0 ! new code only: FFE_NodeScale (convergence studies)
  INTEGER::verbose=0           ! new code only: FFE_Verbose (stage timings)
  LOGICAL::grid_cache=.FALSE.  ! new code only: FFE_Grid_Cache
  CHARACTER(len=400)::grid_cache_dir='./' ! new code only: FFE_Grid_Cache_Dir
  NAMELIST /config/ label,code_label,points_file,out_file,icoll,A1,Z1,A2,Z2,E1,E2,chff,&
       conj_sum,mc_glauber,tag1,tag2,iepsilon,do_grid,do_eval,do_pnohad1,eval_kinds,&
       nlookup,nW,Wmin,Wmax,nWdb,Wdbmin,Wdbmax,algorithm,deltab_grid,dbpoints,legacy_defs,&
       node_scale,verbose,grid_cache,grid_cache_dir
  ! ----------------------------------------------------------------- state --
  CHARACTER(len=400)::arg
  INTEGER::iu,ios,np,i,kind_i
  INTEGER,PARAMETER::NPMAX=20000
  REAL(KIND(1d0)),DIMENSION(NPMAX)::px1,px2,zgrid,zeval,zp1
  INTEGER,DIMENSION(NPMAX)::pkind
  LOGICAL,DIMENSION(NPMAX)::have_eval
  REAL(KIND(1d0))::t_grid,t_grid_wall,t_lookup,t_eval_sum,t_dldw,t_db,z,u1,u2,tc0,tc1
  INTEGER::n_eval
  REAL(KIND(1d0))::t_first_eval
  INTEGER(KIND=8)::c0,c1,crate
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::Wl,dldw,Wdb,db
  REAL(KIND(1d0))::xr1,xr2,acc

  CALL GET_COMMAND_ARGUMENT(1,arg)
  IF(LEN_TRIM(arg).EQ.0)THEN
     WRITE(*,*)"usage: bench_ab <namelist>"
     STOP 1
  ENDIF
  OPEN(NEWUNIT=iu,FILE=TRIM(arg),STATUS='OLD')
  READ(iu,NML=config)
  CLOSE(iu)
  CALL SYSTEM_CLOCK(COUNT_RATE=crate)

  ! beams: ebeam is the energy per nucleon (proton energy for a proton beam)
  ebeam(1)=E1
  ebeam(2)=E2
  IF(A1.EQ.1.AND.Z1.EQ.1)THEN
     nuclearA_beam1=0
     nuclearZ_beam1=0
  ELSE
     nuclearA_beam1=A1
     nuclearZ_beam1=Z1
  ENDIF
  IF(A2.EQ.1.AND.Z2.EQ.1)THEN
     nuclearA_beam2=0
     nuclearZ_beam2=0
  ELSE
     nuclearA_beam2=A2
     nuclearZ_beam2=Z2
  ENDIF
  USE_CHARGEFORMFACTOR4PHOTON=chff
  neutron_tagging(1)=tag1
  neutron_tagging(2)=tag2
  neutron_conj_sum=conj_sum
  use_MC_Glauber=mc_glauber
  IEPSILON_EDFF=iepsilon
  GENERATE_DELTAB_GRID=deltab_grid
#ifdef GUPC_NEW
  IF(algorithm.GE.0)PhotonPhotonFlux_Algorithm=algorithm
  FFE_Legacy_ChFF_Definitions=legacy_defs
  FFE_NodeScale=node_scale
  FFE_Verbose=verbose
  FFE_Grid_Cache=grid_cache
  FFE_Grid_Cache_Dir=grid_cache_dir
#endif

  ! points
  np=0
  IF(LEN_TRIM(points_file).GT.0)THEN
     OPEN(NEWUNIT=iu,FILE=TRIM(points_file),STATUS='OLD')
     DO
        READ(iu,*,IOSTAT=ios)xr1,xr2,kind_i
        IF(ios.NE.0)EXIT
        np=np+1
        IF(np.GT.NPMAX)STOP "too many points"
        px1(np)=xr1
        px2(np)=xr2
        pkind(np)=kind_i
     ENDDO
     CLOSE(iu)
  ENDIF

  ! ----------------------------------------------- first call: grid build --
  t_grid=-1d0
  t_grid_wall=-1d0
  IF(do_grid)THEN
     CALL SYSTEM_CLOCK(c0)
     CALL CPU_TIME(tc0)
     z=flux_grid(1d-3,2d-3)
     CALL CPU_TIME(tc1)
     CALL SYSTEM_CLOCK(c1)
     t_grid=tc1-tc0
     t_grid_wall=DBLE(c1-c0)/DBLE(crate)
     WRITE(*,*)"BENCH: grid build time [s] cpu, wall =",t_grid,t_grid_wall
  ENDIF

  ! ------------------------------------------------------- lookup timing --
  t_lookup=-1d0
  IF(do_grid.AND.nlookup.GT.0)THEN
     acc=0d0
     CALL CPU_TIME(tc0)
     DO i=1,nlookup
        ! deterministic pseudo-random points in [1e-6,1e-1]^2
        u1=MODULO(0.6180339887498949d0*DBLE(i),1d0)
        u2=MODULO(0.7548776662466927d0*DBLE(i),1d0)
        acc=acc+flux_grid(10d0**(-6d0+5d0*u1),10d0**(-6d0+5d0*u2))
     ENDDO
     CALL CPU_TIME(tc1)
     t_lookup=(tc1-tc0)/DBLE(nlookup)
     IF(acc.LT.0d0)WRITE(*,*)"impossible",acc
  ENDIF

  ! --------------------------------- first direct evaluation (setup cost) --
  t_first_eval=-1d0
  IF(do_eval.AND.np.GT.0)THEN
     CALL CPU_TIME(tc0)
     z=flux_eval(1.234d-3,2.345d-3,.FALSE.)
     CALL CPU_TIME(tc1)
     t_first_eval=tc1-tc0
  ENDIF

  ! --------------------------------------------------------- the points --
  t_eval_sum=0d0
  n_eval=0
  DO i=1,np
     have_eval(i)=.FALSE.
     zgrid(i)=-1d0
     zeval(i)=-1d0
     zp1(i)=-1d0
     IF(do_grid)zgrid(i)=flux_grid(px1(i),px2(i))
     IF(do_eval.AND.IAND(eval_kinds,ISHFT(1,pkind(i))).NE.0)THEN
        CALL CPU_TIME(tc0)
        zeval(i)=flux_eval(px1(i),px2(i),.FALSE.)
        CALL CPU_TIME(tc1)
        t_eval_sum=t_eval_sum+(tc1-tc0)
        n_eval=n_eval+1
        have_eval(i)=.TRUE.
     ENDIF
     IF(do_pnohad1)zp1(i)=flux_eval(px1(i),px2(i),.TRUE.)
  ENDDO

  ! --------------------------------------------------------- dL/dW(W) --
  t_dldw=-1d0
  IF(nW.GT.0)THEN
     ALLOCATE(Wl(nW),dldw(nW))
     CALL CPU_TIME(tc0)
     DO i=1,nW
        IF(nW.EQ.1)THEN
           Wl(i)=Wmin
        ELSE
           Wl(i)=Wmin*(Wmax/Wmin)**(DBLE(i-1)/DBLE(nW-1))
        ENDIF
        dldw(i)=dLgammagammadW_UPC(Wl(i),icoll,1)
     ENDDO
     CALL CPU_TIME(tc1)
     t_dldw=(tc1-tc0)/DBLE(nW)
  ENDIF

  ! ---------------------------------------------------- <Delta b>(W) --
  t_db=-1d0
  IF(nWdb.GT.0)THEN
     ALLOCATE(Wdb(nWdb),db(nWdb))
     CALL CPU_TIME(tc0)
     DO i=1,nWdb
        IF(nWdb.EQ.1)THEN
           Wdb(i)=Wdbmin
        ELSE
           Wdb(i)=Wdbmin*(Wdbmax/Wdbmin)**(DBLE(i-1)/DBLE(nWdb-1))
        ENDIF
        db(i)=DeltaB_UPC_at_W(Wdb(i),icoll,1)
     ENDDO
     CALL CPU_TIME(tc1)
     t_db=(tc1-tc0)/DBLE(nWdb)
  ENDIF

  ! ------------------------------------------------------------- output --
  OPEN(NEWUNIT=iu,FILE=TRIM(out_file),STATUS='REPLACE')
  WRITE(iu,'(A)')'{'
  WRITE(iu,'(A)')'  "schema": "gupc-bench-1",'
  WRITE(iu,'(A)')'  "label": "'//TRIM(label)//'",'
  WRITE(iu,'(A)')'  "code": "'//TRIM(code_label)//'",'
  WRITE(iu,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)')'  "config": {"icoll": ',icoll,', "A1": ',A1,', "Z1": ',Z1,&
       ', "A2": ',A2,', "Z2": ',Z2,','
  WRITE(iu,'(A,A,A,A,A)')'    "E1": ',jnum(E1),', "E2": ',jnum(E2),','
  WRITE(iu,'(A,A,A,I0,A,I0,A,A,A,I0,A,A,A)')'    "chff": ',jlog(chff),', "tag1": ',tag1,', "tag2": ',tag2,&
       ', "conj_sum": ',jlog(conj_sum),', "iepsilon": ',iepsilon,', "mc_glauber": ',jlog(mc_glauber),'},'
  WRITE(iu,'(A,A,A,A,A)')'  "timing": {"clock": "cpu", "grid_wall_s": ',jnum(t_grid_wall),','
  WRITE(iu,'(A,A,A,A,A,A,A,I0,A)')'    "grid_s": ',jnum(t_grid),', "lookup_s": ',jnum(t_lookup),&
       ', "eval_mean_s": ',jnum(merge(t_eval_sum/DBLE(MAX(n_eval,1)),-1d0,n_eval.GT.0)),', "n_eval": ',n_eval,','
  WRITE(iu,'(A,A,A,A,A,A,A)')'    "first_eval_s": ',jnum(t_first_eval),', "dldw_mean_s": ',&
       jnum(t_dldw),', "deltab_mean_s": ',jnum(t_db),'},'
  WRITE(iu,'(A)')'  "points": ['
  DO i=1,np
     WRITE(iu,'(A,A,A,A,A,I0,A,A,A,A,A,A,A)')'    {"x1": ',jnum(px1(i)),', "x2": ',jnum(px2(i)),&
          ', "kind": ',pkind(i),', "grid": ',jopt(zgrid(i),do_grid),', "eval": ',jopt(zeval(i),have_eval(i)),&
          ', "pnohad1": ',jopt(zp1(i),do_pnohad1),'}'//TRIM(MERGE(',',' ',i.LT.np))
  ENDDO
  WRITE(iu,'(A)')'  ],'
  WRITE(iu,'(A)')'  "dldw": ['
  DO i=1,nW
     WRITE(iu,'(A,A,A,A,A)')'    [',jnum(Wl(i)),', ',jnum(dldw(i)),']'//TRIM(MERGE(',',' ',i.LT.nW))
  ENDDO
  WRITE(iu,'(A)')'  ],'
  WRITE(iu,'(A)')'  "deltab": ['
  DO i=1,nWdb
     WRITE(iu,'(A,A,A,A,A)')'    [',jnum(Wdb(i)),', ',jnum(db(i)),']'//TRIM(MERGE(',',' ',i.LT.nWdb))
  ENDDO
  WRITE(iu,'(A)')'  ]'
  WRITE(iu,'(A)')'}'
  CLOSE(iu)
  WRITE(*,*)"BENCH: wrote ",TRIM(out_file)

CONTAINS

  FUNCTION flux_grid(x1,x2)
    REAL(KIND(1d0))::flux_grid
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    IF(icoll.EQ.3)THEN
       flux_grid=PhotonPhotonFlux_AB_WoodsSaxon(x1,x2)
    ELSEIF(icoll.EQ.1)THEN
       flux_grid=PhotonPhotonFlux_pp(x1,x2)
    ELSE
       flux_grid=PhotonPhotonFlux_pA_WoodsSaxon(x1,x2)
    ENDIF
  END FUNCTION flux_grid

  FUNCTION flux_eval(x1,x2,pnohad1)
    REAL(KIND(1d0))::flux_eval
    REAL(KIND(1d0)),INTENT(IN)::x1,x2
    LOGICAL,INTENT(IN)::pnohad1
    IF(dbpoints.AND..NOT.pnohad1)THEN
       IF(icoll.EQ.3)THEN
          flux_eval=PhotonPhotonDeltaB_AB_WoodsSaxon_eval(x1,x2)
       ELSEIF(icoll.EQ.1)THEN
          flux_eval=PhotonPhotonDeltaB_pp_eval(x1,x2)
       ELSE
          flux_eval=PhotonPhotonDeltaB_pA_WoodsSaxon_eval(x1,x2)
       ENDIF
       RETURN
    ENDIF
    IF(icoll.EQ.3)THEN
       flux_eval=PhotonPhotonFlux_AB_WoodsSaxon_eval(x1,x2,pnohad1)
    ELSEIF(icoll.EQ.1)THEN
       flux_eval=PhotonPhotonFlux_pp_eval(x1,x2,pnohad1)
    ELSE
       flux_eval=PhotonPhotonFlux_pA_WoodsSaxon_eval(x1,x2,pnohad1)
    ENDIF
  END FUNCTION flux_eval

  FUNCTION jnum(v)
    CHARACTER(len=34)::jnum
    REAL(KIND(1d0)),INTENT(IN)::v
    IF(v.NE.v.OR.ABS(v).GT.HUGE(v))THEN
       jnum='null'
    ELSE
       WRITE(jnum,'(ES26.17E3)')v
       jnum=ADJUSTL(jnum)
    ENDIF
  END FUNCTION jnum

  FUNCTION jopt(v,have)
    CHARACTER(len=34)::jopt
    REAL(KIND(1d0)),INTENT(IN)::v
    LOGICAL,INTENT(IN)::have
    IF(have)THEN
       jopt=jnum(v)
    ELSE
       jopt='null'
    ENDIF
  END FUNCTION jopt

  FUNCTION jlog(l)
    CHARACTER(len=5)::jlog
    LOGICAL,INTENT(IN)::l
    IF(l)THEN
       jlog='true'
    ELSE
       jlog='false'
    ENDIF
  END FUNCTION jlog

END PROGRAM bench_ab
