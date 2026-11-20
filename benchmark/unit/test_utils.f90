! Unit tests of the numerical utilities of gamma-UPC (milestone M1):
!   tbessk.f90      BESSK0/K1/I0/I1, GUPC_BESSK0E/K1E/I0E/I1E, BESSK(N,X)
!   nintlib.f90     gl_rule, composite_gl, gk21_rule, gk21_adapt
!   interpolation   lagr_uniform_weights, lagr_uniform_eval1d, lagr_uniform_eval2d
! Reads bessel_ref.dat and besselkn_ref.dat (mpmath, gen_refs.py) from the
! current directory and writes results_utils.json there. Exit status 1 if a
! test fails.
MODULE ut_report
  IMPLICIT NONE
  INTEGER,PARAMETER::maxres=200
  INTEGER::nres=0,nfail=0
  CHARACTER(len=48),DIMENSION(maxres)::r_id
  CHARACTER(len=200),DIMENSION(maxres)::r_desc
  CHARACTER(len=8),DIMENSION(maxres)::r_status
  REAL(KIND(1d0)),DIMENSION(maxres)::r_err,r_tol
CONTAINS
  SUBROUTINE add_result(id,desc,err,tol)
    ! tol<0: informational entry (status "info", tolerance null)
    CHARACTER(len=*),INTENT(IN)::id,desc
    REAL(KIND(1d0)),INTENT(IN)::err,tol
    nres=nres+1
    r_id(nres)=id
    r_desc(nres)=desc
    r_err(nres)=err
    r_tol(nres)=tol
    IF(tol.LT.0d0)THEN
       r_status(nres)='info'
    ELSEIF(err.LE.tol)THEN   ! NaN fails
       r_status(nres)='pass'
    ELSE
       r_status(nres)='FAIL'
       nfail=nfail+1
    ENDIF
    WRITE(*,'(A4,1X,A32,1X,ES10.3,1X,A,ES9.2,2X,A)')r_status(nres),id,err,'tol ',tol,TRIM(desc)
  END SUBROUTINE add_result

  FUNCTION fmt_num(v) RESULT(s)
    REAL(KIND(1d0)),INTENT(IN)::v
    CHARACTER(len=24)::s
    IF(v.NE.v)THEN
       s='null'
    ELSE
       WRITE(s,'(ES12.4E3)')v
       s=ADJUSTL(s)
    ENDIF
  END FUNCTION fmt_num

  SUBROUTINE write_json(fname)
    CHARACTER(len=*),INTENT(IN)::fname
    INTEGER::iu,i
    CHARACTER(len=24)::tols
    CHARACTER(len=1)::sep
    OPEN(newunit=iu,file=fname,status='replace')
    WRITE(iu,'(A)')'['
    DO i=1,nres
       IF(r_tol(i).LT.0d0)THEN
          tols='null'
       ELSE
          tols=fmt_num(r_tol(i))
       ENDIF
       sep=','
       IF(i.EQ.nres)sep=' '
       WRITE(iu,'(12A)')'  {"id": "',TRIM(r_id(i)),'", "description": "',TRIM(r_desc(i)),&
            '", "max_error": ',TRIM(fmt_num(r_err(i))),', "tolerance": ',TRIM(tols),&
            ', "status": "',TRIM(r_status(i)),'"}',TRIM(sep)
    ENDDO
    WRITE(iu,'(A)')']'
    CLOSE(iu)
  END SUBROUTINE write_json
END MODULE ut_report

MODULE ut_funcs
  ! integrands; the nested test passes the outer variables through module
  ! variables, as the library does
  IMPLICIT NONE
  REAL(KIND(1d0))::xo=0d0,yo=0d0
  INTEGER::nest_ier=0
  REAL(KIND(1d0)),PARAMETER::nest_eps=1d-13
CONTAINS
  FUNCTION f_exp(x)
    REAL(KIND(1d0))::f_exp,x
    f_exp=EXP(x)
  END FUNCTION f_exp
  FUNCTION f_sin(x)
    REAL(KIND(1d0))::f_sin,x
    f_sin=SIN(x)
  END FUNCTION f_sin
  FUNCTION f_log(x)
    REAL(KIND(1d0))::f_log,x
    f_log=LOG(x)
  END FUNCTION f_log
  FUNCTION f_log2(x)
    REAL(KIND(1d0))::f_log2,x
    f_log2=LOG(x)**2
  END FUNCTION f_log2
  FUNCTION f_isqrt(x)
    REAL(KIND(1d0))::f_isqrt,x
    f_isqrt=1d0/SQRT(x)
  END FUNCTION f_isqrt
  FUNCTION f_osc(x)
    REAL(KIND(1d0))::f_osc,x
    f_osc=EXP(-x)*COS(3d1*x)
  END FUNCTION f_osc
  FUNCTION f_cos2(x)
    REAL(KIND(1d0))::f_cos2,x
    f_cos2=COS(1d2*x)**2
  END FUNCTION f_cos2
  FUNCTION f_kink(x)
    REAL(KIND(1d0))::f_kink,x
    f_kink=ABS(x-1d0)*EXP(x)
  END FUNCTION f_kink
  FUNCTION f_logabs(x)
    REAL(KIND(1d0))::f_logabs,x
    f_logabs=LOG(ABS(x-1d0))
  END FUNCTION f_logabs
  FUNCTION f_pow(x)
    ! x**k with k=NINT(xo)
    REAL(KIND(1d0))::f_pow,x
    f_pow=x**NINT(xo)
  END FUNCTION f_pow
  ! nested: int_0^1 dx int_0^x dy int_0^y dz exp(z)*(1+x*y) = 2e-41/8
  FUNCTION f_inner(z)
    REAL(KIND(1d0))::f_inner,z
    f_inner=EXP(z)*(1d0+xo*yo)
  END FUNCTION f_inner
  FUNCTION f_mid(y)
    USE NINTLIB
    REAL(KIND(1d0))::f_mid,y,r,e
    INTEGER::ier
    yo=y
    CALL gk21_adapt(f_inner,0d0,y,0d0,nest_eps,r,e,ier)
    IF(ier.NE.0)nest_ier=ier
    f_mid=r
  END FUNCTION f_mid
  FUNCTION f_outer(x)
    USE NINTLIB
    REAL(KIND(1d0))::f_outer,x,r,e
    INTEGER::ier
    xo=x
    CALL gk21_adapt(f_mid,0d0,x,0d0,nest_eps,r,e,ier)
    IF(ier.NE.0)nest_ier=ier
    f_outer=r
  END FUNCTION f_outer
END MODULE ut_funcs

PROGRAM test_utils
  USE ut_report
  USE ut_funcs
  USE NINTLIB
  USE interpolation
  IMPLICIT NONE
  REAL(KIND(1d0)),PARAMETER::pi=3.14159265358979323846264338328d0
  REAL(KIND(1d0)),PARAMETER::e1=2.71828182845904523536028747135d0

  CALL test_bessel()
  CALL test_gl()
  CALL test_composite()
  CALL test_gk21()
  CALL test_lagrange()
  CALL write_json('results_utils.json')
  WRITE(*,'(A,I0,A,I0,A)')'tests: ',nres,' entries, ',nfail,' failed'
  IF(nfail.GT.0)STOP 1

CONTAINS

  ! ------------------------------------------------------------- Bessel --
  SUBROUTINE test_bessel()
    REAL(KIND(1d0)),EXTERNAL::BESSK0,BESSK1,BESSI0,BESSI1,BESSK
    REAL(KIND(1d0)),EXTERNAL::GUPC_BESSK0E,GUPC_BESSK1E,GUPC_BESSI0E,GUPC_BESSI1E
    REAL(KIND(1d0)),EXTERNAL::LEGACY_BESSK0,LEGACY_BESSK1,LEGACY_BESSI0,LEGACY_BESSI1
    CHARACTER(len=256)::hdr
    CHARACTER(len=200)::desc
    CHARACTER(len=8),DIMENSION(8),PARAMETER::nm=(/'k0      ','k1      ','i0      ','i1      ',&
         'k0e     ','k1e     ','i0e     ','i1e     '/)
    CHARACTER(len=24),DIMENSION(8),PARAMETER::fn=(/'BESSK0                  ',&
         'BESSK1                  ','BESSI0                  ','BESSI1                  ',&
         'GUPC_BESSK0E = e^x K0   ','GUPC_BESSK1E = e^x K1   ','GUPC_BESSI0E = e^-x I0  ',&
         'GUPC_BESSI1E = e^-x I1  '/)
    INTEGER::iu,n,i,k,nn,nskip(8)
    REAL(KIND(1d0))::x,ref(8),v(8),e,emax(8),xat(8),elmax(4),xlat(4),vl(4)
    REAL(KIND(1d0))::dmax(4),xdat(4),dmaxs(4),xdats(4),lx,vn(4)
    REAL(KIND(1d0))::t0,t1,s,tnew,told
    INTEGER::j
    OPEN(newunit=iu,file='bessel_ref.dat',status='old')
    READ(iu,'(A)')hdr
    READ(hdr(2:),*)n
    emax=0d0
    xat=0d0
    elmax=0d0
    xlat=0d0
    nskip=0
    DO i=1,n
       READ(iu,*)x,ref
       v=(/BESSK0(x),BESSK1(x),BESSI0(x),BESSI1(x),GUPC_BESSK0E(x),GUPC_BESSK1E(x),&
            GUPC_BESSI0E(x),GUPC_BESSI1E(x)/)
       vl=(/LEGACY_BESSK0(x),LEGACY_BESSK1(x),LEGACY_BESSI0(x),LEGACY_BESSI1(x)/)
       DO k=1,8
          ! unscaled functions only where they are normal numbers
          IF(ABS(ref(k)).LT.1d-300.OR.ABS(ref(k)).GT.1d300)THEN
             nskip(k)=nskip(k)+1
             CYCLE
          ENDIF
          e=ABS(v(k)-ref(k))/ABS(ref(k))
          IF(.NOT.(e.LE.emax(k)))THEN
             emax(k)=e
             xat(k)=x
          ENDIF
          IF(k.LE.4)THEN
             e=ABS(vl(k)-ref(k))/ABS(ref(k))
             IF(e.GT.elmax(k))THEN
                elmax(k)=e
                xlat(k)=x
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CLOSE(iu)
    DO k=1,8
       WRITE(desc,'(A,A,I0,A,ES9.3,A,I0,A)')TRIM(fn(k)),' vs mpmath at ',n-nskip(k),&
            ' log-spaced x in [1e-12,700], max at x=',xat(k),' (',nskip(k),' skipped)'
       CALL add_result('bessel_'//TRIM(nm(k)),desc,emax(k),2d-15)
    ENDDO
    ! K_n by recursion
    OPEN(newunit=iu,file='besselkn_ref.dat',status='old')
    READ(iu,'(A)')hdr
    READ(hdr(2:),*)n
    e=0d0
    DO i=1,n
       READ(iu,*)nn,x,ref(1)
       e=MAX(e,ABS(BESSK(nn,x)-ref(1))/ABS(ref(1)))
    ENDDO
    CLOSE(iu)
    CALL add_result('bessel_kn_recursion','BESSK(N,X) vs mpmath, N=2..6, x in {1e-3,0.125,1.2097,5,30,200}',&
         e,1d-14)
    ! legacy (A&S) vs mpmath on the same points
    DO k=1,4
       WRITE(desc,'(A,A,ES9.3)')'legacy (A&S fit) '//TRIM(fn(k)),' vs mpmath, 300 points, max at x=',xlat(k)
       CALL add_result('legacy_bessel_'//TRIM(nm(k))//'_error',desc,elmax(k),-1d0)
    ENDDO
    ! legacy vs new, dense log grid: full range and x<=10 (where photon
    ! densities are non-negligible)
    dmax=0d0
    dmaxs=0d0
    xdat=0d0
    xdats=0d0
    DO i=0,200000
       lx=-12d0+(LOG10(700d0)+12d0)*DBLE(i)/2d5
       x=10d0**lx
       vn=(/BESSK0(x),BESSK1(x),BESSI0(x),BESSI1(x)/)
       vl=(/LEGACY_BESSK0(x),LEGACY_BESSK1(x),LEGACY_BESSI0(x),LEGACY_BESSI1(x)/)
       DO k=1,4
          e=ABS(vl(k)-vn(k))/ABS(vn(k))
          IF(e.GT.dmax(k))THEN
             dmax(k)=e
             xdat(k)=x
          ENDIF
          IF(x.LE.10d0.AND.e.GT.dmaxs(k))THEN
             dmaxs(k)=e
             xdats(k)=x
          ENDIF
       ENDDO
    ENDDO
    DO k=1,4
       WRITE(desc,'(A,ES9.3,A,ES9.3,A,ES9.3)')'max |old-new|/new of '//TRIM(fn(k))//&
            ', 200001 log-spaced x in [1e-12,700], at x=',xdat(k),'; for x<=10: ',dmaxs(k),' at x=',xdats(k)
       CALL add_result('legacy_vs_new_'//TRIM(nm(k)),desc,dmax(k),-1d0)
    ENDDO
    ! timing (stdout only)
    s=0d0
    CALL CPU_TIME(t0)
    DO j=1,20
       DO i=1,100000
          x=1d-3+DBLE(i)*1d-4
          s=s+BESSK0(x)+BESSK1(x)
       ENDDO
    ENDDO
    CALL CPU_TIME(t1)
    tnew=(t1-t0)/4d6*1d9
    CALL CPU_TIME(t0)
    DO j=1,20
       DO i=1,100000
          x=1d-3+DBLE(i)*1d-4
          s=s+LEGACY_BESSK0(x)+LEGACY_BESSK1(x)
       ENDDO
    ENDDO
    CALL CPU_TIME(t1)
    told=(t1-t0)/4d6*1d9
    WRITE(*,'(A,F7.1,A,F7.1,A,ES10.3,A)')'timing  K0/K1 per call, x in [1e-3,10]: new ',tnew,&
         ' ns, legacy ',told,' ns  (checksum ',s,')'
  END SUBROUTINE test_bessel

  ! --------------------------------------------------- Gauss-Legendre --
  SUBROUTINE test_gl()
    REAL(KIND(1d0)),DIMENSION(128)::x,w,x2,w2
    REAL(KIND(1d0))::emax,emaxu,q,ex,esym
    INTEGER::n,k,i,nlist(3)
    CHARACTER(len=200)::desc
    emax=0d0
    esym=0d0
    DO n=1,64
       CALL gl_rule(n,x(1:n),w(1:n))
       DO k=0,2*n-1
          q=0d0
          DO i=1,n
             q=q+w(i)*x(i)**k
          ENDDO
          IF(MOD(k,2).EQ.0)THEN
             ex=2d0/DBLE(k+1)
          ELSE
             ex=0d0
          ENDIF
          emax=MAX(emax,ABS(q-ex)/(2d0/DBLE(k+1)))
       ENDDO
       ! ascending, antisymmetric, positive weights; cached copy identical
       DO i=1,n
          IF(x(i).NE.-x(n+1-i).OR.w(i).NE.w(n+1-i).OR.w(i).LE.0d0)esym=1d0
          IF(i.GT.1)THEN
             IF(x(i).LE.x(i-1))esym=1d0
          ENDIF
       ENDDO
       CALL gl_rule(n,x2(1:n),w2(1:n))
       IF(ANY(x2(1:n).NE.x(1:n)).OR.ANY(w2(1:n).NE.w(1:n)))esym=1d0
    ENDDO
    ! the correctly rounded exact rules give 3.6e-15 on this measure (node
    ! rounding: ~(k+1)*eps/2 for x^k), numpy leggauss 4.3e-13
    CALL add_result('gl_exactness','gl_rule n=1..64: max |Q(x^k)-I(x^k)|/I(|x|^k), k=0..2n-1, on [-1,1]',&
         emax,2d-14)
    CALL add_result('gl_structure','gl_rule n=1..64: ascending, antisymmetric nodes, positive symmetric '//&
         'weights, cached rule identical (1 = violation)',esym,0d0)
    nlist=(/65,100,128/)
    emaxu=0d0
    DO i=1,3
       n=nlist(i)
       CALL gl_rule(n,x(1:n),w(1:n))
       DO k=0,2*n-1,2
          q=SUM(w(1:n)*x(1:n)**k)
          emaxu=MAX(emaxu,ABS(q-2d0/DBLE(k+1))/(2d0/DBLE(k+1)))
       ENDDO
    ENDDO
    WRITE(desc,'(A)')'gl_rule n=65,100,128 (not cached): max relative error of Q(x^k), k even <= 2n-1'
    ! node rounding alone gives ~(k+1)*eps/2 for x^k (2.8e-14 at k=255)
    CALL add_result('gl_exactness_uncached',desc,emaxu,5d-14)
  END SUBROUTINE test_gl

  SUBROUTINE test_composite()
    REAL(KIND(1d0)),DIMENSION(200)::xs,ws
    REAL(KIND(1d0))::q,ex,err
    INTEGER::nout,i
    ! scalar npts, kink at a breakpoint
    CALL composite_gl(3,(/0d0,1d0,3d0/),16,200,xs,ws,nout)
    q=0d0
    DO i=1,nout
       q=q+ws(i)*f_kink(xs(i))
    ENDDO
    ex=e1**3+2d0*e1-2d0
    err=ABS(q-ex)/ABS(ex)
    IF(nout.NE.32)err=1d0
    CALL add_result('composite_gl_kink','composite_gl, 16 pts on [0,1],[1,3]: int_0^3 |x-1| e^x dx = e^3+2e-2',&
         err,1d-13)
    ! vector npts
    CALL composite_gl(5,(/0d0,1d0,2.5d0,5d0,1d1/),(/12,16,20,24/),200,xs,ws,nout)
    q=0d0
    DO i=1,nout
       q=q+ws(i)*EXP(-xs(i))*COS(3d0*xs(i))
    ENDDO
    ex=(EXP(-1d1)*(3d0*SIN(3d1)-COS(3d1))+1d0)/1d1
    err=ABS(q-ex)/ABS(ex)
    IF(nout.NE.72)err=1d0
    CALL add_result('composite_gl_vector','composite_gl, npts=(12,16,20,24) on 0,1,2.5,5,10: '//&
         'int e^-x cos(3x) dx',err,1d-13)
    ! zero-length interval skipped, reversed interval gives the oriented integral
    CALL composite_gl(4,(/0d0,1d0,1d0,2d0/),10,200,xs,ws,nout)
    q=0d0
    DO i=1,nout
       q=q+ws(i)*EXP(xs(i))
    ENDDO
    err=ABS(q-(e1**2-1d0))/(e1**2-1d0)
    IF(nout.NE.20)err=1d0
    CALL composite_gl(2,(/1d0,0d0/),10,200,xs,ws,nout)
    q=0d0
    DO i=1,nout
       q=q+ws(i)*EXP(xs(i))
    ENDDO
    err=MAX(err,ABS(q+(e1-1d0))/(e1-1d0))
    CALL add_result('composite_gl_degenerate','composite_gl: zero-length interval skipped (nout=20), '//&
         'reversed interval = -int',err,1d-14)
  END SUBROUTINE test_composite

  ! ------------------------------------------------------ Gauss-Kronrod --
  SUBROUTINE gk_case(id,desc,f,a,b,epsabs,epsrel,ex,tol,ratio,np,pts)
    CHARACTER(len=*),INTENT(IN)::id,desc
    REAL(KIND(1d0)),EXTERNAL::f
    REAL(KIND(1d0)),INTENT(IN)::a,b,epsabs,epsrel,ex,tol
    REAL(KIND(1d0)),INTENT(INOUT)::ratio
    INTEGER,INTENT(IN),OPTIONAL::np
    REAL(KIND(1d0)),DIMENSION(*),INTENT(IN),OPTIONAL::pts
    REAL(KIND(1d0))::res,abserr,err
    INTEGER::ier,nev
    CHARACTER(len=200)::d2
    IF(PRESENT(np))THEN
       CALL gk21_adapt(f,a,b,epsabs,epsrel,res,abserr,ier,np,pts,neval=nev)
    ELSE
       CALL gk21_adapt(f,a,b,epsabs,epsrel,res,abserr,ier,neval=nev)
    ENDIF
    err=ABS(res-ex)/MAX(ABS(ex),1d-300)
    IF(ier.NE.0)THEN
       WRITE(*,*)'  ',id,': ier=',ier
       err=MAX(err,1d0)
    ENDIF
    IF(abserr.GT.0d0)ratio=MAX(ratio,ABS(res-ex)/abserr)
    WRITE(d2,'(A,A,I0,A,ES8.2,A)')desc,' [neval=',nev,', abserr=',abserr,']'
    CALL add_result(id,d2,err,tol)
  END SUBROUTINE gk_case

  SUBROUTINE test_gk21()
    REAL(KIND(1d0))::ratio,r,e,ra,rc,err,ex
    INTEGER::k,ier,nev
    ratio=0d0
    CALL gk_case('gk21_smooth_exp','gk21_adapt int_0^1 e^x = e-1, epsrel=1e-13',f_exp,0d0,1d0,0d0,1d-13,&
         e1-1d0,1d-13,ratio)
    CALL gk_case('gk21_smooth_sin','gk21_adapt int_0^pi sin x = 2, epsrel=1e-13',f_sin,0d0,pi,0d0,1d-13,&
         2d0,1d-13,ratio)
    CALL gk_case('gk21_reversed','gk21_adapt int_1^0 e^x = 1-e, epsrel=1e-13',f_exp,1d0,0d0,0d0,1d-13,&
         1d0-e1,1d-13,ratio)
    CALL gk_case('gk21_log_endpoint','gk21_adapt int_0^1 log x = -1, epsrel=1e-13',f_log,0d0,1d0,0d0,1d-13,&
         -1d0,1d-13,ratio)
    CALL gk_case('gk21_log2_endpoint','gk21_adapt int_0^1 log^2 x = 2, epsrel=1e-13',f_log2,0d0,1d0,0d0,1d-13,&
         2d0,1d-13,ratio)
    CALL gk_case('gk21_invsqrt_endpoint','gk21_adapt int_0^1 x^-1/2 = 2, epsrel=1e-13',f_isqrt,0d0,1d0,0d0,&
         1d-13,2d0,1d-13,ratio)
    ex=(EXP(-1d1)*(3d1*SIN(3d2)-COS(3d2))+1d0)/9.01d2
    ! |I| = 1.1e-3 << int|f| = 0.64: the QUADPACK error floor 50 eps int|f| = 7e-15
    ! needs epsabs above it
    CALL gk_case('gk21_oscillatory','gk21_adapt int_0^10 e^-x cos(30x), epsabs=1e-14, epsrel=1e-13',f_osc,&
         0d0,1d1,1d-14,1d-13,ex,1d-13,ratio)
    CALL gk_case('gk21_oscillatory2','gk21_adapt int_0^2pi cos^2(100x) = pi, epsrel=1e-13',f_cos2,&
         0d0,2d0*pi,0d0,1d-13,pi,1d-13,ratio)
    CALL gk_case('gk21_breakpoint_kink','gk21_adapt int_0^3 |x-1| e^x, points=(1), epsrel=1e-13',f_kink,&
         0d0,3d0,0d0,1d-13,e1**3+2d0*e1-2d0,1d-13,ratio,1,(/1d0/))
    CALL gk_case('gk21_breakpoint_log','gk21_adapt int_0^2 log|x-1| = -2, points=(1,5,-1,1), epsrel=1e-13',&
         f_logabs,0d0,2d0,0d0,1d-13,-2d0,1d-13,ratio,4,(/1d0,5d0,-1d0,1d0/))
    ! nested 3 levels through module variables (needs reentrancy)
    nest_ier=0
    CALL gk21_adapt(f_outer,0d0,1d0,0d0,nest_eps,r,e,ier,neval=nev)
    ex=2d0*e1-41d0/8d0
    err=ABS(r-ex)/ABS(ex)
    IF(ier.NE.0.OR.nest_ier.NE.0)err=MAX(err,1d0)
    CALL add_result('gk21_nested3','gk21_adapt nested 3 levels: int_0^1 dx int_0^x dy int_0^y dz '//&
         'e^z (1+xy) = 2e-41/8, epsrel=1e-13 per level',err,1d-13)
    ! single panel exactness (degree 31) and invalid input
    err=0d0
    DO k=0,31
       xo=DBLE(k)
       CALL gk21_rule(f_pow,-1d0,1d0,r,e,ra,rc)
       IF(MOD(k,2).EQ.0)THEN
          ex=2d0/DBLE(k+1)
       ELSE
          ex=0d0
       ENDIF
       err=MAX(err,ABS(r-ex)/(2d0/DBLE(k+1)))
    ENDDO
    CALL add_result('gk21_rule_exactness','gk21_rule single K21 panel on [-1,1]: max |Q(x^k)-I|/I(|x|^k), '//&
         'k=0..31',err,1d-15)
    CALL gk21_adapt(f_exp,0d0,1d0,0d0,1d-15,r,e,ier)
    err=0d0
    IF(ier.NE.6)err=1d0
    CALL add_result('gk21_invalid_tol','gk21_adapt epsabs=0, epsrel=1e-15 < 50 eps returns ier=6 '//&
         '(1 = not)',err,0d0)
    CALL add_result('gk21_error_estimate','gk21_adapt cases above: max actual error / abserr '//&
         '(<1: estimate is conservative)',ratio,1d0)
  END SUBROUTINE test_gk21

  ! ---------------------------------------------------------- Lagrange --
  FUNCTION lcg(seed) RESULT(r)
    ! deterministic uniform deviates in (0,1) (Park-Miller, no overflow in int64)
    INTEGER(KIND=8),INTENT(INOUT)::seed
    REAL(KIND(1d0))::r
    seed=MOD(16807_8*seed,2147483647_8)
    r=DBLE(seed)/2147483647d0
  END FUNCTION lcg

  FUNCTION poly(c,m,u) RESULT(v)
    INTEGER,INTENT(IN)::m
    REAL(KIND(1d0)),DIMENSION(0:m),INTENT(IN)::c
    REAL(KIND(1d0)),INTENT(IN)::u
    REAL(KIND(1d0))::v
    INTEGER::k
    v=c(m)
    DO k=m-1,0,-1
       v=v*u+c(k)
    ENDDO
  END FUNCTION poly

  SUBROUTINE test_lagrange()
    INTEGER,PARAMETER::n=25,nx=13,ny=17
    REAL(KIND(1d0)),PARAMETER::u0=-0.7d0,h=0.13d0
    REAL(KIND(1d0)),DIMENSION(n)::y
    REAL(KIND(1d0)),DIMENSION(0:9)::c
    REAL(KIND(1d0)),DIMENSION(0:9,0:9)::c2
    REAL(KIND(1d0)),DIMENSION(16)::w
    REAL(KIND(1d0)),DIMENSION(nx,ny)::z
    REAL(KIND(1d0)),DIMENSION(31)::yp
    REAL(KIND(1d0))::u,v,ex,scale,ein,eext,eclip,enode,e2d,eseg,t,ux,uy,v1,v2
    REAL(KIND(1d0))::u0x,hx,u0y,hy
    INTEGER(KIND=8)::seed
    INTEGER::p,i,j,k,l,i0,it,px,py,lo,hi,ipy
    INTEGER,DIMENSION(4),PARAMETER::pylist=(/2,3,5,7/)
    seed=12345_8
    ein=0d0
    eext=0d0
    eclip=0d0
    enode=0d0
    DO p=1,10
       DO k=0,p-1
          c(k)=2d0*lcg(seed)-1d0
       ENDDO
       DO j=1,n
          y(j)=poly(c,p-1,u0+DBLE(j-1)*h)
       ENDDO
       scale=MAXVAL(ABS(y))
       ! interior: reproduce degree p-1
       DO it=1,2000
          u=u0+DBLE(n-1)*h*lcg(seed)
          v=lagr_uniform_eval1d(u,u0,h,n,p,y)
          ein=MAX(ein,ABS(v-poly(c,p-1,u))/scale)
          ! stencil inside [1,n] and nearest-centred unless clipped
          CALL lagr_uniform_weights(u,u0,h,n,p,i0,w(1:p))
          IF(i0.LT.1.OR.i0+p-1.GT.n)eclip=1d0
          t=(u-u0)/h+1d0
          IF(i0.GT.1.AND.i0+p-1.LT.n)THEN
             IF(MOD(p,2).EQ.0)THEN
                IF(t.LT.DBLE(i0+p/2-1)-1d-12.OR.t.GT.DBLE(i0+p/2)+1d-12)eclip=1d0
             ELSE
                IF(ABS(t-DBLE(i0+(p-1)/2)).GT.0.5d0+1d-12)eclip=1d0
             ENDIF
          ENDIF
       ENDDO
       ! extrapolation up to one spacing beyond the ends, clipped stencils
       DO it=1,400
          IF(MOD(it,2).EQ.0)THEN
             u=u0-h*lcg(seed)
          ELSE
             u=u0+DBLE(n-1)*h+h*lcg(seed)
          ENDIF
          v=lagr_uniform_eval1d(u,u0,h,n,p,y)
          eext=MAX(eext,ABS(v-poly(c,p-1,u))/scale)
          CALL lagr_uniform_weights(u,u0,h,n,p,i0,w(1:p))
          IF(MOD(it,2).EQ.0.AND.i0.NE.1)eclip=1d0
          IF(MOD(it,2).EQ.1.AND.i0.NE.n-p+1)eclip=1d0
       ENDDO
       ! far outside and segment-restricted stencils
       DO it=1,400
          u=u0+(DBLE(n+20)*lcg(seed)-10d0)*h
          lo=3+MOD(it,5)
          hi=lo+p-1+MOD(it,7)
          CALL lagr_uniform_weights(u,u0,h,n,p,i0,w(1:p),lo,hi)
          IF(i0.LT.lo.OR.i0+p-1.GT.hi)eclip=1d0
          CALL lagr_uniform_weights(u,u0,h,n,p,i0,w(1:p))
          IF(i0.LT.1.OR.i0+p-1.GT.n)eclip=1d0
       ENDDO
       ! nodes: exact unit vector
       DO j=1,n
          u=u0+DBLE(j-1)*h
          v=lagr_uniform_eval1d(u,u0,h,n,p,y)
          IF(v.NE.y(j))enode=MAX(enode,MAX(ABS(v-y(j)),TINY(1d0)))
       ENDDO
    ENDDO
    CALL add_result('lagr_poly_interior','lagr_uniform_eval1d p=1..10, n=25: max |I(u)-P(u)|/max|P| for '//&
         'deg P = p-1, 2000 random u per p',ein,1d-13)
    CALL add_result('lagr_poly_extrap','lagr_uniform_eval1d p=1..10: same within one spacing outside '//&
         'the grid (clipped end stencils)',eext,1d-12)
    CALL add_result('lagr_stencil_clip','lagr_uniform_weights p=1..10: stencil inside [1,n] or [jlo,jhi], '//&
         'nearest-centred, end stencils when outside (1 = violation)',eclip,0d0)
    CALL add_result('lagr_node_exact','lagr_uniform_eval1d p=1..10 at u=u0+(j-1)h: returns y(j) exactly '//&
         '(max |I-y_j|)',enode,0d0)
    ! piecewise-uniform grid: [0,1] h=0.1 (j=1..11), [1,6] h=0.25 (j=11..31)
    eseg=0d0
    DO p=2,6
       DO k=0,p-1
          c(k)=2d0*lcg(seed)-1d0
       ENDDO
       DO j=1,11
          yp(j)=poly(c,p-1,0.1d0*DBLE(j-1))
       ENDDO
       DO j=12,31
          yp(j)=poly(c,p-1,1d0+0.25d0*DBLE(j-11))
       ENDDO
       scale=MAXVAL(ABS(yp))
       DO it=1,500
          u=6d0*lcg(seed)
          IF(u.LE.1d0)THEN
             v=lagr_uniform_eval1d(u,0d0,0.1d0,31,p,yp,1,11)
          ELSE
             v=lagr_uniform_eval1d(u,-1.5d0,0.25d0,31,p,yp,11,31)
          ENDIF
          eseg=MAX(eseg,ABS(v-poly(c,p-1,u))/scale)
       ENDDO
    ENDDO
    CALL add_result('lagr_piecewise_segments','lagr_uniform_eval1d with jlo,jhi on a two-segment grid '//&
         '(h=0.1, 0.25), p=2..6: max rel. error for deg p-1',eseg,1d-13)
    ! 2D tensor product
    e2d=0d0
    u0x=0.3d0
    hx=0.07d0
    u0y=-2d0
    hy=0.21d0
    DO px=2,6
       DO ipy=1,4
          py=pylist(ipy)
          DO l=0,py-1
             DO k=0,px-1
                c2(k,l)=2d0*lcg(seed)-1d0
             ENDDO
          ENDDO
          DO j=1,ny
             DO i=1,nx
                z(i,j)=p2(c2,px,py,u0x+DBLE(i-1)*hx,u0y+DBLE(j-1)*hy)
             ENDDO
          ENDDO
          scale=MAXVAL(ABS(z))
          DO it=1,1000
             ux=u0x+DBLE(nx-1)*hx*lcg(seed)
             uy=u0y+DBLE(ny-1)*hy*lcg(seed)
             v1=lagr_uniform_eval2d(ux,uy,u0x,hx,nx,u0y,hy,ny,px,py,z)
             e2d=MAX(e2d,ABS(v1-p2(c2,px,py,ux,uy))/scale)
             ! segment bounds x: nodes 2..nx-1, y: nodes 3..ny (points inside them)
             ux=u0x+hx*(1d0+DBLE(nx-3)*lcg(seed))
             uy=u0y+hy*(2d0+DBLE(ny-3)*lcg(seed))
             v2=lagr_uniform_eval2d(ux,uy,u0x,hx,nx,u0y,hy,ny,px,py,z,2,nx-1,3,ny)
             e2d=MAX(e2d,ABS(v2-p2(c2,px,py,ux,uy))/scale)
          ENDDO
       ENDDO
    ENDDO
    CALL add_result('lagr_2d_tensor','lagr_uniform_eval2d px=2..6, py=2,3,5,7: max |I-P|/max|P| for '//&
         'P of degree (px-1,py-1), with/without segment bounds',e2d,1d-13)
  END SUBROUTINE test_lagrange

  FUNCTION p2(c2,px,py,xx,yy)
    ! sum_{k<px,l<py} c2(k,l) x^k y^l
    INTEGER,INTENT(IN)::px,py
    REAL(KIND(1d0)),DIMENSION(0:9,0:9),INTENT(IN)::c2
    REAL(KIND(1d0)),INTENT(IN)::xx,yy
    REAL(KIND(1d0))::p2
    INTEGER::kk,ll
    p2=0d0
    DO ll=0,py-1
       DO kk=0,px-1
          p2=p2+c2(kk,ll)*xx**kk*yy**ll
       ENDDO
    ENDDO
  END FUNCTION p2
END PROGRAM test_utils
