MODULE bessel_integration
  ! Integrals of smooth integrands in kT (GeV), several at once on the same nodes:
  !   res(i) = int_0^kmax dk g_i(k) J_{nj(i)}(b*k)      (nj(i) = -1: no Bessel function)
  ! Gauss-Kronrod 15-point panels: geometric (NPERDEC per decade) from 1e-3*scale,
  ! not longer than one period 2*pi/b of the Bessel functions, with break points at
  ! the edges (kinks of the integrand, e.g. form-factor joins).
  ! A finite kmax is the upper limit. kmax = HUGE(1d0) means infinity: panels up to KRES, then (b > 0) further decades while
  ! the rest predicted from the last two decades exceeds 0.1*RTOL*int|integrand|, then
  ! the tail [K,infinity) with k = K/t, geometric panels in t down to 1e-12.
  ! (as numerical_integrate('gauss') of the Python code, modules/integration_methods.py)
  IMPLICIT NONE
  PRIVATE
  PUBLIC::bessel_integrate,bessel_integrate_finite,sort_unique
  INTEGER,PARAMETER::NPERDEC=10
  REAL(KIND(1d0)),PARAMETER::KRES=1d2,RTOL=1d-8,PANELMAX=2d5
  REAL(KIND(1d0)),PARAMETER::TWOPI=6.28318530717958647692528676656d0
  ! Kronrod 15-point nodes (positive half) and weights
  REAL(KIND(1d0)),DIMENSION(8),PARAMETER::xgk=(/0.991455371120812639206854697526329d0,&
       0.949107912342758524526189684047851d0,0.864864423359769072789712788640926d0,&
       0.741531185599394439863864773280788d0,0.586087235467691130294144845693013d0,&
       0.405845151377397166906606412076961d0,0.207784955007898467600689403773245d0,0d0/)
  REAL(KIND(1d0)),DIMENSION(8),PARAMETER::wgk=(/0.022935322010529224963732008058970d0,&
       0.063092092629978553290700663189204d0,0.104790010322250183839876322541518d0,&
       0.140653259715525918745189590510238d0,0.169004726639267902826583426598550d0,&
       0.190350578064785409913256402421014d0,0.204432940075298892414161999234649d0,&
       0.209482141084727828012999174891714d0/)
  ! Gauss 7-point weights (at the Kronrod nodes 2,4,6,8)
  REAL(KIND(1d0)),DIMENSION(4),PARAMETER::wg7=(/0.129484966168869693270611432679082d0,&
       0.279705391489276667901467771423780d0,0.381830050505118944950369775488975d0,&
       0.417959183673469387755102040816327d0/)
CONTAINS
  SUBROUTINE bessel_integrate(g,n,nj,b,scale,kmax,edges,res)
    ! g(k,n,v): the n integrands at k, without the Bessel function, returned in v(1:n)
    INTEGER,INTENT(IN)::n
    INTEGER,DIMENSION(n),INTENT(IN)::nj
    REAL(KIND(1d0)),INTENT(IN)::b,scale,kmax
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::edges
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::res
    INTERFACE
       SUBROUTINE g(k,n,v)
         REAL(KIND(1d0)),INTENT(IN)::k
         INTEGER,INTENT(IN)::n
         REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::v
       END SUBROUTINE g
    END INTERFACE
    REAL(KIND(1d0))::K,step,a,c
    REAL(KIND(1d0)),DIMENSION(n)::absres,S0,S1,r,rest
    INTEGER::i
    res=0d0; absres=0d0; S0=0d0; S1=0d0
    IF(kmax.LE.0d0)RETURN
    step=HUGE(1d0)
    IF(b.GT.0d0)step=TWOPI/b
    K=kmax ! a finite kmax is the upper limit (also above KRES); KRES only starts the tail of an infinite range
    IF(kmax.GE.HUGE(1d0))K=KRES
    a=MIN(1d-3*scale,K)
    CALL panels(0d0,a)
    CALL decades(a,K)
    IF(kmax.LT.HUGE(1d0))RETURN
    IF(b.GT.0d0)THEN
       DO WHILE(9d0*K/step.LT.PANELMAX)
          r=1d0
          WHERE(S0.GT.0d0)r=S1/S0
          rest=HUGE(1d0)
          WHERE(r.LT.1d0)rest=S1*r/(1d0-r)
          WHERE(S1.EQ.0d0)rest=0d0
          IF(ALL(rest.LE.0.1d0*RTOL*absres))EXIT
          S0=S1; S1=0d0
          K=10d0*K
          CALL decades(K/10d0,K)
       ENDDO
    ENDIF
    ! tail [K,infinity)
    a=1d-12
    DO i=1,12*NPERDEC
       c=a*10d0**(1d0/NPERDEC)
       CALL gk15(a,c,.TRUE.)
       a=c
    ENDDO
  CONTAINS
    SUBROUTINE decades(lo,hi) ! geometric panels on [lo,hi], broken at the edges (in increasing order)
      REAL(KIND(1d0)),INTENT(IN)::lo,hi
      REAL(KIND(1d0))::a,c,e
      a=lo
      DO WHILE(a.LT.hi)
         c=MIN(a*10d0**(1d0/NPERDEC),hi)
         DO
            e=MINVAL(edges,MASK=edges.GT.a.AND.edges.LT.c)
            IF(e.GE.c)EXIT
            CALL panels(a,e)
            a=e
         ENDDO
         CALL panels(a,c)
         a=c
      ENDDO
    END SUBROUTINE decades
    SUBROUTINE panels(lo,hi) ! [lo,hi] in pieces not longer than step
      REAL(KIND(1d0)),INTENT(IN)::lo,hi
      INTEGER::m,j
      REAL(KIND(1d0))::h
      IF(hi.LE.lo)RETURN
      m=MAX(1,CEILING((hi-lo)/step))
      h=(hi-lo)/m
      DO j=1,m
         CALL gk15(lo+(j-1)*h,lo+j*h,.FALSE.)
      ENDDO
    END SUBROUTINE panels
    SUBROUTINE gk15(lo,hi,tail) ! tail: [lo,hi] in t, k = K/t
      REAL(KIND(1d0)),INTENT(IN)::lo,hi
      LOGICAL,INTENT(IN)::tail
      REAL(KIND(1d0))::c0,h,kk,w,v(n),p(n),pabs(n)
      INTEGER::j,s
      c0=0.5d0*(lo+hi); h=0.5d0*(hi-lo)
      p=0d0; pabs=0d0
      DO j=1,8
         DO s=-1,1,2
            IF(j.EQ.8.AND.s.EQ.1)CYCLE
            kk=c0+s*h*xgk(j)
            w=h*wgk(j)
            IF(tail)THEN
               w=w*K/kk**2
               kk=K/kk
            ENDIF
            CALL g(kk,n,v)
            v=w*v*bessel(kk)
            p=p+v
            pabs=pabs+ABS(v)
         ENDDO
      ENDDO
      res=res+p
      absres=absres+pabs
      IF(tail)RETURN
      IF(hi.GT.K/1d2.AND.hi.LE.K/1d1)S0=S0+pabs
      IF(hi.GT.K/1d1)S1=S1+pabs
    END SUBROUTINE gk15
    FUNCTION bessel(kk) ! J_{nj(i)}(b k) for the n integrands
      REAL(KIND(1d0)),INTENT(IN)::kk
      REAL(KIND(1d0)),DIMENSION(n)::bessel
      REAL(KIND(1d0))::z,j0,j1,j2
      INTEGER::i
      z=b*kk
      j0=BESSEL_J0(z); j1=BESSEL_J1(z); j2=BESSEL_JN(2,z)
      DO i=1,n
         SELECT CASE(nj(i))
         CASE(0); bessel(i)=j0
         CASE(1); bessel(i)=j1
         CASE(2); bessel(i)=j2
         CASE DEFAULT; bessel(i)=1d0
         END SELECT
      ENDDO
    END FUNCTION bessel
  END SUBROUTINE bessel_integrate

  SUBROUTINE bessel_integrate_finite(g,n,nj,b,scale,kmax,edges,res,err,absres)
    ! res(i) = int_0^kmax dk g_i(k) J_{nj(i)}(b*k) on a finite range, panels as numerical_integrate('gauss') of the
    ! Python code (panel_edges): geometric (np.geomspace from k_lo = 1e-3*min(scale,kmax), int(NPERDEC*log10(kmax/k_lo))+2
    ! points), every Bessel period 2*pi/b from 0 (b > 0), and the edges inside (0,kmax); Gauss-Kronrod 15 per panel.
    ! err = sum |K15 - G7|, absres = int |integrand|; results below 1e-11*absres are set to zero (guard of the Python).
    INTEGER,INTENT(IN)::n
    INTEGER,DIMENSION(n),INTENT(IN)::nj
    REAL(KIND(1d0)),INTENT(IN)::b,scale,kmax
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::edges
    REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::res,err,absres
    INTERFACE
       SUBROUTINE g(k,n,v)
         REAL(KIND(1d0)),INTENT(IN)::k
         INTEGER,INTENT(IN)::n
         REAL(KIND(1d0)),DIMENSION(n),INTENT(OUT)::v
       END SUBROUTINE g
    END INTERFACE
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::e
    REAL(KIND(1d0))::klo,l0,l1,step,c0,h,kk,v(n),pk(n),pg(n),pa(n),z,j0,j1,j2,bj(n)
    INTEGER::ngeo,nar,ne,i,j,s,m
    res=0d0; err=0d0; absres=0d0
    IF(kmax.LE.0d0)RETURN
    klo=1d-3*MIN(scale,kmax)
    ngeo=INT(NPERDEC*LOG10(kmax/klo))+2
    nar=0
    IF(b.GT.0d0)THEN
       step=TWOPI/b
       IF(kmax.GT.step)nar=CEILING(kmax/step)
    ENDIF
    ALLOCATE(e(ngeo+nar+2+SIZE(edges)))
    l0=LOG10(klo); l1=LOG10(kmax)
    DO i=0,ngeo-1
       e(i+1)=10d0**(l0+i*((l1-l0)/(ngeo-1)))
    ENDDO
    e(1)=klo; e(ngeo)=kmax
    ne=ngeo
    e(ne+1)=0d0; e(ne+2)=kmax; ne=ne+2
    DO i=0,nar-1
       ne=ne+1; e(ne)=i*step
    ENDDO
    DO i=1,SIZE(edges)
       IF(edges(i).GT.0d0.AND.edges(i).LT.kmax)THEN
          ne=ne+1; e(ne)=edges(i)
       ENDIF
    ENDDO
    CALL sort_unique(e,ne,m)
    DO i=1,m-1
       c0=0.5d0*(e(i)+e(i+1)); h=0.5d0*(e(i+1)-e(i))
       pk=0d0; pg=0d0; pa=0d0
       DO j=1,8
          DO s=-1,1,2
             IF(j.EQ.8.AND.s.EQ.1)CYCLE
             kk=c0+s*h*xgk(j)
             CALL g(kk,n,v)
             z=b*kk
             j0=BESSEL_J0(z); j1=BESSEL_J1(z); j2=BESSEL_JN(2,z)
             WHERE(nj.EQ.0)
                bj=j0
             ELSEWHERE(nj.EQ.1)
                bj=j1
             ELSEWHERE(nj.EQ.2)
                bj=j2
             ELSEWHERE
                bj=1d0
             END WHERE
             v=h*v*bj
             pk=pk+wgk(j)*v
             pa=pa+wgk(j)*ABS(v)
             IF(MOD(j,2).EQ.0)pg=pg+wg7(j/2)*v
          ENDDO
       ENDDO
       res=res+pk; err=err+ABS(pk-pg); absres=absres+pa
    ENDDO
    WHERE(ABS(res).LE.1d-11*absres)res=0d0
    DEALLOCATE(e)
  END SUBROUTINE bessel_integrate_finite

  SUBROUTINE sort_unique(a,n,m) ! a(1:n) sorted increasingly (heapsort), exact duplicates removed: a(1:m)
    REAL(KIND(1d0)),DIMENSION(:),INTENT(INOUT)::a
    INTEGER,INTENT(IN)::n
    INTEGER,INTENT(OUT)::m
    INTEGER::i,j,k,l,ir
    REAL(KIND(1d0))::t
    IF(n.GT.1)THEN
       l=n/2+1; ir=n
       DO
          IF(l.GT.1)THEN
             l=l-1; t=a(l)
          ELSE
             t=a(ir); a(ir)=a(1); ir=ir-1
             IF(ir.EQ.1)THEN
                a(1)=t; EXIT
             ENDIF
          ENDIF
          i=l; j=l+l
          DO WHILE(j.LE.ir)
             IF(j.LT.ir)THEN
                IF(a(j).LT.a(j+1))j=j+1
             ENDIF
             IF(t.LT.a(j))THEN
                a(i)=a(j); i=j; j=j+j
             ELSE
                j=ir+1
             ENDIF
          ENDDO
          a(i)=t
       ENDDO
    ENDIF
    m=MIN(n,1)
    DO k=2,n
       IF(a(k).NE.a(m))THEN
          m=m+1; a(m)=a(k)
       ENDIF
    ENDDO
  END SUBROUTINE sort_unique
END MODULE bessel_integration
