! Stubs of the LHAPDF Fortran interface (LHAGlue), linked when gammaUPC is built without LHAPDF (make, not make
! INELASTIC=1): the structure functions at high Q2 (DGLAP part of the inelastic photon) need LHAPDF.
SUBROUTINE lhapdf_initpdfset_byname(nset,name)
  INTEGER,INTENT(IN)::nset
  CHARACTER(LEN=*),INTENT(IN)::name
  WRITE(*,*)"ERROR: the inelastic photon flux needs LHAPDF: build gammaUPC with make INELASTIC=1 (set ",nset,TRIM(name),")"
  STOP
END SUBROUTINE lhapdf_initpdfset_byname
SUBROUTINE lhapdf_getdatapath(s)
  CHARACTER(LEN=*),INTENT(OUT)::s
  s=''
  STOP "ERROR: LHAPDF not linked (make INELASTIC=1)"
END SUBROUTINE lhapdf_getdatapath
SUBROUTINE lhapdf_xfxq2(nset,nmem,pid,x,q2,xf)
  INTEGER,INTENT(IN)::nset,nmem,pid
  REAL(KIND(1d0)),INTENT(IN)::x,q2
  REAL(KIND(1d0)),INTENT(OUT)::xf
  xf=0d0
  WRITE(*,*)"lhapdf_xfxq2 called with",nset,nmem,pid,x,q2
  STOP "ERROR: LHAPDF not linked (make INELASTIC=1)"
END SUBROUTINE lhapdf_xfxq2
SUBROUTINE lhapdf_alphasq2(nset,nmem,q2,as)
  INTEGER,INTENT(IN)::nset,nmem
  REAL(KIND(1d0)),INTENT(IN)::q2
  REAL(KIND(1d0)),INTENT(OUT)::as
  as=0d0
  WRITE(*,*)"lhapdf_alphasq2 called with",nset,nmem,q2
  STOP "ERROR: LHAPDF not linked (make INELASTIC=1)"
END SUBROUTINE lhapdf_alphasq2
