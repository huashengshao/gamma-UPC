! Quick development check of the kernel-factorised engine (FFE):
!   test_ffe <A1> <Z1> <A2> <Z2> <E1> <E2> <chff 0/1> [nodescale] [tag1 tag2]
! prints direct values, PNOHAD=1 values and the grid timing.
PROGRAM test_ffe
  USE OpticalGlauber_Geometry
  USE ElasticPhotonPhotonFlux
  IMPLICIT NONE
  include 'run90.inc'
  CHARACTER(len=64)::arg
  INTEGER::A1,Z1,A2,Z2,ichff,i,nodes,nv1,nv2,nargs
  REAL(KIND(1d0))::E1,E2,z,zp,zg,t0,t1,sc
  REAL(KIND(1d0)),DIMENSION(5)::times
  REAL(KIND(1d0)),DIMENSION(8)::xs=(/1d-6,1d-4,1d-3,3.16227766016838d-3,1d-2,2d-2,5d-2,1d-1/)
  nargs=COMMAND_ARGUMENT_COUNT()
  CALL GET_COMMAND_ARGUMENT(1,arg); READ(arg,*)A1
  CALL GET_COMMAND_ARGUMENT(2,arg); READ(arg,*)Z1
  CALL GET_COMMAND_ARGUMENT(3,arg); READ(arg,*)A2
  CALL GET_COMMAND_ARGUMENT(4,arg); READ(arg,*)Z2
  CALL GET_COMMAND_ARGUMENT(5,arg); READ(arg,*)E1
  CALL GET_COMMAND_ARGUMENT(6,arg); READ(arg,*)E2
  CALL GET_COMMAND_ARGUMENT(7,arg); READ(arg,*)ichff
  sc=1d0
  IF(nargs.GE.8)THEN
     CALL GET_COMMAND_ARGUMENT(8,arg); READ(arg,*)sc
  ENDIF
  IF(nargs.GE.10)THEN
     CALL GET_COMMAND_ARGUMENT(9,arg); READ(arg,*)neutron_tagging(1)
     CALL GET_COMMAND_ARGUMENT(10,arg); READ(arg,*)neutron_tagging(2)
  ENDIF
  FFE_NodeScale=sc
  FFE_Verbose=2
  ebeam(1)=E1
  ebeam(2)=E2
  nuclearA_beam1=A1
  nuclearZ_beam1=Z1
  nuclearA_beam2=A2
  nuclearZ_beam2=Z2
  USE_CHARGEFORMFACTOR4PHOTON=ichff.EQ.1
  CALL CPU_TIME(t0)
  z=FFE_Flux_AB(1d-3,1d-3,1)
  CALL CPU_TIME(t1)
  WRITE(*,*)"first direct call (setup + K) [s]:",t1-t0
  DO i=1,8
     CALL CPU_TIME(t0)
     z=FFE_Flux_AB(xs(i),xs(i),1)
     CALL CPU_TIME(t1)
     zp=FFE_Flux_AB(xs(i),xs(i),3)
     WRITE(*,'(A,ES12.4,A,ES24.16,A,ES24.16,A,F8.5,A,F8.4)')" x=",xs(i)," direct=",z," pnohad1=",zp,&
          " S=",z/zp," t[ms]=",(t1-t0)*1d3
  ENDDO
  CALL CPU_TIME(t0)
  zg=FFE_Flux_AB(1.1d-3,2.3d-3,0)
  CALL CPU_TIME(t1)
  WRITE(*,*)"grid build [s]:",t1-t0
  CALL FFE_Info(times,nodes,nv1,nv2)
  WRITE(*,*)"times setup,P,W,K,grid:",times
  WRITE(*,*)"nodes, nv:",nodes,nv1,nv2
  DO i=1,8
     zg=FFE_Flux_AB(xs(i),xs(i),0)
     z=FFE_Flux_AB(xs(i),xs(i),1)
     WRITE(*,'(A,ES12.4,A,ES24.16,A,ES10.2)')" x=",xs(i)," grid=",zg," grid/direct-1=",zg/z-1d0
  ENDDO
END PROGRAM test_ffe
