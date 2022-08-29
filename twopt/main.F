      program main
**************************************************************
*     Code to compute ITD in 3d SU(2) gauge theory           *
*     Author: Nikhil Karthik Oct 2020                        *
**************************************************************
#include "size.h"
#include "mom.h"
      complex zi, tpizi, gx, gy, gx_back, pauli, sigH(ncs,ncs)
      complex gx2pt, gx2pt_back
      integer ixex(3,nex), ixsl(3,nsl), px, py, qx, qy
      character*2 hstr
      common /param_sqrt/ delitsq, nitsq
*
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
      common /paulimat/ pauli(ncs,ncs,4)
*
      integer v1,v2,v3,iproc,ierr
      character tt1,tt2,tt3,fix(0:nprocs-1)*3
*
      common /quark/ rm
*      
      common /apen/ nlevape
      common /apee/ epsape
      common /boostn/ nfsmear
      common /boostd/ delta
      common /boostmom/ rqm(2)
**************************************************************
*     Source point used                                      *
*     All points are referenced wrt to this                  *
**************************************************************
      common /origin/ iox, iot
**************************************************************
*     [propagator]*[smearing_kernel]                         *
**************************************************************
      common /prop_forward/ gx(mv*ncs*ncs)
      common /prop_backward/ gx_back(mv*ncs*ncs)
**************************************************************
*     [smearing_kernel]*[propagator]*[smearing_kernel]       *
**************************************************************
      common /prop_forward_2pt/ gx2pt(mv*ncs*ncs)
      common /prop_backward_2pt/ gx2pt_back(mv*ncs*ncs)
**************************************************************
*     Sequential propagator                                  *
*     [propagator](x,y')*sigma_H*[smearing](y',z)*exp(p*z)*  *
*     [propagator](z,w)                                      *
**************************************************************
      common /prop_seq/ gy(mv*ncs*ncs)
**************************************************************
*     Correlator outputs                                     *
**************************************************************
*      
*
      myid=0; numprocs=1
c      call MPI_INIT(ierr)
c      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
*
      if(myid.eq.0)then
         do iproc=0,numprocs-1
            v1=iproc/100
            v2=iproc/10
            v3=v2/10
            v2=v2-10*v3
            v3=iproc-v2*10-v1*100
            tt1=char(48+v1)
            tt2=char(48+v2)
            tt3=char(48+v3)
            fix(iproc)=tt1//tt2//tt3
         enddo
      endif
*
c      call MPI_BCAST(fix(0),nprocs*3,MPI_CHARACTER,0,
c     &     MPI_COMM_WORLD,ierr)
*
      open(10,file='gauge'//fix(myid),status='unknown')
      open(27,file='inp'//fix(myid),status='unknown')
      open(14,file='output'//fix(myid),status='unknown')
      open(36,file='c2pt'//fix(myid),status='unknown')
      open(37,file='fpi'//fix(myid),status='unknown')
      open(101,file='c3pt'//fix(myid),status='unknown')
*
      !take care of initializations
**************************************************************
      read(27,*)fm                        ! Quark mass
      read(27,*)hstr                      ! hstr="S0" or "V2" (scalar or
                                         !                    vector)
      read(27,*)ibool3pt                  ! 1 (calculate 3pt), or don't

      read(27,*)nlevape, epsape
      read(27,*)nfsmear, delta
      read(27,*)rkmom

      read(27,*)iox, ioy, iot
**************************************************************
*
      iseed=0; iseed(4)=1
      pi=acos(-1.)
      tpi=2*pi
      zi=(0.,1.)
      tpiz=tpi*zi
*
      rm =-3.0 + fm
      delitsq=1.0E-8 !exact cg stopping criteria
      nitsq=5000
      delitsq=delitsq**2
*
      rqm(1)=2.0*acos(-1.0)*rkmom/mx
      rqm(2)=0.0
******************************************************************
*     Init pauli matrices= (2x2 sigma)(dir prod)(ncxnc identitity*
******************************************************************
      call pauliinit
      if(hstr.eq."S0")then
         sigH=pauli(:,:,4)
      elseif(hstr.eq."V2")then
         sigH=pauli(:,:,2)
      elseif(hstr.eq."V1")then
         sigH=pauli(:,:,1)
      elseif(hstr.eq."V3")then
         sigH=pauli(:,:,3)
      endif
**************************************************************
*     Lattice addressing init                                *
**************************************************************
      call addrc
**************************************************************
*     Input gauge fields and print out plaq                  *
**************************************************************
c      call ordsun
      call readconfig(beta, fm, dt, nitrc, delit, ntraj,
     1           nlen, istart, time, nmd, nstep)
c      call plaqq(es); write(*,*)es
c      call gtransform

      call plaqq(es); write(14,*)'plaquette: ',es
c      call flush

**************************************************************
*     Initialization on gauge links before inversions        *
**************************************************************
      call ferminit        !initialize wilson-dirac; smear,csw 
      call initu           !take care of link initializations, wilson
                           !   ..etc
      call init_boost_smear !do 2d ape etc
**************************************************************
*      Create source point. All point are in its reference.  *
**************************************************************
c      call createsrcpts
**************************************************************
*      Do the necessary inversions                           *
*      compute 2pt and 3pt contractions                      *
**************************************************************
      call invert_and_contract(ibool3pt,sigH)     !invert wilson and contract
**************************************************************
*      Write out 2 pt                                        *
**************************************************************
      call prc2                 !print c2pt
c      call prcfpi               !print c2pt
**************************************************************
*      Write out 3 pt                                        *
**************************************************************
      call prc3                 !print c3pt
*
      close(10)
      close(27)
      close(14)
      close(36)
      close(37)
      close(101)
c      call MPI_FINALIZE(ierr)
      end program main
