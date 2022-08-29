*************************************************************************
*                                                                       *
*     All lattice extents must be even.                                 *
*                                                                       *
*************************************************************************
      include "mpif.h"
c      integer, parameter :: mx=28, my=28, mt=48, nd=3
      integer, parameter :: mx=40, my=20, mt=48, nd=3
      integer, parameter :: nc=2, nc2=nc*nc, ncs=2*nc
      integer, parameter :: mv=mt*mx*my, mb=nd*mv, mv2=(nd-1)*mv
      integer, parameter :: mvd2=mv/2, msp=mx*my
      integer, parameter :: nspr=2,mvs=nspr*mv,mbs=nspr*mb
      integer, parameter :: ntot=nc2*mb
      integer, parameter :: mux=0, muy=mv, mut=2*mv
      integer, parameter :: nvsize=2*nc*mv

      integer, parameter :: nprocs=4000, mprocs=1
      integer myid,numprocs

      real, parameter :: rzeroep=1.0E-10

      
      parameter (nex=1, nsl=0)
      parameter (n2pt=4,n3pt=4)

      parameter (nf=1)
      parameter (nfhe=nf/2, nfho=(nf-1)/2)

      parameter(nmeson=2)

      character*2 ip2str(n2pt)
      character*3 ixij(3)
      data ip2str /'V1','V2','V3','S0'/
      data ixij /'X12','X23','X31'/
