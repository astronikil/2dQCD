*************************************************************************
*                                                                       *
*     All lattice extents must be even.                                 *
*                                                                       *
*************************************************************************
      implicit none
c      include "mpif.h"
      integer, parameter :: mx=4, mt=4, nd=2

c     nc=2 will remain since we want work with SU(2)
c     ncr is the size of the representation of the fermions and n=2j+1 where j is the representation.

      integer, parameter :: nc=2, nc2=nc*nc
      integer, parameter :: twoXj=2, ncr=twoXj+1
c      
      integer, parameter :: nca=nc2-1
      integer, parameter :: mv=mx*mt, mb=nd*mv, mv2=(nd-1)*mv
      integer, parameter :: mvd2=mv/2, msp=mx
      integer, parameter :: nspr=2
      integer, parameter :: ntot=nc2*mb
      integer, parameter :: mux=0, mut=mv
      integer, parameter :: ncs=nspr*ncr, mvs=ncs*mv, mvs2=2*mvs

      integer, parameter :: nf=1
      integer, parameter :: nfhe=nf/2, nfho=(nf-1)/2

      integer, parameter :: nprocs=4000, mprocs=1
      integer myid,numprocs
