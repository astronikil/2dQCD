      !subroutine energy(beta,em,es,ef,ene)
      subroutine energy(beta,em,es,ef,ene,ncgo,ncgw)
************************************************************************
*                                                                      *
*     Calculates total energy E=P^2/2+S(plaq)+S(fermion)               *
*                                                                      *
************************************************************************
#include "size.h"
*
      real*8, parameter :: DACCU1=1.0E-12, DACCU2=1.0E-13
c      real*8, parameter :: DACCU1=1.0E-5, DACCU2=1.0E-7
*
      integer ifl, l, ncgo, ncgw, nitov, nitsq, ic, jc
      real*8 beta, delitov, delitsq, dovhold, dovsq 
      real*8 es, ef, em, ene
      complex*16 rphi, mydot
      complex*16 fo, qf
      complex*16 xsol(mvs)
*
      common /phi/ rphi(mvs,nf)
      common /forces/ fo(ntot), qf(ntot)
*
      common /param_overlap/ delitov, nitov
      common /param_sqrt/ delitsq, nitsq
*
      call plaq(es)
      es=2.0D0*beta*es
*
      dovhold=delitov
      dovsq=delitsq
      delitov=(DACCU1)**2
      delitsq=(DACCU2)**2
*
      ef=0.0D0
      do ifl=1,nf
         call H_overlap_pow_n1by2(rphi(:,ifl),xsol,ncgo,ncgw)
         ef=ef+real(mydot(rphi(:,ifl),xsol))
      enddo
*
      delitov=dovhold
      delitsq=dovsq
*
      em=0.0D0
      do l=1,ntot
         em=em+0.5D0*abs(qf(l))**2
      enddo
*
      ene=em+es+ef
c      ene=em+es
*
      return
      end subroutine energy
