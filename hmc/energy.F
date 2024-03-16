      subroutine energy(beta,em,es,ef,ene,ncgo,ncgw)
************************************************************************
*                                                                      *
*     Calculates total energy E=P^2/2+S(plaq)+S(fermion)               *
*                                                                      *
************************************************************************
#include "size.h"
*
      real*8, parameter :: DACCU1=1.0E-7, DACCU2=1.0E-9
c      real*8, parameter :: DACCU1=1.0E-5, DACCU2=1.0E-7
*
      integer ifl, l, ncgo, ncgw, nitov, nitsq, ic, jc
      real*8 beta, delitov, delitsq, dovhold, dovsq 
      real*8 es, ef, em, ene
      complex*16 rphi, mydot
      complex*16 fo, qf
      complex*16 xsol(mvs)
c      complex*16 x1(mvs),x2(mvs),x3(mvs),x4(mvs) !!!
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
         if(fermtype.eq.'D')then
            call invert_overlap(rphi(:,ifl),xsol,ncgo,ncgw)
         elseif(fermtype.eq.'M')then
            write(*,*)'did M in energy'
            call H_overlap_pow_n1by2(rphi(:,ifl),xsol,ncgo,ncgw)
         endif
         ef=ef+real(mydot(rphi(:,ifl),xsol))
      enddo

c      x1=rphi(:,1)
c      call H_overlap_pow_n1by2(x1,x2,ncgo,ncgw) !x2=H^-0.5 x1
c      call H_overlap_pow_p1by4(x2,x3,ncgo,ncgw) !x2=H^-0.5 x1
c      call H_overlap_pow_n1by2(x3,x2,ncgo,ncgw) !x2=H^-0.5 x1
c      !call invert_overlap(x1,x3,ncgo,ncgw)
c      !call H_overlap(x3,x4,ncgo)
c      !call H_overlap(x4,x2,ncgo)
c      write(*,*)'norms:',mydot(x1,x1),mydot(x2,x2)
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
*
      return
      end subroutine energy
