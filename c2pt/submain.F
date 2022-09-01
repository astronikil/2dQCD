      subroutine meson
#include "size.h"
#include "param.h"
      integer ixex(nd,nex), ixsl(nd,nsl)
      real*8 pi, tpi, es
      complex*16 zi, tpiz, gx, gy
      complex*16 c2ptex
*
      common /consts/ pi, tpi, tpiz, zi
*
      common /prop/ gx(mv*ncs*ncs)
*
      common /corr2ex/c2ptex(mt*n2pt,nex)
*
      call pauliinit
*
      call ferminit        
*
      call createsrcpts(ixex) 
      call invert_and_contract     
      call prc2                 !print c2pt
      contains
      subroutine invert_and_contract
      integer iex, ix, it, it2, it3, isep, js
      real*8 etime, stime
      do iex=1,nex
         ix=ixex(1,iex)
         it=ixex(2,iex)
         call propagator(ix,it,gx(:))
         call correl2pt(ix,it,c2ptex(:,iex))
      enddo
      end
      subroutine prc2
#include "size.h"
#include "param.h"
      integer is, ic, ip, i, iex
      complex*16 c2ptex
      common /corr2ex/ c2ptex(mt,n2pt,nex)
c      write(36,*)'# nexact=',nex
      do ic=1,n2pt
      do i=1,mt
         write(36,*)i,
     1   (c2ptex(i,ic,iex),iex=1,nex),
     1   ip2str(ic)
      enddo
      enddo
      end
      end
