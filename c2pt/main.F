      program main
#include "size.h"
#include "param.h"

      integer iseed
      integer ixex(nd,nex), ixsl(nd,nsl)
      integer nitov, nitsq
      real*8 pi, tpi, rm, fm, es
      real*8 beta, fm2, time
      real*8 delitovex, delitsqex, delitov, delitsq
      complex*16 zi, tpiz, gx, gy
      complex*16 c2ptex
      common /param_overlap/ delitov, nitov
      common /param_sqrt/ delitsq, nitsq
*
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
*
      common /quark/ rm
      common /prop/ gx(mv*ncs*ncs)
*
      common /corr2ex/c2ptex(mt*n2pt,nex)
*
      iseed=0; iseed(4)=1
      pi=acos(-1.0D0)
      tpi=2.0D0*pi
      zi=(0.0D0,1.0D0)
      tpiz=tpi*zi
*      
      fm=1.0D0   !wilson mass
*
      rm = 2.0D0 - fm
*
      delitovex=1.0E-5; delitsqex=1.0E-7 
      nitov=5000
      nitsq=5000
*
      delitovex=delitovex**2
      delitsqex=delitsqex**2
*
      call addrc
      call pauliinit
*
      !call readconfig2(beta, fm2, time) 
*
      call ransun
c      call plaq(es)
c      write(*,*)'before:',es
c      call gtransform
c      call plaq(es)
c      write(*,*)'after:',es
*
      call ferminit        
*
      call createsrcpts(ixex) 
      call invert_and_contract     
      call prc2                 !print c2pt
      stop
      contains
      subroutine invert_and_contract
      integer iex, ix, it, it2, it3, isep, js
      real*8 etime, stime
      call cpu_time(stime)
      delitov=delitovex
      delitsq=delitsqex
      do iex=1,nex
         ix=ixex(1,iex)
         it=ixex(2,iex)
         write(*,*)'Source at',ix,it,'ex'
         call propagator(ix,it,gx(:))
         call correl2pt(ix,it,c2ptex(:,iex))
      enddo
      call cpu_time(etime)
      write(*,*)'exact time=',etime-stime
      end
      subroutine prc2
#include "size.h"
#include "param.h"
      integer is, ic, ip, i, iex
      complex*16 c2ptex
      common /corr2ex/ c2ptex(mt,n2pt,nex)
      write(36,*)'# nexact=',nex
      do ic=1,n2pt
      do i=1,mt
         write(36,*)i,
     1   (c2ptex(i,ic,iex),iex=1,nex),
     1   ip2str(ic)
      enddo
      enddo
      end
      end program main
