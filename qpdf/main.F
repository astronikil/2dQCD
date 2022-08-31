      program main
#include "size.h"
      complex zi, tpizi, gx, gy, xt
      integer ixex(3,nex), ixsl(3,nsl), jtau
      common /param_overlap/ delitov, nitov
      common /param_sqrt/ delitsq, nitsq

      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi

      common /quark/ rm
      common /prop/ gx(mv*ncs*ncs),gy(mv*ncs*ncs)

      common /corr3ex/c3ptex(mt*n3pt,nsep,nex)
      common /corr2ex/c2ptex(mt*3*n2pt,nsep,nex)
      common /taus/jtau(nsep)

      !take care of initializations
      !define fmfac1, fmfac2 and connect to commons

      read(*,*)jtau(1:nsep)

      iseed=0; iseed(4)=1
      pi=acos(-1.0D0)
      tpi=2.0D0*pi
      zi=(0.0D0,1.0D0)
      tpiz=tpi*zi
      xt=0.0D0
      fm=1.0D0
      delitovex=1.0E-5; delitsqex=1.0E-7
      nitov=5000
      nitsq=5000

      delitovex=delitovex**2
      delitsqex=delitsqex**2
*
      rm =2.0 - fm
*
      call addrc
*
      call ordu1
c      call readrajamani
c      call gtransform
      call ferminit
*
      call createsrcpts(ixex,ixsl)
*
*        Exact correlators
*
      call cpu_time(stime)
      delitov=delitovex
      delitsq=delitsqex
      do iex=1,nex
         ix=ixex(1,iex)
         it=ixex(2,iex)
         write(*,*)ix,it,'ex'
         call propagator(ix,it,gx(:))
         do js=1,nsep
            isep=jtau(js)
            it3=mod(it+isep-1,mt)+1
            it2=mod(it+isep/2-1,mt)+1
            call propagator(ix,it3,gy(:))
            call correl3pt(ix,iy,it,it2,it3,c3ptex(:,js,iex))
            call correl2pt(ix,iy,it,it2,it3,c2ptex(:,js,iex))
         enddo
      enddo
      call cpu_time(etime)
      write(*,*)'exact time=',etime-stime
*
      call prc2
      call prc3
      stop
      end program main
      subroutine prc2(q)
#include "size.h"
      common /taus/jtau(nsep)
      common /corr2ex/c2ptex(mt,3,n2pt,nsep,nex)
      common /corr2sl/c2ptsl(mt,3,n2pt,nsep,nsl)
      write(36,*)'# nexact=',nex, 'sloppy=',nsl
      do is=1,nsep
      do ic=1,n2pt
      do ip=1,3
      do i=1,mt
         write(36,*)i,
     1 (c2ptex(i,ip,ic,is,iex),iex=1,nex),
     1 (c2ptsl(i,ip,ic,is,isl),isl=1,nsl),
     1 ixij(ip),' ',
     1 ip2str(ic),' ',
     1 'tau=',jtau(is),'Q=',q
      enddo
      enddo
      enddo
      enddo
      end
      subroutine prc3(q)
#include "size.h"
      common /taus/jtau(nsep)
      common /corr3ex/c3ptex(mt,n3pt,nsep,nex)
      common /corr3sl/c3ptsl(mt,n3pt,nsep,nsl)
      write(37,*)'# nexact=',nex, 'sloppy=',nsl
      do is=1,nsep
      do ic=1,n3pt
      do i=1,mt
         write(37,*)i,
     1 (c3ptex(i,ic,is,iex),iex=1,nex),
     1 (c3ptsl(i,ic,is,isl),isl=1,nsl),
     1 ip3str(ic),' ',
     1 'tau=',jtau(is),'Q=',q
      enddo
      enddo
      enddo
      end
      subroutine createsrcpts(ixex,ixsl)
#include "size.h"
      integer ixex(3,nex), ixsl(3,nsl)
      idx=mx/nex
      ix0=1
      do ix=1,nex
         ixex(1,ix)=ix0 
         ixex(2,ix)=ix0 
         ixex(3,ix)=ix0 
         ix0=mod(ix0+idx-1,mx)+1
      enddo
*
      idx=mx/nsl
      ix0=1
      do ix=1,nsl
         ixsl(1,ix)=ix0 
         ixsl(2,ix)=ix0 
         ixsl(3,ix)=ix0 
         ix0=mod(ix0+idx-1,mx)+1
      enddo
      end
      subroutine propagator(ix,iy,it,gp)
#include "size.h"
      complex gp(mv,2,2), v, w, x
      common /wrkry/ v(2,mv), w(2,mv), x(2,mv)
*
      do ispin=1, 2
         call source(v,ix,iy,it,ispin)
*
         call invert_overlap(v,w,ncg,ncgs)
*
c         write(55,*)ispin, ncg, ncgs
c         call flush
*
         call vmat(w,v,ncg1)
         call vdagmat(w,x,ncg2)
*
c         write(55,*)ispin, ncg1, ncg2
c         call flush
*
         v=v-x
*
         call checktolinear(v,gp(:,:,ispin))
*
      enddo
*
      return
      end subroutine propagator
      subroutine source(v,ix,iy,it,ispin)
#include "size.h"
      complex v(2,mv)
      common /newold/ ino(mv)
*
      v=(0.,0.)
*
      l=ino(ix+mx*(iy-1+my*(it-1)))
      v(ispin,l)=(1.,0.)
*
      return
      end subroutine source 
      subroutine checktolinear(x,y)
#include "size.h"
      complex x(2,mv), y(mv,2)
      common /newold/ ino(mv)
*
      j=0
      do it=1,mt
         do iy=1,my
            do ix=1,mx
               j=j+1
               jc=ino(j)
               y(j,1)=x(1,jc)
               y(j,2)=x(2,jc)
            enddo
         enddo
      enddo
*
      return
      end subroutine checktolinear
