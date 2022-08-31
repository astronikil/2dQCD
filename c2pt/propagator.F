      subroutine propagator(ix,it,gp)
#include "size.h"
      integer ispin, ix, it, ncg, ncgs, ncg1, ncg2, l
      complex*16 gp(mv,ncs,ncs), v, w, x
      common /wrkry/ v(nc,nspr,mv), w(nc,nspr,mv), x(nc,nspr,mv)
*
      do ispin=1, ncs
         call source(v,ix,it,ispin)
*************************************************************************
*     w = (Ho)^{-2}.v                                                   *
*************************************************************************
         call invert_overlap(v,w,ncg,ncgs) 
*
         write(*,*)'ispin=',ispin, 'ncgo=',ncg, 'ncgw=',ncgs
*************************************************************************
*     x = Ho.w = (Ho)^-1.v                                              *
*************************************************************************
         call H_overlap(w, x, ncg1)
*************************************************************************
*     v = V.x = V.(Ho)^-1.v                                             *
*************************************************************************
         call V_overlap(x, v, ncg2)
*************************************************************************
*     v = 0.5(x - v) = 0.5(1.0-V).(Ho)^-1.v                             *
*************************************************************************
         v=(x-v)/2.0D0
*************************************************************************
*     v = \sigma.v = \sigma_3.[0.5(1.0-V)].[(Ho)^-1].v                  *
*************************************************************************
         do l=1, mv
            v(:,2,l)=-v(:,2,l)
         enddo
*
         call checktolinear(v,gp(:,:,ispin))
*
      enddo
*
      return
      end subroutine propagator
      subroutine source(v,ix,it,ispin)
#include "size.h"
      integer l, ino, ix, it, ispin
      complex*16 v(ncs,mv)
      common /newold/ ino(mv)
*
      v=(0.0D0,0.0D0)
*
      l=ino(ix+mx*(it-1))
      v(ispin,l)=(1.0D0,0.0D0)
*
      return
      end subroutine source 
      subroutine createsrcpts(ixex)
#include "size.h"
#include "param.h"
      integer ix, idx, ix0, ixex(nd,nex), ixsl(nd,nsl)
      idx=mx/nex
      ix0=1
      do ix=1,nex
         ixex(1,ix)=ix0 
         ixex(2,ix)=ix0 
         ix0=mod(ix0+idx-1,mx)+1
      enddo
      end
