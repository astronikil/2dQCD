      subroutine gtransform
#include "size.h"
      integer iseed(4), iup, idn, l, ll
      complex*16 u, g, gc, y(nc,nc)
      common /iupidn/ iup(mb), idn(mb)
      common /config/ u(nc,nc,mb)
      common /gomega/ g(nc,nc,mv), gc(nc,nc,mv)
*      
      call rang(g)
      call hconj(g,gc,mv)
*
      do l=1,mb
         ll=mod(l-1,mv)+1
         y(:,:)=matmul(g(:,:,ll),u(:,:,l))
         u(:,:,l)=matmul(y(:,:),gc(:,:,iup(l)))
      enddo
*      
      end
      subroutine rang(g)
*************************************************************************
*                                                                       *
*     Generate random  su(2) matrices                                   *
*                                                                       *
*************************************************************************
#include "size.h"
      integer l
      complex*16 g(nc,nc,mv)
      complex*16 u, ud
      common /config/ u(nc,nc,mb)
      common /donfig/ ud(nc,nc,mb)
*
      ud = u
      call ransun
*
      do l=1, mv
         g(:,:,l)=u(:,:,l)
      enddo
*
      u = ud
      return
      end
