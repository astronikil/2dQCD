      subroutine H_overlap(v,w,ncg)
#include "size.h"
      integer ncg, l
      real*8 r1, r2
      complex*16 v(nc,nspr,mv), w(nc,nspr,mv)
      complex*16 x, y
      common /fermwrk1/ x(nc,nspr,mv), y(nc,nspr,mv)
*
      r1=1.0D0/2.0D0
      r2=1.0D0/2.0D0
      call H_wilson(v,x)
      call sqrtoverlap(x,y,ncg)
*
      do l=1,mv
         w(:,1,l) =  r1*v(:,1,l) + r2*y(:,1,l)
         w(:,2,l) = -r1*v(:,2,l) + r2*y(:,2,l)
      enddo
*
      return
      end
      subroutine sqrtoverlap(v,w,ncg)
#include "size.h"
#include "zolo.h"
      integer ncg, izol
      complex*16 v(mvs), w(mvs), x, ps
      common /wrkforce/ ps(mvs,nzol), x(mvs,nzol)
*
      call invert_wilson_shifted(v,ncg)
*      
      w=rz(0)*v
      do izol=1,nzol
         w=w+rz(izol)*x(:,izol)
      enddo
      end
