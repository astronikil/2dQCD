      subroutine H_overlap(v,w,ncg)
#include "size.h"
      integer ncg, l
      real*8 r1, r2, fmo
      complex*16 v(ncr,nspr,mv), w(ncr,nspr,mv)
      complex*16 x, y
      complex*16 mydot
      common /quark_ov/ fmo
      common /fermwrk1/ x(ncr,nspr,mv), y(ncr,nspr,mv)
*
      r1=0.5D0*(1.0D0+fmo)
      r2=0.5D0*(1.0D0-fmo)
      call H_wilson(v,x)
      call sqrtoverlap(x,y,ncg)

c      write(*,*)'##### eps=',mydot(y,y)/mydot(v,v)
*
      do l=1,mv
         w(:,1,l) =  r1*v(:,1,l) + r2*y(:,1,l)
         w(:,2,l) = -1.0D0*r1*v(:,2,l) + r2*y(:,2,l)
      enddo
*
      return
      end
      subroutine V_overlap(v,w,ncg)
#include "size.h"
      integer ncg, l
      real*8 r1, r2
      complex*16 v(ncr,nspr,mv), w(ncr,nspr,mv)
      complex*16 x, y
      common /fermwrk1/ x(ncr,nspr,mv), y(ncr,nspr,mv)
*
      call H_wilson(v,x)
      call sqrtoverlap(x,w,ncg)
*
      do l=1,mv
         w(:,2,l) = -1.0D0*w(:,2,l)
      enddo
*
      return
      end
      subroutine Vdag_overlap(v,w,ncg)
#include "size.h"
      integer ncg, l
      real*8 r1, r2
      complex*16 v(ncr,nspr,mv), w(ncr,nspr,mv)
      complex*16 x, y
      common /fermwrk1/ x(ncr,nspr,mv), y(ncr,nspr,mv)
*      
      do l=1,mv
         w(:,1,l) =  v(:,1,l)
         w(:,2,l) = -1.0D0*v(:,2,l)
      enddo
*
      call H_wilson(w,x)
      call sqrtoverlap(x,w,ncg)
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
c      subroutine VVdag_overlap(v,w,ncg)
c***** This subroutine not tested. This is a stub ********  
c#include "size.h"
c      integer ncg, l
c      real*8 r1, r2
c      complex*16 v(nc,nspr,mv), w(nc,nspr,mv)
c      complex*16 x, y
c      common /fermwrk2/ x(nc,nspr,mv), y(nc,nspr,mv)
c*      
c      call V_overlap(v, x)
c      call Vdag_overlap(v, y)
c*
c      w = 2.0D0*v + x + y
c*
c      return
c      end
