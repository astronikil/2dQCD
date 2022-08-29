      subroutine writeinit
#include "size.h"
      integer ll, it, iy, ix, l, ino, ic, jc, ifl
      complex*16 r,u,fo,qf
      common /forces/ fo(nc,nc,mv,2), qf(nc,nc,mv,2)
      common /config/ u(nc,nc,mv,2)
      common/rinit/r(nc,2,mv,nf)
      common /newold/ ino(mv)
      ll=0
      do it=1,mt
      do ix=1,mx
         ll=ll+1
         l=ino(ll)
         do ic=1,nc
            do jc=1,nc
               write(80,*)ix,it,ic,jc,u(ic,jc,l,1),
     1         u(ic,jc,l,2)
               write(81,*)ix,it,ic,jc,qf(ic,jc,l,1),
     1         qf(ic,jc,l,2)
            enddo
            write(82,*)ix,it,ic,(r(ic,1,l,ifl),r(ic,2,l,ifl), ifl=1,nf)
         enddo
      enddo
      enddo
      end
