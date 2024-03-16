      subroutine readinit
#include "size.h"
      integer ll, it, iy, ix, l, ino, ic, jc
      integer ix1,it1,ic1,jc1,mu,mu1
      complex*16 u, zz
      common /config/ u(nc,nc,mv,2)
      common /newold/ ino(mv)
      do mu=1,2
      do ix=1,mx
         do it=1,mt
            do ic=1,nc
               do jc=1,nc
                  read(80,*)ix1,it1,ic1,jc1,mu1,zz 
                  l=ino(ix1+mx*(it1-1))
                  u(ic1,jc1,l,mu1)=zz
               enddo
            enddo
          enddo
      enddo
      enddo

c      call unitar
      end
