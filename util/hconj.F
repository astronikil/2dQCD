      subroutine hconj(u,v,nsize)
#include "size.h"
      integer nsize, ic, jc, l
      complex*16 u(nc,nc,nsize), v(nc,nc,nsize)
      do l=1,nsize
         do ic=1,nc
            do jc=1,nc
               v(ic,jc,l)=conjg(u(jc,ic,l))
            enddo
         enddo
      enddo
      return
      end
