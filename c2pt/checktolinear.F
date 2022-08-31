      subroutine checktolinear(x,y)
#include "size.h"
      integer ics, ino, it, ix, j, jc
      complex*16 x(ncs,mv), y(mv,ncs)
      common /newold/ ino(mv)
*
      j=0
      do it=1,mt
         do ix=1,mx
            j=j+1
            jc=ino(j)
            do ics=1,ncs
               y(j,ics)=x(ics,jc)
            enddo
         enddo
      enddo
*
      return
      end subroutine checktolinear
