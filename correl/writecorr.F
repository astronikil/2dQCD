      subroutine writecorrels
#include "size.h"
      parameter (nsize3=mt*n3pt*nsep, nsize2=mt*3*n2pt*nsep)
      common /corr3ex/c3ptex(nsize3,nex)
      common /corr3sl/c3ptsl(nsize3,nex)
      common /corr2ex/c2ptex(nsize2,nsl)
      common /corr2sl/c2ptsl(nsize2,nsl)
      do i=1,nex
         do l=1,nsize3
            write(40,*)c3ptex(l,i)
         enddo
      enddo
      do i=1,nsl
         do l=1,nsize3
            write(40,*)c3ptsl(l,i)
         enddo
      enddo
      do i=1,nex
         do l=1,nsize2
            write(41,*)c2ptex(l,i)
         enddo
      enddo
      do i=1,nsl
         do l=1,nsize2
            write(41,*)c2ptsl(l,i)
         enddo
      enddo
      end
