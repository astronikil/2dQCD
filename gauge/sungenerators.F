      subroutine sungenerators
#include "size.h"
      integer ic, jc, l
      complex*16 tgen, zi
      common /gellmann/ tgen(nc, nc, nca)
*
      zi=(0.0D0,1.0D0)
*
      tgen=0.0D0
*
      l=0
      do ic=1,nc
         do jc=ic+1,nc
            l=l+1
            tgen(ic,jc,l)=1.0D0
            tgen(jc,ic,l)=1.0D0
         enddo
      enddo
*
      do ic=1,nc
         do jc=ic+1,nc
            l=l+1
            tgen(ic,jc,l)=-zi
            tgen(jc,ic,l)=zi
         enddo
      enddo
*
      do ic=1,nc-1
         l=l+1
         do jc=1,ic
            tgen(jc,jc,l)=sqrt(2.0D0/(ic*(ic+1.0D0)))
         enddo
         tgen(ic+1,ic+1,l)=-ic*sqrt(2.0D0/(ic*(ic+1.0D0)))
      enddo
      end
*      
c      program main
cinclude "size.h"
c     integer ic, jc, l
c     complex*16 tgen
c     common /gellmann/ tgen(nc, nc, nca)

c     call sungenerators

c     do l=1,nca
c        do ic=1,nc
c           write(*,*)(tgen(ic,jc,l),jc=1,nc)
c        enddo
c        write(*,*)'###########################'
c     enddo
c     end
