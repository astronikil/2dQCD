      subroutine sungenerators
#include "size.h"
      integer ic, jc, l
      complex*16 tgen, zi
      common /gellmann/ tgen(nc, nc, nca)
*
      zi=(0.0D0,1.0D0)
*
      tgen=(0.0D0,0.0D0)
*
      if(nc.eq.2)then
      
      tgen(1,1,1)=0.0D0; tgen(1,2,1)=1.0D0
      tgen(2,1,1)=1.0D0; tgen(2,2,1)=0.0D0

      tgen(1,1,2)=0.0D0; tgen(1,2,2)=zi
      tgen(2,1,2)=-zi;   tgen(2,2,2)=0.0D0

      tgen(1,1,3)=-1.0D0; tgen(1,2,3)=0.0D0
      tgen(2,1,3)=0.0D0;  tgen(2,2,3)=1.0D0

      else

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
            tgen(jc,jc,l)=sqrt(2.0D0/(
     1                       dfloat(ic)*(dfloat(ic)+1.0D0)))
         enddo
         tgen(ic+1,ic+1,l)=-dfloat(ic)*
     1           sqrt(2.0D0/(dfloat(ic)*(dfloat(ic)+1.0D0)))
      enddo
      endif
      end
c      program main
c#include "size.h"
c      integer ic, jc, l
c      complex*16 tgen
c      common /gellmann/ tgen(nc, nc, nca)
c
c      call sungenerators
c
c      do l=1,nca
c        do ic=1,nc
c           write(*,*)(tgen(ic,jc,l),jc=1,nc)
c        enddo
c        write(*,*)'###########################'
c      enddo
c      end
