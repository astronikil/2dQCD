      subroutine ranalgebra
*************************************************************************
*                                                                       *
*     Generate random  SU(N) matrices                                   *
*                                                                       *
*************************************************************************
#include "size.h"
      integer, parameter :: K2MAX=nca*mb 
      integer j, k, l, iseed, ino, ica
      real*8 pi, tpi, rmag
      real*8 RUNUM(K2MAX)
      complex*16 u, tang, tgen
      common /config/ u(nc2,mb)
      common /expdat/ tang(nc2, mb)
      common /gellmann/ tgen(nc2, nca)
      common /random/ iseed(4)
*
      call sungenerators
      call dlarnv(2,iseed,K2MAX,RUNUM)
*
      tang=(0.0D0,0.0D0)
      j=1
      do l=1,mb
         do ica=1,nca
            tang(:,l)=tang(:,l)+RUNUM(j)*tgen(:,ica)
            j=j+1
         enddo
      enddo
*
      call expsun
*
      u = tang
      write(28,'(2x,a)') 'Initialised to random link configuration from
     1  algbra'
      return
      end
