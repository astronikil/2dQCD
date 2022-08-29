      subroutine ranv(v)
*************************************************************************
*                                                                       *
*     Fill a vector with Gaussian random numbers                        *
*                                                                       *
*************************************************************************
#include "size.h"
      integer, parameter :: K1MAX=2*mvs
      integer j, l, iseed
      real*8 pi, tpi
      real*8 RUNUM(K1MAX)
      complex*16 zi, tpiz, v(mvs)
      common /consts/ pi, tpi, tpiz, zi
      common /random/ iseed(4)
*
      call dlarnv(1,iseed, K1MAX, RUNUM)
*
      j=1
      do l=1,mvs
         v(l)=sqrt(-log(RUNUM(j)))*exp(tpiz*RUNUM(j+1))
         j=j+2
      enddo
*
      return
      end
