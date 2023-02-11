      subroutine z2noise(r)
************************************************************************
*                                                                      *
*     Initialisation of z2 noise vector                                *
*                                                                      *
************************************************************************
#include "size.h"
*
      integer iseed, j, l
      complex*16 r(mvs)
      real*8 fnorm, xre, xim
      real*8 RUMNUM(2*mvs)
      common /random/ iseed(4)
      fnorm=1.0/sqrt(2.0)

*
      call dlarnv(1,iseed,2*mvs, RUMNUM)
*      
      j = 1
      do l=1,mvs
         xre = RUMNUM(j)
         xim = RUMNUM(j+1)
         if(xre.le.0.5)then
            xre = 1.0
         else
            xre = -1.0
         endif

         if(xim.le.0.5)then
            xim = 1.0
         else
            xim = -1.0
         endif

         r(l) = cmplx(xre,xim)*fnorm

         j = j+2
      enddo
      return
      end
