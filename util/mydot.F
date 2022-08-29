      complex*16 function mydot(x,y)
***********************************************************************
*                                                                     *
*     The dot product: x^dagger y.                                    *
*                                                                     *
***********************************************************************
#include "size.h"
      integer l
      complex*16 dotp, x(mvs), y(mvs)
*
      dotp=(0.0D0,0.0D0)
!$OMP parallel do default(shared) reduction(+: dotp)
      do l=1,mvs
         dotp=dotp+conjg(x(l))*y(l)
      enddo
      mydot=dotp
*
      return
      end
