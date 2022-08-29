      complex function mydot(x,y)
***********************************************************************
*                                                                     *
*     The dot product: x^dagger y.                                    *
*                                                                     *
***********************************************************************
#include "size.h"
      complex dotp, x(nvsize), y(nvsize)
*
      dotp=(0.,0.)
!$OMP parallel do default(shared) reduction(+: dotp)
      do l=1,nvsize
         dotp=dotp+conjg(x(l))*y(l)
      enddo
      mydot=dotp
*
      return
      end
