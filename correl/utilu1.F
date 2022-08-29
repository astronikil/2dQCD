      complex function mydot(x,y)
***********************************************************************
*                                                                     *
*     The dot product: x^dagger y.                                    *
*                                                                     *
***********************************************************************
#include "size.h"
      parameter (ns2=2*mv)
      complex dotp, x(ns2), y(ns2)
*
      dotp=(0.,0.)
!$OMP parallel do default(shared) reduction(+: dotp)
      do l=1,ns2
         dotp=dotp+conjg(x(l))*y(l)
      enddo
      mydot=dotp
*
      return
      end
