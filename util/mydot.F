      complex*16 function mydot(x,y)
***********************************************************************
*                                                                     *
*     The dot product: x^dagger y.                                    *
*                                                                     *
***********************************************************************
#include "size.h"
      integer l
      complex*16 dotp, x(mvs), y(mvs)

      mydot=dot_product(x,y)
*
c      dotp=(0.0D0,0.0D0)
c!$OMP parallel do default(shared) reduction(+: dotp)
c      do l=1,mvs
c         dotp=dotp+conjg(x(l))*y(l)
c      enddo
c      mydot=dotp
*
      return
      end
      complex*16 function mychiraldot(x,y)
#include "size.h"
      integer l
      complex*16 dotp, x(ncr,2,mv), y(ncr,2,mv)
*
      dotp=(0.0D0,0.0D0)
      do l=1,mv
         dotp=dotp+dot_product(x(:,2,l),y(:,2,l))
      enddo
      mychiraldot=dotp
*
      return
      end
