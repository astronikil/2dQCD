      subroutine pauliinit
#include "size.h"
      complex sig, zi, sid(nc,nc), s2g(2,2,4)
      common/paulimat/sig(2*nc,2*nc,4)
      sid=0.0
      do ic=1,nc
         sid(ic,ic)=1.0
      enddo
*
      zi=(0.0,1.0)
      s2g(1,1,1)=0.0; s2g(1,2,1)=1.0
      s2g(2,1,1)=1.0; s2g(2,2,1)=0.0
*
      s2g(1,1,2)=0.0; s2g(1,2,2)=-zi
      s2g(2,1,2)=zi; s2g(2,2,2)=0.0
*
      s2g(1,1,3)=1.0; s2g(1,2,3)=0.0
      s2g(2,1,3)=0.0; s2g(2,2,3)=-1.0
*
      s2g(1,1,4)=1.0; s2g(1,2,4)=0.0
      s2g(2,1,4)=0.0; s2g(2,2,4)=1.0
*     
      do is=1,4
         call dirprod(sid,s2g(:,:,is),sig(:,:,is))
      enddo
      return
      end
