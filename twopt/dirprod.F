      subroutine dirprod(a,b,c)
#include "size.h"
      complex a(nc,nc), b(2,2), c(2*nc,2*nc)
      ics=0
      do ic=1,nc
         do is=1,2
            ics=ics+1
            jcs=0
            do jc=1,nc
               do js=1,2
                  jcs=jcs+1
                  c(jcs,ics)=a(jc,ic)*b(js,is)
               enddo
            enddo
         enddo
      enddo
      return
      end
