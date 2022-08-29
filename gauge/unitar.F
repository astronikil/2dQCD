      subroutine unitar
#include "size.h"
      real*8, parameter :: tolu=1.0E-8
      integer l, icount
      real*8  rmag, es1, es2
      complex*16 u, a, b
      common /config/ u(nc,nc,mb)
*
      call plaq(es1)
      icount=0
      do l=1,mb
         a=u(1,1,l)
         b=u(1,2,l)

         rmag=sqrt(abs(a)*abs(a)+abs(b)*abs(b))
         if(abs(rmag-1.0).ge.tolu)then
            icount=icount+1
         endif
         a=a/rmag
         b=b/rmag

         u(1,1,l)=a
         u(1,2,l)=b
         u(2,1,l)=-conjg(b)
         u(2,2,l)=conjg(a)
      enddo
      call plaq(es2)
*
      write(27,*) ' '
      write(27,*) 'unitarized: frac of defect links: ',(1.0D0*icount)/mb
      write(27,*) 'es bf unitar: ',es1/dfloat(mv*nc)/3.0D0
      write(27,*) 'es af unitar: ',es2/dfloat(mv*nc)/3.0D0
      write(27,*) ' '
      return
      end
