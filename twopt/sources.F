      subroutine createsrcpts
      common /origin/ iox,ioy,iot
      iox=1
      ioy=1
      iot=1
      end
      subroutine pointsource(v,ix,iy,it,ispin)
#include "size.h"
      complex v(ncs,mv)
      common /newold/ ino(mv)
*
      v=(0.0,0.0)
*
      l=ino(ix+mx*(iy-1+my*(it-1)))
      v(ispin,l)=(1.0,0.0)
*
      return
      end
      subroutine seqsource(y,
     1                     px,
     1                     py,
     1                     it_sink,
     1                     x)
#include "size.h"
      integer px, py
      complex x(ncs,mv), y(ncs,mv), z(ncs,mv), zfac, zi
      common /newold/ ino(mv)
      common /origin/iox,ioy,iot
      zi=(0.0,1.0)
      tpi=2.0*acos(-1.0)
*
      x=0.0
      it=it_sink
      do iy=1,my
         fac1=tpi*py*(iy-ioy)/my
         do ix=1,mx
            j=ino(ix+mx*((iy-1)+my*(it-1)))
            fac2=tpi*px*(ix-iox)/mx
            zfac=exp(zi*(fac1+fac2))
            do ics=1,ncs
               x(ics,j)=y(ics,j)*zfac
            enddo
         enddo
      enddo
*
      return
      end
      subroutine momsource(v,px,py,it,ispin)
#include "size.h"
      integer px, py
      complex v(ncs,mv), zi
      common /newold/ ino(mv)
      tpi=2.0*acos(-1.0)
      zi=(0.0,1.0)
*
      v=(0.0,0.0)
*
      do ix=1,mx
         fac1=tpi*px*ix/mx
         do iy=1,my
            fac2=tpi*py*iy/my
            l=ino(ix+mx*(iy-1+my*(it-1)))
            v(ispin,l)=exp(zi*(fac1+fac2))
         enddo
      enddo
*
      return
      end
