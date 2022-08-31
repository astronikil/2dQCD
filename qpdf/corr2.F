      subroutine correl2pt(ix,iy,it,it00,it0,c2pt)
#include "size.h"
      complex gx, gy
      real c2pt(mt,3,n2pt)
      common /prop/ gx(mv,2,2), gy(mv,2,2)
*
*                           * x2=(x,y',t00)
*         *x1=(x,y,t)                   * x3=(x,y,t0)
*
      do lt=1,mt
         ! construct G(x0,x1)
         ix3=ix
         iy3=iy
         it3=mod(it+lt-1,mt)+1
         ll=(ix3+mx*((iy3-1)+my*(it3-1)))
         c2pt(lt,1,1)=abs(gx(ll,1,1))**2+abs(gx(ll,2,2))**2
     1               +abs(gx(ll,1,2))**2+abs(gx(ll,2,1))**2
         c2pt(lt,1,2)=real(gx(ll,1,1)*conjg(gx(ll,2,2)))

         ix2=ix
         iy2=mod(iy+lt-1,my)+1
         it2=it00
         ll=(ix2+mx*((iy2-1)+my*(it2-1)))

         c2pt(lt,2,1)=abs(gy(ll,1,1))**2+abs(gy(ll,2,2))**2
     1               +abs(gy(ll,1,2))**2+abs(gy(ll,2,1))**2
         c2pt(lt,2,2)=real(gy(ll,1,1)*conjg(gy(ll,2,2)))


         c2pt(lt,3,1)=abs(gx(ll,1,1))**2+abs(gx(ll,2,2))**2
     1               +abs(gx(ll,1,2))**2+abs(gx(ll,2,1))**2
         c2pt(lt,3,2)=real(gx(ll,1,1)*conjg(gx(ll,2,2)))


      enddo
*
      return
      end
