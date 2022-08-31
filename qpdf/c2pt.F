      subroutine correl2pt(ix,iy,it,it00,it0,c2pt)
#include "size.h"
      complex gx, gy
      complex sig(2,2,4), zi, smat(2,2)
      complex amat(2,2), bmat(2,2), cmat(2,2), dmat(2,2), emat(2,2)
      real c2pt(mt,3,n2pt)
      common /prop/ gx(mv,2,2), gy(mv,2,2)
*
*                           * x2=(x,y',t00)
*         *x1=(x,y,t)                   * x3=(x,y,t0)
*
      zi=(0.0,1.0)
      sig(1,1,1)=0.0; sig(1,2,1)=1.0
      sig(2,1,1)=1.0; sig(2,2,1)=0.0
*
      sig(1,1,2)=0.0; sig(1,2,2)=-zi
      sig(2,1,2)=zi; sig(2,2,2)=0.0
*
      sig(1,1,3)=1.0; sig(1,2,3)=0.0
      sig(2,1,3)=0.0; sig(2,2,3)=-1.0
*
      sig(1,1,4)=1.0; sig(1,2,4)=0.0
      sig(2,1,4)=0.0; sig(2,2,4)=1.0
*
      do lt=1,mt
         ! construct G(x0,x1)
         ix3=ix
         iy3=iy
         it3=mod(it+lt-1,mt)+1
         ll=(ix3+mx*((iy3-1)+my*(it3-1)))
         amat(1,1)=gx(ll,1,1)
         amat(1,2)=gx(ll,1,2)
         amat(2,1)=gx(ll,2,1)
         amat(2,2)=gx(ll,2,2)
         bmat(1,1)=conjg(gx(ll,1,1))
         bmat(1,2)=conjg(gx(ll,2,1))
         bmat(2,1)=conjg(gx(ll,1,2))
         bmat(2,2)=conjg(gx(ll,2,2))
         do i2pt=1,n2pt
            smat(1:2,1:2)=sig(1:2,1:2,i2pt)
            dmat=matmul(smat,amat)
            emat=matmul(dmat,smat)
            dmat=matmul(emat,bmat)
            c2pt(lt,1,i2pt)=real(dmat(1,1)+dmat(2,2))
         enddo

         ix2=ix
         iy2=mod(iy+lt-1,my)+1
         it2=it00
         ll=(ix2+mx*((iy2-1)+my*(it2-1)))

         amat(1,1)=gy(ll,1,1)
         amat(1,2)=gy(ll,1,2)
         amat(2,1)=gy(ll,2,1)
         amat(2,2)=gy(ll,2,2)
         bmat(1,1)=conjg(gy(ll,1,1))
         bmat(1,2)=conjg(gy(ll,2,1))
         bmat(2,1)=conjg(gy(ll,1,2))
         bmat(2,2)=conjg(gy(ll,2,2))
         do i2pt=1,n2pt
            smat(1:2,1:2)=sig(1:2,1:2,i2pt)
            dmat=matmul(smat,amat)
            emat=matmul(dmat,smat)
            dmat=matmul(emat,bmat)
            c2pt(lt,2,i2pt)=real(dmat(1,1)+dmat(2,2))
         enddo

         bmat(1,1)=gx(ll,1,1)
         bmat(1,2)=gx(ll,1,2)
         bmat(2,1)=gx(ll,2,1)
         bmat(2,2)=gx(ll,2,2)
         amat(1,1)=conjg(gx(ll,1,1))
         amat(1,2)=conjg(gx(ll,2,1))
         amat(2,1)=conjg(gx(ll,1,2))
         amat(2,2)=conjg(gx(ll,2,2))
         do i2pt=1,n2pt
            smat(1:2,1:2)=sig(1:2,1:2,i2pt)
            dmat=matmul(smat,amat)
            emat=matmul(dmat,smat)
            dmat=matmul(emat,bmat)
            c2pt(lt,3,i2pt)=real(dmat(1,1)+dmat(2,2))
         enddo

      enddo
*
      return
      end
      subroutine correl3pt(ix,iy,it,it00,it0,c3pt)
#include "size.h"
      complex gx, gy
      complex sig(2,2,4), zi, smat(2,2), s2mat(2,2)
      real c3pt(mt,n3pt)
      complex amat(2,2), bmat(2,2), cmat(2,2), dmat(2,2), emat(2,2)
      common /prop/ gx(mv,2,2), gy(mv,2,2)
*
*                           * x2=(x,y',t00)
*         *x1=(x,y,t)                   * x3=(x,y,t0)
*
*
      zi=(0.0,1.0)
      sig(1,1,1)=0.0; sig(1,2,1)=1.0
      sig(2,1,1)=1.0; sig(2,2,1)=0.0
*
      sig(1,1,2)=0.0; sig(1,2,2)=-zi
      sig(2,1,2)=zi;  sig(2,2,2)=0.0
*
      sig(1,1,3)=1.0; sig(1,2,3)=0.0
      sig(2,1,3)=0.0; sig(2,2,3)=-1.0
*
      sig(1,1,4)=1.0; sig(1,2,4)=0.0
      sig(2,1,4)=0.0; sig(2,2,4)=1.0
*
      do lt=1,mt
         ix3=ix
         iy3=iy
         it3=it0
         l3=(ix3+mx*((iy3-1)+my*(it3-1)))
          
         ix2=ix
         iy2=mod(iy+lt-1,my)+1
         it2=it00
         l2=(ix2+mx*((iy2-1)+my*(it2-1)))
         
         amat(1,1)=gx(l3,1,1) 
         amat(1,2)=gx(l3,1,2) 
         amat(2,1)=gx(l3,2,1) 
         amat(2,2)=gx(l3,2,2) 

         bmat(1,1)=gy(l2,1,1) 
         bmat(1,2)=gy(l2,1,2) 
         bmat(2,1)=gy(l2,2,1) 
         bmat(2,2)=gy(l2,2,2) 

         cmat(1,1)=-conjg(gx(l2,1,1)) 
         cmat(1,2)=-conjg(gx(l2,2,1)) 
         cmat(2,1)=-conjg(gx(l2,1,2)) 
         cmat(2,2)=-conjg(gx(l2,2,2)) 

         do i3pt=1,n3pt
         smat(1:2,1:2)=sig(1:2,1:2,i3pt)
         s2mat(1:2,1:2)=sig(1:2,1:2,3)
         dmat=matmul(smat,amat)
         emat=matmul(dmat,smat)
         dmat=matmul(emat,bmat)
         emat=matmul(dmat,s2mat)
         dmat=matmul(emat,cmat)
         c3pt(lt,i3pt)=real(dmat(1,1)+dmat(2,2))*sqrt(2.0)
         enddo
      enddo
*
      return
      end
