      subroutine correl2pt(ix,it,c2pt)
#include "size.h"
#include "param.h"
      integer ic, jc, ll, ix, it, ix3, it3, lt, i2pt
      complex*16 gx
      complex*16 zi, smat(ncs,ncs), sig
      complex*16 amat(ncs,ncs), bmat(ncs,ncs), cmat(ncs,ncs), 
     1        dmat(ncs,ncs), emat(ncs,ncs)
      complex*16 c2pt(mt,n2pt)
      common /paulimat/ sig(ncs,ncs,4)
      common /prop/ gx(mv,ncs,ncs)
*
      zi=(0.0D0,1.0D0)
*
      do lt=1,mt
         ix3=ix
         it3=mod(it+lt-1,mt)+1
         ll=(ix3+mx*(it3-1))

         do ic=1,ncs
            do jc=1,ncs
               amat(ic,jc)=gx(ll,ic,jc)
               bmat(ic,jc)=conjg(gx(ll,jc,ic))
            enddo
         enddo
*
         do i2pt=1,n2pt
            smat(:,:)=sig(:,:,i2pt)
            dmat=matmul(smat,amat)
            emat=matmul(dmat,smat)
            dmat=matmul(emat,bmat)
            c2pt(lt,i2pt)=(0.0D0,0.0D0)
            do ic=1,ncr
               c2pt(lt,i2pt)=c2pt(lt,i2pt) + dmat(ic,ic)
            enddo
         enddo
*
      enddo
*
      return
      end
