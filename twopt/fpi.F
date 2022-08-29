      subroutine fpi2pt(ix,iy,it,px,py,c2pt)
#include "size.h"
      integer px, py
      complex gx, gy, gx_back
      complex sig, zi, smat(2*nc,2*nc), zfac
      complex amat(2*nc,2*nc), bmat(2*nc,2*nc), cmat(2*nc,2*nc), 
     1        dmat(2*nc,2*nc), emat(2*nc,2*nc)
      complex c2pt(mt,n2pt)
      common /prop_forward/ gx(mv,2*nc,2*nc)
      common /prop_backward/ gx_back(mv,2*nc,2*nc)
      common /prop_seq/ gy(mv,2*nc,2*nc)
      common /paulimat/ sig(2*nc,2*nc,4)
      tpi=2.0*acos(-1.0)
      zi=(0.0,1.0)
*
      c2pt=0.0
      do lt=1,mt
         do iy3=1,my
            fac2=tpi*(iy3-iy)*py/my
            do ix3=1,mx
               fac1=tpi*(ix3-ix)*px/mx
               zfac=exp(zi*(fac1+fac2))
               it3=mod(it+lt-1,mt)+1
               ll=(ix3+mx*((iy3-1)+my*(it3-1)))
               do ics=1,ncs
                  do jcs=1,ncs
                     amat(ics,jcs)=gx(ll,ics,jcs)
                     bmat(ics,jcs)=conjg(gx_back(ll,jcs,ics))
                  enddo
               enddo
               do i2pt=1,n2pt
                  smat(:,:)=sig(:,:,i2pt)
                  dmat=matmul(amat,smat)
                  emat=matmul(dmat,bmat)
                  do ics=1,ncs
                      c2pt(lt,i2pt)=c2pt(lt,i2pt)+
     1                    emat(ics,ics)*zfac
                  enddo
               enddo
            enddo
         enddo
      enddo
*
      return
      end
      subroutine prcfpi
#include "size.h"
#include "mom.h"
      complex f2ptex
      common /fpi2ex/f2ptex(mt,n2pt,nqmom,nmom)
      do imom=1,nmom
         do iqmom=1,nqmom
            do ic=1,n2pt
               do i=1,mt
                  write(37,*)i,
     1             real(f2ptex(i,ic,iqmom,imom)), 
     1             aimag(f2ptex(i,ic,iqmom,imom)), 
     1             ip2str(ic), 
     1             ' px=',phad(1,imom), ' py=',phad(2,imom),
     1             ' qx=',qins(1,iqmom),' qy=',qins(2,iqmom)
               enddo
            enddo
         enddo
      enddo
      end
