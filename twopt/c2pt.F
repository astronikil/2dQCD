      subroutine correl2pt(ix,iy,it,px,py,c2pt)
#include "size.h"
      integer px, py
      complex gx, gy, gx_back
      complex sig, zi, smat(2*nc,2*nc), zfac
      complex amat(2*nc,2*nc), bmat(2*nc,2*nc), cmat(2*nc,2*nc), 
     1        dmat(2*nc,2*nc), emat(2*nc,2*nc)
      complex c2pt(mt,n2pt)
      common /prop_forward_2pt/ gx(mv,2*nc,2*nc)
      common /prop_backward_2pt/ gx_back(mv,2*nc,2*nc)
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
                  dmat=matmul(smat,amat)
                  emat=matmul(dmat,smat)
                  dmat=matmul(emat,bmat)
                  do ics=1,ncs
                      c2pt(lt,i2pt)=c2pt(lt,i2pt)+
     1                    dmat(ics,ics)*zfac
                  enddo
               enddo
            enddo
         enddo
      enddo
*
      return
      end
      subroutine correl3pt(qx,qy,it_source,sigH,c3pt)
#include "size.h"
      integer qx,qy
      complex gx, gy, gx_back, zfac, sigH(2*nc,2*nc)
      complex sig, zi, smat(2*nc,2*nc), op_phase
      complex c3pt(-mx:mx,mt,n3pt)
      complex amat(2*nc,2*nc), bmat(2*nc,2*nc), cmat(2*nc,2*nc), 
     1        dmat(2*nc,2*nc), emat(2*nc,2*nc), wz(2*nc,2*nc)
      common /prop_forward/ gx(mv,2*nc,2*nc)
      common /prop_backward/ gx_back(mv,2*nc,2*nc)
      common /prop_seq/ gy(mv,2*nc,2*nc)
      common /paulimat/ sig(2*nc,2*nc,4)
*
      tpi=2.0*acos(-1.0)
      zi=(0.0,1.0)
*
      c3pt=0.0
      do it_ins=1,mt
         do iy_ins=1,my
            facy=tpi*iy_ins*qy/my
            do ix_ins=1,mx
               facx=tpi*ix_ins*qx/mx

               op_phase=exp(zi*(facx+facy))

               it_ins_abs=mod(it_source+it_ins-1,mt)+1
               l_ins=(ix_ins+mx*((iy_ins-1)+my*(it_ins_abs-1)))

               !seq prop from source to operator psi(z=0) 
               do ics=1,ncs
                  do jcs=1,ncs
                     amat(ics,jcs)=gy(l_ins,ics,jcs)
                  enddo
               enddo

               do iw=-mx,mx    !wilson line loop from x2 to x2+iw

                  call wilsonline(ix_ins,iy_ins,it_ins_abs,iw,wz)

                  if(iw.ge.0)then
                     ix_ins_iw=mod(ix_ins+iw-1,mx)+1
                  else
                     ix_ins_iw=mod(mx+ix_ins+iw-1,mx)+1
                  endif
                  iy_ins_iw=iy_ins
                  it_ins_abs_iw=it_ins_abs
                  l_ins_iw=
     1            (ix_ins_iw+mx*((iy_ins_iw-1)+my*(it_ins_abs_iw-1)))

                  !backward propagator from x2+iw to source
                  do ics=1,ncs
                     do jcs=1,ncs
                        bmat(ics,jcs)=-conjg(gx_back(l_ins_iw,jcs,ics)) 
                     enddo
                  enddo
*
                  do is1=1,4
                     smat(:,:)=sig(:,:,is1)
                     dmat=matmul(smat,amat)
                     emat=matmul(wz,dmat)
                     dmat=matmul(bmat,emat)
                     emat=matmul(sigH,dmat)
                     do ics=1,ncs
                        c3pt(iw,it_ins,is1)=c3pt(iw,it_ins,is1)+
     1                     emat(ics,ics)*op_phase
                     enddo
                  enddo !loop over operators
               enddo !loop over wilson line lengths
            enddo !loop over ix_o
         enddo !loop over iy_o
      enddo  !loop over it_o
*
      return
      end
      subroutine wilsonline(ix0,iy0,it0,nw,wl)
#include "size.h"
      complex u, uc, sid(2,2)
      complex wg(nc,nc), wgd(nc,nc), wl(ncs,ncs)
      common /configsmr_wil/ u(nc,nc,mb)
      common /configsmr_wil_conj/ uc(nc,nc,mb)
      common /newold/ ino(mv)
      ix=ix0
      iy=iy0
      it=it0
*
      sid=0.0
      sid(1,1)=1.0; sid(2,2)=1.0
*      
      wg=0.0
      do ic=1,nc
         wg(ic,ic)=1.0
      enddo
*      
      if(nw.ge.0)then
        do iw=1,nw
           l=ino(ix+mx*((iy-1)+my*(it-1)))
           wg=matmul(wg,u(:,:,l))
           ix=mod(ix,mx)+1
        enddo
      else
        do iw=1,-nw
           ix=mod(mx+ix-2,mx)+1
           l=ino(ix+mx*((iy-1)+my*(it-1)))
           wg=matmul(wg,uc(:,:,l))
        enddo
      endif

      do ic=1,nc
         do jc=1,nc
            wgd(ic,jc)=conjg(wg(jc,ic))
         enddo
      enddo

      call dirprod(wgd,sid,wl)
      end
