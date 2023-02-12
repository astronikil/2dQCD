      subroutine dbydaniso(dF)
************************************************************************
*                                                                      *
*     Derivative of derm determinant wrt anistropy xi                  *
*                                                                      *
************************************************************************
#include "size.h"
#include "improve.h"
#include "zolo.h"
      integer ifl, izol, l, ncgs, ncgs0, ncgo, ncgs1, ncgs2
      integer iup, idn
      real*8 dF
      complex*16 ru, rud, p, r, g, w, wrk, fo, qf
      complex*16 rphi, xs, ps
      complex*16 u, ud, zi
      complex*16 x, t, zsum
      complex*16 y, z, yt, zt
      complex*16 mydot
*      
      common /wrkforce/ ps(mvs,nzol), xs(mvs,nzol)
      common /scratch1/ z(mvs,nzol)
      common /scratch2/ x(mvs), w(mvs)
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /phi/ rphi(mvs,nf)
*
      zi=(0.0D0,1.0D0)
*      
      call ferminit
*      
      dF = 0.0D0
      do ifl=1,nf
         call invert_overlap(rphi(:,ifl), x, ncgo,ncgs0)         
         call invert_wilson_shifted(x,ncgs1)
         z=xs    ! define Y_k

         do izol=1,nzol

            call H_wilson(z(:,izol), w) ! define W_k

            call dfmvbyaniso(w, x)

            dF=dF+rz(izol)*real(mydot(w,x))

            call dfmvbyaniso(z(:,izol), x)

            dF=dF-rz(izol)*pz(izol)*real(mydot(z(:,izol),x))

         enddo
*
      enddo
      return
      end
      subroutine dbydaniso_stoch(dF, er, nvec)
************************************************************************
*                                                                      *
*     Derivative of derm determinant wrt anistropy xi                  *
*                                                                      *
************************************************************************
#include "size.h"
#include "improve.h"
#include "zolo.h"
      integer ifl, izol, l, ncgs, ncgs0, ncgo, ncgs1, ncgs2
      integer iup, idn, ivec, nvec
      real*8 dF, av, er
      complex*16 ru, rud, p, r, g, w, wrk, fo, qf
      complex*16 rphi, xs, ps, rv
      complex*16 u, ud, zi, zR, zL
      complex*16 x, t, zsum, x1, w1
      complex*16 y, z, yt, zt
      complex*16 mydot
*      
      common /wrkforce/ ps(mvs,nzol), xs(mvs,nzol)
      common /scratch12/ zR(mvs,nzol), zL(mvs,nzol)
      common /scratch2/ x(mvs), w(mvs)
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /stochv/ rv(mvs)
*
      zi=(0.0D0,1.0D0)
*      
      call ferminit
*      
      dF = 0.0D0
      er = 0.0D0
      do ivec=1,nvec
         call z2noise(rv)

         call invert_wilson_shifted(rv,ncgs1)
         zR=xs    ! define Y_k

         call invert_overlap(rv, w, ncgo,ncgs0)
         call H_overlap(w,x,ncgs)
         call invert_wilson_shifted(x,ncgs1)
         zL=xs    ! define Y_k

         av=0.0D0
         do izol=1,nzol

            call H_wilson(zR(:,izol), w) ! define W_k

            call dfmvbyaniso(w, x)

            call H_wilson(x, w)

            av=av+rz(izol)*real(mydot(zL(:,izol),w))

            call dfmvbyaniso(zR(:,izol), x)

            av=av-rz(izol)*pz(izol)*real(mydot(zL(:,izol),x))


         enddo
         er=er+av**2
         dF=dF+av
      enddo
*
      er =
     1  nf*sqrt(er/dfloat(nvec)-dF**2/dfloat(nvec*nvec))/sqrt(nvec*1.0) 
      dF = nf*dF/dfloat(nvec)  ! factor 2 due to taking real part
      return
      end
      subroutine dfmvbyaniso(w, v)
#include "size.h"
#include "improve.h"
      integer l, iup, idn, lmu
      real*8 xaniso
      complex*16 zi
      complex*16 x(nc,nc)
      complex*16 u, ud, v(nc,nspr,mv), w(nc,nspr,mv)
      complex*16 vt1(nc), vt2(nc)
      common /configsmr/ u(nc,nc,mb)
      common /donfigsmr/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /anisotropy/ xaniso
      zi=(0.0D0,1.0D0)
*
      do l=1, mv
         v(:,1,l)=-1.0D0/xaniso**3*w(:,1,l)
         v(:,2,l)= 1.0D0/xaniso**3*w(:,2,l)
*
         x(:,:) = u(:,:, l+mut)
         lmu = iup(l+mut)
         vt1 = matmul(x, w(:, 1, lmu))
         vt2 = matmul(x, w(:, 2, lmu))
*
         v(:,1,l)=v(:,1,l)+0.25D0/xaniso**2*(zi*vt2+2.0D0*vt1/xaniso)
         v(:,2,l)=v(:,2,l)+0.25D0/xaniso**2*(zi*vt1-2.0D0*vt2/xaniso)
*
         x(:,:) = ud(:,:, l+mut)
         lmu = idn(l+mut)
         vt1 = matmul(x, w(:, 1, lmu))
         vt2 = matmul(x, w(:, 2, lmu))
*
         v(:,1,l)=v(:,1,l)+0.25D0/xaniso**2*(-zi*vt2+2.0D0*vt1/xaniso)
         v(:,2,l)=v(:,2,l)+0.25D0/xaniso**2*(-zi*vt1-2.0D0*vt2/xaniso)
      enddo
      return
      end
