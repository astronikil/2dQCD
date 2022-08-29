      subroutine H_wilson(w,v)
#include "size.h"
      integer l
      complex*16 v(nc,nspr,mv), w(nc,nspr,mv)
*
      call fmv_wilson(w,v)
*
      do l=1,mv
         v(:,2,l) = -1.0D0*v(:,2,l) 
      enddo
*
      end
      subroutine fmv_wilson(w,v)
*************************************************************************
*                                                                       *
*     The fermion matrix multiplying a vector: v = M w                  *
*                                                                       *
* Mathematica code snippets for explicit expressions in spinor space    *
*                                                                       *
*                                                                       *
* g1 = {{0, 0, 0, -I}, {0, 0, -I, 0}, {0, I, 0, 0}, {I, 0, 0, 0}};      *
* g2 = {{0, 0, 0, -1}, {0, 0, 1, 0}, {0, 1, 0, 0}, {-1, 0, 0, 0}};      *
* g3 = {{0, 0, -I, 0}, {0, 0, 0, I}, {I, 0, 0, 0}, {0, -I, 0, 0}};      *
* g4 = {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, 0}};        *
* id = IdentityMatrix[4];                                               *
*                                                                       *
* q = (g1 - id).vx - (g1 + id).wx + (g2 - id).vy - (g2 + id).wy + (g3 - * 
*     id).vz - (g3 + id).wz + (g4 - id).vt - (g4 + id).wt               *
*                                                                       *
* (g3.g4 - g4.g3).vx/2/I                                                *
*************************************************************************
#include "size.h"
#include "improve.h"
      integer l, ic, jc, iup, idn, lmu
      real*8 rm
      complex*16 zi
      complex*16 x(nc,nc)
      complex*16 u, ud, v(nc,nspr,mv), w(nc,nspr,mv)
      complex*16 vx1(nc), vx2(nc)
      complex*16 vt1(nc), vt2(nc)
      complex*16 wx1(nc), wx2(nc)
      complex*16 wt1(nc), wt2(nc)
      common /configsmr/ u(nc,nc,mb)
      common /donfigsmr/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /quark/ rm
      zi=(0.0D0,1.0D0)
*
*************************************************************************
*                                                                       *
*     Action of the fermion matrix.                                     *
*                                                                       *
*************************************************************************
      v=0.0D0
!$OMP parallel do default(shared)
!$OMP^ private(x)
!$OMP^ private(lmu, l)
!$OMP^ private(vx1, vx2)
!$OMP^ private(vt1, vt2)
!$OMP^ private(wx1, wx2)
!$OMP^ private(wt1, wt2)
      do l=1, mv
         x(:,:) = u(:,:, l+mux)
         lmu = iup(l+mux)
         vx1 = matmul(x, w(:, 1, lmu))
         vx2 = matmul(x, w(:, 2, lmu))

         x(:,:) = ud(:,:, l+mux)
         lmu = idn(l+mux)
         wx1 = matmul(x, w(:, 1, lmu))
         wx2 = matmul(x, w(:, 2, lmu))

         x(:,:) = u(:,:, l+mut)
         lmu = iup(l+mut)
         vt1 = matmul(x, w(:, 1, lmu))
         vt2 = matmul(x, w(:, 2, lmu))

         x(:,:) = ud(:,:, l+mut)
         lmu = idn(l+mut)
         wt1 = matmul(x, w(:, 1, lmu))
         wt2 = matmul(x, w(:, 2, lmu))

         v(:,1,l) = v(:,1,l)
     1    -vt1 - zi * vt2
     1    -vx1 +      vx2
     1    -wt1 + zi * wt2
     1    -wx1 -      wx2

         v(:,2,l) = v(:,2,l)
     1    +zi*vt1 -   vt2
     1    +vx1 -      vx2
     1    -zi*wt1 -   wt2
     1    -wx1 -      wx2


      enddo
!$OMP end parallel do
*
      v = 0.5D0*v + rm*w
*
      return
      end
      subroutine fmtv_wilson(w,v)
*************************************************************************
*                                                                       *
*     The fermion matrix multiplying a vector: v = M w                  *
*                                                                       *
*************************************************************************
#include "size.h"
#include "improve.h"
      integer l, ic, jc, iup, idn, lmu
      real*8 rm
      complex*16 zi
      complex*16 x(nc,nc)
      complex*16 u, ud, v(nc,nspr,mv), w(nc,nspr,mv)
      complex*16 vx1(nc), vx2(nc)
      complex*16 vt1(nc), vt2(nc)
      complex*16 wx1(nc), wx2(nc)
      complex*16 wt1(nc), wt2(nc)
      common /configsmr/ u(nc,nc,mb)
      common /donfigsmr/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /quark/ rm
      zi=(0.0D0,1.0D0)
*
*************************************************************************
*                                                                       *
*     Action of the fermion matrix.                                     *
*                                                                       *
*************************************************************************
      v=0.0
!$OMP parallel do default(shared)
!$OMP^ private(x)
!$OMP^ private(lmu, l)
!$OMP^ private(vx1, vx2)
!$OMP^ private(vt1, vt2)
!$OMP^ private(wx1, wx2)
!$OMP^ private(wt1, wt2)
      do l=1, mv
         x(:,:) = u(:,:, l+mux)
         lmu = iup(l+mux)
         vx1 = matmul(x, w(:, 1, lmu))
         vx2 = matmul(x, w(:, 2, lmu))

         x(:,:) = ud(:,:, l+mux)
         lmu = idn(l+mux)
         wx1 = matmul(x, w(:, 1, lmu))
         wx2 = matmul(x, w(:, 2, lmu))

         x(:,:) = u(:,:, l+mut)
         lmu = iup(l+mut)
         vt1 = matmul(x, w(:, 1, lmu))
         vt2 = matmul(x, w(:, 2, lmu))

         x(:,:) = ud(:,:, l+mut)
         lmu = idn(l+mut)
         wt1 = matmul(x, w(:, 1, lmu))
         wt2 = matmul(x, w(:, 2, lmu))

         v(:,1,l) = v(:,1,l)
     1    -vt1 + zi * vt2
     1    -vx1 -      vx2
     1    -wt1 - zi * wt2
     1    -wx1 +      wx2

         v(:,2,l) = v(:,2,l)
     1   -zi*vt1 -    vt2
     1   -vx1 -       vx2
     1   +zi*wt1 -    wt2
     1   +wx1 -       wx2
      enddo
!$OMP end parallel do
*
      v = 0.5D0*v + rm*w
*
      return
      end
