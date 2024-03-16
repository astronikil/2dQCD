      subroutine H_wilson(w,v)
#include "size.h"
      integer l
      complex*16 v(ncr,nspr,mv), w(ncr,nspr,mv)
*
      if(fermtype.eq.'D')then
         call fmv_wilson(w,v)
      elseif(fermtype.eq.'M')then
         call fmv_wilson_real(w,v)
      endif
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
      real*8 rm, xaniso
      complex*16 zi
      complex*16 x(ncr,ncr)
      complex*16 u, ud, v(ncr,nspr,mv), w(ncr,nspr,mv)
      complex*16 vx1(ncr), vx2(ncr)
      complex*16 vt1(ncr), vt2(ncr)
      complex*16 wx1(ncr), wx2(ncr)
      complex*16 wt1(ncr), wt2(ncr)
      common /configsmr_irrep/ u(ncr,ncr,mb)
      common /donfigsmr_irrep/ ud(ncr,ncr,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /anisotropy/ xaniso
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
     1    -vt1/xaniso**2 - zi * vt2/xaniso
     1    -vx1 +      vx2
     1    -wt1/xaniso**2 + zi * wt2/xaniso
     1    -wx1 -      wx2

         v(:,2,l) = v(:,2,l)
     1    +zi*vt1/xaniso -   vt2/xaniso**2
     1    +vx1 -      vx2
     1    -zi*wt1/xaniso -   wt2/xaniso**2
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
      real*8 rm, xaniso
      complex*16 zi
      complex*16 x(ncr,ncr)
      complex*16 u, ud, v(ncr,nspr,mv), w(ncr,nspr,mv)
      complex*16 vx1(ncr), vx2(ncr)
      complex*16 vt1(ncr), vt2(ncr)
      complex*16 wx1(ncr), wx2(ncr)
      complex*16 wt1(ncr), wt2(ncr)
      common /configsmr_irrep/ u(ncr,ncr,mb)
      common /donfigsmr_irrep/ ud(ncr,ncr,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /anisotropy/ xaniso
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
     1    -vt1/xaniso**2 + zi * vt2/xaniso
     1    -vx1 -      vx2
     1    -wt1/xaniso**2 - zi * wt2/xaniso
     1    -wx1 +      wx2

         v(:,2,l) = v(:,2,l)
     1   -zi*vt1/xaniso -    vt2/xaniso**2
     1   -vx1 -       vx2
     1   +zi*wt1/xaniso -    wt2/xaniso**2
     1   +wx1 -       wx2
      enddo
!$OMP end parallel do
*
      v = 0.5D0*v + rm*w
*
      return
      end
      subroutine fmv_wilson_real(w,v)
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
      real*8 rm, xaniso
      real*8 u, ud
      real*8 x(ncr,ncr)
      real*8 rw(ncr,nspr,mv), zw(ncr,nspr,mv)
      real*8 rvx1(ncr), rvx2(ncr)
      real*8 rvt1(ncr), rvt2(ncr)
      real*8 rwx1(ncr), rwx2(ncr)
      real*8 rwt1(ncr), rwt2(ncr)

      real*8 zvx1(ncr), zvx2(ncr)
      real*8 zvt1(ncr), zvt2(ncr)
      real*8 zwx1(ncr), zwx2(ncr)
      real*8 zwt1(ncr), zwt2(ncr)

      complex*16 zi
      complex*16 v(ncr,nspr,mv), w(ncr,nspr,mv)

      common /real_configsmr_irrep/ u(ncr,ncr,mb)
      common /real_donfigsmr_irrep/ ud(ncr,ncr,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /anisotropy/ xaniso
      common /quark/ rm
      zi=(0.0D0,1.0D0)
*
      rw=real(w)
      zw=aimag(w)
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
         rvx1 = matmul(x, rw(:, 1, lmu))
         rvx2 = matmul(x, rw(:, 2, lmu))
         zvx1 = matmul(x, zw(:, 1, lmu))
         zvx2 = matmul(x, zw(:, 2, lmu))

         x(:,:) = ud(:,:, l+mux)
         lmu = idn(l+mux)
         rwx1 = matmul(x, rw(:, 1, lmu))
         rwx2 = matmul(x, rw(:, 2, lmu))
         zwx1 = matmul(x, zw(:, 1, lmu))
         zwx2 = matmul(x, zw(:, 2, lmu))

         x(:,:) = u(:,:, l+mut)
         lmu = iup(l+mut)
         rvt1 = matmul(x, rw(:, 1, lmu))
         rvt2 = matmul(x, rw(:, 2, lmu))
         zvt1 = matmul(x, zw(:, 1, lmu))
         zvt2 = matmul(x, zw(:, 2, lmu))

         x(:,:) = ud(:,:, l+mut)
         lmu = idn(l+mut)
         rwt1 = matmul(x, rw(:, 1, lmu))
         rwt2 = matmul(x, rw(:, 2, lmu))
         zwt1 = matmul(x, zw(:, 1, lmu))
         zwt2 = matmul(x, zw(:, 2, lmu))

         v(:,1,l) = v(:,1,l)+
     1     zph*(
     1    -(rvt1+zi*zvt1)/xaniso**2 
     1    -zi*(rvt2+zi*zvt2)/xaniso
     1    -(rvx1+zi*zvx1)
     1    +(rvx2+zi*zvx2)
     1          ) + 
     1     zphconj*(
     1    -(rwt1+zi*zwt1)/xaniso**2 
     1    +zi*(rwt2+zi*zwt2)/xaniso
     1    -(rwx1+zi*zwx1) 
     1    -(rwx2+zi*zwx2)
     1              )

         v(:,2,l) = v(:,2,l)+
     1     zph*(
     1     zi*(rvt1+zi*zvt1)/xaniso 
     1    -(rvt2+zi*zvt2)/xaniso**2
     1    +(rvx1+zi*zvx1) 
     1    -(rvx2+zi*zvx2)
     1          ) + 
     1     zphconj*(
     1    -zi*(rwt1+zi*zwt1)/xaniso 
     1    -(rwt2+zi*zwt2)/xaniso**2
     1    -(rwx1+zi*zwx1) 
     1    -(rwx2+zi*zwx2)
     1           )


      enddo
!$OMP end parallel do
*
      v = 0.5D0*v + rm*w
*
      return
      end
