      subroutine fermfo(ncgo,ncgs)
************************************************************************
*                                                                      *
*     Molecular dynamics: fermionic part of the force.                 *
*                                                                      *
************************************************************************
#include "size.h"
#include "improve.h"
#include "zolo.h"
      integer ifl, izol, l, ncgs, ncgs0, ncgo, ncgs1, ncgs2
      integer iup, idn
      complex*16 ru, rud, p, r, g, w, wrk, fo, qf
      complex*16 rphi, xs, ps
      complex*16 u, ud, zi
      complex*16 x, t, zsum
      complex*16 y, z, yt, zt
*      
      !v(nc,nspr,mv)
      common /wrkforce/ ps(mvs,nzol), xs(mvs,nzol)
      common /scratch1/ z(mvs,nzol)
      common /scratch2/ x(mvs), w(mvs)
      common /scratch3/ wrk(ntot)
*
      common /forces/ fo(ntot), qf(ntot)
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /phi/ rphi(mvs,nf)
*
      zi=(0.0D0,1.0D0)
      fo=0.0D0
*************************************************************************
*        Prerequisites for force calculation                            *
*************************************************************************
      call ferminit
*************************************************************************
*        Construct R=M.X                                                *
*************************************************************************
      do ifl=1,nf
         call invert_overlap(rphi(:,ifl), x, ncgo,ncgs0)         
         call invert_wilson_shifted(x,ncgs1)
         z=xs    ! define V

         do izol=1,nzol
            call H_wilson(z(:,izol), w)

            call wdfmv(w, wrk)

            fo=fo+0.25D0*zi*rz(izol)*wrk

            call wdfmv(z(:,izol), wrk)

            fo=fo-0.25D0*zi*rz(izol)*pz(izol)*wrk

         enddo
*
      enddo
      ncgs=max(ncgs0,ncgs1)
*
      return
      end
      subroutine wdfmv(w, v)
#include "size.h"
#include "improve.h"
      integer l, ic, jc, iup, idn, lmu
      real*8 xaniso
      complex*16 zi, rgu, rgd
      complex*16 x(nc,nc)
      complex*16 u, ud, v(nc,nc,mb), w(nc,nspr,mv)
      complex*16 vx1(nc), vx2(nc)
      complex*16 vt1(nc), vt2(nc)
      common /configsmr/ u(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /anisotropy/ xaniso
      zi=(0.0D0,1.0D0)
*
      v=0.0
      do l=1, mv
         x(:,:) = u(:,:, l+mux)
         lmu = iup(l+mux)
         vx1 = matmul(x, w(:, 1, lmu))
         vx2 = matmul(x, w(:, 2, lmu))

         x(:,:) = u(:,:, l+mut)
         lmu = iup(l+mut)
         vt1 = matmul(x, w(:, 1, lmu))
         vt2 = matmul(x, w(:, 2, lmu))

         ! spin trace over (1+\sigma_i).sigma_3.v
         do ic=1, nc
            do jc=1, nc
               v(ic,jc,l+mux)=
     1         conjg(w(jc, 1, l))*(vx1(ic)-vx2(ic)) 
     1        +conjg(w(jc, 2, l))*(vx1(ic)-vx2(ic))

               v(ic,jc,l+mut)=
     1         conjg(w(jc, 1, l))*(vt1(ic)/xaniso**2+zi*vt2(ic)/xaniso) 
     1        +conjg(w(jc, 2, l))*(zi*vt1(ic)/xaniso-vt2(ic)/xaniso**2)
            enddo 
         enddo
      enddo
      return
      end
