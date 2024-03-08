      subroutine fermfo(ncgo,ncgw)
************************************************************************
*                                                                      *
*     Molecular dynamics: fermionic part of the force.                 *
*                                                                      *
************************************************************************
#include "size.h"
#include "improve.h"
#include "zolo.h"
      integer ifl, izol, jzol, l, ncgw, ncgs0, ncgo, ncgs1, ncgs2
      integer iup, idn
      integer ic, jc, ica
      real*8 fmo
      real*8 pzov, rzov
      complex*16 ru, rud, p, r, g, w, wrk, fo, qf, fo_ir
      complex*16 rphi, xs, ps, xsov, psov, tgen
      complex*16 u, ud, zi
      complex*16 x, t, zsum
      complex*16 y, z, yt, zt
      complex*16 a1, a2, a3
      complex*16 T1, T2, T3, Tgen_norm, mydot
*
      common /quark_ov/ fmo

      common /gellmann/ tgen(nc,nc,nca)

      common /wrkforce/ ps(mvs,nzol), xs(mvs,nzol)
      common /wrkforce_ovshift/ psov(mvs,nzolmax), xsov(mvs,nzolmax)
      common /scratch1/ z(mvs,nzol)
      common /scratch2/ x(mvs), w(mvs)
      common /scratch3/ wrk(ncr,ncr,mb)
      common /zolvals/ pzov(nzolmax), rzov(nzolmax)
*
      common /forces/ fo(nc,nc,mb), qf(nc,nc,mb)
      common /force_irrep/ fo_ir(ncr,ncr,mb)
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /phi/ rphi(mvs,nf)
*
      common/algebra2/T1(ncr,ncr), T2(ncr,ncr), T3(ncr,ncr)
*
      zi=(0.0D0,1.0D0)
      pzov(1:nzol_n1by2) = p_n1by2(1:nzol_n1by2)
      rzov(1:nzol_n1by2) = r_n1by2(1:nzol_n1by2)
      fo=(0.0D0,0.0D0)
      fo_ir=(0.0D0,0.0D0)
*************************************************************************
*        Prerequisites for force calculation                            *
*************************************************************************
      call ferminit
*************************************************************************
*        Construct R=M.X                                                *
*************************************************************************
      do ifl=1,nf
         call invert_overlap_shifted(rphi(:,ifl),nzol_n1by2,ncgo,ncgw)
         do jzol=0,nzol_n1by2
            if(jzol.eq.0)then
               call invert_wilson_shifted(rphi(:,ifl),ncgs1)
            else
               call invert_wilson_shifted(xsov(:,jzol),ncgs1)
            endif
            z=r_n1by2(jzol)*xs    ! define V

            do izol=1,nzol
               call H_wilson(z(:,izol), w)

               call wdfmv(w, wrk)

               fo_ir=fo_ir+0.25D0*zi*rz(izol)*wrk

               call wdfmv(z(:,izol), wrk)

               fo_ir=fo_ir-0.25D0*zi*rz(izol)*pz(izol)*wrk

            enddo

         enddo
*
      enddo
*
      fo_ir = fo_ir*(1.0D0-fmo*fmo)
*
      ! Hard coded for su(2) irrep conversion

      do l=1,mb

         a1=(0.0D0,0.0D0)
         a2=(0.0D0,0.0D0)
         a3=(0.0D0,0.0D0)

         do ic=1, ncr
            do jc=1, ncr
               a1 = a1 + T1(ic,jc)*(fo_ir(jc,ic,l) +
     1                        conjg(fo_ir(ic,jc,l)))/2.0D0
               a2 = a2 + T2(ic,jc)*(fo_ir(jc,ic,l) +
     1                        conjg(fo_ir(ic,jc,l)))/2.0D0
               a3 = a3 + T3(ic,jc)*(fo_ir(jc,ic,l) +
     1                        conjg(fo_ir(ic,jc,l)))/2.0D0
            enddo
         enddo

         fo(:,:,l) = a1*tgen(:,:,1) + 
     1               a2*tgen(:,:,2) + 
     1               a3*tgen(:,:,3) 
      enddo
      return
      end
      subroutine wdfmv(w, v)
#include "size.h"
#include "improve.h"
      integer l, ic, jc, iup, idn, lmu
      real*8 xaniso
      complex*16 zi, rgu, rgd
      complex*16 x(ncr,ncr)
      complex*16 u, ud, v(ncr,ncr,mb), w(ncr,nspr,mv)
      complex*16 vx1(ncr), vx2(ncr)
      complex*16 vt1(ncr), vt2(ncr)
      common /configsmr_irrep/ u(ncr,ncr,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /anisotropy/ xaniso
      zi=(0.0D0,1.0D0)
*
      v=(0.0D0,0.0D0)
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
         do ic=1, ncr
            do jc=1, ncr
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
