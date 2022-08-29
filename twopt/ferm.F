      subroutine ferminit
      parameter(nlev_sea=1)
*************************************************************************
*                                                                       *
*     Initializations before fermion computation.                       *
*                                                                       *
*************************************************************************
      call stout(nlev_sea)
      !call clover
      return
      end subroutine ferminit
      subroutine initu
#include "size.h"
      parameter(nlev_wil=1)
      complex*16 u, uw, uwc, utmp(nc,nc,mb)
      common /configsmr/ u(nc,nc,mb)
      common /configsmr_wil/ uw(nc,nc,mb)
      common /configsmr_wil_conj/ uwc(nc,nc,mb)
      call linkc

      utmp=u
      call stout(nlev_wil)
      uw=u
      u=utmp

      call hconj(uw,uwc,mb)
      end 
      subroutine fmv(w,v)
*************************************************************************
*                                                                       *
*     The fermion matrix multiplying a vector: v = M w                  *
*                                                                       *
*************************************************************************
#include "size.h"
      integer l, ic, jc, iup, idn
      real*8 rm
      complex*16 zi
      complex*16 p1,p2
      complex*16 t11,t12,t21,t22,t31,t32
      complex*16 u, ud, v(2,nc,mv), w(2,nc,mv), v1, v2
      common /configsmr/ u(nc,nc,mb)
      common /donfig/ ud(nc,nc,mb)
      common /iupidn/ iup(mv,3), idn(mv,3)
      common /quark/ rm
      zi=(0.0D0,1.0D0)
      v=(0.0D0,0.0D0)
*
      call linkc
*************************************************************************
*                                                                       *
*     Action of the fermion matrix.                                     *
*                                                                       *
*************************************************************************
!$OMP parallel do default(shared)
!$OMP^ private(p1,p2)
!$OMP^ private(t11,t12,t21,t22,t31,t32)
!$OMP^ private(v1,v2)
      do l = 1, mv
         do ic=1, nc
            v(1,ic,l) = rm*w(1,ic,l)
            v(2,ic,l) = rm*w(2,ic,l)
            do jc=1, nc
*
               t11=u(ic,jc,l+mux)*(w(1,jc,iup(l,1))+w(2,jc,iup(l,1)))
*
               t12=ud(ic,jc,l+mux)*(-w(1,jc,idn(l,1))+w(2,jc,idn(l,1)))
*
               t21=u(ic,jc,l+muy)*(w(1,jc,iup(l,2))-zi*w(2,jc,iup(l,2)))
*
               t22=ud(ic,jc,l+muy)*(-w(1,jc,idn(l,2))-
     1             zi*w(2,jc,idn(l,2)))
*
               t31=2.0D0*u(ic,jc,l+mut)*w(1,jc,iup(l,3))
*
               t32=-2.0D0*ud(ic,jc,l+mut)*w(2,jc,idn(l,3))
*
               v(1,ic,l) = v(1,ic,l)+
     *            0.5D0*(
     *                  t11-t12
     *                 +t21-t22
     *                 +t31
     *           )
*
              v(2,ic,l) = v(2,ic,l)+
     *             0.5D0*(
     *                  t11+t12
     *                 +zi*t21+zi*t22
     *                 -t32
     *           )
            enddo
         enddo
      enddo
*
      return
      end
      subroutine fmtv(w,v)
*************************************************************************
*                                                                       *
*     The fermion matrix multiplying a vector: v = M^+ w                *
*                                                                       *
*************************************************************************
#include "size.h"
      integer l, ic, jc, iup, idn
      real*8 rm
      complex*16 zi
      complex*16 p1,p2
      complex*16 t11,t12,t21,t22,t31,t32
      complex*16 u, ud, v(2,nc,mv), w(2,nc,mv), v1, v2
      common /configsmr/ u(nc,nc,mb)
      common /donfig/ ud(nc,nc,mb)
      common /iupidn/ iup(mv,3), idn(mv,3)
      common /quark/ rm
      zi=(0.0D0,1.0D0)
      v=(0.0D0,0.0D0)
*
      call linkc
*************************************************************************
*                                                                       *
*     Action of the fermion matrix.                                     *
*                                                                       *
*************************************************************************
!$OMP parallel do default(shared)
!$OMP^ private(p1,p2)
!$OMP^ private(t11,t12,t21,t22,t31,t32)
!$OMP^ private(v1,v2,l)
      do l=1,mv
         do ic=1, nc
            v(1,ic,l) = rm*w(1,ic,l)
            v(2,ic,l) = rm*w(2,ic,l)
            do jc=1, nc
*
               t11=u(ic,jc,l+mux)*(-w(1,jc,iup(l,1))+w(2,jc,iup(l,1)))
*
               t12=ud(ic,jc,l+mux)*(w(1,jc,idn(l,1))+w(2,jc,idn(l,1)))
*
               t21=u(ic,jc,l+muy)*(-w(1,jc,iup(l,2))-
     1             zi*w(2,jc,iup(l,2)))
*
               t22=ud(ic,jc,l+muy)*(w(1,jc,idn(l,2))-
     1             zi*w(2,jc,idn(l,2)))
*
               t31=-2.0D0*u(ic,jc,l+mut)*w(2,jc,iup(l,3))
*
               t32=2.0D0*ud(ic,jc,l+mut)*w(1,jc,idn(l,3))
*
               v(1,ic,l) = v(1,ic,l)
     *           - 0.5D0*(
     *                  t11-t12
     *                 +t21-t22
     *                     -t32
     *           )
*
               v(2,ic,l) = v(2,ic,l)
     *           - 0.5D0*(
     *                 -t11-t12
     *                 -zi*t21-zi*t22
     *                 +t31
     *           )
            enddo
         enddo
      enddo
*
      return
      end
