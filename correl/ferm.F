      subroutine ferminit
*************************************************************************
*                                                                       *
*     Initializations before fermion computation.                       *
*                                                                       *
*************************************************************************
      call hyp_gauge
      call clover
*
      return
      end subroutine ferminit
      subroutine initcompactu
#include "size.h"
      complex u, zi
      real ua, q, pifac
      integer l
      common /configsmr/ u(mb)
      common /configasmr/ ua(mb)
      common /cftdata/ q, QQ
      zi=(0.0D0,1.0D0)
      pifac=acos(-1.0)/mx !antiperiodic bc in all three directions
      do l=1,mb
         u(l)=exp(zi*(ua(l)*q*QQ+pifac))
      enddo
      call linkc
      end subroutine initcompactu
      subroutine fmv_wilson(w,v)
*************************************************************************
*                                                                       *
*     The fermion matrix multiplying a vector: v = M w                  *
*                                                                       *
*************************************************************************
#include "size.h"
      complex zi
      complex p1,p2
      complex t11,t12,t21,t22,t31,t32
      complex u, ud, v(2,mv), w(2,mv), v1, v2
      real cv, cvi
      common /configsmr/ u(mb)
      common /donfig/ ud(mb)
      common /clvr/cv(mb), cvi(mb)
      common /iupidn/ iup(mv,3), idn(mv,3)
      common /quark/ rm
      zi=(0.,1.)
      v=(0.,0.)
*************************************************************************
*                                                                       *
*     Action of the fermion matrix.                                     *
*                                                                       *
*************************************************************************
!$OMP parallel do default(shared)
!$OMP^ private(p1,p2)
!$OMP^ private(t11,t12,t21,t22,t31,t32)
!$OMP^ private(v1,v2)
      do  l = 1, mv

          v1 = rm*w(1,l)
          v2 = rm*w(2,l)
*
          t11=u(l+mux)*(w(1,iup(l,1)) + w(2,iup(l,1)))
*
          t12=ud(l+mux)*(-w(1,idn(l,1)) + w(2,idn(l,1)))
*
          t21=u(l+muy)*(w(1,iup(l,2)) - zi*w(2,iup(l,2)))
*
          t22=ud(l+muy)*(-w(1,idn(l,2)) - zi*w(2,idn(l,2)))
*
          t31=2*u(l+mut)*w(1,iup(l,3))
*
          t32=-2*ud(l+mut)*w(2,idn(l,3))
*
          v(1,l) = v1
     *           + 0.5*(
     *                  t11-t12
     *                 +t21-t22
     *                 +t31
     *           )
     *           -cv(l)*w(2,l)+zi*cv(l+muy)*w(2,l)-cv(l+mut)*w(1,l)
*
          v(2,l) = v2
     *           + 0.5*(
     *                  t11+t12
     *                 +zi*t21+zi*t22
     *                 -t32
     *           )
     *           -cv(l)*w(1,l)-zi*cv(l+muy)*w(1,l)+cv(l+mut)*w(2,l)

      enddo
*
      return
      end
      subroutine fmtv_wilson(w,v)
*************************************************************************
*                                                                       *
*     The fermion matrix multiplying a vector: v = M^+ w                *
*                                                                       *
*************************************************************************
#include "size.h"
      complex zi
      complex p1,p2
      complex t11,t12,t21,t22,t31,t32
      complex u, ud, v(2,mv), w(2,mv), v1, v2
      real cv, cvi
      common /configsmr/ u(mb)
      common /donfig/ ud(mb)
      common /clvr/cv(mb), cvi(mb)
      common /iupidn/ iup(mv,3), idn(mv,3)
      common /quark/ rm
      zi=(0.,1.)
      v=(0.,0.)
*************************************************************************
*                                                                       *
*     Action of the fermion matrix.                                     *
*                                                                       *
*************************************************************************
!$OMP parallel do default(shared)
!$OMP^ private(p1,p2)
!$OMP^ private(t11,t12,t21,t22,t31,t32)
!$OMP^ private(v1,v2,l)
      do  l = 1, mv
          v1 = rm*w(1,l)
          v2 = rm*w(2,l)
*
          t11=u(l+mux)*(-w(1,iup(l,1)) + w(2,iup(l,1)))
*
          t12=ud(l+mux)*(w(1,idn(l,1)) + w(2,idn(l,1)))
*
          t21=u(l+muy)*(-w(1,iup(l,2)) - zi*w(2,iup(l,2)))
*
          t22=ud(l+muy)*(w(1,idn(l,2)) - zi*w(2,idn(l,2)))
*
          t31=-2*u(l+mut)*w(2,iup(l,3))
*
          t32=2*ud(l+mut)*w(1,idn(l,3))
*
          v(1,l) = v1
     *           - 0.5*(
     *                  t11-t12
     *                 +t21-t22
     *                     -t32
     *           )
     *           -cv(l)*w(2,l)+zi*cv(l+muy)*w(2,l)-cv(l+mut)*w(1,l)
*
          v(2,l) = v2
     *           - 0.5*(
     *                 -t11-t12
     *                 -zi*t21-zi*t22
     *                 +t31
     *           )
     *           -cv(l)*w(1,l)-zi*cv(l+muy)*w(1,l)+cv(l+mut)*w(2,l)
                
      enddo
*
      return
      end
      subroutine fmtmv_overlap(v,w,ncg)
#include "size.h"
      complex v(mvs), w(mvs), x, y
      common /fermwrk1/ x(mvs), y(mvs)
      common /qmass_overlap/fmo,fmfac1,fmfac2
*
      call fmtv_wilson(v,x)
      call sqrtoverlap(x,y,ncg1)
      w=y
      call sqrtoverlap(v,x,ncg2)
      call fmv_wilson(x,y)
      w=w+y
*
      w=fmfac1*v+fmfac2*w
*
      ncg=max(ncg1,ncg2)
      return
      end
      subroutine fmtv_overlap(v,w,ncg)
#include "size.h"
      complex v(mvs), w(mvs), x, y
      common /fermwrk1/ x(mvs), y(mvs)
      common /qmass_overlap/fmo,fmfac1,fmfac2
*
      r1=(1.-fmo)/2.
      r2=(1.+fmo)/2.
      call fmtv_wilson(v,x)
      call sqrtoverlap(x,y,ncg1)
*
      w=r1*v-r2*y
*
      ncg=max(ncg1,ncg1)
      return
      end
      subroutine fmv_overlap(v,w,ncg)
#include "size.h"
      complex v(mvs), w(mvs), x, y
      common /fermwrk1/ x(mvs), y(mvs)
      common /qmass_overlap/fmo,fmfac1,fmfac2
      r1=(1.-fmo)/2.
      r2=(1.+fmo)/2.
*
      call sqrtoverlap(v,x,ncg1)
      call fmv_wilson(x,y)
*
      w=r1*v-r2*y
*
      ncg=max(ncg1,ncg1)
      return
      end
      subroutine sqrtoverlap(v,w,ncg)
#include "size.h"
#include "zolo.h"
      complex v(mvs), w(mvs), x, ps
      common /wrkforce/ ps(mvs,nzol), x(mvs,nzol)

      call invert_wilson_shifted(v,ncg)
      w=0.0
      do izol=1,nzol
         w=w+rz(izol)*x(:,izol)
      enddo
      end
      subroutine vmat(v,w,ncg)
#include "size.h"
      complex v(mvs), w(mvs), x, y
      common /fermwrk1/ x(mvs), y(mvs)
      common /qmass_overlap/fmo,fmfac1,fmfac2
*
      call sqrtoverlap(v,x,ncg)
      call fmv_wilson(x,w)
*
      return
      end
      subroutine vdagmat(v,w,ncg)
#include "size.h"
      complex v(mvs), w(mvs), x, y
      common /fermwrk1/ x(mvs), y(mvs)
      common /qmass_overlap/fmo,fmfac1,fmfac2
*
      call fmtv_wilson(v,x)
      call sqrtoverlap(x,w,ncg)
*
      return
      end
