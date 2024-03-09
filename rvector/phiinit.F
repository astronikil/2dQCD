      subroutine phiinit
*************************************************************************
*                                                                       *
*     Construct the noise vector for matrix inversion.                  *
*                                                                       *
*************************************************************************
#include "size.h"
*
      integer ifl, l, ll, ino, ix, iy, it, ncgo, ncgw
      complex*16 r(mvs), rphi, rchk(mvs), rchk2(mvs)
      complex*16 mydot
      complex*16 rin
      common /phi/ rphi(mvs,nf)
      common /rinit/ rin(mvs,nf)
      common /newold/ ino(mv)
*
      call linkc
*
      rin=(0.0D0,0.0D0)
      do ifl=1,nf
         call ranv(r)
         rin(:,ifl)=r(:)
         if(fermtype.eq.'D')then
            call H_overlap(r,rphi(:,ifl),ncgw)
         elseif(fermtype.eq.'M')then
            call H_overlap_pow_p1by4(r,rphi(:,ifl),ncgo, ncgw)
      enddo
*
      call projchiral
      return
      end
      subroutine projchiral
#include "size.h"
      integer l, ifl
      complex*16 rphi
      common /phi/ rphi(ncr,nspr,mv,nf)
*
      do ifl=1,nf
         do l=1,mv
            rphi(:,2,l,ifl)=(0.0D0,0.0D0)
         enddo
      enddo
*
      end
