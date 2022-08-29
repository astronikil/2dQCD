      subroutine stout(klev)
#include "size.h"
#include "improve.h"
      integer l, ilev, klev
      complex*16 u, us, uin, ud, v,eq, qm, zid
      complex*16 a1, a2, a3, a4
      real*8 tr
      common /config/ u(nc2,mb)
      common /configsmr/ us(nc2,mb)
      common /configsmr_in/ uin(nc2,mb)
      common /expdat/ eq(nc2,mb)
      common /qmat/ qm(nc2,mb)
      common /donfig/ ud(nc2,mb)
      zid=(0.0D0,2.0D0)
*
      us=u
*
      do ilev=1,klev
*
         uin=us
         call rest_in !output in ud
*
         do l=1,mb
            qm(1,l)=-eps4*aimag(ud(1,l))
            qm(4,l)=-eps4*aimag(ud(4,l))
            qm(2,l)=-eps4*(ud(2,l)-conjg(ud(3,l)))/zid
            qm(3,l)=-eps4*(ud(3,l)-conjg(ud(2,l)))/zid
            tr=(qm(1,l)+qm(4,l))/dfloat(nc)
            qm(1,l)=qm(1,l)-cmplx(tr,0.0D0)
            qm(4,l)=qm(4,l)-cmplx(tr,0.0D0)
        enddo
*
        call expsun ! eq= exp(i qm)
*
        do l=1,mb
           us(1,l)=eq(1,l)*uin(1,l)+eq(3,l)*uin(2,l)
           us(2,l)=eq(2,l)*uin(1,l)+eq(4,l)*uin(2,l)
           us(3,l)=eq(1,l)*uin(3,l)+eq(3,l)*uin(4,l)
           us(4,l)=eq(2,l)*uin(3,l)+eq(4,l)*uin(4,l)
        enddo
      enddo
*
      return
      end
