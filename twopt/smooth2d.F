      subroutine stout2d
#include "size.h"
      integer l, ilev, klev
      complex*16 u, us, uin, ud, v,eq, qm, zid
      complex*16 a1, a2, a3, a4
      real*8 tr
      common /config/ u(nc2,mb)
      common /config_y/ us(nc2,mb)
      common /configsmr_in/ uin(nc2,mb)
      common /expdat/ eq(nc2,mb)
      common /qmat/ qm(nc2,mb)
      common /donfig/ ud(nc2,mb)
      common /apen/ nlevsmear
      common /apee/ epssmear
*
      eps4=epssmear/2.0D0
*
      zid=(0.0D0,2.0D0)
*
      us=u
*
      do ilev=1,nlevsmear
*
         uin=us
         call rest2d_in !output in ud
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
      subroutine rest2d_in
************************************************************************
*                                                                      *
*     Calculate the product of gauge fields which interact with a link *
*                                                                      *
************************************************************************
#include "size.h"
      integer l, ll, l1, l2, ii, mu, nu, m, iup, idn
      integer LR
      complex*16 u, ud
      complex*16 a1,a2,a3,a4
      common /configsmr_in/ u(nc2,mb)
      common /donfig/ ud(nc2,mb)
      common /iupidn/ iup(mb), idn(mb)
*
      integer MVA(3,3)
      data MVA /0,1,2, 1,0,2, 2,0,1/
*
************************************************************************
*           Products of 2 sets of 3 links completing plaquettes in this*
************************************************************************
*
      LR  = mb 
      do mu=0, mv, mv
         if(mu.eq.0)then
            nu=mv
         else
            nu=0
         endif
*
      do ll = 1, mv
         l=ll+mu
         ud(:,l) = (0.0D0, 0.0D0)
*
         l1 = ll + nu
         l2 = iup(ll+nu) + mu
         a1=u(1,l1)*u(1,l2)+u(3,l1)*u(2,l2)
         a2=u(2,l1)*u(1,l2)+u(4,l1)*u(2,l2)
         a3=u(1,l1)*u(3,l2)+u(3,l1)*u(4,l2)
         a4=u(2,l1)*u(3,l2)+u(4,l1)*u(4,l2)
*
         l1 = iup(ll+mu) + nu
         ud(1,l)=ud(1,l)+u(1,l1)*conjg(a1)+u(3,l1)*conjg(a3)
         ud(2,l)=ud(2,l)+u(2,l1)*conjg(a1)+u(4,l1)*conjg(a3)
         ud(3,l)=ud(3,l)+u(1,l1)*conjg(a2)+u(3,l1)*conjg(a4)
         ud(4,l)=ud(4,l)+u(2,l1)*conjg(a2)+u(4,l1)*conjg(a4)
*
         l1 = idn(ll+nu) + mu
         l2 = idn(iup(ll+mu) + nu) + nu
         a1=u(1,l1)*u(1,l2)+u(3,l1)*u(2,l2)
         a2=u(2,l1)*u(1,l2)+u(4,l1)*u(2,l2)
         a3=u(1,l1)*u(3,l2)+u(3,l1)*u(4,l2)
         a4=u(2,l1)*u(3,l2)+u(4,l1)*u(4,l2)
*
         l1 = idn(ll+nu) + nu
         ud(1,l)=ud(1,l)+u(1,l1)*conjg(a1)+u(2,l1)*conjg(a2)
         ud(2,l)=ud(2,l)+u(1,l1)*conjg(a3)+u(2,l1)*conjg(a4)
         ud(3,l)=ud(3,l)+u(3,l1)*conjg(a1)+u(4,l1)*conjg(a2)
         ud(4,l)=ud(4,l)+u(3,l1)*conjg(a3)+u(4,l1)*conjg(a4)

         a1=u(1,l)*ud(1,l)+u(3,l)*ud(2,l)
         a2=u(2,l)*ud(1,l)+u(4,l)*ud(2,l)
         a3=u(1,l)*ud(3,l)+u(3,l)*ud(4,l)
         a4=u(2,l)*ud(3,l)+u(4,l)*ud(4,l)
         ud(1,l)=a1
         ud(2,l)=a2
         ud(3,l)=a3
         ud(4,l)=a4
      end do
*
      end do
*
      return
      end
