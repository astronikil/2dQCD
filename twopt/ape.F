      subroutine ape2d
#include "size.h"
      complex u, uin, ud, uy
      common /configsmr/ u(nc2*mv,3)
      common /config_y/ uy(nc2*mv,3)
      common /configsmr_in/ uin(nc2*mv,3)
      common /donfig/ ud(nc2*mv,3)
      common /apen/ nlevape
      common /apee/ epsape
      epsape2=epsape/2.0
*
      uin=u
      do ilev=1, nlevape
         call staple_2d
         uin(:,1)=(1.0-epsape)*uin(:,1)+epsape2*ud(:,1)
         uin(:,2)=(1.0-epsape)*uin(:,2)+epsape2*ud(:,2)
         uin(:,3)=0.0
      enddo  
*
      uy=uin
      end subroutine ape2d
      subroutine staple_2d
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
      complex*16 utmp(nc2,mb)
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
*
      enddo
      enddo
*
      call hconj(ud,utmp,mb)
      ud=utmp
      return
      end
