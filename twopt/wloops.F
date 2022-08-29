      subroutine rest
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
      common /config/ u(nc2,mb)
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
*
      do l = 1, LR
         ii = ((l - 1) / mv) + 1
         ll = mod(l - 1, mv) + 1
         mu = MVA(1, ii) * mv
         ud(:,l) = (0.0D0, 0.0D0)
      do m = 2, 3
         nu = MVA(m, ii) * mv
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
      end do
*
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
      return
      end
      subroutine rest_in
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
*
      do l = 1, LR
         ii = ((l - 1) / mv) + 1
         ll = mod(l - 1, mv) + 1
         mu = MVA(1, ii) * mv
         ud(:,l) = (0.0D0, 0.0D0)
      do m = 2, 3
         nu = MVA(m, ii) * mv
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
      end do
*
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
      return
      end
      subroutine staple_in
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
*
      do l = 1, LR
         ii = ((l - 1) / mv) + 1
         ll = mod(l - 1, mv) + 1
         mu = MVA(1, ii) * mv
         ud(:,l) = (0.0D0, 0.0D0)
      do m = 2, 3
         nu = MVA(m, ii) * mv
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
      end do
*
      end do
*
      return
      end
      subroutine line(pol)
***********************************************************************
*                                                                     *
*     Construct the average Polyakov loop                             *
*                                                                     *
***********************************************************************
#include "size.h"
      integer i, it, l, lpt
      complex*16 u, v, pol
      complex*16 a1, a2, a3, a4
      common /config/ u(nc2,mb)
      common /arrays/ v(nc2,mb)
      common /polptr/ lpt(mt,msp)
***********************************************************************
*                                                                     *
*     Initialise the matrices                                         *
*                                                                     *
***********************************************************************
      do 10 i=1,msp
         l=lpt(1,i)
         v(:,i)=u(:,l)
 10      continue
***********************************************************************
*                                                                     *
*     Construct the Polyakov loops                                    *
*                                                                     *
***********************************************************************
      do 20 it=2,mt
         do 21 i=1,msp
            l=lpt(it,i)
            a1=v(1,i)*u(1,l)+v(3,i)*u(2,l)
            a2=v(2,i)*u(1,l)+v(4,i)*u(2,l)
            a3=v(1,i)*u(3,l)+v(3,i)*u(4,l)
            a4=v(2,i)*u(3,l)+v(4,i)*u(4,l)
            v(1,i)=a1
            v(2,i)=a2
            v(3,i)=a3
            v(4,i)=a4
 21         continue
 20      continue
***********************************************************************
*                                                                     *
*     ... and finally, the average trace.                             *
*                                                                     *
***********************************************************************
      pol=(0.0D0,0.0D0)
      do 30 i=1,msp
         pol=pol+v(1,i)+v(4,i)
 30      continue
      pol=pol/msp/nc
*
      return
      end
      subroutine plaq(es)
***********************************************************************
*                                                                     *
*     Construct the average plaquette values.                         *
*                                                                     *
***********************************************************************
#include "size.h"
      integer iup, idn
      integer mu, nu, muu, nuu, l, l1, l2, l0, ic, jc
      real*8 eser, eter, es
      complex*16 y, xes, xet
      complex*16 plqs,plqt,xeser,xeter
      complex*16 g(nc,nc),h(nc,nc),yc(nc,nc,mb)
      common /config/ y(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
***********************************************************************
*                                                                     *
*     The temporal plaquette value.                                   *
*                                                                     *
***********************************************************************
      do l=1,mb
         do ic=1,nc
            do jc=1,nc
                yc(ic,jc,l)=conjg(y(jc,ic,l))
            enddo
         enddo
      enddo
*
      xet=(0.0D0,0.0D0)
      plqt=(0.0D0,0.0D0)
      mu=3*mv
      do muu=1, 3
         do nuu=muu+1, 3
            mu=(muu-1)*mv
            nu=(nuu-1)*mv
            do l=1,mv
               l1=l+mu
               l2=iup(l+mu)+nu
               g=matmul(y(:,:,l1),y(:,:,l2))
               l2=iup(l+nu)+mu
               h=matmul(g(:,:),yc(:,:,l2))
               l0=l+nu
               g=matmul(h,yc(:,:,l0))
               xet = xet + g(1,1)+g(2,2)
            enddo
         enddo
      enddo
*
      es=real(xet)
      return
      end
      subroutine plaqq(es)
***********************************************************************
*                                                                     *
*     Construct the average plaquette values.                         *
*                                                                     *
***********************************************************************
#include "size.h"
      integer iup, idn
      integer mu, nu, muu, nuu, l, l1, l2, l0, ic, jc
      real*8 eser, eter, es
      complex*16 y, xes, xet
      complex*16 plqs,plqt,xeser,xeter
      complex*16 g(nc,nc),h(nc,nc),yc(nc,nc,mb)
      common /config/ y(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
***********************************************************************
*                                                                     *
*     The temporal plaquette value.                                   *
*                                                                     *
***********************************************************************
      do l=1,mb
         do ic=1,nc
            do jc=1,nc
                yc(ic,jc,l)=conjg(y(jc,ic,l))
            enddo
         enddo
      enddo
*
      xet=(0.0D0,0.0D0)
      plqt=(0.0D0,0.0D0)
      mu=3*mv
      do muu=1, 3
         do nuu=muu+1, 3
            mu=(muu-1)*mv
            nu=(nuu-1)*mv
            do l=1,mv
               l1=l+mu
               l2=iup(l+mu)+nu
               g=matmul(y(:,:,l1),y(:,:,l2))
               l2=iup(l+nu)+mu
               h=matmul(g(:,:),yc(:,:,l2))
               l0=l+nu
               g=matmul(h,yc(:,:,l0))
               xet = xet + g(1,1)+g(2,2)
            enddo
         enddo
      enddo
*
      es=real(xet)/3./mv/2.
      return
      end
      subroutine plaqq2d(es)
***********************************************************************
*                                                                     *
*     Construct the average plaquette values.                         *
*                                                                     *
***********************************************************************
#include "size.h"
      integer iup, idn
      integer mu, nu, muu, nuu, l, l1, l2, l0, ic, jc
      real*8 eser, eter, es
      complex*16 y, xes, xet
      complex*16 plqs,plqt,xeser,xeter
      complex*16 g(nc,nc),h(nc,nc),yc(nc,nc,mb)
      common /config_y/ y(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
***********************************************************************
*                                                                     *
*     The temporal plaquette value.                                   *
*                                                                     *
***********************************************************************
      do l=1,mb
         do ic=1,nc
            do jc=1,nc
                yc(ic,jc,l)=conjg(y(jc,ic,l))
            enddo
         enddo
      enddo
*
      xet=(0.0D0,0.0D0)
      plqt=(0.0D0,0.0D0)
      muu=1
      nuu=2
      mu=(muu-1)*mv
      nu=(nuu-1)*mv
      do l=1,mv
         l1=l+mu
         l2=iup(l+mu)+nu
         g=matmul(y(:,:,l1),y(:,:,l2))
         l2=iup(l+nu)+mu
         h=matmul(g(:,:),yc(:,:,l2))
         l0=l+nu
         g=matmul(h,yc(:,:,l0))
         xet = xet + g(1,1)+g(2,2)
      enddo
*
      es=real(xet)/mv/2.0
      return
      end
