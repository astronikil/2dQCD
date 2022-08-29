      subroutine rest_in
************************************************************************
*                                                                      *
*     Calculate the product of gauge fields which interact with a link *
*                                                                      *
************************************************************************
#include "size.h"
      integer l, ll, l1, l2, ii, mu, nu, m, iup, idn
      integer ic, jc
      integer LR
      complex*16 u, ud
      complex*16 v(nc,nc), w(nc,nc)
      common /config/ u(nc,nc,mb)
      common /donfig/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
*
      integer MVA(2,2)
      data MVA /0,1, 1,0/
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
         ud(:,:,l) = (0.0D0, 0.0D0)
      do m = 2, nd
         nu = MVA(m, ii) * mv
*
         l1 = ll + nu
         l2 = iup(ll+nu) + mu

         v = matmul(u(:,:,l1),u(:,:,l2))

         do ic=1,nc
            do jc=1,nc
               w(jc,ic)=conjg(v(ic,jc))
            enddo
         enddo
*
         l1 = iup(ll+mu) + nu

         ud(:,:,l) = ud(:,:,l)+matmul(u(:,:,l1), w(:,:))
*
         l1 = idn(ll+nu) + mu
         l2 = idn(iup(ll+mu) + nu) + nu

         v = matmul(u(:,:,l1), u(:,:,l2))

         do ic=1,nc
            do jc=1,nc
               w(jc,ic)=conjg(v(ic,jc))
            enddo
         enddo
*
         l1 = idn(ll+nu) + nu

         ud(:,:,l) = ud(:,:,l)+matmul(w(:,:),u(:,:,l1))
*
      end do
*
      v = matmul(u(:,:,l), ud(:,:,l))
      ud(:,:,l) = v
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
      integer ic, jc
      integer LR
      complex*16 u, ud
      complex*16 v(nc,nc), w(nc,nc)
      common /config/ u(nc,nc,mb)
      common /donfig/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
*
      integer MVA(2,2)
      data MVA /0,1, 1,0/
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
         ud(:,:,l) = (0.0D0, 0.0D0)
      do m = 2, nd
         nu = MVA(m, ii) * mv
*
         l1 = ll + nu
         l2 = iup(ll+nu) + mu

         v = matmul(u(:,:,l1),u(:,:,l2))

         do ic=1,nc
            do jc=1,nc
               w(jc,ic)=conjg(v(ic,jc))
            enddo
         enddo
*
         l1 = iup(ll+mu) + nu

         ud(:,:,l) = ud(:,:,l)+matmul(u(:,:,l1), w(:,:))
*
         l1 = idn(ll+nu) + mu
         l2 = idn(iup(ll+mu) + nu) + nu

         v = matmul(u(:,:,l1), u(:,:,l2))

         do ic=1,nc
            do jc=1,nc
               w(jc,ic)=conjg(v(ic,jc))
            enddo
         enddo
*
         l1 = idn(ll+nu) + nu

         ud(:,:,l) = ud(:,:,l)+matmul(w(:,:),u(:,:,l1))
*
      end do
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
      integer i, it, l, lpt, ic, jc
      complex*16 u, v, pol
      complex*16 w(nc,nc)
      common /config/ u(nc,nc,mb)
      common /arrays/ v(nc,nc,mb)
      common /polptr/ lpt(mt,msp)
***********************************************************************
*                                                                     *
*     Initialise the matrices                                         *
*                                                                     *
***********************************************************************
      do i=1,msp
         l=lpt(1,i)
         v(:,:,i)=u(:,:,l)
      enddo
***********************************************************************
*                                                                     *
*     Construct the Polyakov loops                                    *
*                                                                     *
***********************************************************************
      do it=2,mt
         do i=1,msp
            l=lpt(it,i)
            w=matmul(v(:,:,i), u(:,:,l))
            v(:,:,i)=w
         enddo
      enddo
***********************************************************************
*                                                                     *
*     ... and finally, the average trace.                             *
*                                                                     *
***********************************************************************
      pol=(0.0D0,0.0D0)
      do i=1,msp
         do ic=1,nc
            pol=pol+v(ic,ic,i)+v(ic,ic,i)
         enddo
      enddo
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
      integer mu, nu, muu, nuu, l, l1, l2, l0, ic, jc, nplaq
      real*8 eser, eter, es
      complex*16 y, yc, xes, xet
      complex*16 plqs,plqt,xeser,xeter
      complex*16 g(nc,nc),h(nc,nc)
      common /config/ y(nc,nc,mb)
      common /arrays/ yc(nc,nc,mb)
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
      nplaq=0
      do muu=1, nd
         do nuu=muu+1, nd
            mu=(muu-1)*mv
            nu=(nuu-1)*mv
            do l=1,mv
            nplaq=nplaq+1
               l1=l+mu
               l2=iup(l+mu)+nu
               g=matmul(y(:,:,l1),y(:,:,l2))
               l2=iup(l+nu)+mu
               h=matmul(g(:,:),yc(:,:,l2))
               l0=l+nu
               g=matmul(h,yc(:,:,l0))
               do ic=1,nc
                  xet = xet + g(ic,ic)
               enddo
            enddo
         enddo
      enddo
*
      write(*,*)nplaq,mv
      es=dfloat(nc*nplaq)-real(xet) !/nplane/mv/nc
      !es=-real(xet) !/nplane/mv/nc
      return
      end
