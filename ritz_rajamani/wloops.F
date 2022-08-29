      subroutine rest_2
************************************************************************
*                                                                      *
*     Calculate the product of gauge fields which interact with a link *
*                                                                      *
************************************************************************
      include 'global.f'
      complex*16 u
      complex*16 a, b, c
      real*8 ud, ua
      integer iup, idn, lr, l, ii, ll, mu, nu, l1, m
      common /config/ u(mb)
      common /donfiga/ ud(mb)
      common /configa/ ua(mb)
      common /iupidn/ iup(mb), idn(mb)

      integer MVA(3,3)
      data MVA /0,1,2, 1,0,2, 2,0,1/

************************************************************************
*           Products of 2 sets of 3 links completing plaquettes in this*
************************************************************************

      LR  = mb 

      do l = 1, LR
         ii = ((l - 1) / mv) + 1
         ll = mod(l - 1, mv) + 1
         mu = MVA(1, ii) * mv
         ud(l) = (0., 0.)
      do m = 2, 3
         nu = MVA(m, ii) * mv
         l1 = iup(ll+nu) + mu
         b = ua(l1)
         l1 = ll + nu
         a = ua(l1)
         c = a+b 
         a = -c
         l1 = iup(ll+mu) + nu
         b = ua(l1)
         ud(l) = ud(l) + (b+a)
*
         l1 = idn(iup(ll+mu) + nu) + nu
         a = ua(l1)
         l1 = idn(ll+nu) + mu
         b = ua(l1)
         c = b+a 
         a = -c
         l1 = idn(ll+nu) + nu
         b = ua(l1)
         ud(l) = ud(l) + (a+b)
      end do
      end do

      return
      end
