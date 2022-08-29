      subroutine ranu1
*************************************************************************
*                                                                       *
*     Generate random  u(1) matrices                                    *
*                                                                       *
*************************************************************************
#include "size.h"
      complex u, zi, tpiz
      real ua
      real    RUNUM3(mb)
      common /config/ u(mb)
      common /configa/ ua(mb)
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
      common /newold/ ino(mv)
      call dlarnv(2,iseed,mb, RUNUM3)
!$OMP parallel do default(shared)
      do l=1,mb
         ua(l)=pi*RUNUM3(l)
         u(l)=exp(zi*ua(l))
      enddo                  
*
      write(*,'(2x,a)') 'Initialised to random link configuration'
      return
      end
      subroutine ordu1
*************************************************************************
*                                                                       *
*     Generate ordered links                                            *
*                                                                       *
*************************************************************************
#include "size.h"
      complex u, zi
      common /config/ u(mb)
      common /configa/ ua(mb)
      zi=(0.,1.)
      do l=1,mb
         ua(l)=(0.,0.)+acos(-1.)/mx
         u(l)=exp(zi*ua(l))
      enddo
      write(*,'(2x,a)') 'Initialised to ordered link configuration'
      end subroutine ordu1
      subroutine linkc
*************************************************************************
*                                                                       *
*     Fill the complex conjugate of link variables into an array.       *
*                                                                       *
*************************************************************************
#include "size.h"
      complex u, ud
      common /configsmr/ u(mb)
      common /donfig/ ud(mb)
      common /iupidn/ iup(mb), idn(mb)
*
      do 10 mu=1,3
         nu=(mu-1)*mv
         do 10 l=1,mv
            nn=l+nu
            ndo=idn(nn)+nu
 10         ud(nn)=conjg(u(ndo))
*
      return
      end
      subroutine gtransform
#include "size.h"
      integer iseed(4), iup, idn, l, ll
      real g(mv), ua
      common /iupidn/ iup(mb), idn(mb)
      common /configa/ ua(mb)
      iseed=1; iseed(4)=2
      call dlarnv(2,iseed,mv,g)
      do l=1,mb
         ll=mod(l-1,mv)+1
         ua(l)=ua(l)+10.0*(g(iup(l))-g(ll))
      enddo
      end
