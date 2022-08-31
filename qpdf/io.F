      subroutine readconfig(beta, nflav, fm, dt, nitrc, delit, ntraj, 
     1           nmd, istart, time)
************************************************************************
*                                                                      *
*     Write the configuration file                                     *
*                                                                      *
************************************************************************
#include "size.h"
      complex u, zi, tpiz
      real ua
      common /config/ u(mb)
      common /configa/ ua(mb)
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
      read(10,*) beta, nflav, fm, dt, nitrc, delit, ntraj, nmd,
     2         istart, time,
     3         iseed, (ua(l),l=1,mb)
      do l=1,mb
         u(l)=exp(zi*ua(l))
      enddo
      end subroutine readconfig
      subroutine readrajamani
#include "size.h"
      integer i1, i2, i3, id, l, ll, ino
      real theta, pifac
      real ua
      common /configa/ ua(mb)
      common /newold/ ino(mv)
      do l=1,mb
         read(10,*)i1,i2,i3,id,theta
         ll=ino(i1+mx*(i2-1+my*(i3-1)))
         ua(ll+(id-1)*mv)=theta
      enddo
      end
