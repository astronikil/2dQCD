      subroutine ransun
*************************************************************************
*                                                                       *
*     Generate random  SU(N) matrices                                   *
*                                                                       *
*************************************************************************
#include "size.h"
      integer k, l, iseed, ino
      real*8  RUNUM(4*ntot), pi, tpi, rmag
      complex*16 u, zi, tpiz, a, b
      common /config/ u(nc,nc,mb)
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
      common /newold/ ino(mv)
*
      call dlarnv(2,iseed,4*ntot, RUNUM)
      k=1
!$OMP parallel do default(shared)
      do l=1,mb
         a=cmplx(RUNUM(k),RUNUM(k+1))
         b=cmplx(RUNUM(k+2),RUNUM(k+3))
         rmag=sqrt(abs(a)*abs(a)+abs(b)*abs(b))
         a=a/rmag
         b=b/rmag
         k=k+4
         u(1,1,l)=a
         u(1,2,l)=b
         u(2,1,l)=-conjg(b)
         u(2,2,l)=conjg(a)
      enddo                  
*
      write(28,'(2x,a)') 'Initialised to random link configuration'
      return
      end
      subroutine ordsun
*************************************************************************
*                                                                       *
*     Generate ordered links                                            *
*                                                                       *
*************************************************************************
#include "size.h"
      integer l, ic
      real*8 rpi
      complex*16 u, zi
      common /config/ u(nc,nc,mb)

      zi=(0.0D0,1.0D0)
      rpi=acos(-1.0D0)
*
      u=(0.0D0,0.0D0)
*
      do l=1,mb
         u(1,1,l)=(1.0D0,0.0D0)*exp(-zi*rpi/mx)
         u(2,2,l)=(1.0D0,0.0D0)*exp(+zi*rpi/mx)
      enddo
*
      write(28,'(2x,a)') 'Initialised to ordered link configuration'
      return
      end
