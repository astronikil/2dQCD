      subroutine ransun
*************************************************************************
*                                                                       *
*     Generate random  SU(N) matrices                                   *
*                                                                       *
*************************************************************************
#include "size.h"
      integer, parameter :: K2MAX=2*nc2*mb 
      integer j, l, iseed, ic
      real*8 RUNUM(K2MAX)
      complex*16 u, zi
      common /config/ u(nc2,mb)
      common /random/ iseed(4)
*
      zi=(0.0D0, 1.0D0)
*      
      call dlarnv(2,iseed,K2MAX,RUNUM)
*
      j=1
      do l=1,mb
         do ic=1,nc2
            u(ic,l)=RUNUM(j)+zi*RUNUM(j+1)
            j=j+2
         enddo
      enddo
*
      call unitar
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
      integer l, ic, iseed
      real*8 rpi, r(mb)
      complex*16 u, zi
      common /config/ u(nc,nc,mb)
      common /random/ iseed(4)

      call dlarnv(1,iseed,mb,r)

      zi=(0.0D0,1.0D0)
      rpi=acos(-1.0D0)
*
      u=(0.0D0,0.0D0)
*
c      do l=1,mb
c         do ic=1, nc-1
c            u(ic,ic,l)=(1.0D0,0.0D0)*exp(-(rpi*zi/mx))
c         enddo
c         u(nc,nc,l)=(1.0D0,0.0D0)*exp((nc-1)*(rpi*zi/mx))
c      enddo

      do l=1,mb
c         do ic=1, nc
            !u(ic,ic,l)=(1.0D0,0.0D0)*exp((2.0*rpi*zi/mx)/nc)

c            if(r(l).lt.0.5)then
c               u(1,1,l)=(1.0D0,0.0D0) !*exp(rpi*zi/mx)
c               u(2,2,l)=(1.0D0,0.0D0) !*exp(-rpi*zi/mx)
c            else
               u(1,1,l)=(1.0D0,0.0D0)*exp(2.0D0*rpi*zi/mx/(ncr-1.0D0))
               u(2,2,l)=(1.0D0,0.0D0)*exp(-2.0D0*rpi*zi/mx/(ncr-1.0D0))
c            endif
c         enddo
      enddo
*
      write(28,'(2x,a)') 'Initialised to ordered link configuration'
      return
      end
      subroutine ransun_pert
*************************************************************************
*                                                                       *
*     Generate random  SU(N) matrices                                   *
*                                                                       *
*************************************************************************
#include "size.h"
      integer, parameter :: K2MAX=2*nc2*mb 
      integer j, l, iseed, ic, jc
      real*8 RUNUM(K2MAX), scl
      complex*16 u, v, zi
      common /config/ u(nc,nc,mb)
      common /donfig/ v(nc,nc,mb)
      common /random/ iseed(4)

      scl=0.1
*
      zi=(0.0D0, 1.0D0)

      v=(0.0D0,0.0D0)
*
      do l=1,mb
         v(1,1,l)=(1.0D0,0.0D0)
         v(2,2,l)=(1.0D0,0.0D0)
      enddo
*      
      call dlarnv(2,iseed,K2MAX,RUNUM)
*
      j=1
      do l=1,mb
         do ic=1,nc
            do jc=1,nc           
               v(ic,jc,l)=v(ic,jc,l) + scl*(RUNUM(j)+zi*RUNUM(j+1))
               j=j+2
            enddo
         enddo
      enddo

      do l=1,mb
         u(:,:,l)=matmul(u(:,:,l),v(:,:,l))
      enddo
*
      call unitar
*
      write(28,'(2x,a)') 'Initialised to random link configuration'
      return
      end
