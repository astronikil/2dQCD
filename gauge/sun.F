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
         do ic=1, nc-1
            u(ic,ic,l)=(1.0D0,0.0D0)*exp(-(rpi*zi/mx))
         enddo
         u(nc,nc,l)=(1.0D0,0.0D0)*exp((nc-1)*(rpi*zi/mx))
      enddo
*
      write(28,'(2x,a)') 'Initialised to ordered link configuration'
      return
      end
