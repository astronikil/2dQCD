      subroutine gnoise
************************************************************************
*                                                                      *
*     Initialisation of Gaussian random momenta.                       *
*                                                                      *
************************************************************************
#include "size.h"
      integer, parameter :: K2MAX=2*(nc*nc-1)*mb 
      real*8, parameter :: fac=1.0D0/sqrt(dfloat(nc))
*
      integer iseed, ino, j, l
      real*8 RUNUM(K2MAX)
      real*8 pi ,tpi, a1, a2, a3
      complex*16 qf, fo
      complex*16 tpiz, zi
*
      common /consts/ pi, tpi, tpiz, zi
      common /forces/ fo(nc2,mb), qf(nc2,mb)
      common /random/ iseed(4)
      common /newold/ ino(mv)
*
      call dlarnv(1,iseed,K2MAX, RUNUM)
************************************************************************
*     Langevin step: set up the gaussian noise (momenta).              *
************************************************************************
      j=1
      do l=1,mb
         a1=sqrt(-2.0D0*log(RUNUM(j)))*cos(tpi*RUNUM(j+1))
         j=j+2
         a2=sqrt(-2.0D0*log(RUNUM(j)))*cos(tpi*RUNUM(j+1))
         j=j+2
         a3=sqrt(-2.0D0*log(RUNUM(j)))*cos(tpi*RUNUM(j+1))
         j=j+2
*         
         qf(1,l)=a3
         qf(2,l)=a1+zi*a2
         qf(3,l)=a1-zi*a2
         qf(4,l)=-a3
      enddo
*
      qf=qf*fac
*
      return
      end
