      subroutine gnoise
************************************************************************
*                                                                      *
*     Initialisation of Gaussian random momenta.                       *
*                                                                      *
************************************************************************
#include "size.h"
      integer, parameter :: K2MAX=2*nca*mb 
      real*8, parameter :: fac=1.0D0/sqrt(dfloat(2.0D0))
*
      integer iseed, ino, ica, j, l
      real*8 RUNUM(K2MAX)
      real*8 pi ,tpi, a1
      complex*16 qf, fo, tgen
      complex*16 tpiz, zi
*
      common /consts/ pi, tpi, tpiz, zi
      common /forces/ fo(nc2,mb), qf(nc2,mb)
      common /gellmann/ tgen(nc2, nca)
      common /random/ iseed(4)
      common /newold/ ino(mv)
*
      call dlarnv(1,iseed,K2MAX, RUNUM)
************************************************************************
*     Langevin step: set up the gaussian noise (momenta).              *
************************************************************************
      qf=(0.0D0,0.0D0)
      j=1
      do l=1,mb
         do ica=1,nca
            a1=sqrt(-2.0D0*log(RUNUM(j)))*cos(tpi*RUNUM(j+1))
            qf(:,l)=qf(:,l)+a1*tgen(:,ica)
            j=j+2
         enddo
      enddo
*
      !correct for normalization Tr(g_1 g_2)=2 delta_12 for generators
      qf=qf*fac 
*
      return
      end
