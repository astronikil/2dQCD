      subroutine metro(ene0,enef,nacc)
#include "size.h"
************************************************************************
*                                                                      *
*     Metropolis accept/reject step                                    *
*                                                                      *
************************************************************************
      integer nacc, iseed
      real*8 r(1), ene0, enef, fac
      common /random/ iseed(4)
*
      if(enef.le.ene0)then
         nacc=1
      else
         call dlarnv(1,iseed,1,r)
         fac=exp(ene0-enef)
         if(r(1).lt.fac)then
           nacc=1
         else
           nacc=0
         endif 
      endif
      end subroutine metro
