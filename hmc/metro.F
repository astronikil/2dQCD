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
      subroutine adapttime(ene0,enef,nin,nout)
#include "size.h"
      integer, parameter :: nmin=4
      real*8, parameter :: dmax=0.4, dmin=0.1, dmaxx=0.8
      integer nin, nout
      real*8 ene0, enef, de
*
      if(enef.le.ene0)then
         nout=nin
      else
         de=abs(enef-ene0)
*
         if(de.gt.dmax)nout=nin+1
         if(de.gt.dmaxx)nout=nin+2
         if(de.lt.dmin)nout=nin-1
         if(de.ge.dmin.and.de.le.dmax)nout=nin
*
         if(nout.lt.nmin)nout=nmin
      endif
*
      return
      end
      subroutine adaptstep(ene0,enef,nmd,nmdp)
#include "size.h"
      integer, parameter :: nmin=10, nmax=80
      real*8, parameter :: dmax=0.4, dmin=0.1, dmaxx=0.8
      integer nmd, nmdp
      real*8 ene0, enef, de
*
      if(enef.le.ene0)then
         nmdp=nmd
      else
         de=abs(enef-ene0)
         if(de.gt.dmax)nmdp=nmd+2
         if(de.gt.dmaxx)nmdp=nmdp+4
         if(de.lt.dmin)nmdp=nmd-2
         if(de.ge.dmin.and.de.le.dmax)nmdp=nmd
         if(nmdp.lt.nmin)nmdp=nmin
         if(nmdp.gt.nmax)nmdp=nmax
      endif
*        
      return
      end subroutine adaptstep
