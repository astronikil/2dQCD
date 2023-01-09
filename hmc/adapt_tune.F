      subroutine adapttime(ene0,enef,nin,nout)
#include "size.h"
      integer, parameter :: nmin=2
      real*8, parameter :: dmax=0.4D0, dmin=0.1D0, dmaxx=0.8D0
      integer nin, nout
      real*8 ene0, enef, de
*
      de=abs(enef-ene0)
*
      if(de.gt.dmax)nout=nin+1
      if(de.gt.dmaxx)nout=nin+1
      if(de.lt.dmin)nout=nin-1
      if(de.ge.dmin.and.de.le.dmax)nout=nin
*
      if(nout.lt.nmin)nout=nmin
*
      return
      end
      subroutine adaptstep(ene0,enef,nmd,nmdp)
#include "size.h"
      integer, parameter :: nmin=2, nmax=10
      real*8, parameter :: dmax=0.4D0, dmin=0.1D0, dmaxx=0.8D0
      integer nmd, nmdp
      real*8 ene0, enef, de
*
      de=abs(enef-ene0)
      if(de.gt.dmax)nmdp=nmd+1
      if(de.lt.dmin)nmdp=nmd-1
      if(de.ge.dmin.and.de.le.dmax)nmdp=nmd
      if(nmdp.lt.nmin)nmdp=nmin
      if(nmdp.gt.nmax)nmdp=nmax
*        
      return
      end subroutine adaptstep
