      subroutine addrc
***********************************************************************
*                                                                     *
*     Loads the addresses required during the update                  *
*                                                                     *
***********************************************************************
#include "size.h"
      integer nshf, ino, iup, idn, li, ic, i1, i2, i, j, mu, j1, j2, j4
      integer jxup(mx), jxdn(mx),
     1        jtup(mt), jtdn(mt), lpt, ix, iy, it, l
      parameter (nshf=mv/2)
*
      common /newold/ ino(mv)
      common /iupidn/ iup(mv,nd), idn(mv,nd)
      common /polptr/ lpt(mt,msp)
***********************************************************************
*                                                                     *
*     Some pre-computation of modular arithmetic                      *
*                                                                     *
***********************************************************************
      do 10 i=1,mx
         jxup(i)=mod(i,mx)+1
         jxdn(i)=mod(i+mx-2,mx)+1
 10      continue
*
      do 11 i=1,mt
         jtup(i)=mod(i,mt)+1
         jtdn(i)=mod(i+mt-2,mt)+1
 11      continue
*
***********************************************************************
*                                                                     *
*     Computation of forward and backward pointers and chequerboard   *
*                                                                     *
***********************************************************************
      i1=0
      i2=nshf
      li=0
         do j2=1,mt
            do j1=1,mx
               li=li+1
               if (mod(j1+j2,2).eq.0) then
                  i1=i1+1
                  ic=i1
               else
                  i2=i2+1
                  ic=i2
               endif
               ino(li)=ic
               iup(ic,1)=jxup(j1)+mx*(j2-1)
               iup(ic,2)=j1+mx*(jtup(j2)-1)
               idn(ic,1)=jxdn(j1)+mx*(j2-1)
               idn(ic,2)=j1+mx*(jtdn(j2)-1)
            enddo
         enddo
***********************************************************************
*                                                                     *
*     Conversion from linear to chequerboard coordinates.             *
*                                                                     *
***********************************************************************
      do 30 mu=1,nd
         do 31 i=1,mv
            j=iup(i,mu)
            iup(i,mu)=ino(j)
            j=idn(i,mu)
            idn(i,mu)=ino(j)
 31         continue
 30      continue
***********************************************************************
*                                                                     *
*     Pointers for the Polyakov loop along t-direction                *
*                                                                     *
***********************************************************************
      it=1   !compute the addresses on the slice perp to y-axis
      l=0
      do ix=1,mx
         l=l+1
         lpt(it,l)=ino(ix+mx*(it-1))
      enddo
*
      do it=2,mt
         do i=1,msp
            lpt(it,i)=iup(lpt(it-1,i),2)
         enddo
      enddo
*
      mu=mv
      do it=1,mt
         do i=1,msp
            lpt(it,i)=lpt(it,i)+mu
         enddo
      enddo
      return
      end
