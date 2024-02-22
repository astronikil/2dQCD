      subroutine pauliinit
*************************************************************
*                                                           *
* Pauli matrices in 2d                                      *
* This subroutine is hard-coded for 2d gamma matrices       *
*                                                           *
*************************************************************
#include "size.h"
      integer ic, is
      complex*16 sig, zi, sid(ncr,ncr), s2g(nspr,nspr,4)
      common /paulimat/ sig(nspr*ncr,nspr*ncr,4)
      zi=(0.0,1.0)
*      
      sid=0.0D0
      do ic=1,ncr
         sid(ic,ic)=1.0D0
      enddo
*
*     1
      s2g(1,1,1)=1.0D0; s2g(1,2,1)=0.0D0
      s2g(2,1,1)=0.0D0; s2g(2,2,1)=1.0D0
*
*     Gamma_x
      s2g(1,1,2)=0.0D0; s2g(1,2,2)=1.0D0
      s2g(2,1,2)=1.0D0; s2g(2,2,2)=0.0D0
*     Gamma_t
      s2g(1,1,3)=0.0D0; s2g(1,2,3)=-zi
      s2g(2,1,3)=zi;    s2g(2,2,3)=0.0D0
*
*     Gamma_5
      s2g(1,1,4)=1.0D0; s2g(1,2,4)=0.0D0
      s2g(2,1,4)=0.0D0; s2g(2,2,4)=-1.0D0
*
      do is=1,4
         call dirprod(sid,s2g(:,:,is),sig(:,:,is))
      enddo
      return
      end
      subroutine dirprod(a,b,c)
#include "size.h"
      integer ics, is, ic, jcs, js, jc
      complex*16 a(ncr,ncr), b(nspr,nspr), c(ncs,ncs)
      ics=0
      do is=1,nspr
         do ic=1,ncr
            ics=ics+1
            jcs=0
            do js=1,nspr
               do jc=1,ncr
                  jcs=jcs+1
                  c(jcs,ics)=a(jc,ic)*b(js,is)
               enddo
            enddo
         enddo
      enddo
      return
      end
