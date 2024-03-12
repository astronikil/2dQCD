      subroutine expsun_stout
#include "size.h"
      real*8, parameter :: tny=1.0E-13
      integer l
      complex*16 eq, qm, zi
      real*8 a, b, q
      common /qmat/ qm(nc2,mb)
      common /expdat/ eq(nc2,mb)
      zi=(0.0D0,1.0D0)
*
!$OMP parallel do default(shared)
!$OMP^ private(a,b,q)
      do l=1,mb
         a=qm(1,l)*qm(1,l)+qm(3,l)*qm(2,l)
         b=qm(2,l)*qm(3,l)+qm(4,l)*qm(4,l)
         q=sqrt(0.5D0*(a+b))
         if(q.gt.tny)then
            eq(1,l)=zi*sin(q)/q*qm(1,l)
            eq(2,l)=zi*sin(q)/q*qm(2,l)
            eq(3,l)=zi*sin(q)/q*qm(3,l)
            eq(4,l)=zi*sin(q)/q*qm(4,l)
         else
            eq(1,l)=zi*qm(1,l)
            eq(2,l)=zi*qm(2,l)
            eq(3,l)=zi*qm(3,l)
            eq(4,l)=zi*qm(4,l)
         endif
         eq(1,l)=eq(1,l)+cmplx(cos(q),0.0D0)
         eq(4,l)=eq(4,l)+cmplx(cos(q),0.0D0)
      enddo
*
      return
      end
