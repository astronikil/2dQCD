      subroutine expsun
#include "size.h"
      real*8, parameter :: tny=1.0E-13
      integer l
      integer info, i, j, k
      complex*16 eq, qm, zi, v(nc,nc), u(nc,nc), zw(nc)
      real*8 a, b, q
      real*8 w(nc)
      common /expdat/ qm(nc, nc, mb)
      zi=(0.0D0,1.0D0)
*
      do l=1,mb
         v(:,:) = qm(:,:,l)
         call hermdiag(v,w,u,info)
         do i=1,nc
            zw(i)=exp(zi*w(i))
         enddo
*
         do i=1,nc
            do j=1,nc
               qm(i,j,l)=0.0D0
               do k=1,nc
                  qm(i,j,l)=qm(i,j,l)+zw(k)*u(i,k)*conjg(u(j,k))
               enddo
            enddo
         enddo
      enddo
*
      return
      end
      subroutine hermdiag(h,w,u,info)
************************************************************************************
*     w(i) = U^dag_{i,j) H_{j,k} U_{k,i}                                           *
************************************************************************************
#include "size.h"
*
      integer LDA,LWORK,LRW
      parameter (LDA=nc, LWORK=2*nc-1, LRW=3*nc-2)
      integer INFO
      complex*16 h(nc,nc), u(nc,nc)
      real*8   w(nc), RWORK(LRW)
      complex*16 WORK(LWORK)
*
      u=h
*
      call ZHEEV( 'V', 'U', nc, u,
     1      nc, w, WORK, LWORK, RWORK,
     1       INFO )
*
c      if(INFO.ne.0)then
c        write(*,*)'something went wrong in dsyev', INFO
c      endif
*
      return
      end
