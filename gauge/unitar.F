      subroutine unitar
#include "size.h"
      real*8, parameter :: tolu=1.0E-8
      integer l, icount, ic, jc
      real*8  rmag, es1, es2
      complex*16 u, a, b, det, v(nc,nc)
      common /config/ u(nc,nc,mb)
*
      call plaq(es1)
      icount=0
      do l=1,mb
*         
         rmag=0.0D0
         do ic=1, nc
            do jc=1, nc
               rmag = rmag + abs(u(ic,jc,l))**2
            enddo
         enddo
*
         if(abs(rmag-1.0D0*nc).ge.tolu)then
            icount=icount+1
         endif
*
         v = u(:,:,l)
*
         call lapack_unitar(v)
*
         u(:,:,l) = v
*
         call gauss(v, nc, det)
*
         u(nc, :, l) = u(nc, :, l)/det
*
      enddo
      call plaq(es2)
*
      write(27,*) ' '
      write(27,*) 'unitarized: frac of defect links: ',icount
      write(27,*) 'es bf unitar: ',es1/dfloat(mv*nc)/1.0D0
      write(27,*) 'es af unitar: ',es2/dfloat(mv*nc)/1.0D0
      write(27,*) ' '
      return
      end
      subroutine lapack_unitar(v)
#include "size.h"
      real*8, parameter :: tny=1.0E-13
      integer l
      integer info, i, j, k
      complex*16 x(nc,nc)
      complex*16 v(nc,nc), u(nc,nc), zw(nc)
      real*8 w(nc)
*
      do i=1,nc
         do j=1,nc
            x(i,j)=conjg(v(j,i))
         enddo
      enddo
*
      x = matmul(x,v)
*
      call hermdiag(x,w,u,info)
*
      do i=1,nc
         zw(i)=1.0D0/sqrt(w(i))
      enddo
*
      do i=1,nc
         do j=1,nc
            x(i,j)=0.0D0
            do k=1,nc
               x(i,j)=x(i,j)+zw(k)*u(i,k)*conjg(u(j,k))
            enddo
         enddo
      enddo
*
      v = matmul(v, x)
      return
      end
c      program main
c#include "size.h"
c      integer iseed
c      complex*16 u
c      common /random/ iseed(4)
c      common /config/ u(nc,nc,mb)
c      iseed=0; iseed(4)=1
c      call addrc
c      call ransun
c
c      u(3,2,3)=1.4*u(3,2,3)
c
c      call unitar
c      call unitar
c      end
