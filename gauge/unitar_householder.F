      program main
#include "size.h"
      integer iseed
      complex*16 u
      common /random/ iseed(4)
      common /config/ u(nc,nc,mb)
      iseed=0; iseed(4)=1
      call addrc
      call ransun

      call unitar
      call unitar
      end
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
         
         rmag=0.0D0
         do ic=1, nc
            rmag = rmag + abs(u(1,ic,l))**2
         enddo


         !write(*,*)rmag;call flush

         if(abs(rmag-1.0).ge.tolu)then
            icount=icount+1
         endif

         v = u(:,:,l)
         call gauss(v, nc, det)
         write(*,*)'##',det

         v = u(:,:,l)

         call lapack_unitar(v)

         u(:,:,l) = v

         call gauss(v, nc, det)

         !write(*,*)det

         u(nc, :, l) = u(nc, :, l)/det
         !u(:,:,l) = u(:,:,l)/det

         !v = u(:,:,l)
         !call gauss(v, nc, det)
         !write(*,*)'##',det


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
      subroutine lapack_unitar(u)
************************************************************************************
*     w(i) = U^dag_{i,j) H_{j,k} U_{k,i}                                           *
************************************************************************************
#include "size.h"
*
      integer LDA,LWORK,LRW
      parameter (LDA=nc, LWORK=2*nc-1, LRW=3*nc-2)
      integer INFO
      complex*16 h(nc,nc), u(nc,nc)
      complex*16 WORK(LWORK), TAU(nc)
*
      call zgeqrf(nc, nc, u, LDA, TAU, WORK, LWORK, INFO)
*
      call zungqr(nc, nc, nc, u, LDA, TAU, WORK, LWORK, INFO) 
*
      return
      end
