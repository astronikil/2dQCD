      subroutine ransun
*************************************************************************
*                                                                       *
*     Generate random  u(1) matrices                                    *
*                                                                       *
*************************************************************************
#include "size.h"
      integer k, l, iseed, ino
      real*8  RUNUM(4*ntot), pi, tpi, rmag
      complex*16 u, zi, tpiz, a, b
      common /config/ u(nc,nc,mb)
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
      common /newold/ ino(mv)
*
      call dlarnv(2,iseed,4*ntot, RUNUM)
!$OMP parallel do default(shared)
      k=1
      do l=1,mb
         a=cmplx(RUNUM(k),RUNUM(k+1))
         b=cmplx(RUNUM(k+2),RUNUM(k+3))
         rmag=sqrt(abs(a)*abs(a)+abs(b)*abs(b))
         a=a/rmag
         b=b/rmag
         k=k+4
         u(1,1,l)=a
         u(1,2,l)=b
         u(2,1,l)=-conjg(b)
         u(2,2,l)=conjg(a)
      enddo                  
*
c      call unitarize
*
      write(14,'(2x,a)') 'Initialised to random link configuration'
      return
      end
      subroutine ordsun
*************************************************************************
*                                                                       *
*     Generate ordered links                                            *
*                                                                       *
*************************************************************************
#include "size.h"
      integer l, ic
      complex*16 u, zi, zfac
      common /config/ u(nc,nc,mb)
      zi=(0.0,1.0)
      zfac=exp(zi*acos(-1.0)/mx)
*
      u=(0.0D0,0.0D0)
      do l=1,2*mv
         u(1,1,l)=zfac !(1.0D0,0.0D0)
         u(2,2,l)=conjg(zfac) !(1.0D0,0.0D0)
         !u(1,1,l)=(1.0D0,0.0D0)
         !u(2,2,l)=(1.0D0,0.0D0)
      enddo
      do l=2*mv+1,mb
         u(1,1,l)=(1.0D0,0.0D0)
         u(2,2,l)=(1.0D0,0.0D0)
      enddo
      write(14,'(2x,a)') 'Initialised to ordered link configuration'
      return
      end
      subroutine linkc
*************************************************************************
*                                                                       *
*     Fill the complex conjugate of link variables into an array.       *
*                                                                       *
*************************************************************************
#include "size.h"
      integer iup, idn, mu, nu, l ,nn, ndo, ic, jc
      complex*16 u, ud
      common /configsmr/ u(nc,nc,mb)
      common /donfig/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
*
      do mu=1,3
         nu=(mu-1)*mv
         do l=1,mv
            nn=l+nu
            ndo=idn(nn)+nu
            do ic=1,nc
               do jc=1,nc
                  ud(ic,jc,nn)=conjg(u(jc,ic,ndo))
               enddo
            enddo
         enddo
      enddo
*
      return
      end
      subroutine rang(g)
*************************************************************************
*                                                                       *
*     Generate random  su(2) matrices                                   *
*                                                                       *
*************************************************************************
#include "size.h"
      integer k, l, iseed
      real*8  RUNUM(4*ntot), pi, tpi, rmag
      complex*16 zi, tpiz, a, b, g(nc,nc,mv)
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
*
      call dlarnv(2,iseed,4*ntot, RUNUM)
!$OMP parallel do default(shared)
      k=1
      do l=1,mv
         a=cmplx(RUNUM(k),RUNUM(k+1))
         b=cmplx(RUNUM(k+2),RUNUM(k+3))
         rmag=sqrt(abs(a)*abs(a)+abs(b)*abs(b))
         a=a/rmag
         b=b/rmag
         k=k+4
         g(1,1,l)=a
         g(1,2,l)=b
         g(2,1,l)=-conjg(b)
         g(2,2,l)=conjg(a)
      enddo                  
*
c      call unitarize
*
      write(14,'(2x,a)') 'Initialised to random link configuration'
      return
      end
      subroutine gtransform
#include "size.h"
      integer iseed(4), iup, idn, l, ll
      complex u, g(nc,nc,mv), gc(nc,nc,mv), y(nc,nc)
      common /iupidn/ iup(mb), idn(mb)
      common /config/ u(nc,nc,mb)
*      
      call rang(g)
      call hconj(g,gc,mv)
*
      do l=1,mb
         ll=mod(l-1,mv)+1
         y(:,:)=matmul(g(:,:,ll),u(:,:,l))
         u(:,:,l)=matmul(y(:,:),gc(:,:,iup(l)))
      enddo
*      
      end
      subroutine hconj(u,v,nsize)
#include "size.h"
      complex u(nc,nc,nsize), v(nc,nc,nsize)
      do l=1,nsize
         do ic=1,nc
            do jc=1,nc
               v(ic,jc,l)=conjg(u(jc,ic,l))
            enddo
         enddo
      enddo
      return
      end
