      subroutine init_boost_smear
#include "size.h"
      complex u, uc
      common /config_y/ u(nc2,mb)
      common /config_y_conj/ uc(nc2,mb)
c      call ape2d
      call stout2d
      call hconj(u,uc,mb)
      call plaqq2d(es); write(14,*)'2d plaquette: ',es
      return
      end
      subroutine boost_smear(v,w,sgn)
#include "size.h"
      real fcn, fcd
      complex px1(nc), px2(nc), px3(nc), px4(nc), px5(nc), px6(nc)
      Complex u, uc, v(2,nc,mv), x(2,nc,mv), w(2,nc,mv), zi, mydot
      character*1 sgn
      common /config_y/ u(nc,nc,mb)
      common /config_y_conj/ uc(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
      common /boostn/ nfsmear
      common /boostd/ delta
      common /boostmom/ rqm(2)
*************************************************************************
*                                                                       *
*     The matrix multiplication                                         *
*                                                                       *
*************************************************************************
      zi=(0.0,1.0)
*
      if(sgn.eq."+")then
         isgn=1
      else
         isgn=-1
      endif
*         
      fcd=1.0/(1.0+4.0*delta)
      fcn=fcd*delta
*
      x=v
!$OMP parallel do default(shared)
!$OMP+      private(nnu,px1,px2,px3,px4,px5,px6)
      do ifsmear=1,nfsmear
*
         w=0.0
*
         do l=1,mv

            do ic=1,nc
               w(1,ic,l)=fcd*x(1,ic,l)
               w(2,ic,l)=fcd*x(2,ic,l)
            enddo
*
            do nuu=1,2,1
               nu=(nuu-1)*mv
*
               nnu=l+nu
*
               inu=iup(nnu)
               do ic=1,nc
                  px1(ic)=x(1,ic,inu)
                  px2(ic)=x(2,ic,inu)
               enddo
               inu=idn(nnu)
               lnu=inu+nu
               do ic=1,nc
                  px4(ic)=x(1,ic,inu)
                  px5(ic)=x(2,ic,inu)
               enddo
*
               do ic=1,nc
                  do jc=1,nc
                     w(1,ic,l)=w(1,ic,l)+fcn*(
     1                 u(ic,jc,nnu)*px1(jc)*exp(zi*isgn*rqm(nuu))+
     1                 uc(ic,jc,lnu)*px4(jc)*exp(-zi*isgn*rqm(nuu))
     1                   )
                     w(2,ic,l)=w(2,ic,l)+fcn*(
     1                 u(ic,jc,nnu)*px2(jc)*exp(zi*isgn*rqm(nuu))+
     1                 uc(ic,jc,lnu)*px5(jc)*exp(-zi*isgn*rqm(nuu))
     1                   )
                  enddo
               enddo
*
            enddo
         end do
         x=w
      enddo
*
      w=x
      return
      end subroutine
