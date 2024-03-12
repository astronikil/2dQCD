      subroutine smearforce
#include "size.h"
#include "improve.h"
*
      integer ifl, l, iup, idn, nitcg, klev
      real*8 es
      complex*16 fo, qf
      complex*16 x, ru, rud, p, r, g
      complex*16 rphi, eq, sig
      complex*16 u, ud, v, us, zi
      complex*16 a1, a2, a3, a4
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /forces/ fo(nc2,mb), qf(nc2,mb)
      common /expdat/ eq(nc2,mb)
      common /sigmat/ sig(nc2,mb)
      common /config/ u(nc2,mb)
      common /configsmr/ us(nc2,mb)
*
      zi=(0.0D0,1.0D0)
*************************************************************************
*        Prerequisites for force calculation                            *
*************************************************************************
      call smear(nlev)
c      call plaq_smear(es)
c      write(*,*)'smeared plaq= ',es

      do l=1,mb
         sig(1,l)=conjg(us(1,l))*fo(1,l)+conjg(us(2,l))*fo(2,l)
         sig(2,l)=conjg(us(3,l))*fo(1,l)+conjg(us(4,l))*fo(2,l)
         sig(3,l)=conjg(us(1,l))*fo(3,l)+conjg(us(2,l))*fo(4,l)
         sig(4,l)=conjg(us(3,l))*fo(3,l)+conjg(us(4,l))*fo(4,l)
      enddo

      sig=2.0D0*sig/zi ! convert to notation of sigma in 3d SU2 code
*      
      do klev=nlev-1,0,-1
         do l=1,mb
            fo(1,l)=sig(1,l)*eq(1,l)+sig(3,l)*eq(2,l)
            fo(2,l)=sig(2,l)*eq(1,l)+sig(4,l)*eq(2,l)
            fo(3,l)=sig(1,l)*eq(3,l)+sig(3,l)*eq(4,l)
            fo(4,l)=sig(2,l)*eq(3,l)+sig(4,l)*eq(4,l)
         enddo
*
         call calclambda
         call sumoverplaq
*
         do l=1,mb
            sig(1,l)=fo(1,l)+zi*eps4*sig(1,l)
            sig(2,l)=fo(2,l)+zi*eps4*sig(2,l)
            sig(3,l)=fo(3,l)+zi*eps4*sig(3,l)
            sig(4,l)=fo(4,l)+zi*eps4*sig(4,l)
         enddo
*
         call smear(klev)
      enddo
*
!$OMP parallel do default(shared)
      do l=1,mb
         fo(1,l)=(u(1,l)*sig(1,l)+u(3,l)*sig(2,l))
         fo(2,l)=(u(2,l)*sig(1,l)+u(4,l)*sig(2,l))
         fo(3,l)=(u(1,l)*sig(3,l)+u(3,l)*sig(4,l))
         fo(4,l)=(u(2,l)*sig(3,l)+u(4,l)*sig(4,l))
      enddo

      fo = zi*fo/2.0D0 ! convert back to sigma to 2d QCD force def
*
      return
      end
      subroutine calclambda
#include "size.h"
      real*8, parameter :: tny=1.0E-12
      integer l, iup, idn
      complex*16 h, sig, qm
      complex*16 u, ud, uin, zi
      complex*16 tr1, tr2
      complex*16 eq(nc2)
      complex*16 a,b
      complex*16 a1, a2, a3, a4
      real*8 tr, ar, br, sfac, cfac, csfac, q
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /qmat/ qm(nc2,mb)
      common /configsmr_in/ u(nc2,mb)
      common /sigmat/ sig(nc2,mb)
      common /lamb/ h(nc2,mb)
*
      zi=(0.0D0,1.0D0)
*
!$OMP parallel do default(shared)
!$OMP^ private(a1,a2,a3,a4)
      do l=1, mb
         a1=u(1,l)*sig(1,l)+u(3,l)*sig(2,l)
         a2=u(2,l)*sig(1,l)+u(4,l)*sig(2,l)
         a3=u(1,l)*sig(3,l)+u(3,l)*sig(4,l)
         a4=u(2,l)*sig(3,l)+u(4,l)*sig(4,l)
         sig(1,l)=a1
         sig(2,l)=a2
         sig(3,l)=a3
         sig(4,l)=a4
      enddo
*
!$OMP parallel do default(shared)
!$OMP^ private(ar,br,q,sfac,cfac,csfac,a,b,tr2,tr)
      do l=1,mb
         ar=qm(1,l)*qm(1,l)+qm(3,l)*qm(2,l)
         br=qm(2,l)*qm(3,l)+qm(4,l)*qm(4,l)
*
         q=sqrt(0.5d0*(ar+br))
         if(q.gt.tny)then
            sfac=sin(q)/q
            cfac=cos(q)
            csfac=(cfac-sfac)/q/q
         else
            sfac=1.0D0
            cfac=cos(q)
            csfac=-1.0D0/3.0D0
         endif
* 
*
         tr1=sig(1,l)+sig(4,l)
*
         a=sig(1,l)*qm(1,l)+sig(3,l)*qm(2,l)
         b=sig(2,l)*qm(3,l)+sig(4,l)*qm(4,l)
         tr2=a+b
*          
         h(1,l)=0.5D0*(-tr1*sfac+zi*tr2*csfac)*qm(1,l)
     1             +zi*sfac*sig(1,l)
         h(2,l)=0.5D0*(-tr1*sfac+zi*tr2*csfac)*qm(2,l)
     1             +zi*sfac*sig(2,l)
         h(3,l)=0.5D0*(-tr1*sfac+zi*tr2*csfac)*qm(3,l)
     1             +zi*sfac*sig(3,l)
         h(4,l)=0.5D0*(-tr1*sfac+zi*tr2*csfac)*qm(4,l)
     1             +zi*sfac*sig(4,l)
*
         tr=(h(1,l)+h(4,l))/2.0D0
*
         h(1,l)=real(h(1,l))-cmplx(tr,0.0D0)
         h(4,l)=real(h(4,l))-cmplx(tr,0.0D0)
         h(3,l)=0.5D0*(h(3,l)+conjg(h(2,l)))
         h(2,l)=conjg(h(3,l))
      enddo
*
      return
      end
      subroutine sumoverplaq
#include "size.h"
      real*8, parameter :: tny=1.0E-12
      integer l, ii , ll, l1, l2, mu, nu, m
      integer iup, idn
      complex*16 h, s
      complex*16 u, ud, zi
      complex*16 a1, a2, a3, a4
      complex*16 b1, b2, b3, b4
*
      common /iupidn/ iup(mb), idn(mb)
*
      common /configsmr_in/ u(nc2,mb)
      common /donfig/ ud(nc2,mb)
      common /sigmat/ s(nc2,mb)
      common /lamb/ h(nc2,mb)
*
      integer MVA(2,2)
      data MVA /0,1, 1,0/
*
      call staple_smear
* 
!$OMP parallel do default(shared)
!$OMP^ private(a1,a2,a3,a4)
      do l=1,mb    ! ud =\Lambda. u and h=u^\dag \Lambda
         a1=ud(1,l)*h(1,l)+ud(3,l)*h(2,l)
         a2=ud(2,l)*h(1,l)+ud(4,l)*h(2,l)
         a3=ud(1,l)*h(3,l)+ud(3,l)*h(4,l)
         a4=ud(2,l)*h(3,l)+ud(4,l)*h(4,l)
         ud(1,l)=a1
         ud(2,l)=a2
         ud(3,l)=a3
         ud(4,l)=a4

         a1=h(1,l)*u(1,l)+h(3,l)*u(2,l)
         a2=h(2,l)*u(1,l)+h(4,l)*u(2,l)
         a3=h(1,l)*u(3,l)+h(3,l)*u(4,l)
         a4=h(2,l)*u(3,l)+h(4,l)*u(4,l)
         h(1,l)=a1
         h(2,l)=a2
         h(3,l)=a3
         h(4,l)=a4
      enddo
*
!$OMP parallel do default(shared)
!$OMP^ private(ii,ll,mu,nu,m,l1,l2,a1,a2,a3,a4)
      do l=1,mb
         ii=((l-1)/mv)+1
         ll=mod(l-1,mv)+1
         mu=MVA(1,ii)*mv
         s(1,l)=(0.0D0, 0.0D0)
         s(2,l)=(0.0D0, 0.0D0)
         s(3,l)=(0.0D0, 0.0D0)
         s(4,l)=(0.0D0, 0.0D0)
         do m=2,2
            nu=MVA(m, ii)*mv
*         
*           Term 1
*
            l1=idn(ll+nu)+mu
            l2=iup(idn(ll+nu)+mu)+nu
            a1=u(1,l1)*u(1,l2)+u(3,l1)*u(2,l2)
            a2=u(2,l1)*u(1,l2)+u(4,l1)*u(2,l2)
            a3=u(1,l1)*u(3,l2)+u(3,l1)*u(4,l2)
            a4=u(2,l1)*u(3,l2)+u(4,l1)*u(4,l2)
*
            l1=idn(ll+nu)+nu
*
            s(1,l)=s(1,l)+conjg(a1)*h(1,l1)+conjg(a2)*h(2,l1)
            s(2,l)=s(2,l)+conjg(a3)*h(1,l1)+conjg(a4)*h(2,l1)
            s(3,l)=s(3,l)+conjg(a1)*h(3,l1)+conjg(a2)*h(4,l1)
            s(4,l)=s(4,l)+conjg(a3)*h(3,l1)+conjg(a4)*h(4,l1)
*         
*           Term 2
*
            l1=idn(ll+nu)+mu
            l2=iup(idn(ll+nu)+mu)+nu
            a1=h(1,l1)*u(1,l2)+h(3,l1)*u(2,l2)
            a2=h(2,l1)*u(1,l2)+h(4,l1)*u(2,l2)
            a3=h(1,l1)*u(3,l2)+h(3,l1)*u(4,l2)
            a4=h(2,l1)*u(3,l2)+h(4,l1)*u(4,l2)
*
            l1=idn(ll+nu)+nu
*
            s(1,l)=s(1,l)-conjg(a1)*u(1,l1)-conjg(a2)*u(2,l1)
            s(2,l)=s(2,l)-conjg(a3)*u(1,l1)-conjg(a4)*u(2,l1)
            s(3,l)=s(3,l)-conjg(a1)*u(3,l1)-conjg(a2)*u(4,l1)
            s(4,l)=s(4,l)-conjg(a3)*u(3,l1)-conjg(a4)*u(4,l1)
*         
*           Term 3
*
            l1=ll+nu
            l2=iup(ll+nu)+mu
            a1=h(1,l1)*u(1,l2)+h(3,l1)*u(2,l2)
            a2=h(2,l1)*u(1,l2)+h(4,l1)*u(2,l2)
            a3=h(1,l1)*u(3,l2)+h(3,l1)*u(4,l2)
            a4=h(2,l1)*u(3,l2)+h(4,l1)*u(4,l2)
*
            l1=iup(ll+mu)+nu
*
            s(1,l)=s(1,l)-u(1,l1)*conjg(a1)-u(3,l1)*conjg(a3)
            s(2,l)=s(2,l)-u(2,l1)*conjg(a1)-u(4,l1)*conjg(a3)
            s(3,l)=s(3,l)-u(1,l1)*conjg(a2)-u(3,l1)*conjg(a4)
            s(4,l)=s(4,l)-u(2,l1)*conjg(a2)-u(4,l1)*conjg(a4)
*         
*           Term 4
*
            l1=idn(ll+nu)+mu
            l2=idn(iup(ll+mu)+nu)+nu
            a1=u(1,l1)*h(1,l2)+u(3,l1)*h(2,l2)
            a2=u(2,l1)*h(1,l2)+u(4,l1)*h(2,l2)
            a3=u(1,l1)*h(3,l2)+u(3,l1)*h(4,l2)
            a4=u(2,l1)*h(3,l2)+u(4,l1)*h(4,l2)
*
            l1=idn(ll+nu)+nu
*
            s(1,l)=s(1,l)-conjg(a1)*u(1,l1)-conjg(a2)*u(2,l1)
            s(2,l)=s(2,l)-conjg(a3)*u(1,l1)-conjg(a4)*u(2,l1)
            s(3,l)=s(3,l)-conjg(a1)*u(3,l1)-conjg(a2)*u(4,l1)
            s(4,l)=s(4,l)-conjg(a3)*u(3,l1)-conjg(a4)*u(4,l1)
*         
*           Term 5
*
            l1=ll+nu
            l2=iup(ll+nu)+mu
            a1=u(1,l1)*h(1,l2)+u(3,l1)*h(2,l2)
            a2=u(2,l1)*h(1,l2)+u(4,l1)*h(2,l2)
            a3=u(1,l1)*h(3,l2)+u(3,l1)*h(4,l2)
            a4=u(2,l1)*h(3,l2)+u(4,l1)*h(4,l2)
*
            l1=iup(ll+mu)+nu
            s(1,l)=s(1,l)-u(1,l1)*conjg(a1)-u(3,l1)*conjg(a3)
            s(2,l)=s(2,l)-u(2,l1)*conjg(a1)-u(4,l1)*conjg(a3)
            s(3,l)=s(3,l)-u(1,l1)*conjg(a2)-u(3,l1)*conjg(a4)
            s(4,l)=s(4,l)-u(2,l1)*conjg(a2)-u(4,l1)*conjg(a4)
*         
*           Term 6
*
            l1=ll+nu
            l2=iup(ll+nu)+mu
            a1=u(1,l1)*u(1,l2)+u(3,l1)*u(2,l2)
            a2=u(2,l1)*u(1,l2)+u(4,l1)*u(2,l2)
            a3=u(1,l1)*u(3,l2)+u(3,l1)*u(4,l2)
            a4=u(2,l1)*u(3,l2)+u(4,l1)*u(4,l2)
*
            l1=iup(ll+mu)+nu
*
            s(1,l)=s(1,l)+h(1,l1)*conjg(a1)+h(3,l1)*conjg(a3)
            s(2,l)=s(2,l)+h(2,l1)*conjg(a1)+h(4,l1)*conjg(a3)
            s(3,l)=s(3,l)+h(1,l1)*conjg(a2)+h(3,l1)*conjg(a4)
            s(4,l)=s(4,l)+h(2,l1)*conjg(a2)+h(4,l1)*conjg(a4)
         enddo   
         a1=s(1,l)+ud(1,l)
         a2=s(2,l)+ud(2,l)
         a3=s(3,l)+ud(3,l)
         a4=s(4,l)+ud(4,l)
         s(1,l)=a1
         s(2,l)=a2
         s(3,l)=a3
         s(4,l)=a4
      enddo
*
      return
      end
