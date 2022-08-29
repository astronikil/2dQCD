      subroutine hyp(eps1,eps2,nlev)
      include 'global.f'
      integer nlev, l, ilev
      complex*16 u, zi
      real*8 ua, ud, y, ys, eps1, eps2
      common /configa/ ua(mb)
      common /donfiga/ ud(mb)
      common /config/ u(mb)
      common /smear/ ys(mb)
      zi=(0.,1.)
      
!$OMP parallel do default(shared)
      do l=1, mb
         ud(l)=ua(l)
      enddo
*
      do ilev=1, nlev
         call hypsmear_y(eps1,eps2)
      enddo
*
!$OMP parallel do default(shared)
      do l=1, mb
         u(l)=exp(zi*ud(l))
      enddo
*
      return
      end subroutine hyp
      subroutine hypsmear_y(eps1,eps2)
      include 'global.f'
      real*8 y, ys, eps1, eps2, ep4, epsq8
      integer l, ll, ii, mu, nu, nrho
      integer MVA(3,3), iup, idn
      common /donfiga/ y(mb)
      common /smear/ ys(mb)
      common /iupidn/ iup(mb), idn(mb)
      data MVA /0,1,2, 1,2,0, 2,0,1/
*
      ep4=eps1/4.
      epsq8=eps2*eps2/8.
*
      do l=1, mb
         ll = mod(l - 1, mv) + 1
         ii = ((l - 1) / mv) + 1
         mu = MVA(1, ii) * mv
         nu = MVA(2, ii) * mv
         nrho = MVA(3, ii) * mv
*
         ys(l)=(1.0-eps1)*y(l)
         ys(l)=ys(l)+ep4*(
     1         y(mu+idn(ll+nu))+y(mu+iup(ll+nu))
     2        -y(nu+iup(ll+mu))+y(nu+iup(idn(ll+nu)+mu))
     3        +y(mu+idn(ll+nrho))+y(mu+iup(ll+nrho))
     4        -y(nrho+iup(ll+mu))+y(nrho+iup(idn(ll+nrho)+mu))
     5        +y(nu+ll)-y(nu+idn(ll+nu))+y(nrho+ll)-y(nrho+idn(ll+nrho))
     6        )
*
         ys(l)=ys(l)+epsq8*(
     1     2*(y(mu+idn(idn(ll+nu)+nrho))+y(mu+idn(iup(ll+nu)+nrho))
     2       +y(mu+idn(iup(ll+nrho)+nu))+y(mu+iup(iup(ll+nu)+nrho)) 
     3       -y(mu+idn(ll+nu))-y(mu+iup(ll+nu))
     4       -y(mu+idn(ll+nrho))-y(mu+iup(ll+nrho)))
     5       -y(nu+iup(idn(ll+nrho)+mu))
     6       +y(nu+iup(idn(idn(ll+nrho)+nu)+mu))
     7       -y(nu+iup(iup(ll+mu)+nrho))
     8       +y(nu+idn(iup(iup(ll+mu)+nrho)+nu))
     9       -y(nrho+iup(idn(ll+nu)+mu))-y(nrho+iup(iup(ll+mu)+nu))
     1       +y(nrho+iup(idn(idn(ll+nu)+nrho)+mu))
     2       +y(nrho+idn(iup(iup(ll+mu)+nu)+nrho))
     3       +y(nu+idn(ll+nrho))-y(nu+idn(idn(ll+nu)+nrho))
     4       +y(nu+iup(ll+nrho))-y(nu+iup(idn(ll+nu)+nrho))
     5       +y(nrho+idn(ll+nu))+y(nrho+iup(ll+nu))
     6       -y(nrho+idn(idn(ll+nu)+nrho))-y(nrho+iup(idn(ll+nrho)+nu))
     7       )
 
      enddo 
* 
      y=ys
*
      return
      end subroutine hypsmear_y
