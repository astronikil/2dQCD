      subroutine hyp(eps1,eps2,nlev)
************************************************************************
*                                                                      *
*     HYP smear the links. Smear parameters should be set in           *
*     improve.h                                                        *
*                                                                      *
*     Dependence: A call to plaquette_gauge(ua) before the call to hyp *
*                                                                      *
************************************************************************
      include 'global.f'
      complex*16 u, zi
      real*8 ua, uas, ud, y, ys, pl, plg, eps1, eps2
      integer l, nlev, ilev, nq
************************************************************************
*     Inputs                                                           *
************************************************************************
      common /configa/ ua(mb)
      common /plaq_gauge/plg(mb)
************************************************************************
*     outputs                                                          *
************************************************************************
      common /config/ u(mb)
      common /configasmr/ uas(mb)
************************************************************************
*     Workspace                                                        *
************************************************************************
      common /plaq_wrk/ pl(mb)
      common /smear/ ys(mb)
      zi=(0.,1.)
*
      do ilev=1,nlev
         call plaquette_gauge
*
         do l=1,mb
            pl(l)=plg(l)
         enddo
*
         call hypsmear_y(eps1,eps2)

         do l=1, mb
            ua(l)=ua(l)-ys(l) !this step is compatible only with csw=0
         enddo
      enddo
*
*
      return
      end subroutine hyp
      subroutine hypsmear_y(eps1,eps2)
************************************************************************
*                                                                      *
*     General subroutine to compute the plaquette part of HYP.         *
*     Smears whatever is in plaq_wrk.                                  *
*                                                                      *
************************************************************************
      include 'global.f'
      real*8 ys, pl, sg, eps1, eps2, ep4, epsq8, w(mv)
      integer MVA(3,3), sgn(3,3), iup, idn
      integer l, mu, nu, nrho
      common /iupidn/ iup(mv,nd), idn(mv,nd)
************************************************************************
*     Inputs                                                           *
************************************************************************
      common /plaq_wrk/ pl(mv,nd)
************************************************************************
*     Output                                                           *
************************************************************************
      common /smear/ ys(mv,nd)
************************************************************************
*     Plaquettes are stored as P_rho=eps_{rho mu nu} F_{mu nu}         *
*     MVA(mu,nrho) gives mu perpendicular to  rho in an order          *
*     consistent with evaluation of P_rho.                             *
*                                                                      *
*     sgn(nu,mu) gives the sign to recover F_munu from P_rho           *
************************************************************************
      data MVA /1,2,3, 2,3,1, 3,1,2/
      data sgn /0,1,-1, -1,0,1, 1,-1,0/
*
      ep4=eps1/4.
      epsq8=eps2*eps2/8.
*
      do mu=1,nd
         nu = MVA(2, mu)
         nrho = MVA(3, mu)

         sg=sgn(nu,mu)

************************************************************************
*     Del_nu P_rho                                                     *
************************************************************************
         do l=1,mv
            w(l)=pl(l,nrho)-pl(idn(l,nu),nrho)
         enddo
*
         do l=1,mv
            ys(l,mu)=sg*(ep4*w(l)+
     1               epsq8*(w(iup(l,nrho))+w(idn(l,nrho))))
         enddo

         sg=sgn(nrho,mu)
         
************************************************************************
*     Del_rho P_nu                                                     *
************************************************************************
         do l=1,mv
            w(l)=pl(l,nu)-pl(idn(l,nrho),nu)
         enddo

         do l=1,mv
            ys(l,mu)=ys(l,mu)+sg*(ep4*w(l)+
     1               epsq8*(w(iup(l,nu))+w(idn(l,nu))))
         enddo
      enddo
* 
      return
      end subroutine hypsmear_y
      subroutine plaquette(v)
************************************************************************
*                                                                      *
*     Computes all the plaquettes formed out of v                      *
*                                                                      *
************************************************************************
      include 'global.f'
      real*8 v(mb), a, b, c, pl
      integer MVA(3,3)
      integer l, ll, ii, mu, nu, l1, iup, idn
      common /plaq_wrk/ pl(mb)
      common /iupidn/ iup(mb), idn(mb)
      data MVA /0,1,2, 1,2,0, 2,0,1/
*
!$OMP parallel do default(shared) private(ll,ii,mu,nu,l1,a,b,c)
      do l = 1, mb
         ll = mod(l - 1, mv) + 1
         ii = ((l - 1) / mv) + 1
         mu = MVA(2, ii) * mv
         nu = MVA(3, ii) * mv
         l1 = iup(ll+nu) + mu
         b = v(l1)
         l1 = ll + nu
         a = v(l1)
         c = a+b 
         a = -c
         l1 = iup(ll+mu) + nu
         b = v(l1)
         pl(l) = b+a+v(ll+mu)
      end do
*
      return
      end subroutine plaquette
      subroutine plaquette_gauge
      include 'global.f'
      real*8 ua, pl, plg
      common /configa/ ua(mb)
      common /plaq_gauge/plg(mb)
      common /plaq_wrk/ pl(mb)
*
      call plaquette(ua)
      plg=pl
      return
      end subroutine plaquette_gauge
