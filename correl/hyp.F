      subroutine hyp_gauge
************************************************************************
*                                                                      *
*     HYP smear the links. Smear parameters should be set in           *
*     improve.h                                                        *
*                                                                      *
*     Dependence: A call to plaquette_gauge(ua) before the call to hyp *
*                                                                      *
************************************************************************
#include "improve.h"
#include "size.h" 
      complex u, us, zi
      real ua, ud, y, ys
************************************************************************
*     Inputs                                                           *
************************************************************************
      common /configa/ ua(mb)
      common /plaq_gauge/plg(mb)
************************************************************************
*     outputs                                                          *
************************************************************************
      common /configsmr/ us(mb)
      common /configasmr/ uas(mb)
************************************************************************
*     Workspace                                                        *
************************************************************************
      common /plaq_wrk/ pl(mb)
      common /smear/ ys(mb)
      common /cftdata/ q, QQ
      zi=(0.,1.)
*
      call plaquette_gauge
*
      do l=1,mb
         pl(l)=plg(l)
      enddo
*
      call hypsmear_y
*
      pifac=acos(-1.0)/mx
!$OMP parallel do default(shared)
      do l=1, mb
         uas(l)=(ua(l)-ys(l)) !*q*QQ+pifac
         !us(l)=exp(zi*uas(l))
      enddo
*
      return
      end subroutine hyp_gauge
      subroutine hyp_force(w)
************************************************************************
*                                                                      *
*     HYP smear the force. Smear parameters should be set in           *
*     improve.h                                                        *
*                                                                      *
*     Dependence: A call to plaquette_force(ua) before the call to hyp *
*                                                                      *
************************************************************************
#include "improve.h"
#include "size.h" 
      real w(mb)
************************************************************************
*     Inputs/Output                                                    *
************************************************************************
      common /plaq_force/plf(mb)
      common /plaq_gauge/plg(mb)
************************************************************************
*     Workspace                                                        *
************************************************************************
      common /plaq_wrk/ pl(mb)
      common /smear/ ys(mb)
*
      call plaquette_gauge
      call plaquette_force(w)
*
c      do l=1,mb
c         pl(l)=pl(l)*cos(plg(l))
c      enddo
*
      call hypsmear_y
*
!$OMP parallel do default(shared)
      do l=1, mb
         w(l)=w(l)-ys(l)
      enddo
*
      return
      end subroutine hyp_force
      subroutine hypsmear_y
************************************************************************
*                                                                      *
*     General subroutine to compute the plaquette part of HYP.         *
*     Smears whatever is in plaq_wrk.                                  *
*                                                                      *
************************************************************************
#include "improve.h"
#include "size.h" 
      real ys
      real w(mv)
      integer MVA(3,3), sgn(3,3)
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
