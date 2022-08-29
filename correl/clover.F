      subroutine clover
************************************************************************
*                                                                      *
*     Computes the compact clover term with thin-links.                *
*                                                                      *
*     Requires a call to plaquette before calling clover.              *
*                                                                      *
************************************************************************
#include "improve.h"
#include "size.h"
      real cvr, cvi
      integer MVA(3,3)
      common /configa/ ua(mb)
      common /plaq_gauge/plg(mb)
      common /clvr/ cvr(mb), cvi(mb)
      common /iupidn/ iup(mb), idn(mb)
      data MVA /0,1,2, 1,2,0, 2,0,1/
*
      call plaquette_gauge
*
!$OMP parallel do default(shared) private(ll,ii,nrho,mu,nu)
      do l=1, mb
         ll = mod(l - 1, mv) + 1
         ii = ((l - 1) / mv) + 1
         nrho = MVA(1,ii)*mv
         mu = MVA(2, ii) * mv
         nu = MVA(3, ii) * mv
         cvr(ll+nrho)=0.25*kcs*(sin(plg(ll+nrho))+
     1                sin(plg(idn(ll+mu)+nrho))+
     1                sin(plg(idn(ll+nu)+nrho))+
     1                sin(plg(idn(idn(ll+nu)+mu)+nrho)))
      enddo
*
!$OMP parallel do default(shared)
      do l=1, mb
         cvi(l)=kcs*cos(plg(l))
      enddo
*
      return
      end subroutine clover
