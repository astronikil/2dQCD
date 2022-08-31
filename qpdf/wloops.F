      subroutine plaquette(v)
************************************************************************
*                                                                      *
*     Computes all the plaquettes formed out of v                      *
*                                                                      *
************************************************************************
#include "size.h"
      real v(mb), a, b, c
      integer MVA(3,3)
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
#include "size.h"
      common /configa/ ua(mb)
      common /plaq_gauge/plg(mb)
      common /plaq_wrk/ pl(mb)
*
      call plaquette(ua)
      plg=pl
      return
      end subroutine plaquette_gauge
      subroutine plaquette_force(w)
#include "size.h"
      real w(mb)
      common /plaq_force/plf(mb)
      common /plaq_wrk/ pl(mb)
*
      call plaquette(w)
      plf=pl
      return
      end subroutine plaquette_force
      subroutine plaq(es)
***********************************************************************
*                                                                     *
*     Construct the average plaquette values.                         *
*                                                                     *
***********************************************************************
#include "size.h"
      common /plaq_gauge/plg(mb)
***********************************************************************
*                                                                     *
*     The temporal plaquette value.                                   *
*                                                                     *
***********************************************************************
      call plaquette_gauge
      es=0.
      do l=1,mb
         es=es+0.5*plg(l)**2
      enddo
*
      return
      end
      subroutine plaqcompact(es)
***********************************************************************
*                                                                     *
*     Construct the average plaquette values.                         *
*                                                                     *
***********************************************************************
#include "size.h"
      common /plaq_gauge/plg(mb)
***********************************************************************
*                                                                     *
*     The temporal plaquette value.                                   *
*                                                                     *
***********************************************************************
      call plaquette_gauge
      es=0.
      do l=1,mb
         es=es+1.-cos(plg(l))
      enddo
*
      return
      end
