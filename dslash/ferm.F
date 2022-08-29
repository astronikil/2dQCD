      subroutine ferminit
#include "size.h"
#include "improve.h"
*************************************************************************
*                                                                       *
*     Initializations before fermion computation.                       *
*                                                                       *
*************************************************************************
      call smear(nlev)
      call linkc
c      call clover
*
      return
      end subroutine ferminit
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
      common /donfigsmr/ ud(nc,nc,mb)
      common /iupidn/ iup(mb), idn(mb)
*
      do mu=1,nd
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
