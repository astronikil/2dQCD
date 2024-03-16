      subroutine invert_overlap(rphi,xsol,nitcg,ncgs)
*************************************************************************
*                                                                       *
*     The conjugate gradient routine for solving Ax=b.                  *
*     Finds vectors X and R in fermion force.                           *
*                                                                       *
*************************************************************************
#include "size.h"
*
      integer iup, idn, nitrc, nitcg, l, ncgs, nx
      complex*16 rphi(mvs), xsol(mvs), rphinorm(mvs)
      complex*16 dtp, ddtp, p, r, x
      complex*16 mydot, mychiraldot
      real*8 delit, relit, alphan, alphad, betan, betad, beta, alpha
      real*8 rnorm, rchinorm, rcheck
*
      common /param_overlap/ delit, nitrc
      common /arrays_overlap/ r(mvs), p(mvs), x(mvs),
     1                dtp(mvs), ddtp(mvs)
      common /iupidn/ iup(mb), idn(mb)
*************************************************************************
*        Find U^+ used in D and D^+                                     *
*************************************************************************
      rnorm=sqrt(real(mydot(rphi,rphi)))
      rchinorm=real(mychiraldot(rphi,rphi))/rnorm/rnorm
      rphinorm=rphi/rnorm
*
      x=(0.0D0,0.0D0)
      r=rphinorm
      p=rphinorm
*
      alphan=real(mydot(r,r))
      betad=alphan
      relit=delit
      nitcg=0
*************************************************************************
*                                                                       *
*     The main conjugate gradient loop.                                 *
*                                                                       *
*************************************************************************
      do 20 nx=1,nitrc
         nitcg=nitcg+1
*************************************************************************
*        Multiplying by the matrix: M^+ M                               *
*************************************************************************
         if(rchinorm.lt.relit)then
c            write(*,*)'switch to overlap_square'
            call H_overlap_square(p,ddtp,ncgs)
         else
c            write(*,*)'switch to exact overlap2'
            call H_overlap(p,dtp,ncgs)
            call H_overlap(dtp,ddtp,ncgs)
         endif
*************************************************************************
*        Update the solution and residual                               *
*************************************************************************
         alphad=real(mydot(p,ddtp))
         alpha=alphan/alphad
         x=x+alpha*p
         r=r-alpha*ddtp
*************************************************************************
*        Test the exit condition                                        *
*************************************************************************
         betan=real(mydot(r,r))
         if (betan.lt.relit) go to 30
*************************************************************************
*        Update the conjugate vector                                    *
*************************************************************************
         beta=betan/betad
         betad=betan
         alphan=betan
         p=r+beta*p
 20      continue
*
 30   continue
*
      r=r*rnorm
      p=p*rnorm
      x=x*rnorm
      dtp=dtp*rnorm
      ddtp=ddtp*rnorm

      xsol=x

c      rcheck=real(mydot(rphi,xsol))
c      if(rcheck.lt.0.0)then
c         write(*,*)'!!!!!!!!!!!-ve norm alert!!!!!!!!'
c         write(*,*)rcheck,nitcg,ncgs
c         stop
c      endif
*
      return
      end
