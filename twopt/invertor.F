      subroutine invert_wilson(rphi,xsol,nitcg,ncgs)
*************************************************************************
*                                                                       *
*     The conjugate gradient routine for solving Ax=b.                  *
*     Finds vectors X and R in fermion force.                           *
*                                                                       *
*************************************************************************
#include "size.h"
      integer, parameter :: nsd=2*nvsize
      integer l, nx, nitrc, nitcg, ncgs
*
      complex*16 rphi(nvsize), xsol(nvsize), rphinorm(nvsize)
      complex*16 dtp, ddtp, p, r, x
      complex*16 mydot
      real*8 rt(nsd)
      real*8 alpha, alphan, alphad, beta, betan, betad, delit, relit,
     1       rnorm
*
      common /param_sqrt/ delit, nitrc
      common /arrays_wil/ r(nvsize), p(nvsize), x(nvsize),
     1                dtp(nvsize), ddtp(nvsize)
      equivalence (rt(1),r(1))
*************************************************************************
*        Find U^+ used in D and D^+                                     *
*************************************************************************
      rnorm=sqrt(real(mydot(rphi,rphi)))
      rphinorm=rphi/rnorm
!$OMP parallel do default(shared)
      do l=1,nvsize
         x(l)=(0.0D0,0.0D0)
         r(l)=rphinorm(l)
         p(l)=rphinorm(l)
      end do
*
      alphan=0.0D0
!$OMP parallel do default(shared) reduction(+: alphan)
      do l=1,nsd
         alphan=alphan+rt(l)*rt(l)
      end do
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
         call fmtv(p,dtp)
         call fmv(dtp,ddtp)
*************************************************************************
*        Update the solution and residual                               *
*************************************************************************
         alphad=real(mydot(p,ddtp))
         alpha=alphan/alphad
!$OMP parallel do default(shared)
         do l=1,nvsize
            x(l)=x(l)+alpha*p(l)
            r(l)=r(l)-alpha*ddtp(l)
         end do
*************************************************************************
*        Test the exit condition                                        *
*************************************************************************
         betan=0.0D0
!$OMP parallel do default(shared) reduction(+: betan)
         do l=1,nsd
            betan=betan+rt(l)*rt(l)
         end do
         if (betan.lt.relit) go to 30
*************************************************************************
*        Update the conjugate vector                                    *
*************************************************************************
         beta=betan/betad
         betad=betan
         alphan=betan
         
!$OMP parallel do default(shared)
         do l=1,nvsize
            p(l)=r(l)+beta*p(l)
         end do
 20      continue
*
 30   continue
*
      r=r*rnorm
      p=p*rnorm
      x=x*rnorm
      dtp=dtp*rnorm
      ddtp=ddtp*rnorm
*
      xsol=x

      ncgs=0
      return
      end subroutine
