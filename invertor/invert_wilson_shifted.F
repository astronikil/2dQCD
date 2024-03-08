      subroutine invert_wilson_shifted(rphi,nitcg)
************************************************************************
*                                                                      *
*     Finds the the inverse of H^2+p_i for multiple p_i.               *
*                                                                      *
************************************************************************
#include "size.h"
#include "zolo.h"
      integer kzol, l, nitrc, nitcg, nx, izol
      complex*16 ps, xs
      complex*16 p, r, x, dtp, ddtp
      complex*16 mydot
      complex*16 rphi(mvs)
      real*8 xinew(nzol), xicurr(nzol), xiold(nzol)
      real*8 betas(nzol), alphas(nzol)
      real*8 rt(mvs2)
      real*8 delit, relit
      real*8 alphan, alphad, betan, betad, beta, alpha, alphold
      real*8 rnorm
      complex*16 rphinorm(mvs)
*
      common /param_sqrt/ delit, nitrc
*
      common /arrays_wil/ r(mvs), p(mvs), x(mvs),
     1                dtp(mvs), ddtp(mvs)
*
      common /wrkforce/ ps(mvs,nzol), xs(mvs,nzol)
*
      equivalence (rt(1),r(1))
*
      rnorm=sqrt(real(mydot(rphi,rphi)))
      rphinorm=rphi/rnorm
*
      kzol=nzol
*************************************************************************
*        Find U^+ used in D and D^+                                     *
*************************************************************************
      do l=1,mvs
         x(l)=(0.0D0,0.0D0) 
         r(l)=rphinorm(l)
         p(l)=rphinorm(l)
         xs(l,1:nzol)=(0.0D0,0.0D0)
         ps(l,1:nzol)=rphinorm(l)
      end do
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
      alphold=1.0D0
      xiold=1.0D0
      xicurr=1.0D0
      betas=0.0D0
      beta=0.0D0
      do 20 nx=1,nitrc
         nitcg=nitcg+1
*
*************************************************************************
*        Multiplying by the matrix: M^+ M                               *
*************************************************************************
         call H_wilson(p,dtp)
         call H_wilson(dtp,ddtp)
*************************************************************************
*        Update the solution and residual                               *
*************************************************************************
         alphad=real(mydot(p,ddtp))
         alpha=alphan/alphad
*
         do izol=1, kzol
            xinew(izol)=xicurr(izol)*xiold(izol)*alphold
     1          /(alpha*beta*(xiold(izol)-xicurr(izol))
     2          +xiold(izol)*alphold*(1.0D0+pz(izol)*alpha))
            alphas(izol)=alpha*xinew(izol)/xicurr(izol)
         enddo
*
         do l=1,mvs
            r(l)=r(l)-alpha*ddtp(l)+tny
*
            do izol=1,kzol
               xs(l,izol)=xs(l,izol)+alphas(izol)*ps(l,izol)
            enddo
         end do
*************************************************************************
*        Test the exit condition                                        *
*************************************************************************
         betan=real(mydot(r,r))
         if (betan*xinew(1)*xinew(1) .lt. relit) go to 30
         if (betan*xinew(kzol)*xinew(kzol) .lt. relit) kzol=kzol-1
*************************************************************************
*        Update the conjugate vector                                    *
*************************************************************************
         beta=betan/betad
         betad=betan
         alphan=betan
*
         !calculate betas for shifted system
         do izol=1,kzol
            betas(izol)=beta*xinew(izol)*alphas(izol)
     1                  /xicurr(izol)/alpha
         enddo
*
         do l=1,mvs
            p(l)=r(l)+beta*p(l)+tny
            do izol=1, kzol
               ps(l,izol)=xinew(izol)*r(l)+betas(izol)*ps(l,izol)+tny
            enddo
         end do
*
         alphold=alpha
         xiold=xicurr
         xicurr=xinew
 20      continue
*
 30   continue
*
      r=r*rnorm
      p=p*rnorm
      x=x*rnorm
      dtp=dtp*rnorm
      ddtp=ddtp*rnorm
      xs=xs*rnorm
      ps=ps*rnorm

      return
      end subroutine  invert_wilson_shifted
