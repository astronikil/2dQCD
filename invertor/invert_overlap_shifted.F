      subroutine invert_overlap_shifted(rphi,nzolin,nitcg,ncgs)
************************************************************************
*                                                                      *
*     Finds the the inverse of H^2+p_i for multiple p_i.               *
*                                                                      *
************************************************************************
#include "size.h"
c#include "zolo.h"
      integer, parameter :: nzolmax=30
      real*8, parameter :: tny=1.0E-15
      integer kzol, l, nitrc, nitcg, ncgs, nx, izol, nzolin
      real*8 pz, rz
      complex*16 ps, xs
      complex*16 p, r, x, dtp, ddtp
      complex*16 mydot, mychiraldot
      complex*16 rphi(mvs)
      real*8 xinew(nzolmax), xicurr(nzolmax), xiold(nzolmax)
      real*8 betas(nzolmax), alphas(nzolmax)
      real*8 rt(mvs2)
      real*8 delit, relit
      real*8 alphan, alphad, betan, betad, beta, alpha, alphold
      real*8 rnorm, rchinorm
      complex*16 rphinorm(mvs)
*
      common /param_overlap/ delit, nitrc
*
      common /arrays_ovshift/ r(mvs), p(mvs), x(mvs),
     1                dtp(mvs), ddtp(mvs)
*
      common /wrkforce_ovshift/ ps(mvs,nzolmax), xs(mvs,nzolmax)
*      
      common /zolvals/ pz(nzolmax), rz(nzolmax)
*
      equivalence (rt(1),r(1))
*
      rnorm=sqrt(real(mydot(rphi,rphi)))
      rchinorm=real(mychiraldot(rphi,rphi))/rnorm/rnorm
      rphinorm=rphi/rnorm
*
      kzol=nzolin
*************************************************************************
*        Find U^+ used in D and D^+                                     *
*************************************************************************
      do l=1,mvs
         x(l)=(0.0D0,0.0D0) 
         r(l)=rphinorm(l)
         p(l)=rphinorm(l)
         xs(l,1:kzol)=(0.0D0,0.0D0)
         ps(l,1:kzol)=rphinorm(l)
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
         if(rchinorm.lt.relit)then
c            write(*,*)'switch to overlap_square:',rchinorm
            call H_overlap_square(p,ddtp,ncgs)
         else
c            write(*,*)'switch to exact overlap2:',rchinorm
            call H_overlap(p,dtp,ncgs)
            call H_overlap(dtp,ddtp,ncgs)
         endif
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
      end subroutine  invert_overlap_shifted
