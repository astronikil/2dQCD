      subroutine invert_wilson_shifted(rphi,nitcg)
************************************************************************
*                                                                      *
*     Finds the the inverse of H^2+p_i for multiple p_i.               *
*                                                                      *
************************************************************************
#include "size.h"
#include "zolo.h"
      complex ps, xs
      complex p, r, x, dtp, ddtp
      complex mydot
      complex rphi(mvs)
      real xinew(nzol), xicurr(nzol), xiold(nzol)
      real betas(nzol), alphas(nzol)
      real rt(mvs2)
*
      integer idone(nzol)
*
      common /param_sqrt/ delit, nitrc

      common /arrays_wil/ r(mvs), p(mvs), x(mvs),
     1                dtp(mvs), ddtp(mvs)

      common /wrkforce/ ps(mvs,nzol), xs(mvs,nzol)

      equivalence (rt(1),r(1))
      kzol=nzol
      call linkc
*************************************************************************
*        Find U^+ used in D and D^+                                     *
*************************************************************************
!$OMP parallel default(shared)
!$OMP do 
      do l=1,mvs
         x(l)=(0.,0.) 
         r(l)=rphi(l)
         p(l)=rphi(l)
         xs(l,1:nzol)=(0.,0.)
         ps(l,1:nzol)=rphi(l)
      end do
!$OMP end do 
*
      alphan=0.
!$OMP BARRIER
!$OMP do reduction(+:alphan)
      do l=1,mvs2
         alphan=alphan+rt(l)*rt(l)
      end do
!$OMP end do 
!$OMP end parallel
      betad=alphan
      relit=delit
      nitcg=0
*************************************************************************
*                                                                       *
*     The main conjugate gradient loop.                                 *
*                                                                       *
*************************************************************************
      alphold=1.
      xiold=1.
      xicurr=1.
      betas=0.
      beta=0.
      do 20 nx=1,nitrc
         nitcg=nitcg+1
*
*************************************************************************
*        Multiplying by the matrix: M^+ M                               *
*************************************************************************
         call fmv_wilson(p,dtp)
         call fmtv_wilson(dtp,ddtp)
*************************************************************************
*        Update the solution and residual                               *
*************************************************************************
         alphad=0.
!$OMP parallel default(shared)
!$OMP do reduction(+: alphad)
         do l=1, mvs
            alphad=alphad+conjg(p(l))*ddtp(l)
         enddo
!$OMP end do 
!$OMP end parallel
         alpha=alphan/alphad
*
         do izol=1, kzol
            xinew(izol)=xicurr(izol)*xiold(izol)*alphold
     1          /(alpha*beta*(xiold(izol)-xicurr(izol))
     2          +xiold(izol)*alphold*(1+pz(izol)*alpha))
            alphas(izol)=alpha*xinew(izol)/xicurr(izol)
         enddo
*
!$OMP parallel default(shared) private(izol)
!$OMP do 
         do l=1,mvs
            r(l)=r(l)-alpha*ddtp(l)+tny
*
            do izol=1,kzol
               xs(l,izol)=xs(l,izol)+alphas(izol)*ps(l,izol)
            enddo
         end do
!$OMP end do 
!$OMP end parallel
*************************************************************************
*        Test the exit condition                                        *
*************************************************************************
         betan=0.
!$OMP parallel default(shared)
!$OMP do reduction(+: betan)
         do l=1,mvs2
            betan=betan+rt(l)*rt(l)
         enddo
!$OMP end do 
!$OMP end parallel
c         write(99,*)nx,betan*xinew(1)*xinew(1),alpha
c         call flush
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
!$OMP parallel default(shared) private(izol)
!$OMP do 
         do l=1,mvs
            p(l)=r(l)+beta*p(l)+tny
            do izol=1, kzol
               ps(l,izol)=xinew(izol)*r(l)+betas(izol)*ps(l,izol)+tny
            enddo
         end do
!$OMP end do 
!$OMP end parallel
*
         alphold=alpha
         xiold=xicurr
         xicurr=xinew
 20      continue
*
 30   continue
c         write(99,*)'##########'
c         write(91,*)mx
         call flush
*
      return
      end subroutine  invert_wilson_shifted
      subroutine invert_overlap(rphi,xsol,nitcg,ncgs)
*************************************************************************
*                                                                       *
*     The conjugate gradient routine for solving Ax=b.                  *
*     Finds vectors X and R in fermion force.                           *
*                                                                       *
*************************************************************************
#include "size.h"
*
      complex rphi(2,mv), xsol(2,mv)
      complex dtp, ddtp, p, r, x
      complex mydot
      real rt(mvs2)
*
      common /param_overlap/ delit, nitrc
      common /arrays_overlap/ r(2,mv), p(2,mv), x(2,mv),
     1                dtp(2,mv), ddtp(2,mv)
      common /iupidn/ iup(mb), idn(mb)
      equivalence (rt(1),r(1,1))
*************************************************************************
*        Find U^+ used in D and D^+                                     *
*************************************************************************
      call linkc
!$OMP parallel do default(shared)
      do l=1,mv
         x(1,l)=(0.,0.)
         x(2,l)=(0.,0.)
         r(1,l)=rphi(1,l)
         r(2,l)=rphi(2,l)
         p(1,l)=rphi(1,l)
         p(2,l)=rphi(2,l)
      end do
*
      alphan=0.
!$OMP parallel do default(shared) reduction(+: alphan)
      do l=1,mvs2
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
         call fmtmv_overlap(p,ddtp,ncgs)
*************************************************************************
*        Update the solution and residual                               *
*************************************************************************
         alphad=real(mydot(p,ddtp))
         alpha=alphan/alphad
!$OMP parallel do default(shared)
         do l=1,mv
            x(1,l)=x(1,l)+alpha*p(1,l)
            x(2,l)=x(2,l)+alpha*p(2,l)
            r(1,l)=r(1,l)-alpha*ddtp(1,l)
            r(2,l)=r(2,l)-alpha*ddtp(2,l)
         end do
*************************************************************************
*        Test the exit condition                                        *
*************************************************************************
         betan=0.
!$OMP parallel do default(shared) reduction(+: betan)
         do l=1,mvs2
            betan=betan+rt(l)*rt(l)
         end do
c         write(*,*)betan
         if (betan.lt.relit) go to 30
*************************************************************************
*        Update the conjugate vector                                    *
*************************************************************************
         beta=betan/betad
         betad=betan
         alphan=betan
         
!$OMP parallel do default(shared)
         do l=1,mv
            p(1,l)=r(1,l)+beta*p(1,l)
            p(2,l)=r(2,l)+beta*p(2,l)
         end do
 20      continue
*
 30   continue
      xsol=x
*
      call flush
      return
      end
