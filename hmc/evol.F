      subroutine pevol(dt)
************************************************************************
*                                                                      *
*     Molecular dynamics: evolve the canonical momenta.                *
*                                                                      *
************************************************************************
#include "size.h"
      integer l, ic, jc
      real*8 dt
      complex*16 fo, qf, tr, v(nc,nc)
      common /forces/ fo(nc,nc,mb), qf(nc,nc,mb)
************************************************************************
*     evolve the momenta.                                              *
************************************************************************
      do l=1, mb
         tr = (0.0D0,0.0D0)
         v = (0.0D0,0.0D0)
         do ic=1,nc
            do jc=ic,nc
               v(ic,jc)=fo(ic,jc,l)+conjg(fo(jc,ic,l))
               v(jc,ic)=conjg(v(ic,jc))
            enddo
            tr = tr + v(ic,ic)
         enddo 
         fo(:,:,l) = v(:,:)
*
         do ic=1,nc
            fo(ic,ic,l) = fo(ic,ic,l) - tr/dfloat(nc)
         enddo
      enddo
*
      qf = fo*dt + qf
*
      return
      end
      subroutine uevol(dt)
************************************************************************
*                                                                      *
*     Molecular dynamics: evolve the link matrices.                    *
*                                                                      *
************************************************************************
#include "size.h"
      integer l
      real*8 dt
      complex*16 u, zi
      complex*16 qf, fo, eq
      common /config/ u(nc,nc,mb)
      common /forces/ fo(nc,nc,mb), qf(nc,nc,mb)
      common /expdat/ eq(nc,nc,mb)
      zi=(0.0D0,1.0D0)
************************************************************************
*     Exponentiate the force.                                          *
************************************************************************
      eq=dt*qf
      call expsun !eq=exp(i eq)
************************************************************************
*     Evolve by unitary evolution                                      *
************************************************************************
      do l=1,mb
         u(:,:,l) = matmul(eq(:,:,l), u(:,:,l))
      enddo

c      call unitar
*
      return
      end
