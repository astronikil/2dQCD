      subroutine choose(y,i)
*
      INCLUDE 'global.f'
*
      complex*16 y(nvsz)
*
      integer i,ic,is
*
!$OMP parallel do default(shared)
      do ic=1,nvsz
         y(ic)=ur(ic,i) 
      enddo
*
      return
      end
      subroutine replace(y,i)
*
      INCLUDE 'global.f'
*
      complex*16 y(nvsz)
*
      integer i,ic,is
*
!$OMP parallel do default(shared)
      do ic=1,nvsz
         ur(ic,i)=y(ic)
      enddo
*
      return
      end
      subroutine norm(x,res)
*
      INCLUDE 'global.f'
*
      integer ic
      real*8 x(mvs2)
      real*8 res
*
      res=0.d0
!$OMP parallel do default(shared) reduction(+: res)
      do ic=1,mvs2
         res=res+x(ic)**2
      enddo
*
      return
      end
      subroutine inner(x,y,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(mvs),y(mvs),res
      integer ic
*
      res=0.d0
!$OMP parallel do default(shared) reduction(+: res)
      do ic=1,mvs
         res=res+dconjg(x(ic))*y(ic)
      enddo
*
      return
      end
      subroutine scale(x,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(mvs)
      real*8 res
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,mvs
            x(ic)=res*x(ic)
      enddo
*
      return
      end
      subroutine transn(x,y,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(mvs),y(mvs),res
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,mvs
         x(ic)=x(ic)-res*y(ic)
      enddo
*
      return
      end
      subroutine transp(x,y,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(mvs),y(mvs),res
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,mvs
         x(ic)=x(ic)+res*y(ic)
      enddo
*
      return
      end
      subroutine linear(x,y,a,b)
*
      INCLUDE 'global.f'
*
      complex*16 x(mvs),y(mvs)
      real*8 a,b
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,mvs
         x(ic)=a*x(ic)+b*y(ic)
      enddo
*
      return
      end
      subroutine linearn(x,y,z,a)
*
      INCLUDE 'global.f'
*
      complex*16 x(mvs),y(mvs),z(mvs)
      real*8 a
      integer ic,is
*
!$OMP parallel do default(shared)
      do ic=1,mvs
         z(ic)=x(ic)-a*y(ic)
      enddo
*
      return
      end
