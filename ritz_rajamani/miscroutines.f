      subroutine choose(y,i)
*
      INCLUDE 'global.f'
*
      complex*16 y(nv,2)
*
      integer i,ic,is
*
!$OMP parallel do default(shared)
      do ic=1,nv
            y(ic,1)=ur(ic,1,i) 
            y(ic,2)=ur(ic,2,i)
      enddo
*
      return
      end
      subroutine replace(y,i)
*
      INCLUDE 'global.f'
*
      complex*16 y(nv,2)
*
      integer i,ic,is
*
!$OMP parallel do default(shared)
      do ic=1,nv
         ur(ic,1,i)=y(ic,1)
         ur(ic,2,i)=y(ic,2)
      enddo
*
      return
      end
      subroutine norm(x,res)
*
      INCLUDE 'global.f'
*
      integer ic
      real*8 x(nv4)
      real*8 res
*
      res=0.d0
!$OMP parallel do default(shared) reduction(+: res)
      do ic=1,nv4
         res=res+x(ic)**2
      enddo
*
      return
      end
      subroutine inner(x,y,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(nv2),y(nv2),res
      integer ic
*
      res=0.d0
!$OMP parallel do default(shared) reduction(+: res)
      do ic=1,nv2
         res=res+dconjg(x(ic))*y(ic)
      enddo
*
      return
      end
      subroutine scale(x,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(nv2)
      real*8 res
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,nv2
            x(ic)=res*x(ic)
      enddo
*
      return
      end
      subroutine transn(x,y,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(nv2),y(nv2),res
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,nv2
         x(ic)=x(ic)-res*y(ic)
      enddo
*
      return
      end
      subroutine transp(x,y,res)
*
      INCLUDE 'global.f'
*
      complex*16 x(nv2),y(nv2),res
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,nv2
         x(ic)=x(ic)+res*y(ic)
      enddo
*
      return
      end
      subroutine linear(x,y,a,b)
*
      INCLUDE 'global.f'
*
      complex*16 x(nv2),y(nv2)
      real*8 a,b
      integer ic
*
!$OMP parallel do default(shared)
      do ic=1,nv2
         x(ic)=a*x(ic)+b*y(ic)
      enddo
*
      return
      end
      subroutine linearn(x,y,z,a)
*
      INCLUDE 'global.f'
*
      complex*16 x(nv2),y(nv2),z(nv2)
      real*8 a
      integer ic,is
*
!$OMP parallel do default(shared)
      do ic=1,nv2
         z(ic)=x(ic)-a*y(ic)
      enddo
*
      return
      end
