      subroutine gstart(u)

      include 'global.f'
      
      integer ic,jc,kc
      integer lwork,info
      parameter (lwork=2*nc-1)
      character*1 jobz,uplo
      complex*16 h(nc,nc),work(lwork),sum
      real*8 w(nc),rwork(3*nc-2)
      real*8 rr,th,sranf,trace,pi
      complex*16 ai,c,u(nc,nc)


      ai=(0.d0,1.d0)
      pi=dacos(-1.d0)
      
      jobz='V'
      uplo='U'

      trace=0.d0
      do ic=1,nc
         do jc=ic,nc
            rr=sranf()
            th=2.d0*pi*sranf()
            c=dsqrt(-dlog(rr))*cdexp(ai*th)
            if(ic.ne.jc) then
               h(ic,jc)=c
               h(jc,ic)=dconjg(c)
            else
               h(ic,jc)=dreal(c)
               trace=trace+dreal(c)
            endif
         enddo
      enddo
      trace=trace/dfloat(nc)
      do ic=1,nc
         h(ic,ic)=h(ic,ic)-trace
      enddo
      call zheev(jobz,uplo,nc,h,nc,w,work,lwork,rwork,info)
      do ic=1,nc
         do jc=1,nc
            sum=0.d0
            do kc=1,nc
               sum=sum+(h(ic,kc)*dconjg(h(jc,kc))
     &              *cdexp(ai*w(kc)))
            enddo
            u(ic,jc)=sum
         enddo
      enddo

      return
      end
      
