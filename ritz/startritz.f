      SUBROUTINE startritz

      INCLUDE 'global.f'

      integer iev,i

      real*8 sranf,r1,r2
      complex*16 ai

      ai=(0.d0,1.d0)

      do iev=1,nev
          do i=1,nvsz
             r1=sranf()
             r2=sranf()
             ur(i,iev)=r1+ai*r2
         enddo
      enddo

      return
      end
