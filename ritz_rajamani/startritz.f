      SUBROUTINE startritz

      INCLUDE 'global.f'

      integer iev,is,iv

      real*8 sranf,r1,r2
      complex*16 ai

      ai=(0.d0,1.d0)

      do iev=1,nev
         do is=1,2
            do iv=1,nv
               r1=sranf()
               r2=sranf()
               ur(iv,is,iev)=r1+ai*r2
            enddo
         enddo
      enddo

      return
      end

