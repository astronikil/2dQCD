      SUBROUTINE GS(x,iv)

      INCLUDE 'global.f'

      COMPLEX*16 x(nv,2),y(nv,2)

      COMPLEX*16 overc
      INTEGER i,is,ic,iv

      DO i=1,iv-1
         call choose(y,i)
         call inner(y,x,overc)
         call transn(x,y,overc)
      enddo

      RETURN
      END
