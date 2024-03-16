      subroutine su2alg

      include 'globalgen.f'

      integer i1,i2,i3,jj,mm
      real*8 j,m
      complex*16 ai, ss
      complex*16 h1(ncr,ncr), h2(ncr,ncr), h3(ncr,ncr)
      complex*16 r(ncr,ncr), lp(ncr,ncr), lm(ncr,ncr),
     1           l3(ncr,ncr)

      ai=(0.d0,1.d0)
      
c     The three generators in the |j,m> representation are formed
      
      do i1=1,ncr
         do i2=1,ncr
            l3(i1,i2)=0.d0
            lp(i1,i2)=0.d0
            lm(i1,i2)=0.d0
         enddo
      enddo

      j=0.5d0*dfloat(ncr-1)
      do i1=1,ncr
         m=dfloat(i1-1)-j
         l3(i1,i1)=m
         if(i1.lt.ncr) lp(i1+1,i1)=sqrt(j*(j+1.d0)
     1                 -m*m-m)
         if(i1.gt.1)
     1    lm(i1-1,i1)=sqrt(j*(j+1.d0)-
     1                (m*m)+(m))
      enddo

      if(mod(ncr,2).eq.0) then
         do i1=1,ncr
            do i2=1,ncr
               r(i1,i2)=0.d0
            enddo
            r(i1,i1)=1.d0
         enddo
      else
         do i1=1,ncr
            do i2=1,ncr
               r(i1,i2)=0.d0
            enddo
         enddo
         jj=(ncr-1)/2
         do mm=1,jj
            r(mm+jj+1,mm+jj+1)=1.d0/sqrt(2.d0)
            r(-mm+jj+1,mm+jj+1)=ai/sqrt(2.d0)
            if(mod(mm,2).eq.0) then
               r(mm+jj+1,-mm+jj+1)=1.d0/sqrt(2.d0)
               r(-mm+jj+1,-mm+jj+1)=-ai/sqrt(2.d0)
            else
               r(mm+jj+1,-mm+jj+1)=-1.d0/sqrt(2.d0)
               r(-mm+jj+1,-mm+jj+1)=ai/sqrt(2.d0)
            endif
         enddo
         r(jj+1,jj+1)=1.d0
      endif

c     Rotate it by r to make in real for an integer representation.

c      T3 = l3
c      lp = T1 + i T2;  T1 = (lp + lm)/2.
c      lm = T1 - i T3   T2 = (lp - lm)/2i

      T1 = 0.5D0*(lp + lm)
      T2 = -0.5D0*ai*(lp - lm)
      T3 = l3

      if(mod(ncr,2).ne.0) then
         do i1=1,ncr
            do i2=1,ncr
               h1(i1,i2)=0.d0
               h2(i1,i2)=0.d0
               h3(i1,i2)=0.d0
               do i3=1,ncr
                  h1(i1,i2)=h1(i1,i2)+T1(i1,i3)*dconjg(r(i2,i3))
                  h2(i1,i2)=h2(i1,i2)+T2(i1,i3)*dconjg(r(i2,i3))
                  h3(i1,i2)=h3(i1,i2)+T3(i1,i3)*dconjg(r(i2,i3))
               enddo
            enddo
         enddo
         do i1=1,ncr
            do i2=1,ncr
               T1(i1,i2)=0.d0
               T2(i1,i2)=0.d0
               T3(i1,i2)=0.d0
               do i3=1,ncr
                  T1(i1,i2)=T1(i1,i2)+r(i1,i3)*h1(i3,i2)
                  T2(i1,i2)=T2(i1,i2)+r(i1,i3)*h2(i3,i2)
                  T3(i1,i2)=T3(i1,i2)+r(i1,i3)*h3(i3,i2)
               enddo
            enddo
         enddo
      endif

c      h1=0.0D0
c      ss=0.0D0
c      do i1=1,ncr
c         do i2=1,ncr
c            ss = ss + T1(i1,i2)*T2(i2,i1)
c         enddo
c      enddo
c      write(*,*)'trace=',ss
c
c         do i1=1,ncr
c            do i2=1,ncr
c               h1(i1,i2)=0.d0
c               h2(i1,i2)=0.d0
c               h3(i1,i2)=0.d0
c               do i3=1,ncr
c        h1(i1,i2)=h1(i1,i2)+T1(i1,i3)*T2(i3,i2)-T2(i1,i3)*T1(i3,i2)
c        h2(i1,i2)=h2(i1,i2)+T2(i1,i3)*T3(i3,i2)-T3(i1,i3)*T2(i3,i2)
c        h3(i1,i2)=h3(i1,i2)+T3(i1,i3)*T1(i3,i2)-T1(i1,i3)*T3(i3,i2)
c               enddo
c            enddo
c         enddo
c
c        do i1=1,ncr
c           !write(*,*)(h1(i1,i2)-ai*T3(i1,i2),i2=1,ncr)
c           write(*,*)(T3(i1,i2),i2=1,ncr)
c        enddo
c        write(*,*)'###########################'
c        do i1=1,ncr
c           !write(*,*)(h2(i1,i2)-ai*T1(i1,i2),i2=1,ncr)
c           write(*,*)(T1(i1,i2),i2=1,ncr)
c        enddo
c        write(*,*)'###########################'
c        do i1=1,ncr
c           !write(*,*)(h3(i1,i2)-ai*T2(i1,i2),i2=1,ncr)
c           write(*,*)(T2(i1,i2),i2=1,ncr)
c        enddo

      return
      end
      
