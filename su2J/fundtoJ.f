      subroutine fundtoJ

      include 'globalgen.f'
      integer i1,i2,i3,m, ll, muu
      real*8 theta,snt,th1,th2,th3,j

      complex*16 uin, uout, zphase(nd)
      complex*16 u(nc,nc), v(ncr,ncr),thp,thm,ai
      integer lwork,info
      parameter (lwork=2*ncr-1)
      character*1 jobz,uplo
      complex*16 h(ncr,ncr),work(lwork),sum
      real*8 w(ncr),rwork(3*ncr-2)
      common /configsmr/ uin(nc,nc,mb)
      common /configsmr_irrep/ uout(ncr,ncr,mb)

c     zheev is the LAPACK routine used for eigenvalue and eigenvector computations

      ai=(0.d0,1.d0)
      zphase(1)=exp(ai*acos(-1.0D0)/mx)
      zphase(2)=exp(ai*acos(-1.0D0)/mt)
      jobz='V'
      uplo='U'

c     computation of algebra angles from u

      do ll=1,mb

      muu=(ll-1)/mv+1
      
      u(:,:)=uin(:,:,ll)
      
      theta=dacos(dreal(u(1,1)))
      if(theta.eq.0.d0) then
         do i1=1,ncr
            do i2=1,ncr
               v(i1,i2)=0.d0
            enddo
            v(i1,i1)=1.d0
         enddo
         return
      endif
      snt=dsin(theta)/theta
      th1=dimag(u(1,2))/snt
      th2=dreal(u(1,2))/snt
      thp=0.5d0*(th1-ai*th2)
      thm=0.5d0*(th1+ai*th2)
      th3=dimag(u(1,1))/snt


c     Form the H in U=e^{-iH} and H is the nxn matrix      
c      do i1=1,ncr
c         do i2=1,ncr
c            h(i1,i2)=2.0D0*(thp*lp(i1,i2)+thm*lm(i1,i2)+th3*l3(i1,i2))
c         enddo
c      enddo

      h = 2.0D0*(th1*T1 + th2*T2 + th3*T3)
      !h = 2.0D0*(th1*T2 + th2*T3 + th3*T1)

c     Diagonalize H
      call zheev(jobz,uplo,ncr,h,ncr,w,work,lwork,rwork,info)

c     Form the nxn special unitary matrix in the |j,m> basis.
      
      do i1=1,ncr
         do i2=1,ncr
            sum=0.d0
            do i3=1,ncr
               sum=sum+(h(i1,i3)*dconjg(h(i2,i3))
     &              *cdexp(ai*w(i3)))
            enddo
            v(i1,i2)=sum
         enddo
      enddo

c     Rotate it by r to make in real for an integer representation.

c      if(mod(ncr,2).ne.0) then
c         do i1=1,ncr
c            do i2=1,ncr
c               h(i1,i2)=0.d0
c               do i3=1,ncr
c                  h(i1,i2)=h(i1,i2)+v(i1,i3)*dconjg(r(i2,i3))
c               enddo
c            enddo
c         enddo
c         do i1=1,ncr
c            do i2=1,ncr
c               v(i1,i2)=0.d0
c               do i3=1,ncr
c                  v(i1,i2)=v(i1,i2)+r(i1,i3)*h(i3,i2)
c               enddo
c            enddo
c         enddo
c      endif

      uout(:,:,ll) = v(:,:)*zphase(muu)

c      print*, 'The nxn special unitary matrix'
c      do i1=1,ncr
c         do i2=1,ncr
c            print*, i1,i2,v(i1,i2)
c         enddo
c      enddo
c      print*
c      print*, 'check orthogonality'
c      do i1=1,ncr
c         do i2=1,ncr
c            sum=0.d0
c            do i3=1,ncr
c               sum=sum+v(i1,i3)*dconjg(v(i2,i3))
c            enddo
c            print*, i1,i2,sum
c         enddo
c      enddo

      enddo  !end loop over lattice sites
               
      return
      end
      
