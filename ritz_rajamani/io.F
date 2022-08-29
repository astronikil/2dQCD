      subroutine readrajamani
      include 'global.f'
      integer i1, i2, i3, id, l, ll, ino
      real*8 theta, pifac
      real*8 ua
      common /configa/ ua(mb)
      common /newold/ ino(mv)
      pifac=acos(-1.0)/mx
      do l=1,mb
         read(10,*)i1,i2,i3,id,theta
         ll=ino(i1+mx*(i2-1+my*(i3-1)))
         ua(ll+(id-1)*mv)=theta
      enddo
      !do l=1,mb
      !   ua(l)=ua(l)+pifac
      !enddo
      end
