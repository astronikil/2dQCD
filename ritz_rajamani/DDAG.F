      subroutine DDAG(v,w)
      include 'global.f'
c      integer, parameter :: ncgmax = 5000
c      real*8, parameter :: DACCU1=1.0E-4, DACCU2=1.0E-6
      integer nitov, nitsq, ncgo, l
c      real*8 delitov, delitsq
      complex*16 v(nc,2,nv), w(nc,2,nv), x(nc,2,nv)
c      common /param_overlap/ delitov, nitov
c      common /param_sqrt/ delitsq, nitsq
*       
c      delitov=(DACCU1)**2
c      delitsq=(DACCU2)**2
c      nitov = ncgmax
c      nitsq = ncgmax
*
c      call linkc
c      call fmv_wilson(v,x)
c      call fmtv_wilson(x,w)
      do l=1,nv
         v(:,2,l)=(0.0D0,0.0D0)
      enddo

      call H_overlap(v,x,ncgo)
      call H_overlap(x,w,ncgo)

      do l=1,nv
         w(:,2,l)=(0.0D0,0.0D0)
      enddo
*
      return
      end subroutine DDAG
