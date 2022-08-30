      subroutine DDAG(v,w)
      include 'global.f'
      integer, parameter :: ncgmax = 5000
      real*8, parameter :: DACCU1=1.0E-4, DACCU2=1.0E-6
      integer nitov, nitsq, ncgo
      real*8 delitov, delitsq
      complex*16 v(nc,2,nv), w(nc,2,nv), x(nc,2,nv)
      common /param_overlap/ delitov, nitov
      common /param_sqrt/ delitsq, nitsq
*       
      delitov=(DACCU1)**2
      delitsq=(DACCU2)**2
      nitov = ncgmax
      nitsq = ncgmax
*
c      call linkc
c      call fmv_wilson(v,x)
c      call fmtv_wilson(x,w)
      call H_overlap(v,x,ncgo)
      call H_overlap(x,w,ncgo)
*
      return
      end subroutine DDAG
