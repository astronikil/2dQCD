      subroutine DDAG(v,w)
      include 'global.f'
      integer ncgo
      complex*16 v(nc,2,nv), w(nc,2,nv), x(nc,2,nv)
*
      call linkc
c      call fmv_wilson(v,x)
c      call fmtv_wilson(x,w)
      call H_overlap(v,x,ncgo)
      call H_overlap(x,w,ncgo)
*
      return
      end subroutine DDAG
