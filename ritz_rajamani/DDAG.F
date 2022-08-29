      subroutine DDAG(v,w)
      include 'global.f'
      integer ncg
      complex*16 v(nv,2), w(nv,2)
*
      call linkc
      call fmtmv_overlap(v,w,ncg)
*
      return
      end subroutine DDAG
