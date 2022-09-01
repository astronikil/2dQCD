      subroutine findeigs
      include 'global.f'
      integer counter, lseed, iev
      real*8 beta, fm, time, rm, rmv, csw, es
      REAL*8 lambda(nev), ipr(nev)
      complex*16 x, y, t, p, wl

      common /wrkry/ x(nvsz), y(nvsz), t(nvsz), p(nvsz)
*
c      open(55,file='eig'//fix(myid),access='append',status='old')
c      open(88,file='back'//fix(myid),status='unknown')
c      open(99,file='hist'//fix(myid),status='unknown')
*
      lseed=137
      call setrn(lseed)
      
      call plaq(es)
      call line(wl)

      call ferminit

      call startritz

      call spectrum(lambda)

      call partratio(ipr)
*
      do iev=1, nev
         write(55,6001)iev,lambda(iev), ipr(iev), mx, es, real(wl)
      enddo
*
 6001 format(1x, i3, 2x, e16.8,2x,e16.8,2x,i3,2x,e16.8,2x,e16.8)
      end
