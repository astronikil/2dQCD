      program findeigs
      include 'global.f'
      integer counter, lseed, iev, iseed
      integer nitsq,  nitov
      real*8 beta, fm, time, rm, rmv, csw, es
      REAL*8 lambda(nev), ipr(nev)
      real*8 xaniso, delitov, delitsq
      complex*16 x, y, t, p, wl

      common /wrkry/ x(nvsz), y(nvsz), t(nvsz), p(nvsz)
      common /random/ iseed(4)
      common /param_overlap/ delitov, nitov
      common /param_sqrt/ delitsq, nitsq
      common /anisotropy/ xaniso
      common /quark/ rm
      iseed=0; iseed(4)=1
      fm=1.0D0
      xaniso=1.0D0
      delitov=1.0E-4
      delitsq=1.0E-6
      nitov=5000
      nitsq=5000

      read(*,*)fm

      rm = 1.0D0 + 1.0D0/xaniso**2 - fm
      delitov=delitov*delitov
      delitsq=delitsq*delitsq
*
c      open(55,file='eig'//fix(myid),access='append',status='old')
c      open(88,file='back'//fix(myid),status='unknown')
c      open(99,file='hist'//fix(myid),status='unknown')
************************************************************************
*     Address arithmetic.                                              *
************************************************************************
      call addrc
************************************************************************
*     setup sun generators                                             *
************************************************************************
      call sungenerators
      call su2alg
      call ransun
      call ransun_pert
************************************************************************
*     initialise the configuration                                     *
************************************************************************
*
      lseed=137
      call setrn(lseed)
      
      call plaq(es)

      write(*,*)'plaq=',es;call flush
      call line(wl)

      call ferminit

      call startritz

      call spectrum(lambda)

      call partratio(ipr)
*
      do iev=1, nev
         write(55,6001)iev,lambda(iev), fm !ipr(iev), mx, es, real(wl)
      enddo
*
 6001 format(1x, i3, 2x, e16.8,2x,e16.8,2x,i3,2x,e16.8,2x,e16.8)
      end
