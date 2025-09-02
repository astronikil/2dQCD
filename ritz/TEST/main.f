      include 'global.f'
      integer counter, lseed, iev
      real*8 beta, fm, time, rm, rmv, csw, es
      REAL*8 lambda(nev), ipr(nev)
      complex*16 x, y, t, p, wl

      integer v1,v2,v3,iproc,ierr
      character tt1,tt2,tt3,fix(0:nprocs-1)*3
      common /quark/ rm
      common /wrkry/ x(nvsz), y(nvsz), t(nvsz), p(nvsz)

      myid=0
      numprocs=1
*
c      call MPI_INIT(ierr)
c      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
*
      if(myid.eq.0)then
         read(*,*)rmv
         do iproc=0,numprocs-1
            v1=iproc/100
            v2=iproc/10
            v3=v2/10
            v2=v2-10*v3
            v3=iproc-v2*10-v1*100
            tt1=char(48+v1)
            tt2=char(48+v2)
            tt3=char(48+v3)
            fix(iproc)=tt1//tt2//tt3
         enddo
      endif
*
c      call MPI_BCAST(fix(0),nprocs*3,MPI_CHARACTER,0,
c     &     MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(rmv,1,MPI_DOUBLE_PRECISION,0,
c     &     MPI_COMM_WORLD,ierr)
*
      open(55,file='eig'//fix(myid),access='append',status='old')
      open(88,file='back'//fix(myid),status='unknown')
      open(99,file='hist'//fix(myid),status='unknown')
      open(10,file='gauge'//fix(myid),status='unknown')
*
      rm=2.0-rmv
*
      lseed=137
      call setrn(lseed)
      
      call addrc

      call ordsun

      !call readconfig2(beta, fm, time) 

      call plaq(es)
      call line(wl)

      call ferminit

      call startritz

      call spectrum(lambda)

      call partratio(ipr)
*
      do iev=1, nev
         write(55,6001)iev,lambda(iev), ipr(iev), mx, es, real(wl),rmv
      enddo
*
      close(10)
      close(55)
      close(88)
      close(99)
c      call MPI_FINALIZE(ierr)
      stop
 6001 format(1x, i3, 2x, e16.8,2x,e16.8,2x,i3,2x,e16.8,2x,e16.8,2x,
     1       e16.8)
      end
