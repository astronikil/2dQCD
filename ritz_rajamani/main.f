      include 'global.f'
      integer counter, lseed, iev
      integer nitov, nitsq
      integer ibegin, iend, i
      integer icharge, ncharge, rstart
      real*8 beta, fm, time, rm, rmv, csw
      REAL*8 lambda(nev), ipr(nev), rmag(nev), rang(nev)
      real*8 fmo, fmfac1, fmfac2, delitov, delitsq
      real*8 q
      complex*16 x, y, t, p

      integer v1,v2,v3,iproc,ierr
      character tt1,tt2,tt3,fix(0:nprocs-1)*3
      common /quark/ rm
      common /wrkry/ x(nv,2), y(nv,2), t(nv,2), p(nv,2)

      common /qmass_overlap/fmo,fmfac1,fmfac2
      common /param_sqrt/ delitsq, nitsq

      parameter (ncharge=6)
      real*8 qferm(ncharge)
      data qferm /32.0,16.0,8.0,4.0,2.0,1.0/
      !data qferm /1.0,2.0,4.0,8.0,16.0,32.0/
*
        myid=0
        numprocs=1
        read(*,*)rmv, fmo, csw
        read(*,*)delitsq, nitsq
        read(*,*)ibegin, iend
        read(*,*)rstart
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
*
      open(55,file='eig'//fix(myid),access='append',status='old')
      open(77,file='ritzev'//fix(myid),status='unknown')
      open(88,file='back'//fix(myid),status='unknown')
      open(99,file='hist'//fix(myid),status='unknown')
      open(10,file='gauge'//fix(myid),status='unknown')
*

      fmfac1=(1.0+fmo*fmo)/2.
      fmfac2=(1.0-fmo*fmo)/4.

      delitsq=mv*delitsq*delitsq
      
      rm=-3.0+rmv
      
      call addrc

      !call readconfig(beta, fm, time)
      do i=1,ibegin-1
         call readrajamani
      enddo

      do i=ibegin,iend
      call readrajamani

      lseed=137
      call setrn(lseed)

c     define a variable in global called rstart.
c     rstart=0 (generates a random set of starting vectors)
c     rstart \ne 0 (reads from a file. I defined 77 above)
      
      if(rstart.eq.0) then
         call startritz
      else
         read(77,*) ur
      endif

c     start the do loop over charge. I am assuming charge is an input to fermit
      call ferminit

      do icharge=1,ncharge

       !  call startritz

         q=qferm(icharge)

         call initcompactu(q)

         call spectrum(lambda)

         write(77,*) ur
       
         call partratio(ipr)

         call veigcheck(rmag,rang)
 
c      call writeigs

c     write the value of the charge also in 55
         do iev=1, nev
            write(55,*)iev,lambda(iev),ipr(iev),fmo,rmv,rmag(iev),
     1                 rang(iev),' Q=',q
         enddo
      enddo
      enddo
*
c     close the do loop over charge
      
      close(10)
      close(77)
      close(55)
      close(88)
      close(99)
      stop
      end
