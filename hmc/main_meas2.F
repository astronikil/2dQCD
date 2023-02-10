       program hmc
************************************************************************
*                                                                      *
*      Hybrid Monte Carlo for SU(Nc) theory coupled to massless        *
*      overlap fermions in 2d                                          *
*                                                                      *
************************************************************************
#include "size.h"
#include "vers.h"
#include "improve.h"
*
      integer iseed, ncnf, nmd, nmp, nitov, nitsq, istart,
     1        nodat, ncgomn, ncgomx, ncgsmn, ncgsmx, nmt, nmdp, nr,
     1        iter, ntraj, nitrc, nstep, nadaptstep, nadaptnmd, nlen,
     1        ilen, istep, l, nmd_min, nstep_new, nmd_new, initconf,
     1        iconf, nconf, itune,  ncgoemn, ncgoemx, ncgsemn,
     1        ncgsemx, ncgoe, ncgse, nmon, inewseed, ntop3, nstring,
     1        iboolwrite, ntherm, nconfnotune
        
      real*8 ua, uas, uao, rm, fm, fmo, betag, dt, delitov, delitsq,
     1       pi, tpi, fmfac1, fmfac2, beta, time, esp, racc, delit,
     1       ene0, enef, tmon, tmonsmr, rmo, start, finish, avd, erd,
     1       startcur, finishcur, h1, zeta, xaniso

      complex*16 u, ud, uo, tpiz, wl, zi, sact, sact2, ppb

      integer v1,v2,v3,iproc,ierr
      character tt1,tt2,tt3,fix(0:nprocs-1)*3
      character*19 gstring
      character*12 sio
      character*18 eeigstring
      character*19 ehiststring
      character*19 ebackstring
      character*19 messtring
*
      common /param_overlap/ delitov, nitov
      common /param_sqrt/ delitsq, nitsq
*
      common /config/ u(nc,nc,mb)
      common /oldconf/uo(nc,nc,mb)
*
      common /random/ iseed(4)
      common /consts/ pi, tpi, tpiz, zi
      common /anisotropy/ xaniso
      common /quark/ rm
      common /cnfno/ncnf
************************************************************************
*                                                                      *
*     Read in parameters.                                              *
*        beta		  coupling constant                            *
*        nflav		  number of flavours                           *
*        fm		  fermion mass                                 *
*        dt		  time step in hmc trajectory                  *
*        nitrc		  number of congrad iterations                 *
*        delit		  congrad stopping criterion                   *
*        ntraj		  number of trajectories                       *
*        nmd		  number of md steps per trajectory            *
*        istart		  <0 old conf; =0 cold; >0 hot                 *
*        nodat		  <=0 old data file, >0 new data file          *
*        nmp              =1 metropolis step;!=0 no metropolis         *
*                                                                      *
************************************************************************
*
      myid=0
      numprocs=1
c      call MPI_INIT(ierr)
c      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
*
      if(myid.eq.0)then
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
*
      open(27,file='inp'//fix(myid),status='unknown')
      read(27,*)betag, fm, xaniso
      read(27,*)nitov, delitov, nitsq, delitsq
      read(27,*)istart, nodat, iboolwrite
      read(27,*)nmd, nconf, ntraj, nmp, nadaptnmd
      read(27,*)ncnf
      read(27,*)ntherm
      read(27,*)nconfnotune
      close(27)
*
      write(sio,6102)fix(myid)
      open(31,file=sio,status='unknown')
      read(31,*)initconf
      close(31)
*      
      beta = dfloat(mx*mt)/betag/xaniso
      rm = 1.0D0 + 1.0D0/xaniso - fm
*
      delitov=delitov*delitov
      delitsq=delitsq*delitsq
*
      !openfiles
      if(nodat.gt.0)then
         open(14,file='data'//fix(myid),status='unknown')
      else
         open(14,file='data'//fix(myid),access='append',status='old')
      endif
*
      if(istart.eq.-1.or.istart.eq.-3)then
         open(10,file='gaugeold'//fix(myid),status='old')
      elseif(istart.eq.-2)then
         write(gstring,6101)fix(myid),initconf
         open(10,file=gstring,status='old')
      endif   
      open(28,file='outrun'//fix(myid),access='append',status='old')
*     
      !Initialize..read files.. set params...call addressing...
      call init(istart,time,dt,nmd,nstep)
*
      nmd=nmd*nstep
      nstep=1
      dt=1.0D0/dfloat(nmd)
*
      write(28,6000) version
      write(28,6001) mx,mt,nf,fm,betag,ntraj,nmd,dt,delit,eps
************************************************************************
*                                                                      *
*     Main loop: hybrid monte carlo update                             *
*                                                                      *
************************************************************************
      nmt=0
*
      do iconf=initconf+1, nconf
*
         call unitar
*
         nr=0
         do iter=1,ntraj
            ncgomn=nitov
            ncgomx=0
            ncgsmn=nitsq
            ncgsmx=0
            ncgoemn=nitov
            ncgoemx=0
            ncgsemn=nitsq
            ncgsemx=0

            ilen=0
            nmd_min=nmd
*
            uo=u      !Copy old configutation in case of mc rejection
*
            dt=1.0D0/dfloat(nmd) !MC step size for
                                                 !trajectory of size nmd which might have changed in previous
                                                 !step
            ilen=ilen+nmd !Keep track of actual number of steps taken to
                          !complete t=1.0
*
            call traj(beta,dt,ncgomn,ncgomx,ncgsmn,ncgsmx,nmd, 
     1                ncgoe,ncgse,ene0, enef, nmt)    !Take nmd molecular steps
*
            ncgoemn=min(ncgoemn,ncgoe)
            ncgoemx=max(ncgoemx,ncgoe)
            ncgsemn=min(ncgsemn,ncgse)
            ncgsemx=max(ncgsemx,ncgse)
*

            if(iconf.ge.nconfnotune)then
               call adaptstep(ene0,enef,nmd,nmd_new) !Tune nmd
*
               if(nmt.eq.0.and.nmp.eq.1)then 
                  u=uo    !Copy back old config if rejected
                  nr=nr+1 !increase reject counter
               endif 
*
               if(nadaptnmd.eq.1)then !tune nmd
                  nmd=nmd_new
               endif
            endif
*
            nlen=ilen  !Total number of dt steps
            time=time+1.0D0
************************************************************************
*        Measurments                                                   *
************************************************************************
            call line(wl)  !Polyakov loop
            call plaq(esp)
            sact=1.0-real(esp)/dfloat(mv*nc)/1.0D0
************************************************************************
*        Write on the data file.                                       *
************************************************************************
            write(14,*)real(sact), real(wl), aimag(wl)
************************************************************************
*        Write on the output file.                                     *
************************************************************************
            write(28,6002) iter, time, 
     1                  ncgomn, ncgomx, ncgsmn, ncgsmx, 
     1                  ncgoemn, ncgoemx, ncgsemn, ncgsemx, 
     1                  dt, nlen, sact, wl, ene0, enef
            call flush
         enddo
         racc=1.0D0-dfloat(nr)/dfloat(ntraj)
         write(28,6003)racc
************************************************************************
*                                                                      *
*        Measurement of fermionic observables           	       *
*                                                                      *
************************************************************************

         if(iconf.gt.ntherm)then

         write(eeigstring,6103)fix(myid),iconf
         write(ebackstring,6104)fix(myid),iconf
         write(ehiststring,6105)fix(myid),iconf
         open(55,file=eeigstring,status='new')
         open(88,file=ebackstring,status='new')
         open(99,file=ehiststring,status='new')
*
         call findeigs
*
         close(55)
         close(88)
         close(99)
*         
         write(messtring,6106)fix(myid),iconf
         open(36,file=messtring,status='new')
*
         call meson
*
         close(36)

         endif
************************************************************************
*                                                                      *
*        Write out the final configuration 			       *
*                                                                      *
************************************************************************
         write(sio,6102)fix(myid)
         open(31,file=sio,status='replace')
         write(31,*)iconf
         close(31)
*      
         if(iboolwrite.eq.1)then
            write(gstring,6101)fix(myid),iconf
            open(12,file=gstring,status='new')
            call writeconfig(beta, fm, dt, nitrc, delit, ntraj, 
     1           nstep, istart, time, nmd)
            close(12)
         else
            open(12,file='gaugeold'//fix(myid),status='replace')
            call writeconfig(beta, fm, dt, nitrc, delit, ntraj, 
     1           nstep, istart, time, nmd)
            close(12)
         endif
         call flush

      enddo !end loop over iconf
*      
      close(14)
      close(28)
      close(63)
      close(22)
c      call MPI_FINALIZE(ierr)
*
      stop
************************************************************************
*     Format statements.                                               *
************************************************************************
 6000 format(//2x,'HMC: program version ',f8.3)

 6001 format(2x,'Lattice size: nx=',i2,' ny=',i2/
     1 2x,i1,' flavours of fermions with mass=',e12.4,2x,'Beta=',e12.4/
     1 2x,'Run',i3,' trajectories with nmd=',i3,' and dt=',e12.4/
     1 2x,'Conjugate gradient stopping criterion=',e12.4/
     1 2x,'Stout smearing eps1=',f8.3)

 6002 format(/2x,'Trajectory',i6,t46,'MD time=',f16.2,t46,/
     1 2x,'ncgo=',t11,2i6/
     1 2x,'ncgs=',t11,2i6/
     1 2x,'ncgoE=',t11,2i6/
     1 2x,'ncgsE=',t11,2i6/
     1 2x,'dt=',t11,e12.4,t46,'nmd=',t56,i6/
     1 2x,'eb=',t11,2e12.4, t46,'wl=',2e12.4/
     1 2x,'ene=',t11,2e12.4)

 6003 format(//2x,'Acceptance rate= ',f8.3)

 6101 format('./GAUGE/GG',A3,'.',i5.5)     
 6102 format('inputfile',A3)     
 6103 format('./EIG/eig',A3,'.',i5.5)     
 6104 format('./EIG/back',A3,'.',i5.5)     
 6105 format('./EIG/hist',A3,'.',i5.5)     
 6106 format('./MES/corr',A3,'.',i5.5)     

 7000 format(e14.6)
      end program
