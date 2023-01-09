      subroutine init(istart,time,dt,nmd,nstep)
************************************************************************
*     initialises constants                                            *
*     call addrc to initialise pointers and fermion phases             *
*     initialise gauge field configuration                             *
*        istart < 0  read old configuration                            *
*        istart = 0  ordered start                                     *
*        istart > 0  random start                                      *
************************************************************************
#include "size.h"
      integer istart, istrt, nmd, nflav, nitrc, ntraj, nlen, nstep,
     1  nmin, nmax
      real*8 time, dt, pi, tpi, fm, beta, delit
      complex*16 r1, r2, tpiz, zi
      common /consts/ pi, tpi, tpiz, zi
************************************************************************
*     calculate constants                                              *
************************************************************************
      zi=(0.0D0,1.0D0)
      pi=acos(-1.0D0)
      tpi=2.0D0*pi
      tpiz=tpi*zi
************************************************************************
*     Address arithmetic.                                              *
************************************************************************
      call addrc
************************************************************************
*     setup sun generators                                             *
************************************************************************
      call sungenerators
************************************************************************
*     initialise the configuration                                     *
************************************************************************
      if (istart.gt.0) then			! hot start
         time=0.0D0
         call setseed
         call ransun
      else if (istart.eq.0) then		! cold start
         time=0.0D0
         call setseed
         call ordsun
      else if (istart.eq.-1. or. istart.eq.-2) then				! old start
         call readconfig(beta, fm, dt, nitrc, delit, ntraj,
     1           nstep, istrt, time, nmd)

         if(istart.eq.-1)then
           write(28,'(5x,a)') 'Initialised to gaugeold'
         endif
         if(istart.eq.-2)then
           write(28,'(5x,a)') 'Initialised to end of previous run'
         endif
c
      else if (istart.eq.-3) then				! old start
         call readconfig(beta, fm, dt, nitrc, delit, ntraj,
     1           nstep, istrt, time, nmd)
         call setseed
         write(28,'(5x,a)') 
     1     'Initialised to gaugeold and new seed'
      endif
*
      return
 6001 format(2x,' trajectories, MD time=',f14.2)
      end
      subroutine setseed
      integer, parameter :: nmax=1000000, nrand=4095
      integer iseed, ncnf, pos
      real*8 r(nmax)
      common /random/ iseed(4)
      common /cnfno/ ncnf
      iseed=0; iseed(4)=1
      call dlarnv(1,iseed,nmax,r)
      lpos=4*(ncnf-1)+1
      iseed(1)=r(lpos)*nrand
      iseed(2)=r(lpos+1)*nrand
      iseed(3)=r(lpos+3)*nrand
      iseed(4)=r(lpos+4)*nrand
      if(mod(iseed(4),2).eq.0)then
        iseed(4)=iseed(4)+1
      endif
      end subroutine setseed
