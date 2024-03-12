      program tunestout
#include "size.h"
#include "improve.h"
      integer counter, lseed, iev, iseed
      integer nitsq,  nitov
      real*8 beta, fm, time, rm, rmv, csw, es1, es2
      real*8 xaniso, delitov, delitsq, fmo
      complex*16 x, y, t, p, wl, ztw
      integer istart, istrt, nmd, nflav, nitrc, ntraj, nlen, nstep,
     1  nmin, nmax
      real*8 dt, pi, tpi,delit,zeta
      common /random/ iseed(4)
      common /twisted/ztw(mv,-nd:nd ,-nd:nd)
      iseed=0; iseed(4)=1
      fm=1.0D0
      fmo=0.0D0
      delitov=1.0E-4
      delitsq=1.0E-6
      nitov=5000
      nitsq=5000
      ztw=1.0D0
*
************************************************************************
*     Address arithmetic.                                              *
************************************************************************
      call addrc
************************************************************************
*     setup sun generators                                             *
************************************************************************
      call readconfig(beta, fm, dt, nitrc, delit, ntraj,
     1           nstep, istrt, time, nmd)
************************************************************************
*     initialise the configuration                                     *
************************************************************************
*
      call plaq(es1)
      es1=1.0-real(es1)/dfloat(mv*nc)/1.0D0

      call smear(nlev)

      call plaq_smear(es2)

      write(*,*)eps4,es1,es2, nc, mx,mt
*
      end
