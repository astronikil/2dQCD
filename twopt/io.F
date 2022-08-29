      subroutine readconfig(beta, fm, dt, nitrc, delit, ntraj, 
     1           nlen, istart, time, nmd, nstep)
************************************************************************
*                                                                      *
*     Write the configuration file                                     *
*                                                                      *
************************************************************************
#include "size.h"
      integer nflav, nitrc, ntraj, nmd, istart, iseed, nstep, nlen
      real*8 beta, fm, dt, delit, time
      complex*16 u
      common /config/ u(ntot)
      common /random/ iseed(4)
      read(10,*) beta, nflav, fm, dt, nitrc, delit, ntraj, nmd,
     2         istart, iseed, time, u
      end subroutine readconfig
