      subroutine writeconfig(beta, fm, dt, nitrc, delit, ntraj, 
     1           nlen, istart, time, nmd)
************************************************************************
*                                                                      *
*     Write the configuration file                                     *
*                                                                      *
************************************************************************
#include "size.h"
      integer nitrc, ntraj, nmd, istart, iseed, nlen
      real*8 beta, fm, dt, delit, time
      complex*16 u
      common /config/ u(ntot)
      common /random/ iseed(4)
      write(12,*) beta, nf, fm, dt, nitrc, delit, ntraj, nmd,
     2         istart, iseed, time, u
      write(12,*)nmd
      end subroutine writeconfig
      subroutine readconfig(beta, fm, dt, nitrc, delit, ntraj, 
     1           nlen, istart, time, nmd)
************************************************************************
*                                                                      *
*     Read  the configuration file                                     *
*                                                                      *
************************************************************************
#include "size.h"
      integer nflav, nitrc, ntraj, nmd, istart, iseed, nlen
      real*8 beta, fm, dt, delit, time
      complex*16 u
      common /config/ u(ntot)
      common /random/ iseed(4)
      read(10,*) beta, nflav, fm, dt, nitrc, delit, ntraj, nmd,
     2         istart, iseed, time, u
      read(10,*)nmd
      close(10)
      end subroutine readconfig
      subroutine readconfig2(beta, fm, time)
************************************************************************
*                                                                      *
*     Read  the configuration file                                     *
*                                                                      *
************************************************************************
#include "size.h"
      integer nflav, nitrc, ntraj, nmd, istart, iseed, nlen
      real*8 beta, fm, dt, delit, time
      complex*16 u
      common /config/ u(ntot)
      common /random/ iseed(4)
      read(10,*) beta, nflav, fm, dt, nitrc, delit, ntraj, nmd,
     2         istart, iseed, time, u
      read(10,*)nmd
      end subroutine readconfig2
