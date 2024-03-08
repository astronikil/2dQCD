      subroutine bosefo(beta)
************************************************************************
*     Molecular dynamics: bosonic part of the force.                   *
************************************************************************
#include "size.h"
*
      integer l
      complex*16 fo, qf, zi
      complex*16 st
      real*8 beta
*
      common /forces/ fo(nc2,mb), qf(nc2,mb)
      common /donfig/ st(nc2,mb)
      zi=(0.0D0,1.0D0)
************************************************************************
*     Construct the plaquette values.                                  *
************************************************************************
      call rest_in
************************************************************************
*     Construct the bosonic part of the driving force.                 *
************************************************************************
      do l=1,mb
         fo(:,l)=fo(:,l)+zi*beta*st(:,l)
      enddo
*
      return
      end
