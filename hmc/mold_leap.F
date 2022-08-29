      subroutine traj(beta,dt,ncgomn,ncgomx,ncgsmn,ncgsmx,nmd,
     1  ncgoe, ncgse, ene0,enef,naccept)
#include "size.h"
      integer nmd,naccept,imd
      integer ncgomn,ncgomx,ncgsmn,ncgsmx, ncgo, ncgs, ncgoe, ncgse
      real*8 beta, dt, em0, es0, ef0, ene0, em, es, ef, enef,ebr, et
      real*8 start, finish
************************************************************************
*                                                                      *
*     The molecular dynamics trajectory. In this version of the        *
*     program, a first order leap-frog is used with the momenta        *
*     at integer time steps and links at half integer steps.           *
*                                                                      *
************************************************************************
************************************************************************
*     Initialisation of the molecular dynamics trajectory.             *
************************************************************************
      call gnoise                               ! set the gaussian momenta
      call ferminit
      call phiinit                              ! initialise the phi vector
c      call writeinit
      call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
      call uevol(dt*0.5D0)                        ! evolve the gauge fields
************************************************************************
*     The rest of the steps of the trajectory.                         *
************************************************************************
      do 110 imd=1,nmd
         call fermfo(ncgo,ncgs)                     ! the fermionic force
         ncgomn=min(ncgomn,ncgo)
         ncgomx=max(ncgomx,ncgo)
         ncgsmn=min(ncgsmn,ncgs)
         ncgsmx=max(ncgsmx,ncgs)
         call bosefo(beta,ebr,et)               ! the bosonic force
         call pevol(dt)                         ! evolve the momenta
         call uevol(dt)                         ! evolve the gauge fields
 110     continue
************************************************************************
*     The last step of the trajectory                                  *
************************************************************************
      call fermfo(ncgo,ncgs)                     ! the fermionic force
      ncgomn=min(ncgomn,ncgo)
      ncgomx=max(ncgomx,ncgo)
      ncgsmn=min(ncgsmn,ncgs)
      ncgsmx=max(ncgsmx,ncgs)
      call bosefo(beta,ebr,et)                  ! the bosonic force
      call pevol(dt)
      call uevol(dt*0.5D0)                        ! evolve the gauge fields
      call ferminit
      call energy(beta,em,es,ef,enef,ncgoe,ncgse)
      call metro(ene0,enef,naccept)             ! metropolis step
      write(*,*)ene0,em0,es0,ef0
      write(*,*)enef,em,es,ef
      write(*,*)' '
*
      return
      end
