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
      call writeinit
      call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
c      write(*,*)'T0    :',ene0,em0,es0,ef0
      call uevol(dt*0.5D0)                        ! evolve the gauge fields
c      call ferminit
c      call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
c      write(*,*)'U-half:',ene0,em0,es0,ef0
************************************************************************
*     The rest of the steps of the trajectory.                         *
************************************************************************
      do 110 imd=1,nmd-1
c         write(*,*)dt,nmd
         call fermfo(ncgo,ncgs)                     ! the fermionic force
         ncgomn=min(ncgomn,ncgo)
         ncgomx=max(ncgomx,ncgo)
         ncgsmn=min(ncgsmn,ncgs)
         ncgsmx=max(ncgsmx,ncgs)
         call bosefo(beta)                      ! the bosonic force
         call pevol(dt)                         ! evolve the momenta
c         call ferminit
c         call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
c         write(*,*)'P-full:',ene0,em0,es0,ef0
         call uevol(dt)                         ! evolve the gauge fields
c         call ferminit
c         call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
c         write(*,*)'U-full:',ene0,em0,es0,ef0
 110     continue
************************************************************************
*     The last step of the trajectory                                  *
************************************************************************
      call fermfo(ncgo,ncgs)                     ! the fermionic force
      ncgomn=min(ncgomn,ncgo)
      ncgomx=max(ncgomx,ncgo)
      ncgsmn=min(ncgsmn,ncgs)
      ncgsmx=max(ncgsmx,ncgs)
      call bosefo(beta)                          ! the bosonic force
      call pevol(dt)
c      call ferminit
c      call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
c      write(*,*)'P-full:',ene0,em0,es0,ef0
      call uevol(dt*0.5D0)                        ! evolve the gauge fields
c      call ferminit
c      call energy(beta,em0,es0,ef0,ene0,ncgoe,ncgse)
c      write(*,*)'U-half:',ene0,em0,es0,ef0
      call ferminit
      call energy(beta,em,es,ef,enef,ncgoe,ncgse)
      call metro(ene0,enef,naccept)             ! metropolis step
c      write(*,*)'E0= ',ene0,em0,es0,ef0
c      write(*,*)'EF= ',enef,em,es,ef
c      write(*,*)' '
*
      return
      end
