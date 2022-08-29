      subroutine invert_and_contract(ibool3pt,sigH)
#include "size.h"
#include "tsep.h"
#include "mom.h"
      complex gx, gy, gx_back, pauli, sigH(ncs,ncs)
      complex gx2pt, gx2pt_back
      complex c3ptex, c2ptex, f2ptex
      integer px, py, qx, qy
      character*2 hstr
      common /param_sqrt/ delitsq, nitsq
*
      common /random/ iseed(4)
      common /paulimat/ pauli(ncs,ncs,4)
*
      common /quark/ rm
*      
      common /origin/ iox, ioy, iot
      common /prop_forward/ gx(mv*ncs*ncs)
      common /prop_backward/ gx_back(mv*ncs*ncs)
      common /prop_forward_2pt/ gx2pt(mv*ncs*ncs)
      common /prop_backward_2pt/ gx2pt_back(mv*ncs*ncs)
      common /prop_seq/ gy(mv*ncs*ncs)
      common /corr3ex/ c3ptex(n3pt*       ! ubar_x sigma_i W_x,y u_y, i=1,n3pt
     1                      (2*mx+1)*    ! number of |x-y|
     1                       mt,         ! # of op insertions
     1                       nsep,       ! # of source-sink seps
     1                       nqmom,      ! # of operator mom insertions
     1                       nmom)       ! # of hadron mom 
      common /corr2ex/ c2ptex(mt*
     1                       n2pt,
     1                       nqmom,
     1                       nmom)
      common /fpi2ex/ f2ptex(mt*n2pt,
     1                       nqmom,
     1                       nmom)
*      
      call cpu_time(stime)

c      write(14,*)'Source at',iox,ioy,iot,'ex'

      call propagator(iox,ioy,iot,0,0,"frw","boost",gx(:), gx2pt(:))

      call propagator(iox,ioy,iot,0,0,"bck","boost",gx_back(:),
     1                gx2pt_back(:))
*
      do imom=1, nmom
*
         px=phad(1,imom)
         py=phad(2,imom)
*
         do iqmom=1, nqmom

            qx=qins(1,iqmom)
            qy=qins(2,iqmom)

c            write(14,*)'Computing c2pt:', px+qx, py+qy

            call correl2pt(iox,ioy,iot,px+qx,py+qy,c2ptex(:,iqmom,imom))

c            call fpi2pt(iox,ioy,iot,px+qx,py+qy,f2ptex(:,iqmom,imom))

         enddo
*
         if (ibool3pt .eq. 1) then
           do jsep=1, nsep

              isep=jtau(jsep)
              it3=mod(iot+isep-1,mt)+1
c              write(14,*)'Computing seqprop: dt=', isep

              call seq_propagator(px,py,it3,sigH,gx(:),gy(:))

              do iqmom=1, nqmom
                 qx=qins(1,iqmom)
                 qy=qins(2,iqmom)
c                 write(14,*)'Computing contraction: qx,qy=',qx,qy

                 call correl3pt(qx,qy,iot,sigH,
     1                          c3ptex(:,jsep,iqmom,imom))

              enddo
           enddo
         endif 
      enddo
      call cpu_time(etime)
c      write(14,*)'exact time=',etime-stime
      end
