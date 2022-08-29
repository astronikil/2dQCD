      subroutine seq_propagator(px,      !sink mom
     1                          py,      
     1                          it_sink, !sink t-slice
     1                          sigH,    !sink hadron sigma
     1                          gi,      !forward prop from source
     1                          gp)      !seq prop output
#include "size.h"
      integer px, py
      complex gp(mv,ncs,ncs), gi(mv,ncs,ncs), gf(mv,ncs,ncs)
      complex sigH(ncs,ncs)
      complex v, w, x
      common /wrkry/ v(ncs,mv), w(ncs,mv), x(ncs,mv)
*
      do l=1,mv
         gf(l,:,:)=matmul(sigH(:,:),gi(l,:,:))            
      enddo
*      
      do ispin= 1, ncs
         call lineartocheck(gf(:,:,ispin), v)

         call boost_smear(v,w,"+")

         call seqsource(w,px,py,it_sink,v)

         call boost_smear(v,w,"-")
*         
         call fmv(w,v)
*
         call invert_wilson(v,w,ncg,ncgs)
*
c         write(14,*)'ispin=',ispin, 'ncgo=',ncg, 'ncgw=',ncgs
*
         call checktolinear(w,gp(:,:,ispin))
      enddo
      return
      end
