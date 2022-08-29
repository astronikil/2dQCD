      subroutine propagator(ix,            !source pos
     1                      iy,
     1                      it,
     1                      px,            !source spatial mom
     1                      py,
     1                      qtype,         !frw or bck
     1                      sourcetype,    !smearing S; point or mproj or boost
     1                      gp,            !M^-1.S
     1                      gp2pt)         !S.M^-1.S
#include "size.h"
      integer px, py
      character*5 sourcetype
      character*3 qtype
      complex gp(mv,ncs,ncs), gp2pt(mv,ncs,ncs), v, w, x, gx
      common /wrkry/ v(ncs,mv), w(ncs,mv), x(ncs,mv)
      common /prop_forward/ gx(mv,ncs,ncs)
      common /boostmom/ rqm(2)
*
**************************************************************
*     Don't recompute backward prop if not boost smearing    *
**************************************************************
      if(sourcetype.eq."point".and.qtype.eq."bck")then
        gp=gx
        gp2pt=gp
        goto 20
      endif
*
      if(sourcetype.eq."mproj".and.qtype.eq."bck")then
        gp=gx
        gp2pt=gp
        goto 20
      endif
*
c      if(sourcetype.eq."boost".and.qtype.eq."bck".and.
c     1  abs(rqm(1)).le.rzeroep)then
c        gp=gx
c        gp2pt=gp
c        write(*,*)'rkfac found to consistent with zero.So I didnot do
c     1     back'
c        goto 20
c      endif
*      
      do ispin=1, ncs
**************************************************************
*     Choose source                                          *
**************************************************************
         if(sourcetype.eq."point")then
            call pointsource(v,ix,iy,it,ispin)

         elseif(sourcetype.eq."mproj")then
            call momsource(v,px,py,it,ispin)

         elseif(sourcetype.eq."boost")then
            call pointsource(w,ix,iy,it,ispin)

            if(qtype.eq."frw")then
              call boost_smear(w,v,"+")

            else
              call boost_smear(w,v,"-")

            endif

         endif
**************************************************************
*     Invert Wilson on source                                *
**************************************************************
         call invert_wilson(v,w,ncg,ncgs)
*
c         write(14,*)'ispin=',ispin, 'ncg=',ncg
*
         call fmtv(w,v)
*
**************************************************************
*     create prop for smeared-source point-sink; used for fpi*
**************************************************************
         call checktolinear(v,gp(:,:,ispin))
*
**************************************************************
*     Hit prop with another [Smear] for 2pt contractions     *
**************************************************************

         if(sourcetype.eq."boost")then

            if(qtype.eq."frw")then
              call boost_smear(v,w,"+")

            else
              call boost_smear(v,w,"-")

            endif
         else
            w=v
         endif
*
         call checktolinear(w,gp2pt(:,:,ispin))
      enddo
*
 20   continue
      return
      end subroutine propagator
