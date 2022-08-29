      subroutine prc2
#include "size.h"
#include "mom.h"
      complex c2ptex
      common /corr2ex/c2ptex(mt,n2pt,nqmom,nmom)
      do imom=1,nmom
         do iqmom=1,nqmom
            do ic=1,n2pt
               do i=1,mt
                  write(36,*)i,
     1             real(c2ptex(i,ic,iqmom,imom)), 
     1             aimag(c2ptex(i,ic,iqmom,imom)), 
     1             ip2str(ic), 
     1             ' px=',phad(1,imom), ' py=',phad(2,imom),
     1             ' qx=',qins(1,iqmom),' qy=',qins(2,iqmom)
               enddo
            enddo
         enddo
      enddo
      end
      subroutine prc3
#include "size.h"
#include "mom.h"
#include "tsep.h"
      complex c3ptex
      common /corr3ex/c3ptex(-mx:mx,mt,n3pt,nsep,nqmom,nmom)
      idat=0
      do imom=1,1 !nmom
         do iqmom=1,1 !nqmom
            idat=idat+1
c            write(99,*)idat, phad(1,imom), phad(2,imom), qins(1,imom),
c     1        qins(2,imom)
            do iobs=3,3 !1,n3pt
               do its=1,nsep
                  do itau=1,jtau(its)-1 !mt
                     do iz=-mx/2,mx/2
                         write(101,*)iz, 
     1                    -real(c3ptex(iz,itau,iobs,its,iqmom,imom)),
     1                    -aimag(c3ptex(iz,itau,iobs,its,iqmom,imom)),
     1                    ip2str(iobs), 
     1                    ' Ts=', jtau(its),
     1                    ' tau=', itau
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end
