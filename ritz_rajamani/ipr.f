      subroutine partratio(ipr)
      include 'global.f'
      integer l, iev
      REAL*8 ipr(nev), rnorm(nev)
      ipr=0.
      rnorm=0.0
      do l=1, nv
         do iev=1, nev
            ipr(iev)=ipr(iev)+abs(ur(l,1,iev))**4+abs(ur(l,2,iev))**4
            rnorm(iev)=rnorm(iev)+
     1      abs(ur(l,1,iev))**2+abs(ur(l,2,iev))**2
         enddo
      enddo

      do iev=1,nev
         ipr(iev)=ipr(iev)/rnorm(iev)**2
      enddo
      return
      end subroutine partratio
      subroutine veigcheck(rmag,rang)
      include 'global.f'
      REAL*8 rmag(nev), rang(nev)
      complex*16 zi, zval, v(nv,2)
      integer iev
      zi=(0.,1.)

      do iev=1,nev
         call vmat(ur(:,:,iev),v)
         call inner(ur(:,:,iev),v,zval)
         rmag(iev)=abs(zval)
         rang(iev)=-zi*log(zval/rmag(iev))
      enddo
      end subroutine veigcheck
c      subroutine writeigs
c      include 'global.f'
c
c      write(92,*)ur(:,:,1)
c      write(92,*)ur(:,:,2)
c      write(92,*)ur(:,:,3)
c      write(92,*)ur(:,:,4)
c
cc      end
