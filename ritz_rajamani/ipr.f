      subroutine partratio(ipr)
      include 'global.f'
      integer l, iev, is, ic, ll
      REAL*8 ipr(nev), rnorm(nev), ss(nev)
      ipr=0.
      rnorm=0.0
      l=0
      do ll=1, nv
         ss=0.0
         do is=1,2
            do ic=1,nc
               l=l+1
               do iev=1, nev
                  ss(iev)=ss(iev)+abs(ur(l,iev))**2+abs(ur(l,iev))**2
               enddo
            enddo
         enddo
         do iev=1,nev
            ipr(iev)=ipr(iev)+ss(iev)**2
            rnorm(iev)=rnorm(iev)+ss(iev)
         enddo
      enddo

      do iev=1,nev
         ipr(iev)=ipr(iev)/rnorm(iev)**2
      enddo
      return
      end subroutine partratio
