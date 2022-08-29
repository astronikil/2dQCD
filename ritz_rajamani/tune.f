include 'global.f'
      real*8 beta, fm, time, rmi, rmf, rtune, errl ,errm, eps
      real*8 lamb, csw
      integer counter,lseed, nlev
      parameter(eps=0.6, nlev=1)
      
      read(*,*)rmi,rmf, csw
      
      call addrc

      call readconfig(beta, fm, time)

      call ferminit(csw)

      lseed=137
      call setrn(lseed)

      call tune(rmi,rmf,rtune,errl,errm,time,lamb)

      write(65,*)rtune, errm, errl, lamb

      stop

      end
      subroutine tune(rmi,rmf,rmm,rat,dw,time,lamb)
      include 'global.f'
      real*8 fm, rmin, rmi, rmf, rmm, dw, rat, rm, tol, err,time
      real*8 lambdai(nev), lambdaf(nev), lambdam(nev), ln, lo, errl
      real*8 tolw, tolw2, lamb
      integer itmax, iter
      parameter (itmax=40, tol=1.0E-2, tolw=1.0E-3, tolw2=1.0E-4)
      common /quark/ rm
*
      rmm=(rmi+rmf)/2.0
      dw=(rmf-rmi)/2.0
*
      do iter=1, itmax 
*
         call lambda3pts(rmi,rmm,rmf,lambdai,lambdam,lambdaf)
*
         call parabola(rmi,rmm,rmf,lambdai(1),lambdam(1),lambdaf(1),
     1        rmin)
*
         if(lambdai(1)-lambdam(1).gt.0.0.and.
     1      lambdaf(1)-lambdam(1).gt.0.0)then
            rat=abs(0.5*(lambdaf(1)+lambdai(1))/lambdam(1)-1.)
           if((rat.le.tol).and.(dw.le.tolw).or.(dw.le.tolw2))then
              goto 40
           endif
           dw=dw/2.
           rmm=rmin
           rmi=rmm-dw
           rmf=rmm+dw
         elseif(lambdai(1)-lambdam(1).lt.0.0)then
             rmi=rmi-dw
             dw=min(2.*dw,0.1)
         else
             rmf=rmf+dw
             dw=min(2.*dw,0.1)
         endif
*
      enddo
 40   continue
*
      lamb=lambdam(1)
      return
      end subroutine tune
      subroutine parabola(x1,x2,x3,y1,y2,y3,rmin)
      real*8 x1, x2, x3, y1, y2, y3, rmin, b, c
*
      b=(x3**2*(-y1 + y2) + x2**2*(y1 - y3) + x1**2*(-y2 + y3))
     1   /((x1 - x2)*(x1 - x3)*(-x2 + x3))
      c=(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))
     1   /((x1 - x2)*(x1 - x3)*(x2 - x3))
*
      rmin=-0.5*b/c
*
      return
      end subroutine parabola
      subroutine lambda3pts(rmi,rmm,rmf,lambdai,lambdam,lambdaf)
      include 'global.f'
      real*8 lambdai(nev), lambdaf(nev), lambdam(nev)
      real*8 rmi, rmm, rmf, rm
      common /quark/ rm
*
      rm=-3.+rmi
      call startritz
      call spectrum(lambdai)
*
      rm=-3.+rmm
      call startritz
      call spectrum(lambdam)
*
      rm=-3.+rmf
      call startritz
      call spectrum(lambdaf)
*
      return
      end subroutine lambda3pts
