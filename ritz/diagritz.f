      SUBROUTINE diagritz(lambda,iflag,ks)
      
      INCLUDE 'global.f'
      COMPLEX*16 hm(nev,nev),p(nev),ee(nev)
      REAL*8 d(nev),ep(nev)
      REAL*8 lambda(nev)
      REAL*8 temp
      REAL*8 diff,w, rip
      COMPLEX*16 tempc
      INTEGER ks,i,j,k,iflag
      INTEGER is,ic,ll,l
      complex*16 x,y,t,pp
      common /wrkry/ x(nvsz), y(nvsz), t(nvsz), pp(nvsz)
      common /ipr/ rip(nev)


      DO j=1,nev

         call choose(x,j)

         CALL DDAG(x,y)      

         DO i=1,j
            call choose(t,i)
            call inner(t,y,tempc)
            hm(i,j)=tempc
            hm(j,i)=dconjg(hm(i,j))

         enddo
      ENDDO

      call TRIDIAG (hm,nev,d,ee,ep,p)
      call TQLI(d,ep,nev,hm)      
      
      DO i=1,nev-1
         DO j=i+1,nev
            if(d(j).lt.d(i)) then
               temp=d(i)
               d(i)=d(j)
               d(j)=temp
               DO k=1,nev
                  tempc=hm(k,i)
                  hm(k,i)=hm(k,j)
                  hm(k,j)=tempc
               ENDDO
            endif
         ENDDO
      ENDDO
      iflag=1
      DO i=1,nev
         if(ks.gt.1) then
            diff=dabs(d(i)-lambda(i))
            w=10.d0*w2*dabs(lambda(i))
	    if(10.d0*w1.gt.w) w=10.d0*w1
	    if(diff.gt.w) iflag=0
         else
            iflag=0
         endif
         lambda(i)=d(i)
         write(88,*) 'hov2', i,lambda(i)
      ENDDO


!$OMP parallel do default(shared)
!$OMP^ private(l,k,j)
      do ll=1,nvszev
         l=mod(ll-1,nvsz)+1
         k=(ll-1)/nvsz+1
         up(l,k)=0.d0
         do j=1,nev
            up(l,k)=up(l,k)+dconjg(hm(j,k))*ur(l,j)
         enddo
      enddo 
      ur=up
*
      RETURN
      END
