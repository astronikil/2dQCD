      SUBROUTINE spectrum(lambda)

      INCLUDE 'global.f'
      
      complex*16 x,y,t,p
      REAL*8 lambda(nev), rm
      common /quark/ rm
      common /wrkry/ x(nv,2), y(nv,2), t(nv,2), p(nv,2)

      COMPLEX*16 overc
      REAL*8 over,muritz,pnorm,gnorm,gstart,w,overi
      INTEGER ks,iev,i,iflag
      INTEGER is,ic, ll, l, k

      iflag=0
      
      DO ks=1,nks

         DO iev=1,nev
            
!$OMP parallel do default(shared)
            do ic=1,nv
                x(ic,1)=ur(ic,1,iev) 
                x(ic,2)=ur(ic,2,iev)
            enddo
*
            CALL GS(x,iev) 

            call norm(x,over)

            if(over.lt.1.e-18) then
               write(99,*) 'zero vector over at start',over,ks,iev
            endif
            over=1.0/sqrt(over)
*
!$OMP parallel do default(shared)
            do ic=1,nv
               x(ic,1)=over*x(ic,1)
               x(ic,2)=over*x(ic,2)
               ur(ic,1,iev)=x(ic,1)
               ur(ic,2,iev)=x(ic,2)
            enddo
*
            CALL DDAG(x,y) 

            call inner(x,y,overc)

            muritz=overc

            call linearn(y,x,p,muritz)

            CALL GS(p,iev)  

            call norm(p,pnorm)
            pnorm=dsqrt(pnorm)
            if(pnorm.lt.1.e-18) then
               write(99,*) 'zero search at start',pnorm,ks,iev
            endif
            gnorm=pnorm 
	    gstart=gnorm


	    DO i=2,nrt
               CALL ritz(iev,muritz,pnorm,gnorm)
               
	       w=w2*muritz
	       if(w1.gt.w) w=w1
               if(i.gt.10) then
                  if(gnorm.lt.w) go to 2
                  if(gnorm.lt.gstart*cv_fac) go to 2
               endif
               if(mod(i,10).eq.0) then
                  CALL GS(x,iev)  
                  call norm(x,over)
                  if(over.lt.1.e-18) then
                     write(99,*) 'zero vector over',over,ks,iev,i
                  endif
                  over=dsqrt(over)
                  overi=1.d0/over
*
!$OMP parallel do default(shared)
                  do ic=1,nv
                     p(ic,1)=over*p(ic,1)
                     p(ic,2)=over*p(ic,2)
                     x(ic,1)=overi*x(ic,1)
                     x(ic,2)=overi*x(ic,2)
                  enddo
                  over=overi
*
                  CALL DDAG(x,y) 
                  call inner(x,y,overc)
                  muritz=overc
                  call linearn(y,x,t,muritz)
                  call norm(t,gnorm)
                  gnorm=dsqrt(gnorm)

                  call inner(x,p,overc)
                  call transn(p,x,overc)
                  call inner(y,p,overc)
                  overc=(overc/(gnorm*gnorm))-1.d0
                  call transn(p,t,overc)
                  CALL GS(p,iev)  
                  call norm(p,pnorm)
                  pnorm=dsqrt(pnorm)

               endif
               if(pnorm.lt.1.e-18) then
                  write(99,*) 'zero search in Ritz',pnorm,ks,iev,i
               endif
               call flush
            ENDDO

            
 2          continue


            if(i.gt.nrt) then
               write(88,*) 'No of Ritz steps=', nrt
            else
               write(88,*) 'No of Ritz steps=', i
            endif
*
!$OMP parallel do default(shared)
            do ic=1,nv
               ur(ic,1,iev)=x(ic,1)
               ur(ic,2,iev)=x(ic,2)
            enddo
*

         ENDDO


     
         CALL diagritz(lambda,iflag,ks)

         if(iflag.eq.1) go to 1

      ENDDO

 1    continue
*
*
      if(ks.gt.nks) then
         write(88,*) 'No of KS steps= ', nks
      else
         write(88,*) 'No of KS steps= ', ks
      endif

      return
      end
