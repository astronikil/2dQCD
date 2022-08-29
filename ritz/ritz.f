      SUBROUTINE ritz(iev,muritz,pnorm,gnorm)

      
      INCLUDE 'global.f'

      COMPLEX*16 zz(nvsz)

      complex*16 x,y,t,p

      COMPLEX*16 e1c
      REAL*8 muritz,pnorm,gnorm,tnorm,e1
      REAL*8 s1,s2,s3,a,cosdel,sindel,costh,sinth,betar

      INTEGER iev
      INTEGER is,ic

      common /wrkry/ x(nvsz), y(nvsz), t(nvsz), p(nvsz)

      CALL DDAG(p,zz)

      call inner(p,zz,e1c)
      e1=e1c

      s1=0.5d0*(muritz+(e1/(pnorm*pnorm)))
      s2=0.5d0*(muritz-(e1/(pnorm*pnorm)))
      s3=(gnorm*gnorm)/pnorm
      if(dabs(s2).ge.s3) then
         a=dabs(s2)*dsqrt(1.d0+(s3/s2)**2)
      else
         a=s3*dsqrt(1.d0+(s2/s3)**2)
      endif
      cosdel=s2/a
      sindel=s3/a
      if(cosdel.gt.0) then
         sinth=-sqrt(0.5d0*(1.d0+cosdel))
         costh=-0.5d0*sindel/sinth
      else
         costh=sqrt(0.5d0*(1.d0-cosdel))
         sinth=-0.5d0*sindel/costh
      endif

      muritz=muritz-2.d0*a*sinth*sinth

      sinth=sinth/pnorm

      call linear(x,p,costh,sinth)
      call linear(y,zz,costh,sinth)
      call linearn(y,x,zz,muritz)
      call norm(zz,tnorm)

      tnorm=dsqrt(tnorm)

      betar=(costh*tnorm*tnorm)/(gnorm*gnorm)
      gnorm=tnorm
      if((betar*pnorm*costh).gt.(20.d0*gnorm))betar=0.0

      call inner(x,p,e1c)
      call transn(p,x,e1c)
      call linear(p,zz,betar,1.d0)

      
      CALL GS(p,iev)  

      call norm(p,pnorm)
      
      pnorm=DSQRT(pnorm)

      RETURN
      END
