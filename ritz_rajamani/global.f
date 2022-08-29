      include 'size.h'
      integer nv,nev,nks,nrt, nvev
      integer nv2,nv4
      real*8 w1,w2,cv_fac
      parameter (nks=100,nrt=100)
      parameter (nev=8)
      parameter (w1=1.d-5,w2=10.0*w1,cv_fac=0.1d0)
*
      parameter (nv=mv)
      parameter (nv2=mv2, nv4=4*nv)
      parameter (nvev=nv*nev)
*
      complex*16 ur(nv,2,nc,nev)
      complex*16 up(nv,2,nc,nev)
      common/evec/ur
      common/evect/up
