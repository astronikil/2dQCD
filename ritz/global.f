      include 'size.h'
      integer nv,nev,nks,nrt, nvsz, nvszev
      integer nv2,nv4
      real*8 w1,w2,cv_fac
      parameter (nv=mv)
      parameter (nks=100,nrt=100)
      parameter (nev=20)
      parameter (nvsz=2*nc*nv, nvszev=nvsz*nev)
      parameter (w1=1.d-6,w2=1.d-7,cv_fac=0.1d0)
      complex*16 ur(nvsz,nev)
      complex*16 up(nvsz,nev)
      common/evec/ur
      common/evect/up
