      include 'size.h'

      
c     l3,lp,lm are the three su(2) algebra generators in the |j,m> basis      
c      real*8 l3(ncr,ncr),lp(ncr,ncr),lm(ncr,ncr)
      complex*16 T1(ncr,ncr), T2(ncr,ncr), T3(ncr,ncr)
      complex*16 Tgen_norm

c     r is the rotation matrix that makes the nxn special unitary matrix real when j is integer      
c      complex*16 r(ncr,ncr)
      
c      common/algebra/l3,lp,lm,r
      common/algebra2/T1, T2, T3
      common/tgennorm/Tgen_norm
