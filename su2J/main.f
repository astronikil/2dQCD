      include 'global.f'

      integer counter,lseed
      common /vsranfom2/ counter
      save /vsranfom2/
      complex*16 u(nc,nc),v(n,n)

c     This is just for test purposes.

c     A random 2X2 SU(2) matrix is generated.      
      lseed=1
      call setrn(lseed)
      call gstart(u)

c     The 2X2 matrix is converted to a special unitary matrix in the nxn representation
c     n=2j+1. If j is an integer, the matrix is real.

c     su2alg generates l3,lp,lm and r. r is set to identity if the representation is even dimensional.
      call su2alg
      
c     This routine converts u to v.
      call fundtoJ(u,v)


      stop
      end
      
