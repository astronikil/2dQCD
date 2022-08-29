*************************************************************************
*                                                                       *
*     All lattice extents must be even.                                 *
*                                                                       *
*************************************************************************
      parameter (mx=32, my=mx, mt=mx, nd=3)
      parameter (mv=mt*mx*my, mb=nd*mv, mv2=(nd-1)*mv)
      parameter (mvd2=mv/2, msp=mx*my)
      parameter (nspr=2,mvs=nspr*mv,mvs2=2*mvs,mbs=nspr*mb)
      parameter (mux=0, muy=mv, mut=2*mv)
      
      parameter (nex=1, nsl=1)
      parameter (nsep=1)
      parameter (n2pt=4,n3pt=4)

      parameter (nf=1)
      parameter (nfhe=nf/2, nfho=(nf-1)/2)

      parameter(nmeson=2)

      character*2 ip2str(n2pt)
      character*3 ixij(3)
      data ip2str /'V1','V2','V3','S0'/
      data ixij /'x12','x23','x31'/

      character*6 ip3str(n3pt)
      data ip3str /'V1V3V1','V2V3V2','V3V3V3','S0V3S0'/
