*************************************************************************
*                                                                       *
*     All lattice extents must be even.                                 *
*                                                                       *
*************************************************************************
      parameter (mx=LLL, my=mx, mt=mx, nd=3)
      parameter (mv=mt*mx*my, mb=nd*mv, mv2=(nd-1)*mv)
      parameter (mvd2=mv/2, msp=mx*my)
      parameter (nspr=2,mvs=nspr*mv,mvs2=2*mvs,mbs=nspr*mb)
      parameter (mux=0, muy=mv, mut=2*mv)

      parameter (nf=1)
      parameter (nfhe=nf/2, nfho=(nf-1)/2)

      parameter(nmeson=2)
