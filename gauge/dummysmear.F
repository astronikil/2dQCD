      subroutine smear(nlev)
#include "size.h"
      integer nlev
      complex*16 u, us
      common /config/ u(nc,nc,mb)
      common /configsmr/ us(nc,nc,mb)
      us=u
      end
