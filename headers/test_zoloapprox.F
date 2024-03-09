      program main
#include "zolo.h"
      do ix=1,100
         x=(ix-1)*0.005
         y=r_p1by4(0)
         do izol=1,nzol_n1by2
            y=y+r_p1by4(izol)/(x+p_p1by4(izol))
         enddo
         write(*,*)x, y, sqrt(sqrt(x))
      enddo
      end
