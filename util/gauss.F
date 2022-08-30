      subroutine gauss(am,nsize,det)
************************************************************************
*                                                                      *
*   Using Gauss elimination, calculates determinant det of matrix g    *
*                        
************************************************************************
      integer i1, i2, k, nsize
      complex*16 am(nsize,nsize)
      complex*16 det,x,y
      det=(1.0D0,0.0D0)
      do i1=1,nsize
         x=am(i1,i1)
         det=det*x
         do i2=i1+1,nsize
            y=am(i2,i1)/x
            do k=i1,nsize
               am(i2,k)=am(i2,k)-y*am(i1,k)
            enddo
         enddo
      enddo
      return
      end
