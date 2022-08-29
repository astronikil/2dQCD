      subroutine pikart(n,arr)
      implicit none
      integer n
      real*8 arr(n)
      integer i,j
      real*8 a
      do j=2,n
         a=arr(j)
         do i=j-1,1,-1
            if(arr(i).le.a) go to 10
            arr(i+1)=arr(i)
         enddo
         i=0
 10      arr(i+1)=a
      enddo
      return
      end
      subroutine pikart2(n,arr,mat)
      implicit none
      integer n
      real*8 arr(n)
      complex*16 mat(n,n),vec(n)
      integer i,j,k
      real*8 a
      do j=2,n
         a=arr(j)
         do k=1,n
            vec(k)=mat(k,j)
         enddo
         do i=j-1,1,-1
            if(arr(i).le.a) go to 10
            arr(i+1)=arr(i)
            do k=1,n
               mat(k,i+1)=mat(k,i)
            enddo
         enddo
         i=0
 10      arr(i+1)=a
         do k=1,n
            mat(k,i+1)=vec(k)
         enddo
      enddo
      return
      end
