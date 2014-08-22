      subroutine lubksb(aatmp,n,np,indx,bb) bind(c,name='lubksb')
!---------------------------------------------------------------------
!
! LU back substitution routine.
!
!---------------------------------------------------------------------
      use iso_c_binding
      implicit none

      ! input - LU decomposition
      real (c_double), dimension(*) :: aatmp

      ! input - size of original matrix
      integer (c_int) :: n    

      ! input - size of decomposed matrix 
      integer (c_int) :: np   

      ! input - row permutation 
      integer (c_int), dimension(*) :: indx              

      ! input/output - solution vector
      real (c_double), dimension(*) :: bb

      real (c_double), dimension(50,50) :: aa
      integer idx, i, j, ii, ll
      real sum

!.... convert the vector to a matrix
!
      idx = 1
      do i=1,n
         do j=1,n
            aa(i,j) = aatmp(idx)
            idx = idx+1
         enddo
      enddo
      
      ii=0
      do i=1,n
         ll=indx(i)
         sum=bb(ll)
         bb(ll)=bb(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-aa(i,j)*bb(j)
            enddo
         else if (sum.ne.0.0) then
            ii=i
         endif
         bb(i)=sum
      enddo
      do i=n,1,-1
         sum=bb(i)
         do j=i+1,n
            sum=sum-aa(i,j)*bb(j)
         enddo
         bb(i)=sum/aa(i,i)
      enddo
      
      return
      end
      
            
      
