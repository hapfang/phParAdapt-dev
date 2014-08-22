     subroutine ludcmp(aatmp,n,np,indx,d) bind(c,name='ludcmp')
!---------------------------------------------------------------------
!
!  find the LU decomposition of a matrix. 
!
!  Numerical Recipes: p. 38
!
!---------------------------------------------------------------------
      use iso_c_binding
      implicit none

      ! input/output - original matrix
      real (c_double), dimension(*) :: aatmp

      ! input - size of original matrix
      integer (c_int) :: n    

      ! input - size of decomposed matrix 
      integer (c_int) :: np   

      ! output - The vector which records the 
      !   row permutation effected by the partial 
      !   pivoting. The values returned for index 
      !   are needed in the call to LUBKSB.
      integer (c_int), dimension(*) :: indx              

      ! output - An indicator of the number of row interchanges
      real (c_double) :: d           

      real (c_double), dimension(50,50) :: aa
      real (c_double), dimension(50) :: vv

      integer idx, i, j, k, imax
      real aamax, sum, dum
      
!
!.... convert the vector to a matrix
!
      idx = 1
      do i=1,n
         do j=1,n
            aa(i,j) = aatmp(idx)
            idx = idx+1
         enddo
      enddo

      d=1.0
      do i=1,n
         aamax=0.0
         do j=1,n
            if (abs(aa(i,j)).gt.aamax) aamax=abs(aa(i,j))
         enddo
         vv(i)=1.0/aamax
      enddo
      
      do j=1,n
         do i=1,j-1
            sum=aa(i,j)
            do k=1,i-1
               sum=sum-aa(i,k)*aa(k,j)
            enddo
            aa(i,j)=sum
         enddo
         aamax = 0.0
         do i=j,n
            sum=aa(i,j)
            do k=1,j-1
               sum=sum-aa(i,k)*aa(k,j)
            enddo
            aa(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax) then
            do k=1,n
               dum=aa(imax,k)
               aa(imax,k)=aa(j,k)
               aa(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(j.ne.n)then
            dum=1.0/aa(j,j)
            do i=j+1,n
               aa(i,j)=aa(i,j)*dum
            enddo
         endif
      enddo

!
!.... convert the matrix back to a vector
!
      idx = 1
      do i=1,n
         do j=1,n
            aatmp(idx) = aa(i,j)
            idx = idx+1
         enddo
      enddo
      
      return
      end subroutine ludcmp
