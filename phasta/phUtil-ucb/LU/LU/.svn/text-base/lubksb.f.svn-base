      subroutine lubksb(aatmp,n,np,indx,bb)
c---------------------------------------------------------------------
c
c LU back substitution routine.
c
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      
      dimension indx(n),aa(50,50),bb(n),aatmp(np*np)
c
c.... convert the vector to a matrix
c
      index = 1
      do i=1,n
         do j=1,n
            aa(i,j) = aatmp(index)
            index = index+1
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
      
            
      
