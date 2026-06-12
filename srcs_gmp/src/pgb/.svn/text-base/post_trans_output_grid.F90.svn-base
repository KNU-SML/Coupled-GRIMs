!
   subroutine post_trans_output_grid(io,jo,a,lot)
!-------------------------------------------------------------------------------
   integer              ::  io,jo,lot
   real                 ::  a(2*io,(jo+1)/2,lot)
   real                 ::  b(io,jo)
!-------------------------------------------------------------------------------
!
   do k = 1,lot
     do j = 1,(jo+1)/2
       do i = 1,io
         b(i,j)=a(i,j,k)
         b(i,jo-j+1)=a(io+i,j,k)
       enddo
     enddo
     do ij = 1,io*jo
       a(ij,1,k)=b(ij,1)
     enddo
   enddo
!
   return
!
   end subroutine post_trans_output_grid
