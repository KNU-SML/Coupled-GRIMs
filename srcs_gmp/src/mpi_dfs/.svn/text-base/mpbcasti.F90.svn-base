!
   subroutine mpbcasti(n,len)
!-------------------------------------------------------------------------------
!
!  subprogram documentation block
!
! subprogram:    mpbcasti
!            
! usage:   mpbcasti(n,len)
!
!    input argument lists:
!   n   - integer (len) array from master pe 
!   len   - integer length of array n
!
!    output argument list:
!   n   - integer (len) array from master pe 
! 
! subprograms called:
!   mpi_bcast   - to broadcast to all pe in the comm
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,master,comm_world,int_type
!-------------------------------------------------------------------------------
   integer  ::  len
   integer  ::  n(len)
!
   call mpi_bcast(n,len,int_type,0,comm_world,ierr)
   if (ierr.ne.0) then
     write(6,*)'mpbcasti failed'
     call mpabort
   endif
!
   return
   end subroutine mpbcasti
!-------------------------------------------------------------------------------
