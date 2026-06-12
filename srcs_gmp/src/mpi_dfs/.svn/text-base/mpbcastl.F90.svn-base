!
   subroutine mpbcastl(n,len)
!-------------------------------------------------------------------------------
!
!  subprogram documentation block
!
! subprogram:    mpbcasti
!            
! abstract:  broadcast integer array to all pes
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
   use commpi, only : comm_world
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   integer  ::  len
   logical  ::  n(len)
!
   call mpi_bcast(n,len,mpi_logical,0,comm_world,ierr)
!
   return
   end subroutine mpbcastl
!-------------------------------------------------------------------------------
