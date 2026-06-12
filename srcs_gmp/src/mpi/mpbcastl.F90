!
   subroutine mpbcastl(n,len)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpbcasti
!            
! abstract:  broadcast integer array to all pes
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
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
!soojin_couple
!   use commpi, only : mpi_comm_world, mpi_logical
   use commpi, only : mpi_comm_private, mpi_logical
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer     ::  len, ierr
   logical     ::  n(len)
!
!soojin_couple
!   call mpi_bcast(n,len,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(n,len,mpi_logical,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   return
   end subroutine mpbcastl
!-------------------------------------------------------------------------------
