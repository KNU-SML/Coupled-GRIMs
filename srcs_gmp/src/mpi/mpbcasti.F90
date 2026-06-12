!
   subroutine mpbcasti(n,len)
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
   use commpi
!-------------------------------------------------------------------------------
   integer     :: len
   integer     :: n(len)
!
!soojin_couple
!   call mpi_bcast(n,len,mpi_integer,0,mpi_comm_world,ierr)
   call mpi_bcast(n,len,mpi_integer,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   return
   end subroutine mpbcasti
!-------------------------------------------------------------------------------
