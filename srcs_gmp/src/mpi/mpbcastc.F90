!
   subroutine mpbcastc(n,len)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpbcastc
!            
! abstract: broadcast character array to all pes
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   mpbcastc(n,len)
!
!    input argument lists:
!   n   - character (len) charater array from master pe 
!   len   - integer length of array n
!
!    output argument list:
!   n   - character (len) charater array from master pe 
! 
! subprograms called:
!   mpi_bcast   - to broadcast to all pe in the comm
!
!-------------------------------------------------------------------------------
   use commpi
!-------------------------------------------------------------------------------
   integer           ::  len,ierr
   character(len=1)  ::  n(len)
!
!soojin_couple
!   call mpi_bcast(n,len,mpi_character,0,mpi_comm_world,ierr)
   call mpi_bcast(n,len,mpi_character,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   return
   end subroutine mpbcastc
!-------------------------------------------------------------------------------
