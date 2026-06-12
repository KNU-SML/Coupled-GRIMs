!
   subroutine mpbcastc(n,len)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
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
   use commpi, only : mype,comm_world
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   include 'mpif.h'
   integer  ::  len,ierr
   character*1 n(len)
!
   call mpi_bcast(n,len,mpi_character,0,comm_world,ierr)
   if (ierr.ne.0) then
     write(6,*)'mpbcastc failed'
     call mpabort
   endif
!
   return
   end subroutine mpbcastc
!-------------------------------------------------------------------------------
