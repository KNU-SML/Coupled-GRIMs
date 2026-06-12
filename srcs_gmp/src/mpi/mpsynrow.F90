!
   subroutine mpsynrow
!-------------------------------------------------------------------------------
! subprogram documentation block
!
! subprogram:    mpsynrow
!            
! abstract: set mpi barrier for comm_row nodes
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpsynrow
!
! subprograms called:
!   mpi_barrier  - set barrier for comm
!
!-------------------------------------------------------------------------------
   use commpi, only : mype, comm_row
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ierr

   call mpi_barrier(comm_row,ierr)
   if( ierr.ne.0 ) then
     print *,'PE',mype,':***** Error stop in mpsynrow ******** '
     print *,'PE',mype,':error code from mpi_barrier = ',ierr
     call mpabort
   endif
!
   return
   end subroutine mpsynrow
!-------------------------------------------------------------------------------
