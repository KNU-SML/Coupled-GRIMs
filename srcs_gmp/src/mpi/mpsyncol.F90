!
   subroutine mpsyncol
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpsyncol
!            
! abstract:  set mpi barrier for comm_column nodes
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpsyncol
!
! subprograms called:
!   mpi_barrier  - set barrier for comm
!
!-------------------------------------------------------------------------------
   use commpi, only : comm_column, mype
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ierr
!
   call mpi_barrier(comm_column,ierr)
   if( ierr.ne.0 ) then
     print *,'PE',mype,':****** Error stop in mpsyncol ****** '
     print *,'PE',mype,':error code from mpi_barrier = ',ierr
     call mpabort
   endif
!
   return
   end subroutine mpsyncol
!-------------------------------------------------------------------------------
