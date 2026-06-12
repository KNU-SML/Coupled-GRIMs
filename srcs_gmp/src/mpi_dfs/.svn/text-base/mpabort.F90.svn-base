!
   subroutine mpabort
!-------------------------------------------------------------------------------
! subprogram documentation block
!
! subprogram:    mpabort
!            
! abstract: abort the mpi by any pe.
!
! usage:   call mpfine
!
!    input argument lists:
!
!    output argument list:
! 
! subprograms called:
!   mpi_abort   - to abort mpi
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,comm_world
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  info
!
   call mpi_barrier(comm_world,info)
   call mpi_abort(comm_world,1,info)
   if( info.ne.0 ) then
     print *,'PE',mype,': ********* Error stop in mpabort ****** '
     print *,'PE',mype,': error code from mpi_abort =',info
     call exit(1)
   endif
!
   return
   end subroutine mpabort
!-------------------------------------------------------------------------------
