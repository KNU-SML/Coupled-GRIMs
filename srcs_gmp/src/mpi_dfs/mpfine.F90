!
   subroutine mpfine(endwtime)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpfine
!            
! abstract: finalizing the mpi by each pe.
!
! usage:   call mpfine
!
!    input argument lists:
!
!    output argument list:
! 
! subprograms called:
!   mpi_finalize   - to end of mpi
!
!-------------------------------------------------------------------------------
   use commpi, only : mype
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   real*8   ::  endwtime
   integer  ::  info
!
   endwtime=mpi_wtime()
   call mpi_barrier(mpi_comm_world,info)
   call mpi_finalize(info)
!
   if( info.ne.0 ) then
      print *,'PE',mype,': ********* Error stop in mpfine ******* '
      print *,'PE',mype,': error code from mpi_finalize =',info
      stop 
   endif
!
   return
   end subroutine mpfine
!-------------------------------------------------------------------------------
