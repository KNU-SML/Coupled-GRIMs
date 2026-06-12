!
   subroutine mpsynall
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpsynall
!            
! abstract:  set mpi barrier for all nodes
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpsynall
!
! subprograms called:
!   mpi_barrier  - set barrier for comm
!
!-------------------------------------------------------------------------------
!soojin
!   use commpi, only : mpi_comm_world, mype
   use commpi, only : mpi_comm_private, mype
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ierr
!soojin_couple
!   call mpi_barrier(mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if( ierr.ne.0 ) then
     print *,'PE',mype,':***** Error stop in mpsynall ******** '
     print *,'PE',mype,':error code from mpi_barrier = ',ierr
     call mpabort
   endif
!
   return
   end subroutine mpsynall
!-------------------------------------------------------------------------------
