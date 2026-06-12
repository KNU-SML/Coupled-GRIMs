#include <define.h>
   subroutine mpgetspd(spdmax)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgetspd
!            
! abstract:  get spdmax from all pe
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpgetspd(spdmax)
!
!    input argument lists:
!   spdmax   - real (levs) array of maximal speed for each pe
!
!    output argument list:
!   spdmax   - real (levs) array of maximal speed from all pe for
!                 master pe.
! 
! subprograms called:
!   mpi_allreduce - to gather message from all pe  to master
!
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!   use commpi, only : MPI_REAL, mpi_max, mpi_comm_world
!soojin_couple
   use commpi, only : MPI_REAL, mpi_max, mpi_comm_private
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                          ::  k,ierr
   real                             ::  spdmax(levs_)
   real(_mpi_real_), allocatable    ::  spdsnt(:), spdrcv(:)
!
   allocate(spdsnt(levs_))
   allocate(spdrcv(levs_))
!
   do k = 1,levs_
     spdsnt(k) = spdmax(k)
   enddo
!soojin_couple
!   call mpi_allreduce(spdsnt,spdrcv,levs_,MPI_REAL,                            &
!                      mpi_max,mpi_comm_world,ierr)
   call mpi_allreduce(spdsnt,spdrcv,levs_,MPI_REAL,                            &
                      mpi_max,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do k = 1,levs_
     spdmax(k) = spdrcv(k)
   enddo
!
   deallocate(spdrcv)
   deallocate(spdsnt)
!
   return
   end subroutine mpgetspd
!-------------------------------------------------------------------------------
