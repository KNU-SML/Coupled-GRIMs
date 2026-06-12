#include "define.h"
   subroutine mpgetspd(spdmax)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgetspd
!            
! abstract:  get spdmax from all pe
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
   use dfsvar, only : levs
   use commpi, only : comm_world,real_type
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   integer  ::  k,ierr
   real     ::  spdmax(levs)
   real(_mpi_real_),allocatable  ::  spdsnt(:),spdrcv(:)
!
   allocate(spdsnt(levs))
   allocate(spdrcv(levs))
!
   do k = 1,levs
     spdsnt(k)=spdmax(k)
   enddo
!
   call mpi_allreduce(spdsnt,spdrcv,levs,real_type,mpi_max,comm_world,ierr)
!
   do k = 1,levs
     spdmax(k)=spdrcv(k)
   enddo
!
   deallocate(spdrcv)
   deallocate(spdsnt)
!
   return
   end subroutine mpgetspd
!-------------------------------------------------------------------------------
