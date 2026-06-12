#include "define.h"
   subroutine mpsumrow(a,ij)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgetspd
!            
! abstract:  get sum from row comm.
!
! usage:   call mpsumrow(spdmax)
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
   use commpi, only : comm_row,real_type
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   integer                       ::  ij,k,ierr
   real                          ::  a(ij)
   real(_mpi_real_),allocatable  ::  snt(:),rcv(:)
!
   allocate(snt(ij))
   allocate(rcv(ij))
!
   do k = 1,ij
     snt(k)=a(k)
   enddo
!
   call mpi_allreduce(snt,rcv,ij,real_type,mpi_sum,comm_row,ierr)
!
   do k = 1,ij
     a(k)=rcv(k)
   enddo
!
   deallocate(rcv)
   deallocate(snt)
!
   return
   end subroutine mpsumrow
!-------------------------------------------------------------------------------
