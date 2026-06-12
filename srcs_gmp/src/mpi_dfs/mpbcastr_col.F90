#include "define.h"
   subroutine mpbcastr_col(a,len)
!-------------------------------------------------------------------------------
!
!  subprogram documentation block
!
! subprogram:    mpbcastr
!
! usage:   mpbcastr(a,len)
!
!    input argument lists:
!   a   - real (len) array from master pe
!   len   - integer length of array a
!
!    output argument list:
!   a   - real (len) array from master pe
!
! subprograms called:
!   mpi_bcast   - to broadcast to all pe in the comm
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,myrow,comm_col,real_type
!-------------------------------------------------------------------------------
   integer  ::  len,i
   real     ::  a(len)
   real(_mpi_real_),allocatable  ::  tmpa(:)
!
   allocate(tmpa(len))
!
   if( myrow.eq.0 ) then
     do i = 1,len
       tmpa(i) = a(i)
     enddo
   endif
!
   call mpi_bcast(tmpa,len,real_type,0,comm_col,ierr)
!
   if (ierr.ne.0) then
     write(6,*)'mpbcastr_col failed at pe=',mype
     call mpabort
   endif
!
   do i = 1,len
     a(i) = tmpa(i)
   enddo
   deallocate(tmpa)   ! SCC 12/20/05
!
   return
   end subroutine mpbcastr_col
!-------------------------------------------------------------------------------
