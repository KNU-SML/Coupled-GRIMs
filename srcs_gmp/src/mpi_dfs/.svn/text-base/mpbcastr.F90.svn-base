#include "define.h"
   subroutine mpbcastr(a,len)
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
   use commpi, only : mype,master,comm_world,real_type
!-------------------------------------------------------------------------------
   integer  ::  len,i
   real     ::  a(len)
!
   real(_mpi_real_),allocatable  ::  tmpa(:)
!
   allocate(tmpa(len))
!
   if( mype.eq.master ) then
     do i = 1,len
       tmpa(i) = a(i)
     enddo
   endif
!
   call mpi_bcast(tmpa,len,real_type,0,comm_world,ierr)
!
   if (ierr.ne.0) then
     write(6,*)'mpbcastr failed'
     call mpabort
   endif
!
   do i = 1,len
     a(i) = tmpa(i)
   enddo
   deallocate(tmpa) ! SCC 12/20/05
!
   return
   end subroutine mpbcastr
!-------------------------------------------------------------------------------
