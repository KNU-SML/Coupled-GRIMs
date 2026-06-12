#include <define.h>
   subroutine mpfine(endwtime)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpfine
!            
! abstract: finalizing the mpi by each pe.
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
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
   real*8               ::  endwtime, mpi_wtime
   integer              ::  info
!
   endwtime = mpi_wtime()
#ifdef MPFINE
   call mpi_finalize(info)
   if( info.ne.0 ) then
      print *,'PE',mype,': ********* Error stop in mpfine ******* '
      print *,'PE',mype,': error code from mpi_finalize =',info
      call mpabort 
   endif
#endif
!
   return
   end subroutine mpfine
