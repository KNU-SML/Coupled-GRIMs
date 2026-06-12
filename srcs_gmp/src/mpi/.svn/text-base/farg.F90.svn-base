#include "define.h"
   integer function mpir_iargc()
!-------------------------------------------------------------------------------
#ifdef X1E
   mpir_iargc = command_argument_count()
#else
   mpir_iargc = iargc()
#endif
!
   return
   end function mpir_iargc
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine mpir_getarg( i, s )
!-------------------------------------------------------------------------------
   integer       ::  i
   character*(*) ::  s
!
#ifdef X1E
   call get_command_argument(i,s)
#else
   call getarg(i,s)
#endif
!
   return
   end subroutine mpir_getarg
!-------------------------------------------------------------------------------
