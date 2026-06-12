#include "define.h"
   module commpi
!-------------------------------------------------------------------------------
   integer  ::  npes,ncol,nrow,mype,master,msgtag,myrow,mycol,comm_row,       &
                comm_col,comm_world
   integer  ::  real_type,int_type,istat_size
   integer, allocatable, save,dimension(:)  ::                                &
                levstr,levlen,lonstr,lonlen,latstr,latlen,lerlen,lerstr,      &
                lmstr,lmlen,lnstr,lnlen
!
   contains
!-------------------------------------------------------------------------------
   subroutine mpivarini(npes)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  npes,jmax
   integer  ::  ios
!
   allocate(levstr(0:npes-1),levlen(0:npes-1),lonstr(0:npes-1),lonlen(0:npes-1)&
           ,latstr(0:npes-1),latlen(0:npes-1),lerstr(0:npes-1),lerlen(0:npes-1)&
           ,lmstr(0:npes-1),lmlen(0:npes-1),lnstr(0:npes-1),lnlen(0:npes-1)    &
           ,stat=ios)
!
   if (ios /= 0) then
     print*,'mem allocate error: ',ios
     call mpabort
   endif
!
   return
   end subroutine mpivarini
!-------------------------------------------------------------------------------
   end module commpi
