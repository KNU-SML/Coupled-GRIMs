#include <define.h>
   module comreduce
!-------------------------------------------------------------------------------
   use paramodel, only  : latg_
#ifdef MP
   use paramodel, only  : npes_,latg2p_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private              :: latg_
#ifdef MP
   private              :: npes_,latg2p_
#endif
!
   integer, allocatable :: lcapd (:)
#ifdef MP
   integer, allocatable :: lonfd (:)
   integer, allocatable :: lcapdp(:,:)
   integer, allocatable :: lonfdp(:,:)
   integer, allocatable :: lonfds(:,:)
#else
   integer, allocatable :: lonfd (:)
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine comreduce_init
!-------------------------------------------------------------------------------
   allocate(lcapd (latg_/2))
#ifdef MP
   allocate(lonfd (latg_/2))
   allocate(lcapdp(latg_/2,0:npes_-1))
   allocate(lonfdp(latg2p_,0:npes_-1))
   allocate(lonfds(latg2p_,0:npes_-1))
#else
   allocate(lonfd (latg_/2))
#endif
   end subroutine comreduce_init
!-------------------------------------------------------------------------------
   end module comreduce
