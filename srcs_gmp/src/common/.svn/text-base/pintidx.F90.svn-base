#include <define.h>
   module pintidx
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  IGRD12S,border_
!-------------------------------------------------------------------------------
   private             ::  IGRD12S,border_
!
   integer             ::  ncntmax
   integer             ::  ncount
   integer, allocatable  ::  irra(:,:), ibba(:,:),        &
                             iwa (:,:), ina (:,:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine pintidx_init
!-------------------------------------------------------------------------------
   ncntmax=border_*2
   allocate( irra(IGRD12S,ncntmax), ibba(IGRD12S,ncntmax),&
             iwa (IGRD12S,ncntmax), ina (IGRD12S,ncntmax) )
!
   end subroutine pintidx_init
!-------------------------------------------------------------------------------
#endif
   end module pintidx
