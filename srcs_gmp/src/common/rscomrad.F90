#include <define.h>
   module rscomrad
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  IGRD12S,JGRD12S
!-------------------------------------------------------------------------------
   private             ::  IGRD12S,JGRD12S
!
   real   , allocatable, dimension(:,:)    ::  cvavg
   real   , allocatable, dimension(:,:,:)  ::  fluxr
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomrad_init
!-------------------------------------------------------------------------------
   allocate(   cvavg     (IGRD12S,JGRD12S)                                    ,&
               fluxr     (IGRD12S,JGRD12S,26)                                  )
!
   end subroutine rscomrad_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscomrad
