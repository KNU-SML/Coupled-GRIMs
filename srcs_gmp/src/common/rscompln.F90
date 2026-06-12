#include <define.h>
   module rscompln
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  jwav1_,jgrd12_
!-------------------------------------------------------------------------------
   private             ::  jwav1_,jgrd12_
!
! do not consider MP dimension here 
!
   real   , allocatable, dimension(:,:)   ::  ccosg                           ,&
                                              csing                           ,&
                                              gcosc                           ,&
                                              gsinc
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscompln_init
!-------------------------------------------------------------------------------
   allocate(    ccosg  (jwav1_,jgrd12_)                                       ,&
                csing  (jwav1_,jgrd12_)                                       ,&
                gcosc  (jwav1_,jgrd12_)                                       ,&
                gsinc  (jwav1_,jgrd12_)                                        )
!
   end subroutine rscompln_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscompln
