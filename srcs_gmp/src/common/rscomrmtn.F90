#include <define.h>
   module rscomrmtn
!-------------------------------------------------------------------------------
   use paramodel, only : imn_, jmn_
!-------------------------------------------------------------------------------
#ifndef GTOPO30
   integer, allocatable, dimension(:,:),target    ::  zavg
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomrmtn_init
!-------------------------------------------------------------------------------
#ifndef GTOPO30
   allocate(   zavg    (imn_,jmn_)    )
#endif
!
   end subroutine rscomrmtn_init
!-------------------------------------------------------------------------------
   end module rscomrmtn

