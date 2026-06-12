#include <define.h>
   module rsparltb
!-------------------------------------------------------------------------------
#ifdef RMP
! 
!   begin parltb 
!
   use paramodel, only  :  IGRD1S,JGRD12S,bgf_,border_
!-------------------------------------------------------------------------------
   private             ::  IGRD1S,JGRD12S,bgf_,border_
!
   integer             ::  ibgd1, jbgd1, lngrdb
   integer             ::  ib1,ib2,jb1,jb2,jb3,jb4,jbx,istr,ilen
!
   contains
!-------------------------------------------------------------------------------
   subroutine rsparltb_init
!-------------------------------------------------------------------------------
   ibgd1=(IGRD1S-1)/bgf_+2*border_+1
   jbgd1=(JGRD12S-1)*2/bgf_+4*border_+2
   lngrdb=ibgd1*jbgd1
!
   end subroutine rsparltb_init
!-------------------------------------------------------------------------------
#endif
   end module rsparltb
