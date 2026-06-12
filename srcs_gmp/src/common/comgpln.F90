#include <define.h>
   module comgpln
!-------------------------------------------------------------------------------
#ifdef MP
   use paramodel, only  : lln2p_,TWOJ1S,latg2_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private             :: lln2p_,TWOJ1S,latg2_
#else
   use paramodel, only  : lnt2_,TWOJ1S,latg2_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private             :: lnt2_,TWOJ1S,latg2_
#endif
!
   real, allocatable  ::  qtt(:,:),qvv(:,:),qdd(:,:),qww(:,:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine comgpln_init
!-------------------------------------------------------------------------------
#ifdef MP
#define LNT2S  lln2p_
#else
#define LNT2S  lnt2_
#endif
   allocate(              qtt(LNT2S,latg2_),qvv(TWOJ1S,latg2_)                ,&
                          qdd(LNT2S,latg2_),qww(LNT2S,latg2_)                  )
!
   end subroutine comgpln_init
!-------------------------------------------------------------------------------
#undef LNT2S
#undef TWOJ1S
   end module comgpln
