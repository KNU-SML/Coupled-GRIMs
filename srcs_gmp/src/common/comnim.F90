#include <define.h>
   module comnim
#ifdef NIM
   use paramodel, only : LATG2S, LONF2S, levs_
!
   real*8, allocatable, dimension(:,:)    :: wor1d
   real*8, allocatable, dimension(:,:,:)  :: wor2d
   real*8, allocatable, dimension(:)      :: wor1d1,wor1d2,wor1d3

   contains

   subroutine comnim_init
   allocate ( wor1d(LONF2S,100)         )
   allocate ( wor2d(LONF2S,levs_,100)   )
   allocate ( wor1d1(LONF2S)            )
   allocate ( wor1d2(LONF2S)            )
   allocate ( wor1d3(LONF2S)            )

   wor1d=0.  ; wor2d=0.
   wor1d1=0. ; wor1d2=0. ; wor1d3=0.
   

   end subroutine comnim_init
#endif
   end module comnim
