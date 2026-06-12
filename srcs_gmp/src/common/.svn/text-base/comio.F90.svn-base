#include <define.h>
   module comio
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  LNWAVS,lnt2_
!-------------------------------------------------------------------------------
   private             ::  LNWAVS,lnt2_
#else
   use paramodel, only  :  lnt2_
!-------------------------------------------------------------------------------
   private             ::  lnt2_
#endif
!
#ifdef NIM
   integer             :: ncpu,ndbg1,kdbg1
#endif
   integer             :: ifin,icen,igen,icen2,ienst,iensi,ienss
   real                :: runid,usrid
#ifndef DFS
   integer,allocatable :: ndex(:)
   real, allocatable   :: snnp1(:)
   logical             :: iope
#endif
#ifdef RMP
   real                :: delx, dely, rnnp1max
   real, allocatable   :: rsnnp1(:),rnnp1(:),epsx(:),epsy(:)
#endif
   character(len=32)   :: lab,labs
!
   contains
!-------------------------------------------------------------------------------
   subroutine comio_init
!-------------------------------------------------------------------------------
#ifndef DFS
   allocate(ndex(lnt2_))
   allocate(snnp1(lnt2_))
#endif
#ifdef RMP
   allocate(rsnnp1(LNWAVS),rnnp1(LNWAVS),epsx(LNWAVS),epsy(LNWAVS))
#endif
!
   end subroutine comio_init
!-------------------------------------------------------------------------------
   end module comio
