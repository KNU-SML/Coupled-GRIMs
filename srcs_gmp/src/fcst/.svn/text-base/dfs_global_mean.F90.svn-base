#include "define.h"
   SUBROUTINE dfs_global_mean ( AJX,sum,mt,nls,nle,levs,MAKE_0)  
                                                 ! MAKE_0 =+1 : Make 0
!-------------------------------------------------------------------------------
   use dfsvar, only : mls,ngs,nge,iope
!-------------------------------------------------------------------------------
!                                                                    
! abstract :    To make GLOBAL AVERAGE zero ( When MAKE_0 = +1 )               
!                                                                    
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer              ::  mt,nls,nle,levs,MAKE_0
#ifdef MP
   real                 ::  AJX(mt,nls:nle,levs),sum(levs)
#else
   real                 ::  AJX(mt,ngs:nge,levs),sum(levs)
#endif
!
   integer              ::  mj,mzr,k
#ifdef MP
   real                 ::  aj(nls:nle,levs)
#endif
   real                 ::  ajg(ngs:nge,levs)
   real                 ::  dd
!
   if (mls.eq.0) then
     mzr=mt/2+1
#ifdef MP
     if (nls.ne.ngs.or.nle.ne.nge) then
       do k = 1,levs
         aj(nls:nle,k)=AJX(mzr,nls:nle,k)
       enddo
       call mpmn2m(aj,nle-nls+1,ajg,nge-ngs+1,1,levs)
     else
       do k = 1,levs
         ajg(nls:nle,k)=AJX(mzr,nls:nle,k)
       enddo
     endif
#else
     do k = 1,levs
       ajg(ngs:nge,k)=AJX(mzr,ngs:nge,k)
     enddo
#endif
!
!     global mean
!
      do k = 1,levs
        sum(k)= 0.
        do mj = ngs,nge,2
          dd= mj*mj-1.0
          !dd= DFLOAT(mj)**2 - 1.0
          sum(k)= sum(k) - ajg(mj,k)/dd
        enddo
      enddo
   endif
!
#ifdef MP
   if (nls.eq.ngs.and.nle.eq.nge.and.MAKE_0.ne.0) then
     call mpbcastr_col(sum,levs)
   else
     call mpbcastr(sum,levs)
   endif
#endif
!
   if (MAKE_0.eq.+1.and.mls.eq.0) then
     if (nls.eq.0) then
       do k=1,levs
!
!        global-mean subtracted
!
         AJX(mzr,nls,k)= AJX(mzr,nls,k) - sum(k)
       enddo
     endif
   endif
!
   RETURN
   END subroutine dfs_global_mean
