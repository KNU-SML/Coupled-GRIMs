!#include <define.h>
   subroutine dfs_cut_alias(anm,mt,ns,ne,lot)
!-------------------------------------------------------------------------------
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
#ifdef ALIASED
   use dfsvar, only : nmxa,iope,mta
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  mt,ns,ne,lot
   real                 ::  anm(mt,ns:ne,lot)
!
   integer m,l,k,ncs,nend,jtg,mx
!
   if (ne.gt.mta) then
     do k = 1,lot
       do m = 1,mt
         ncs=nmxa(m)
         if (ncs.gt.0) anm(m,ncs:ne,k)=0.0
       enddo
     enddo
   endif
#endif
!
   return
   end subroutine dfs_cut_alias
