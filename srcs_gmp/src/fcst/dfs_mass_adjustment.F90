#include "define.h"
   subroutine dfs_mass_adjustment(q,avprs0)
!-------------------------------------------------------------------------------
!
! subroutine: dfs_mass_adjustment
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   use dfsvar, only : ib,jb,jbw,mt,jl,ngs,nge,jlg,coslat,iope
#ifdef MP
   use commpi, only : mype
#endif
!-------------------------------------------------------------------------------
   implicit none
!
   real                 ::  avprs0
   real                 ::  q(mt,jlg)
   real                 ::  prs(ib,jbw)
   real                 ::  avprs
   integer              ::  lat,i
!-------------------------------------------------------------------------------
   CALL dfs_global_mean( q, avprs, mt,ngs,nge,1,0 )
#ifdef DBG
   if (iope) then
     write(6,'(2(A,F15.10))')' avprs=',avprs,' avprs0=',avprs0
   endif
#endif
!
!  pressure correction loop
!
   if(avprs0.gt.0) then
     call dfs_fft_driver(-1,prs,ib,jbw,1,q,mt,jlg,1,jlg,1,1,0,coslat,1)
#ifdef OPENMP
!$omp parallel do private(lat,i)
#endif
     do lat = 1,jbw
       do i = 1,ib
         prs(i,lat)=prs(i,lat)-avprs+avprs0
       enddo
     enddo
     call dfs_fft_driver(1,prs,ib,jbw,1,q,mt,jlg,1,jlg,1,1,0,coslat,1)
   else
     avprs0=avprs
   endif
!
   return
   end subroutine dfs_mass_adjustment
