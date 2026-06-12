#include "define.h"
   subroutine dfs_wind_damping(div,vor,tem,q3,deltim,spdmax)
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
   use constant, only : rerth_
   use dfsvar, only : mt,jl,levs,levh,ntotal,nls,midx,iope
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                          ::  deltim
   real, dimension(mt,jl,levs)  ::  div,vor,tem
   real, dimension(mt,jl,levh)  ::  q3
   real                          ::  spdmax(levs)
!
! local
!
   integer                       ::  m,n,k,mx,jj,kk,it
   real                          ::  cinv,rncrit,coef
!
! local save
!
   real, allocatable, save       ::  arn(:,:)
   real, save                    ::  alfa,beta,alfadt
   logical                       ::  first
   data first/.true./
!
   if (first) then
     first=.false.
     allocate(arn(mt,jl))
     do n = 1,jl
       jj=nls+n-1
       do m = 1,mt
         mx=abs(midx(m))
         arn(m,n)=sqrt((mx+jj)*(mx+jj+1)+0.25e0)-0.5D0
       enddo
     enddo
     alfa=2.5e0
     beta=rerth_*1.009e0/deltim
     alfadt=alfa*deltim/rerth_
   endif
!
   do k = 1,levs
     rncrit=beta/spdmax(k)
     coef=alfadt*spdmax(k)
     do n = 1,jl
       do m = 1,mt
         if (arn(m,n).gt.rncrit.and.midx(m).ne.0) then
           cinv=1./(1.+(arn(m,n)-rncrit)*coef)
           div(m,n,k) =div(m,n,k)*cinv
           vor(m,n,k) =vor(m,n,k)*cinv
           tem(m,n,k) =tem(m,n,k)*cinv
           do it = 1,ntotal
             kk=(it-1)*levs+k
             tem(m,n,kk) =tem(m,n,kk)*cinv
           enddo
         endif
       enddo
     enddo
     do it = 1,ntotal
       kk=(it-1)*levs+k
       do n = 1,jl
         do m = 1,mt
           if (arn(m,n).gt.rncrit.and.midx(m).ne.0) then
             cinv=1./(1.+(arn(m,n)-rncrit)*coef)
             q3(m,n,kk) =q3(m,n,kk)*cinv
           endif
         enddo
       enddo
     enddo
   enddo         ! k=1,levs
!
   return
   end subroutine dfs_wind_damping
