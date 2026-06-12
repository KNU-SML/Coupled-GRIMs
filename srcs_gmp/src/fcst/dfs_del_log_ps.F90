#include "define.h"
   subroutine dfs_del_log_ps(psl,mt,jlg,ps,xps,yps,ib,jbw,coslat)
!-------------------------------------------------------------------------------
   use dfsvar, only : jla,jgs,jge
!-------------------------------------------------------------------------------
   implicit none
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
   integer                        ::  mt,jlg,ib,jbw
   real   , dimension(mt,jlg)     ::  psl
   real   , dimension(ib,jbw)     ::  ps,xps,yps
   real   , dimension(jgs:jge,3)  ::  coslat
!
! local 
!
   integer                        ::  i,j
   real   , dimension(mt,jlg)     ::  sanmu
   real   , dimension(mt,jlg+1)   ::  sanmv
!-------------------------------------------------------------------------------
   sanmv=0.
   sanmu=0.
!
   call dfs_psi2wind( psl,sanmv,sanmu,mt,jla,1 )
   call dfs_fft_driver(-1, ps,ib,jbw,1,psl  ,mt,jlg  ,1,jlg  ,1,1,0,coslat,1)
   call dfs_fft_driver(-1,xps,ib,jbw,1,sanmu,mt,jlg  ,1,jlg  ,1,1,0,coslat,1)
   call dfs_fft_driver(-1,yps,ib,jbw,1,sanmv,mt,jlg+1,1,jlg+1,1,1,0,coslat,1)
!
   forall(i=1:ib,j=1:jbw) yps(i,j)=-yps(i,j)
!     
   return
   end subroutine dfs_del_log_ps
