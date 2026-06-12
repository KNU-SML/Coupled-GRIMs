#include "define.h"
   subroutine dyn_dfs_driver
!-------------------------------------------------------------------------------
!
! subroutine: dyn_dfs_driver         
!
! abstract: 
! - computes dynamic non-linear tendency terms of temp. div. ln(ps)
! - omputes predicted values of vorticity and moisture
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   use dfsvar, only : VOR,DIO,TAI,QAI,XPSL,YPSL,APSN,coslat,spdmax,levs       ,&
                      V1,D1,VJ,DJ,TJ,QJ,PSJ                                   ,&
#ifdef HDSZ
                      PSL2                                                    ,&
#endif
#ifdef BIN_DIAG
                      jl,jlg,levsp,ntotal                                     ,&
#endif
#ifdef NISLQ
                      levh                                                    ,&
#endif
#if defined(BIN_DIAG) ||  defined(NISLQ) ||  defined(STOCH)
                      ib,jbw,vxc,vyc                                          ,&
#endif
#if defined(NISLQ) ||  defined(HYBRID)
                      prs                                                     ,&
#endif
                      iope
   use comfver, only : kdt
#ifdef HYBRID
   use comfver, only : ak5,bk5,deltim
#endif
#ifdef BIN_DIAG
   use module_file_write, only : file_write_bin2 
#endif
#ifdef DG3
   use diag_3d_module, only    : diag_3d_get, diag_3d_arrange 
#endif
#ifdef NISLQ
   use nislq         , only    : slq_q2,slq_psfc2,slq_u2,slq_v2,slq_w2
#ifndef HYBRID
   use dfsvar        , only    : sigdot
#endif
#endif /* NISLQ end */
#ifdef STOCH
   use stochastic
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! local array
!
   integer                              ::  i,j,k
#ifdef BIN_DBG
   integer                              ::  nvar,lotg,lotl
#endif
#include "abort.h"
!
! calculate t grid values of vxc,vyc,sigdot,VOR,DIO,TAI,QAI
! and XPSL,YPSL,APSN will be used in adv_3dim
!
   CALL dfs_wave2grid_dynamics
#ifdef BIN_DBG
   call file_write_bin2(334,dio,ib,jbw,(3+ntotal)*levs+1,0)   !  D,V,T,Q,PS
   call file_write_bin2(334,vxc,ib,jbw,2*levs  ,0)   !  U,V
#endif
#ifdef HDSZ
   ! Held-Suarez(1994, BAMS) & Williamson(1998, MWR)
   CALL T_RELAX ( PSL2 ) ! T_REF, T_VIS, Z_VIS
#endif
#ifdef STOCH
!
! update grid-point values for u,v,t,q
!
   do k = 1,levs
     do j = 1,jbw
       do i = 1,ib
         dfsg_u1(i,j,k)=dfsg_u2(i,j,k)
         dfsg_v1(i,j,k)=dfsg_v2(i,j,k)
         dfsg_t1(i,j,k)=dfsg_t2(i,j,k)
         dfsg_q1(i,j,k)=dfsg_q2(i,j,k)
         dfsg_u2(i,j,k)=vxc(i,j,k)
         dfsg_v2(i,j,k)=vyc(i,j,k)
         dfsg_t2(i,j,k)=tai(i,j,k)
         dfsg_q2(i,j,k)=qai(i,j,k)
       enddo
     enddo
   enddo
#endif /* STOCH end */
!
! compute dynamic tendency in VJ,DJ,TJ,QJ,PSJ
!
#ifdef NISLQ
#define QAI slq_q2
#endif
   call dfs_dynamics_advection                                                 &
#ifdef HYBRID
                               (V1,D1,VOR,DIO,TAI,QAI,PRS,XPSL,YPSL,APSN      ,&
                                spdmax,VJ,DJ,TJ,QJ,PSJ,ak5,bk5,kdt            ,&
#ifndef NISLQ
                                deltim)
#else                                 
                                deltim,slq_w2)
#endif
#else
                               (V1,D1,VOR,DIO,TAI,QAI,XPSL,YPSL,APSN,          &
                                spdmax,VJ,DJ,TJ,QJ,PSJ,kdt )
#endif /* HYBRID end */
#ifdef NISLQ
#undef QAI
!
! dyn_grid for nislq
!
   do j = 1,jbw
     do i = 1,ib
       slq_psfc2(i,j)=exp(prs(i,j))  ! Psfc
     enddo
   enddo
!
   do k = 1,levs
     do j = 1,jbw
       do i = 1,ib
         slq_u2(i,j,k)=vxc(i,j,levs+1-k)    ! u-wind
         slq_v2(i,j,k)=vyc(i,j,levs+1-k)    ! v-wind
       enddo
     enddo
   enddo
#ifndef HYBRID
!
   do k = 1,levs+1
     do j = 1,jbw
       do i = 1,ib
         slq_w2(i,j,k)=sigdot(i,j,levs+2-k)*slq_psfc2(i,j) !sigdot:b2t (-1)
       enddo
     enddo
   enddo
#endif
#endif /* NISLQ end */
!
#ifdef BIN_DBG
   nvar=3+ntotal        ! D,V,T, and Q    
   lotl=nvar*levsp+1    ! +PS
   lotg=nvar*levs+1     ! +PS
   call fft_dfs_235_sc(-1,dio,ib,jbw,lotg,dj,mtg,jl,lotl,jlg                  ,&
                      levsp,levs,nvar,coslat,1)
   call file_write_bin2(334,dio,ib,jbw,lotg,0)
   call MPABORT
#endif
   do k = 1,levs
     spdmax(k)=sqrt(spdmax(k))
   enddo
#ifdef MP
   call mpgetspd(spdmax)
#endif
   if (iope) then
     write(6,100)(spdmax(k),k=1,levs)
   endif
100   format(' global checked speed maxima for all layers ',                   &
        :/' spdmx(001:010)=',10f5.0,:/' spdmx(011:020)=',10f5.0,               &
        :/' spdmx(021:030)=',10f5.0,:/' spdmx(031:040)=',10f5.0,               &   
        :/' spdmx(041:050)=',10f5.0,:/' spdmx(051:060)=',10f5.0,               &   
        :/' spdmx(061:070)=',10f5.0,:/' spdmx(071:080)=',10f5.0,               &   
        :/' spdmx(081:090)=',10f5.0,:/' spdmx(091:100)=',10f5.0,               &
        :/' spdmx(101:110)=',10f5.0,:/' spdmx(111:120)=',10f5.0,               &
        :/' spdmx(121:130)=',10f5.0,:/' spdmx(131:140)=',10f5.0,               &               
        :/' spdmx(141:150)=',10f5.0,:/' spdmx(151:160)=',10f5.0,               &
        :/' spdmx(161:170)=',10f5.0,:/' spdmx(171:180)=',10f5.0,               &
        :/' spdmx(181:190)=',10f5.0,:/' spdmx(191:200)=',10f5.0)
!
#ifndef DCMIP
   do k = 1,levs
     if(spdmax(k).eq.0.) then
       if( iope ) then
         write(6,*)'run failure.  spdmax=0.'
         call flush(6)
       endif
       call MPABORT
     endif
   enddo
#endif
!
   return
   end subroutine dyn_dfs_driver
