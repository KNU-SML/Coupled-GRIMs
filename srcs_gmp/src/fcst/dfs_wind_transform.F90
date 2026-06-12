#include <define.h>
   subroutine dfs_wind_transform
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!    [dfs_wind_transform]
!      |
!      |--- [dfs_wave2grid_dynamics] *
!      |--- [dfs_wave2grid_physics] *
!      |--- [dfs_grid2wave_physics] *
!      |--- [dfs_wind2divor] *
!      |--- [dfs_divor2wind] *
!      |--- [dfs_divor2wind_kinetic] *
!      |--- [dfs_psi2wind] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   end subroutine dfs_wind_transform
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_wave2grid_dynamics
!-------------------------------------------------------------------------------
   use dfsvar, only : v2,d2,t2,q2,psl2,xps2,vor,dio,vxc,ypsl                  ,&
#ifndef HYBRID
                      xpsl,apsn,vyc                                           ,&
                      latdef,jls,sigdot,dwR,wR,apsk,divn,adt1,tsigma,delsig   ,&
#endif
                      mt,jlg,levsp,jla,jlha,jgs,levs,ib,jbw,levs,ntotal       ,&
                      amatm,AMATm_,dmatm,DMATm_,amsquar,jcol2js               ,&
                      coslat,mls,tlevs,tlevsp
!-------------------------------------------------------------------------------
!
!     3-DIMENSIONAL velocity field : VXC,VYC,SIGDOT                  
!        dwR,wR,apsk,divn,adt1, XPSL, YPSL, and APSN 
!                  are used in dfs_dynamics_advection                
!       V2, D2, T2, PSL2 are input wave
!       VOR,DIO,TAI, XPSL, YPSL are grid value
!      SIGDOT, VXC,VYC, APSK
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! local-variable
!
   integer                                       ::  k,mx,mj,i,j,mm,jj,mzr
   integer                                       ::  nvar,lot,jb,js,lotg,lotl
   real   , dimension(ib)                        ::  wxy,wxz
   real   , dimension(mt,jlg  ,2*levsp), target  ::  p2, sanm1
   real   , dimension(mt,jlg+1,2*levsp), target  ::             sanmu
   real   , dimension(:,:,:)           , pointer ::      sanm2, sanmv
   real   , dimension(mt,jlg+1,2*levsp+1)        ::  worka
!-------------------------------------------------------------------------------
!
! ***  wave to grid ***
! worka for vxc,vyc,ypsl
!
   jb=jbw/2
!
   if(mls.eq.0) then
     mzr=mt/2+1
     do k = 1,2*levsp
       p2(mzr,:,k)=0.
     enddo
   endif
   sanm2=>sanm1(1:mt,1:jlg  ,levsp+1:2*levsp)
   sanmv=>sanmu(1:mt,1:jlg+1,levsp+1:2*levsp)
!
! v2, d2
!
   CALL dfs_divergence_laplacian1( p2, v2, amatm,AMATm_,dmatm,DMATm_,          &
                     amsquar,jcol2js,mt,jlha,2*levsp,1 )
   CALL dfs_psi2wind( p2, sanmu, sanm1 ,mt, jla,2*levsp )
!
! 1/a(-cos(lat)d/dlat, d/dlon)
!
   do k = 1,levsp
     do mj = 1,jlg
       do mx = 1,mt
         worka(mx,mj,k      )= SANMU(mx,mj,k)+SANM2(mx,mj,k)
         worka(mx,mj,levsp+k)=-SANMV(mx,mj,k)+SANM1(mx,mj,k)
       enddo
     enddo
     do mx = 1,mt
       worka(mx,jlg+1,k      )= SANMU(mx,jlg+1,k)
       worka(mx,jlg+1,levsp+k)=-SANMV(mx,jlg+1,k)
     enddo
   enddo
!
   CALL dfs_psi2wind( psl2, worka(1,1,2*levsp+1), xps2, mt, jla, 1)
!
! (3D) vor,dio,tai,qai (2D) prs,xpsl <- (3D) d3,t3,q3 (2D) psl3,xps3
!
#ifdef NISLQ
   nvar=3         ! vor,dio,tai
#else
   nvar=3+ntotal  ! vor,dio,tai,qai
#endif
   lotg=tlevs+2
   lotl=tlevsp+2
   call dfs_fft_driver(-1,vor,ib,jbw,lotg,v2,mt,jlg,lotl,jlg,                  &
                      levsp,levs,nvar,coslat,1)
!
! (3D) vxc,vyc (2D) ypsl <- (3D) v2,d2 (2D) psl2
!  ... ypsl is multipied by (-1.), previous:ypsl=-coslat*(dp/dlat)
!
   nvar=2         ! vxc, vyc
   lotg=nvar*levs+1
   lotl=nvar*levsp+1
   call dfs_fft_driver(-1,vxc,ib,jbw,lotg,worka,mt,jlg+1,lotl,jlg+1,           &
                      levsp,levs,nvar,coslat,1)
   forall(i=1:ib,j=1:jbw) ypsl(i,j)=-ypsl(i,j)
!
!--beg-( d SIGMA / d T )---------------------------------------------!
!
#ifndef HYBRID
   mm=1
   do J = 1,jb
     jj=latdef(j+jls-jgs)
     js=jbw-j+1
     do I = 1,IB
       wxy(i)=(VXC(i,j ,mm)*xpsL(i,j )+VYC(i,j ,mm)*ypsL(i,j ))*COSLAT(JJ,3)
       dwR(i,j,mm)=(DIO(i,j ,mm)+wxy(i))*delsig(mm)
       wR(i,j,mm)=wxy(i) *delsig(mm)
       APSK(i,j,mm)=wxy(i) 
!
! S.H
!
       wxz(i)=(VXC(i,js,mm)*xpsL(i,js)+VYC(i,js,mm)*ypsL(i,js))*COSLAT(JJ,3)
       dwR(i,js,mm)=(DIO(i,js,mm)+wxz(i))*delsig(mm)
       wR(i,js,mm)=wxz(i) *delsig(mm)
       APSK(i,js,mm)=wxz(i) 
     enddo
   enddo
!
   do mm = 2,levs
     do J = 1,jb
       jj=latdef(j+jls-jgs)
       js=jbw-j+1
       do I = 1,ib
         wxy(i)=(VXC(i,j ,mm)*xpsL(i,j )+VYC(i,j ,mm)*ypsL(i,j ))*COSLAT(jj,3)
         dwR(i,j ,mm)=(DIO(i,j ,mm)+wxy(i))*delsig(mm) + dwR(i,j ,mm-1)
         wR(i,j ,mm)=wxy(i) *delsig(mm) +  wR(i,j ,mm-1)
         APSK(i,j ,mm)=wxy(i)         ! adv-ps-k-level
!
! S.H
!
         wxz(i)=(VXC(i,js,mm)*xpsL(i,js)+VYC(i,js,mm)*ypsL(i,js))*COSLAT(jj,3)
         dwR(i,js,mm)=(DIO(i,js,mm)+wxz(i))*delsig(mm) + dwR(i,js,mm-1)
         wR(i,js,mm)=wxz(i) *delsig(mm) +  wR(i,js,mm-1)
         APSK(i,js,mm)=wxz(i)         ! adv-ps-k-level
       enddo
     enddo
   enddo ! mm=2,levs
!
   do J = 1,jbw
     do I = 1,ib
       DIVN(I,J)=  dwR(i,j,levs)-wR(i,j,levs)       
       APSN(I,J)=  -wR(i,j,levs)                    
     enddo
   enddo
!
!============================================  for level r  ====!
!
   do k = 1,levs-1
     do j = 1,jbw
       do i = 1,ib
         !*(-1.) ! (-1) up K
         sigdot(i,j,k+1)= -(tsigma(k)*dwR(i,j,levs)-dwR(i,j,k))
         ADT1(i,j,k+1)= -(tsigma(k)* wR(i,j,levs)- wR(i,j,k))
       enddo
     enddo
   enddo
!
   do j = 1,jbw
     do i = 1,ib
       sigdot(i,j,  1   )= 0
       ADT1(i,j,  1   )= 0
       sigdot(i,j,levs+1)= 0
       ADT1(i,j,levs+1)= 0
     enddo
   enddo
#endif /* ~HYBRID end */
!
!============================================  for level r  ====!
!
!--end-( d SIGMA / d T )---------------------------------------------!
!
   return
   end subroutine dfs_wave2grid_dynamics
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_wave2grid_physics
!-------------------------------------------------------------------------------
   use dfsvar, only : v3,d3,t3,q3,psl3,xps3,lap3,dio,tai,ypsl,vxc             ,&
                      mt,jlg,levsp,levhp,jla,jlha,levs,ib,jbw,levs,levh,ntotal,&
                      amatm,AMATm_,dmatm,DMATm_,amsquar,jcol2js               ,&
                      coslat,mls,tlevs,tlevsp
!-------------------------------------------------------------------------------
!
! wave to grid transform after physics_main_solver
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                       ::  k,mx,mj,i,j,mm,jj,mzr
   integer                                       ::  nvar,lot,jb,js,lotg,lotl
   real   , dimension(ib)                        ::  wxy,wxz
   real   , dimension(mt,jlg  ,2*levsp), target  ::  p3, sanm1
   real   , dimension(mt,jlg+1,2*levsp), target  ::             sanmu
   real   , dimension(:,:,:)           , pointer ::      sanm2, sanmv
   real   , dimension(mt,jlg+1,2*levsp+1)        ::  worka
!
! ***  wave to grid ***
! worka for vxc,vyc,ypsl
!
   jb=jbw/2
!
   if(mls.eq.0) then
     mzr=mt/2+1
     do k = 1,2*levsp
       p3(mzr,:,k)=0.
     enddo
   endif
   sanm2=>sanm1(1:mt,1:jlg  ,levsp+1:2*levsp)
   sanmv=>sanmu(1:mt,1:jlg+1,levsp+1:2*levsp)
!
! v3,d3 input
!
   CALL dfs_divergence_laplacian1( p3, v3, amatm,AMATm_,dmatm,DMATm_,          &
                     amsquar,jcol2js,mt,jlha,2*levsp,1 )
   CALL dfs_psi2wind( p3, sanmu, sanm1 ,mt, jla,2*levsp )
!
! 1/a(-cos(lat)d/dlat, d/dlon)
!
   do k = 1,levsp
     do mj = 1,jlg
       do mx = 1,mt
         worka(mx,mj,k      )= SANMU(mx,mj,k)+SANM2(mx,mj,k)
         worka(mx,mj,levsp+k)=-SANMV(mx,mj,k)+SANM1(mx,mj,k)
       enddo
     enddo
     do mx = 1,mt
       worka(mx,jlg+1,k      )= SANMU(mx,jlg+1,k)
       worka(mx,jlg+1,levsp+k)=-SANMV(mx,jlg+1,k)
     enddo
   enddo
!
   CALL dfs_psi2wind( psl3, worka(1,1,2*levsp+1), xps3, mt, jla, 1)
   call dfs_divergence_laplacian1 ( psl3, lap3 ,amatm,AMATm_,dmatm,DMATm_,     &
                             amsquar,jcol2js,mt,jlha, 1,0 )
!
! (3D) dio,tai,qai (2D) prs,xpsl,plap <- (3D) d3,t3,q3 (2D) psl3,xps3,lap3
!
#ifdef NISLQ
   nvar=2          ! dio,tai
#else
   nvar=2+ntotal   ! dio,tai,qai
#endif
   lotg=nvar*levs+3
   lotl=nvar*levsp+3
   call dfs_fft_driver(-1,dio,ib,jbw,lotg,d3,mt,jlg,lotl,jlg,                  &
                      levsp,levs,nvar,coslat,1)
!
! (3D) vxc,vyc (2D) ypsl <- (3D) v3,d3 (2D) psl3
!  ... ypsl is multipied by (-1.), previous:ypsl=-coslat*(dp/dlat)
!
   nvar=2          ! vxc,vyc
   lotg=nvar*levs+1
   lotl=nvar*levsp+1
   call dfs_fft_driver(-1,vxc,ib,jbw,lotg,worka,mt,jlg+1,lotl,jlg+1,           &
                      levsp,levs,nvar,coslat,1)
   forall(i=1:ib,j=1:jbw) ypsl(i,j)=-ypsl(i,j)
!
   return
   end subroutine dfs_wave2grid_physics
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dfs_grid2wave_physics
!-------------------------------------------------------------------------------
   use constant, only : rrerth_
   use dfsvar, only   : vxc,vyc,Vj,Dj,Tj,coslat,tai                           ,&
                        ib,jbw,levs,mt,jlg,levsp,ntotal                       ,&
                        jgs,jge,jls,jla,midx,mop,ngs,nge,latdef
!-------------------------------------------------------------------------------
!
! grid to wave transform after physics_main_solver
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer , parameter                               ::  isuvc=1
   integer                                           ::  k,j,i,mj,mx,jj,m     ,&
                                                         jb,js,nvar,lotl,lotg
   real    , target   , dimension(ib,jbw  ,2*levs)   ::  SORK1
   real    , target   , dimension(mt,jlg+1,2*levsp)  ::  SANMU
   real    , target   , dimension(mt,jlg  ,2*levsp)  ::  SANM1
   real    , pointer  , dimension(:,:,:)             ::  SORK2,SANMV,SANM2
   real               , dimension(2*levsp)           ::  xmat0
!
   SORK2=>SORK1(1:ib,1:jbw,1+levs:2*levs)
   jb=jbw/2
!
! virtual wind
!
   do k = 1,levs
     do j = 1,jb
       jj=latdef(jls+j-1)
       js=jbw-j+1
       do i=1,ib
         if (isuvc.ne.0) then
           SORK1(I,J ,k)= VXC(I,J ,k)*coslat(jj,3)              ! 1/cos^2(lat)
           SORK2(I,J ,k)= VYC(I,J ,k)*coslat(jj,3)
           SORK1(I,js,k)= VXC(I,js,k)*coslat(jj,3)              ! 1/cos^2(lat)
           SORK2(I,js,k)= VYC(I,js,k)*coslat(jj,3)
         else
           SORK1(I,J ,k)= VXC(I,J ,k)*coslat(jj,2)              ! 1/cos(lat)
           SORK2(I,J ,k)= VYC(I,J ,k)*coslat(jj,2)
           SORK1(I,js,k)= VXC(I,js,k)*coslat(jj,2)              ! 1/cos(lat)
           SORK2(I,js,k)= VYC(I,js,k)*coslat(jj,2)
         endif
       enddo  ! ib
     enddo    ! jb
   enddo      ! levs
   call dfs_fft_driver(1,SORK1,ib,jbw,2*levs,SANMU,mt,jlg+1,2*levsp,jlg+1,     &
                        levsp,levs,2,coslat,1)
!
! Vj,Dj
!
   SANMV=>SANMU(1:mt,1:jlg+1,1+levsp:2*levsp)
   SANM2=>SANM1(1:mt,1:jlg  ,1+levsp:2*levsp)
   CALL dfs_advection_y( SANMU, SANM1, MT, JLA,2*levsp  )
!
   do k = 1,levsp
     do mj = 1,jlg
       do m = 1,mt      
         mx=midx(m)
         Vj(m,mj,k)= SANMV(mop(m),mj,k)*(-mx)*rrerth_ - SANM1(m,mj,k)
         Dj(m,mj,k)= SANMU(mop(m),mj,k)*(-mx)*rrerth_ + SANM2(m,mj,k)
       enddo
     enddo
   end do
   call dfs_global_mean( Vj,XMAT0, mt,ngs,nge,2*levsp,+1 )
!
! Tj,Qj
!
#ifdef NISLQ
   nvar=1         ! tai
#else
   nvar=1+ntotal  ! tai,qai
#endif
   lotg=nvar*levs
   lotl=nvar*levsp
   call dfs_fft_driver( 1,tai,ib,jbw,lotg,Tj,mt,jlg,lotl,jlg,                  &
                         levsp,levs,nvar,coslat,1)
!
   return
   end subroutine dfs_grid2wave_physics
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_wind2divor(VXC,VYC,ib,jbw,levs,V2,D2,mt,jlg,levsp,           &
                       COSLAT,isuvc)
!-------------------------------------------------------------------------------
   use constant, only : rrerth_
   use dfsvar, only   : jgs,jge,jls,jla,midx,mop,ngs,nge,latdef
!-------------------------------------------------------------------------------
!                                                                    
!     VORTICITY, DIVERGENCE from U V : U=u*COSLAT, V=v*COSLAT        
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer,intent(in)               :: ib,jbw,mt,jlg,levsp,levs,isuvc
   real,intent(in)                  :: vxc(IB,JBw,levs)
   real,intent(in)                  :: vyc(IB,JBw,levs)
   real,intent(out)                 :: V2(mt,jlg,levsp),D2(mt,jlg,levsp)
   real,intent(in)                  :: COSLAT(jgs:jge,3)
!
! local-variable
!
   real,target                      ::  SORK1(IB,JBW,2*levs)
   real,target                      ::  SANMU(mt,jlg+1,2*levsp)
   real,target                      ::  SANM1(mt,jlg  ,2*levsp)
   real,pointer,dimension(:,:,:)    ::  SORK2,SANMV,SANM2
   real                             ::  xmat0(2*levsp)
   integer                          ::  k,j,i,mj,mx,jj,m,jb,js
!
   SORK2=>SORK1(1:ib,1:jbw,1+levs:2*levs)
   jb=jbw/2
!
   do k = 1,levs
     do j = 1,jb
       jj=latdef(jls+j-1)
       js=jbw-j+1
       do i = 1,ib
         if (isuvc.ne.0) then
           SORK1(I,J ,k)= VXC(I,J ,k)*coslat(jj,3)           ! 1/cos^2(lat)
           SORK2(I,J ,k)= VYC(I,J ,k)*coslat(jj,3)
           SORK1(I,js,k)= VXC(I,js,k)*coslat(jj,3)           ! 1/cos^2(lat)
           SORK2(I,js,k)= VYC(I,js,k)*coslat(jj,3)
         else
           SORK1(I,J ,k)= VXC(I,J ,k)*coslat(jj,2)           ! 1/cos(lat)
           SORK2(I,J ,k)= VYC(I,J ,k)*coslat(jj,2)
           SORK1(I,js,k)= VXC(I,js,k)*coslat(jj,2)           ! 1/cos(lat)
           SORK2(I,js,k)= VYC(I,js,k)*coslat(jj,2)
         endif
       enddo        ! ib
     enddo          ! jb
   enddo            ! levs
!
   call dfs_fft_driver(1,SORK1,ib,jbw,2*levs,SANMU,mt,jlg+1,2*levsp,jlg+1,     &
                        levsp,levs,2,coslat,1)
   SANMV=>SANMU(1:mt,1:jlg+1,1+levsp:2*levsp)
   SANM2=>SANM1(1:mt,1:jlg  ,1+levsp:2*levsp)
   CALL dfs_advection_y( SANMU, SANM1, MT, JLA,2*levsp  )
!
   do k = 1,levsp
     do mj = 1,jlg
       do m = 1,MT      
         mx=midx(m)
         V2(m,mj,k)= SANMV(mop(m),mj,k)*(-mx)*rrerth_ - SANM1(m,mj,k)
         D2(m,mj,k)= SANMU(mop(m),mj,k)*(-mx)*rrerth_ + SANM2(m,mj,k)
       enddo
     enddo
   enddo
!
! D2 & V2
!
   CALL dfs_global_mean( V2, XMAT0(1)      , mt,ngs,nge,levsp,+1 )
   CALL dfs_global_mean( D2, XMAT0(1+levsp), mt,ngs,nge,levsp,+1 )
!
   RETURN
   END SUBROUTINE dfs_wind2divor
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dfs_divor2wind(v2,d2,mt,jlg,levsp,                                &
                         AMATm,DMATm,AMATm_,DMATm_,aMSQUAR,jcol2js,            &
                         vo,di,u,v,ib,jbw,coslat)
!-------------------------------------------------------------------------------
   use dfsvar, only : jlha,jla,mls,iope,jbwa,levs
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                       ::  mt,jlg,levsp,ib,jbw
   real                          ::  v2(mt,jlg,levsp),d2(mt,jlg,levsp)
   real                          ::  AMATm (MT,3,0:jlha,2),DMATm (MT,3,0:jlha,2)
   real                          ::  DMATm_(MT,3,0:jlha,2),AMATm_(MT,3,0:jlha,2)
   real                          ::  aMSQUAR(MT)
   integer                       ::  jcol2js(0:jlha,2)
!
! u & v and di & vo allocated as continously
!
   real                          ::  u(ib,jbw,levs), v(ib,jbw,levs)
   real                          ::  vo(ib,jbw,levs), di(ib,jbw,levs)
   real                          ::  coslat(jbwa,3)
!
! local 
!
   real, target                     ::  P2(MT,JLg,2*levsp)
   real, target                     ::  SANM1(MT,JLg,2*levsp)
   real, target                     ::  SANMU(MT,JLg+1,2*levsp)
   real, pointer, dimension(:,:,:)  ::  X2,SANM2,SANMV
   integer                          ::  i,j,k,mj,mx,mzr,lotg,lotl,nvar
!
   if (mls.eq.0) then
     mzr=mt/2+1
     do k = 1,2*levsp
       P2(mzr,:,k)=0.
     enddo
   endif
!
   X2=>   P2   (1:MT,1:JLg  ,1+levsp:2*levsp)
   SANM2=>SANM1(1:MT,1:JLg  ,1+levsp:2*levsp)
   SANMV=>SANMU(1:MT,1:JLg+1,1+levsp:2*levsp)
!
!--beg-( U,V )-------------------------------------------------------!
!  ->Psi : -> P1
!
   CALL dfs_divergence_laplacian1( P2, V2, amatm,AMATm_,dmatm,DMATm_,          &
                     amsquar,jcol2js,mt,jlha,2*levsp,1 )
!
! 1/a(-cos(lat)d/dlat, d/dlon)
!
   CALL dfs_psi2wind( P2, SANMU, SANM1 ,MT, JLA,2*levsp )

   do k = 1,levsp
     do mj = 1,JLG
       do mx = 1,MT
         sanmu(mx,mj,k)= SANMU(mx,mj,k)+SANM2(mx,mj,k)
         sanmv(mx,mj,k)=-SANMV(mx,mj,k)+SANM1(mx,mj,k)
       enddo
     enddo
     do mx = 1,MT
       sanmu(mx,jlg+1,k)= SANMU(mx,jlg+1,k)                
       sanmv(mx,jlg+1,k)=-SANMV(mx,jlg+1,k)                
     enddo
   enddo
!
! di & vo
!
   nvar=2
   lotg=levs*nvar
   lotl=levsp*nvar
   call dfs_fft_driver(-1,vo,ib,jbw,lotg,V2   ,mt,jlg  ,lotl,jlg,              &
                        levsp,levs,nvar,coslat,1)
!
! u & v
!
   call dfs_fft_driver(-1,u ,ib,jbw,lotg,SANMU,mt,jlg+1,lotl,jlg+1,            &
                        levsp,levs,nvar,coslat,1)
   return
   end subroutine dfs_divor2wind
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dfs_divor2wind_kinetic(d2,v2,mt,jlg,levsp,                       &
                        AMATm,DMATm,AMATm_,DMATm_,aMSQUAR,jcol2js,             &
                        di,vo,u,v,ib,jbw,coslat)
!-------------------------------------------------------------------------------
   use dfsvar, only : jlha,jla,mls,jgs,jge,jls,latdef,levs
!-------------------------------------------------------------------------------
   implicit none
!
   integer                       ::  mt,jlg,levsp,ib,jbw
   real, dimension(mt,jlg,levsp) ::  d2,v2
   real                          ::  AMATm (MT,3,0:jlha,2),DMATm (MT,3,0:jlha,2)
   real                          ::  DMATm_(MT,3,0:jlha,2),AMATm_(MT,3,0:jlha,2)
   real                          ::  aMSQUAR(MT)
   integer                       ::  jcol2js(0:jlha,2)
!
! u & v and di & vo allocated as continously
!
   real                          ::  u(ib,jbw,levs),v(ib,jbw,levs)
   real                          ::  di(ib,jbw,levs),vo(ib,jbw,levs)
   real                          ::  coslat(jgs:jge,3)
!
! local 
!
   real, target                     ::  P2(MT,JLg,2*levsp)
   real, target                     ::  SANM1(MT,JLg,2*levsp)
   real, target                     ::  SANMU(MT,JLg+1,2*levsp)
   real, pointer, dimension(:,:,:)  ::  X2,SANM2,SANMV
   integer                          ::  i,j,k,mj,mx,mzr,jj,jb,js
!
   jb=jbw/2
   if (mls.eq.0) then
     mzr=mt/2+1
     do k = 1,2*levsp
       P2(mzr,:,k)=0.
     enddo
   endif
!
   X2   =>P2   (1:MT,1:JLg  ,1+levsp:2*levsp)
   SANM2=>SANM1(1:MT,1:JLg  ,1+levsp:2*levsp)
   SANMV=>SANMU(1:MT,1:JLg+1,1+levsp:2*levsp)
!
!--beg-( U,V )-------------------------------------------------------!
!
! D2 & V2
!
   CALL dfs_divergence_laplacian1( X2, D2, amatm,AMATm_,dmatm,DMATm_,          &
                     amsquar,jcol2js,mt,jlha,levsp,1 ) ! ->Psi : -> P1
   CALL dfs_divergence_laplacian1( P2, V2, amatm,AMATm_,dmatm,DMATm_,          &
                     amsquar,jcol2js,mt,jlha,levsp,1 ) ! ->Psi : -> P1
!
! 1/a(-cos(lat)d/dlat, d/dlon)
!
   CALL dfs_psi2wind( P2, SANMU, SANM1 ,MT, JLA,2*levsp )
!
   do k = 1,levsp
     do mj = 1,JLg
       do mx = 1,MT
         sanmu(mx,mj,k)= SANMU(mx,mj,k)+SANM2(mx,mj,k)
         sanmv(mx,mj,k)=-SANMV(mx,mj,k)+SANM1(mx,mj,k)
       enddo
     enddo
     do mx = 1,MT
       sanmu(mx,jlg+1,k)= SANMU(mx,jlg+1,k)                
       sanmv(mx,jlg+1,k)=-SANMV(mx,jlg+1,k)                
     enddo
   enddo
!
   call dfs_fft_driver(-1,di,ib,jbw,2*levs,D2,mt,jlg,2*levsp,jlg,              &
                        levsp,levs,2,coslat,1)
   call dfs_fft_driver(-1,u,ib,jbw,2*levs,SANMU,mt,jlg+1,2*levsp,jlg+1,        &
                        levsp,levs,2,coslat,1)
!
   do k = 1,levs
     do j = 1,jb
       jj=latdef(jls+j-1)
       js=jbw-j+1
       do i = 1,ib
         u(i,j ,k)=u(i,j ,k)*coslat(jj,2)
         v(i,j ,k)=v(i,j ,k)*coslat(jj,2)
         u(i,js,k)=u(i,js,k)*coslat(jj,2)
         v(i,js,k)=v(i,js,k)*coslat(jj,2)
       enddo
     enddo
   enddo
!
   return
   end subroutine dfs_divor2wind_kinetic
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_psi2wind                                                     &
                                                         ! SANMU : - of Y-diff.
                      ( P2, SANMU,SANM2, MT, JLA,lot)    ! SANM2 : + of X-diff.
!-------------------------------------------------------------------------------
   use constant, only : rrerth_
   use dfsvar, only   : nod,nev,mods,mevs,midx,mls,mop
!-------------------------------------------------------------------------------
!                                                                    
!     U*COS(LAT), V*COS(LAT) FROM STREAMFUNCTION                     
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jla,lot
   real     ::  P2(MT,0:JLA,lot)
   real     ::  SANMU(MT,0:(JLA+1),lot),SANM2(MT,0:JLA,lot)
   integer  ::  mj,mjp1,mjm1,mx,k,m
!
   do k = 1,lot   
     SANMU(1:MT,0,k)=0.
     SANM2(1:MT,0,k)=0.
!
!-psi:m=1,3,...  U : y-dif
!  u*COS(Lat)=
!
     do mj = 1,JLA-1                          ! -(Y-dif OF PSI)*COS(Lat)
       mjp1=mj+1
       mjm1=mj-1
       Do m = 1,nod
         mx=mods(m)
         SANMU(mx,mj,k)=-((mjm1)*P2(mx,mjm1,k)-(mjp1)*P2(mx,mjp1,k))*0.5
       ENDDo ! par
     enddo
!
! +1 level, P2(JLA+1,mx)= 0
!
     mj=JLA
     mjm1=mj-1
     Do m = 1,nod
       mx=mods(m)
       SANMU(mx,mj,k)= -(mjm1)*P2(mx,mjm1,k) *0.5
     ENDDo ! par
!
! +1 level, P2(JLA+2,mx)= 0
!
     mj=JLA+1
     mjm1=mj-1
     Do m = 1,nod
       mx=mods(m)
       SANMU(mx,mj,k)= -(mjm1)*P2(mx,mjm1,k) *0.5
     ENDDo ! par
!
!      V : x-dif,  X-dif of PSI =  v*COS(LAT)
!
     do mj = 1,JLA                            !!
       Do m = 1,nod
         mx=mods(m)
         SANM2(mx,mj,k)=  P2(mop(mx),mj,k)*(-midx(mx))
       ENDDo ! par
     enddo
!
!-psi:m=2,4,...  u : y-dif                  
!
     do mj = 1,JLA-1                      !  u=
       mjp1=mj+1                          ! -(Y-dif OF ([PSI^*]*[COSLAT]))
       mjm1=mj-1
       Do m = 1,nev
         mx=mevs(m)
         SANMU(mx,mj,k)= -(mj  )*(P2(mx,mjm1,k)-P2(mx,mjp1,k))*0.5
       ENDDo ! par
     enddo
!
     mj=JLA
     mjm1=mj-1
     Do m = 1,nev
       mx=mevs(m)
       SANMU(mx,mj,k)= -(mj  )* P2(mx,mjm1,k) *0.5   !! +1 level
     ENDDo ! par
!
     mj=JLA+1
     mjm1=mj-1
     Do m = 1,nev
       mx=mevs(m)
       SANMU(mx,mj,k)= -(mj  )* P2(mx,mjm1,k) *0.5   !! +1 level
     ENDDo ! par
!
!    V : x-dif 
!
     do mj = 0,JLA                            
       Do m = 1,nev
         mx=mevs(m)
         SANM2(mx,mj,k)=  P2(mop(mx),mj,k)*(-midx(mx))    ! X-dif of [PSI^*] 
       ENDDo ! par
     enddo
!
!-psi:m=0        U : y-dif
!
     if (mls.eq.0) then
       mx=mt/2+1                               !  u*COS(Lat)=
       do mj = 1,JLA-1                         ! -(Y-dif OF PSI)*COS(Lat)
         SANMU(mx,mj,k)=-((mj-1)*P2(mx,mj-1,k)-                                &
                         (mj+1)*P2(mx,mj+1,k))*0.5
       enddo                                     
!
! +1 level, P2(JLA+1,mx)= 0
!
       mj=JLA
       SANMU(mx,mj,k)= -(mj-1)*P2(mx,mj-1,k) *0.5
! +1 level, P2(JLA+2,mx)= 0
       mj=JLA+1
       SANMU(mx,mj,k)= -(mj-1)*P2(mx,mj-1,k) *0.5
       SANMU(mx,0 ,k)=  ( 0+1)*P2(mx,0+1 ,k) *0.5 
     endif
!
!                V : x-dif = 0
!
! [1/radius] Dimensionalize by Earth's parameter
! (-) for N.P starting
!
     do mj = 0,JLA
       do mx = 1,mt
         SANMU(mx,mj,k)= -SANMU(mx,mj,k)*rrerth_
         SANM2(mx,mj,k)=  SANM2(mx,mj,k)*rrerth_
       enddo
     enddo
     mj=JLA+1
     do mx = 1,mt
       SANMU(mx,mj,k)= -SANMU(mx,mj,k)*rrerth_
     enddo
   enddo      ! k
!
   RETURN
   END SUBROUTINE dfs_psi2wind
!-------------------------------------------------------------------------------
