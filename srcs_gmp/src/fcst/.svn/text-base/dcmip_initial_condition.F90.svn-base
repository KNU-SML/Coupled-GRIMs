#include <define.h>
   subroutine dcmip_initial_condition
!-------------------------------------------------------------------------------
!
!  [Abstract] 
!    To set up the initial condition for DCMIP in GRIMs
!
!  [History]
!    2014.07.01 : Test case 1-1, 1-2    SPH
!    2014.07.11 :                       DFS (also semi-Lagrangian)
!    2014.07.14 : Test case 1-3
!
!-------------------------------------------------------------------------------
   end subroutine dcmip_initial_condition
!-------------------------------------------------------------------------------
#ifdef DCMIP
!
!-------------------------------------------------------------------------------
#ifdef NISLQ
   subroutine dcmip_update(utime,gz,q,di,ze,te,rq,z00,slq)
#else
   subroutine dcmip_update(utime,gz,q,di,ze,te,rq,z00)
#endif
!-------------------------------------------------------------------------------
   use paramodel  , only : levs_,levh_,ntotal_
   use constant   , only : fv_,rerth_
   use dcmip_grims, only : icase,cfv,grav,hybrid_eta,hyam,hybm                ,&
                           pl,zl,u,v,us,vs,pwl,t,tv,phis,ps,rho               ,&
                           pi,zi,          pwi      ,hyai,hybi                ,&
                           q0,q1,q2,q3,q4,psfc,hsfc,ptop,p0
#ifndef DFS
   use dcmip_grims, only : gphf
#endif
   use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation,&
                                                   test1_advection_hadley     ,&
                                                   test1_advection_orography  ,&
                                                   test2_steady_state_mountain,&
                                                   test3_gravity_wave
   use dcmip_initial_conditions_test_4,     only : test4_baroclinic_wave
   use dcmip_initial_conditions_test_5,     only : test5_tropical_cyclone
#ifdef DFS
   use dfsvar     , only : iope,mt,jlg,ib,jbw,levsp,levhp,latdef,jls,coslat,xlon,xlat
#else /* SPH */
   use paramodel  , only : LONF2S,LATG2S,LNT22S,jcap_,lnt2_,lnt22_
   use comfgrid   , only : rbs2
   use comfphys   , only : xlon,xlat
   use comio      , only : iope,snnp1
#ifdef MP
   use paramodel  , only : lnt22p_
   use commpi     , only : lwvdef
#endif
#endif /* DFS end */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   real                                   , intent(in)      ::  utime
#ifdef DFS
   real   , dimension(mt,jlg)             , intent(out)     ::  gz,q
   real   , dimension(mt,jlg,levsp)       , intent(out)     ::  ze,di,te
   real   , dimension(mt,jlg,levhp)       , intent(out)     ::  rq
#ifdef NISLQ
   real   , dimension(ib,jbw,levh_)       , intent(out)     ::  slq
#endif
#else /* SPH end */
   real   , dimension(LNT22S)             , intent(out)     ::  gz,q
   real   , dimension(LNT22S,levs_)       , intent(out)     ::  ze,di,te
   real   , dimension(LNT22S,levh_)       , intent(out)     ::  rq
#ifdef NISLQ
   real   , dimension(LONF2S,levh_,LATG2S), intent(out)     ::  slq
#endif
#endif /* DFS end */
   real                                   , intent(out)     ::  z00
!
! local variables
!
   integer                                                  ::  i,j,k,moist
   real                                                     ::  dum
#ifdef DFS
   integer                                                  ::  jj
   real   , dimension(ib,jbw,levh_)                         ::  qq
#else /* SPH */
   real   , dimension(LONF2S,levh_,LATG2S)                  ::  qq
   real   , dimension(lnt22_)                               ::  spec
   real   , dimension(LNT22S)                               ::  gph
#endif /* DFS end */
!-------------------------------------------------------------------------------
#ifdef DFS
#define LONF2S ib
#define LATG2S jbw
#define IJK i,j,k
#define IJKm1 i,j,k-1
#else
#define IJK i,k,j
#define IJKm1 i,k-1,j
#endif
!
! initialize 
!
   phis=0. ; ps=0. ; rho=0. ; u=0. ;  v=0. ; us=0. ; vs=0. ; pwl=0.
   t=0. ; q0=0. ;  q1=0. ;  q2=0. ;  q3=0. ;  q4=0.
!
! test cases
!
   if (icase.eq.11) then
     do j = 1,LATG2S
       do i = 1,LONF2S
         do k = 1,levs_
           call test1_advection_deformation                                    &
                               ( utime,xlon(i,j),xlat(i,j),pl(IJK),zl(IJK),0  ,&
                                 u(IJK),v(IJK),pwl(IJK),t(IJK)                ,&
                                 phis(i,j),ps(i,j),rho(IJK)                   ,&
                                 q0(IJK),q1(IJK),q2(IJK),q3(IJK),q4(IJK) )
         enddo
         do k = 1,levs_+1
           call test1_advection_deformation                                    &
                               ( utime,xlon(i,j),xlat(i,j),pi(IJK),zi(IJK),0  ,&
                                 dum,dum,pwi(IJK),dum                         ,&
                                 dum,dum,dum                                  ,&
                                 dum,dum,dum,dum,dum )
         enddo
       enddo
     enddo
   else if (icase.eq.12) then
     do j = 1,LATG2S
       do i = 1,LONF2S
         do k = 1,levs_
           call test1_advection_hadley                                         &
                               ( utime,xlon(i,j),xlat(i,j),pl(IJK),zl(IJK),0  ,&
                                 u(IJK),v(IJK),pwl(IJK),t(IJK)                ,&
                                 phis(i,j),ps(i,j),rho(IJK)                   ,&
                                 q0(IJK),q1(IJK) )
         enddo
         do k = 1,levs_+1
           call test1_advection_hadley                                         &
                               ( utime,xlon(i,j),xlat(i,j),pi(IJK),zi(IJK),0  ,&
                                 dum,dum,pwi(IJK),dum                         ,&
                                 dum,dum,dum                                  ,&
                                 dum,dum )         
         enddo
       enddo
     enddo
   else if (icase.eq.13) then
     do j = 1,LATG2S
       do i = 1,LONF2S
         do k = 1,levs_
           call test1_advection_orography                                      &
                               (       xlon(i,j),xlat(i,j),pl(IJK),zl(IJK),0  ,&
                                 cfv,hybrid_eta,hyam(k),hybm(k),dum           ,&
                                 u(IJK),v(IJK),pwl(IJK),t(IJK)                ,&
                                 phis(i,j),ps(i,j),rho(IJK)                   ,&
                                 q0(IJK),q1(IJK),q2(IJK),q3(IJK),q4(IJK) )
         enddo
         do k = 1,levs_+1
           call test1_advection_orography                                      &
                               (       xlon(i,j),xlat(i,j),pi(IJK),zi(IJK),0  ,&
                                 cfv,hybrid_eta,hyai(k),hybi(k),dum           ,&
                                 dum,dum,pwi(IJK),dum                         ,&
                                 dum,dum,dum                                  ,&
                                 dum,dum,dum,dum,dum )                    
         enddo
       enddo
     enddo
   else if (icase.eq.200) then
     do j = 1,LATG2S
       do i = 1,LONF2S
         do k = 1,levs_
           call test2_steady_state_mountain                                    &
                               (       xlon(i,j),xlat(i,j),pl(IJK),zl(IJK),0  ,&
                                 hybrid_eta,hyam(k),hybm(k)                   ,&
                                 u(IJK),v(IJK),pwl(IJK),t(IJK)                ,&
                                 phis(i,j),ps(i,j),rho(IJK)                   ,&
                                 q0(IJK) )
         enddo
         do k = 1,levs_+1
           call test2_steady_state_mountain                                    &
                               (       xlon(i,j),xlat(i,j),pi(IJK),zi(IJK),0  ,&
                                 hybrid_eta,hyai(k),hybi(k)                   ,&
                                 dum,dum,pwi(IJK),dum                         ,&
                                 dum,dum,dum                                  ,&
                                 dum )
         enddo
       enddo
     enddo
   else if (icase.eq.410 .or. icase.eq.42 .or. icase.eq.43) then
     if (icase.eq.410) then
       moist=0
     else
       moist=1
     endif
     do j = 1,LATG2S
       do i = 1,LONF2S
         do k = 1,levs_
           call test4_baroclinic_wave                                          &
                            ( moist,1.,xlon(i,j),xlat(i,j),pl(IJK),zl(IJK),0  ,&
                              u(IJK),v(IJK),pwl(IJK),t(IJK)                   ,&
                              phis(i,j),ps(i,j),rho(IJK)                      ,&
                              q0(IJK),q1(IJK),q2(IJK) )
         enddo
         do k = 1,levs_+1
           call test4_baroclinic_wave                                          &
                            ( moist,1.,xlon(i,j),xlat(i,j),pi(IJK),zi(IJK),0  ,&
                              dum,dum,pwi(IJK),dum                            ,&
                              dum,dum,dum                                     ,&
                              dum,dum,dum )
         enddo
       enddo
     enddo
   else if (icase.eq.51) then
     do j = 1,LATG2S
       do i = 1,LONF2S
         do k = 1,levs_
           call test5_tropical_cyclone                                         &
                               (       xlon(i,j),xlat(i,j),pl(IJK),zl(IJK),0  ,&
                                 hybrid_eta,hyam(k),hybm(k)                   ,&
                                 u(IJK),v(IJK),pwl(IJK),t(IJK)                ,&
                                 phis(i,j),ps(i,j),rho(IJK)                   ,&
                                 q0(IJK) )
         enddo
         do k = 1,levs_+1
           call test5_tropical_cyclone                                         &
                               (       xlon(i,j),xlat(i,j),pi(IJK),zi(IJK),0  ,&
                                 hybrid_eta,hyai(k),hybi(k)                   ,&
                                 dum,dum,pwi(IJK),dum                         ,&
                                 dum,dum,dum                                  ,&
                                 dum )
         enddo
       enddo
     enddo
   endif
#undef IJK
#undef IJKm1
!
! 3D
!
   do j = 1,LATG2S
#ifdef DFS
     if(j.le.jbw/2) then
       jj=latdef(jls+j-1)
     else
       jj=latdef(jls+jbw-j)
     endif
#endif
     do k = 1,levs_
       do i = 1,LONF2S
#ifdef DFS
         us(i,j,k)        =u(i,j,k)*coslat(jj,1)         ! u
         vs(i,j,k)        =v(i,j,k)*coslat(jj,1)         ! v
         tv(i,j,k)        =t(i,j,k)*(1.+fv_*q0(i,j,k))-300.   ! t
         qq(i,j,        k)=q0(i,j,k)                     ! spfh
         qq(i,j,  levs_+k)=0.                            ! o3
         qq(i,j,2*levs_+k)=0.                            ! tke
         qq(i,j,3*levs_+k)=q1(i,j,k)                     ! q1
         qq(i,j,4*levs_+k)=q2(i,j,k)                     ! q2
         qq(i,j,5*levs_+k)=q3(i,j,k)                     ! q3
         qq(i,j,6*levs_+k)=q4(i,j,k)                     ! q4
#else
         us(i,k,j)        =u(i,k,j)*sqrt(rbs2(j))        ! u
         vs(i,k,j)        =v(i,k,j)*sqrt(rbs2(j))        ! v
         tv(i,k,j)        =t(i,k,j)*(1.+fv_*q0(i,k,j))   ! t
         qq(i,        k,j)=q0(i,k,j)                     ! spfh
         qq(i,  levs_+k,j)=0.                            ! o3
         qq(i,2*levs_+k,j)=0.                            ! tke
         qq(i,3*levs_+k,j)=q1(i,k,j)                     ! q1
         qq(i,4*levs_+k,j)=q2(i,k,j)                     ! q2
         qq(i,5*levs_+k,j)=q3(i,k,j)                     ! q3
         qq(i,6*levs_+k,j)=q4(i,k,j)                     ! q4
#endif
       enddo
     enddo
   enddo
!
! 2D
!
   do j = 1,LATG2S
     do i = 1,LONF2S
       psfc(i,j)=log(ps(i,j)*1.e-3)  ! Pa -> cb (log pressure)
       hsfc(i,j)=phis(i,j)/grav      ! geopotential -> geopotential height
     enddo
   enddo
#ifdef NISLQ
   slq=qq
#endif
!
! wave to grid
!
#ifdef DFS
#undef LONF2S
#undef LATG2S
   call dfs_wind2divor(us,vs,ib,jbw,levs_,ze,di,mt,jlg,levsp,coslat,1)                       ! ze,di
   call dfs_fft_driver(1,tv  ,ib,jbw,levs_,te,mt,jlg,levsp,jlg,levsp,levs_,1,coslat,1)       ! te
   call dfs_fft_driver(1,qq  ,ib,jbw,levh_,rq,mt,jlg,levhp,jlg,levsp,levs_,ntotal_,coslat,1) ! rq
   call dfs_fft_driver(1,psfc,ib,jbw,1    ,q ,mt,jlg,1    ,jlg,1    ,1    ,0,coslat,1)       ! q
   call dfs_fft_driver(1,hsfc,ib,jbw,1    ,gz,mt,jlg,1    ,jlg,1    ,1    ,0,coslat,1)       ! gz
   if(iope) then
     z00=0.
   endif
   call mpbcastr(z00,1)
#else
   call sph_uv2dize(us,vs,di,ze)                 ! update di,ze
   call sph_fft_driver(1,tv,te,1)                ! update te
   call sph_fft_driver(1,qq,rq,ntotal_)          ! update rq
   call sph_fft_driver_2d(1,psfc,q  ,1)          ! update q
   call sph_fft_driver_2d(1,hsfc,gph,1)
!
! laplacian of geopotential
!
#ifdef MP
   call mpsp2f(gph,lnt22p_,gphf,lnt22_,1)
   spec=gphf
#else
   spec=gph
#endif
   if(iope) then
     z00=spec(1)
     call spcshfli(spec,lnt22_,1,jcap_,lwvdef)
     do j = 1,lnt2_
       spec(j)=spec(j)*snnp1(j)*grav/(rerth_*rerth_)
     enddo
     call spcshflo(spec,lnt22_,1,jcap_,lwvdef)
   endif
#ifdef MP
   call mpbcastr(z00,1)
   call mpsf2p(spec,lnt22_,gz,lnt22p_,1)         ! update gz
#endif
#endif /* DFS end */
!
   return
   end subroutine dcmip_update
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dcmip_write
!-------------------------------------------------------------------------------
   use paramodel         , only : levs_
   use module_file_write , only : file_write_bin
   use dcmip_grims       , only : phis,ps,u,v,pwl,t,q0,q1,q2,q3,q4,rho,zl,pl
#ifdef DFS
   use dfsvar            , only : ib,jbw,xlon,xlat,iope
#define LONF2S ib
#define LATG2S jbw
#else
   use paramodel         , only : LONF2S,LATG2S
   use comfphys          , only : xlon,xlat
   use comio             , only : iope
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   integer, parameter                            ::  nn=102
!
! local variables
!
   integer, parameter                            ::  n3d=12,n2d=4
   integer                                       ::  i,j,k,n
   real   , dimension(LONF2S,LATG2S,n3d*levs_)   ::  temp
   real   , dimension(LONF2S,LATG2S,n2d)         ::  temp2
   logical, save                                 ::  lfirst
   data lfirst/.true./
!-------------------------------------------------------------------------------
#ifdef DFS
#define IJK i,j,k
#else
#define IJK i,k,j
#endif
!
! initialize
!
   temp=0.
   temp2=0.
!
! save 3D
!
   do k = 1,levs_
     do j = 1,LATG2S
       do i = 1,LONF2S
         temp(i,j,         k)=u(IJK)
         temp(i,j, 1*levs_+k)=v(IJK)
         temp(i,j, 2*levs_+k)=pwl(IJK)
         temp(i,j, 3*levs_+k)=t(IJK)
         temp(i,j, 4*levs_+k)=q0(IJK)
         temp(i,j, 5*levs_+k)=q1(IJK)
         temp(i,j, 6*levs_+k)=q2(IJK)
         temp(i,j, 7*levs_+k)=q3(IJK)
         temp(i,j, 8*levs_+k)=q4(IJK)
         temp(i,j, 9*levs_+k)=rho(IJK)
         temp(i,j,10*levs_+k)=pl(IJK)
         temp(i,j,11*levs_+k)=zl(IJK)
       enddo
     enddo
   enddo
#undef IJK
!
! save 2D
!
   do j = 1,LATG2S
     do i = 1,LONF2S
       temp2(i,j,1)=xlon(i,j)
       temp2(i,j,2)=xlat(i,j)
       temp2(i,j,3)=phis(i,j)
       temp2(i,j,4)=ps(i,j)
     enddo
   enddo
!
! write
!
   if(lfirst) then
     lfirst=.false.
     if(iope) open(nn,file='dcmip.bin',form='unformatted')
   endif
   do n = 1,n3d
     k = (n-1)*levs_+1
     call file_write_bin(nn,temp(1,1,k) ,LONF2S,LATG2S,levs_,1)
   enddo
   do n = 1,n2d
     call file_write_bin(nn,temp2(1,1,n),LONF2S,LATG2S,1    ,1)
   enddo
#ifdef DFS
#undef LONF2S
#undef LATG2S
#endif
!
   return
   end subroutine dcmip_write
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dcmip_normalized_error(rq)
!-------------------------------------------------------------------------------
   use paramodel  , only : levs_,levh_
#ifdef MP
   use paramodel  , only : lonf2p_,latg2p_
#endif
#ifdef DFS
   use dfsvar     , only : iope,mt,jlg,ib,jbw,levsp,iba,jbwa,coslat,levhp
#else
   use paramodel  , only : LONF2S,LATG2S,LNT22S,lonf2_,latg2_
   use comio      , only : iope
#endif
   use dcmip_grims, only : icase,rho,q1,q2,q3,q4
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
#ifdef NISLQ
#ifdef DFS
   real   , dimension(ib,jbw,levh_)             ::  rq
#else
   real   , dimension(LONF2S,levh_,LATG2S)      ::  rq
#endif
#else /* ~NISLQ */
#ifdef DFS
   real   , dimension(mt,jlg,levhp)             ::  rq
#else
   real   , dimension(LNT22S,levh_)             ::  rq
#endif
#endif /* NISLQ end */
!
! local variables
!
   integer, parameter                           ::  nt=5
   integer                                      ::  n,i,j,k,kk
   real                                         ::  vol,q,qt                  ,&
                                                    l1_u,l2_u                 ,&
                                                    l1_l,l2_l
#ifdef DFS
   real   , dimension(iba,jbwa,nt*levs_)        ::  model,refer
   real   , dimension(iba,jbwa,levs_)           ::  den,l0_u,l0_l
   real   , dimension(ib ,jbw ,levs_)           ::  model_q1,model_q2         ,&
                                                    model_q3,model_q4         ,&
                                                    model_q5,q5
#else
   real   , dimension(lonf2_ ,latg2_ ,nt*levs_) ::  model,refer
   real   , dimension(lonf2_ ,latg2_ ,levs_)    ::  den,l0_u,l0_l
   real   , dimension(LONF2S ,levs_  ,LATG2S)   ::  model_q1,model_q2         ,&
                                                    model_q3,model_q4         ,&
                                                    model_q5,q5
#ifdef MP
   real   , dimension(lonf2p_,latg2p_,nt*levs_) ::  model_qq,refer_qq
   real   , dimension(lonf2p_,latg2p_,levs_)    ::  denp
#endif
#endif /* DFS end */
   real   , dimension(nt)                       ::  l1,l2,l0
!-------------------------------------------------------------------------------
   model=0.
   refer=0.
!
! wave to grid
!
#ifdef DFS
#ifdef NISLQ
   do k = 1,levs_
     do j = 1,jbw
       do i = 1,ib
         model_q1(i,j,k)=rq(i,j,3*levs_+k)
         model_q2(i,j,k)=rq(i,j,4*levs_+k)
         model_q3(i,j,k)=rq(i,j,5*levs_+k)
         model_q4(i,j,k)=rq(i,j,6*levs_+k)
       enddo
     enddo
   enddo
#else
   call dfs_fft_driver(-1,model_q1,ib,jbw,levs_,rq(1,1,3*levsp+1),mt,jlg,levsp,jlg,levsp,levs_,1,coslat,1)
   call dfs_fft_driver(-1,model_q2,ib,jbw,levs_,rq(1,1,4*levsp+1),mt,jlg,levsp,jlg,levsp,levs_,1,coslat,1)
   call dfs_fft_driver(-1,model_q3,ib,jbw,levs_,rq(1,1,5*levsp+1),mt,jlg,levsp,jlg,levsp,levs_,1,coslat,1)
   call dfs_fft_driver(-1,model_q4,ib,jbw,levs_,rq(1,1,6*levsp+1),mt,jlg,levsp,jlg,levsp,levs_,1,coslat,1)
#endif /* NISLQ end */
   do k = 1,levs_
     do j = 1,jbw
       do i = 1,ib
         model_q5(i,j,k) = 0.3*(model_q1(i,j,k) + model_q2(i,j,k) + model_q3(i,j,k) ) + model_q4(i,j,k)
               q5(i,j,k) = 1.
       enddo
     enddo
   enddo
#ifdef MP
   call mpgp2fd(model_q1,ib,jbw,model(1,1,        1),iba,jbwa,levs_)
   call mpgp2fd(model_q2,ib,jbw,model(1,1,  levs_+1),iba,jbwa,levs_)
   call mpgp2fd(model_q3,ib,jbw,model(1,1,2*levs_+1),iba,jbwa,levs_)
   call mpgp2fd(model_q4,ib,jbw,model(1,1,3*levs_+1),iba,jbwa,levs_)
   call mpgp2fd(model_q5,ib,jbw,model(1,1,4*levs_+1),iba,jbwa,levs_)
   call mpgp2fd(      q1,ib,jbw,refer(1,1,        1),iba,jbwa,levs_)
   call mpgp2fd(      q2,ib,jbw,refer(1,1,  levs_+1),iba,jbwa,levs_)
   call mpgp2fd(      q3,ib,jbw,refer(1,1,2*levs_+1),iba,jbwa,levs_)
   call mpgp2fd(      q4,ib,jbw,refer(1,1,3*levs_+1),iba,jbwa,levs_)
   call mpgp2fd(      q5,ib,jbw,refer(1,1,4*levs_+1),iba,jbwa,levs_)
   call mpgp2fd(     rho,ib,jbw,den                 ,iba,jbwa,levs_)
#endif
#define lonf2_ iba
#define latg2_ jbwa
#else /* SPH */
#ifndef NISLQ
   call sph_fft_driver(-1,model_q1,rq(1,3*levs_+1),1)
   call sph_fft_driver(-1,model_q2,rq(1,4*levs_+1),1)
   call sph_fft_driver(-1,model_q3,rq(1,5*levs_+1),1)
   call sph_fft_driver(-1,model_q4,rq(1,6*levs_+1),1)
#endif /* ~NISLQ end */
#ifdef MP
#define MODEL model_qq
#define REFER refer_qq
#define DEN denp
#endif
   do k = 1,levs_
     do j = 1,LATG2S
       do i = 1,LONF2S
#ifdef NISLQ
         MODEL(i,j,        k)=rq(i,3*levs_+k,j)
         MODEL(i,j,  levs_+k)=rq(i,4*levs_+k,j)
         MODEL(i,j,2*levs_+k)=rq(i,5*levs_+k,j)
         MODEL(i,j,3*levs_+k)=rq(i,6*levs_+k,j)
#else
         MODEL(i,j,        k)=model_q1(i,k,j)
         MODEL(i,j,  levs_+k)=model_q2(i,k,j)
         MODEL(i,j,2*levs_+k)=model_q3(i,k,j)
         MODEL(i,j,3*levs_+k)=model_q4(i,k,j)
#endif
         MODEL(i,j,4*levs_+k)=0.3*( MODEL(i,j,k)+MODEL(i,j,levs_+k)+MODEL(i,j,2*levs_+k) ) + MODEL(i,j,3*levs_+k)
         REFER(i,j,        k)=q1(i,k,j)
         REFER(i,j,  levs_+k)=q2(i,k,j)
         REFER(i,j,2*levs_+k)=q3(i,k,j)
         REFER(i,j,3*levs_+k)=q4(i,k,j)
         REFER(i,j,4*levs_+k)=1.
         DEN(i,j,k)=rho(i,k,j)
       enddo
     enddo
   enddo
#ifdef MP
#undef MODEL
#undef DCMIP
#undef DEN
!
! partial to full
!
   call mpgp2f(model_qq,lonf2p_,latg2p_,model,lonf2_,latg2_,nt*levs_)
   call mpgp2f(refer_qq,lonf2p_,latg2p_,refer,lonf2_,latg2_,nt*levs_)
   call mpgp2f(denp    ,lonf2p_,latg2p_,den  ,lonf2_,latg2_,levs_  )
#endif /* MP end */
#endif /* DFS end */
!
! open file
!
   open(313,file='dcmip_norm.txt',form='formatted')
!
! Normalized error norms
!
   if(iope) then
     do n = 1,nt
       l1_u=0.  ; l2_u=0. ; l0_u=0.
       l1_l=0.  ; l2_l=0. ; l0_l=0.
       do k = 1,levs_
         do j = 1,latg2_
           do i = 1,lonf2_
             kk=(n-1)*levs_+k
             vol=den(i,j,k)
             q=model(i,j,kk)
             qt=refer(i,j,kk)
             !
             l0_u(i,j,k)=abs(q-qt)
             l0_l(i,j,k)=abs(qt)
             l1_u=l1_u+l0_u(i,j,k)*vol
             l1_l=l1_l+l0_l(i,j,k)*vol
             l2_u=l2_u+(q-qt)**2*vol
             l2_l=l2_l+qt**2*vol
           enddo
         enddo
       enddo
       l1(n)=l1_u/l1_l
       l2(n)=sqrt(l2_u/l2_l)
       l0(n)=maxval(l0_u)/maxval(l0_l)
     enddo
!
     write(313,201) 'for Q1, L1= ',l1(1),'  L2= ',l2(1),'  L0= ',l0(1)
     write(313,201) 'for Q2, L1= ',l1(2),'  L2= ',l2(2),'  L0= ',l0(2)
     write(313,201) 'for Q3, L1= ',l1(3),'  L2= ',l2(3),'  L0= ',l0(3)
     write(313,201) 'for Q4, L1= ',l1(4),'  L2= ',l2(4),'  L0= ',l0(4)
     write(313,201) 'for Q5, L1= ',l1(5),'  L2= ',l2(5),'  L0= ',l0(5)
   endif
201 format(a,f12.9,a,f12.9,a,f12.9)
#ifdef DFS
#undef lonf2_
#undef latg2_
#endif
!
   return
   end subroutine dcmip_normalized_error
!-------------------------------------------------------------------------------
#ifndef DFS
!
!-------------------------------------------------------------------------------
   subroutine sph_uv2dize(ug,vg,di,ze)
!-------------------------------------------------------------------------------
   use paramodel  , only : LONF2S,LATG2S,latg2_,levs_,LEVSS                   ,&
                           LNT22S,LCAP22S,JCAP1S
   use comfgrid   , only : wgt
   use comgpln    , only : qvv,qww
#ifdef MP
   use paramodel  , only : lonf22p_,lonf22_,latg2p_,levsp_                    ,&
                           lln22p_,lnt22p_,lcap22p_
   use commpi     , only : latdef,lwvdef,latstr,latlen,lwvstr,lwvlen,mype
#else
   use comfgrid   , only : latdef,lwvdef
#endif
#ifdef REDUCE_GRID
   use comreduce
#else
   use paramodel  , only : lonf_,jcap1_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   real   , dimension(LONF2S,levs_,LATG2S), intent(in)     ::  ug,vg
   real   , dimension(LNT22S,levs_)       , intent(out)    ::  di,ze
!
! local variables
!
   integer                                                 ::  lat,lat1,lat2,&
                                                               latdon,lcapf ,&
                                                               lonff,llstr  ,&
                                                               llens,llensd ,&
                                                               lan,i,j,k    ,&
                                                               jstr,jend
#ifdef MP
   real   , dimension(lonf22p_,2*levs_ ,latg2p_)           ::  gra
   real   , dimension(lonf22_ ,2*levsp_,latg2p_)           ::  anf
   real   , dimension(lln22p_ ,2*levsp_)        , target   ::  uua,dia
   real   , dimension(:,:)                      , pointer  ::  vva,zea
#else
   real   , dimension(LNT22S,levs_)             , target   ::  uu
   real   , dimension(:,:)                      , pointer  ::  vv
#endif
   real   , dimension(LCAP22S ,2*LEVSS ,latg2_ )           ::  anl
   real   , dimension(2,JCAP1S,2*LEVSS ,latg2_ )           ::  flp,flm
   real   , dimension(2,JCAP1S,2*LEVSS)                    ::  anltop
!-------------------------------------------------------------------------------
!
! initialize
!
   di=0.
   ze=0.
#ifdef MP
   gra=0.
   anf=0.
   uua=0.
   dia=0.
#else
   uu=0.
#endif
   anl=0.
   flp=0.
   flm=0.
   anltop=0.
!
! pointer
!
#ifdef MP
   vva=>uua(1:lln22p_,levsp_+1:2*levsp_)
   zea=>dia(1:lln22p_,levsp_+1:2*levsp_)
#else
   vv=>uu(1:LNT22S,levs_+1:2*levs_)
#endif
!
! index for wave and grid
!
#ifdef MP
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
#define ANL gra
#else
   llstr=0
   llens=jcap1_
#endif
!
   do j = 1,LATG2S
     do k = 1,levs_
       do i = 1,LONF2S
         ANL(i,      k,j)=ug(i,k,j)
         ANL(i,levs_+k,j)=vg(i,k,j)
       enddo
     enddo
   enddo
#ifdef MP
#undef ANL
   call mpnx2nk(gra,lonf22p_,2*levs_,anf,lonf22_,2*levsp_,latg2p_,levs_,levsp_,1,1,2)
#define ANL anf
#endif
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
     lcapf=lcapd(latdef(lat))
     lonff=lonfd(latdef(lat))
#else
     lcapf=jcap1_
     lonff=lonf_
#endif
     call sph_wave2grid(ANL(1,1,lan),ANL(1,1,lan),4*LEVSS,lcapf,lonff,latdef(lat),-1)
   enddo
#ifdef MP
#undef ANL
   call mpny2nl(anf,lonf22_,latg2p_,anl,lcap22p_,latg2_,2*levsp_,1,2*levsp_)
#endif
!
   lat1=1
   lat2=latg2_
   latdon=0
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif /* REDUCE_GRID end */
     call sph_grid2wave(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),llensd,2*LEVSS)
   enddo
!
#ifdef MP
#define UU uua
#define VV vva
#define DI dia
#define ZE zea
#endif
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif
     call sph_comp_coeff2(   flp(1,1,1,lan),flm(1,1,1,lan),UU,qww(1,lat)      ,&
                             llstr,llensd,lwvdef,2*LEVSS)
     call sph_sum_wind_coeff(flp(1,1,      1,lan),flm(1,1,      1,lan)        ,&
                             flp(1,1,LEVSS+1,lan),flm(1,1,LEVSS+1,lan)        ,&
                             anltop(1,1,1)       ,anltop(1,1,1+LEVSS)         ,&
                             qvv(1,lat),wgt(lat),llstr,llensd,lwvdef,LEVSS)
   enddo
!
   call sph_wind2divor(UU,VV,DI,ZE,anltop(1,1,1),anltop(1,1,1+LEVSS),llstr,llens,lwvdef)
!
#ifdef MP
#undef UU
#undef VV
#undef DI
#undef ZE
   call mpnk2nn(dia,lln22p_,levsp_,di,lnt22p_,levs_,1)
   call mpnk2nn(zea,lln22p_,levsp_,ze,lnt22p_,levs_,1)
#endif
!
   return
   end subroutine sph_uv2dize
!-------------------------------------------------------------------------------
#endif /* ~DFS end */
#endif /* DCMIP end */
