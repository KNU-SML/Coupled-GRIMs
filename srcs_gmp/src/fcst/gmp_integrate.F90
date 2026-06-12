#include "define.h"
   subroutine gmp_integrate
!-------------------------------------------------------------------------------
!
! subroutine: gmp_integrate   
!
! abstract: make global forecast 
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! references : 
!   hong et al. (2013, apjas): global/regional integrated model system (grims)
!   park et al. (2013, mwr): dfs dynamical core
!   byun and hong (2007, j. climate): single-column model (SMP)
!   kanamitsu et al. (2002, bams): ncep dynamical seasonal forecast system 2000
!
! ::: structure ::: 
!  
!    [gmp_integrate]
!      |
!      |--- [phys_main_driver]  
!      |
!      |---- A-1) [dyn_sph_driver]  
!      | |-- A-2) [dyn_dfs_driver]  
!      |
!      |---- B-1) [sph_horizontal_diffusion]  
!      | |-- B-2) [dfs_diffusion_driver]  
!      |
!      |--- [dyn_partial_time_filter] *
!      |--- [dyn_full_time_filter] *
!      |--- [rdscmbd] *
!
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S,LNT22S,levs_,levh_
#ifndef DFS
   use paramodel, only : lnt22_
#endif
   use varsfc,only     : lalbd_
   use comsfc          ! sfcftyp,sfcfcs,sfcv
#ifdef DFS
   use dfsvar, only    : EIGVAL,EIGVEC,AEIGVEC                                ,&
                         d1,v1,t1,q1,psl1,sf1                                 ,&
#ifdef PERT
                         d0                                                   ,&
#endif
#ifdef BARO_TEST
                         dio,prs                                              ,&
#endif
                         d2,v2,t2,q2,psl2                                     ,&
                         d3,v3,t3,q3,psl3                                     ,&
                         PSJ,DJ,TJ,vj,qj,tavexy                               ,&
                         gamma_v,gamma_d,gamma_t,gamma_q,lap_dim              ,&
                         damp2,coslat,spdmax                                  ,&
                         sl=>sigma,si=>sigmafull,del=>delsig                  ,&
                         levsp,levs,levhp,levh,ntotal,kls,ngs,nge             ,&
                         lnt2,mt,mtg,jl,jlg,ib,jbw,iope,mta
#else
#ifndef HYBRID
   use module_sph_semi_implicit,only : sph_semi_hydro
#else
   use module_sph_semi_implicit,only : sph_semi_hydro_hybrid
#endif
#ifndef SMP
   use module_sph_semi_implicit,only : sph_semi_gwave
#endif
#endif
#ifdef MP
   use commpi, only                  : mype,npes, master
#endif
#ifdef DG3
   use diag_3d_module, only          : diag_3d_zero_out
#endif
   use constant, only                : pi_, rhoh2o_
   use comfibm
   use comcon
   use comgpd
   use comfgsm
#ifdef LFM
   use comlfm
#endif
#ifdef RIVER
   use comfrivh
#endif
#ifdef OMLWRF
   use module_ocean_slab_wrf
#endif
!
#ifdef SMP
   use comfcst, only  : dtbdy,curtime,vvel,hour1,hour2,hours,houre,ftim1,ftim2
#ifdef CLM_CWF
   use comfcst, only  : wdiv,hadq
#endif
   use comfcst, only  : uo1,vo1,to1,qo1,wo1,ps1                               ,&
                        uo2,vo2,to2,qo2,wo2,ps2                               ,&
#ifdef REV_FRC
                        ac1, ac2                                              ,&
#endif
#ifdef RLX_FRC
                        tscl1, tscl2                                          ,&
#endif
                        at1,aq1,dv1,at2,aq2,dv2
#ifdef SMP_NUDGING_SFC
#ifdef SMP_NUDGING_FLX
   use comfcst, only  : scmsh,scmlh
#endif
#ifdef SMP_NUDGING_RAD
   use comfcst, only  : psfc,tsfc,qsfc,prec,scmlwup,scmlwdn,scmswup,scmswdn
#endif
#endif
#endif       /* SMP */
!
#ifdef SAS_DIAG
   use comfcst, only : dcu,dcv,dct,dcq,dch,fcu,fcd                            ,&
                       deltb,delqb,delhb,cbmf,dlt,dlq,dlh,cvbp,cvtp
#endif
#ifdef EXPLICIT_CLOUDINESS
   use comfcst, only : qcicps,qrscps,taucld,cldwp,cldip
#endif
   use module_file_write, only : file_write_bin,   file_write_cps_bin,         &
                                 file_write_point, file_write_lfm
#ifdef NISLQ
   use nislq, only : slq_q1,slq_q2,slq_q3,slq_psfc2,slq_u2,slq_v2,slq_w2
#endif
#ifdef DCMIP
   use dcmip_grims, only : icase
#endif
!soojin_couple
#ifdef AOMG
   use coupling
   use module_trans
   use mod_oasis
   use couple_oasis, only : cpl_recv
   use comsfc,    only : tsea, snwdph, slmsk
   use paramodel, only : lonf_, latg_, lonfp_, latgp_, lonf2p_, latg2p_, lonf2_,latg2_
#endif /*AOMG*/
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include <abort.h>
#ifdef DFS
#define LNT22S lnt2
   integer                          ::  lotg,kk,lapdim,kvar
#endif
   integer ncho
   character(len=128)               ::  fno
   character(len=4)                 ::  sfcftypin
   data sfcftypin /'dumy'/
#ifdef DFS
   real                             ::  WNM2D1(mt,jlg)
#endif
!
! local variables
!
   integer                          ::  i,j,m,k,k1,k2,iday
   real                             ::  cvmod
#ifdef TIME_CHECK
   integer                          ::  hz,clock0,clock1
   real                             ::  t
#endif
#ifdef DCMIP
   real                                     ::  dum,utime
#ifdef DFS
   real   , dimension(mt,jlg,levhp)         ::  dum3
#ifdef NISLQ
   real   , dimension(ib,jbw,levh)          ::  dumq
#endif
#else /* SPH */
   real   , dimension(LNT22S,levh_)         ::  dum3
#ifdef NISLQ
   real   , dimension(LONF2S,levh_,LATG2S)  ::  dumq
#endif
#endif /* DFS end */
#endif /* DCMIP end */
!soojin_couple
#ifdef AOMG
   integer              :: len, ierror, intsec, jj
   real, allocatable    :: work(:),work2(:,:)
   real, allocatable    :: exchp(:,:)
   real                 :: msk(LONF2S,LATG2S),tice(LONF2S,LATG2S), snowdp(LONF2S,LATG2S)
   real                 :: tseatmp(LONF2S,LATG2S) !soojin couple
#endif /*AOMG*/
!
!  start one time step
!
   if( stepone ) then
     if (jdt.gt.limlow) inistp=0
#ifdef DGP
     if (jdt.gt.limlow) isave = 0
#endif
     kdt = jdt - limlow + 1
     lastep=jdt.eq.1
   else
     kdt=jdt
     lastep=jdt.eq.maxstp
#ifdef DGP
     if(npoint.gt.0.and.itnum.gt.nstken.and.iope) then
#ifndef NOPRINT
       if(iope) write(6,*)'DGP points disabled - time levels exceed ',nstken
#endif
       npoint=0
     endif
#endif
   endif
! 
!  The following print is kept
!
#ifndef NOPRINT
   if(iope) then
     if(stepone) then
       write(6,81)kdt
     else
       write(6,82)kdt
     endif
   endif
#else
   if(iope) then
     if(mod(kdt-1,10).eq.0) then
       if(stepone) then
         write(6,81)kdt
       else
         write(6,82)kdt
       endif
     endif
   endif
#endif
81 format(1h ,'forward step: kdt in gmp_integrate=',i5)
82 format(1h ,'full step   : kdt in gmp_integrate=',i5)
! 
   if( .not. stepone ) then
!
!  set switch for saving kuo data (for interactive clouds)..
!
     cvmod= mod (solhr+dthr,dtcvav)
     if(cvmod.lt.hdthr.or.cvmod.ge.dtcvav-hdthr) then
       clstp=min(dtcvav,(shour+deltim)/3600.)
     elseif(clstp.gt.0.) then
       clstp=0.
     else
       clstp=-10.
     endif
   endif
#ifdef ONELOOP
!
!  one loop includes rad_driver, dyn_sph_driver, and phys_main_driver.
!
   call dyn_driver_onestep
#ifdef DBG
!
   if(iope) then
     if(stepone) then
       write(6,*)' forward time step:  after dyn_driver_onestep '
     else
       write(6,*)' full time step: after dyn_driver_onestep '
     endif
   endif
#ifdef DFS
   call dyn_comp_rms(PSJ,DJ,TJ,vj,del,QJ)
#else
   call dyn_comp_rms(z ,x ,y ,w ,del,rt)
#endif
#endif /* DBG end */
#else /* ~ONELOOP */
#ifndef NO_PHYSICS
#ifdef SMP
!
! smp off_line test
!
   call rdscmbd
#endif
!
#ifdef PSPLIT
!
! process splitting physics tendency initialization
!
   if( jdt.eq.limlow) then
     call phys_main_driver
   endif
#endif
#endif /* ~NO_PHYSICS end */
!
!  dynamic loop
!
#ifdef TIME_CHECK
   if (.not.stepone.and.kdt.eq.2) then
     call system_clock(count_rate=hz)
     call system_clock(count=clock0)
   endif
#endif
!
!soojin_couple!!!!!!!!!!!!!
#ifdef AOMG
!receive parts
   intsec = ( int(fhour)*3600 ) + int(shour) -1200.
   print*, "intsec_receive =", intsec
   ! first step intsec -> skip exchange











!!!!!!!!!!!!!!!!!end receive variable
     endif !-- end for mod intsec
    endif!for skip exchange forward
#endif /* AOMG */
!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef DCMIP
   call rad_read_prepare(fhour,idate)
#endif /* DCMIP end */
!
#ifdef DFS
   call dyn_dfs_driver
#else
   call dyn_sph_driver
#endif
#ifdef DBG
   if(iope) then
     if(stepone) then
       write(6,*)' forward time step:  after dyn_sph_driver'
     else
       write(6,*)' full time step: after dyn_sph_driver'
     endif
   endif
#ifdef DFS
   call dyn_comp_rms(PSJ,DJ,TJ,vj,del,QJ)
#else
   call dyn_comp_rms(z ,x  ,y  ,w  ,del,rt)
#endif
#endif /* DBG end */
#endif /* ONELOOP end */
!
! advance shour (second) and thour (total hour)
!
   if(stepone) then
     shour=deltim
     dtpost=deltim
   endif
   shour=shour+deltim
   dtpost=dtpost+deltim
   if(.not.stepone) thour=fhour+ifix(shour/3600.+0.5e0)
#ifdef SMP
   write(6,*) 'fhour=',fhour,' shour=',shour
   curtime = fhour + shour/3600.
#endif
!
#ifndef SMP
!  semi-implicit
!
#ifdef DFS
#ifdef HYBRID
   call dfs_divergence_hybrid(deltim,SF1,D1,T1,PSL1,DJ,TJ,PSJ,                 &
               EIGVAL,EIGVEC,AEIGVEC,D2,T2,PSL2,D3)
#else
   call dfs_divergence(deltim,SF1,D1,T1,PSL1,DJ,TJ,PSJ,                        &
               EIGVAL,EIGVEC,AEIGVEC,D3)
#endif
!
! integrate vorticity and specific humidity
!
   call dfs_diffusion(deltim,D1,V1,Q1,T1,PSL1,VJ,QJ,TJ,PSJ,                    &
               D3,V3,Q3,T3,PSL3)
#else /* SPH */
#ifndef HYBRID
   call sph_semi_hydro(dim,tem,qm,x,y,z,uln,vln,lnts2,lnoffset)
#else
   call sph_semi_hydro_hybrid(dim,tem,qm,x,y,z,di,te,q,uln,vln,lnts2,lnoffset)
#endif
#endif /* DFS end */
#ifdef DBG
   if(iope) then
     if(stepone) then
       write(6,*)' forward time step:  after sph_semi_hydro '
     else
       write(6,*)' full time step: after sph_semi_hydro '
     endif
   endif
#ifdef DFS
   call dyn_comp_rms(psl3,d3,t3,v3,del,q3)
#else
   call dyn_comp_rms(z ,x  ,y  ,w  ,del,rt)
#endif
#endif /* DBG end */
!
!  horizontal diffusion
!
#ifdef DFS
#ifdef PERT
! 
! pertubation diffusion
!
   do k = 1, (1+ntotal)*levs
     t3(1:mt, 1:jl, k) = t3(1:mt, 1:jl, k) - t1(1:mt, 1:jl, k)
   enddo
#endif
   do k = 1,levsp/2
     k1=k*2-1
     k2=k*2
     kk=kls+k2-1
     CALL dfs_diffusion_driver(v3(1,1,k1),v3(1,1,k2),gamma_v,LAP_DIM,sl(kk),1)
     CALL dfs_diffusion_driver(d3(1,1,k1),d3(1,1,k2),gamma_d,LAP_DIM,sl(kk),2)
     CALL dfs_diffusion_driver(t3(1,1,k1),t3(1,1,k2),gamma_t,LAP_DIM,sl(kk),3)
   end do ! k=1,levsp/2
   WNM2D1(1:mt,1:jlg)=0.0
   IF(MOD(levsp,2).eq.1) then
     k1=levsp
     kk=kls+k1-1
     CALL dfs_diffusion_driver(v3(1,1,k1),WNM2D1,gamma_v,LAP_DIM,sl(kk),1)
     CALL dfs_diffusion_driver(d3(1,1,k1),WNM2D1,gamma_d,LAP_DIM,sl(kk),2)
     CALL dfs_diffusion_driver(t3(1,1,k1),WNM2D1,gamma_t,LAP_DIM,sl(kk),3)
   ENDIF
#ifndef NISLQ
   do kvar = 1,ntotal
     do k = 1,levsp/2
       k1=k*2-1+(kvar-1)*levsp
       k2=k1+1
       kk=kls+2*k-1
       CALL dfs_diffusion_driver(q3(1,1,k1),q3(1,1,k2),gamma_q,LAP_DIM,sl(kk),4)
     enddo
     IF(MOD(levsp,2).eq.1) then
       k1=levsp+(kvar-1)*levsp
       kk=kls+levsp-1
       CALL dfs_diffusion_driver(q3(1,1,k1),WNM2D1,gamma_q,LAP_DIM,sl(kk),4)
     ENDIF
   enddo
#endif /* ~NISLQ end */
#ifdef PERT
   ! 
   ! pertubation diffusion
   !
   do k = 1,(1+ntotal)*levs
     t3(1:mt, 1:jl, k) = t3(1:mt, 1:jl, k) + t1(1:mt, 1:jl, k)
   end do
#endif
#else /* SPH */
   call sph_horizontal_diffusion(rt,w,deltim,qm,sl,x,y,lnts2,lnoffset)
#endif /* DFS end */
#ifdef DCMIP
!
! dcmip update
!
   !if(mod(kdt,max(ksig,1)).eq.0) call dcmip_write
   utime = dble(fhour*3600.+shour)
   if(iope) print *,'kdt,deltim,shour,utime',kdt,deltim,shour,utime
#ifdef DFS
#define DGPH sf1
#define DSLP psl3
#define DDIV d3
#define DVOR v3
#define DTEM t3
#else
#define DGPH gz
#define DSLP z
#define DDIV x
#define DVOR w
#define DTEM y
#endif /* DFS end */
   if (icase.eq.11 .or. icase.eq.12 .or. icase.eq.13) then
#ifdef NISLQ
     call dcmip_update(utime,DGPH,DSLP,DDIV,DVOR,DTEM,dum3,dum,dumq)
#else
     call dcmip_update(utime,DGPH,DSLP,DDIV,DVOR,DTEM,dum3,dum)
#endif /* NISLQ end */
   endif
#undef DGPH
#undef DSLP
#undef DDIV
#undef DVOR
#undef DTEM
#endif /* DCMIP end */
#ifdef DBG
   if(iope) then
     if(stepone) then
       write(6,*)' forward time step:  after sph_horizontal_diffusion '
     else
       write(6,*)' full time step: after sph_horizontal_diffusion '
     endif
   endif
#ifndef DFS
   call dyn_comp_rms(z ,x  ,y  ,w  ,del,rt)
#endif
#endif /* DBG end */
#ifdef NISLQ
!
! forward semi-Largrangian advection for specific humidity (t+dt)
!
   call nislq_mono_advect(deltim,slq_psfc2,slq_u2,slq_v2,slq_w2,slq_q1,slq_q3)
#ifdef CHEM
!
! forward semi-Largrangian advection for chemical tracers (t+dt)
!
   call transport_chem
#endif
#endif /* NISLQ end */
#ifdef NISLQ_MONOMASS
!
! 3D mass conservation for nislq
!
   if(kdt.ne.1) call nislq_mass_adjustment(slq_q1,slq_q3,qm,z)
#endif /* NISLQ_MONOMASS */
#endif /* SMP end */
!
!  partial time filter to n and save in n-1
!
   if( .not. stepone ) then
#ifdef DFS
#define TEM t1
#define TE t2
#define DIM d1
#define DI d2
#define ZEM v1
#define ZE v2
#define RM q1
#define RQ q2
#define QM psl1
#define Q psl2
#define LEVS levsp
#define LEVH levhp
#else
#define LEVS levs_
#define LEVH levh_
#endif
#ifdef NISLQ
#define RM slq_q1
#define RQ slq_q2
#endif
     call dyn_partial_time_filter(TEM,TE,DIM,DI,ZEM,ZE,RM,RQ,QM,Q,             &
                                filta,LNT22S,lnts2,LEVS,LEVH)
#ifdef DBG
     if(iope) print *,' full time step: after dyn_time_filter t-dt'
     call dyn_comp_rms(QM,DIM,TEM,ZEM,del,rm)
#endif
#undef TEM
#undef TE
#undef DIM
#undef DI
#undef ZEM
#undef ZE
#undef RM
#undef RQ
#undef QM
#undef Q
#undef LEVS
#undef LEVH
   endif
#ifdef DFS
#ifdef ALIASED
   call dfs_cut_alias(v3,mt,ngs,nge,(3+ntotal)*levsp+1)
#endif
   psl2(1:mt,1:jlg)=psl3(1:mt,1:jlg)
#else /* SPH */
   q(1:lnts2)=z(1:lnts2,1)
#endif /* DFS end */
!
#ifdef DCMIP
!
! simple physics for dcmip
!
   if(icase.eq.42 .or. icase.eq.51) call phys_main_driver_dcmip
#endif
#ifndef NO_PHYSICS
#ifndef ONELOOP
#ifndef PSPLIT
!
!  physics loop
!
#ifdef OMLWRF
   if(mod(kdt-1,max(komlx,1)).eq.0.or.kdt.eq.1) then
           call ocean_ml_wrf_init(tsea,tml,t0ml,                               &
                                hml,h0ml,hml0,huml,hvml,tmoml,1,LONF2S,1,LATG2S)
#ifdef DBG
     if(iope) then
       write(6,*)' ocean mixed layer is initialized at kdt = ',kdt
       print*,1,LONF2S,1,LATG2S
       print*,' tsea ',tsea(LONF2S/2,LATG2S/2)
       print*,' t0ml ',t0ml(LONF2S/2,LATG2S/2)
       print*,' h0ml ',hml0(LONF2S/2,LATG2S/2)
     endif
#endif /* DBG end */
   endif
#endif /* OMLWRF end */
!
   call phys_main_driver
#ifdef DBG
!
   if(iope) then
     if(stepone) then
       write(6,*)' forward time step:  after phys_main_driver t+dt'
     else
       write(6,*)' full time step: after phys_main_driver t+dt'
     endif
   endif
#ifdef DFS
   call dyn_comp_rms(psl2,d3  ,t3  ,v3  ,del,q3)
#else
   call dyn_comp_rms(q ,x  ,y  ,w  ,del,rt)
#endif /* DFS end */
#endif /* DBG end */
#endif /* ~PSPLIT end */
#endif /* ~ONELOOP end */
#endif /* ~NO_PHYSICS end */
#ifdef RIVER
!
!  river discharge
!
   call phys_river_main_driver
#endif
#ifndef SMP
#ifndef DFS
!
!  zonal wind damping
!
   call sph_wind_damping(x,w,y,rt,deltim,uln,vln,spdmax,snnp1,lnts2,lnoffset)
#endif /* ~DFS end */
#ifdef DBG
   if(iope) then
     if(stepone) then
       write(6,*)' forward time step:  after sph_wind_damping t+dt'
     else
       write(6,*)' full time step: after sph_wind_damping t+dt'
     endif
   endif
#ifdef DFS
   call dyn_comp_rms(psl3,d3,t3,v3,del,q3)
#else
   call dyn_comp_rms(q ,x  ,y  ,w  ,del,rt)
#endif
#endif /* DBG end */
#endif /* ~SMP end */
!
!  update n by n+1
!
#ifdef DFS
#ifdef ALIASED
   call dfs_cut_alias(v3,mt,ngs,nge,(3+ntotal)*levsp+1)
#endif
   do k = 1,levsp
     do j = 1,jlg
       do m = 1,mt
         d2(m,j,k)=d3(m,j,k)
         v2(m,j,k)=v3(m,j,k)
         t2(m,j,k)=t3(m,j,k)
       enddo
     enddo
   enddo
#ifdef NISLQ
   slq_q2(1:ib,1:jbw,1:levh)=slq_q3(1:ib,1:jbw,1:levh)
#else
   do k = 1,levhp
     do j = 1,jlg
       do m = 1,mt
         q2(m,j,k)=q3(m,j,k)
       enddo
     enddo
   enddo
#endif
#else /* SPH */
   do k = 1,levs_
     do j = 1,lnts2
       di(j,k)=x(j,k)
       ze(j,k)=w(j,k)
       te(j,k)=y(j,k)
     enddo
   enddo
#ifdef NISLQ
   slq_q2(1:LONF2S,1:levh_,1:LATG2S)=slq_q3(1:LONF2S,1:levh_,1:LATG2S)
#else
   do k = 1,levh_
     do j = 1,lnts2
       rq(j,k)=rt(j,k)
     enddo
   enddo
#endif
#endif /* DFS end */
#ifdef DCMIP
!
! dcmip normalized error (L1,L2,L00)
!
!mskoo   if(lastep .and. .not. stepone) then
   if(lastep ) then
     if(icase.eq.11 .or. icase.eq.12 .or. icase.eq.13) then
#ifdef NISLQ
       call dcmip_normalized_error(slq_q2)
#else
#ifdef DFS
       call dcmip_normalized_error(q2)
#else
       call dcmip_normalized_error(rq)
#endif
#endif /* NISLQ end */
     endif
   endif
#endif /* DCMIP end */
#ifdef TRCON
!
! check global conservation
!
#ifdef DFS
#define QM psl1
#define Q psl2
#define RM q1 
#define RQ q2 
#endif
   if(kdt.eq.1) then 
#ifdef NISLQ
     call dyn_3d_contest(QM,slq_q1)
#else
     call dyn_3d_contest(QM,RQ)
#endif
   else
#ifdef NISLQ
     call dyn_3d_contest(Q,slq_q2)
#else
     call dyn_3d_contest(Q,RQ)
#endif
   endif
#undef QM
#undef Q
#undef RM
#undef RQ
#endif /* TRCON end */
!
!  complete time filter for n and save to n-1
!
   if( .not. stepone ) then
#ifdef DFS
#define TEM t1
#define TE t2
#define DIM d1
#define DI d2
#define ZEM v1
#define ZE v2
#define RM q1
#define RQ q2
#define QM psl1
#define Q psl2
#define LEVS levsp
#define LEVH levhp
#else
#define LEVS levs_
#define LEVH levh_
#endif
#ifdef NISLQ
#define RM slq_q1
#define RQ slq_q2
#endif
     call dyn_full_time_filter(TEM,TE,DIM,DI,ZEM,ZE,RM,RQ,                     &
                          filta,LNT22S,lnts2,LEVS,LEVH)
#ifdef DBG
     if(iope) write(6,*)' full time step: after dyn_full_time_filter'
     call dyn_comp_rms(Q,DI,TE,ZE,del,rq)
#endif /* DBG end */
#undef TEM
#undef TE
#undef DIM
#undef DI
#undef ZEM
#undef ZE
#undef RM
#undef RQ
#undef QM
#undef Q
#undef LEVS
#undef LEVH
   endif
#ifdef LFM
!
!  lfm filtering
!
   if( .not. stepone ) then
     ifstep=kdt+ipstep
     call dyn_low_freq_model(ifstep,thour)
   endif
#endif
!                    
#ifdef MP            
#ifdef GMPDAMP       
   if(iope) print *,' start dyn_global_nudging'
   call dyn_global_nudging
#endif 
#endif
!
#ifdef CHEM
   call integrate_chem
#endif
!
! update the matrix in case of changing deltim within stepone
!
   if( stepone ) then
     deltim=deltim*2.e0
#ifndef SMP
#ifdef DFS
     call dfs_diffusion_coef(deltim,gamma_v,mta,1)
     call dfs_diffusion_coef(deltim,gamma_d,mta,2)
     call dfs_diffusion_coef(deltim,gamma_t,mta,3)
#ifndef NISLQ
     call dfs_diffusion_coef(deltim,gamma_q,mta,4)
#endif
#else /* SPH */
     call sph_semi_gwave(deltim,am,bm,gv,sv,cm)
#ifdef DBG
     if( iope ) print *,' forward time: sph_semi_gwave for deltim= ',deltim
#endif
#endif /* DFS end */
#endif /* ~SMP end */
   endif
#ifdef DFI
! dyn_digital_filter : 
!        call digital filter initialization every step if con(3).gt.0.0
!
   if(con(3).ne.0.0) call dyn_digital_filter(1,n1,con(3),shour/3600.,solsec)
#endif
!
!  advance solhr
!
   if( .not. stepone ) then
     solsec=solsec+deltim
     solhr=solsec/3600.e0
     iday=solhr/24.e0
     solhr=solhr-iday*24.e0
   endif
#ifdef DGP
!
!   for grid point diag advance itnum, if proper time, and set isave
!
   isave = 0
#endif
#ifdef DFI
   if(numsum.lt.0) then
#endif
#ifdef DGP
     if (ikfreq.gt.1) then
       imodk = mod(jdt,ikfreq)
       if (imodk.eq.0) then
         isave = 1
         itnum = itnum + 1
       end if
     else
       isave = 1
       itnum = itnum + 1
     end if
#endif
#ifdef DFI
   end if
#endif
#ifdef RMP
!
! RMP : call main routine of regional forecast
!
   call rmp_integrate(fhour,shour,gz,q,te,di,ze,rq)
#endif
!
#ifndef BARO_TEST
#ifndef SMP
#ifdef DFS
#define GLOOPP dfs_mass_adjustment
#define PSL2 psl2
#define PSL1 psl1
#else
#define GLOOPP sph_mass_adjustment
#define PSL2 q
#define PSL1 qm
#endif
!
!  mass conservation, surface merge and output files
!
   if(mod(kdt,max(kpfix,1)).eq.0.or.lastep) then
     if( .not. stepone ) then
       call GLOOPP(PSL2,avprs0)
       call GLOOPP(PSL1,avprs0)
     endif
     avprs0=0.
     call GLOOPP(PSL2,avprs0)
   endif
#ifdef DBG
   if(iope) write(6,*)' full time step: after driver_sfcp '
#ifdef DFS
   call dyn_comp_rms(psl2,d2,t2,v2,del,q2)
#else
   call dyn_comp_rms(q ,di  ,te  ,ze  ,del,rq)
#endif
#endif /* DBG end */
#undef GLOOPP
#undef PSL2
#undef PSL1
#endif
#endif         /* BARO_TEST */
! 
#ifdef TIME_CHECK
   if (kdt.eq.40) then
     call system_clock(count=clock1)
     t=real(clock1-clock0)
#ifdef MP
     write(6,'(A,I4,F10.3)')'mype,40 sec=',mype,t/real(hz)
     call mpabort
#else
     write(6,'(A,F10.3)')'40 sec=',t/real(hz)
     stop
#endif
   endif
#endif
! write sigma and surface files
!
   nstep1: if( .not. stepone ) then
!
! write sigit for diagnostics and RMP boundary conditions
!
   if(mod(kdt,max(ksig,1)).eq.0.or.lastep) then
#ifdef BARO_TEST
     call dfs_fft_driver(-1,dio,ib,jbw,3*levs,D1,mt,jlg,3*levsp,jlg,           &
                     levsp,levs,3,coslat,1)
     call dfs_fft_driver(-1,prs,ib,jbw,1,PSL1,mt,jlg,1,jlg,                    &
                     levsp,levs,0,coslat,1)
     call file_write_bin(333,PRS,ib,jbw,1      ,0)
     call file_write_bin(333,dio,ib,jbw,3*levs,0)
#else

#ifdef DFS
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,sl,si,sf1,tavexy,   &
#else
     call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,ak5,bk5,sf1,tavexy, &
#endif
#else /* SPH */
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,sl,si,gz,z00,      &
#else
     call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,    &
#endif
#endif
#ifdef SMP
                    vvel,                                                      &
#endif
                    1)
#ifdef DBG
!
!  monitor sigma file
!
     if( iope )                                                                &
             write(6,*)'predicted full values at the end of forecst segment'
#ifdef DFS
       call dyn_comp_rms(psl2,d2,t2,v2,del,q2)
#else
       call dyn_comp_rms(q,di,te,ze,del,rq)
#endif
#endif
#endif          /* not BARO_TEST */
     endif
!
!  change surface field by merging 
!  (not executed at the beginning of time loop).
!
#if !defined(BARO_TEST) && !defined(NO_PHYSICS)
#ifndef SKIPSFCMRG
     if(mod(kdt,max(ksfcx,1)).eq.0) then
       call sfc_cycle_driver                                                   &
                  (sfcftypin,n1,idate(4),idate(2),idate(3),idate(1),thour)
     endif
#endif
!
!  write surface file
!
     if(mod(kdt,max(ksfc,1)).eq.0.or.lastep) then
       call file_name('sfc',3,thour,fno,ncho)
       call sfc_read_file_driver(n1,fno,sfcftyp,                               &
                  labs,idate(4),idate(2),idate(3),idate(1),thour,              &
                  sfcfcs,LONF2S,LATG2S,1)
     endif
#endif            /* BARO_TEST */
   endif nstep1
#ifndef BARO_TEST
!
!  write restart files
!
   if(mod(kdt,max(krestart,1)).eq.0.or.lastep) then
#ifdef DFS
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,sl,si,sf1,tavexy,   &
#else
     call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,ak5,bk5,sf1,tavexy, &
#endif
#else /* SPH */
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,sl,si,gz,z00,      &
#else
     call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,    &
#endif
#endif
#ifdef SMP
               vvel,                                                           &
#endif
               4)
#ifdef DFS
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,psl2,t2,d2,v2,q2,sl,si,sf1,tavexy,   &
#else
     call file_write_sigma(n1,thour,idate,psl2,t2,d2,v2,q2,ak5,bk5,sf1,tavexy, &
#endif
#else
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,q ,te ,di ,ze ,rq,sl,si,gz,z00,      &
#else
     call file_write_sigma(n1,thour,idate,q ,te ,di ,ze ,rq,ak5,bk5,gz,z00,    &
#endif
#endif
#ifdef SMP
               vvel,                                                           &
#endif
               5)
#ifndef NO_PHYSICS
     fno='sfci'
     call sfc_read_file_driver(n1,fno,sfcftyp,                                 &
                  labs,idate(4),idate(2),idate(3),idate(1),thour,              &
                  sfcfcs,LONF2S,LATG2S,1)
#ifdef RIVER
!
!  river file
!
     if (iope) then
       close(42)
       open(42,file='rivi',form='unformatted')
       print*,'writing rivi'
       write(42)gdriv
       close(42)
     endif
#endif
#endif
   endif
!
!  write diagnostic files
!
   if(mod(kdt,max(kpost,1)).eq.0.or.lastep) then
! not used
!        if( stepone ) then
!          xhour=shour/3600.
!        else
!          xhour=thour
!        endif
#if !defined(MRG_POST) && !defined(NO_PHYSICS)
     call file_write_flux(n1)
#elif !defined(NO_PHYSICS)
#ifndef HYBRID
     call post_pgb_2d(thour,idate,sl,si,                                    &
#else
     call post_pgb_2d(thour,idate,ak5,bk5,                                  &
#endif
                    gz,q,te,di,ze,rq,colrad,dummy,n1)
#endif
#ifdef SAS_DIAG
     do i = 1,LONF2S
       do j = 1,LATG2S
         if (cvb(i,j).ne.100.) then
           kbcnv = int(cvb(i,j))
           cvbp(i,j) = sl(kbcnv) * psurf(i,j) * 1e+3
         else
           cvbp(i,j) =  -999.0
         endif
         if (cvt(i,j).ne.0.) then
           ktcnv = int(cvt(i,j))
           cvtp(i,j) = sl(ktcnv) * psurf(i,j) * 1e+3
         else
           cvtp(i,j) =  -999.0
         endif
       enddo
     enddo
     call file_write_cps_bin(thour,LONF2S,LATG2S,levs_,                        &
                   dcu,dcv,dct,dcq,dch,fcu,fcd,deltb,delqb,delhb,              &
                   cbmf,cvbp,cvtp,dlt,dlq,dlh)
#endif
!
!  write dg3 file
!
#ifdef EXPLICIT_CLOUDINESS
#ifndef KSAS
     call file_write_cld_grib (fhour,thour,idate,sl,colrad,                    &
                      fluxr,cvavg,qcicps,qrscps,                               &
                      taucld,cldwp,cldip,                                      &
                      n1,lastep)
#endif
#endif
#ifdef DG3
     call file_write_3d (slmsk,n1)
     if( .not. stepone ) call diag_3d_zero_out
#endif
#ifdef DGP
!
! write DGP
!
#ifdef MP
     call mpgetken
#endif
     call file_write_point(npoint,ikfreq,imodk,itnum,svdata,                   &
                 lab,thour,idate,si,sl,nvrken,nptken,nstken,n1)
#endif
     if( .not. stepone ) then
         call dyn_flux_zero_out(0)
         dtpost=0.
     endif
     call dyn_tmax_zero_out
   endif
#endif /* BARO_TEST end */
!
#ifdef LFM
   if(mod(kdt,max(klfm,1)).eq.0.or.lastep) then
      call file_write_lfm(ifstep,thour,n1)
   endif
#endif
#ifdef RMP
!
!  save RMP out
!
   if(mod(kdt,max(krsm,1)).eq.0.or.lastep) then
      call rmp_end(fhour)
   endif
#endif
!
!  reset timers
!
   if(.not.stepone ) then
     if(mod(kdt,max(min(ksfcx,kpost),1)).eq.0) then
       fhour=thour
       solhr=fhour+idate(1)
       iday=solhr/24.e0
       solhr=solhr-iday*24.e0
       solsec=solhr*3600.
       shour=0.
#ifndef NOPRINT
       if(iope) write(6,*)'reset fhour=',fhour
#endif
     endif
   endif
!
! end of stepone
!
   if( lastep .and. stepone ) then
!
! re-write sig.ft00 with added tracers
!
#ifdef DFS
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,sl,si,sf1,tavexy,   &
#else
     call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,ak5,bk5,sf1,tavexy, &
#endif
#else /* SPH */
#ifndef HYBRID
     call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,sl,si,gz,z00,      &
#else
     call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,    &
#endif
#endif /* DFS end */
#ifdef SMP
               vvel,                                                           &
#endif
               1)
!
     stepone=.false.
     solsec=solsec+deltim
     solhr=solsec/3600.0
     deltim=con(1)
#ifdef DG3
     call diag_3d_zero_out
#endif
#ifdef DGP
     call dyn_tmax_zero_out
#endif
#ifdef LFM
     ifstep=2
     call dyn_low_freq_model(ifstep,fhour)
     ipstep=1
#endif
#ifdef DBG
     if(iope) print *,' forward step done. '
#ifdef DFS
     call dyn_comp_rms(psl2,d2,t2,v2,del,q2)
#else
     call dyn_comp_rms(q, di, te, ze, del, rq)
#endif
#endif /* DBG end */
!
   endif
!
!  end  time  step
!
   return
   end subroutine gmp_integrate
!
!-------------------------------------------------------------------------------
   subroutine dyn_partial_time_filter(tem,te,dim,di,zem,ze,rm,rq,qm,q,         &
                                    filta,lnt22,lnt2,levs,levh)              
#ifdef NISLQ
#ifdef DFS
   use dfsvar, only : ib,jbw
#else
   use paramodel, only : LONF2S,LATG2S
#endif
#endif /* NISLQ end */
#if defined(RAW) || defined(HORA)
   use comdyn, only :   workt,workd,workz,workr,coef1,coef2,coef3
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                 ::  lnt22,lnt2,levs,levh,i,j,k
   real                                    ::  filta,filtb
   real   , dimension(lnt22)               ::  qm            , q
   real   , dimension(lnt22,levs)          ::  tem, dim, zem , te , di , ze
#ifdef NISLQ
#ifdef DFS
   real   , dimension(ib    ,jbw ,levh  )  ::  rm            , rq
#else
   real   , dimension(LONF2S,levh,LATG2S)  ::  rm            , rq
#endif
#else
   real   , dimension(lnt22,levh)          ::  rm            , rq
#endif /* NISLQ end */
#ifdef HORA
   real   , dimension(lnt22,levs)          ::  wkt, wkd, wkz
   real   , dimension(lnt22,levh)          ::  wkr         
#endif
!
! filter coefficient
!
#if defined(RAW)
   filta = coef2
   filtb = coef1
#elif defined(HORA)
   filta = coef3
   filtb = coef2
#else
   filtb = (1.-filta)*0.5                                                    
#endif
#ifdef HORA
!
! filtered value at n-2 (high-order time filter)
!
   forall(j=1:lnt2,k=1:levs)
     wkt(j,k)   = workt(j,k)
     wkd(j,k)   = workd(j,k)
     wkz(j,k)   = workz(j,k)
   end forall
   forall(j=1:lnt2,k=1:levh)
     wkr(j,k)   = workr(j,k)
   end forall
  !
   forall(j=1:lnt2,k=1:levs)
     workt(j,k) = tem(j,k)
     workd(j,k) = dim(j,k)
     workz(j,k) = zem(j,k)
   end forall
   forall(j=1:lnt2,k=1:levh)
     workr(j,k) = rm(j,k)
   end forall
#endif
!
! surface pressure
!
   do j = 1,lnt2
     qm(j)=q(j)
   enddo
!
! temperature, divergence, vorticity
!
#ifdef RAW
   ! raw time filter (n+1)
   forall(j=1:lnt2,k=1:levs)
     workt(j,k)=coef3*(0.5*tem(j,k)-te(j,k))
     workd(j,k)=coef3*(0.5*dim(j,k)-di(j,k))
     workz(j,k)=coef3*(0.5*zem(j,k)-ze(j,k))
   end forall
#endif
   do k = 1,levs
     do j = 1,lnt2
       tem(j,k)=filtb*tem(j,k)+filta*te(j,k)                                     
       dim(j,k)=filtb*dim(j,k)+filta*di(j,k)                                     
       zem(j,k)=filtb*zem(j,k)+filta*ze(j,k)                                     
     enddo                                                                    
   enddo
#ifdef HORA
   ! high-order ra time filter
   forall(j=1:lnt2,k=1:levs)
     tem(j,k)=tem(j,k)+coef1*wkt(j,k)
     dim(j,k)=dim(j,k)+coef1*wkd(j,k)
     zem(j,k)=zem(j,k)+coef1*wkz(j,k)
   end forall
#endif
!
! moisture
!
#ifdef NISLQ
#ifdef DFS
   forall(i=1:ib,j=1:jbw,k=1:levh) rm(i,j,k)=filtb*rm(i,j,k)+filta*rq(i,j,k)
#else
   forall(i=1:LONF2S,k=1:levh,j=1:LATG2S) rm(i,k,j)=filtb*rm(i,k,j)+filta*rq(i,k,j)
#endif
#else /* ~NISLQ */
#ifdef RAW
   ! raw time filter (n+1)
   forall(j=1:lnt2,k=1:levh)
     workr(j,k)=coef3*(0.5*rm(j,k)-rq(j,k))
   end forall
#endif
   do k = 1,levh                                                         
     do j = 1,lnt2                                                         
       rm(j,k)=filtb*rm(j,k)+filta*rq(j,k)                                       
     enddo                                                               
   enddo
#ifdef HORA
   ! high-order ra time filter
   forall(j=1:lnt2,k=1:levh) 
     rm(j,k)=rm(j,k)+coef1*wkr(j,k)                                       
   end forall
#endif
#endif /* NISLQ end */
!
   return                                                                    
   end subroutine dyn_partial_time_filter
!
!-------------------------------------------------------------------------------
   subroutine dyn_full_time_filter(tem,te,dim,di,zem,ze,rm,rq,                 &
#ifdef AAA
                     qm,q,                                                     &
#endif
                     filta,lnt22,lnt2,levs,levh)              
#ifdef NISLQ
#ifdef DFS
   use dfsvar, only : ib,jbw
#else
   use paramodel, only : LONF2S,LATG2S
#endif
#endif /* NISLQ end */
#if defined(RAW) || defined(HORA)
   use comdyn, only : workt,workd,workz,workr,coef4
#ifdef RAW
   use comdyn, only : coef1
#endif
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                 ::  lnt22,lnt2,levs,levh,i,j,k
   real                                    ::  filta,filtb
#ifdef DFS
   real   , dimension(lnt22)               ::  qm            , q
#endif
   real   , dimension(lnt22,levs)          ::  tem, dim, zem , te , di , ze
#ifdef NISLQ
#ifdef DFS
   real   , dimension(ib    ,jbw ,levh  )  ::  rm            , rq
#else
   real   , dimension(LONF2S,levh,LATG2S)  ::  rm            , rq
#endif
#else
   real   , dimension(lnt22,levh)          ::  rm            , rq
#endif
!
! filter coefficient
!
#if defined(RAW)
   filtb = coef1
#elif defined(HORA)
   filtb = coef4
#else
   filtb = (1.-filta)* 0.5                                                   
#endif /* RAW end */
!
#ifdef AAA
   do  j = 1,lnt2
     qm(j)=qm(j)+filtb*q(j)
   enddo
#endif
!
! temperature, divergence, vorticity
!
   do k = 1,levs                                                         
     do j = 1,lnt2                                                         
       tem(j,k)=tem(j,k)+filtb*te(j,k)                                           
       dim(j,k)=dim(j,k)+filtb*di(j,k)                                           
       zem(j,k)=zem(j,k)+filtb*ze(j,k)                                           
     enddo                                                                    
   enddo
#ifdef RAW
   ! raw time filter (n+1)
   forall(j=1:lnt2,k=1:levs)
     te(j,k)=workt(j,k)+coef4*te(j,k)
     di(j,k)=workd(j,k)+coef4*di(j,k)
     ze(j,k)=workz(j,k)+coef4*ze(j,k)
   end forall
#endif
!
! moisture
!
#ifdef NISLQ
#ifdef DFS
   forall(i=1:ib,j=1:jbw,k=1:levh) rm(i,j,k)=rm(i,j,k)+filtb*rq(i,j,k)
#else
   forall(i=1:LONF2S,k=1:levh,j=1:LATG2S) rm(i,k,j)=rm(i,k,j)+filtb*rq(i,k,j)
#endif
#else /* ~NISLQ end */
   do k = 1,levh                                                         
     do j = 1,lnt2                                                         
       rm(j,k)=rm(j,k)+filtb*rq(j,k)                                             
     enddo                                                                    
   enddo
#ifdef RAW
   ! raw time filter (n+1)
   forall(j=1:lnt2,k=1:levh)
     rq(j,k)=workr(j,k)+coef4*rq(j,k)
   end forall
#endif
#endif /* ~NISLQ end */
!
   return                                                                    
   end subroutine dyn_full_time_filter
!
!--------------------------------------------------------------------------------
   subroutine rdscmbd
!--------------------------------------------------------------------------------
#ifdef SMP
   use paramodel, only : JCAP1S,LNT2S,LATG2S,LEVSS,LEVHS,LCAPS,LCAP22S,         &
                         levs_,levh_,lnt22_,lnt2_,ncpus_
#ifndef DYNAMIC_ALLOC
   use paramodel, only : ncpus=>ncpus_
#endif
   use constant , only : cp_,g_,rd_,rv_,rerth_,hfus_,hvap_,qmin_
   use varsfc   , only : msub_,nsoil_,lsoil_,lalbd_
   use comfibm
   use comfcst  , only : dtbdy,curtime,vvel,hour1,hour2,hours,houre,ftim1,ftim2
#ifdef CLM_CWF
   use comfcst  , only : wdiv, hadq
#endif
#ifdef SMP_NUDGING_SFC
#ifdef SMP_NUDGING_FLX
   use comfcst, only : scmsh, scmlh
#endif
#ifdef SMP_NUDGING_RAD
   use comfcst, only : psfc,tsfc,qsfc,prec,scmlwup,scmlwdn,scmswup,scmswdn
#endif
#endif  /* SMP_NUDGING_SFC end */
   use comfcst, only : ps1,uo1,vo1,to1,qo1,wo1     ,&
                       ps2,uo2,vo2,to2,qo2,wo2     ,&
#ifdef REV_FRC
                       ac1,ac2                     ,&
#endif
#ifdef RLX_FRC
                       tscl1, tscl2                ,&
#endif
                       at1,aq1,dv1                 ,&
                       at2,aq2,dv2
!--------------------------------------------------------------------------------
#include "abort.h"
!
#ifdef SMP_NUDGING_SFC
#ifdef SMP_NUDGING_FLX
   real,save,allocatable::  scmsh1(:),scmlh1(:)
   real,save,allocatable::  scmsh2(:),scmlh2(:)
#endif
#ifdef SMP_NUDGING_RAD
   real,save,allocatable::  psfc1(:),tsfc1(:)
   real,save,allocatable::  qsfc1(:),prec1(:)
   real,save,allocatable::  psfc2(:),tsfc2(:)
   real,save,allocatable::  qsfc2(:),prec2(:)
!
   real,save,allocatable::  scmlwup1(:),scmlwdn1(:)
   real,save,allocatable::  scmswup1(:),scmswdn1(:)
   real,save,allocatable::  scmlwup2(:),scmlwdn2(:)
   real,save,allocatable::  scmswup2(:),scmswdn2(:)
#endif
   integer              ::  iy,im,id,ih
#endif  /* SMP_NUDGING_SFC end */
!
#ifdef DBG
#ifndef CLM_CWF
   real                 ::  hadq(lnt22_,levs_)
#endif
   real                 ::  hadt(lnt22_,levs_)
   real                 ::  vadt(lnt22_,levs_), vadq(lnt22_,levh_)
   real                 ::  adbc(lnt22_,levs_)
   real                 ::  rlxt(lnt22_,levs_), rlxq(lnt22_,levh_)
#endif
   real                 ::  time_scm
   character(len=80)    ::  fno
   integer              ::  idateg(4)
   logical,save         ::  first
   data first/.true./
!
#ifdef SMP_NUDGING_SFC
#ifdef SMP_NUDGING_FLX
   if(.not.allocated(scmsh1))  allocate( scmsh1(lnt22_) )
   if(.not.allocated(scmlh1))  allocate( scmlh1(lnt22_) )
   if(.not.allocated(scmsh2))  allocate( scmsh2(lnt22_) )
   if(.not.allocated(scmlh2))  allocate( scmlh2(lnt22_) )
#endif
#ifdef SMP_NUDGING_RAD
   if(.not.allocated(psfc1))  allocate( psfc1(lnt22_) )
   if(.not.allocated(tsfc1))  allocate( tsfc1(lnt22_) )
   if(.not.allocated(qsfc1))  allocate( qsfc1(lnt22_) )
   if(.not.allocated(prec1))  allocate( prec1(lnt22_) )
   if(.not.allocated(psfc2))  allocate( psfc2(lnt22_) )
   if(.not.allocated(tsfc2))  allocate( tsfc2(lnt22_) )
   if(.not.allocated(qsfc2))  allocate( qsfc2(lnt22_) )
   if(.not.allocated(prec2))  allocate( prec2(lnt22_) )
!
   if(.not.allocated(scmlwup1))  allocate( scmlwup1(lnt22_) )
   if(.not.allocated(scmlwdn1))  allocate( scmlwdn1(lnt22_) )
   if(.not.allocated(scmswup1))  allocate( scmswup1(lnt22_) )
   if(.not.allocated(scmswdn1))  allocate( scmswdn1(lnt22_) )
   if(.not.allocated(scmlwup2))  allocate( scmlwup2(lnt22_) )
   if(.not.allocated(scmlwdn2))  allocate( scmlwdn2(lnt22_) )
   if(.not.allocated(scmswup2))  allocate( scmswup2(lnt22_) )
   if(.not.allocated(scmswdn2))  allocate( scmswdn2(lnt22_) )
#endif
#endif  /* SMP_NUDGING_SFC end */
!=========================================================================
! SMP : read dynamic tendency terms (advection)
!=========================================================================
   if (first) then
     print*,'rdscmbd first'
     hour1 = thour
     hour2 = hour1 + dtbdy
   endif
!  fv = rv_/rd_ - 1.0
!
   time_scm = curtime
   write(6,*) 'rdscmbd check time'
   write(6,'(A,4f8.2)') 'fhour shour thour curtime',fhour,shour,thour,curtime
   write(6,'(A,2f8.2)') 'h1  h2', hour1,hour2
!
   if (time_scm.eq.hour1) then
     nng = 31
     call file_name('scmbdy',6,thour,fno,ncho)
     open(nng,file=fno(1:ncho),form='unformatted',status='old')
     read(nng) thourg, idateg, dtv
     write(6,*) 'SCMBDY rdscmbd ',fno(1:ncho)
#ifdef REV_FRC
     read(nng) uo1, vo1, at1, ac1, aq1, wo1, dv1, ps1
     read(nng) uo2, vo2, at2, ac2, aq2, wo2, dv2, ps2
#else
#ifdef RLX_FRC
#ifndef TOGA
     read(nng) uo1, vo1, to1, qo1, at1, aq1, wo1, dv1, ps1, tscl1
     read(nng) uo2, vo2, to2, qo2, at2, aq2, wo2, dv2, ps2, tscl2
#else
     read(nng) uo1, vo1, to1, qo1, at1, aq1, wo1, dv1, ps1
     read(nng) uo2, vo2, to2, qo2, at2, aq2, wo2, dv2, ps2
#endif
#else
     read(nng) uo1, vo1, at1, aq1, wo1, dv1, ps1
     read(nng) uo2, vo2, at2, aq2, wo2, dv2, ps2
#endif
#endif
     close(nng)
#ifdef SMP_NUDGING_SFC
!
     nnf = 32
     call file_name('SMP_FLX',7,thour,fno,ncho)
     open(nnf,file=fno(1:ncho),form='unformatted',status='old')
     read(nnf) iy,im,id,ih
     write(6,*) 'SMP_NUDGING_FLX rdscmbd ',fno(1:ncho)
     write(6,'(4i5)') iy,im,id,ih
#ifdef SMP_NUDGING_FLX
     read(nnf) scmsh1,scmlh1
     read(nnf) scmsh2,scmlh2
     write(6,'(a11,2f8.3)') ' SH   LH   ',scmsh1(1),scmlh1(1)
#endif
#ifdef SMP_NUDGING_RAD
     read(nnf) psfc1,tsfc1,qsfc1,prec1,scmlwup1, scmlwdn1, scmswup1, scmswdn1
     read(nnf) psfc2,tsfc2,qsfc2,prec2,scmlwup2, scmlwdn2, scmswup2, scmswdn2
     print*,'rdscmbd lnt22',lnt22_
     write(6,'(a11,2f8.3)') ' psfc', psfc1(1),psfc2(1)
     write(6,'(a11,2f8.3)') ' tsfc', tsfc1(1),tsfc2(1)
     write(6,'(a11,2f8.3)') ' qsfc', qsfc1(1),qsfc2(1)
     write(6,'(a11,2f8.3)') ' prec', prec1(1),prec2(1)
     write(6,'(a11,2f8.3)') ' lwup', scmlwup1(1),scmlwup2(1)
     write(6,'(a11,2f8.3)') ' lwdn', scmlwdn1(1),scmlwdn2(1)
     write(6,'(a11,2f8.3)') ' swup', scmswup1(1),scmswup2(1)
     write(6,'(a11,2f8.3)') ' swdn', scmswdn1(1),scmswdn2(1)

#endif
     close(nnf)
#endif
!
! ... check boundary ...
!
     do k = 1,levs_
       if (k.lt.3) then
         write(6,'(a11,i2)') 'rdscmbd k =',k
         write(6,220) ' uo vo to qo ',uo1(1,k),vo1(1,k),to1(1,k),qo1(1,k)
         write(6,220) ' at aq wo dv ',at1(1,k),aq1(1,k),wo1(1,k),dv1(1,k)
         write(6,230) ' ps          ',ps1(1)
#ifdef RLX_FRC
         write(6,230) ' tscl        ',tscl1(1,k)
#endif
       endif
     enddo
220  format(a13,4e15.7)
230  format(a13,2e15.7)
     hours = hour1
     houre = hour2
     hour1 = hour2
     hour2 = hour1 + dtbdy
   endif
   ftim1 = (time_scm-hours)/dtbdy
   ftim2 = (houre-time_scm)/dtbdy
   write(6,'(A,1f8.2)')'rdscmbd cur thour ', time_scm
   write(6,'(A,2f8.2)')'rdscmbd h1  h2    ', hour1,hour2
   write(6,'(A,2f8.2)')'rdscmbd hs  he    ', hours,houre
   write(6,'(A,2f8.2)')'rdscmbd ft1 ft2   ', ftim1,ftim2
   write(6,'(A,3f8.2)')'rdscmbd', real(kdt),deltim,dtbdy
!
! ... set surface values ...
!
#ifdef SMP_NUDGING_SFC
   do l = 1, lnt2_
#ifdef SMP_NUDGING_FLX
! no interpolation
     if (time_scm.eq.0) then
       scmsh(l) = scmsh1(l)
       scmlh(l) = scmlh1(l)
     else
       scmsh(l) = scmsh2(l)
       scmlh(l) = scmlh2(l)
     endif
! temporal interpolation
!     scmsh(l) = ftim2*scmsh1(l) + ftim1*scmsh2(l)
!     scmlh(l) = ftim2*scmlh1(l) + ftim1*scmlh2(l)
     write(6,'(a11,2f8.3)') ' SH   LH   ',scmsh(1),scmlh(1)
#endif
#ifdef SMP_NUDGING_RAD
!no interpolation
!     if (time_scm.eq.0.or.deltim.lt.600.) then
!       psfc(l) = psfc1(l)
!       tsfc(l) = tsfc1(l)
!       qsfc(l) = qsfc1(l)
!     else
!       psfc(l) = psfc2(l)
!       tsfc(l) = tsfc2(l)
!       qsfc(l) = qsfc2(l)
!     endif
!temporal interpolation
     psfc(l) = ftim2*psfc1(l) + ftim1*psfc2(l)
     tsfc(l) = ftim2*tsfc1(l) + ftim1*tsfc2(l)
     qsfc(l) = ftim2*qsfc1(l) + ftim1*qsfc2(l)
! 
     prec(l) = (deltim)*(prec2(l)/1000./3600.)
!
     if (l.eq.1) then
       write(6,*) 'SMP_NUDGING_FLX temporal interpolation'
       write(6,'(a11,2f8.3)') ' Psfc tsfc ',psfc(1),tsfc(1)
       write(6,'(a11,2f8.5)') ' qsfc prec ',qsfc(1),prec(1)
     endif
!no interpolation
! dt=600 sec
!     if (time_scm.eq.0.or.deltim.lt.600.) then
! dt=3 hr
!     if (time_scm.eq.0.or.deltim.lt.10800.) then
!       scmlwup(l) = scmlwup1(l)
!       scmlwdn(l) = scmlwdn1(l)
!       scmswup(l) = scmswup1(l)
!       scmswdn(l) = scmswdn1(l)
!     else
!       scmlwup(l) = scmlwup2(l)
!       scmlwdn(l) = scmlwdn2(l)
!       scmswup(l) = scmswup2(l)
!       scmswdn(l) = scmswdn2(l)
!     endif
!temporal interpolation
     scmlwup(l) = ftim2*scmlwup1(l) + ftim1*scmlwup2(l)
     scmlwdn(l) = ftim2*scmlwdn1(l) + ftim1*scmlwdn2(l)
     scmswup(l) = ftim2*scmswup1(l) + ftim1*scmswup2(l)
     scmswdn(l) = ftim2*scmswdn1(l) + ftim1*scmswdn2(l)
#endif
   enddo
!
   write(6,*) 'SMP_NUDGING_RAD temporal interpolation'
   write(6,'(a11,6f8.3)') ' lwup lwdn ',scmlwup,scmlwdn
   write(6,'(a11,6f8.3)') ' swup swdn ',scmswup,scmswdn
#endif
!
   first=.false.
!
   return
#endif
   end subroutine rdscmbd
