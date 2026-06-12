#include <define.h>
   subroutine phys_rad_solver(lons2,lats2,                                     &
#ifndef NONHYD
                     uugrs,vvgrs,plgr,tvgrs,rqgrs,xxgrs,                       &
#ifndef SMP
                     pllamgr,plphigr,                                          &
#endif
#else /* else of NONHYD */
                     plgr, tvgrs, rqgrs, wwgrs,                                &
#ifdef NIM_DIAG
                     o3_nim,rflux,                                             &
#endif
#endif /* endif NONHYD */
#ifdef DG3
                     gda,                                                      &
#endif
#ifdef DG
                     cldt,clcv,                                                &
#endif
#ifdef EXPLICIT_CLOUDINESS
                     qcicps,qrscps,                                            &
                     taucld,cldwp,cldip,                                       &
#endif
#ifdef NIM
                     ips,ipe,lat,latrue,nimlat,nimlon,si1d,sl1d,prsi_nim,prsl_nim)
#else
                     lat,latrue)
#endif
!-------------------------------------------------------------------------------
!
! subroutine: phys_rad_solver
!
! program history log:
!   1993-09-01  ken campana            mrf implementation
!   1998-01-01  yu-tai ho              chou's scheme
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!  RMP variables:
!
!    ddgrd = xxgrs
!    pxgrd = pllamgr
!    pygrd = plphigr
!    uugrd = uugrs
!    vvgrd = vvgrs
!    ttgrd = tvgrs
!    qqgrd = rqgrs
!    ppgrd = plgr 
!    
!-------------------------------------------------------------------------------
   use paramodel, only      : LONF2S,LONF22S,LATG2S,                           &
                              ngases_,ntotal_,nwmass_,npes_,levs_,levh_,levp1_,&
                              kgases_,ncloud_,icloud_
   use comsfc
   use comio
#ifdef SWRMDC
   use aerparm
#endif
#ifndef RMP
   use paramodel, only      : LATGS
   use comznl
   use comfgrid
   use comfphys
   use comfver
   use comgrad
   use radiag
#ifdef MP
#ifndef DFS
   use commpi
#endif
#endif
#else   /* RMP */
#ifdef MP
   use paramodel, only      : igrd12p_,jgrd12p_
   use commpi
#else
   use paramodel, only      : igrd12_
#endif
   use rscomf_rerun
   use rscomltb
   use rscommap
   use rscomgrad
#endif  /* RMP end */
   use constant, only       : pi_, qmin_
#ifdef DFS
   use dfsvar, only         : iope
#ifdef MP
   use commpi, only         : mype
#endif
#endif
#ifdef DG3
   use comgda
   use diag_3d_module, only : diag_3d_archive
#endif
!-------------------------------------------------------------------------------
#ifndef RMP
   implicit none
!-------------------------------------------------------------------------------
#endif
#include <abort.h>
!
#ifndef RMP
#define LONLENS lons2
#else /* RMP */
#ifdef MP
#ifdef RMPVECTORIZE
#define LONLENS igrd12p_*jgrd12p_
#else
#define LONLENS lonlen(mype)*2
#endif
#else /* ~MP */
#undef RMPVECTORIZE
#define LONLENS igrd12_
#endif /* MP end */
#endif /* ~RMP end */
   integer              ::  lat,latrue
#ifndef NONHYD
   real                 ::  pllamgr(LONF22S)
   real                 ::  plphigr(LONF22S)
   real                 ::  uugrs(LONF22S,levs_)
   real                 ::  vvgrs(LONF22S,levs_)
   real                 ::  xxgrs(LONF22S,levs_)
#endif /* NONHYD */
   real                 ::  tvgrs(LONF22S,levs_)
   real                 ::  rqgrs(LONF22S,levh_)
   real                 ::  plgr(LONF22S)
#ifdef NONHYD
   real                 ::  wwgrs(LONF22S,levp1_)
#endif /* NONHYD */
!
! local array
!
#ifndef NONHYD
   real                 ::  plamgr(LONF2S)
   real                 ::  pphigr(LONF2S)
   real                 ::  ugrs(LONF2S,levs_)
   real                 ::  vgrs(LONF2S,levs_)
   real                 ::  xgrs(LONF2S,levs_)
#endif /* NONHYD */
   real                 ::  tgrs(LONF2S,levs_)
   real                 ::  qgrs(LONF2S,levh_)
   real                 ::  pgr(LONF2S)
   real                 ::  qcis(LONF2S,levs_)
#ifdef NONHYD
   real                 ::  wgrs(LONF2S,levp1_)
#endif /* NONHYD */
#ifdef NIM
   real                 ::  nimlat(LONF22S),nimlon(LONF22S)
   real                 ::  prsi_nim(LONF22S,levp1_),prsl_nim(LONF22S,levs_)
   real                 ::  si1d(levp1_),sl1d(levs_)
#ifdef NIM_DIAG
   real                 ::  o3_nim(LONF22S,levs_)
   real                 ::  rflux(LONF22S,6)
#endif
!
   real                 ::  fnet(LONF2S)
#endif
!
#ifdef DG3
   real                 ::  gda(nwgda,kdgda)
#endif
#ifdef DG
   real                 ::  cldt(LONF2S,levs_),clcv(LONF2S,levs_)
#endif
#ifdef EXPLICIT_CLOUDINESS
   real                 ::  qcicps(LONF2S,levs_), qrscps(LONF2S,levs_)
   real                 ::  taucld(LONF2S,levs_)
   real                 ::  cldwp(LONF2S,levs_), cldip(LONF2S,levs_)
   real                 ::  qcis_cps(LONF2S,levs_)
#endif
!
   real                 ::  workr(LONF2S,levs_)
   real                 ::  w2(LONF2S,levs_)
   real                 ::  tsear(LONF2S),shelgr(LONF2S)
   real                 ::  albdoa(LONF2S)
   real                 ::  cldary(LONF2S,levs_)
   real                 ::  cldsa(LONF2S,4)
   real                 ::  cldtot(LONF2S,levs_)
   real                 ::  cldcnv(LONF2S,levs_)
   integer              ::  mtopa(LONF2S,3)
   integer              ::  mbota(LONF2S,3)
   real                 ::  swhr(LONF2S,levs_)
   real                 ::  hlwr(LONF2S,levs_)
   real                 ::  sfnswr(LONF2S),sfdlwr(LONF2S)
#ifdef VIC
   real                 ::  sfdswr(LONF2S),sfuswr(LONF2S)
#endif
   real                 ::  tsflwr(LONF2S)
   real                 ::  tgr(LONF2S)
   real                 ::  gdfvbr(LONF2S),gdfnbr(LONF2S)
   real                 ::  gdfvdr(LONF2S),gdfndr(LONF2S)
#ifdef SWRMDC
!
! new aerosol scheme
!
   real                 ::  rhrh(LONF2S,levs_)
   integer              ::  idxc(nxc,LONF2S)
   real                 ::  cmix(nxc,LONF2S)
   real                 ::  denn(ndn,LONF2S)
   integer              ::  kprf(LONF2S)
   real                 ::  rdg,hdlt,tmp1,tmp2
   integer              ::  i2,j2,i1,j1
#endif
#ifdef NIM
   integer              ::  ips,ipe
#endif
   integer              ::  lons2,lats2,latco
   integer              ::  i,k,ic,kc
#ifdef DBG
   integer              ::  numlevs
#endif
#ifndef MP
   integer,parameter    ::  mype=0
#endif
!
! initialize local variables
!
   tgrs=0.   ; qgrs=0.   ; pgr=0.     ; qcis=0.
#ifndef NONHYD
   plamgr=0. ; pphigr=0. ; ugrs=0.    ; vgrs=0.   ; xgrs=0.
#else
   wgrs=0.
#endif /* ~NONHYD end */
   workr=0.  ; w2=0.     ; tsear=0.   ; shelgr=0. ; albdoa=0.
   cldary=0. ; cldsa=0.  ; cldtot=0.  ; cldcnv=0. ; mtopa=0.
   mbota=0.  ; swhr=0.   ; hlwr=0.    ; sfnswr=0. ; sfdlwr=0.
#ifdef VIC
   sfdswr=0. ; sfuswr=0.
#endif
   tsflwr=0. ; tgr=0.
   gdfvbr=0. ; gdfnbr=0. ; gdfvdr=0.  ; gdfndr=0.
#ifdef SWRMDC
   rhrh=0.   ; idxc=0    ; cmix=0.    ; denn=0.   ;  kprf=0
#endif
#ifdef EXPLICIT_CLOUDINESS
   qcis_cps=0.
#endif
!
#ifdef DG
   dtacc = min(dtswav,dtlwav)
   dtacc = dtacc*3600.e0
#endif
!
#ifndef RMP
#ifdef NIM
   latco=1     ! temporary value
#else
   latco=LATGS+1-lat
#endif
#endif
!
   do i = 1,LONLENS
#ifndef RMP
#ifdef NIM
     pgr(i) = plgr(i)/1000.
#else
     pgr(i) = exp(plgr(i))
#endif
#else
     pgr(i) = plgr(i)
#endif
#ifndef NONHYD
#ifndef SMP
     plamgr(i)=pllamgr(i)
     pphigr(i)=plphigr(i)
#endif
#endif /* NONHYD */
   enddo
!....
   do k = 1,levs_
     do i = 1,LONLENS
#ifndef NONHYD
       ugrs(i,k)=uugrs(i,k)
       vgrs(i,k)=vvgrs(i,k)
       xgrs(i,k)=xxgrs(i,k)
#endif /* NONHYD */
       tgrs(i,k)=tvgrs(i,k)
     enddo
   enddo
!
   do k = 1,levh_
     do i = 1,LONLENS
       qgrs(i,k)=rqgrs(i,k)
     enddo
   enddo
!
#ifdef NONHYD
  do k = 1,levp1_
    do i = 1,LONLENS
      wgrs(i,k)=wwgrs(i,k)
    enddo
  enddo
#endif
!
!    convert virt. temp to thermodynamic temp.
!
   do k = 1,levs_
     do i = 1,LONLENS
       if(qgrs(i,k).le.0.0) qgrs(i,k)=1.0e-10
       w2(i,k)=1.+0.6*qgrs(i,k)
       qcis(i,k) = 0.
       if(ncloud_.ge.1) then
         do ic = icloud_,nwmass_
           kc = (ic-1)*levs_ + k
           qcis(i,k) = max(qgrs(i,kc),qmin_) + qcis(i,k)
         enddo
         w2(i,k)=w2(i,k) - qcis(i,k)
       endif
#ifdef EXPLICIT_CLOUDINESS
       qcis_cps(i,k) = 0.
       qcis_cps(i,k) = qcicps(i,k) + qrscps(i,k)
#endif
     enddo
   enddo
!
   do k = 1,levs_
     do i = 1,LONLENS
       tgrs(i,k)=tgrs(i,k)/w2(i,k)
     enddo
   enddo
#ifdef NIM
!
! calculate qcis
!
   do k = 1,levs_
     do i = 1,LONLENS
       if(qgrs(i,k).le.0.0) qgrs(i,k)=1.0e-10
       qcis(i,k) = 0.
       if(ncloud_.ge.1) then
         do ic = icloud_,nwmass_
           kc = (ic-1)*levs_ + k
           qcis(i,k) = max(qgrs(i,kc),qmin_) + qcis(i,k)
         enddo
       endif
     enddo
   enddo
#endif
!
   do i = 1,LONLENS
     tsear(i)=tsea(i,lat)
     shelgr(i)=snoweq(i,lat)
     tgr(i)=stc(i,lat,1)
   enddo
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifdef MP
   if (iope) then
#endif
     write(6,*)'maxmin from phys_rad_solver at lat=',lat
#ifndef NONHYD
     call print_maxmin_seven(ugrs,LONLENS,LONF2S,levs_,1,levs_,'ugrs')
     call print_maxmin_seven(vgrs,LONLENS,LONF2S,levs_,1,levs_,'vgrs')
     call print_maxmin_seven(xgrs,LONLENS,LONF2S,levs_,1,levs_,'xgrs')
#endif
     call print_maxmin_seven(tgrs,LONLENS,LONF2S,levs_,1,levs_,'tgrs')
     call print_maxmin_seven(qgrs,LONLENS,LONF2S,levh_,1,levh_,'q')
     call print_maxmin_seven(tsear(1),LONLENS,LONLENS,1,1,1,'tsear')
     call print_maxmin_seven(shelgr(1),LONLENS,LONLENS,1,1,1,'shelgr')
     call print_maxmin_seven(tgr(1),LONLENS,LONLENS,1,1,1,'tgr')
#ifdef RMP
     call print_maxmin_seven(albed (1,lat),LONLENS,LONF2S,1,1,1,'albed')
#endif
#ifdef MP
   endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
   call rad_cloud_prepare(LONLENS,                                             &
#ifndef NONHYD
#ifdef SMP
               xgrs(1,1),                                                      &
#else
               xgrs(1,1),plamgr(1),pphigr(1),                                  &
#endif
               ugrs(1,1),vgrs(1,1),                                            &
#endif /* NONHYD */
               tgrs(1,1),qgrs(1,1),pgr(1),                                     &
#ifdef NONHYD
               wgrs(1,1),                                                      &
#endif /* NONHYD */
#ifndef RMP
               albedr(1,lat),slmskr(1,lat),                                    &
#ifndef NIM
               xlon(1,lat),xlat(1,lat),                                        &
#else
               nimlon(1),nimlat(1),                                            &
#endif
#else
               albed (1,lat),slmsk (1,lat),                                    &
               rlon(1,lat),rlat(1,lat),                                        &
#endif
               tsear(1),shelgr(1),tgr(1),                                      &
               cvr(1,lat),cvtr(1,lat),cvbr(1,lat),rhcl,                        &
               ozon(1,1,lat),albdoa(1),cldary(1,1),                            &
#ifndef SWRMDC
               cldtot(1,1),cldcnv(1,1),                                        &
#else
               cldtot(1,1),cldcnv(1,1),rhrh(1,1),                              &
#endif
               cldsa(1,1),mtopa(1,1),mbota(1,1),                               &
#if defined(RMP) || defined(NIM)
               rrs2,lat,latco,latrue,istrat,                                   &
#else
               rbs2,lat,latco,latrue,istrat,                                   &
#endif
#ifdef NIM
               ips,ipe,si1d(1),sl1d(1),prsi_nim(1,1),prsl_nim(1,1),            &
#endif
               kalb,jo3,slag,rsin1,rcos1,rcos2,                                &
               fjd,dlt,jsno,workr,levs_                                        &
#ifdef DG
               ,cldt(1,1),clcv(1,1)                                            &
#endif
#ifdef EXPLICIT_CLOUDINESS
               ,qcis_cps(1,1)                                                  &
               ,qcicps(1,1),qrscps(1,1)                                        &
#endif
               ,icwp,qcis(1,1))
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifdef MP
   if (iope) then
#endif
     call print_maxmin_seven(slmsk (1,lat),LONLENS,LONF2S,1,1,1,'slmsk')
#ifdef RMP
     call print_maxmin_seven(rlat (1,lat),LONLENS,LONF2S,1,1,1,'rlat')
     call print_maxmin_seven(rlon (1,lat),LONLENS,LONF2S,1,1,1,'rlon')
#else
#ifndef NIM
     call print_maxmin_seven(xlat (1,lat),LONLENS,LONF2S,1,1,1,'rlat')
     call print_maxmin_seven(xlon (1,lat),LONLENS,LONF2S,1,1,1,'rlon')
#endif
#endif
     call print_maxmin_seven(cvr (1,lat),LONLENS,LONF2S,1,1,1,'cvr')
     call print_maxmin_seven(cvtr (1,lat),LONLENS,LONF2S,1,1,1,'cvtr')
     call print_maxmin_seven(cvbr (1,lat),LONLENS,LONF2S,1,1,1,'cvbr')
     call print_maxmin_seven(albdoa(1),LONLENS,LONF2S,1,1,1,' albdoa')
     write(6,*)' after rad_cloud_prepare',LONF2S,                              &
#ifndef NONHYD
#ifdef SMP
               xgrs(1,1),                                                      &
#else
               xgrs(1,1),plamgr(1),pphigr(1),                                  &
#endif
               ugrs(1,1),vgrs(1,1),                                            &
#endif /* NONHYD */
               tgrs(1,1),qgrs(1,1),pgr(1),                                     &
#ifdef NONHYD
               wgrs(1,1),                                                      &
#endif /* NONHYD */
#ifdef RMP
               albed (1,lat),slmsk (1,lat),                                    &
               rlon(1,lat),rlat(1,lat),                                        &
#else
#ifdef NIM
               albedr (1,lat),slmsk (1,lat),                                   &
               nimlon(1),nimlat(1),                                            &
#else
               albedr (1,lat),slmsk (1,lat),                                   &
               xlon(1,lat),xlat(1,lat),                                        &
#endif
#endif
               tsear(1),shelgr(1),tgr(1),                                      &
               cvr(1,lat),cvtr(1,lat),cvbr(1,lat),rhcl(1,1,1,1,1),             &
               ozon(1,1,lat),albdoa(1),cldary(1,1),                            &
               cldtot(1,1),cldcnv(1,1),                                        &
               cldsa(1,1),mtopa(1,1),mbota(1,1),                               &
#if defined(RMP) || defined(NIM)
               rrs2,                                                           &
#else
               rbs2,                                                           &
#endif
               lat,latco,latrue,istrat,                                        &
               kalb,jo3,slag,rsin1,rcos1,rcos2,                                &
               fjd,dlt,jsno,workr(1,1),levs_
#ifdef MP
   endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif /* DBG end */
!
!  use forecasted ozone if possible
!
    if (ngases_.ge.1.and.thour.gt.0.) then
     do k = 1,levs_
       do i = 1,LONLENS
#ifdef RMPVECTORUZE
         ozon(i,lat,k) = qgrs(i,k+kgases_-1)
#else
         ozon(i,k,lat) = qgrs(i,k+kgases_-1)
#endif
       enddo
     enddo
   endif
#ifdef SWRMDC
!
!  map global aerosol data to model grid
!
   if (iswsrc(1) .gt. 0) then
     do i = 1,LONLENS
       kprf(i) = kprfg(i,lat)
       do k = 1,nxc
         idxc(k,i) = idxcg(k,i,lat)
         cmix(k,i) = cmixg(k,i,lat)
       enddo
       do k = 1,ndn
         denn(k,i) = denng(k,i,lat)
       enddo
     enddo
   endif
!
#endif
   call rad_main_driver(LONLENS,                                               &
               tgrs(1,1),qgrs(1,1),qcis(1,1),pgr(1),                           &
#ifndef SWRMDC
               paerr(1,1,lat),ozon(1,1,lat),albdoa(1),                         &
#else
               ozon(1,1,lat),albdoa(1),                                        &
#endif
#ifndef RMP
               slmskr(1,lat),coszer(1,lat),coszdg(1,lat),                      &
#ifndef NIM
               xlat(1,lat),tsear(1),                                           &
#else
               nimlat(1), tsear(1),                                            &
               si1d(1), sl1d(1), prsi_nim(1,1), prsl_nim(1,1),                 &
#ifdef NIM_DIAG
               fnet,                                                           &
#endif
#endif
#else
               slmsk (1,lat),coszer(1,lat),coszdg(1,lat),                      &
               rlat(1,lat),tsear(1),                                           &
#endif
#if defined (RRTMGSW) || defined (RRTMGLW)
               snoweq (1,lat),vtype(1,lat),                                    &
#endif
               alvbr(1,lat),alnbr(1,lat),alvdr(1,lat),alndr(1,lat),            &
               cldary(1,1),cldtot(1,1),cldcnv(1,1),cldsa(1,1),                 &
               mtopa(1,1),mbota(1,1),                                          &
#ifdef SWRMDC
               rhrh(1,1),kprf(1),idxc(1,1),cmix(1,1),denn(1,1),                &
#endif
#ifndef RMP
               lat,latco,sdec,solc,rsin1,rcos1,rcos2,                          &
#else
               lat,sdec,solc,rsin1,rcos1,rcos2,                                &
#endif
               raddt,dtlw,itimsw,itimlw,kalb,                                  &
               iswsrc,ibnd,ko3,icwp,icfc,                                      &
               swhr(1,1),hlwr(1,1),                                            &
#ifdef EXPLICIT_CLOUDINESS
               qcis_cps(1,1),                                                  &
               taucld(1,1),                                                    &
               cldwp(1,1), cldip(1,1),                                         &
#endif
               sfnswr(1),sfdlwr(1),tsflwr(1),                                  &
#ifdef VIC
               sfdswr(1),sfuswr(1),                                            &
#endif
               gdfvdr(1),gdfndr(1),                                            &
               gdfvbr(1),gdfnbr(1))
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifdef MP
   if (iope) then
#endif
     write(6,*)' after rad_main_driver i',LONF2S,                              &
               tgrs(1,1),qgrs(1,1),qcis(1,1),pgr(1),                           &
#ifndef SWRMDC
               paerr(1,1,lat),ozon(1,1,lat),albdoa(1),                         &
#else
               ozon(1,1,lat),albdoa(1),                                        &
#endif
               slmsk (1,lat),coszer(1,lat),coszdg(1,lat),                      &
#ifdef RMP
               rlat(1,lat),tsear(1),                                           &
#else
#ifdef NIM
               nimlat(1),tsear(1),                                             &
#else
               xlat(1,lat),tsear(1),                                           &
#endif
#endif
               alvbr(1,lat),alnbr(1,lat),alvdr(1,lat),alndr(1,lat),            &
               cldary(1,1),cldtot(1,1),cldcnv(1,1),cldsa(1,1),                 &
               mtopa(1,1),mbota(1,1),                                          &
               lat,sdec,solc,rsin1,rcos1,rcos2,                                &
               raddt,dtlw,itimsw,itimlw,kalb,                                  &
               iswsrc,ibnd,ko3,icwp,icfc,                                      &
               swhr(1,1),hlwr(1,1),                                            &
               sfnswr(1),sfdlwr(1),tsflwr(1),                                  &
               gdfvdr(1),gdfndr(1),                                            &
               gdfvbr(1),gdfnbr(1)
#ifdef MP
   endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif /* DBG end */
#ifdef DG
   do k = 1,levs_
     do i = 1,LONLENS
       if (cldt(i,k).lt.0.e0) then
         cldt(i,k) = 0.e0
       else if (cldt(i,k).gt.1.e0) then
         cldt(i,k) = 100.e0
       else
         cldt(i,k) = cldt(i,k) * 100.e0
       endif
     enddo
   enddo
!
   do k = 1,levs_
     do i = 1,LONLENS
       if (clcv(i,k).lt.0.e0) then
         clcv(i,k) = 0.e0
       else if (clcv(i,k).gt.1.e0) then
         clcv(i,k) = 100.e0
       else
         clcv(i,k) = clcv(i,k) * 100.e0
       endif
     enddo
   enddo
#endif
#ifdef DG3
   call diag_3d_archive(LONLENS,LONF2S,cldt(1,1),dtacc,kdtcld,gda(1,1))
   call diag_3d_archive(LONLENS,LONF2S,clcv(1,1),dtacc,kdtccv,gda(1,1))
#endif
!
   if(itimsw.eq.1) then
     do k = 1,levs_
       do i = 1,LONLENS
#ifdef RMPVECTORIZE
         swh(i,lat,k)=swhr(i,k)
#else
         swh(i,k,lat)=swhr(i,k)
#endif
       enddo
     enddo
     do i = 1,LONLENS
       sfcnsw(i,lat)=sfnswr(i)
#ifdef VIC
       sfcdsw(i,lat)=sfdswr(i)
       sfcusw(i,lat)=sfuswr(i)
#endif
     enddo
   endif
!
   if(itimlw.eq.1) then
     do k = 1,levs_
       do i = 1,LONLENS
#ifdef RMPVECTORIZE
         hlw(i,lat,k)=hlwr(i,k)
#else
         hlw(i,k,lat)=hlwr(i,k)
#endif
       enddo
     enddo
     do i = 1,LONLENS
       sfcdlw(i,lat)=sfdlwr(i)
       tsflw(i,lat)=tsflwr(i)
     enddo
   endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifdef MP
   if(iope) then
#endif
   numlevs=1
   call print_maxmin_seven(swh(1,1,lat),LONLENS,LONF2S,levs_,1,numlevs,'swh')
   write(6,*)'lwh in phys_rad_solver'
   call print_maxmin_seven(hlw(1,1,lat),LONLENS,LONF2S,levs_,1,numlevs,'lwh')
#ifdef MP
   endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif /* DBG end */
!
   do i = 1,LONLENS
     cvavg(i,lat) = cvavg(i,lat) + raddt * cvr(i,lat)
   enddo
!
#ifdef NIM_DIAG
   do k = 1,levs_
     do i = 1,LONLENS
       o3_nim(i,k)  =ozon(i,k,lat)
     enddo
   enddo
!
!  do i=1,lons2
!    write(42,9000) &
!    its,lat,nimlat*180./3.141592,nimlon*180./3.141592,(o3_nim(i,k),k=1,levs_)
!    print*, 'o3_nim',lat,nimlat*180./3.141592,nimlon*180./3.141592,o3_nim(i,25)
!  enddo
!9000 format(' 2007',I6,10x,2(f15.6,1x),F20.15)
!
   do i = 1,LONLENS
     if(itimlw.eq.1) rflux(i,1)=fluxr(i,lat,1)/dtlw   ! up_lwr at the top
     if(itimsw.eq.1) rflux(i,2)=fluxr(i,lat,2)/raddt  ! up_swr at the top
     if(itimsw.eq.1) rflux(i,3)=fluxr(i,lat,3)/raddt  ! up_swr at the sfc
     if(itimsw.eq.1) rflux(i,4)=fluxr(i,lat,4)/raddt  ! down_swr at the sfc
     if(itimsw.eq.1) rflux(i,5)=fluxr(i,lat,18)/raddt ! down_swr at the top
     if(itimsw.eq.1) rflux(i,6)=fnet(i)               ! swr_net at the top
   enddo
#endif
!
#if 0
   call rmp_fill_array(swh,levs_)
   call rmp_fill_array(hlw,levs_)
   call rmp_fill_array(sfcnsw,1)
#ifdef VIC
   call rmp_fill_array(sfcdsw,1)
   call rmp_fill_array(sfcusw,1)
#endif
   call rmp_fill_array(sfcdlw,1)
   call rmp_fill_array(tsflw,1)
   call rmp_fill_array(cvavg,1)
#endif
!
   return
   end subroutine phys_rad_solver
