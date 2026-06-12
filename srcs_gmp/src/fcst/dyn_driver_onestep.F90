#include "define.h"
   subroutine dyn_driver_onestep
!-------------------------------------------------------------------------------
!
! subprogram: dyn_driver_onestep 
!    computes dynamic non-linear tendency terms of temp. div. ln(ps)
!    computes predicted values of vorticity and moisture                           
!                                                                               
! abstract:                                                                     
!   program  starts with spectral  coefficients temp.                           
!   of vorticity, divergence, specific humidity, and                            
!   ln((psfc).  converts them to the gaussian grid at each                      
!   latitude and calls fidi,  for the northern and southern                     
!   hemispheres at the same time.  after return from fidi                       
!   sr.  completes calculation of tendencies of temp. div. and lnps.            
!   specific humidity, and vorticity are predicted by sr. sigvor                
!   all input/output  is via commons.                                           
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
! usage:    call dyn_driver_onestep                                                          
!                                                                               
!-------------------------------------------------------------------------------
! syn(1, 0*levs_+1, lan)  ze                                            
! syn(1, 1*levs_+1, lan)  di                                            
! syn(1, 2*levs_+1, lan)  te                                            
! syn(1, 3*levs_+1, lan)  rq                                            
! syn(1, 4*levs_+1, lan)  uln                                           
! syn(1, 5*levs_+1, lan)  vln                                           
! syn(1, 6*levs_+1, lan)  dpdphi                                        
! syn(1, 6*levs_+2, lan)  dpdlam                                        
! syn(1, 6*levs_+3, lan)  q                                           
! syn(1, 6*levs_+4, lan)  qlap
!
! dyn(1, 0*levs_+1, lan)  d(t)/d(phi)                                   
! dyn(1, 1*levs_+1, lan)  d(rq)/d(phi)                                  
! dyn(1, 2*levs_+1, lan)  d(t)/d(lam)                                   
! dyn(1, 3*levs_+1, lan)  d(rq)/d(lam)                                  
! dyn(1, 4*levs_+1, lan)  d(u)/d(lam)                                   
! dyn(1, 5*levs_+1, lan)  d(v)/d(lam)                                   
! dyn(1, 6*levs_+1, lan)  d(u)/d(phi)                                   
! dyn(1, 7*levs_+1, lan)  d(v)/d(phi)                                   
!
! anl(1, 0*levs_+1, lan)  z     dqdt                                    
! anl(1, 0*levs_+2, lan)  uu    dvdt                                    
! anl(1, 1*levs_+2, lan)  vv    dtdt                                    
! anl(1, 2*levs_+2, lan)  y     dtdt                                    
! anl(1, 3*levs_+2, lan)  rt    drdt                                    
!-------------------------------------------------------------------------------
   use paramodel, only : JCAP1S,LNT2S,LONF2S,LONF22S,LATG2S,LEVSS,LCAPS,       &
                         LCAP22S,jcap_,latg2_,jcap1_,lonf2_,levs_,lnt2_,       &
                         lmx=>levs_,ntotal_ 
   use constant, only  : cp_,g_,hvap_,rd_,rerth_
#ifdef DFS
   use dfsvar, only    : iope
#endif
#ifdef MP
   use paramodel, only : lln22p_,levsp_,lonf22p_,latg2p_,lonf22_,lnt22p_
   use commpi
   use compspec
#endif
#ifdef REDUCE_GRID
   use comreduce
#else
   use paramodel, only : lonf_
#endif
#ifdef DG3
   use diag_3d_module, only   : diag_3d_get, diag_3d_arrange
#endif
#ifdef RAS
   use phys_ras_module, only  : ras_setup
#endif
   use comfibm
   use comgrad
   use comgda
   use comfcst, only   :  spdlat
#ifndef DYNAMIC_ALLOC
   use paramodel, only : ncpus=>ncpus_
#endif
   integer    ::  lots    
   integer    ::  lotst   
   integer    ::  ksz     
   integer    ::  ksd     
   integer    ::  kst     
   integer    ::  ksr     
   integer    ::  ksu     
   integer    ::  kstb    
   integer    ::  ksv     
   integer    ::  kspphi  
   integer    ::  ksplam  
   integer    ::  ksp     
   integer    ::  ksplap  
   integer    ::  lotd    
   integer    ::  kdtphi  
   integer    ::  kdrphi  
   integer    ::  kdtlam  
   integer    ::  kdrlam  
   integer    ::  kdulam  
   integer    ::  kdvlam  
   integer    ::  kduphi  
   integer    ::  kdvphi  
   integer    ::  lota    
   integer    ::  lotat   
   integer    ::  kap     
   integer    ::  kau     
   integer    ::  kav     
   integer    ::  kat     
   integer    ::  kar     
!
   integer    ::  lotss   
   integer    ::  lotsts  
   integer    ::  kszs    
   integer    ::  ksds    
   integer    ::  ksts    
   integer    ::  ksrs    
   integer    ::  ksus    
   integer    ::  kstbs   
   integer    ::  ksvs    
   integer    ::  kspphis 
   integer    ::  ksplams 
   integer    ::  ksps    
   integer    ::  ksplaps 
   integer    ::  lotds   
   integer    ::  kdtphis 
   integer    ::  kdrphis 
   integer    ::  kdtlams 
   integer    ::  kdrlams 
   integer    ::  kdulams 
   integer    ::  kdvlams 
   integer    ::  kduphis 
   integer    ::  kdvphis 
   integer    ::  lotas   
   integer    ::  lotats  
   integer    ::  kaps    
   integer    ::  kaus    
   integer    ::  kavs    
   integer    ::  kats    
   integer    ::  kars    
!
   integer    ::  lotf    
   integer    ::  kfu     
   integer    ::  kfv     
   integer    ::  kft     
   integer    ::  kfr     
!
   real, allocatable   ::  fnl(:,:,:)
#ifdef MP
   real, allocatable    ::  syf(:,:,:)
   real, allocatable    ::  grs(:,:,:)
   real, allocatable    ::  dyf(:,:,:)
   real, allocatable    ::  dgr(:,:,:)
   real, allocatable    ::  anf(:,:,:)
   real, allocatable    ::  gra(:,:,:)
#endif
!
   real, allocatable    ::  scs2(:)
   real, allocatable    ::  syn(:,:,:)
   real, allocatable    ::  syntop(:,:,:)
   real, allocatable    ::  dyn(:,:,:)
   real, allocatable    ::  anl(:,:,:)
   real, allocatable    ::  anltop(:,:,:)
   real, allocatable    ::  flp(:,:,:,:)
   real, allocatable    ::  flm(:,:,:,:)
!
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
#ifdef DG
   real                 ::  cldt(LONF2S,levs_,NCPUSS)
   real                 ::  clcv(LONF2S,levs_,NCPUSS)
   real                 ::  tgmxl(NCPUSS)
   real                 ::  tgmnl(NCPUSS)
   integer              ::  igmxl(NCPUSS),kgmxl(NCPUSS)
   integer              ::  igmnl(NCPUSS),kgmnl(NCPUSS)
#endif
#ifdef DG3
   real                 ::  gda(nwgda,kdgda,NCPUSS)
#endif
#undef NCPUSS
!
   logical,save         ::  fnl1st
   data                     fnl1st/.true./
   logical              ::  dorad
#ifdef RAS
   logical,parameter    ::  ras=.true.
   integer,parameter    ::  nsphys=1
   real   ,parameter    ::  cp=cp_, alhl=hvap_, grav=g_, rgas=rd_
   real                 ::  sig(lmx+1),prj(lmx+1), prh(lmx)
   real                 ::  fpk(lmx),  hpk(lmx)
   real                 ::  sgb(lmx),  ods(lmx),rasal(lmx), prns(lmx/2)
   real                 ::  rannum(200)
!
   call ras_setup(lmx, si, sl, del, cp, rgas, deltim                           &
               ,nsphys, thour                                                  &
               ,sig, sgb, prh, prj, hpk, fpk, ods, prns                        &
               ,rasal, lm, krmin, krmax, nstrp                                 &
               ,ncrnd, rannum, afac, ufac)
#endif
!
#ifdef MP
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   lons2=lonlen(mype)*2
   lats2=latlen(mype)
#else
   llstr=0
   llens=jcap1_
   lons2=lonf2_
   lats2=latg2_
#endif
!
   rdt=0.5/deltim
   lots    =6*levs_+4
   lotst   =2*levs_+1
   ksz     =1
   ksd     =1*levs_+1
   kst     =2*levs_+1
   ksr     =3*levs_+1
   ksu     =4*levs_+1
   kstb    =4*levs_+1
   ksv     =5*levs_+1
   kspphi  =6*levs_+1
   ksplam  =6*levs_+2
   ksp     =6*levs_+3
   ksplap  =6*levs_+4
   lotd    =8*levs_
   kdtphi  =0*levs_+1
   kdrphi  =1*levs_+1
   kdtlam  =2*levs_+1
   kdrlam  =3*levs_+1
   kdulam  =4*levs_+1
   kdvlam  =5*levs_+1
   kduphi  =6*levs_+1
   kdvphi  =7*levs_+1
   lota    =4*levs_+1
   lotat   =2*levs_
   kap     =1
   kau     =2
   kav     =1*levs_+2
   kat     =2*levs_+2
   kar     =3*levs_+2
!
   lotss   =6*LEVSS+4
   lotsts  =2*LEVSS+1
   kszs    =1
   ksds    =1*LEVSS+1
   ksts    =2*LEVSS+1
   ksrs    =3*LEVSS+1
   ksus    =4*LEVSS+1
   kstbs   =4*LEVSS+1
   ksvs    =5*LEVSS+1
   kspphis =6*LEVSS+1
   ksplams =6*LEVSS+2
   ksps    =6*LEVSS+3
   ksplaps =6*LEVSS+4
   lotds   =8*LEVSS
   kdtphis =0*LEVSS+1
   kdrphis =1*LEVSS+1
   kdtlams =2*LEVSS+1
   kdrlams =3*LEVSS+1
   kdulams =4*LEVSS+1
   kdvlams =5*LEVSS+1
   kduphis =6*LEVSS+1
   kdvphis =7*LEVSS+1
   lotas   =4*LEVSS+1
   lotats  =2*LEVSS
   kaps    =1
   kaus    =2
   kavs    =1*LEVSS+2
   kats    =2*LEVSS+2
   kars    =3*LEVSS+2
!
   lotf    =4*levs_
   kfu     =1
   kfv     =kfu+levs_
   kft     =kfv+levs_
   kfr     =kft+levs_
!
   allocate(fnl(LONF22S,lotf,LATG2S))
#ifdef MP
   allocate(syf(lonf22_,lotss,latg2p_))
   allocate(grs(lonf22p_,lots,latg2p_))
   allocate(dyf(lonf22_,lotds,latg2p_))
   allocate(dgr(lonf22p_,lotd,latg2p_))
   allocate(anf(lonf22_,lotas,latg2p_))
   allocate(gra(lonf22p_,lota,latg2p_))
#endif
!
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus
#endif
   allocate(scs2(JCAP1S))
   allocate(syn(LCAP22S,lotss,NCPUSS))
   allocate(syntop(2,JCAP1S,lotsts))
   allocate(dyn(LCAP22S,lotds,NCPUSS))
   allocate(anl(LCAP22S,lotas,NCPUSS))
   allocate(anltop(2,JCAP1S,lotats))
   allocate(flp(2,JCAP1S,lotas,NCPUSS))
   allocate(flm(2,JCAP1S,lotas,NCPUSS))
#undef NCPUSS
!
   call rad_prepare(idate,dorad)
!
   if(dorad.and.iope)                                                          &
       write(6,*) ' run radiation at solar hour=',solhr
!
   if(fnl1st) then
      if(iope)                                                                 &
         write(6,*) ' run t-1 initial physics forcing '
!
#ifdef MP
     call mpnn2nk(zem ,lnt22p_,levs_, zea,lln22p_,levsp_,1)
     call mpnn2nk(dim ,lnt22p_,levs_, dia,lln22p_,levsp_,1)
     call mpnn2nk(tem ,lnt22p_,levs_, tea,lln22p_,levsp_,1)
     call mpnn2nk(rm  ,lnt22p_,levs_, rqa,lln22p_,levsp_,ntotal_)
     call mpnn2n (qm,lnt22p_, qa,lln22p_,1)
#define QMS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define DIMS dia
#define ZEMS zea
#define ULNS ulna
#define VLNS vlna
#else
#define QMS qm
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define DIMS dim
#define ZEMS zem
#define ULNS uln
#define VLNS vln
#endif
     call sph_del_log_ps(QMS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,            &
               llstr,llens,lwvdef)
     call sph_solve_laplacian(QMS,QLAPS,llstr,llens,lwvdef)
     call sph_divor2wind(DIMS,ZEMS,ULNS,VLNS,syntop(1,1,1),                    &
               syntop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef QMS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIMS
#undef ZEMS
#undef ULNS
#undef VLNS
!
#ifndef MP
     last=mod(latg2_,ncpus)
     nggs=(latg2_-last)/ncpus
     if(last.ne.0)nggs=nggs+1
     inclat=ncpus
     lat1=1-ncpus
     lat2=0
     latdon=0
     do ngg = 1,nggs
       if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
       lat1=lat1+ncpus
       lat2=lat2+inclat
#endif
!
#ifdef MP
       lat1=1
       lat2=latg2_
       latdon=0
#define QMS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define DIMS dia
#define TEMS tea
#define RMS  rqa
#define ULNS ulna
#define VLNS vlna
#else
#define QMS qm
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define DIMS dim
#define TEMS tem
#define RMS  rm 
#define ULNS uln
#define VLNS vln
#endif
#ifdef ORIGIN_THREAD
!$doacross share(qtt,qvv,colrad,dim,dia,syn,uln,ulna,qm,qa,qlap,qlapa,
!$&              syntop,lat1,lat2,latdon,llstr,llens,lwvdef,
!$&              lcapd,lcapdp),
!$&        local(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
       do lat = lat1,lat2
         lan=lat-latdon
#ifdef MP
#ifdef REDUCE_GRID
         llensd=lcapdp(lat,mype)
#else
         llensd=llens
#endif
         call sph_sum_coeff(dia    ,syn(1,ksds,lan),                           &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(tea    ,syn(1,ksts,lan),                           &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(rqa    ,syn(1,ksrs,lan),                           &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEHSS)
         call sph_sum_coeff(ulna   ,syn(1,ksus,lan),                           &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(vlna   ,syn(1,ksvs,lan),                           &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(dpdphia,syn(1,kspphis,lan),                        &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
         call sph_sum_coeff(dpdlama,syn(1,ksplams,lan),                        &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
         call sph_sum_coeff(qa     ,syn(1,ksplaps,lan),                        &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
         call sph_sum_coeff(qlapa  ,syn(1,ksps,lan)   ,                        &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
#else
#ifdef REDUCE_GRID
         llensd=lcapd(lat)
#else
         llensd=llens
#endif
         call sph_sum_coeff(dim   ,syn(1,ksds,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(tem   ,syn(1,ksts,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(rm    ,syn(1,ksrs,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(uln   ,syn(1,ksus,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(vln   ,syn(1,ksvs,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff(dpdphi,syn(1,kspphis,lan),                         &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
         call sph_sum_coeff(dpdlam,syn(1,ksplams,lan),                         &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
         call sph_sum_coeff(qm    ,syn(1,ksplaps,lan),                         &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
         call sph_sum_coeff(qlap  ,syn(1,ksps,lan)   ,                         &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
#endif
         call sph_sum_coeff_top(syn(1,ksus,lan)   ,syntop(1,1,1)        ,      &
                                          qvv(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff_top(syn(1,ksvs,lan)   ,syntop(1,1,LEVSS+1)  ,      &
                                          qvv(1,lat),llstr,llensd,lwvdef,LEVSS)
         call sph_sum_coeff_top(syn(1,kspphis,lan),syntop(1,1,2*LEVSS+1),      &
                                          qvv(1,lat),llstr,llensd,lwvdef,1)
       enddo
#undef QMS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIMS
#undef TEMS
#undef RMS
#undef ULNS
#undef VLNS
!
#ifdef MP
       call mpnl2ny(syn,lcap22p_,latg2_,                                       &
                    syf,lonf22_,latg2p_,lotss,ksds,5*levsp_+4)
       lat1=jstr
       lat2=jend
       latdon=jstr-1
#define SYNS syf
#else
#define SYNS syn
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(syn,syf,lat1,lat2,latdon,lcapd,lonfd,lstdef),
!$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syn,syf,lat1,lat2,latdon,lcapd,lonfd,latdef)
!mic$& private(lat,lan,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
       do lat = lat1,lat2
         lan=lat-latdon
#ifdef REDUCE_GRID
         lcapf=lcapd(latdef(lat))
         lonff=lonfd(latdef(lat))
#else
         lcapf=jcap1_
         lonff=lonf_
#endif
         call sph_wave2grid(SYNS(1,ksds,lan),SYNS(1,ksds,lan),                 &
                      10*LEVSS+8,lcapf,lonff,latdef(lat),1)
       enddo
#undef SYNS
!
#ifdef MP
       call mpnk2nx(syf,lonf22_,lotss,                                         &
                 grs,lonf22p_,lots,latg2p_,levsp_,levs_,ksds,ksd,5)
       call mpx2nx (syf,lonf22_,lotss,                                         &
                 grs,lonf22p_,lots,latg2p_,kspphis,kspphi,4)
       lat1=jstr
       lat2=jend
       latdon=jstr-1
#define SYNS grs
#define LATX lan
#else
#define SYNS syn
#define LATX lat
#endif
!
       if(dorad) then
!
#ifdef ORIGIN_THREAD
!$doacross share(syn,grs,lat1,lat2,latdon,gda,cldt,clcv,sfcp,
!$&              lons2,lats2,lonfdp,mype)
!$&              ,local(lat,lan,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syn,grs,lat1,lat2,latdon,gda,cldt,clcv,sfcp)
!mic$& shared(lons2,lats2,lonfdp,mype)
!mic$& private(lat,lan,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lonsd2)
#endif
         do lat = lat1,lat2
           latrue=latdef(lat)
           lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
           lonsd2=lonfdp(lan,mype)*2
#else
           lonsd2=lonfd(latdef(lat))*2
#endif
#else
           lonsd2=lons2
#endif
#ifdef MP
           if( lonsd2.gt.0 ) then
#endif
             call phys_rad_solver(lonsd2,lats2,                                &
                     SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),                    &
                     SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),          &
                     SYNS(1,kst,lan),SYNS(1,ksr,lan),SYNS(1,ksd,lan),          &
#ifdef DG3
                     gda(1,1,lan),                                             &
#endif
#ifdef DG
                     cldt(1,1,lan),clcv(1,1,lan),                              &
#endif
                     LATX,latrue)
#ifdef MP
           endif
#endif
         enddo
       endif
!
#ifdef ORIGIN_THREAD
!$doacross share(syn,grs,fnl,lat1,lat2,latdon,rdt,gda,lons2,
!$&              tmgxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,ras,
!$&              lmx,cp,alhl,grav,rgas,sig,sgb,prh,prj,hpk,fpk,ods,
!$&              prns,rasal,lm,krmin,krmax,nstrp,ncrnd,rannum,afac,
!$&              ufac,lonfdp,mype),local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syn,grs,fnl,lat1,lat2,latdon,rdt,gda,lons2)
!mic$& shared(tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl)
!mic$& shared(ras,lmx,cp,alhl,grav,rgas)
!mic$& shared(sig, sgb, prh, prj, hpk, fpk, ods, prns)
!mic$& shared(rasal, lm, krmin, krmax, nstrp)
!mic$& shared(ncrnd, rannum, afac, ufac,lonfdp,mype)
!mic$& private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
       do lat = lat1,lat2
         lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
         lonsd2=lonfdp(lan,mype)*2
#else
         lonsd2=lonfd(latdef(lat))*2
#endif
#else
         lonsd2=lons2
#endif
!
#ifdef MP
         if( lonsd2.gt.0 ) then
#endif
           call phys_main_solver(lonsd2,                                       &
                  SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),                       &
                  SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),             &
                  SYNS(1,kst,lan),SYNS(1,ksr,lan),SYNS(1,ksd,lan),             &
                  SYNS(1,ksplap,lan),                                          &
                  fnl(1,kft,LATX),fnl(1,kfr,LATX),                             &
                  fnl(1,kfu,LATX),fnl(1,kfv,LATX),                             &
#ifdef DG
                  tgmxl(lan),igmxl(lan),kgmxl(lan),                            &
                  tgmnl(lan),igmnl(lan),kgmnl(lan),                            &
#endif
#ifdef DG3
                  gda(1,1,lan),                                                &
#endif
#ifdef RAS
                  ras,lmx,cp,alhl,grav,rgas,                                   &
                  sig, sgb, prh, prj, hpk, fpk, ods, prns,                     &
                  rasal, lm, krmin, krmax, nstrp,                              &
                  ncrnd, rannum, afac, ufac,                                   &
#endif
                  LATX,0.)
#ifdef MP
         endif
#endif
!
         do k = 1,levs_
           do j = 1,lonsd2
             fnl(j,kfu-1+k,LATX)=                                              &
                     (fnl(j,kfu-1+k,LATX)-SYNS(j,ksu-1+k,lan))*rdt
             fnl(j,kfv-1+k,LATX)=                                              &
                     (fnl(j,kfv-1+k,LATX)-SYNS(j,ksv-1+k,lan))*rdt
             fnl(j,kft-1+k,LATX)=                                              &
                     (fnl(j,kft-1+k,LATX)-SYNS(j,kst-1+k,lan))*rdt
             fnl(j,kfr-1+k,LATX)=                                              &
                     (fnl(j,kfr-1+k,LATX)-SYNS(j,ksr-1+k,lan))*rdt
           enddo
         enddo
       enddo
#undef LATX
#undef SYNS
!
#ifndef MP
         latdon=latdon+(lat2-lat1+1)
     enddo      ! end of do ngg
#endif
      fnl1st=.false.
      dorad=.false.
   endif
!
#ifdef MP
   call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,1)
   call mpnn2nk(di ,lnt22p_,levs_, dia,lln22p_,levsp_,1)
   call mpnn2nk(te ,lnt22p_,levs_, tea,lln22p_,levsp_,1)
   call mpnn2nk(rq ,lnt22p_,levs_, rqa,lln22p_,levsp_,ntotal_)
   call mpnn2n(q ,lnt22p_, qa,lln22p_,1)
#define QS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define DIS dia
#define ZES zea
#define ULNS ulna
#define VLNS vlna
#else
#define QS q
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define DIS di
#define ZES ze
#define ULNS uln
#define VLNS vln
#endif
   call sph_del_log_ps(QS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,               &
               llstr,llens,lwvdef)
   call sph_solve_laplacian(QS,QLAPS,llstr,llens,lwvdef)
   call sph_divor2wind(DIS,ZES,ULNS,VLNS,syntop(1,1,1),                        &
               syntop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
!
   do n = 1,2                                                                  
     do j = 1,JCAP1S
       do l = 1,lotats
         anltop(n,j,l)=0.0                                                   
       enddo                                                                 
     enddo                                                                   
   enddo                                                                     
#ifdef DG
   tgmx=-1.e20
   tgmn= 1.e20
#endif
!
#ifdef MP
#define LNT2X lln22p_
#define UUS uua
#define VVS vva
#define YS ya
#define RTS rta
#define ZS za
#else
#define LNT2X lnt2_
#define UUS uu
#define VVS vv
#define YS y
#define RTS rt
#define ZS z
#endif
!
#ifdef ORIGIN_THREAD
!
!$doacross share(uu,uua,vv,vva,y,ya,rt,rta),local(j,k)  
!
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(uu,uua,vv,vva,y,ya,rt,rta)        
!mic$1 private(j,k)                                                             
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
   do k = 1,LEVSS 
     do j = 1,LNT2X
       UUS(j,k)=0.0
       VVS(j,k)=0.0
       YS(j,k)=0.0
       RTS(j,k)=0.0
     enddo
   enddo
   do j = 1,LNT2X 
     ZS(j,1)=0.e0
   enddo
#undef LNT2X
#undef UUS
#undef VVS
#undef YS
#undef RTS
#undef ZS
!
#ifndef MP
! compute latitude band limits                                                  
!
   last=mod(latg2_,ncpus)
   nggs=(latg2_-last)/ncpus
   if(last.ne.0)nggs=nggs+1
   inclat=ncpus
   lat1=1-ncpus
   lat2=0
   latdon=0
!                                                                              
   do ngg = 1,nggs
     if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
     lat1=lat1+ncpus
     lat2=lat2+inclat
#endif
!
#ifdef DG3
#ifdef MP
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define LATX lan
#else
#define LATX lat
#endif
     do lat = lat1,lat2
       lan=lat-latdon
       call diag_3d_get(LATX,nwgda*kdgda,gda(1,1,lan))
     enddo
#undef LATX
#endif
#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
     do k = 1,NCPUSS
       tgmxl(k)=tgmx
       tgmnl(k)=tgmn
     enddo
#undef NCPUSS
#endif
!
#ifdef MP
     lat1=1
     lat2=latg2_
     latdon=0
#define ZES zea
#define DIS dia
#define TES tea
#define RQS rqa
#define ULNS ulna
#define VLNS vlna
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QS qa
#define QLAPS qlapa
#else
#define ZES ze
#define DIS di
#define TES te
#define RQS rq
#define ULNS uln
#define VLNS vln
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QS q
#define QLAPS qlap
#endif
!
! first lat loop                                                                
!
#ifdef ORIGIN_THREAD
!$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
!$&        colrad,ze,zea,llstr,llens,lwvdef,lcapdp,lcapd,mype)
!$&        local(lat,lan,llensd)                                                       
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(syntop,syn,qtt,qvv,lat1,lat2,latdon)
!mic$1 shared(colrad,ze,zea,llstr,llens,lwvdef,lcapdp,lcapd,mype)
!mic$1 private(lat,lan,llensd)                                                         
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
!                                                                               
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
       call sph_sum_coeff(ZES    ,syn(1,kszs,lan),                             &
                                     qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(DIS    ,syn(1,ksds,lan),                             &
                                     qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(TES    ,syn(1,ksts,lan),                             &
                                     qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(RQS    ,syn(1,ksrs,lan),                             &
                                     qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(ULNS   ,syn(1,ksus,lan),                             &
                                     qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(VLNS   ,syn(1,ksvs,lan),                             &
                                     qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(DPDPHIS,syn(1,kspphis,lan),                          &
                                     qtt(1,lat),llstr,llensd,lwvdef,1)
       call sph_sum_coeff(DPDLAMS,syn(1,ksplams,lan),                          &
                                     qtt(1,lat),llstr,llensd,lwvdef,1)
       call sph_sum_coeff(QS     ,syn(1,ksps   ,lan),                          &
                                     qtt(1,lat),llstr,llensd,lwvdef,1)
       call sph_sum_coeff(QLAPS  ,syn(1,ksplaps,lan),                          &
                                     qtt(1,lat),llstr,llensd,lwvdef,1)
       call sph_sum_coeff_top(syn(1,ksus,lan)   ,syntop(1,1,1),                &
                                     qvv(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff_top(syn(1,ksvs,lan)   ,syntop(1,1,LEVSS+1),          &
                                     qvv(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff_top(syn(1,kspphis,lan),syntop(1,1,2*LEVSS+1),        &
                                     qvv(1,lat),llstr,llensd,lwvdef,1)
     enddo
!                                                                               
! compute merid. derivs. of temp. and moisture using qdd.                       
!                                                                               
#ifdef ORIGIN_THREAD
!$doacross share(dyn,qdd,lat1,lat2,latdon,                
!$&        llens,te,tea,llstr,llens,lwvdef,lcapd,lcapdp,mype),
!$&        local(lat,lan,i,k,llensd)
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(dyn,qdd)                                                      
!mic$1 shared(lat1,lat2,latdon,llens)
!mic$1 shared(te,tea,llstr,llens,lwvdef,lcapd,lcapdp,mype)
!mic$1 private(lat,lan,i,k,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k,llensd)
#endif
!                                                                               
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
       call sph_sum_coeff(TES,dyn(1,kdtphis,lan),                              &
                                     qdd(1,lat),llstr,llensd,lwvdef,LEVSS)
       call sph_sum_coeff(RQS,dyn(1,kdrphis,lan),                              &
                                     qdd(1,lat),llstr,llensd,lwvdef,LEVSS)
!
!...... d(t)/d(phi)  d(rq)/d(phi) in s. hemi.
!
       do k = 1,LEVSS
         do i = 1,llensd*2
           dyn(i+LCAPS,kdtphis-1+k,lan)=-dyn(i+LCAPS,kdtphis-1+k,lan)
           dyn(i+LCAPS,kdrphis-1+k,lan)=-dyn(i+LCAPS,kdrphis-1+k,lan)
         enddo
       enddo
     enddo
!
#ifdef ORIGIN_THREAD
!$doacross share(dyn,rcs2,syn,lat1,lat2,latdon,llens,
!$&              llstr,lwvdef,lcapd,lcapdp,mype),
!$&              local(lat,lan,j,k,l,kdd,kss,scs2,llensd) 
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(dyn,rcs2,llstr,llens,lwvdef)
!mic$1 shared(syn,lat1,lat2,latdon,lcapd,lcapdp,mype)
!mic$1 private(lat,lan,j,k,l,kdd,kss,scs2,llensd)         
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
!                                                                               
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
!                                                                               
!   calculate t rq u v zonal derivs. by multiplication with i*l                 
!                                                                               
       do l = 1,llensd
         j=lwvdef(llstr+l)
         scs2(l)=float(j)/rerth_ *rcs2(lat)
       enddo
!
!   calculate t rq u v zonal derivs. by multiplication with i*l
!   note scs2=rcs2*l/rerth_
!
! ....... d(t)/d(lam)  d(rq)/d(lam)  d(u)/d(lam)  d(v)/d(lam) .......
!
       do k = 1,4*LEVSS
         kdd=kdtlams-1+k
         kss=ksts-1+k
         do j = 1,llensd
           dyn(       2*j-1,kdd,lan)=-syn(       2*j  ,kss,lan)*scs2(j)
           dyn(       2*j  ,kdd,lan)= syn(       2*j-1,kss,lan)*scs2(j)
           dyn(LCAPS+2*j-1,kdd,lan)=-syn(LCAPS+2*j  ,kss,lan)*scs2(j)
           dyn(LCAPS+2*j  ,kdd,lan)= syn(LCAPS+2*j-1,kss,lan)*scs2(j)
         enddo
       enddo
     enddo
!
#undef ZES
#undef DIS
#undef TES
#undef RQS
#undef ULNS
#undef VLNS
#undef DPDPHIS
#undef DPDLAMS
#undef QS
#undef QLAPS
!
#ifdef MP
     call mpnl2ny(syn,lcap22p_,latg2_,                                         &
                syf,lonf22_,latg2p_,lotss,1,lotss)
     call mpnl2ny(dyn,lcap22p_,latg2_,                                         &
                dyf,lonf22_,latg2p_,lotds,1,6*levsp_)
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define SYNS syf
#define DYNS dyf
#else
#define SYNS syn
#define DYNS dyn
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(syn,syf,dyn,dyf,lat1,lat2,latdon,lcapd,lonfd,latdef),
!$&        local(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syn,syf,dyn,dyf)
!mic$& shared(lat1,lat2,latdon,lcapd,lonfd,latdef)
!mic$& private(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,k,j,lcapf,lonff)
#endif
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,                 &
                   lcapf,lonff,latdef(lat),1)
!
       call sph_wave2grid(DYNS(1,1,lan),DYNS(1,1,lan),12*LEVSS,                &
                   lcapf,lonff,latdef(lat),1)
       do k = 1,LEVSS
         do j = 1,lonff*2
           DYNS(j,kduphis-1+k,lan)=                                            &
           DYNS(j,kdvlams-1+k,lan)-SYNS(j,kszs-1+k,lan)
           DYNS(j,kdvphis-1+k,lan)=                                            &
           -DYNS(j,kdulams-1+k,lan)+SYNS(j,ksds-1+k,lan)
         enddo
       enddo
     enddo
!
#undef SYNS
#undef DYNS
!
#ifdef MP
     call mpnk2nx(syf,lonf22_,lotss,                                           &
                 grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,6)
     call mpx2nx (syf,lonf22_,lotss,                                           &
                 grs,lonf22p_,lots,latg2p_,kspphis,kspphi,4)
     call mpnk2nx(dyf,lonf22_,lotds,                                           &
                 dgr,lonf22p_,lotd,latg2p_,levsp_,levs_,1,1,8)
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define SYNS grs
#define DYNS dgr
#define ANLS gra
#define LATX lan
#else
#define SYNS syn
#define DYNS dyn
#define ANLS anl
#define LATX lat
#endif
!
     if( dorad ) then
!
#ifdef ORIGIN_THREAD
!$doacross share(syn,grs,lat1,lat2,latdon,gda,
!$&              cldt,clcv,sfcp,lons2,lats2,lonfdp,mype),
!$&        local(lat,lan,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syn,grs,lat1,lat2,latdon,gda)
!mic$& shared(cldt,clcv,sfcp,lons2,lats2,lonfdp,mype)
!mic$& private(lat,lan,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lonsd2)
#endif
!
     do lat = lat1,lat2
       latrue=latdef(lat)
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       lonsd2=lonfdp(lan,mype)*2
#else
       lonsd2=lonfd(latdef(lat))*2
#endif
#else
       lonsd2=lons2
#endif
#ifdef MP
       if( lonsd2.gt.0 ) then
#endif
         call phys_rad_solver(lonsd2,lats2,                                    &
               SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),                          &
               SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),                &
               SYNS(1,kst,lan),SYNS(1,ksr,lan),SYNS(1,ksd,lan),                &
#ifdef DG3
               gda(1,1,lan),                                                   &
#endif
#ifdef DG
               cldt(1,1,lan),clcv(1,1,lan),                                    &
#endif
               LATX,latrue)
#ifdef MP
         endif
#endif
       enddo
     endif
!
#ifdef ORIGIN_THREAD
!$doacross share(syn,grs,anl,gra,fnl,lat1,lat2,latdon,gda,rdt,lons2,
!$&              tmgx,igmx,kgmx,tgmn,igmn,kgmn,
!$&              tmgxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,ras,
!$&              lmx,cp,alhl,grav,rgas,sig,sgb,prh,prj,hpk,fpk,ods,
!$&              prns,rasal,lm,krmin,krmax,nstrp,ncrnd,rannum,afac,
!$&              ufac,lonfdp,mype),local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syn,grs,anl,gra,fnl,lat1,lat2,latdon,gda,rdt,lons2)
!mic$& shared(tgmx ,igmx ,kgmx ,tgmn ,igmn ,kgmn )
!mic$& shared(tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl)
!mic$& shared(ras,lmx,cp,alhl,grav,rgas)
!mic$& shared(sig, sgb, prh, prj, hpk, fpk, ods, prns)
!mic$& shared(rasal, lm, krmin, krmax, nstrp)
!mic$& shared(ncrnd, rannum, afac, ufac,lonfdp,mype)
!mic$& private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       lonsd2=lonfdp(lan,mype)*2
#else
       lonsd2=lonfd(latdef(lat))*2
#endif
#else
       lonsd2=lons2
#endif
!
! ------------------------- subgrid forcing -----------------------
!
#ifdef MP
       if( lonsd2.gt.0 ) then
#endif
         do k = 1,levs_
           do j = 1,lonsd2
             ANLS(j,kau-1+k,lan)=fnl(j,kfu-1+k,LATX)
             ANLS(j,kav-1+k,lan)=fnl(j,kfv-1+k,LATX)
             ANLS(j,kat-1+k,lan)=fnl(j,kft-1+k,LATX)
             ANLS(j,kar-1+k,lan)=fnl(j,kfr-1+k,LATX)
           enddo
         enddo
!
         call phys_main_solver(lonsd2,                                         &
               SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),                          &
               SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),                &
               SYNS(1,kst,lan),SYNS(1,ksr,lan),SYNS(1,ksd,lan),                &
               SYNS(1,ksplap,lan),                                             &
               fnl(1,kft,LATX),fnl(1,kfr,LATX),                                &
               fnl(1,kfu,LATX),fnl(1,kfv,LATX),                                &
#ifdef DG
               tgmxl(lan),igmxl(lan),kgmxl(lan),                               &
               tgmnl(lan),igmnl(lan),kgmnl(lan),                               &
#endif
#ifdef DG3
               gda(1,1,lan),                                                   &
#endif
#ifdef RAS
               ras,lmx,cp,alhl,grav,rgas,                                      &
               sig, sgb, prh, prj, hpk, fpk, ods, prns,                        &
               rasal, lm, krmin, krmax, nstrp,                                 &
               ncrnd, rannum, afac, ufac,                                      &
#endif
               LATX,1.)
!
         do k = 1,levs_
           do j = 1,lonsd2
             fnl(j,kfu-1+k,LATX)=                                              &
                     (fnl(j,kfu-1+k,LATX)-SYNS(j,ksu-1+k,lan))*rdt
             fnl(j,kfv-1+k,LATX)=                                              &
                     (fnl(j,kfv-1+k,LATX)-SYNS(j,ksv-1+k,lan))*rdt
             fnl(j,kft-1+k,LATX)=                                              &
                     (fnl(j,kft-1+k,LATX)-SYNS(j,kst-1+k,lan))*rdt
             fnl(j,kfr-1+k,LATX)=                                              &
                     (fnl(j,kfr-1+k,LATX)-SYNS(j,ksr-1+k,lan))*rdt
           enddo
         enddo
#ifdef DG
         if(tgmxl(lan).gt.tgmx) then
           tgmx=tgmxl(lan)
           igmx=igmxl(lan)
           kgmx=kgmxl(lan)
           jgmx=lat
         else if(tgmnl(lan).lt.tgmn) then
           tgmn=tgmnl(lan)
           igmn=igmnl(lan)
           kgmn=kgmnl(lan)
           jgmn=lat
         endif
#endif
#ifdef MP
       endif
#endif
#undef LATX
!
     enddo
!
#ifdef DG3
#ifdef MP
#define LATX lan
#else
#define LATX lat
#endif
     do lat = lat1,lat2
       lan=lat-latdon
       call diag_3d_arrange(LATX,nwgda*kdgda,gda(1,1,lan))
     enddo
#undef LATX
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(dyn,dgr,rcs2,syn,grs,anl,gra,lat1,lat2,latdon,spdlat,
!$&              del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype),
!$&              local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(dyn,dgr,rcs2)
!mic$1 shared(syn,grs,anl,gra,lat1,lat2,latdon,spdlat)
!mic$1 shared(del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype)
!mic$1 private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
!                                                                               
     do lat = lat1,lat2
       lan=lat-latdon
!                                                                               
#ifdef REDUCE_GRID
#ifdef MP
#define LATX lan
       lonsd2=lonfdp(lan,mype)*2
#else
#define LATX lat
       lonsd2=lonfd(latdef(lat))*2
#endif
#else
#define LATX lat
       lonsd2=lons2
#endif
#ifdef MP
       if( lonsd2.gt.0 ) then
#endif
         do k = 1,levs_ 
           do j = 1,lonsd2
             SYNS(j,kst-1+k,lan)=SYNS(j,kst-1+k,lan)-tov(k)
           enddo      
         enddo      
!                                                                               
#ifndef HYBRID
         call sph_nonlinear_tend(lonsd2,                                       &
               SYNS(1,ksd,lan),SYNS(1,kst,lan),                                &
               SYNS(1,ksz,lan),SYNS(1,ksu,lan),                                &
               SYNS(1,ksv,lan),SYNS(1,ksr,lan),                                &
               SYNS(1,kspphi,lan),SYNS(1,ksplam,lan),                          &
               rbs2(LATX),del,rdel2,ci,p1,p2,h1,h2,tov,spdlat(1,LATX),         &
               DYNS(1,kdtphi,lan),DYNS(1,kdtlam,lan),                          &
               DYNS(1,kdrphi,lan),DYNS(1,kdrlam,lan),                          &
               DYNS(1,kdulam,lan),DYNS(1,kdvlam,lan),                          &
               DYNS(1,kduphi,lan),DYNS(1,kdvphi,lan),                          &
               ANLS(1,kap,lan),ANLS(1,kat,lan),                                &
               ANLS(1,kar,lan),ANLS(1,kau,lan),                                &
               ANLS(1,kav,lan))
#endif
#ifdef MP
       endif
#endif
#undef LATX
!
     enddo
#undef SYNS
#undef DYNS
#undef ANLS
!
#ifdef MP
     call mpnx2nk(gra,lonf22p_,lota,                                           &
                anf,lonf22_ ,lotas,latg2p_,levs_,levsp_,kau,kaus,4)
     call mpnx2x (gra,lonf22p_,lota,                                           &
                anf,lonf22_ ,lotas,latg2p_,kap,kaps,1)
     if( levlen(mype).lt.levsp_ ) then
       do lan = 1,lats2
         do i = 1,lonf22_
               anf(i,kaus+levsp_-1,lan)=0.0
               anf(i,kavs+levsp_-1,lan)=0.0
               anf(i,kats+levsp_-1,lan)=0.0
               anf(i,kars+levsp_-1,lan)=0.0
         enddo
       enddo
     endif
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define ANLS anf
#else
#define ANLS anl
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(latdon,lat1,lat2,anl,anf,lcapd,lonfd,latdef),                                 
!$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(latdon,lat1,lat2,anl,anf,lcapd,lonfd,latdef) 
!mic$1 private(lat,lan,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(ANLS(1,1,lan),ANLS(1,1,lan),2*lotas,                 &
                  lcapf,lonff,latdef(lat),-1)
     enddo
#undef ANLS
!
#ifdef MP
     call mpny2nl(anf,lonf22_ ,latg2p_,                                        &
                anl,lcap22p_,latg2_ ,lotas,kaps,lotas)
!
     lat1=1
     lat2=latg2_
     latdon=0
#define ZS za
#define UUS uua
#define VVS vva
#define YS ya
#define RTS rta
#else
#define ZS z
#define UUS uu
#define VVS vv
#define YS y
#define RTS rt
#endif
!                                                                               
#ifdef ORIGIN_THREAD
!$doacross share(latdon,lat1,lat2,anl,flp,flm,llens,lcapdp,mype),
!$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(latdon,lat1,lat2,anl,flp,flm,llens,lcapdp,mype)
!mic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
!                                                                               
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
!
       call sph_grid2wave(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),          &
               llensd,lotas)
     enddo
!
! -----------------------------------------------------------------             
! ******** no multi-threads should be given in following loop ****!
! ******** otherwise non-produciable will occur. *****************!
! -----------------------------------------------------------------             
!
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
!                                                                               
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
       call sph_grid2wave(flp(1,1,kaps,lan),flm(1,1,kaps,lan),ZS,qww(1,lat),   &
                                           llstr,llensd,lwvdef,1)
       call sph_grid2wave(flp(1,1,kaus,lan),flm(1,1,kaus,lan),UUS,qww(1,lat),  &
                                           llstr,llensd,lwvdef,LEVSS)
       call sph_grid2wave(flp(1,1,kavs,lan),flm(1,1,kavs,lan),VVS,qww(1,lat),  &
                                           llstr,llensd,lwvdef,LEVSS)
       call sph_grid2wave(flp(1,1,kats,lan),flm(1,1,kats,lan),YS,qww(1,lat),   &
                                           llstr,llensd,lwvdef,LEVSS)
       call sph_grid2wave(flp(1,1,kars,lan),flm(1,1,kars,lan),RTS,qww(1,lat),  &
                                           llstr,llensd,lwvdef,LEVSS)
#endif
#ifdef DEFAULT
       call sph_comp_coeff2(flp(1,1,kaps,lan),flm(1,1,kaps,lan),ZS,qww(1,lat), &
                                           llstr,llensd,lwvdef,1)
       call sph_comp_coeff2(flp(1,1,kaus,lan),flm(1,1,kaus,lan),UUS,qww(1,lat),&
                                           llstr,llensd,lwvdef,LEVSS)
       call sph_comp_coeff2(flp(1,1,kavs,lan),flm(1,1,kavs,lan),VVS,qww(1,lat),&
                                           llstr,llensd,lwvdef,LEVSS)
       call sph_comp_coeff2(flp(1,1,kats,lan),flm(1,1,kats,lan),YS,qww(1,lat), &
                                           llstr,llensd,lwvdef,LEVSS)
       call sph_comp_coeff2(flp(1,1,kars,lan),flm(1,1,kars,lan),RTS,qww(1,lat),&
                                           llstr,llensd,lwvdef,LEVSS)
#endif
!                                                                               
       call sph_sum_wind_coeff(flp(1,1,kaus,lan),flm(1,1,kaus,lan),            &
                 flp(1,1,kavs,lan),flm(1,1,kavs,lan),                          &
                 anltop(1,1,1),anltop(1,1,LEVSS+1),                            &
                 qvv(1,lat),wgt(lat),                                          &
                 llstr,llensd,lwvdef,LEVSS)
     enddo
#undef ZS
#undef UUS
#undef VVS
#undef YS
#undef RTS
!
#ifndef MP
     latdon=latdon+(lat2-lat1+1)
   enddo
!
#endif
!                                                                               
#ifdef ORIGIN_THREAD
!$doacross share(spdmax,spdlat,lats2),local(k,lat)
#endif
#ifdef CRAY_THREAD
!mic$ do all                                                                    
!mic$1 shared(spdmax,spdlat,lats2)
!mic$1 private(k,lat)                                                           
#endif
#ifdef OPENMP
!$omp parallel do private(k,lat)
#endif
   do k = 1,levs_
     spdmax(k) = 0.0
     do lat = 1,lats2
       spdmax(k)=max(spdmax(k),spdlat(k,lat))
     enddo                                  
     spdmax(k)=sqrt(spdmax(k))             
   enddo                                  
#ifdef MP
   call mpgetspd(spdmax)
#endif
   if( iope ) then
     write(6,100)(spdmax(k),k=1,levs_)
100  format(' global checked speed maxima for all layers ',                    &
         :/' spdmx(01:10)=',10f5.0,:/' spdmx(11:20)=',10f5.0,                  &
         :/' spdmx(21:30)=',10f5.0,:/' spdmx(31:40)=',10f5.0,                  &
         :/' spdmx(41:50)=',10f5.0,:/' spdmx(51:60)=',10f5.0,                  &
         :/' spdmx(61:70)=',10f5.0,:/' spdmx(71:80)=',10f5.0,                  &
         :/' spdmx(81:90)=',10f5.0,:/' spdmx(91:00)=',10f5.0)
   endif
!                                                                               
!     input : w=d(u)/d(t) x=d(v)/d(t)                                           
!     output: uln=d(di)/d(t) vln=d(ze)/d(t)                                     
!                                                                               
#ifdef MP
#define UUS uua
#define VVS vva
#define XS xa
#define WS wa
#else
#define UUS uu
#define VVS vv
#define XS x
#define WS w
#endif
   call sph_wind2divor(UUS,VVS,XS,WS,anltop(1,1,1),                            &
               anltop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef UUS
#undef VVS
#undef XS
#undef WS
!
#ifdef MP
   call mpnk2nn(ya ,lln22p_,levsp_, y ,lnt22p_,levs_,1)
   call mpnk2nn(rta,lln22p_,levsp_, rt,lnt22p_,levs_,ntotal_)
   call mpnk2nn(xa ,lln22p_,levsp_, x ,lnt22p_,levs_,1)
   call mpnk2nn(wa ,lln22p_,levsp_, w ,lnt22p_,levs_,1)
   call mpn2nn(za,lln22p_,z,lnt22p_,1)
#endif
!
!  subtract off linear dependence on divergence                                 
!
   do k = 1,levs_ 
     do j = 1,levs_ 
       do i = 1,LNT2S
         y(i,k)=y(i,k)-bm(k,j)*di(i,j)
       enddo
     enddo
   enddo
!
! move div tendency into x and add topog. contrib.
! integrate vorticity amd moisture in time
! remember uln is old x
! remember vln is old w
!
   do k = 1,levs_
     do i = 1,LNT2S
       x(i,k)=x(i,k)+gz(i)
       w(i,k)=zem(i,k)+2.*deltim*w(i,k)
     enddo
   enddo
   do k = 1,levs_
     do i = 1,LNT2S 
       rt(i,k)= rm(i,k)+2.*deltim* rt(i,k)
     enddo
   enddo
!
   if( iope ) then
     do k = 1,levs_
       w(1,k)=0.e0
       w(2,k)=0.e0
     enddo
   endif
!
! ---------------- for phys_main_solver ------------------
!
#ifdef DG
   write(6,'(" dyn_driver_onestep t range ",2(4x,f6.1," @i,k,lat ",3i4))')     &
      tgmx,igmx,kgmx,jgmx,tgmn,igmn,kgmn,jgmn
#endif
!
! ---------------- for phys_rad_solver -------------------
!
   if(inistp.ne.0) then
     runrad=.false.
   else
     runrad=.true.
   endif
!
   deallocate(fnl)
#ifdef MP
   deallocate(syf, grs, dyf, dgr, anf, gra)
#endif
   deallocate(scs2, syn, syntop, dyn, anl, anltop, flp, flm)
!
   return                                                                    
   end subroutine dyn_driver_onestep
!
!-------------------------------------------------------------------------------
