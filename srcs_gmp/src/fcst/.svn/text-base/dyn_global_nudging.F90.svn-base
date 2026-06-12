#include <define.h>
   subroutine dyn_global_nudging
!-------------------------------------------------------------------------------
#ifdef MP
   use paramodel, only : JCAP1S,LATG2S,LEVSS,LEVHS,LCAPS,LCAP22S,              &
                         LONF22S,LNT2S,levs_,levh_,lnt22_,ncpus_,              &
                         jcap1_,lonf22_,lonf2_,latg2_,ntotal_
   use paramodel, only : lonf22p_,levhp_,lnt2p_,jcap1p_,lcapp_,                &
                         lnt22p_,lln22p_,levsp_,lcap22p_,latg2p_
#ifdef GMPDAMP
   use paramodel, only : gdamp_crlat_
#endif
   use constant, only  : cp_,pi_,g_,rd_,rv_,rerth_
   use comfibm
   use commpi
   use compspec
   use comfcst, only :  spdlat
#ifdef REDUCE_GRID
   use comreduce
#endif
#undef GDAMP09
!-------------------------------------------------------------------------------
!
! subprogram:    dyn_global_nudging  (only for MPI) 
!  - computes dynamic non-linear tendency terms
!    of temp. div. ln(ps) (for global dynamical downscaling)
!  - computes predicted values of vorticity and moisture
!
! abstract:
!   program nudges the Fourier coefficients of forecasted run ("fcst data")
!   with those of reanalysis (or other external) run (called "base data").
!   program starts with the conversion of spectrum base data to 
!   fourier series. this conversion is taken place in every 6 hours
!   (nesthour). the base data is linearly interpolated for the time
!   in between. the spectra of fcst data are also converted to
!   fourier series.
!   two kinds of fourier series are nudged with a weighting constant as,
!    fnew=(1/(a+1))*(fold+a*fbase)
!   for those with longer wavelength than lnest (e.g. 1000km)
!   (for those with shorter, a=0).
!   a can be variable, too (attenuating function).
!   all of the variables, t,q,u,v,ps,gz can be nudged, but practically,
!   ps,u,v and t should be nudged. instead, zonal q (the first coefficient
!   of the fourier series) is replacable to that of the base data.
!
! default setting:
!   u,v,t,ps : nudging
!   q        : nothing (used to be zonal average replace)
!   n-scale  : 2000 km
!   n-weight : 0.9 (indaa=0)
!              used to be attenuating function (indaa=1)
!
! program history log:
!   07-07-24  kei yoshimura             modified from dyn_sph_driver
!
! references :
!   yoshimura, k. and m. kanamitsu (2008, mon. wea. rev.)
!
!  ::: structure :::
!
!    [dyn_global_nudging] --- [rmp_base_filename]
!
!-------------------------------------------------------------------------------
!
!  version with stacked transforms
!
!  sof include
!
! syn(1, 0*levs_+1, lan)  ze
! syn(1, 1*levs_+1, lan)  di
! syn(1, 2*levs_+1, lan)  te
! syn(1, 3*levs_+1, lan)  rq
! syn(1, 4*levs_+1, lan)  uln
! syn(1, 5*levs_+1, lan)  vln
! syn(1, 6*levs_+1, lan)  dpdphi
! syn(1, 6*levs_+2, lan)  dpdlam
! syn(1, 6*levs_+3, lan)  q
!
   integer              ::  lots    ,lotst
   integer              ::  ksz     
   integer              ::  ksd     
   integer              ::  kst     
   integer              ::  ksr     
   integer              ::  ksu     ,kstb
   integer              ::  ksv     
   integer              ::  kspphi  
   integer              ::  ksplam  
   integer              ::  ksp     
   integer              ::  ksplap  
!
   integer              ::  lotd    
   integer              ::  kdtphi  
   integer              ::  kdrphi  
   integer              ::  kdtlam  
   integer              ::  kdrlam  
   integer              ::  kdulam  
   integer              ::  kdvlam  
   integer              ::  kduphi  
   integer              ::  kdvphi  
   integer              ::  lota    ,lotat
   integer              ::  kap     
   integer              ::  kau     
   integer              ::  kav     
   integer              ::  kat     
   integer              ::  kar     
!
   integer              ::  lotss   ,lotsts
   integer              ::  kszs    
   integer              ::  ksds    
   integer              ::  ksts    
   integer              ::  ksrs    
   integer              ::  ksus    ,kstbs
   integer              ::  ksvs    
   integer              ::  kspphis 
   integer              ::  ksplams 
   integer              ::  ksps    
   integer              ::  ksplaps 
!
   integer              ::  lotds   
   integer              ::  kdtphis 
   integer              ::  kdrphis 
   integer              ::  kdtlams 
   integer              ::  kdrlams 
   integer              ::  kdulams 
   integer              ::  kdvlams 
   integer              ::  kduphis 
   integer              ::  kdvphis 
!
#ifdef GDAMP09
   integer              ::  lotas   ,lotats
#else
   integer              ::  lotas   ,lotats
#endif
   integer              ::  kaps    
   integer              ::  kaus    
   integer              ::  kavs    
   integer              ::  kats    
   integer              ::  kars    
#ifdef GDAMP09
   integer              ::  kads    
   integer              ::  kazs    
#endif
!
! local array
!
   real, allocatable    ::  syf(:,:,:),grs(:,:,:), anf(:,:,:)
#define NCPUSS latg2_
   real, allocatable    ::  scs2(:) ,                                          &
                            syn(:,:,:),syntop(:,:,:),                          &
                            anl(:,:,:),anltop(:,:,:),                          &
                            flp(:,:,:,:),flm(:,:,:,:)
!
! by Kei
!
   character(len=20)    ::  bfname,bfname2
   real                 ::  fnxx,sshour,nestsec,sshour2,fnxx2,fnxb
   data    sshour,sshour2 /0.,0./
   data    fnxb/-1./
   real                 ::  fhour1,fhour2,nestsec2
   integer              ::  nchar,idate1(4),idate2(4),nbase1,nbase2
   data    nbase1,nbase2/91,92/
   data    nchar/15/
   integer              ::  ifirst
   data    ifirst/0/
   integer              ::  m,stm,edm,ii,jj1,jj2
!
   integer,allocatable,save :: jlim(:)
   integer,allocatable,save :: jlimt(:)
   real   ,allocatable,save :: nlat(:)
   real   ,allocatable,save :: nlatb(:)
!
   save    ifirst,sshour,fnxb,bfname2
   real                 ::  lnest,nesthour, aa   !! nudging scale, interval, weight
   data    lnest,nesthour, aa /2000.e3,  6.  ,0.9/
#ifdef GDAMP09
   real                 ::  lnestt               !! nudging scale for temp.
   data    lnestt /10000.e3/
#endif
   real                 ::  aa2
   integer              ::  jnest,indaa    !! basedata resolution, use of AF for aa
   data    jnest,indaa /   62,    0/     !! indaa=0: use constant value.
!
   real     ::  gz1     (lnt22p_),                                             &
                 q1     (lnt22p_),                                             &
                ze1     (lnt22p_,levs_),                                       &
                di1     (lnt22p_,levs_),                                       &
                te1     (lnt22p_,levs_),                                       &
                rq1     (lnt22p_,levh_),                                       &
                z001
   real     ::  gz2     (lnt22p_),                                             &
                 q2     (lnt22p_),                                             &
                ze2     (lnt22p_,levs_),                                       &
                di2     (lnt22p_,levs_),                                       &
                te2     (lnt22p_,levs_),                                       &
                rq2     (lnt22p_,levh_),                                       &
                z002
   real     ::  bzea1    (lln22p_,levsp_),                                     &
                bdia1    (lln22p_,levsp_),                                     &
                btea1    (lln22p_,levsp_),                                     &
                brqa1    (lln22p_,levhp_),                                     &
                bulna1   (lln22p_,levsp_),                                     &
                bvlna1   (lln22p_,levsp_),                                     &
                bdpdphia1(lln22p_),                                            &
                bdpdlama1(lln22p_),                                            &
                bqa1     (lln22p_),                                            &
                bqlapa1  (lln22p_)
   real     ::  bzea2    (lln22p_,levsp_),                                     &
                bdia2    (lln22p_,levsp_),                                     &
                btea2    (lln22p_,levsp_),                                     &
                brqa2    (lln22p_,levhp_),                                     &
                bulna2   (lln22p_,levsp_),                                     &
                bvlna2   (lln22p_,levsp_),                                     &
                bdpdphia2(lln22p_),                                            &
                bdpdlama2(lln22p_),                                            &
                bqa2     (lln22p_),                                            &
                bqlapa2  (lln22p_)
!
   real, allocatable, dimension(:,:,:)         ::  byf
   real, allocatable, save, dimension(:,:,:)   ::  byf1,byf2
   real, allocatable, dimension(:,:,:)         ::  byn1,byn2
   real, allocatable, dimension(:,:,:)         ::  byntop1, byntop2
!
   real             ::  ditmp(lnt22p_,levs_)
   real             ::  ek(lonf22p_,levs_)      
   real             ::  crcolat
!-------------------------------------------------------------------------------
   lots    =5*levs_+levh_+4
   lotst   =2*levs_+1
   ksz     =1
   ksd     =1*levs_+1
   kst     =2*levs_+1
   ksr     =3*levs_+1
   ksu     =3*levs_+levh_+1
   kstb    =3*levs_+levh_+1
   ksv     =4*levs_+levh_+1
   kspphi  =5*levs_+levh_+1
   ksplam  =5*levs_+levh_+2
   ksp     =5*levs_+levh_+3
   ksplap  =5*levs_+levh_+4
!
   lotd    =6*levs_+2*levh_
   kdtphi  =0*levs_+1
   kdrphi  =1*levs_+1
   kdtlam  =1*levs_+levh_+1
   kdrlam  =2*levs_+levh_+1
   kdulam  =2*levs_+2*levh_+1
   kdvlam  =3*levs_+2*levh_+1
   kduphi  =4*levs_+2*levh_+1
   kdvphi  =5*levs_+2*levh_+1
   lota    =3*levs_+levh_+1
   lotat   =2*levs_
   kap     =1
   kau     =2
   kav     =1*levs_+2
   kat     =2*levs_+2
   kar     =3*levs_+2
!
   lotss   =5*LEVSS+LEVHS+4
   lotsts  =2*LEVSS+1
   kszs    =1
   ksds    =1*LEVSS+1
   ksts    =2*LEVSS+1
   ksrs    =3*LEVSS+1
   ksus    =3*LEVSS+LEVHS+1
   kstbs   =3*LEVSS+LEVHS+1
   ksvs    =4*LEVSS+LEVHS+1
   kspphis =5*LEVSS+LEVHS+1
   ksplams =5*LEVSS+LEVHS+2
   ksps    =5*LEVSS+LEVHS+3
   ksplaps =5*LEVSS+LEVHS+4
!
   lotds   =6*LEVSS+2*LEVHS
   kdtphis =0*LEVSS+1
   kdrphis =1*LEVSS+1
   kdtlams =1*LEVSS+LEVHS+1
   kdrlams =2*LEVSS+LEVHS+1
   kdulams =2*LEVSS+2*LEVHS+1
   kdvlams =3*LEVSS+2*LEVHS+1
   kduphis =4*LEVSS+2*LEVHS+1
   kdvphis =5*LEVSS+2*LEVHS+1
!
#ifdef GDAMP09
   lotas   =5*LEVSS+LEVHS+1
   lotats  =2*LEVSS
#else
   lotas   =3*LEVSS+LEVHS+1
   lotats  =2*LEVSS
#endif
   kaps    =1
   kaus    =2
   kavs    =1*LEVSS+2
   kats    =2*LEVSS+2
   kars    =3*LEVSS+2
#ifdef GDAMP09
   kads    =3*LEVSS+LEVHS+2
   kazs    =4*LEVSS+LEVHS+2
#endif
   allocate(                                                                   &
      syf(lonf22_,lotss,latg2p_),grs(lonf22p_,lots,latg2p_),                   &
      anf(lonf22_,lotas,latg2p_) )
#define NCPUSS latg2_
   allocate(                                                                   &
      scs2(JCAP1S),                                                            &
      syn(LCAP22S,lotss,NCPUSS),syntop(2,JCAP1S,lotsts),                       &
      anl(LCAP22S,lotas,NCPUSS),anltop(2,JCAP1S,lotats),                       &
      flp(2,JCAP1S,lotas,NCPUSS),flm(2,JCAP1S,lotas,NCPUSS) )
!
   allocate(byf(lonf22_,lotss,latg2p_))
   if(.not.allocated(byf1)) allocate(byf1(lonf22_,lotss,latg2p_) )
   if(.not.allocated(byf2)) allocate(byf2(lonf22_,lotss,latg2p_) )
!
   allocate(byn1(lcap22p_,lotss,latg2_),byn2(lcap22p_,lotss,latg2_) )
   allocate(byntop1(2,jcap1p_,lotsts), byntop2(2,jcap1p_,lotsts))
!
   if(.not.allocated(jlim))  allocate(jlim(latg2p_))
   if(.not.allocated(jlimt)) allocate(jlimt(latg2p_))
   if(.not.allocated(nlat))  allocate(nlat(latg2_))
   if(.not.allocated(nlatb)) allocate(nlatb(latg2p_))
!
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   lons2=lonlen(mype)*2
   lats2=latlen(mype)
!
!=========================================================================
! GMP : calculate dynamic tendency terms
!=========================================================================
!
!  KEI    
!  nudging by base field
!  read base
!
   crcolat = pi_/2.-gdamp_crlat_*pi_/180.
   print*, 'iope', iope
   nestsec=nesthour*60.*60.
   sshour2=sshour
   if (ifirst.eq.0) then
     if (iope) then
       print*,'dyn_global_nudging, lnest,aa,indaa',lnest,aa,indaa
       print*,' critical lat (deg) colat (rad) ',gdamp_crlat_,crcolat
     endif
     if (shour.gt.1800.) then
       sshour=shour+(thour-1.)*60.*60.
       sshour2=(thour-1.)*60.*60.
     else
       sshour=shour+(thour)*60.*60.
       sshour2=(thour)*60.*60.
     endif
!
     do j = 1,latg2_
!
! only tropics (latitudinal weighting)
!
       if (colrad(j).ge.crcolat) then
         nlat(j)=1.
       else
         nlat(j)=0.
       endif
     enddo
!
     do j = 1,latlen(mype)
       jj=latdef(latstr(mype)+j-1)
       nlatb(j)=nlat(jj)
!
! with latitudinal weighting
!
       jlim(j)=min(int((2*pi_*rerth_*sin(colrad(jj))                           &
                /lnest)*nlatb(j)),jnest)
!       jlim(j)=min(int(2*pi_*rerth_*sin(colrad(jj))                            &
!                /lnest),jnest)
#ifdef GDAMP09
       jlimt(j)=min(int(2*pi_*rerth_*sin(colrad(jj))                           &
                /lnestt),jnest)
#endif
     enddo
     lat1=jstr
     lat2=jend
     latdon=jstr-1
     do lat = lat1,lat2
       lan=lat-latdon
       if(iope) then
         print 104,'JJ',mype,lat,jlim(lan),nlatb(lan),colrab(lan)
       endif
     enddo
!
     fnxx=int(sshour2/nestsec)*nesthour
     if(iope) then
       print*, 'call rmp_base_filename1' 
       call rmp_base_filename(idate(4),idate(2),idate(3),idate(1),             &
              fnxx,                                                            &
              bfname,nchar)
       close(nbase1)
       open(unit=nbase1,file=bfname(1:nchar),status='unknown',                 &
              form='unformatted')
     endif
!
     call file_read_sigma(nbase1,fhour1,idate1,gz1,q1,te1,di1,ze1,rq1,         & 
#ifndef HYBRID
             sl,si,z001)
#else
             ak5,bk5,z001)
#endif
     print*, 'after file_read_sigma1'
!
! ... without divergence component ...
!
!        do j = 1,lnt22p_
!          do k = 1,levs_
!            di1(j,k)=0.
!          enddo
!        enddo
!
   else
     sshour=sshour+deltim
   endif
!
104 format (a5,3i4,f6.2,f15.7)
   fnxx=int(sshour2/nestsec)*nesthour
   fnxx2=fnxx+nesthour
   if (iope) then
     print*, 'fnxb,fnxx', fnxb, fnxx
   endif
!
!   read basefiles (sigma-level spectral coef.)
!
   if (fnxb.ne.fnxx) then
!
     if (iope) then
       print*, 'call rmp_base_filename' 
       call rmp_base_filename(idate(4),idate(2),idate(3),idate(1),             &
           fnxx2,                                                              &
           bfname2,nchar)
       close(nbase2)
       open(unit=nbase2,file=bfname2(1:nchar),status='unknown',                &
           form='unformatted')
     endif
!
!        call rdsig2(nbase2,fhour2,idate2,gz2,q2,te2,di2,ze2,rq2,
     call file_read_sigma(nbase2,fhour2,idate2,gz2,q2,te2,di2,ze2,rq2,         &
#ifndef HYBRID
             sl,si,z002)
#else
             ak5,bk5,z002)
#endif
     print*, 'after file_read_sigma2'
!
! ... without divergence component ...
!
!        do j = 1,lnt22p_
!          do k = 1,levs_
!            di2(j,k)=0.
!          enddo
!        enddo
   endif
!     if (iope) then
   print 113,'Kei',fnxx,fnxx2,sshour,sshour2,bfname2,thour,fnxb
!     endif
113 format(a5,2f8.2,2f10.0,' ',a15,2f8.2)
!
! ... without divergence component ...
!
!      do j = 1,lnt22p_
!         do k = 1,levs_
!            ditmp(j,k)=di(j,k)
!            di(j,k)=0.
!         enddo
!      enddo
!
   call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,1)
   call mpnn2nk(di ,lnt22p_,levs_, dia,lln22p_,levsp_,1)
   call mpnn2nk(te ,lnt22p_,levs_, tea,lln22p_,levsp_,1)
   call mpnn2nk(rq ,lnt22p_,levs_, rqa,lln22p_,levsp_,ntotal_)
   call mpnn2n(q ,lnt22p_, qa,lln22p_,1)
   call sph_del_log_ps(qa,dpdphia,syntop(1,1,2*levsp_+1),dpdlama,              &
                llstr,llens,lwvdef)
   call sph_solve_laplacian(qa,qlapa,llstr,llens,lwvdef)
   call sph_divor2wind(dia,zea,ulna,vlna,syntop(1,1,1),                        &
                syntop(1,1,levsp_+1),llstr,llens,lwvdef)
!
!  for base cc
!
   if (fnxb.ne.fnxx) then
!
     if (ifirst.eq.0) then
       call mpnn2nk(ze1 ,lnt22p_,levs_, bzea1,lln22p_,levsp_,1)
       call mpnn2nk(di1 ,lnt22p_,levs_, bdia1,lln22p_,levsp_,1)
       call mpnn2nk(te1 ,lnt22p_,levs_, btea1,lln22p_,levsp_,1)
       call mpnn2nk(rq1 ,lnt22p_,levs_, brqa1,lln22p_,levsp_,ntotal_)
       call mpnn2n(q1 ,lnt22p_, bqa1,lln22p_,1)
       call sph_del_log_ps(bqa1,bdpdphia1,byntop1(1,1,2*levsp_+1),bdpdlama1,   &
                llstr,llens,lwvdef)
       call sph_solve_laplacian(bqa1,bqlapa1,llstr,llens,lwvdef)
       call sph_divor2wind(bdia1,bzea1,bulna1,bvlna1,byntop1(1,1,1),           &
                byntop1(1,1,levsp_+1),llstr,llens,lwvdef)
     endif
!
     call mpnn2nk(ze2 ,lnt22p_,levs_, bzea2,lln22p_,levsp_,1)
     call mpnn2nk(di2 ,lnt22p_,levs_, bdia2,lln22p_,levsp_,1)
     call mpnn2nk(te2 ,lnt22p_,levs_, btea2,lln22p_,levsp_,1)
     call mpnn2nk(rq2 ,lnt22p_,levs_, brqa2,lln22p_,levsp_,ntotal_)
     call mpnn2n(q2 ,lnt22p_, bqa2,lln22p_,1)
     call sph_del_log_ps(bqa2,bdpdphia2,byntop2(1,1,2*levsp_+1),bdpdlama2,     &
                 llstr,llens,lwvdef)
     call sph_solve_laplacian(bqa2,bqlapa2,llstr,llens,lwvdef)
     call sph_divor2wind(bdia2,bzea2,bulna2,bvlna2,byntop2(1,1,1),             &
                 byntop2(1,1,levsp_+1),llstr,llens,lwvdef)
   endif
!
   do n = 1,2
     do j = 1,jcap1p_
       do l = 1,lotats
         anltop(n,j,l)=0.0
       enddo
     enddo
   enddo
!
   lat1=1
   lat2=latg2_
   latdon=0
!
! first lat loop
!
!$omp parallel do private(lat,lan,llensd)
!
   do 1000 lat = lat1,lat2
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
     call sph_sum_coeff(zea    ,syn(1,kszs,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff(dia    ,syn(1,ksds,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff(tea    ,syn(1,ksts,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff(rqa    ,syn(1,ksrs,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVHS)
     call sph_sum_coeff(ulna   ,syn(1,ksus,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff(vlna   ,syn(1,ksvs,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff(dpdphia,syn(1,kspphis,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
     call sph_sum_coeff(dpdlama,syn(1,ksplams,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
     call sph_sum_coeff(qa     ,syn(1,ksps,lan)   ,                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
     call sph_sum_coeff(qlapa  ,syn(1,ksplaps,lan),                            &
                                          qtt(1,lat),llstr,llensd,lwvdef,1)
     call sph_sum_coeff_top(syn(1,ksus,lan)    ,syntop(1,1,1)        ,         &
                                          qvv(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff_top(syn(1,ksvs,lan)    ,syntop(1,1,LEVSS+1)  ,         &
                                          qvv(1,lat),llstr,llensd,lwvdef,LEVSS)
     call sph_sum_coeff_top(syn(1,kspphis,lan),syntop(1,1,2*LEVSS+1),          &
                                          qvv(1,lat),llstr,llensd,lwvdef,1)
!
! for base cc
!
     if (fnxb.ne.fnxx) then
!
       if (ifirst.eq.0) then
!
         call sph_sum_coeff(bzea1,byn1(1,kszs,lan),qtt(1,lat),                 &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff(bdia1,byn1(1,ksds,lan),qtt(1,lat),                 &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff(btea1,byn1(1,ksts,lan),qtt(1,lat),                 &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff(brqa1,byn1(1,ksrs,lan),qtt(1,lat),                 &
                     llstr,llensd,lwvdef,levhp_)
         call sph_sum_coeff(bulna1,byn1(1,ksus,lan),qtt(1,lat),                &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff(bvlna1,byn1(1,ksvs,lan),qtt(1,lat),                &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff(bqa1,byn1(1,ksps,lan),qtt(1,lat),                  &
                     llstr,llensd,lwvdef,1)
         call sph_sum_coeff_top(byn1(1,ksus,lan),byntop1(1,1,1),qvv(1,lat),    &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff_top(byn1(1,ksvs,lan),byntop1(1,1,1),qvv(1,lat),    &
                     llstr,llensd,lwvdef,levsp_)
         call sph_sum_coeff_top(byn1(1,kspphis,lan),byntop1(1,1,1),qvv(1,lat), &
                       llstr,llensd,lwvdef,1)
       endif
!
       call sph_sum_coeff(bzea2,byn2(1,kszs,lan),qtt(1,lat),                   &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff(bdia2,byn2(1,ksds,lan),qtt(1,lat),                   &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff(btea2,byn2(1,ksts,lan),qtt(1,lat),                   &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff(brqa2,byn2(1,ksrs,lan),qtt(1,lat),                   &
                   llstr,llensd,lwvdef,levhp_)
       call sph_sum_coeff(bulna2,byn2(1,ksus,lan),qtt(1,lat),                  &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff(bvlna2,byn2(1,ksvs,lan),qtt(1,lat),                  &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff(bqa2,byn2(1,ksps,lan),qtt(1,lat),                    &
                   llstr,llensd,lwvdef,1)
       call sph_sum_coeff_top(byn2(1,ksus,lan),byntop2(1,1,1),qvv(1,lat),      &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff_top(byn2(1,ksvs,lan),byntop2(1,1,1),qvv(1,lat),      &
                   llstr,llensd,lwvdef,levsp_)
       call sph_sum_coeff_top(byn2(1,kspphis,lan),byntop2(1,1,1),qvv(1,lat),   &
                   llstr,llensd,lwvdef,1)
     endif
!       
1000 continue
!
! compute merid. derivs. of temp. and moisture using qdd.
!
   call mpnl2ny(syn,lcap22p_,latg2_,syf,lonf22_,latg2p_,lotss,1,lotss)
!
! for base cc
!
   if (fnxb.ne.fnxx) then
     if (ifirst.eq.0) then
       ifirst=1
       call mpnl2ny(byn1,lcap22p_,latg2_,byf1,lonf22_,latg2p_,lotss,1,lotss)
     else
       do i = 1,lonf22_
         do k = 1,lotss
           do j = 1,latg2p_
             byf1(i,k,j)=byf2(i,k,j)
           enddo
         enddo
       enddo
     endif
     call mpnl2ny(byn2,lcap22p_,latg2_,byf2,lonf22_,latg2p_,lotss,1,lotss)
   endif
!
   do i = 1,lonf22_
     do k = 1,lotss
       do j = 1,latg2p_
         byf(i,k,j)                                                            &
              =byf1(i,k,j)*(1.-(sshour-fnxx*3600.)/nestsec)                    &
              +byf2(i,k,j)*    (sshour-fnxx*3600.)/nestsec
       enddo
     enddo
   enddo
!
! nudging cc
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
!
   do lat = lat1,lat2
     lan=lat-latdon
     do n = 1,2
       do i = 1,(jlim(lan)+1)*2
         ii=(n-1)*(lonf22_-2)/2+i
         if (indaa.eq.0) then
           aa2=aa/(1.+aa)
         else
           aa2=(1.-real(int((i-1)/2)/max(real(jlim(lan)),1.)))**0.5
!                                 alpha infinity->0 sqrt (YK07)
!              aa2=(cos(int((i-1)/2)/max(jlim(lan),1)*pi_)+1.)/2.
!                                 alpha Infinity->0 cosine
         endif
!                                 aa2 is transformed in aa/(1+aa)
         do k = 1,levsp_
#ifdef GDAMP09  
           syf(ii,kszs+k-1,lan)                                                &
               =syf(ii,kszs+k-1,lan)*(1.-aa2)                                  &
               +byf(ii,kszs+k-1,lan)*(   aa2) !! ze
           syf(ii,ksds+k-1,lan)                                                &
               =syf(ii,ksds+k-1,lan)*(1.-aa2)                                  &
               +byf(ii,ksds+k-1,lan)*(   aa2) !! di
#else
           syf(ii,ksus+k-1,lan)                                                &
               =syf(ii,ksus+k-1,lan)*(1.-aa2)                                  &
               +byf(ii,ksus+k-1,lan)*(   aa2) !! u
           syf(ii,ksvs+k-1,lan)                                                &
               =syf(ii,ksvs+k-1,lan)*(1.-aa2)                                  &
               +byf(ii,ksvs+k-1,lan)*(   aa2) !! v
           syf(ii,ksts+k-1,lan)                                                & 
               =syf(ii,ksts+k-1,lan)*(1.-aa2)                                  &
               +byf(ii,ksts+k-1,lan)*(   aa2) !! t
#endif
         enddo
#ifndef GDAMP09  
         syf(ii,ksps,lan)                                                      &
             =syf(ii,ksps,lan)*(1.-aa2)                                        &
             +byf(ii,ksps,lan)*(   aa2) !! q (sfP)
#endif
       enddo
#ifdef GDAMP09  
       do i = 1,(jlimt(lan)+1)*2
         ii=(n-1)*(lonf22_-2)/2+i
         aa2=aa/(1.+aa)
         do k = 1,levsp_
           syf(ii,ksts+k-1,lan)                                                &
               =syf(ii,ksts+k-1,lan)*(1.-aa2)                                  &
               +byf(ii,ksts+k-1,lan)*(   aa2) !! t
         enddo
       enddo
#endif
!          do i = 1,(jlim(lan)+1)*2 !! for same as UV-SSBC
       do i = 1,2              !! for zonal mean nudging
         ii=(n-1)*(lonf22_-2)/2+i
!            do k = 1,levsp_
!              syf(ii,ksts+k-1,lan)=byf(ii,ksts+k-1,lan) !! te(tmp)
!            enddo
!            do k = 1,levhp_
!              syf(ii,ksrs+k-1,lan)=byf(ii,ksrs+k-1,lan) !! rq(vpr)
!            enddo
!            syf(ii,ksps,lan)=byf(ii,ksps,lan)           !! q (sfP)
       enddo
     enddo
   enddo
!
#ifdef DBG
!      if (mype.eq.master) then
!         do i = 1,lonf22_
!            print 102,'syf2',i,ksps,syf(i,ksps,10),byf(i,ksps,10)
!         enddo
!         do i = 1,lonf22_
!            print 102,'syf2',i,ksus,syf(i,ksus,10),byf(i,ksus,10)
!         enddo
!         do i = 1,lonf22_
!            print 102,'syf2',i,ksvs,syf(i,ksvs,10),byf(i,ksvs,10)
!         enddo
!         do i = 1,lonf22_
!            print 102,'syf2',i,ksts,syf(i,ksts,10),byf(i,ksts,10)
!         enddo
!         do i = 1,lonf22_
!            print 102,'syf2',i,ksrs,syf(i,ksrs,10),byf(i,ksrs,10)
!         enddo
!      endif
102  format (a5,2i5,2e15.7)
#endif
!
!$omp parallel do private(j,k)
!
   do k = 1,levsp_
     do j = 1,lln22p_
       uua(j,k)=0.0
       vva(j,k)=0.0
       ya(j,k)=0.0
#ifdef GDAMP09
       xa(j,k)=0.0
       wa(j,k)=0.0
#endif
     enddo
   enddo
!
   do k = 1,levhp_
     do j = 1,lln22p_
       rta(j,k)=0.0
     enddo
   enddo
   do j = 1,lln22p_
     za(j,1)=0.e0
   enddo
!
!$omp parallel do private(i,j,k)
!
   do i = 1,lonf22_
     do j = 1,latg2p_
       anf(i,kaps,j)=syf(i,ksps,j)            !! q (sfP)
       do k = 1,levsp_
         anf(i,kaus-1+k,j)=syf(i,ksus-1+k,j) !! u
         anf(i,kavs-1+k,j)=syf(i,ksvs-1+k,j) !! v
         anf(i,kats-1+k,j)=syf(i,ksts-1+k,j) !! te (tmp)
#ifdef GDAMP09
         anf(i,kads-1+k,j)=syf(i,ksds-1+k,j) !! di
         anf(i,kazs-1+k,j)=syf(i,kszs-1+k,j) !! ze
#endif
       enddo
       do k = 1,levhp_
         anf(i,kars-1+k,j)=syf(i,ksrs-1+k,j) !! rt (vpr)
       enddo
     enddo
   enddo
!
   call mpny2nl(anf,lonf22_ ,latg2p_,                                          &
                   anl,lcap22p_,latg2_ ,lotas,kaps,lotas)
!
   lat1=1
   lat2=latg2_
   latdon=0
!
!$omp parallel do private(lat,lan,llensd)
!
   do 2500 lat = lat1,lat2
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
     do k = 1,levsp_*2
       kaa=kaus-1+k
       do j = 1,llensd*2
         anl(       j,kaa,lan)=anl(       j,kaa,lan)*rcs2(lat)
         anl(lcapp_+j,kaa,lan)=anl(lcapp_+j,kaa,lan)*rcs2(lat)
       enddo
     enddo
!
     call sph_grid2wave(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),            &
                  llensd,lotas)
   2500 continue
! 
!   no multi-threads should be given in following loop
!   otherwise results will be non-produceable
!
   do 3000 lat = lat1,lat2
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
     call sph_comp_coeff2(flp(1,1,kaps,lan),flm(1,1,kaps,lan),za ,qww(1,lat),  &
                                  llstr,llensd,lwvdef,1)
     call sph_comp_coeff2(flp(1,1,kaus,lan),flm(1,1,kaus,lan),uua,qww(1,lat),  &
                                  llstr,llensd,lwvdef,LEVSS)
     call sph_comp_coeff2(flp(1,1,kavs,lan),flm(1,1,kavs,lan),vva,qww(1,lat),  &
                                  llstr,llensd,lwvdef,LEVSS)
     call sph_comp_coeff2(flp(1,1,kats,lan),flm(1,1,kats,lan),ya ,qww(1,lat),  &
                                  llstr,llensd,lwvdef,LEVSS)
     call sph_comp_coeff2(flp(1,1,kars,lan),flm(1,1,kars,lan),rta,qww(1,lat),  &
                                  llstr,llensd,lwvdef,LEVHS)
#ifdef GDAMP09
     call sph_comp_coeff2(flp(1,1,kads,lan),flm(1,1,kads,lan),xa ,qww(1,lat),  &
                                  llstr,llensd,lwvdef,LEVSS)
     call sph_comp_coeff2(flp(1,1,kazs,lan),flm(1,1,kazs,lan),wa ,qww(1,lat),  &
                                  llstr,llensd,lwvdef,LEVSS)
#else
     call sph_sum_wind_coeff(flp(1,1,kaus,lan),flm(1,1,kaus,lan),              &
                                  flp(1,1,kavs,lan),flm(1,1,kavs,lan),         &
                                  anltop(1,1,1),anltop(1,1,levsp_+1),          &
                                  qvv(1,lat),wgt(lat),                         &
                                  llstr,llensd,lwvdef,levsp_)
#endif
!
   3000 continue
!
#ifdef DBG
!   Check UV (Unnecessary) !!!
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
!
   do 1300 lat = lat1,lat2
#ifdef REDUCE_GRID
     lcapf=lcapd(latdef(lat))
     lonff=lonfd(latdef(lat))
#endif
     lan=lat-latdon
     call sph_wave2grid(syf(1,1,lan),syf(1,1,lan),lotss*2,                     &
                lcapf,lonff,latdef(lat),1)
   1300 continue
!
   call mpnk2nx(syf,lonf22_,lotss,                                             &
                grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,                    &
                5+ntotal_)
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
!
   do 2000 lat = lat1,lat2
     lan=lat-latdon
#ifdef REDICE_GRID
     lonsd2=lonfdp(lan,mype)*2
#endif
!
     do k = 1,levs_
       spdlat(k,lan)=0.
     enddo
!
     do 140 k = 1,levs_
       do 140 j = 1,lonsd2
         ek(j,k)=(grs(j,ksu-1+k,lan)*grs(j,ksu-1+k,lan)                        &
                 +grs(j,ksv-1+k,lan)*grs(j,ksv-1+k,lan))*rbs2(lan)
         if (ek(j,k) .gt. spdlat(k,lan))  spdlat(k,lan)=ek(j,k)
140  continue
!
2000 continue
!
   do k = 1,levs_
     spdmax(k) = 0.0
     do lat = 1,lats2
       spdmax(k)=max(spdmax(k),spdlat(k,lat))
     enddo
     spdmax(k)=sqrt(spdmax(k))
   enddo
   call mpgetspd(spdmax)
   if( mype.eq.master ) then
     print 100,(spdmax(k),k=1,levs_)
100  format(' global checked speed maxima for all layers 2',                   &
          :/' spdmx(01:10)=',10f5.0,:/' spdmx(11:20)=',10f5.0,                 &
          :/' spdmx(21:30)=',10f5.0,:/' spdmx(31:40)=',10f5.0,                 &
          :/' spdmx(41:50)=',10f5.0,:/' spdmx(51:60)=',10f5.0,                 &
          :/' spdmx(61:70)=',10f5.0,:/' spdmx(71:80)=',10f5.0,                 &
          :/' spdmx(81:90)=',10f5.0,:/' spdmx(91:00)=',10f5.0)
   endif
!
   if( mype.eq.master ) then
     do k = 1,levs_
       if(spdmax(k).eq.0.) then
         print *,'run failure.  spdmax=0.'
         call abort
       endif
     enddo
   endif
!
!   End of UV check !!!
#endif
!
#ifndef GDAMP09
   call sph_wind2divor(uua,vva,xa,wa,anltop(1,1,1),                            &
                  anltop(1,1,levsp_+1),llstr,llens,lwvdef)
#endif
   call mpnk2nn(ya,lln22p_,levsp_, te,lnt22p_,levs_,1)
   call mpnk2nn(rta,lln22p_,levsp_, rq,lnt22p_,levs_,ntotal_)
   call mpnk2nn(xa,lln22p_,levsp_, di,lnt22p_,levs_,1)
!
! ... without divergence component ...
!
!      do j = 1,lnt22p_
!        do k = 1,levs_
!          di(j,k)=ditmp(j,K)
!        enddo
!      enddo
   call mpnk2nn(wa,lln22p_,levsp_, ze,lnt22p_,levs_,1)
   call mpn2nn(za,lln22p_,q,lnt22p_,1)
!
   fnxb=fnxx
!
#endif
   return
   end subroutine dyn_global_nudging 
!
!-------------------------------------------------------------------------------
   subroutine rmp_base_filename(iy,im,id,ih,fh,fnam,nchar)
!-------------------------------------------------------------------------------
!
!  returns base field file name with verifying date as
!  a suffix (sig.yyyymmddhh).
!
!-------------------------------------------------------------------------------
   implicit none
   integer            ::  iy,im,id,ih,jy,jm,jd,jh,nchar
   real               ::  fh,rjday
   character(len=4)   ::  cy
   character(len=2)   ::  cm,cd,ch
   character(len=20)  ::  fnam
!
   nchar=15
!
   call sfc_interp_date(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
!
   if(jy.lt.1000) then
     print *,'year less than 1000'
     call abort
   endif
!
   100 format('0',i1)
   101 format(i2)
   102 format(i4)
   write(cy,102) jy
   if(jm.lt.10) then
     write(cm,100) jm
   else
     write(cm,101) jm
   endif
   if(jd.lt.10) then
     write(cd,100) jd
   else
     write(cd,101) jd
   endif
   if(jh.lt.10) then
     write(ch,100) jh
   else
     write(ch,101) jh
   endif
!
   fnam='base.'//cy//cm//cd//ch
!
   return
   end subroutine rmp_base_filename
!
!-------------------------------------------------------------------------------
