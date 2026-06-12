#include <define.h>
   subroutine post_sgb_solver(fhour,idate,nsig,si,sl,                          &
                   nsgb,ncpus,ids,icen,icen2,igen,igases,iwater,keid)
!-------------------------------------------------------------------------------
!
! subprogram: post_sgb_solver           transforms sigma to sigma grib
!
! abstract: transforms a sigma spectral file to sigma grib1.
!   one latitude slice at a time, first sigma grid data
!   is transformed from sigma spectral coefficients.
!   the input data consists of vorticity, divergence, wind components,
!   temperature and specific humidity on the sigma surfaces as well as
!   surface pressure and orography and their horizontal gradients,
!   plus the potential of the bottom pressure gradient force.
!   relative humidity, vertical velocity and geopotential heights
!   are computed on the sigma surfaces and then interpolated to sigma
!   along with wind and temperature.  post_sgb_2d_diag fields are also computed.
!   the output data is quarterpacked and transposed to horizontal fields
!   which are then interpolated to the output grid and rounded
!   and packed into grib messages and written to the sigma grib file.
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
! usage:    call post_sgb_solver(fhour,idate,nsig,si,sl,
!    &                nsgb,ncpus,ids,pot,icen,icen2,igen)
!   input arguments:
!     fhour        real forecast hour
!     idate        integer (4) date
!     nsig         integer unit from which to read sigma file
!     si           real (levs+1) sigma interface values
!     sl           real (levs) sigma full level values
!     nsgb         integer unit to which to write grib messages
!     ncpus        integer number of cpus over which to distribute work
!     ids          integer (255) decimal scaling
!     icen         integer forecast center identifier
!     icen2        integer forecast subcenter identifier
!     igen         integer generating model identifier
!
! subprograms called:
!   post_read_coeff         read sigma coefficients
!   post_sgb_2d_param       set parameters for post_sgb_2d_diag fields
!   post_pole_vector        set defaults for field parameter identifiers
!   post_wave2grid          transform sigma coefficients
!   post_get_rh             compute relative humidity
!   post_get_omega          compute vertical velocity
!   post_get_height         compute geopotential heights
!   post_sgb_2d_diag        compute post_sgb_2d_diag fields
!   ptranw                  quarterpack and transpose data
!   ptranr                  unpack quarterpacked transposed data
!   post_pole_extrap        extrapolate pole values
!   post_trans_output_grid  separate hemispheric rows on grid
!   gribit                  create grib message
!   file_write_byte         write data by bytes
!
!  ::: structure ::: This file contains ... 
!
!    [post_sgb_solver] * ----- [post_sgb_2d_param] *
!                          |-- [post_sgb_2d_diag] *
!
!-------------------------------------------------------------------------------
   use module_parm_ptr
   use constant, only          : omega_,rd_,rv_,akapa_
   use paramodel, only         : ngases_,nwater_,ntotal_,                      &
                                 io=>lonf_,jo=>latg_,ko=>levs_,                &
                                 jcap=>jcap_,levs=>levs_,                      &
#ifdef DFS
                                 io2=>lonf_ ,io22=>lonf_ ,johf=>latg_
#else
                                 io2=>lonf2_,io22=>siodd_,johf=>latg2_,        &
                                 nctop => twoj1_
#endif
#ifdef DFS
   use dfsvar, only            : mtg,jlg
#endif
   use module_file_write, only : file_write_bin
#ifdef DCMIP
   use dcmip_grims, only       : icase,dcmip_grims_init,H,ptop,p0
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!
   integer, dimension(4)    ::  idate(4)
   integer, dimension(255)  ::  ids(255)
   integer                  ::  nsig,nsgb,ncpus,icen,icen2,igen,igases,iwater,mgrid
!
!  in this routine wave number is an input model resoultion
!
   integer, parameter       ::  mxbit=32
   real                     ::  fhour
#ifdef DCMIP
   integer, parameter  ::  nupa=22
#ifdef HYBRID
   real                ::  dum1(levs+1),dum2(levs+1)
#else
   real                ::  dum1(levs+1),dum2(levs)
#endif
#else
   integer, parameter  ::  nupa=18
#endif /* DCMIP end */
   integer, parameter  ::  nsun=9   !! 5 +  qc qr qi qs qg accumulation
   integer, parameter  ::  lenpds=28,lengds=32
!
   integer             ::  keid
   integer             ::  nc,ncho,jfhour,nfldp,nflds                 ,&
                           n,i,j,k,k1,k2,kan,lbm,ierr,is,lpds,lgds    ,&
                           kpmo,kpmr,kpmu,kpmv,kpmt,kpmz,kpma         ,&
                           kpmq,kpmqc,kpmqr,kpmqi,kpmqs,kpmqg         ,&
                           kpmnccn,kpmnc,kpmnr,kpmo3,kpmtke           ,&
#ifdef DCMIP
                           kpmq1,kpmq2,kpmq3,kpmq4                    ,&
                           ksq1,ksq2,ksq3,ksq4                        ,&
                           kpq1,kpq2,kpq3,kpq4                        ,&
#endif
                           kstke                                      ,&
                           ksnccn,ksnc,ksnr                           ,&
                           ksqg,kpnccn,kpnc,kpnr,kptke,ienst,iensi    ,&
                           lat,idir
#ifdef DCMIP
   real                ::  ps,p
#endif
   real                ::  rmaxsiin,rmaxsi,rk1,rkinv,dif
   real                ::  si(levs+1)
   real                ::  sl(levs)
   real                ::  sihyb(io22,levs+1),slhyb(io22,levs)   
   logical             ::  gauss_lat
   real                ::  colat1
#ifndef DFS
   integer             ::  ifax(100)
   real                ::  trig(io2)
   real, allocatable   ::  eps(:)
   real                ::  epstop(nctop/2)
   real                ::  clat(johf),slat(johf),wlat(johf)
   real, allocatable   ::  ss(:,:),sstop(:,:)
#else
   real                ::  clat(jo,3),slat(jo)
   real, allocatable   ::  fxys(:,:,:)
#endif
   real, allocatable   ::  fxs(:,:)
   real                ::  oxs(io22,levs),rxs(io22,levs)
   integer             ::  ipnum
   integer, allocatable, save  ::  ipo(:), npo(:)
   real                ::  zxs(io22,levs),zxi(io22,levs)
   real                ::  fxp(io22,nupa*ko+nsun)
   real                ::  fxy(io2,johf,nupa*ko+nsun)
!
#define DEFAULT
#ifdef DYNAMIC_ALLOC
#undef DEFAULT
   integer             ::  lgrib(ncpus)
   character grib(30+lenpds+lengds+io*jo*(mxbit+1)/8,ncpus)
#endif
#ifdef DEFAULT
   integer             ::  lgrib(1)
   character grib(30+lenpds+lengds+io*jo*(mxbit+1)/8,1)
#endif
   integer             ::  ipu(nupa*ko+nsun)
   integer, allocatable, save  ::  itl(:),il1(:),il2(:)
   integer             ::  mpf(0:255)
   integer, parameter                          ::  ipuu=33                    ,&
                                                   ipuv=34                    ,&
                                                   ipuo=39                    ,&
                                                   ipuz=7                     ,&
                                                   iput=11                    ,&
                                                   ipur=52                    ,&
                                                   ipus=175                   ,&
                                                   ipuq=51                    ,&
                                                   ipuqc=153                  ,&
                                                   ipuqr=152                  ,&
                                                   ipuqi=151                  ,&
                                                   ipuqs=150                  ,&
                                                   ipuqg=255                  ,&
                                                   ipunccn=181                ,&
                                                   ipunc=182                  ,&
                                                   ipunr=183                  ,&
                                                   ipuo3=154                  ,&
                                                   iputke=158                 ,&
#ifdef DCMIP
                                                   ipuq1=196                  ,&
                                                   ipuq2=197                  ,&
                                                   ipuq3=198                  ,&
                                                   ipuq4=199                  ,&
#endif
                                                   ipua=41

   integer             ::  ipusun(nsun)
   integer             ::  itlsun(nsun),il1sun(nsun),il2sun(nsun)
   integer             ::  ip1sun(nsun),ip2sun(nsun)
   integer             ::  kslp(2)
!
   integer             ::  iens(5)
!
   integer             ::  ipudef(nsun),itldef(nsun)
   data ipudef/001,054,007,001,76,77,67,68,69/
   data itldef/001,200,001,102,200,200,200,200,200/
   logical lpqc,lpqr,lpqi,lpqs,lpqg,lpo3,lptke
   logical lpnccn,lpnc,lpnr
#ifdef DCMIP
   logical lpq1,lpq2,lpq3,lpq4
#endif
!
   character(LEN=3), parameter  ::  fni='sgb'
   integer, parameter           ::  nchi=3
   character(LEN=80)            ::  fno
#ifdef ASSIGN
   character(LEN=80)            ::  asgnstr
#endif
   logical first, ipinfo
   data first/.true./, ipinfo/.true./
!-------------------------------------------------------------------------------
#ifdef DFS
   mgrid=0
#else
   mgrid=1
#endif
   if(mgrid.eq.1) then   ! gaussian
     idir=4
     gauss_lat=.true.
   else                  ! uniform
     idir=0
     gauss_lat=.false.
     !
     ! there is a problem in valid digit for lat-lon grid!
     ! simple remedy but not exact
     !
     colat1=90.-0.5*180./jo
   endif
   if(.not.allocated(ipo)) allocate( ipo(ko) )
   if(.not.allocated(npo)) allocate( npo(ko) )
   if(.not.allocated(itl)) allocate( itl(nupa*ko+nsun) )
   if(.not.allocated(il1)) allocate( il1(nupa*ko+nsun) )
   if(.not.allocated(il2)) allocate( il2(nupa*ko+nsun) )
#ifdef DFS
   nc=mtg*jlg
#else
   nc=(jcap+1)*(jcap+2)+1
#endif /* DFS end */
   nflds=5*levs+ntotal_*levs+6
!
#ifndef DFS
   allocate(eps(nc/2)  )
   allocate(ss(nc,nflds),sstop(nctop,nflds)  )
#else
   allocate(fxys(io,jo,nflds)  )
#endif
   allocate(fxs(io22,nflds)  )
!
   if (first) then
     first=.false.
     call ini_parm_ptr(levs,ko)
   endif
!
   call file_name(fni,nchi,fhour,fno,ncho)
#ifdef ASSIGN
   write(asgnstr,'(23hassign -s unblocked  u:,I2,)') nsgb
   call assign('assign -R')
   call assign(asgnstr)
#endif
   open(unit=nsgb,file=fno(1:ncho),form='unformatted',err=900)
   go to 901
900 continue
   write(6,*) ' error in opening file ',fno(1:ncho)
   call abort
901 continue
   write(6,*) ' file ',fno(1:ncho),' opened. unit=',nsgb
!
!  set some parameters
!
#ifdef DFS
   call post_read_coeff(nsig,keid,nflds,igases,iwater,nc,jcap,fhour,si,fxys   ,&
                        io,jo,io2,io22,johf)
#else
   call post_read_coeff(nsig,keid,nflds,igases,iwater,nc,nctop,jcap,fhour     ,&
             si,clat,slat,wlat,trig,ifax,eps,epstop,ss,sstop                  ,&
             io,jo,io2,io22,johf,gauss_lat)
   if(gauss_lat) colat1=acos(slat(1))
#endif
!
   lpqc=.false.
   lpqr=.false.
   lpqi=.false.
   lpqs=.false.
   lpqg=.false.
   lpnccn=.false.
   lpnc=.false.
   lpnr=.false.
   lpo3=.false.
   lptke=.false.
#ifdef DCMIP
   lpq1=.false.
   lpq2=.false.
   lpq3=.false.
   lpq4=.false.
#endif
   lpqc=iwater.ge.2
   lpqr=iwater.ge.3
   lpqi=iwater.ge.4
   lpqs=iwater.ge.5
   lpqg=iwater.ge.6
   if (iwater.eq.8) then
     lpqg=.false.
   endif
   if (iwater.ge.8) then
     lpnccn=.true.
     lpnc=.true.
     lpnr=.true.
   endif
!   lpnccn=iwater.ge.7
!   lpnc=iwater.ge.8
!   lpnr=iwater.ge.9
   lpo3=igases.ge.1
   lptke=igases.ge.2
#ifdef DCMIP
   lpq1=igases.ge.3
   lpq2=igases.ge.4
   lpq3=igases.ge.5
   lpq4=igases.ge.6
#endif
   print *,' return from post_read_coeff is ',nflds
   print *,' iwater igases from post_read_coeff = ',iwater,igases
   print *,' lpqc lpqr lpqi lpqs lpqg lpo3 = ',lpqc,lpqr,lpqi,lpqs,lpqg,lpo3
   print *,' lpnccn lpnc lpnr = ',lpnccn,lpnc,lpnr
   print *,' lpo3 lptke = ',lpo3, lptke
#ifdef DCMIP
   print *,' lpq1 lpq2 lpq3 lpq4 = ',lpq1,lpq2,lpq3,lpq4
#endif
!
   jfhour=nint(fhour)
   nfldp=nupa*ko+nsun
!
   kpmo=ko
   kpmr=ko
   kpmu=ko
   kpmv=ko
   kpmt=ko
   kpmz=ko
   kpma=ko
   kpmq=ko
   kpmqc=ko
   kpmqr=ko
   kpmqi=ko
   kpmqs=ko
   kpmqg=ko
   kpmnccn=ko
   kpmnc=ko
   kpmnr=ko
   kpmo3=ko
   kpmtke=ko
#ifdef DCMIP
   kpmq1=ko
   kpmq2=ko
   kpmq3=ko
   kpmq4=ko
#endif
!
   call post_sgb_2d_param(ko,sl,lpqc,lpqr,lpqi,lpqs,lpqg,nsun,                 &
               ipusun,itlsun,ip1sun,ip2sun,kslp)
!
   do n = 1,nsun
     ipusun(n)=ipudef(n)
     itlsun(n)=itldef(n)
   enddo
! 
!  set both input and output indices
! 
   kszs=1
   kszsx=2
   kszsy=3
   ksps=4
   kspsx=5
   kspsy=6
   ksd=7
   ksz=levs+7
   ksu=2*levs+7
   ksv=3*levs+7
   kst=4*levs+7
!
! -- wsm 1
!
   ksq=5*levs+7
!
! -- wsm 2
!
   if(lpqc) then
     ksqc = ksq + levs
   else
     ksqc = ksq
   endif
!
! --  wsm 3
!
   if(lpqr) then
     ksqr = ksqc + levs
   else
     ksqr = ksqc
   endif
!
! --  wsm 5
!
   if(lpqi) then
     ksqi = ksqr + levs
     ksqs = ksqi + levs
   else
     ksqi = ksqr
     ksqs = ksqi
   endif
!
! --  wsm 6
!
   if(lpqg) then
     ksqg = ksqs + levs
   else
     ksqg = ksqs
   endif
!
! --  wdm 5 and 6
!
   if(lpnccn) then
     ksnccn = ksqg + levs
     ksnc = ksnccn + levs
     ksnr = ksnc + levs
   else
     ksnccn = ksqg
     ksnc = ksqg
     ksnr = ksqg
   endif
!
! --  o3
!
   if(lpo3) then
     kso3 = ksnr + levs
   else
     kso3 = ksnr
   endif
!
! --  tke
!
   if(lptke) then
     kstke = kso3 + levs
   else
     kstke = kso3
   endif
#ifdef DCMIP
!
! - dcmip
!
   if(lpq1) then
     ksq1= kstke + levs
   else
     ksq1= kstke
   endif
   if(lpq2) then
     ksq2= ksq1 + levs
   else
     ksq2= ksq1
   endif
   if(lpq3) then
     ksq3= ksq2 + levs
   else
     ksq3= ksq2
   endif
   if(lpq4) then
     ksq4= ksq3 + levs
   else
     ksq4= ksq3
   endif
#endif
#ifdef DBG
   print*,' ksqc ksqr ksqi ksqs ksqg ksnccn ksnc ksnr kso3 '
   print*, ksqc,ksqr,ksqi,ksqs,ksqg,ksnccn,ksnc,ksnr,kso3,kstke
#ifdef DCMIP
   print*,' ksq1 ksq2 ksq3 ksq4'
   print*, ksq1,ksq2,ksq3,ksq4
#endif
#endif
   kpz=1
   kpu=ko+1
   kpv=2*ko+1
   kpr=3*ko+1
   kpq=4*ko+1
   kpt=5*ko+1
   kpo=6*ko+1
   kpa=7*ko+1
   kpqc=8*ko+1
   kpqr=9*ko+1
   kpqi=10*ko+1
   kpqs=11*ko+1
   kpqg=12*ko+1
   kpnccn=13*ko+1
   kpnc=14*ko+1
   kpnr=15*ko+1
   kpo3=16*ko+1
   kptke=17*ko+1
#ifndef DCMIP
   kpsun=18*ko+1
#else
   kpq1=18*ko+1
   kpq2=19*ko+1
   kpq3=20*ko+1
   kpq4=21*ko+1
   kpsun=22*ko+1
#endif
!
!  set some grib parameters
!
   do i = 1,nupa*ko
     ipu(i)=1
   enddo
   do i = kpu,kpu+kpmu-1
     ipu(i)=ipuu
   enddo
   do i = kpv,kpv+kpmv-1
     ipu(i)=ipuv
   enddo
   do i = kpo,kpo+kpmo-1
     ipu(i)=ipuo
   enddo
   do i = kpz,kpz+kpmz-1
     ipu(i)=ipuz
   enddo
   do i = kpt,kpt+kpmt-1
     ipu(i)=iput
   enddo
   do i = kpa,kpa+kpma-1
     ipu(i)=ipua
   enddo
   do i = kpr,kpr+kpmr-1
     ipu(i)=ipur
   enddo
   do i = kpq,kpq+kpmq-1
     ipu(i)=ipuq
   enddo
   do i = kpqc,kpqc+kpmqc-1
     if(lpqc) then
       ipu(i)=ipuqc
     else
       ipu(i)=0
     endif
   enddo
   do i = kpqr,kpqr+kpmqr-1
     if(lpqr) then
       ipu(i)=ipuqr
     else
       ipu(i)=0
     endif
   enddo
   do i = kpqi,kpqi+kpmqi-1
     if(lpqi) then
       ipu(i)=ipuqi
     else
       ipu(i)=0
     endif
   enddo
   do i = kpqs,kpqs+kpmqs-1
     if(lpqs) then
       ipu(i)=ipuqs
     else
       ipu(i)=0
     endif
   enddo
   do i = kpqg,kpqg+kpmqg-1
     if(lpqg) then
       ipu(i)=ipuqg
     else
       ipu(i)=0
     endif
   enddo
   do i = kpnccn,kpnccn+kpmnccn-1
     if(lpnccn) then
       ipu(i)=ipunccn
     else
       ipu(i)=0
     endif
   enddo
   do i = kpnc,kpnc+kpmnc-1
     if(lpnc) then
       ipu(i)=ipunc
     else
       ipu(i)=0
     endif
   enddo
   do i = kpnr,kpnr+kpmnr-1
     if(lpnr) then
       ipu(i)=ipunr
     else
       ipu(i)=0
     endif
   enddo
   do i = kpo3,kpo3+kpmo3-1
     ipu(i)=ipuo3
   enddo
   do i = kptke,kptke+kpmtke-1
     ipu(i)=iputke
   enddo
!
#ifdef DCMIP
   do i = kpq1,kpq1+kpmq1-1
     ipu(i)=ipuq1
   enddo
   do i = kpq2,kpq2+kpmq2-1
     ipu(i)=ipuq2
   enddo
   do i = kpq3,kpq3+kpmq3-1
     ipu(i)=ipuq3
   enddo
   do i = kpq4,kpq4+kpmq4-1
     ipu(i)=ipuq4
   enddo
#endif
   do i = kpsun,kpsun+nsun-1
     ipu(i)=ipusun(i-kpsun+1)
   enddo
!
   ienst=0
   iensi=0
   iens(1)=1
   iens(2)=ienst
   iens(3)=iensi
   iens(4)=1
   iens(5)=255
!
   call post_pole_vector(2,mpf)
!
#ifdef DCMIP
   call dcmip_grims_init(dum1,dum2)    ! for icase,ptop,H
   ids(007)=2   ! geopotential height (m)
   ids(158)=6   ! tracer 1
   ids(196)=6   ! tracer 1
   ids(197)=6   ! tracer 2
   ids(198)=6   ! tracer 3
   ids(199)=6   ! tracer 4
#endif /* DCMIP end */
#ifdef SMP
   loop_lat: do lat=1,johf
#else
#ifdef DFS
   loop_lat: do lat=1,jo
#else
   loop_lat: do lat=2,johf
#endif /* DFS end */
#endif /* SMP end */
     do n = 1,nfldp
       do i = 1,io22
         fxp(i,n)=0.0
       enddo
     enddo
#ifndef DFS
#ifdef SMP
     call post_wave2grid(nflds,igases,iwater,ss,fxs,                           &
                       io2,io22,johf)
#else
     call post_wave2grid(nflds,igases,iwater,trig,ifax,eps,epstop,             &
                       ss,sstop, clat(lat),slat(lat),fxs(1,1),                 &
                       io2,io22,johf)
#endif
#else /* DFS */
     do n = 1,nflds
       do i = 1,io22
         fxs(i,n)=fxys(i,lat,n)
       enddo
     enddo
#endif /* ~DFS end */
#ifdef DBG
     do k = 1,levs
       print 1960, k,fxs(1,kszs), fxs(1,ksps),                                 &
                  fxs(1,ksd+k-1), fxs(1,ksz+k-1), fxs(1,ksu+k-1),              &
                  fxs(1,kst+k-1), fxs(1,ksq+k-1)
1960   format(i5,7e13.5)
     enddo
#endif
!
     rmaxsiin=si(1)
     do k = 1,levs+1
       rmaxsi=max(rmaxsi,si(k))
     enddo
!
     if (rmaxsi.gt.1. .or. rmaxsi.eq.0.) then
!
       rk1 = akapa_ + 1.
       rkinv=1./akapa_
!
       do i = 1,io2
         do k = 1,levs
           sihyb(i,k)=si(k)/fxs(i,ksps)/1000.+sl(k)
         enddo
         sihyb(i,levs+1)=0.
         do k = 1,levs
           dif = sihyb(i,k)**rk1 - sihyb(i,k+1)**rk1
           dif = dif / (rk1*(sihyb(i,k)-sihyb(i,k+1)))
           slhyb(i,k) = dif**rkinv
         enddo
       enddo
       ipnum=109    ! hybrid
     else
       do i = 1,io2
         do k = 1,levs+1
           sihyb(i,k)=si(k)
         enddo
         do k = 1,levs
           slhyb(i,k)=sl(k)
         enddo
       enddo
       ipnum=107    ! sigma
     endif
#ifdef DCMIP
     if(icase.eq.11 .or. icase.eq.12) ipnum=105      ! height
#endif
!
! get rh, omega, height
!
     call post_get_rh(io2,io22,levs,slhyb,                                     &
              fxs(1,ksps    ),fxs(1,ksq    ),fxs(1,kst    ),                   &
              rxs)
!
     call post_get_omega(io,io2,io22,levs,sihyb,slhyb,                         &
#ifdef SMP
              fxs(1,ksu    ),                                                  &
#else
              fxs(1,ksps    ),fxs(1,kspsx    ),fxs(1,kspsy    ),               &
              fxs(1,ksd    ),fxs(1,ksu    ),fxs(1,ksv    ),                    &
#endif
              oxs)
!
     call post_get_height(io2,io22,levs,sihyb,slhyb,                           &
              fxs(1,kszs    ),fxs(1,kst    ),fxs(1,ksq    ),                   &
              zxs,zxi)
#ifdef DCMIP
!
! DCMIP pressure,height
!
     if(icase.eq.11 .or. icase.eq.12) then
       do k = 1,ko
         do i = 1,io2
           ps       = fxs(i,ksps)*1.e3
#ifdef HYBRID
           p        = 0.5*(si(k)+si(k+1))+0.5*(sl(k)+sl(k+1))*ps
           zxs(i,k) = H*log(p0/p)+1.
#else
           p        = ptop+slhyb(i,k)*(ps-ptop)
           zxs(i,k) = H*log(p0/p)
#endif
         enddo
       enddo
     endif
#endif /* DCMIP end */
!
     if(ipinfo) then
       ipinfo=.false.
!
       do k = 1,ko
         ipo(k)=ipnum
         npo(k)=nint(slhyb(1,k)*1.e4)
       enddo
#ifdef DCMIP
       if(icase.eq.11 .or. icase.eq.12) then
         do k = 1,ko
           npo(k)=nint(zxs(1,k))
         enddo
       endif
#ifdef HYBRID
       if(icase.eq.13 .or. icase.eq.200 .or. icase.eq.410 .or. icase.eq.42 ) then
         do k = 1,ko
           npo(k)=nint( ( 0.5*(si(k)+si(k+1))/p0+0.5*(sl(k)+sl(k+1)) )*10000. )
         enddo
       endif
#endif
#endif /* DCMIP end */
!
! itl
!
       do i = kpu,kpu+ko-1
         itl(i)=ipo(i-kpu+1)
       enddo
       do i = kpv,kpv+ko-1
         itl(i)=ipo(i-kpv+1)
       enddo
       do i = kpo,kpo+ko-1
         itl(i)=ipo(i-kpo+1)
       enddo
       do i = kpz,kpz+ko-1
         itl(i)=ipo(i-kpz+1)
       enddo
       do i = kpt,kpt+ko-1
         itl(i)=ipo(i-kpt+1)
       enddo
       do i = kpa,kpa+ko-1
         itl(i)=ipo(i-kpa+1)
       enddo
       do i = kpr,kpr+ko-1
         itl(i)=ipo(i-kpr+1)
       enddo
       do i = kpq,kpq+ko-1
         itl(i)=ipo(i-kpq+1)
       enddo
       do i = kpo3,kpo3+ko-1
         itl(i)=ipo(i-kpo3+1)
       enddo
       do i = kptke,kptke+ko-1
         itl(i)=ipo(i-kptke+1)
       enddo
#ifdef DCMIP
       do i = kpq1,kpq1+ko-1
         itl(i)=ipo(i-kpq1+1)
       enddo
       do i = kpq2,kpq2+ko-1
         itl(i)=ipo(i-kpq2+1)
       enddo
       do i = kpq3,kpq3+ko-1
         itl(i)=ipo(i-kpq3+1)
       enddo
       do i = kpq4,kpq4+ko-1
         itl(i)=ipo(i-kpq4+1)
       enddo
#endif
       do i = kpqc,kpqc+ko-1
         itl(i)=ipo(i-kpqc+1)
       enddo
       do i = kpqr,kpqr+ko-1
         itl(i)=ipo(i-kpqr+1)
       enddo
       do i = kpqi,kpqi+ko-1
         itl(i)=ipo(i-kpqi+1)
       enddo
       do i = kpqs,kpqs+ko-1
         itl(i)=ipo(i-kpqs+1)
       enddo
       do i = kpqg,kpqg+ko-1
         itl(i)=ipo(i-kpqg+1)
       enddo
       do i = kpsun,kpsun+nsun-1
         itl(i)=itlsun(i-kpsun+1)
       enddo
!
! il1
!
       do i = 1,nupa*ko+nsun
         il1(i)=0
       enddo
       do i = kpsun,kpsun+nsun-1
!        il1(i)=il1sun(i-kpsun+1)
         il1(i)=0
       enddo
!
       do i = kpu,kpu+ko-1
         il2(i)=npo(i-kpu+1)
       enddo
       do i = kpv,kpv+ko-1
         il2(i)=npo(i-kpv+1)
       enddo
       do i = kpo,kpo+ko-1
         il2(i)=npo(i-kpo+1)
       enddo
       do i = kpz,kpz+ko-1
         il2(i)=npo(i-kpz+1)
       enddo
       do i = kpt,kpt+ko-1
         il2(i)=npo(i-kpt+1)
       enddo
       do i = kpa,kpa+ko-1
         il2(i)=npo(i-kpa+1)
       enddo
       do i = kpr,kpr+ko-1
         il2(i)=npo(i-kpr+1)
       enddo
       do i = kpq,kpq+ko-1
         il2(i)=npo(i-kpq+1)
       enddo
       do i = kpqc,kpqc+ko-1
         il2(i)=npo(i-kpqc+1)
       enddo
       do i = kpqr,kpqr+ko-1
         il2(i)=npo(i-kpqr+1)
       enddo
       do i = kpqi,kpqi+ko-1
         il2(i)=npo(i-kpqi+1)
       enddo
       do i = kpqs,kpqs+ko-1
         il2(i)=npo(i-kpqs+1)
       enddo
       do i = kpqg,kpqg+ko-1
         il2(i)=npo(i-kpqg+1)
       enddo
       do i = kpo3,kpo3+ko-1
         il2(i)=npo(i-kpo3+1)
       enddo
       do i = kptke,kptke+ko-1
         il2(i)=npo(i-kptke+1)
       enddo
#ifdef DCMIP
       do i = kpq1,kpq1+ko-1
         il2(i)=npo(i-kpq1+1)
       enddo
       do i = kpq2,kpq2+ko-1
         il2(i)=npo(i-kpq2+1)
       enddo
       do i = kpq3,kpq3+ko-1
         il2(i)=npo(i-kpq3+1)
       enddo
       do i = kpq4,kpq4+ko-1
         il2(i)=npo(i-kpq4+1)
       enddo
#endif

       do i = kpsun,kpsun+nsun-1
         il2(i)=0
       enddo
     endif
!
! fxp
!
     do k = 1,ko
#ifdef SMP
       fxp(:,kpu+k-1)=fxs(:,ksd+k-1)
       fxp(:,kpv+k-1)=fxs(:,ksz+k-1)
#else
       fxp(:,kpu+k-1)=fxs(:,ksu+k-1)
       fxp(:,kpv+k-1)=fxs(:,ksv+k-1)
#endif
       fxp(:,kpa+k-1)=fxs(:,ksz+k-1)
       fxp(:,kpo+k-1)=oxs(:,k)
       fxp(:,kpz+k-1)=zxs(:,k)
       fxp(:,kpt+k-1)=fxs(:,kst+k-1)
       fxp(:,kpr+k-1)=rxs(:,k)
       fxp(:,kpq+k-1)=fxs(:,ksq+k-1)
       if(lpqc) fxp(:,kpqc+k-1)=fxs(:,ksqc+k-1)
       if(lpqr) fxp(:,kpqr+k-1)=fxs(:,ksqr+k-1)
       if(lpqi) fxp(:,kpqi+k-1)=fxs(:,ksqi+k-1)
       if(lpqs) fxp(:,kpqs+k-1)=fxs(:,ksqs+k-1)
       if(lpqg) fxp(:,kpqg+k-1)=fxs(:,ksqg+k-1)
       if(lpo3) fxp(:,kpo3+k-1)=fxs(:,kso3+k-1)
       if(lptke) fxp(:,kptke+k-1)=fxs(:,kstke+k-1)
#ifdef DCMIP
       if(lpq1) fxp(:,kpq1+k-1)=fxs(:,ksq1+k-1)
       if(lpq2) fxp(:,kpq2+k-1)=fxs(:,ksq2+k-1)
       if(lpq3) fxp(:,kpq3+k-1)=fxs(:,ksq3+k-1)
       if(lpq4) fxp(:,kpq4+k-1)=fxs(:,ksq4+k-1)
#endif
     enddo
#ifdef DBG
     do k = 1, levs
       print 1964, k,fxp(1,kpu+k-1),fxp(1,kpv+k-1),fxp(1,kpo+k-1),             &
                  fxp(1,kpz+k-1), fxp(1,kpt+k-1), fxp(1,kpr+k-1),              &
                  fxp(1,kpq+k-1), fxs(1,ksps)
1964   format('SGB ',i5,8e13.5)
     enddo
#endif
!
     call post_sgb_2d_diag(io2,io22,levs,nsun,kslp,sihyb,                      &
               fxs(1,kszs    ),fxs(1,ksps    ),                                &
               fxs(1,ksq     ),fxp(1,kpz     ),                                &
               fxp(1,kpsun    ),                                               &
               lpqc,lpqr,lpqi,lpqs,lpqg,                                       &
               fxs(1,ksqc),fxs(1,ksqr),                                        &
               fxs(1,ksqi),fxs(1,ksqs),fxs(1,ksqg))
!
     do n = 1,nfldp
       do i = 1,io2
         fxy(i,lat,n)=fxp(i,n)
       enddo
     enddo
!
   enddo loop_lat
!
#ifndef SMP
!
!  fill pole points with extrapolation
!
   do k = 1,nfldp
#ifdef DFS
     call post_pole_extrap(mpf(ipu(k)),io,fxy(1,2,k),fxy(1,jo-1,k),            &
                    fxy(1,1,k),fxy(1,jo,k))
#else
     call post_pole_extrap(mpf(ipu(k)),io,fxy(1,2,k),fxy(1+io,2,k),            &
                      fxy(1,1,k),fxy(1+io,1,k))
#endif
   enddo
#endif /* ~ SMP end */
!
   call print_maxmin_six(fxy(1,1,kpa),io2*johf,ko,1,ko,'REL vor')
#ifdef DFS
   call post_equall_lat(jo,slat,clat,0)
#endif
   do k = 0,ko-1
     do j = 1,johf
       do i = 1,io2
#ifdef DCMIP
         if(icase.eq.200) then
            fxy(i,j,kpa+k)=fxy(i,j,kpa+k)
         else
            fxy(i,j,kpa+k)=fxy(i,j,kpa+k)+2.*omega_*slat(j)
         endif
#else
            fxy(i,j,kpa+k)=fxy(i,j,kpa+k)+2.*omega_*slat(j)
#endif
       enddo
     enddo
   enddo
   call print_maxmin_six(fxy(1,1,kpa),io2*johf,ko,1,ko,'ABS vor')
!
!  loop over groups of horizontal fields
!
   do k1 = 1,nfldp,ncpus
     k2=min(k1+ncpus-1,nfldp)
!
!  and round to the number of bits and engrib the field in parallel
!
     do k = k1,k2
       kan=k-k1+1
       lgrib(kan)=0
       if(ipu(k).gt.0) then
#ifndef DFS
         call post_trans_output_grid(io,jo,fxy(1,1,k),1)
#endif
!
!        if(k.eq.1) then
!          call post_quick_print(fxy(1,1,k),io,jo,1./100.)
!        endif
!        call file_write_bin(222,fxy(1,1,k),io,jo,1,0)
         call file_make_grib(fxy(1,1,k),lbm,idir,io,jo,mxbit,colat1,           &
                        lenpds,2,icen,igen,0,                                  &
                        ipu(k),itl(k),il1(k),il2(k),                           &
                        idate(4),idate(2),idate(3),idate(1),                   &
                        1,jfhour,0,10,0,0,icen2,ids(ipu(k)),iens,              &
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,                         &
                        grib(1,kan),lgrib(kan),ierr)
!
       endif
     enddo
!
!  write out grib messages sequentially
!
     do k = k1,k2
       kan=k-k1+1
       if((.not.lpqc).and.k.ge.kpqc.and.k.lt.kpqr) lgrib(kan) = 0
       if((.not.lpqr).and.k.ge.kpqr.and.k.lt.kpqi) lgrib(kan) = 0
       if((.not.lpqi).and.k.ge.kpqi.and.k.lt.kpqs) lgrib(kan) = 0
       if((.not.lpqs).and.k.ge.kpqs.and.k.lt.kpqg) lgrib(kan) = 0
       if((.not.lpqg).and.k.ge.kpqg.and.k.lt.kpo3) lgrib(kan) = 0
       if((.not.lpo3).and.k.ge.kpo3.and.k.lt.kptke) lgrib(kan) = 0
#ifndef DCMIP
       if((.not.lptke).and.k.ge.kptke.and.k.lt.kpsun) lgrib(kan) = 0
#else
       if((.not.lptke).and.k.ge.kptke.and.k.lt.kpq1) lgrib(kan) = 0
       if((.not.lpq1).and.k.ge.kpq1.and.k.lt.kpq2) lgrib(kan) = 0
       if((.not.lpq2).and.k.ge.kpq2.and.k.lt.kpq3) lgrib(kan) = 0
       if((.not.lpq3).and.k.ge.kpq3.and.k.lt.kpq4) lgrib(kan) = 0
       if((.not.lpq4).and.k.ge.kpq4.and.k.lt.kpsun) lgrib(kan) = 0
#endif
       if(lgrib(kan).gt.0) then
         is=8
         lpds=ichar(grib(is+1,kan))*65536                                      &
              +ichar(grib(is+2,kan))*256                                       &
              +ichar(grib(is+3,kan))
         is=8+lpds
         lgds=ichar(grib(is+1,kan))*65536                                      &
              +ichar(grib(is+2,kan))*256                                       &
              +ichar(grib(is+3,kan))
         call file_write_byte(nsgb,lgrib(kan),grib(1,kan))
         print *,' grib1 written to ',nsgb,' of length ',lgrib(kan)
       endif
     enddo
   enddo
! 
#ifndef DFS
   deallocate(eps )
   deallocate(ss,sstop  )
#else
   deallocate(fxys )
#endif
   deallocate(fxs )
!
   return
   end subroutine post_sgb_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_sgb_2d_param(ko,po,lpqc,lpqr,lpqi,lpqs,lpqg,nsun,           &
                     ipusun,itlsun,ip1sun,ip2sun,kslp)
!-------------------------------------------------------------------------------
!
! subprogram:  post_sgb_2d_param     set parameters for post_sgb_2d_diag fields
!
! abstract: sets parameters for the post_sgb_2d_diag fields.
!   parameters returned are parameter indicator, level type indicator,
!   two level numbers and decimal scaling all required for the pds
!   section of the grib1 message as well relevant sigma layer numbers
!   for the three lower level relative humidity fields.
!   the current nsun=22 post_sgb_2d_diag fields are:
!     1) surface pressure
!     2) precipitable water
!     3) surface orography
!     4) sea level pressure
!
! program history log:
!   92-10-31  iredell
!
! usage:    call post_sgb_2d_param(ko,po,lpqc,lpqr,lpqi,lpqs,lpqg,nsun,
!    &                  ipusun,itlsun,ip1sun,ip2sun,kslp)
!
!   input argument list:
!     km       - integer number of levels
!     ko       - integer number of pressure levels
!     po       - real (ko) sigma
!
!   output argument list:
!     ipusun   - integer (nsun) parameter indicators
!     itlsun   - integer (nsun) level type indicators
!     ip1sun   - integer (nsun) first level numbers
!     ip2sun   - integer (nsun) second level numbers
!     kslp     - integer (2) relevant pressure levels for slp
!
! subprograms called:
!   isrchflex - find first value in an array le target value
!   isrcheqx  - find first value in an array equal to target value
!
!-------------------------------------------------------------------------------
   real     ::  po(ko)
   real     ::  pslp(2)
   integer  ::  ipusun(nsun)
   integer  ::  itlsun(nsun),ip1sun(nsun),ip2sun(nsun)
   integer  ::  kslp(2)
   data pslp/1.,0.5/
!-------------------------------------------------------------------------------
   kslp(1)=mod(isrcheqx(ko,po,1,pslp(1)),ko+1)
   kslp(2)=mod(isrcheqx(ko,po,1,pslp(2)),ko+1)
   do n = 1,nsun
     ip1sun(n)=0
     ip2sun(n)=0
   enddo
! 
   return
   end subroutine post_sgb_2d_param
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_sgb_2d_diag(im,ix,km,nsun,kslp,si,zs,ps,q,zm,sun,           &
                     lpqc,lpqr,lpqi,lpqs,lpqg,qc,qr,qi,qs,qg)
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_
!-------------------------------------------------------------------------------
!
! subprogram:    post_sgb_2d_diag      compute post_sgb_2d_diag fields
!
! abstract: computes post_sgb_2d_diag fields from sigma level winds, omega,
!   temperature, moisture, and relative humidity.
!   relative humidity is perversely averaged across sigma layers.
!   the current nsun=4 post_sgb_2d_diag fields are:
!     1) surface pressure
!     2) precipitable water
!     3) surface orography
!     4) sea level pressure
!
! subprograms called:
!   post_sigma2tropopause    interpolate sigma to tropopause level
!   post_max_wind            interpolate sigma to maxwind level
!   post_lifting_index       compute best lifted index
!
! program history log:
!   92-10-31  mccalla,iredell
!   96-01-18  kanamistu  greatly simplified
!
! usage:    call post_sgb_2d_diag(im,ix,km,kslp,rnhour,si,
!    &                  zs,ps,q,zm,sun,
!    +                  lpqc,lpqr,lpqi,lpqs,qc,qr,qi,qs,qg)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of levels
!     kslp     - integer (2) relevant pressure levels for slp
!     rnhour   - real hours rain accumulated
!     si       - real (km) sigma interfaces
!     zs       - real (im) surface orography in m
!     ps       - real (im) surface pressure in kpa
!     q        - real (ix,km) specific humidity in kg/kg
!     zm       - real (ix,*) height on sigma surface in m
!
!   output argument list:
!     sun      - real (ix,nsun) post_sgb_2d_diag fields given above
!
!-------------------------------------------------------------------------------
   integer          ::  kslp(2)
   real             ::  si(ix,km+1)
   real             ::  zs(im),ps(im)
!
   real             ::  q(ix,km)
   real             ::  zm(ix,km)
   real             ::  sun(ix,nsun)
   real             ::  qc(ix,km), qr(ix,km), qi(ix,km), qs(ix,km), qg(ix,km)
   real             ::  sk(2*lonf_),dum(2*lonf_)
   real, parameter  ::  g= 9.8000e+0 ,cp= 1.0046e+3 
   real, parameter  ::  rd= 2.8705e+2 ,rv= 4.6150e+2 
   real, parameter  ::  pm1=1.e5,tm1=287.45,zm1=113.,zm2=5572.
   real, parameter  ::  rk=rd/cp,fslp=g*(zm2-zm1)/(rd*tm1)
!
   logical lpqc,lpqr,lpqi,lpqs,lpqg
!-------------------------------------------------------------------------------
!
!  surface pressure
!
   do i = 1,im
     sun(i,1)=ps(i)*1.e3
   enddo
#ifdef DBG
   call print_maxmin_six(sun(1,1),im,1,1,1,'ps')
#endif
!
!  precipitable water
!
   do i = 1,im
     sun(i,2)=0.
   enddo
   do k = 1,km
     do i = 1,im
       ds=si(i,k)-si(i,k+1)
       sun(i,2)=sun(i,2)+q(i,k)*ds
     enddo
   enddo
   do i = 1,im
     sun(i,2)=sun(i,2)*ps(i)*1.e3/g
   enddo
#ifdef DBG
   call print_maxmin_six(sun(1,2),im,1,1,1,'precipitable water')
#endif
! 
!  surface orography
!
   do i = 1,im
     sun(i,3)=zs(i)
   enddo
#ifdef DBG
   call print_maxmin_six(sun(1,3),im,1,1,1,'surface orography')
#endif
! 
!  sea level pressure
!
   if(kslp(1).gt.0.and.kslp(2).gt.0) then
     k1=kslp(1)
     k2=kslp(2)
     do i = 1,im
       sun(i,4)=pm1*exp(fslp*zm(i,k1)/(zm(i,k2)-zm(i,k1)))
     enddo
   else
     do i = 1,im
       sun(i,4)=0.
     enddo
   endif
#ifdef DBG
   call print_maxmin_six(sun(1,4),im,1,1,1,'sea level pressure')
#endif
! 
!  accumulated qc
!
   if (lpqc) then
     do i = 1,im
       sun(i,5)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,5)=sun(i,5)+qc(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,5)=sun(i,5)*ps(i)*1.e3/g
     enddo
#ifdef DBG
     call print_maxmin_six(sun(1,5),im,1,1,1,'qc')
#endif
   endif
!
!  accumulated qr
!
   if (lpqr) then
     do i = 1,im
       sun(i,6)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,6)=sun(i,6)+qr(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,6)=sun(i,6)*ps(i)*1.e3/g
     enddo
#ifdef DBG
     call print_maxmin_six(sun(1,6),im,1,1,1,'qr')
#endif
   endif
! 
!  accumulated qi
!
   if (lpqi) then
     do i = 1,im
       sun(i,7)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,7)=sun(i,7)+qi(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,7)=sun(i,7)*ps(i)*1.e3/g
     enddo
#ifdef DBG
    call print_maxmin_six(sun(1,7),im,1,1,1,'qi')
#endif
   endif
!
!  accumulated qs
!
   if (lpqs) then
     do i = 1,im
       sun(i,8)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,8)=sun(i,8)+qs(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,8)=sun(i,8)*ps(i)*1.e3/g
     enddo
#ifdef DBG
     call print_maxmin_six(sun(1,8),im,1,1,1,'qs')
#endif
   endif
!
!  accumulated qg
!
   if (lpqg) then
     do i = 1,im
       sun(i,9)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,9)=sun(i,9)+qg(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,9)=sun(i,9)*ps(i)*1.e3/g
     enddo
#ifdef DBG
      call print_maxmin_six(sun(1,9),im,1,1,1,'qg')
#endif
   endif
!
   return
   end subroutine post_sgb_2d_diag
!-------------------------------------------------------------------------------
