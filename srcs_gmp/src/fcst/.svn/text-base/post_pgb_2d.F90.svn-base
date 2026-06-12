#include <define.h>
#ifndef RMP
#undef STDAMP_P
#endif
#ifndef HYBRID
   subroutine post_pgb_2d(fhour,idate,sl,si,                                   &
#else
   subroutine post_pgb_2d(fhour,idate,ak5,bk5,                                 &
#endif
                     fgz,fq,fte,fuu,fvv,frq,flat,flon,n1                       &
#ifdef RMP
                    ,lbaseout,fm2)
#else
                     )
#endif
!-------------------------------------------------------------------------------
#ifdef MRG_POST
!
! subprogram: post_pgb_2d        transforms a sigma to pressure grib
!
! abstract: transforms a sigma spectral file to pressure grib1
!
! usage:    call post_pgb_2d(fhour,idate,sl,si,
!    &                  fgz,fq,fte,fuu,fvv,frq,flat,flon,n1)
!   input arguments:
!     fhour        real forecast hour
!     idate        integer (4) date
!     si           real (levr+1) sigma interface values
!     sl           real (levr) sigma full level values
!     fgz
!     fq
!     fte
!     fuu
!     fvv
!     frq
!     flat
!     flon
!     n1
!
! subprograms called:
!   idsdef       set ids
!   post_pgb_2d_param       set parameters for post_pgb_2d_diag fields
!   rpost_wave2grid         convert coef to grid
!   post_get_rh             compute relative humidity
!   post_get_omega          compute vertical velocity
!   post_get_height         compute geopotential heights
!   post_sigma2pressure     interpolate sigma to pressure
!   post_pgb_2d_diag        compute post_pgb_2d_diag fields
!   ptranw                  quarterpack and transpose data
!   ptranr                  unpack quarterpacked transposed data
!   gribit                  create grib message
!   file_write_byte         write data by bytes
!
!  define MPIGRIB will gather all the arrays to master processor then
!  convert to grib and writes on master processor.  this is most inefficient.
!
!  define PGBGATHER will conver to grib on each processor, then gather
!  gribbed arrays to master processor and then master processor writes.
!  this is second efficient, but i/o is serial.
!
!  undef MPIGRIB, undef PGBGATHER will do the coversion to grib on
!  each processor and then each processor writes.  this is most efficient
!  since computation and i/o are both parallel.
!
!  when MPIGRIB is defined, PGBGATHER should not be defined.
!  when PGBGATHER is defined, PGBGRIB should not be defined.
!
!  ::: structure ::: This file contains ... 
!
!    [post_pgb_2d] *   
!         |
!         |-- [post_pgb_2d_param] *
!         |-- [post_pgb_2d_diag] *
!         |-- [post_lifting_index] *
!         |-- [post_max_wind] *
!                  |
!             (post_interp_spline)
!                  |-- [post_spline_derivative]
!                  |     : compute 2nd derivatives for cubic spline
!                  |-- [post_spline_max]
!                        : determine maximum value of cubic spline
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use constant, only  : g_,cp_,rd_,rv_
   use paramodel, only : levs_,npes_,ngases_,nwmass_,ko_,kt_
#ifdef RMP
   use paramodel, only : igrd1_,igrd12_
#ifdef MP
   use paramodel, only : igrd12p_,jgrd12p_
#else
   use paramodel, only : lnwav_
#endif
   use rscomloc
   use rscomltb
#else   /* not RMP */
   use paramodel, only : lonf_,latg_,lonf2_,latg2_,io_,jo_
#ifdef MP
   use paramodel, only : lnt22p_,lonf2p_,latg2p_
#else
   use paramodel, only : lnt22_
#endif
#ifdef REDUCE_GRID
   use comfgrid
   use comreduce
#endif
   use module_trans, only      :  rmp_trans2output_grid
#endif  /* RMP */
#ifdef MP
   use commpi
#endif
   use module_trans, only      :  dyn_trans2output_grid
#ifdef MP
#ifdef LINUX
#define PGBGATHER
#endif
#endif
#ifndef RMP
#undef STDAMP_P
#endif
!-------------------------------------------------------------------------------
   integer,parameter    ::  l2=2
   integer,parameter    ::  nupa=11,nupt=4,nsun=29,nupas=8
   integer,parameter    ::  lenpds=28,lengds=32
   real,parameter       ::  fvirt=rv_/rd_-1.
   real,parameter       ::  gor2=g_/rd_*l2,eps=rd_/rv_,epsm1=rd_/rv_-1.
!
   real                 ::  si(levs_+1),sl(levs_),po(ko_)
   integer              ::  idate(4)
   integer              ::  ids(255)
   integer              ::  ipo(ko_),npo(ko_),pokpa(ko_)
   integer              ::  npt1(kt_),npt2(kt_),iz(ko_)
#ifdef RMP
   logical              ::  lbaseout
#endif
!
#include <postplevs.h>
!
   real                 ::  pot(255)
   data pot/255*0./
!-------------------------------------------------------------------------------
#ifdef RMP
#define IDIMF  igrd1_
#define JDIMF  jgrd1_
#define IDIMF2 igrd12_
#define JDIMF2 jgrd12_
#ifdef MP
#define LNWAVS lnwavp_
#define IDIM2 igrd12p_
#define JDIM2 jgrd12p_
#else
#define LNWAVS lnwav_
#define IDIM2 igrd1_
#define JDIM2 jgrd1_
#endif
#define IDIMO igrd1_
#define JDIMO jgrd1_
#else   
#ifdef REDUCE_GRID
#endif   /* RMP */
#define IDIMF  lonf_
#define JDIMF  latg_
#define IDIMF2 lonf2_
#define JDIMF2 latg2_
#ifdef MP
#define LNWAVS lnt22p_
#define IDIM2 lonf2p_
#define JDIM2 latg2p_
#else
#define LNWAVS lnt22_
#define IDIM2 lonf2_
#define JDIM2 latg2_
#endif
#define IDIMO io_
#define JDIMO jo_
#endif
!
   real                 ::  fgz(LNWAVS),fq(LNWAVS)
   real                 ::  fte(LNWAVS,levs_)
   real                 ::  fuu(LNWAVS,levs_),fvv(LNWAVS,levs_)
   real                 ::  frq(LNWAVS,levs_)
#ifdef RMP
   real                 ::  flat(IDIM2,JDIM2),flon(IDIM2,JDIM2)
#else
   real                 ::  flat(1,JDIMF)
#endif
!
   integer,parameter    ::  ntlen=nupa*ko_+nupt*kt_+nsun
   real                 ::  oxs(IDIM2,levs_),osxs(IDIM2,levs_)
   real                 ::  rxs(IDIM2,levs_),qsxs(IDIM2,levs_)
   real                 ::  zxs(IDIM2,levs_),zxi (IDIM2,levs_)
   real                 ::  qxp(IDIM2,kt_)
   real                 ::  fxp(IDIM2,ntlen)
!
   integer              ::  ipu(ntlen)
   integer              ::  itl(ntlen)
   integer              ::  il1(ntlen)
   integer              ::  il2(ntlen)
!
   integer,parameter    ::  ipuu=33,ipuv=34,ipuo=39,ipuz=7,iput=11,ipur=52,ipua=41
   integer,parameter    ::  ipuq=51,ipucl=153,ipuo3=154,ipupr=152
!
   logical              ::  lbm(IDIMF*JDIMF)
   integer              ::  ipusun(nsun)
   integer              ::  itlsun(nsun),il1sun(nsun),il2sun(nsun)
   integer              ::  kslp(2)
   integer              ::  iens(5)
!
   save iens
   data iens/1,0,0,1,255/
!
   logical              ::  lppr,lpcl,lpo3
   real, allocatable    ::  fxs(:,:,:),fxy(:,:,:)
#ifdef MP
   real, allocatable    ::  fullat(:,:),fullon(:,:)
#endif
#ifdef STDAMP_P
   real                 ::  fm2(*)
#endif
!
#ifdef MP
   integer              ::  nstr(0:npes_-1),nend(0:npes_-1)
#endif
   integer              ::  lgrib(ntlen)
#ifdef RMP
   character            ::  grib(30+lenpds+lengds+IDIMF*JDIMF*(16+1)/8,ntlen)
#else
   real                 ::  gout(IDIMO,JDIMO)
   character            ::  grib(30+lenpds+lengds+IDIMO*JDIMO*(16+1)/8,ntlen)
#endif
!
   integer              ::  iens2(5)
!
   real,parameter       ::  pt=30.,mxbit=16,icen=7,icen2=0,igen=99
#ifdef MP
!
   real, allocatable    ::  work(:,:,:)
#endif
#ifdef HYBRID
!
   real                 ::  ak5(levs_+1),bk5(levs_+1)
#endif
   real                 ::  sihyb(IDIM2,levs_+1),slhyb(IDIM2,levs_)
!
   character*6 fni
   character*80 fno
#ifdef MP
   character*3 cpe
#endif
!-------------------------------------------------------------------------------
   lpo3=.false.
   lpcl=.false.
   lppr=.false.
   lpo3=ngases_.ge.1
   lpcl=nwmass_.ge.2
   lppr=nwmass_.ge.3
!
!  set constants
!
   call idsdef(2,ids)
!
#ifdef MP
#ifndef MPIGRIB
   if(ntlen.ge.npes_) then
     nblock=max(ntlen/npes_,1)
     modnb=mod(ntlen,npes_)
     nlenmx=0
     nn=1
     do n = 0,npes_-1
       nstr(n)=nn
       if(n.le.modnb) then
         nblen=nblock+1
       else
         nblen=nblock
       endif
       nend(n)=min(nn+nblen-1,ntlen)
       nlen=nend(n)-nstr(n)+1
       nlenmx=max(nlenmx,nlen)
       nn=nn+nlen
     enddo
   else
     nn=1
     do n = 0,ntlen-1
       nstr(n)=nn
       nend(n)=nn
       nn=nn+1
     enddo
     do n = ntlen,npes_-1
       nstr(n)=0
       nend(n)=-1
     enddo
     nlenmx=1
   endif
   allocate (work(IDIMF2,JDIMF2,nlenmx))
#else
   allocate (work(IDIMF2,JDIMF2,1))
#endif
#endif
!
#ifdef RMP
   if(.not.lbaseout) then
     fni='r_pgb'
     nchi=5
   else
     fni='r_bpgb'
     nchi=6
   endif
#else
   fni='pgb'
   nchi=3
#endif
   call file_name(fni,nchi,fhour,fno,ncho)
#ifdef MP
#ifndef MPIGRIB
#ifndef PGBGATHER
   if(mype.lt.10) then
     ndig=2
     write(cpe,'(1h0,i1)') mype
   elseif(mype.lt.100) then
     ndig=2
     write(cpe,'(i2)') mype
   elseif(mype.lt.1000) then
     ndig=3
     write(cpe,'(i3)') mype
   endif
   fno=fno(1:ncho)//'_pe'//cpe
   ncho=ncho+3+ndig
#endif
#endif
#endif
#ifdef MP
#ifdef MPIGRIB
   if(mype.eq.master) then
#endif
#ifdef PGBGATHER
   if(mype.eq.master) then
#endif
#endif
   open(unit=n1,file=fno(1:ncho),form='unformatted',err=900)
   go to 901
  900 continue
   print *,' error in opening file ',fno(1:ncho)
#ifdef MP
#ifdef RMP
   call rmpabort
#else
   call mpabort
#endif
#else
   call abort
#endif
  901 continue 
   print *,' opening file ',fno(1:ncho)
#ifdef MP
#ifdef MPIGRIB
   endif
#endif
#ifdef PGBGATHER
   endif
#endif
#endif
!
!  set both input and output indices
!
#ifdef MP
   if(mype.eq.master) then
#endif
   print *,' set both input and output indices '
#ifdef MP
   endif
#endif
   ksz=1
   ksd=1+levs_
   kst=1+2*levs_
   ksq=1+3*levs_
   kscl=1+4*levs_
   kspr=1+5*levs_
   kso3=1+6*levs_
   kspsx=1+7*levs_
   kspsy=2+7*levs_
   ksu=3+7*levs_
   ksv=3+8*levs_
   ksps=3+9*levs_
   ksgz=4+9*levs_
   ksgzx=5+9*levs_
   ksgzy=6+9*levs_
   nflds=6+9*levs_
!
   kpz=1
   kpu=1+ko_
   kpv=1+2*ko_
   kpt=1+3*ko_
   kpo=1+4*ko_
   kpr=1+5*ko_
   kpq=1+6*ko_
   kpa=1+7*ko_
   kpo3=1+8*ko_
   kpcl=1+9*ko_
   kppr=1+10*ko_
   kptu=11*ko_+1
   kptv=11*ko_+kt_+1
   kptt=11*ko_+2*kt_+1
   kptr=11*ko_+3*kt_+1
   kpsun=11*ko_+4*kt_+1
!
   allocate (fxs(IDIM2,JDIM2,nflds))
   allocate (fxy(IDIM2,JDIM2,ntlen))
!
!  set some parameters
!
   jfhour=nint(fhour)
   nfldp=nupa*ko_+nupt*kt_+nsun
   nfldps=nupas*ko_
   do k = 1,ko_
     pokpa(k)=po(k)/10.
     if(float(nint(po(k))).eq.po(k).or.po(k).gt.655.) then
       ipo(k)=100
       npo(k)=nint(po(k))
     else
       ipo(k)=120
       npo(k)=nint(po(k)*100.)
     endif
   enddo
   ptkpa=pt/10.
   do k = 1,kt_
     npt1(k)=k*pt
     npt2(k)=(k-1)*pt
   enddo
   kpmu=isrchflt(ko_,po,1,pot(ipuu))
   kpmv=isrchflt(ko_,po,1,pot(ipuv))
   kpmo=isrchflt(ko_,po,1,pot(ipuo))
   kpmz=isrchflt(ko_,po,1,pot(ipuz))
   kpmt=isrchflt(ko_,po,1,pot(iput))
   kpmr=isrchflt(ko_,po,1,pot(ipur))
   kpmq=isrchflt(ko_,po,1,pot(ipuq))
   kpma=isrchflt(ko_,po,1,pot(ipua))
   kpmcl=isrchflt(ko_,po,1,pot(ipucl))
   kpmpr=isrchflt(ko_,po,1,pot(ipupr))
   kpmo3=isrchflt(ko_,po,1,pot(ipuo3))
   call post_pgb_2d_param(ko_,po,lpcl,lppr,ipusun,itlsun,il1sun,il2sun,        &
               kslp,kli)
!
!  set some grib parameters
!
   do nn = 1,nfldp
     ipu(nn)=0
   enddo
   do nn = kpu,kpu+kpmu-1
     ipu(nn)=ipuu
   enddo
   do nn = kpv,kpv+kpmv-1
     ipu(nn)=ipuv
   enddo
   do nn = kpo,kpo+kpmo-1
     ipu(nn)=ipuo
   enddo
   do nn = kpz,kpz+kpmz-1
     ipu(nn)=ipuz
   enddo
   do nn = kpt,kpt+kpmt-1
     ipu(nn)=iput
   enddo
   do nn = kpr,kpr+kpmr-1
     ipu(nn)=ipur
   enddo
   do nn = kpq,kpq+kpmq-1
     ipu(nn)=ipuq
   enddo
   do nn = kpa,kpa+kpma-1
     ipu(nn)=ipua
   enddo
   do nn = kpcl,kpcl+kpmcl-1
     if(lpcl) then
       ipu(nn)=ipucl
     else
       ipu(nn)=0
     endif
   enddo
   do nn = kppr,kppr+kpmpr-1
     if(lpcl) then
       ipu(nn)=ipupr
     else
       ipu(nn)=0
     endif
   enddo
   do nn = kpo3,kpo3+kpmo3-1
     if(lpcl) then
       ipu(nn)=ipuo3
     else
       ipu(nn)=0
     endif
   enddo
   do nn = kptu,kptu+kt_-1
     ipu(nn)=ipuu
   enddo
   do nn = kptv,kptv+kt_-1
     ipu(nn)=ipuv
   enddo
   do nn = kptt,kptt+kt_-1
     ipu(nn)=iput
   enddo
   do nn = kptr,kptr+kt_-1
     ipu(nn)=ipur
   enddo
   do nn = kpsun,kpsun+nsun-1
     ipu(nn)=ipusun(nn-kpsun+1)
   enddo
!
   do nn = kpu,kpu+ko_-1
     itl(nn)=ipo(nn-kpu+1)
   enddo
   do nn = kpv,kpv+ko_-1
     itl(nn)=ipo(nn-kpv+1)
   enddo
   do nn = kpo,kpo+ko_-1
     itl(nn)=ipo(nn-kpo+1)
   enddo
   do nn = kpz,kpz+ko_-1
     itl(nn)=ipo(nn-kpz+1)
   enddo
   do nn = kpt,kpt+ko_-1
     itl(nn)=ipo(nn-kpt+1)
   enddo
   do nn = kpr,kpr+ko_-1
     itl(nn)=ipo(nn-kpr+1)
   enddo
   do nn = kpq,kpq+ko_-1
     itl(nn)=ipo(nn-kpq+1)
   enddo
   do nn = kpa,kpa+ko_-1
     itl(nn)=ipo(nn-kpa+1)
   enddo
   do nn = kpo3,kpo3+ko_-1
     itl(nn)=ipo(nn-kpo3+1)
   enddo
   do nn = kpcl,kpcl+ko_-1
     itl(nn)=ipo(nn-kpcl+1)
   enddo
   do nn = kppr,kppr+ko_-1
     itl(nn)=ipo(nn-kppr+1)
   enddo
   do nn = nupa*ko_+1,nupa*ko_+nupt*kt_
     itl(nn)=116
   enddo
   do nn = kpsun,kpsun+nsun-1
     itl(nn)=itlsun(nn-kpsun+1)
   enddo
   do nn = 1,nupa*ko_+nupt*kt_+nsun
     il1(nn)=0
   enddo
   do nn = kptu,kptu+kt_-1
     il1(nn)=npt1(nn-kptu+1)
   enddo
   do nn = kptv,kptv+kt_-1
     il1(nn)=npt1(nn-kptv+1)
   enddo
   do nn = kptt,kptt+kt_-1
     il1(nn)=npt1(nn-kptt+1)
   enddo
   do nn = kptr,kptr+kt_-1
     il1(nn)=npt1(nn-kptr+1)
   enddo
   do nn = kpsun,kpsun+nsun-1
     il1(nn)=il1sun(nn-kpsun+1)
   enddo
   do nn = kpu,kpu+ko_-1
     il2(nn)=npo(nn-kpu+1)
   enddo
   do nn = kpv,kpv+ko_-1
     il2(nn)=npo(nn-kpv+1)
   enddo
   do nn = kpo,kpo+ko_-1
     il2(nn)=npo(nn-kpo+1)
   enddo
   do nn = kpz,kpz+ko_-1
     il2(nn)=npo(nn-kpz+1)
   enddo
   do nn = kpt,kpt+ko_-1
     il2(nn)=npo(nn-kpt+1)
   enddo
   do nn = kpr,kpr+ko_-1
     il2(nn)=npo(nn-kpr+1)
   enddo
   do nn = kpq,kpq+ko_-1
     il2(nn)=npo(nn-kpq+1)
   enddo
   do nn = kpa,kpa+ko_-1
     il2(nn)=npo(nn-kpa+1)
   enddo
   do nn = kpo3,kpo3+ko_-1
     il2(nn)=npo(nn-kpo3+1)
   enddo
   do nn = kpcl,kpcl+ko_-1
     il2(nn)=npo(nn-kpcl+1)
   enddo
   do nn = kppr,kppr+ko_-1
     il2(nn)=npo(nn-kppr+1)
   enddo
   do nn = kptu,kptu+kt_-1
     il2(nn)=npt2(nn-kptu+1)
   enddo
   do nn = kptv,kptv+kt_-1
     il2(nn)=npt2(nn-kptv+1)
   enddo
   do nn = kptt,kptt+kt_-1
     il2(nn)=npt2(nn-kptt+1)
   enddo
   do nn = kptr,kptr+kt_-1
     il2(nn)=npt2(nn-kptr+1)
   enddo
   do nn = kpsun,kpsun+nsun-1
     il2(nn)=il2sun(nn-kpsun+1)
   enddo
!
!  convert wave to grid
!
#ifdef RMP
#define COF2GRD rpost_wave2grid
#else
#define COF2GRD post_wave2grid_fcst
#endif
#ifdef RMP
   call COF2GRD(fgz,fq,fte,fuu,fvv,frq,fm2,                                    &
#else
   call COF2GRD(fgz,fq,fte,fuu,fvv,frq,                                        &
#endif
                fxs(1,1,ksgz),fxs(1,1,ksgzx),fxs(1,1,ksgzy),                   &
                fxs(1,1,ksps),fxs(1,1,kspsx),fxs(1,1,kspsy),                   &
                fxs(1,1,kst ),                                                 &
                fxs(1,1,ksu ),fxs(1,1,ksv  ),                                  &
                fxs(1,1,ksd ),fxs(1,1,ksz  ),                                  &
                fxs(1,1,ksq )                                                  &
#ifdef RMP
               ,lbaseout)
#else
                )
#endif
!
!  convert virtual temp to real temp
!
   do k = 1,levs_
     do j = 1,JDIM2
       do i = 1,IDIM2
         fxs(i,j,kst+k-1)=fxs(i,j,kst+k-1)/                                    &
                         (1.+fvirt*fxs(i,j,ksq+k-1))
       enddo
     enddo
   enddo
!
!  convert ln of ps to real ps
!
   do j = 1,JDIM2
     do i = 1,IDIM2
       fxs(i,j,ksps)=exp(fxs(i,j,ksps))
     enddo
   enddo
#ifdef DBG
   call print_maxmin_six(fxs(1,1,ksgz ),IDIM2*JDIM2, 1,1, 1,'gz')
   call print_maxmin_six(fxs(1,1,ksgzx),IDIM2*JDIM2, 1,1, 1,'dgz/dx')
   call print_maxmin_six(fxs(1,1,ksgzy),IDIM2*JDIM2, 1,1, 1,'dgz/dy')
   call print_maxmin_six(fxs(1,1,ksps ),IDIM2*JDIM2, 1,1, 1,'ps')
   call print_maxmin_six(fxs(1,1,kspsx),IDIM2*JDIM2, 1,1, 1,'dps/dx')
   call print_maxmin_six(fxs(1,1,kspsy),IDIM2*JDIM2, 1,1, 1,'dps/dy')
   call print_maxmin_six(fxs(1,1,kst  ),IDIM2*JDIM2,levs_,1,levs_,'tmp')
   call print_maxmin_six(fxs(1,1,ksu  ),IDIM2*JDIM2,levs_,1,levs_,'u')
   call print_maxmin_six(fxs(1,1,ksv  ),IDIM2*JDIM2,levs_,1,levs_,'v')
   call print_maxmin_six(fxs(1,1,ksd  ),IDIM2*JDIM2,levs_,1,levs_,'div')
   call print_maxmin_six(fxs(1,1,ksz  ),IDIM2*JDIM2,levs_,1,levs_,'vot')
   call print_maxmin_six(fxs(1,1,ksq  ),IDIM2*JDIM2,levs_,1,levs_,'q')
#endif
!
!  loop over groups of latitudes
!
!  compute auxiliary quantities on sigma and interpolate to pressure
!  and compute post_pgb_2d_diag fields and pack for transpose in parallel
!
#ifdef MP
!
!  fill undefined partial array points
!
   do n = 1,nflds
     do j = 1,JDIM2
       if(j.gt.latlen(mype)) then
         do i = 1,IDIM2
           fxs(i,j,n)=fxs(1,1,n)
         enddo
       endif
     enddo
     do i = 1,IDIM2
       if(i.gt.lonlen(mype)*2) then
         do j = 1,JDIM2
           fxs(i,j,n)=fxs(1,1,n)
         enddo
       endif
     enddo
   enddo
#ifdef RMP
   do j = 1,JDIM2
     if(j.gt.latlen(mype)) then
       do i = 1,IDIM2
         flat(i,j)=flat(1,1)
       enddo
     endif
   enddo
   do i = 1,IDIM2
     if(i.gt.lonlen(mype)*2) then
       do j = 1,JDIM2
         flat(i,j)=flat(1,1)
       enddo
     endif
   enddo
#endif
#endif
!
   j1=1
   j2=JDIM2
   ijlen=IDIM2*JDIM2
   do j = j1,j2
     jj=j-j1+1
     do i = 1,IDIM2
       do k = 1,levs_+1
         sihyb(i,k)=ak5(levs_-k+2)/fxs(i,jj,ksps)+bk5(levs_-k+2)
       enddo
       do k = 1,levs_
         slhyb(i,k)=0.5*(sihyb(i,k)+sihyb(i,k+1))
       enddo
     enddo
     do i = 1,IDIM2
       do k = 1,levs_+1
         sihyb(i,k)=si(k)
       enddo
       do k = 1,levs_
         slhyb(i,k)=sl(k)
       enddo
     enddo
!
     call post_get_rh(IDIM2,ijlen,levs_,slhyb,                                 &
              fxs(1,jj,ksps),fxs(1,jj,ksq),fxs(1,jj,kst),                      &
              qsxs,rxs)
     call post_get_omega(IDIM2,ijlen,levs_,sihyb,slhyb,                        &
              fxs(1,jj,ksps),fxs(1,jj,kspsx),fxs(1,jj,kspsy),                  &
              fxs(1,jj,ksd ),fxs(1,jj,ksu  ),fxs(1,jj,ksv  ),                  &
              oxs,osxs)
#ifndef RMP
     call post_get_height(IDIM2,ijlen,levs_,sihyb,slhyb,                       &
#else
     call rpost_hydro(IDIM2,ijlen,levs_,sihyb,slhyb,                           &
#endif
              fxs(1,jj,ksgz),fxs(1,jj,kst),fxs(1,jj,ksq),                      &
              zxs,zxi)
     call post_sigma2pressure(IDIM2,ijlen,levs_,sihyb,slhyb,fxs(1,jj,ksps),    &
              fxs(1,jj,ksu),fxs(1,jj,ksv),oxs,                                 &
              zxs,zxi,fxs(1,jj,kst),rxs,fxs(1,jj,ksq),                         &
              fxs(1,jj,ksz),                                                   &
              fxs(1,jj,kscl),fxs(1,jj,kspr),fxs(1,jj,kso3),                    &
              ko_,pokpa,lpcl,lppr,lpo3,                                        &
              fxp(1,kpu),fxp(1,kpv),fxp(1,kpo),                                &
              fxp(1,kpz),fxp(1,kpt),fxp(1,kpr),                                &
              fxp(1,kpq),fxp(1,kpa),fxp(1,kpcl),fxp(1,kppr),                   &
              fxp(1,kpo3))
     call post_sigma2p_avg(IDIM2,ijlen,levs_,sihyb,slhyb,fxs(1,jj,ksps),       &
               fxs(1,jj,ksu),fxs(1,jj,ksv),                                    &
               fxs(1,jj,kst),fxs(1,jj,ksq),qsxs,                               &
               kt_,ptkpa,                                                      &
               fxp(1,kptu),fxp(1,kptv),                                        &
               fxp(1,kptt),qxp,fxp(1,kptr))
     call post_pgb_2d_diag(IDIM2,ijlen,levs_,kslp,kli,lpcl,lppr,               &
               flat(1,j),sihyb,slhyb,kt_,ptkpa,                                &
               fxs(1,jj,ksgz),fxs(1,jj,ksps),osxs,                             &
               fxs(1,jj,ksu),fxs(1,jj,ksv),oxs,                                &
               fxs(1,jj,kst),rxs,fxs(1,jj,ksq),qsxs,                           &
               fxp(1,kptt),qxp,                                                &
               fxp(1,kpz),fxp(1,kpt),                                          &
               fxs(1,jj,kscl),fxs(1,jj,kspr),                                  &
               fxp(1,kpsun))
!
! compute absolute vorticity
!
     aomega=2.0*acos(-1.0)/(24.*60.*60.)
     do nn = 1,nfldp
       do i = 1,IDIM2 
         fxy(i,jj,nn)=fxp(i,nn) 
       enddo
     enddo
     do k = 1,ko_
       do i = 1,IDIM2
          fxy(i,jj,kpa-1+k)= fxy(i,jj,kpa-1+k)                                 &
#ifdef RMP
                        +2.0*aomega*sin(flat(i,j))
#else
                        +2.0*aomega*sin(flat(1,j))
#endif
        enddo
     enddo
!
#ifdef COUPLE_ROP
     do i = 1,IDIM2
        romssgz(i,jj)= fxs(i,jj,ksgz)
     enddo
#endif
!
   enddo
!
#ifdef RMP
   proj=rproj
   delx=rdelx
   dely=rdely
#ifdef MP
   allocate (fullat(IDIMF,JDIMF),fullon(IDIMF,JDIMF))
   call rmpgp2f(flat,igrd12p_,jgrd12p_,fullat,igrd12_,jgrd12_,1)
   call rmpgp2f(flon,igrd12p_,jgrd12p_,fullon,igrd12_,jgrd12_,1)
#define FLAT fullat
#define FLON fullon
#else
#define FLAT flat
#define FLON flon
#endif
   call rmp_trans2output_grid(FLAT,1)
   call rmp_trans2output_grid(FLON,1)
   rlat1=FLAT(1,1)
   rlat2=FLAT(IDIMF,JDIMF)
   rlon1=FLON(1,1)
   rlon2=FLON(IDIMF,JDIMF)
#ifndef MP
   call rmp_trans2model_grid(FLAT,1)
   call rmp_trans2model_grid(FLON,1)
#else
   deallocate (fullat,fullon)
#endif
   if( proj.eq. 0. ) then
     idrt=1                    ! mercater
     ortru=rtruth
#ifdef MP
     if(mype.eq.master) then
#endif
       print *,' mercater projection.'
#ifdef MP
     endif
#endif
   elseif( abs(proj).eq.1. ) then
     idrt=5                    ! polar projection
     ortru=rorient
#ifdef MP
     if(mype.eq.master) then
#endif
       print *,' polar projection.'
#ifdef MP
     endif
#endif
   else
     idrt=0
#ifdef MP
     if(mype.eq.master) then
#endif
       print *,' undefine map projection.'
#ifdef MP
     endif
#endif
   endif
#else
   proj=0.
   delx=0.
   dely=0.
   ortru=0.
   rlat1=0.
   rlon1=0.
   rlat2=0.
   rlon2=0.
#endif
#ifndef MP
!
!  loop over groups of horizontal fields
!
   do k1=1,nfldp,ntlen
     k2=min(k1+ntlen-1,nfldp)
!
!  unpack transposed fields and interpolate to output grid
!  and round to the number of bits and engrib the field in parallel
!
#ifdef ORIGIN
!$doacross share(k1,k2,pot,nfldp,fxy,
!$&       ipu,itl,il1,il2,icen,igen,idate,lenpds,jfhour,ids,mxbit,
!$&       grib,lgrib,itruth,idelx,idely,
!$&       ilat1,ilon1,ilat2,ilon2,idrt,
!$&       rmax,rmin,icen2,rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj),
!$& local(k,kan,ierr,pok,poktop,i,j)
!
#endif
#ifdef OPENMP
!$omp parallel do private(k,kan,ierr,pok,poktop,i,j)
!
#endif
#endif
#ifdef MP
#define FXY work
#else
#define FXY fxy
#endif
!
#ifndef MP
     do k = k1,k2
       mype=0
       kan=k-k1+1
#else
#ifdef MPIGRIB
     do k = 1,ntlen
       kan=1
#ifndef RMP
       call mpgp2f(fxy(1,1,k),lonf2p_,latg2p_,work,lonf2_,latg2_,1)
#else
       call rmpgp2f(fxy(1,1,k),igrd12p_,jgrd12p_,work,                         &
                    igrd12_,jgrd12_,1)
#endif
       if(mype.eq.master) then
#else  /* ifndef MPIGRIB */
!
#ifdef RMP
#define MPGPFK2FPK rmpgpfk2fpk
#else
#define MPGPFK2FPK mpgpfk2fpk
#endif
       call MPGPFK2FPK(fxy,IDIM2,JDIM2,nfldp,                                  &
                       work,IDIMF2,JDIMF2,nlenmx,                              &
                       nstr,nend)
       do k = nstr(mype),nend(mype)
         kan=k-nstr(mype)+1
#endif  /* ifdef MPIGRIB */
!
#endif  /* ifdef MP */
         lgrib(kan)=0
         if(ipu(k).gt.0.and.itl(k).ne.107.and.itl(k).ne.116) then
#ifdef DBG
           rmax=FXY(1,1,kan)
           rmin=FXY(1,1,kan)
#ifdef MP
           do j = 1,JDIMF2
#ifdef REDUCE_GRID
             lonsd2=lonfdp(j,mype)*2
#else
             lonsd2=IDIMF2
#endif /* REDUCE_GRID end */
             do i = 1,lonsd2
#else 
           do j = 1,JDIM2
#ifdef REDUCE_GRID
             lonsd2=lonfd(latdef(j))*2
#else
             lonsd2=IDIM2
#endif
             do i = 1,lonsd2
#endif  /* #else ifdef MP  */
               rmax=max(FXY(i,j,kan),rmax)
               rmin=min(FXY(i,j,kan),rmin)
             enddo
           enddo
#endif
#ifndef RMP
           call dyn_trans2output_grid(FXY(1,1,kan),1)
           call post_gauss2latlon(FXY(1,1,kan),IDIMF,JDIMF,                    &
                       0.,90.,360./float(IDIMO),180./float(JDIMO-1),           &
                       gout,IDIMO,JDIMO)
           call file_make_grib(gout        ,lbm,   0,io_,jo_,mxbit,90.,        &
#else
           call rmp_trans2output_grid(FXY(1,1,kan),1)
           call file_make_grib(FXY(1,1,kan),lbm,idrt,igrd1_,jgrd1_,mxbit,0.,   &
#endif
                     lenpds,2,icen,igen,0,                                     &
                     ipu(k),itl(k),il1(k),il2(k),                              &
                     idate(4),idate(2),idate(3),idate(1),                      &
                     1,jfhour,0,10,                                            &
                     0,0,icen2,ids(ipu(k)),iens,                               &
                     rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,             &
                     truth,cotru,                                              &
                     grib(1,kan),lgrib(kan),ierr)
#ifdef DBG
           write(*,'("mype=",i3," k=",i4," max=",                              &
                 e17.8," min=",e17.8," lgrib=",i8)')                           &
                   mype,k,rmax,rmin,lgrib(kan)
#endif
#ifdef MPIGRIB 
           if((.not.lpo3).and.k.ge.kpo3.and.k.lt.kpcl) lgrib(kan) = 0
           if((.not.lpcl).and.k.ge.kpcl.and.k.lt.kppr) lgrib(kan) = 0
           if((.not.lppr).and.k.ge.kppr.and.k.lt.kptu) lgrib(kan) = 0
           if(lgrib(kan).gt.0) then
             call file_write_byte(n1,lgrib(kan),grib(1,kan))
#ifdef DBG
             print *,' grib1 written to ',n1,'length=',lgrib(kan)
#endif
           endif
#endif  /* ifdef MPIGRIB */
         endif
#ifdef MP
#ifdef MPIGRIB
         endif
#endif  /* ifdef MPIGRIB  */
#endif  /* ifdef MP */
       enddo
!
#ifdef MP
#ifdef PGBGATHER
       call mpgfpk2fk(grib,30+lenpds+lengds+IDIMO*JDIMO*(16+1)/8,              &
                      lgrib,nstr,nend,grib,ntlen)
!
#endif
#endif
#ifndef MPIGRIB
#ifdef MP 
#ifndef PGBGATHER
     do k = nstr(mype),nend(mype)
       kan=k-nstr(mype)+1
#else
     if(mype.eq.master) then
     do k = 1,ntlen
       kan=k
#endif
#else
     do k = k1,k2
       kan=k-k1+1
#endif
       if((.not.lpo3).and.k.ge.kpo3.and.k.lt.kpcl) lgrib(kan) = 0
       if((.not.lpcl).and.k.ge.kpcl.and.k.lt.kppr) lgrib(kan) = 0
       if((.not.lppr).and.k.ge.kppr.and.k.lt.kptu) lgrib(kan) = 0
       if(lgrib(kan).gt.0) then
         call file_write_byte(n1,lgrib(kan),grib(1,kan))
#ifdef DBG
         print *,' grib1 written to ',n1,' of length ',lgrib(kan)
#endif
       endif
     enddo
#endif  /* ifndef MPIGRIB */
#ifdef MP
#ifdef PGBGATHER
     endif
#endif
#endif
#ifndef MP
   enddo
#endif
   deallocate (fxs,fxy)
!
!  flx fields
!
!  for global model, separate pgb and flx files
!
#ifndef RMP
#ifdef MP
#ifdef MPIGRIB
   if(mype.eq.master) then
#endif
#ifdef PGBGATHER
   if(mype.eq.master) then
#endif
#endif
     close (n1)
#ifdef MP
#ifdef MPIGRIB
   endif
#endif
#ifdef PGBGATHER
   endif
#endif
#endif
#endif
!
#ifdef MP
#ifndef MPIGRIB
   deallocate (work)
#endif
#endif
!
#ifndef COUPLE_ROP
   call file_write_flux(n1)
#else
   call file_write_flux(n1,nnp)
#endif
!
#endif
   return
   end subroutine post_pgb_2d
!
!-------------------------------------------------------------------------------
   subroutine post_pgb_2d_param(ko,po,lpcl,lppr,                               &
                     ipusun,itlsun,ip1sun,ip2sun,kslp,kli)
!-------------------------------------------------------------------------------
!
! subprogram:    post_pgb_2d_param      set parameters for sundry fields
!
! abstract: sets parameters for the sundry fields.
!   parameters returned are parameter indicator, level type indicator
!   and two level numbers.
!   the current nsun=29 sundry fields are:
!     1) surface pressure
!     2) surface pressure tendency
!     3) column precipitable water
!     4) column relative humidity
!     5) tropopause temperature
!     6) tropopause pressure
!     7) tropopause zonal wind
!     8) tropopause meridional wind
!     9) tropopause vertical wind speed shear
!    10) surface lifted index
!    11) best lifted index
!    12) maximum wind level temperature
!    13) maximum wind level pressure
!    14) maximum wind level zonal wind
!    15) maximum wind level meridional wind
!    16) surface orography
!    17) sea level pressure
!    18) relative humidity in sigma range (0.44,1.00)
!    19) relative humidity in sigma range (0.72,0.94)
!    20) relative humidity in sigma range (0.44,0.72)
!    21) potential temperature at sigma 0.9950
!    22) temperature at sigma 0.9950
!    23) pressure vertical velocity at sigma 0.9950
!    24) relative humidity at sigma 0.9950
!    25) zonal wind at sigma 0.9950
!    26) meridional wind at sigma 0.9950
!    27) specific humidity at sigma 0.9950
!    28) total cloud water
!    29) total precipitation water
!
! program history log:
!   1992-01-01  mark iredell           initial mrf
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call post_pgb_2d_param(ko,po,ipusun,itlsun,ip1sun,ip2sun,kslp,kli)
!
!   input argument list:
!     ko       - integer number of pressure levels
!     po       - real (ko) pressure in millibars
!
!   output argument list:
!     ipusun   - integer (nsun) parameter indicators
!     itlsun   - integer (nsun) level type indicators
!     ip1sun   - integer (nsun) first level numbers
!     ip2sun   - integer (nsun) second level numbers
!     kslp     - integer (2) relevant pressure levels for slp
!     kli      - integer relevant pressure level for lifted index
!
! subprograms called:
!   isrcheq  - find first value in an array equal to target value
!
!-------------------------------------------------------------------------------
   real                ::  po(ko)
   logical             ::  lpcl,lppr,lpo3
   integer, parameter  ::  nsun=29
   integer, parameter  ::  nps  = 1,npst = 2,ntpw = 3,ntrh = 4,ntpt = 5,       &
                           ntpp = 6,ntpu = 7,ntpv = 8,ntpsh= 9,nsli =10,       &
                           nbli =11,nmwt =12,nmwp =13,nmwu =14,nmwv =15,       &
                           nzs  =16,nslp =17,nrh1 =18,nrh2 =19,nrh3 =20,       &
                           ns1th=21,ns1t =22,ns1o =23,ns1r =24,ns1u =25,       &
                           ns1v =26,ns1q =27,ntcl =28,ntpr =29
   integer             ::  ipusun(nsun),itlsun(nsun)
   integer             ::  ip1sun(nsun),ip2sun(nsun)
   integer             ::  kslp(2)
   integer             ::  ipudef(nsun),itldef(nsun)
   integer             ::  ip1def(nsun),ip2def(nsun)
   real                ::  pslp(2)
   data ipudef/001,003,054,052,011,001,033,034,136,131,                        &
               132,011,001,033,034,007,002,052,052,052,                        &
               013,011,039,052,033,034,51,76,77/
   data itldef/001,001,200,200,007,007,007,007,007,001,                        &
               001,006,006,006,006,001,102,108,108,108,                        &
               107,107,107,107,107,107,107,200,200/
   data ip1def/000,000,000,000,000,000,000,000,000,000,                        &
               000,000,000,000,000,000,000,044,072,044,                        &
               00000,00000,00000,00000,00000,00000,00000,                      &
               000,000/
   data ip2def/000,000,000,000,000,000,000,000,000,000,                        &
               000,000,000,000,000,000,000,100,094,072,                        &
               09950,09950,09950,09950,09950,09950,09950,                      &
               000,000/
   data pslp/1000.,500./,pli/500./
!-------------------------------------------------------------------------------
! 
   ipusun=ipudef
   itlsun=itldef
   ip1sun=ip1def
   ip2sun=ip2def
   kslp(1)=mod(isrcheq(ko,po,1,pslp(1)),ko+1)
   kslp(2)=mod(isrcheq(ko,po,1,pslp(2)),ko+1)
   kli=mod(isrcheq(ko,po,1,pli),ko+1)
   if(kslp(1).eq.0.or.kslp(2).eq.0) ipusun(nslp)=0
   if(kli.eq.0) ipusun(nsli)=0
   if(kli.eq.0) ipusun(nbli)=0
   if(.not.lpcl) ipusun(ntcl)=0
   if(.not.lppr) ipusun(ntpr)=0
!
!  do not produce low-layer rh because of a difficulty in grib 
!
   ipusun(nrh1)=0
   ipusun(nrh2)=0
   ipusun(nrh3)=0
!
   return
   end subroutine post_pgb_2d_param
!
!-------------------------------------------------------------------------------
   subroutine post_pgb_2d_diag(im,ix,km,kslp,kli,lpcl,lppr,                    &
                     clat,si,sl,kt,pt,                                         &
                     zs,ps,os,u,v,o,t,r,q,qs,tpt,qpt,zm,tm,                    &
                     cl,pr,sun)
!-------------------------------------------------------------------------------
!
! abstract: computes post_pgb_2d_diag fields.
!   the current nsun=29 post_pgb_2d_diag fields are:
!     1) surface pressure
!     2) surface pressure tendency
!     3) column precipitable water
!     4) column relative humidity
!     5) tropopause temperature
!     6) tropopause pressure
!     7) tropopause zonal wind
!     8) tropopause meridional wind
!     9) tropopause vertical wind speed shear
!    10) surface lifted index
!    11) best lifted index
!    12) maximum wind level temperature
!    13) maximum wind level pressure
!    14) maximum wind level zonal wind
!    15) maximum wind level meridional wind
!    16) surface orography
!    17) sea level pressure
!    18) relative humidity in sigma range (0.44,1.00)
!    19) relative humidity in sigma range (0.72,0.94)
!    20) relative humidity in sigma range (0.44,0.72)
!    21) potential temperature at sigma 0.9950
!    22) temperature at sigma 0.9950
!    23) pressure vertical velocity at sigma 0.9950
!    24) relative humidity at sigma 0.9950
!    25) zonal wind at sigma 0.9950
!    26) meridional wind at sigma 0.9950
!    27) specific humidity at sigma 0.9950
!    28) total cloud water
!    29) total ozone
!
! subprograms called:
!   post_sigma2tropopause    interpolate sigma to tropopause level
!   post_lifting_index       compute best lifted index
!   post_max_wind            interpolate sigma to maxwind level
!
! program history log:
!
! usage:    call post_pgb_2d_diag(im,ix,km,kslp,kli,clat,si,sl,kt,pt,
!    &                  zs,ps,os,u,v,o,t,r,q,qs,tpt,qpt,zm,tm,sun)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of levels
!     kslp     - integer (2) relevant pressure levels for slp
!     kli      - integer relevant pressure level for lifted index
!     clat     - real (im) cosine of latitude
!     si       - real (km) sigma interfaces
!     sl       - real (km) sigma values
!     kt       - integer number of pressure thickness layers
!     pt       - real pressure thickness in kpa
!     zs       - real (im) surface orography in m
!     ps       - real (im) surface pressure in kpa
!     os       - real (im) surface pressure tendency in pa/s
!     u        - real (ix,km) zonal wind in m/s
!     v        - real (ix,km) meridional wind in m/s
!     o        - real (im,km) vertical velocity in pa/s
!     t        - real (ix,km) temperature in k
!     r        - real (im,km) relative humidity in percent
!     q        - real (ix,km) specific humidity in kg/kg
!     qs       - real (im,km) saturated specific humidity in kg/kg
!     tpt      - real (im,kt) temperature in k
!     qpt      - real (im,kt) specific humidity in kg/kg
!     zm       - real (im,*) height on pressure surface in m
!     tm       - real (im,*) temperature on pressure surface in k
!     cl       - real (ix,*) cloud water in kg/kg
!     pr       - real (ix,*) precipitation water in kg/kg
!
!   output argument list:
!     sun      - real (im,nsun) post_pgb_2d_diag fields given above
!
!-------------------------------------------------------------------------------
   use constant, only : g=>g_, rd=>rd_,cp=>cp_
!-------------------------------------------------------------------------------
   integer,parameter    ::  nsun=29
   integer,parameter    ::  nps  = 1,npst = 2,ntpw = 3,ntrh = 4,ntpt = 5
   integer,parameter    ::  ntpp = 6,ntpu = 7,ntpv = 8,ntpsh= 9,nsli =10
   integer,parameter    ::  nbli =11,nmwt =12,nmwp =13,nmwu =14,nmwv =15
   integer,parameter    ::  nzs  =16,nslp =17,nrh1 =18,nrh2 =19,nrh3 =20
   integer,parameter    ::  ns1th=21,ns1t =22,ns1o =23,ns1r =24,ns1u =25
   integer,parameter    ::  ns1v =26,ns1q =27,ntcl =28,ntpr=29
!
! os,o,r,qs
!
   integer              ::  kslp(2)
   logical              ::  lpcl,lppr,lpo3
   real                 ::  si(im,km+1),sl(im,km)
   real                 ::  zs(im),ps(im),os(im),clat(im)
   real                 ::  u(ix,km),v(ix,km),o(im,km)
   real                 ::  t(ix,km),r(im,km),q(ix,km),qs(im,km)
   real                 ::  cl(ix,km),pr(ix,km)
   real                 ::  tpt(im,kt),qpt(im,kt)
   real                 ::  zm(im,*),tm(im,*)
   real                 ::  sun(im,nsun)
   real                 ::  wrk(im)
   real,parameter       ::  rocp=rd/cp
   real,parameter       ::  pm1=1.e5,tm1=287.45,zm1=113.,zm2=5572.
   real,parameter       ::  fslp=g*(zm2-zm1)/(rd*tm1)
   real,parameter       ::  strh1=0.44,strh2=0.72,strh3=0.44
   real,parameter       ::  sbrh1=1.00,sbrh2=0.94,sbrh3=0.72
   real,parameter       ::  sl1=0.9950
!-------------------------------------------------------------------------------
!
!  surface orography, surface pressure and surface pressure tendency
!
   do i = 1,im
     sun(i,nzs)=zs(i)
     sun(i,nps)=ps(i)*1.e3
     sun(i,npst)=os(i)
   enddo
! 
!  column precipitable water and relative humidity
!
   do i = 1,im
     sun(i,ntpw)=0.
     wrk(i)=0.
   enddo
   do k = 1,km
     do i = 1,im
       ds=si(i,k)-si(i,k+1)
       sun(i,ntpw)=sun(i,ntpw)+q(i,k)*ds
       wrk(i)=wrk(i)+qs(i,k)*ds
     enddo
   enddo
   do i = 1,im
     sun(i,ntrh)=min(max(sun(i,ntpw)/wrk(i),0.),1.)*100.
     sun(i,ntpw)=sun(i,ntpw)*ps(i)*1.e3/g
   enddo
!
!  post_pgb_2d_diag tropopause fields
!
   call post_sigma2tropopause(im,ix,km,clat,sl,ps,u,v,t,                       &
               sun(1,ntpp),sun(1,ntpu),sun(1,ntpv),                            &
               sun(1,ntpt),sun(1,ntpsh))
   do i = 1,im
     sun(i,ntpp)=sun(i,ntpp)*1.e3
   enddo
! 
!  lifted index
!
   if(kli.gt.0) then
     call post_lifting_index(im,ix,kt,pt,ps,tpt,qpt,tm(1,kli),                 &
                 sun(1,nsli),sun(1,nbli))
   else
     do i = 1,im
       sun(i,nsli)=0.
       sun(i,nbli)=0.
     enddo
   endif
!
!  post_pgb_2d_diag maxwind fields
!
   call post_max_wind(im,ix,km,sl,ps,u,v,t,                                    &
               sun(1,nmwp),sun(1,nmwu),sun(1,nmwv),sun(1,nmwt))
   do i = 1,im
     sun(i,nmwp)=sun(i,nmwp)*1.e3
   enddo
! 
!  sea level pressure
!
   if(kslp(1).gt.0.and.kslp(2).gt.0) then
     k1=kslp(1)
     k2=kslp(2)
     do i = 1,im
       sun(i,nslp)=pm1*exp(fslp*zm(i,k1)/(zm(i,k2)-zm(i,k1)))
     enddo
   else
     do i = 1,im
       sun(i,nslp)=0.
     enddo
   endif
!
!  average relative humidity 1
!
   do i = 1,im
     sun(i,nrh1)=0.
     wrk(i)=0.
   enddo
   do k = 1,km
     do i = 1,im
       ds=min(si(i,k),sbrh1)-max(si(i,k+1),strh1)
       if(ds.gt.0.) then
         sun(i,nrh1)=sun(i,nrh1)+q(i,k)*ds
         wrk(i)=wrk(i)+qs(i,k)*ds
       endif
     enddo
   enddo
   do i = 1,im
     sun(i,nrh1)=min(max(sun(i,nrh1)/wrk(i),0.),1.)*100.
   enddo
! 
!  average relative humidity 2
!
   do i = 1,im
     sun(i,nrh2)=0.
     wrk(i)=0.
   enddo
   do k = 1,km
     do i = 1,im
       ds=min(si(i,k),sbrh2)-max(si(i,k+1),strh2)
       if(ds.gt.0.) then
         sun(i,nrh2)=sun(i,nrh2)+q(i,k)*ds
         wrk(i)=wrk(i)+qs(i,k)*ds
       endif
     enddo
   enddo
   do i = 1,im
     sun(i,nrh2)=min(max(sun(i,nrh2)/wrk(i),0.),1.)*100.
   enddo
!
!  average relative humidity 3
!
   do i = 1,im
     sun(i,nrh3)=0.
     wrk(i)=0.
   enddo
   do k = 1,km
     do i = 1,im
       ds=min(si(i,k),sbrh3)-max(si(i,k+1),strh3)
       if(ds.gt.0.) then
         sun(i,nrh3)=sun(i,nrh3)+q(i,k)*ds
         wrk(i)=wrk(i)+qs(i,k)*ds
       endif
     enddo
   enddo
   do i = 1,im
     sun(i,nrh3)=min(max(sun(i,nrh3)/wrk(i),0.),1.)*100.
   enddo
!
!  bottom sigma fields interpolated from first two model layers
!
   do i = 1,im
     f2=log(sl(i,1)/sl1)/log(sl(i,1)/sl(i,2))
     sl1k=(sl1*1.e-2)**(-rocp)
     sun(i,ns1t)=t(i,1)+f2*(t(i,2)-t(i,1))
     sun(i,ns1th)=sun(i,ns1t)*sl1k*ps(i)**(-rocp)
     sun(i,ns1o)=o(i,1)+f2*(o(i,2)-o(i,1))
     sun(i,ns1r)=r(i,1)+f2*(r(i,2)-r(i,1))
     sun(i,ns1u)=u(i,1)+f2*(u(i,2)-u(i,1))
     sun(i,ns1v)=v(i,1)+f2*(v(i,2)-v(i,1))
     sun(i,ns1q)=q(i,1)+f2*(q(i,2)-q(i,1))
   enddo
! 
!  total cloud water
!
   if(lpcl) then
     do i = 1,im
       sun(i,ntcl)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,ntcl)=sun(i,ntcl)+cl(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,ntcl)=max(sun(i,ntcl),0.)*ps(i)*1.e3/g
     enddo
   endif
!
!  total precipatation water
!
   if(lppr) then
     do i = 1,im
       sun(i,ntpr)=0.
     enddo
     do k = 1,km
       do i = 1,im
         ds=si(i,k)-si(i,k+1)
         sun(i,ntpr)=sun(i,ntpr)+pr(i,k)*ds
       enddo
     enddo
     do i = 1,im
       sun(i,ntpr)=max(sun(i,ntpr),0.)*ps(i)*1.e3/g
     enddo
   endif
!
   return
   end subroutine post_pgb_2d_diag
!
!-------------------------------------------------------------------------------
   subroutine post_lifting_index(im,ix,kt,pt,ps,t,q,tm,sli,bli)
!-------------------------------------------------------------------------------
!
! abstract
!   - compute lifted indices from sigma
!   - computes both the surface lifted index and best lifted index
!   from profiles in constant pressure thickness layers above ground.
!   the surface lifted index is computed by raising the lowest layer
!   to 500 mb and subtracting its parcel temperature
!   from the environment temperature.
!   the best lifted index is computed by finding the parcel
!   with the warmest equivalent potential temperature,
!   then raising it to 500 mb and subtracting its parcel temperature
!   from the environment temperature.
!
! program history log:
!   1992-10-31  iredell                development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!
! usage:    call post_lifting_index(im,ix,kt,pt,ps,t,q,tm,sli,bli)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     kt       - integer number of layers in profile
!     pt       - real pressure thickness in kpa
!     ps       - real (im) surface pressure in kpa
!     t        - real (ix,kt) temperature in k
!     q        - real (ix,kt) specific humidity in kg/kg
!     tm       - real (im) 500 mb temperature in k
!
!   output argument list:
!     sli      - real (ix) surface lifted index in k
!     bli      - real (ix) best lifted index in k
!
! subprograms called:
!   (fpkap)   - function to compute pressure to the kappa
!   (ftdp)    - function to compute dewpoint temperature
!   (ftlcl)   - function to compute lifting condensation level
!   (fthe)    - function to compute equivalent potential temperature
!   (ftma)    - function to compute moist adiabat temperature
!
!fpp$ expand(fpkap,ftdp,ftlcl,fthe,ftma)
!-------------------------------------------------------------------------------
   use constant, only : cp_,rd_,rv_
   real,parameter  ::  rk=rd_/cp_,eps=rd_/rv_,epsm1=rd_/rv_-1.
   real,parameter  ::  plift=50.
   real            ::  ps(im),t(im,kt),q(im,kt),tm(im),sli(im),bli(im)
   real            ::  p2kmas(im),themas(im),p2kmab(im),themab(im)
!-------------------------------------------------------------------------------
!
!  select the warmest equivalent potential temperature
!
   do k = 1,kt
     do i = 1,im
       p=ps(i)-(k-0.5)*pt
       pv=p*q(i,k)/(eps-epsm1*q(i,k))
       tdpd=max(t(i,k)-ftdp(pv),0.)
       tlcl=ftlcl(t(i,k),tdpd)
       p2klcl=fpkap(p)*tlcl/t(i,k)
       thelcl=fthe(tlcl,p2klcl)
       if(k.eq.1) then
         p2kmas(i)=p2klcl
         themas(i)=thelcl
         p2kmab(i)=p2klcl
         themab(i)=thelcl
       elseif(thelcl.gt.themab(i)) then
         p2kmab(i)=p2klcl
         themab(i)=thelcl
       endif
     enddo
   enddo
! 
!  lift the parcel to 500 mb along a dry adiabat below the lcl
!  or along a moist adiabat above the lcl.
!  the lifted index is the environment minus parcel temperature.
!
   pliftk=(plift/100.)**rk
   do i = 1,im
     if(ps(i).gt.plift) then
       p2ks=min(pliftk,p2kmas(i))
       sli(i)=tm(i)-pliftk/p2ks*ftma(themas(i),p2ks,qma)
       p2kb=min(pliftk,p2kmab(i))
       bli(i)=tm(i)-pliftk/p2kb*ftma(themab(i),p2kb,qma)
     else
       sli(i)=0.
       bli(i)=0.
     endif
   enddo
!
   return
   end subroutine post_lifting_index
!
!-------------------------------------------------------------------------------
   subroutine post_max_wind(im,ix,km,sl,ps,u,v,t,pmw,umw,vmw,tmw)                   
!-------------------------------------------------------------------------------
!
! abstract 
!   - sigma to maxwind interpolation
!   - locates the maximum wind speed level (maxwind level) and            
!   returns the wind speed, components and pressure at that level.              
!   the maxwind level is restricted to be between 50kpa and 7kpa.               
!   the maxwind level is identified by cubic spline interpolation               
!   of the wind speeds in log pressure.                                         
!                                                                               
! usage:    call post_max_wind(im,ix,km,sl,ps,u,v,t,pmw,umw,vmw,tmw)  
!                                                                               
!   input argument list:                                                        
!     im       - integer number of points                                       
!     ix       - integer first dimension of upper air data                      
!     km       - integer number of sigma levels                                 
!     sl       - real (km) sigma values                                         
!     ps       - real (im) surface pressure in kpa                              
!     u        - real (ix,km) zonal wind in m/s                                 
!     v        - real (ix,km) merid wind in m/s                                 
!     t        - real (ix,km) temperature in k                                  
!                                                                               
!   output argument list:                                                       
!     pmw      - real (im) maxwind pressure in kpa                              
!     umw      - real (im) maxwind zonal wind in m/s                            
!     vmw      - real (im) maxwind merid wind in m/s                            
!     tmw      - real (im) maxwind temperature in k                             
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
   real, parameter  ::  pmwbot=500.e-1,pmwtop=70.e-1                                   
!
   real             ::  sl(im,km),ps(im)                                                   
   real             ::  u(ix,km),v(ix,km),t(ix,km)                                      
   real             ::  pmw(im),umw(im),vmw(im),tmw(im)                                 
   real             ::  spdmw(im),smw(im),s(levs_)                              
   real             ::  spd(im,levs_),d2spd(im,levs_)                           
!-------------------------------------------------------------------------------
!
!  fix vertical coordinate proportional to log pressure                         
!  and calculate wind speeds between pmwbot and pmwtop                          
!
   do k = 1,km                                                                 
     do i = 1,im                                                               
       s(k)=-log(sl(i,k))                                                        
       p=sl(i,k)*ps(i)                                                         
       if(p.le.pmwbot.and.p.ge.pmwtop) then                                  
         spd(i,k)=sqrt(u(i,k)**2+v(i,k)**2)                                  
       else                                                                  
         spd(i,k)=0.                                                         
       endif                                                                 
     enddo                                                                   
   enddo                                                                     
!
!  use spline routines to determine maxwind level and wind speed                
!
   call post_spline_derivative(im,km,s,spd,d2spd)         
   call post_spline_max(im,km,s,spd,d2spd,smw,pmw,spdmw)  
!
!  compute maxwind pressure and wind components                                 
!
   do i = 1,im                                                                 
     pmw(i)=exp(-pmw(i))*ps(i)                                               
     k=int(smw(i))                                                           
     if(float(k).eq.smw(i)) then                                             
       ub=u(i,k)                                                             
       vb=v(i,k)                                                             
       tb=t(i,k)                                                             
     else                                                                    
       ub=(k+1-smw(i))*u(i,k)+(smw(i)-k)*u(i,k+1)                            
       vb=(k+1-smw(i))*v(i,k)+(smw(i)-k)*v(i,k+1)                            
       tb=(k+1-smw(i))*t(i,k)+(smw(i)-k)*t(i,k+1)                            
     endif                                                                   
     spdb=sqrt(ub**2+vb**2)                                                  
     umw(i)=ub*spdmw(i)/spdb                                                 
     vmw(i)=vb*spdmw(i)/spdb                                                 
     tmw(i)=tb                                                               
   enddo                                                                     
!
   return                                                                    
   end subroutine post_max_wind
!
!-------------------------------------------------------------------------------
