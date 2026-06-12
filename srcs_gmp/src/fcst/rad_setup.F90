#include "define.h"
   subroutine rad_setup
!-------------------------------------------------------------------------------
!
!  ::: structure ::: 
!
!    [rad_initialize]
!        |---  [rad_cloud_read] *
!        |---  [rad_cloudiness_init] *
!        |---  [rad_ozone_co2_init] ----- [rad_constant_init] *
!        |                            |-- [rad_ozone_interp_step1] *
!        |                            |-- [rad_ozone_interp_step1] *
!        |                            |-- [rad_co2_read] - [rad_lw_table] *
!        |                            |-- for NIM --- [rad_co2_read1]
!        |                                         |- [rad_co2_read2]
!        |--- A-1) [rad_aeros_init] *
!        |    A-2) [rad_aeros_init_ice] *
!        |--- B)   [rad_aeros_init_nasa] *
!        |
!        |--- [rad_albedo_aeros_read] *
!
!    [dyn_sph_driver] or [dyn_dfs_driver]
!        |
!        |--- [rad_prepare] *
!                 |---  [rad_albedo_aeros_read] *
!                   |-  [rmp_albedo_aeros_read] *
!                 |---  [rad_julian_day] *
!                 |---  [rad_fcst_tim] *
!                 |---  [rad_solar_time] *
!                 |---  [rad_cos_zenith] *
!                 |---  [rad_calc_time] *
!                 |---  [rad_print_time] *
!                 |---  [rad_albedo_aerosol] *
!                 |
!                 |---  [rad_cloud_interp] * - [rad_cloud_interp_sub] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine rad_setup
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#ifndef NIM
   subroutine rad_initialize
#else
   subroutine rad_initialize(swhr_int,lwhr_int,nvl,si1d,sl1d)
#endif
!-------------------------------------------------------------------------------
   use varsfc, only    : lalbd_
   use paramodel, only : ncloud_,ngases_,latg_
#ifdef NIM
   use paramodel, only : levs_,levp1_,itsbeg,LONF2S,LATG2S
#endif
   use comsfc
#ifdef DFS
   use dfsvar,only     : si=>sigmafull,sl=>sigma,iope
#endif
#include "abort.h"
#ifndef RMP
   use comfgrid
   use comfphys
   use comfver
   use comio
   use comgrad
   use radiag
#else
   use paramodel, only : IGRD12S,JGRD12S
   use rdparm
   use rscomf_rerun
   use rscomltb
   use rscommap
   use rscomgrad
#endif
#ifdef MP
   use commpi
#endif
!-------------------------------------------------------------------------------
#ifdef NIM
   real*8,intent(in)    ::  si1d(levp1_),sl1d(levs_)
   real*8,intent(in)    ::  swhr_int,lwhr_int
#endif
   integer              ::  ifunit
   logical,save         ::  rad1st
   data rad1st/.true./
!-------------------------------------------------------------------------------
   if(rad1st) then
#ifdef NIM
     dtswav=swhr_int
     dtlwav=lwhr_int    ! org : dtlwav=3.
     cvmint = dtswav
     dtcvav =  min (cvmint,  max (dtswav,dtlwav))
     hprime(:,:,:)=0.   ! org : read(24) hprime
#ifdef LANDSEA
     hprime(:,:,:)=hprime2d(:,:,:)
#endif
#endif
#if defined(RMP) || defined(NIM)
     rrs2=1.0
#endif
     runrad = .true.
#ifndef SWRMDC
     nalaer=49
     kalb=0
     istrat=1
     ko3       = 1 ! nasa o3
     ibnd      = 1 ! =1:use one nir band, =2:use three nir bands
     iswsrc(1) = 1 ! aerosol
     iswsrc(2) = 0 ! o2
     iswsrc(3) = 0 ! co2
     iswsrc(4) = 1 ! water vapor
     iswsrc(5) = 1 ! o3
#ifdef ICECLOUD
     icfc   = 1
     icwp   = 0
     if (ncloud_.ge. 1) icwp = 1
#else
     icfc = 0
     icwp = 0
#endif
#else		/* SWRMDC */
     nfalb =49
     nfaer =47
     kalb=1
     istrat=1
     ko3       = 1 ! nasa o3
     ibnd      = 2 ! =1:use one nir band, =2:use three nir bands
     iswsrc(1) = 1 ! aerosol
     iswsrc(2) = 1 ! o2
     iswsrc(3) = 1 ! co2
     iswsrc(4) = 1 ! water vapor
     iswsrc(5) = 1 ! o3
#ifdef ICECLOUD
     icfc   = 1
     icwp   = 0
     if (ncloud_.ge. 1) icwp = 1
#else
     icfc = 0
     icwp = 0
#endif
#endif
     if(lalbd_.gt.1) kalb = 1
#ifdef AQUA_PLANET
     istrat=0
     ibnd      = 1 ! =1:use one nir band, =2:use three nir bands
     iswsrc(1) = 0 ! aerosol
     iswsrc(2) = 1 ! o2
     iswsrc(3) = 1 ! co2
     iswsrc(4) = 0 ! water vapor
     iswsrc(5) = 1 ! o3
#endif
#ifdef RMP
!
     do j = 1,JGRD12S
       do i = 1,IGRD12S
         sinlar(i,j)=sin(rlat(i,j))
         coslar(i,j)=sqrt(1.e0 - sinlar(i,j)*sinlar(i,j))
       enddo
     enddo
#endif
!
#ifdef NIM
     do k = 1,levs_
       sl(k)=sl1d(k)
     enddo
     do k = 1,levp1_
       si(k)=si1d(k)
     enddo
#endif
     iunco2 = 15
     call rad_cloudiness_init(si)
!
!yh95.. add new data initialization routines for rad and aerosols
!
     call rad_ozone_co2_init(sl,iunco2)
#ifndef SWRMDC
#ifndef ICECLOUD
     call rad_aeros_init(si,sl)
#else
     call rad_aeros_init_ice(si,sl)
#endif
#else
     call rad_aeros_init_nasa(si,sl)
#endif
!
!   specify the latitude where permanent snow resides poleward
!     jsno=latitude closest to 70 deg n-indicating extent of perm
!     snow cover. will change as a fcn. of latitude structure
!
#if defined (RMP) || defined (NIM)
     jsno=0
#else
#ifdef SMP
     jsno=(180+1)/9
#else
     jsno=(latg_+1)/9
#endif
#endif			/* RMP */
!
!  get interval (hrs) between short-wave radiation calls.....
!
     raddt = 3600. * dtswav
     dtlw = 3600. * dtlwav
!
!  paerf when SWRMDC =.false.
!  alvsf,alnsf,alvwf,alnwf from the rad_albedo_aerosol.snl will be discarded
!
     rad1st=.false.
   endif
!
   return
   end subroutine rad_initialize
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_cloud_read(rhcl,ier)
!-------------------------------------------------------------------------------
!
!..  cld-rh relations obtained from mitchell-hahn procedure, here read
!     cld/rh tuning tables for day 0,1,...,5 and merge into 1 file..
!                         .............k.a.c.   mar 93
!     use only one table (day 1) for all fcst hrs....k.a.c. feb 94
!...    4 cld types .... kac  feb96
!    output:
!        rhcl - tuning tables for all forecast days
!        ier  - =1 if tables available.. =-1 if no tables
!
!     this routine is called only by master in mp mode
!
!-------------------------------------------------------------------------------
#ifdef DFS
   use dfsvar, only : iope
#else
   use comio, only  : iope
#endif
!-------------------------------------------------------------------------------
#include "abort.h"
   integer,parameter    ::  mcld=3,nseal=2,ida=1
   integer,parameter    ::  nbin=100,nlon=2,nlat=4
   real                 ::  rhfd(nbin,nlon,nlat,mcld,nseal)
   real                 ::  rrhfd(nbin,nlon,nlat,mcld,nseal)
   real                 ::  rtnffd(nbin,nlon,nlat,mcld,nseal)
   real                 ::  rrnffd(nbin,nlon,nlat,mcld,nseal)
   real                 ::  rhcf(nbin,nlon,nlat,mcld,nseal)
   real                 ::  rtnfcf(nbin,nlon,nlat,mcld,nseal)
   integer              ::  kpts(nlon,nlat,mcld,nseal)
   integer              ::  kkpts(nlon,nlat,mcld,nseal)
   real                 ::  rhc(nlon,nlat,mcld,nseal)
   real                 ::  rhcl (nbin,nlon,nlat,mcld,nseal,ida)
   real                 ::  rhcla(nbin,nlon,nlat,mcld,nseal)
   integer              ::  icdays(15),idate(4)
   integer              ::  ioerr,k
!
   ier = 1
   do itim=1,ida
     icfq = 43 + itim-1
     rewind icfq
!mcl3       nclds=1,2,3 (l,m,h)..jsl=1,2 (land,sea)
!mcl4       mcld=1,2,3,4 (bl,l,m,h)
     binscl = 1./nbin
     do m=1,nseal
       do l=1,mcld
         do k=1,nlat
           do j=1,nlon
             do i=1,nbin
               rrhfd(i,j,k,l,m) = 0.
               rrnffd(i,j,k,l,m) = 0.
             enddo
           enddo
         enddo
       enddo
     enddo
!
     kkpts(:,:,:,:) = 0
!
!....  read the data off the rotating file
!
     read (icfq,err=898,end=899) nbdayi,icdays
#ifdef DBG 
     if (iope) write(6,11)nbdayi
#endif
     go to 123
898  if (iope) write(6,988)itim
     call MPABORT
899  if (iope) write(6,989)itim
     call MPABORT
!
123  continue
!
     do ld = 1,nbdayi
       id = icdays(ld) / 10000
       im = (icdays(ld)-id*10000) / 100
       iy = icdays(ld)-id*10000-im*100
#ifndef NOPRINT
       if (iope) write(6,51)id,im,iy
#endif
     enddo
     read (icfq,err=998,end=999) fhour,idate
#ifdef DBG
     if (iope) write(6,3003)idate,fhour,itim
#endif
!
     go to 223
!
998  if (iope) write(6,988)itim
     call MPABORT
999  if (iope) write(6,989)itim
     call MPABORT
!
223   continue
!
     do kd = 1,nbdayi
       read (icfq) rhfd
       read (icfq) rtnffd
       read (icfq) kpts
       do m = 1,nseal
         do l = 1,mcld
           do k = 1,nlat
             do j = 1,nlon
               do i = 1,nbin
                 rrhfd(i,j,k,l,m) = rrhfd(i,j,k,l,m) + rhfd(i,j,k,l,m)
                 rrnffd(i,j,k,l,m) = rrnffd(i,j,k,l,m)+rtnffd(i,j,k,l,m)
               enddo
             enddo
           enddo
         enddo
       enddo
       do m = 1,nseal
         do l = 1,mcld
           do k = 1,nlat
             do j = 1,nlon
               kkpts(j,k,l,m) = kkpts(j,k,l,m) + kpts(j,k,l,m)
             enddo
           enddo
         enddo
       enddo
!
     enddo
!
     do m = 1,nseal
       do l = 1,mcld
         do k = 1,nlat
           do j = 1,nlon
             do i = 1,nbin
               rhcf(i,j,k,l,m) = rrhfd(i,j,k,l,m)
               rtnfcf(i,j,k,l,m) = rrnffd(i,j,k,l,m)
             enddo
           enddo
         enddo
       enddo
     enddo
     do m = 1,nseal
       do l = 1,mcld
         do k = 1,nlat
           do j = 1,nlon
             kpts(j,k,l,m) = kkpts(j,k,l,m)
           enddo
         enddo
       enddo
     enddo
!
!.....  compute the cumulative frequency distribution..
!
     do n = 1,nseal
       do k = 1,mcld
         do l = 1,nlat
           do j = 1,nlon
             do i = 2,nbin
               rhcf(i,j,l,k,n) = rhcf(i-1,j,l,k,n) + rhcf(i,j,l,k,n)
               rtnfcf(i,j,l,k,n)=rtnfcf(i-1,j,l,k,n) + rtnfcf(i,j,l,k,n)
             enddo
           enddo
         enddo
       enddo
     enddo
     do n = 1,nseal
       do l = 1,nlat
         do j = 1,nlon
           do k = 1,mcld
             do i = 1,nbin
               if (kpts(j,l,k,n).gt.0) then
                 rhcf(i,j,l,k,n) = rhcf(i,j,l,k,n) / kpts(j,l,k,n)
                 rtnfcf(i,j,l,k,n) = rtnfcf(i,j,l,k,n) / kpts(j,l,k,n)
!
!...  cause we mix calculations of rh retune with cray and ibm words
!      the last value of rhcf is close to but ne 1.0,
!      so we reset it in order that the 360 loop gives compleat tabl
!...  rtnfcf caused couple of problems, seems to be ever so slightly
!      gt 1.0
!
                 if (i.eq.nbin) then
                   rhcf(i,j,l,k,n) = 1.0
                 endif
                 if (rtnfcf(i,j,l,k,n).ge.1.0) then
                   rtnfcf(i,j,l,k,n) = 1.0
                 endif
               else
                 rhcf(i,j,l,k,n) = -0.1
                 rtnfcf(i,j,l,k,n) = -0.1
               endif
             enddo
           enddo
         enddo
       enddo
     enddo
#ifdef DBG
     do nsl = 1,nseal
       do kcl = 1,mcld
         if (iope) write(6,264)kcl,nsl
         if (iope) write(6,265)((kpts(i,l,kcl,nsl),i=1,nlon),l=1,nlat)
       enddo
     enddo
#endif
     do nsl = 1,nseal
       do k = 1,mcld
         do l = 1,nlat
           do j = 1,nlon
             if (kpts(j,l,k,nsl).le.0) go to 317
             do i = 1,nbin
               icrit = i
               if (rhcf(i,j,l,k,nsl).ge.rtnfcf(1,j,l,k,nsl)) go to 350
             enddo
!
!... no critical rh
!
317          icrit=-1
!
#ifdef DBG
             if (iope) write(6,210)l,j,nsl
             if (iope) write(6,202)
             do i = 1,nbin
               if (iope) write(6,203)rhcf(i,j,l,k,nsl),rtnfcf(i,j,l,k,nsl)
             enddo
#endif
!
350          rhc(j,l,k,nsl) = icrit * binscl
!
           enddo
         enddo
       enddo
     enddo
!
#ifdef DBG
     do nsl = 1,nseal
       do k = 1,mcld
         if (iope) write(6,1221)k,nsl
         do l = 1,nlat
           if (iope) write(6,211)(rhc(j,l,k,nsl),j=1,nlon)
         enddo
       enddo
     enddo
#endif
     do nsl = 1,nseal
       do ken = 1,mcld
         do l = 1,nlat
           do jl = 1,nlon
             do i = 1,nbin
               rhcl(i,jl,l,ken,nsl,itim) = -0.1
             enddo
           enddo
         enddo
       enddo
     enddo
     do 751 nsl = 1,nseal
     do 751 ken = 1,mcld
     do 751 l = 1,nlat
     do 751 jl = 1,nlon
       if (kpts(jl,l,ken,nsl).le.0) go to 751
       do 753 i = 1,nbin
         do 755 j = 1,nbin
           if (rhcf(j,jl,l,ken,nsl).ge.rtnfcf(i,jl,l,ken,nsl)) then
             rhcl(i,jl,l,ken,nsl,itim) = j*binscl
!
             go to 753
!
           endif
755      continue
753    continue
751  continue
     do lon = 1,nlon
       do lat = 1,nlat
         do nc = 1,mcld
           do nsl = 1,nseal
             isat = 0
             do it = 1,nbin
               cfrac = binscl * (it-1)
               if (rhcl(it,lon,lat,nc,nsl,itim).lt.0.) then
                 if (iope) write(6,1941)it,nsl,nc,lat,lon
               endif
               if (it.lt.nbin.and.rtnfcf(it,lon,lat,nc,nsl).ge.1.) then
                 if (isat.le.0) then
                   isat = it
                   rhsat = rhcl(it,lon,lat,nc,nsl,itim)
                   clsat = cfrac
                 endif
                 rhcl(it,lon,lat,nc,nsl,itim) =                                &
                           rhsat + (1.-rhsat)*(cfrac-clsat)/(1.-clsat)
               endif
               if (it.eq.nbin) rhcl(it,lon,lat,nc,nsl,itim) = 1.
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
!
   do ken = 1,ida
     icfq = 42 + ken
     rewind icfq
   enddo
!
   return
  11 format(1x,' days on file =',i5)
  51 format(1x,' tuned cloud archive data from dd,mm,yy = ',3i4)
 202 format(1x,' model rh ',' obs rtcld')
 203 format(2f10.2)
 210 format(1x,' no crit rh for lat=',i3,' and lon band = ',i3,                &
              ' land(=1) sea=',i3)
 211 format(1x,15f6.2)
 264 format(1x,' number of gg points used in each area..by latitude',          &
              '..for cloud type=',i4,'sealand = ',i2)
 265 format(1x,15i8)
 988 format(1x,'....error reading tables for time = ',i4)
 989 format(1x,'....e.o.f reading tables for time = ',i4)
1221 format(1x,' critical rh for lon,lat arrays for cld type = ',i3,           &
              ' land(=1) sea = ',i3)
1941 format(1x,' neg rhcl for it,nsl,nc,lat,lon = ',5i4,'...stoppp..')
3003 format(5x,'...last date/time and current itim',/,10x,                     &
            4i15,f7.1,i6)
!
   end subroutine rad_cloud_read
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_cloudiness_init (si)
!-------------------------------------------------------------------------------
   use paramodel, only : kdim=>levs_, kdimp=>levp1_
#ifdef DFS
   use dfsvar, only    : iope
#endif
   use constant, only  : cp_,rd_
   use comio
   use comcd1
!-------------------------------------------------------------------------------
!
! --- pressure limits for sfc and top of each cloud domain (l,m,h)
!     in mb, model layers for cld tops are l=7,m=11,h=15 at low
!     latitudes and l= ,m= ,h=  , at pole region.
!
!-------------------------------------------------------------------------------
   real                 ::  si(kdimp), ppptop(4,2)
!
!....     ptop above h changed from 150 to 100, cause
!     data ppptop /1050.,642.,350.,150., 1050.,750.,500.,150./
!          code was truncating tops of convective clouds
!
   data ppptop /1050.,642.,350.,100., 1050.,750.,500.,100./
!
   rocp = rd_ / cp_
!
! --- inverson type cld critical value-istrat=0
!
   clapse = -0.06e0
!
! --- inverson type cld critical value-istrat=1
!
   clapkc = -0.05e0
!
!....critical dtheta/dp for ocean stratus(wgt varies 0 to 1
!                linearly from clapse to clpse)
!
   dclps = -0.01e0
   clpse = clapkc + dclps
   cvtop = 400.0e0
   pstrt = 800.0e0
!
! --- low cld bottom (at sigma=0.95) and top sigma level
!
   do k = 1,kdim
     kk=k
     if (si(kk) .le. 0.95e0) exit
   enddo
   klowb = kk - 1
   silow = ppptop(2,1) * 1.0e-3
   do k = 1,kdim
     kk=k
     if (si(kk) .lt. silow) exit
   enddo
   klowt = kk
!
! --- presure limit at sfc and at top of cloud domains (l,m,h) in mb
!
   do j = 1,2
     do i = 1,4
       ptopc(i,j) = ppptop(i,j)
     enddo
   enddo
!
! --- l cld vertical vel adj boundaries
!
!     vvcld(1) =  0.0003e0
!     vvcld(2) = -0.0005e0
!
!  changed by mk 5/5/98
!
!     vvcld(1) =  0.0006
!     vvcld(2) =  0.0006
!  turned off vertical motion check
!
   vvcld(1) =  100.
   vvcld(2) =  100.
#ifdef vvadj
!
!  turned on vertical motion check (suhee, 2002)
!
   vvcld(1) =  0.0003e0
   vvcld(2) = -0.0005e0
#endif
!
   crhrh = 0.60e0
!
!--- compute llyr--which is topmost non cld(low) layer, for stratiform
!
   xthk = 0.e0
!
!....   default llyr
!
   kl = kdimp
!
!....   topmost noncloud layer will be the one at or above lowest
!         0.1 of the atmosphere..
!
   do k = 1,kdim
     kl = k
     if (si(k).lt.0.9e0) exit
   enddo
!
   llyr = kl-1
#ifndef NOPRINT
   if(iope) write(6,205) llyr,klowb
205 format(1h ,'-------llyr,klowb =',2i5)
#endif
!
   return
   end subroutine rad_cloudiness_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_ozone_co2_init(sigl,nfile)
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm
   use constant, only : pi_,degrad_,daysec_,sbc_
   use comfcst, only : nl, nlp1,                                               &
                       degrad, hsigma, daysec, rco2,                           &
                       dduo3n, ddo3n2, ddo3n3, ddo3n4,                         &
                       xduo3n, xdo3n2, xdo3n3, xdo3n4, prgfdl
!
   real                 ::  sigl(l),pstd(nlp1)
   real, pointer        ::  rad1(:,:), rad2(:,:), rad3(:,:), rad4(:,:)
   real, pointer        ::  xrad1(:,:),xrad2(:,:),xrad3(:,:),xrad4(:,:)
!-------------------------------------------------------------------------------
   rad1=>dduo3n(:,:)
   rad2=>ddo3n2(:,:)
   rad3=>ddo3n3(:,:)
   rad4=>ddo3n4(:,:)
   xrad1=>xduo3n(:,:)
   xrad2=>xdo3n2(:,:)
   xrad3=>xdo3n3(:,:)
   xrad4=>xdo3n4(:,:)
!
!   *******************************************************!
!   *      one time computation of necessary quantities    !
!   *******************************************************!
!===> ... initialize arrays,get constants,etc...
!
   degrad=degrad_
   hsigma=sbc_*1.0e3
   daysec=daysec_
!
!===> ... atmosperic carbon dioxide concentration is now read by rad_co2_read
!         but it defaults to 330 ppm for backward compatibility.
!
   rco2=3.3e-4
!
!===> ... interpolate climo o3 to the current vertical coordinate
!         need layer sigma, get from psfc and layer p for i=1
!
   call rad_ozone_interp_step1(dduo3n,ddo3n2,ddo3n3,ddo3n4,sigl)
!
!===> ... compute detailed o3 profile from the original gfdl pressures
!         where output from rad_ozone_interp_step1 (pstd) is top down in mb*1.e3
!         and psfc=1013.25 mb    ......k.a.c. dec94
!
   call rad_ozone_interp_step2(xduo3n,xdo3n2,xdo3n3,xdo3n4,pstd)
!
!===> ... average climatological valus of o3 from 5 deg lat means,
!         so that time and space interpolation will work
!
   do k = 1,l
     do i = 1,37
       avg=.25e0*(rad1(i,k)+rad2(i,k)+rad3(i,k)+rad4(i,k))
       a1=.5e0*(rad2(i,k)-rad4(i,k))
       b1=.5e0*(rad1(i,k)-rad3(i,k))
       b2=.25e0*((rad1(i,k)+rad3(i,k))-(rad2(i,k)+rad4(i,k)))
       rad1(i,k)=avg
       rad2(i,k)=a1
       rad3(i,k)=b1
       rad4(i,k)=b2
     enddo
   enddo
!
!===> ... average climatological valus of o3 from 5 deg lat means,
!         so that time and space interpolation will work
!        (see subprogram rad_ozone_gfdl)
!
   do k = 1,nl
     do i = 1,37
       avg=.25e0*(xrad1(i,k)+xrad2(i,k)+xrad3(i,k)+xrad4(i,k))
       a1=.5e0*(xrad2(i,k)-xrad4(i,k))
       b1=.5e0*(xrad1(i,k)-xrad3(i,k))
       b2=.25e0*((xrad1(i,k)+xrad3(i,k))-(xrad2(i,k)+xrad4(i,k)))
       xrad1(i,k)=avg
       xrad2(i,k)=a1
       xrad3(i,k)=b1
       xrad4(i,k)=b2
     enddo
   enddo
!
!===> ... get gfdl pressure in cb (flip vertical coordinate)
!
   do n = 1,nl
     prgfdl(n) = pstd(nl+1-n)*1.e-4
   enddo
!
   return
   end subroutine rad_ozone_co2_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_constant_init
!-------------------------------------------------------------------------------
!                                                                               
!     subroutine rad_constant_init defines variables to represent floating- 
!       point constants.                                                        
!                                                                               
!     comdeck hcon contains the common block for these floating-                
!       point constants.                                                        
!                                                                               
!     the naming conventions for the floating-point variables are               
!       as follows:                                                             
!                                                                               
!   1) physical and mathematical constants will be given names                  
!     relevant to their meaning                                                 
!   2) other constants will be given names relevant to their value              
!      and adhering to the following conventions:                               
!       a) the first letter will be represented with an 'h' except              
!          for i) and j) below                                                  
!       b) a decimal point will be represented with a 'p'                       
!       c) there will be no embedded '0'(zero); all 0s will                    
!          be represented with a 'z'                                            
!       d) a minus sign will be represented with an 'm'                         
!       e) the decimal point is assumed after the first digit for               
!          numbers with exponents                                               
!       f) positive exponents are indicated with 'e';negative                   
!          exponents with 'm'                                                   
!       g) digits are truncated in order to have no more than 8                 
!          characters per name                                                  
!       h) numbers less than 0.1 and greater than 10. will be                   
!          represented in exponent format (except a few special cases)          
!       i) the whole numbers from 0.0 through 10.,and 20.,30.,40.,50.,          
!          60.,70.,80.,90.,100.,will be spelled out                             
!       j) good judgment will prevail over all conventions                      
!                                                                               
!       examples                                                                
!     constant           variable name             convention                   
!      600.                 lheatc                  1)                          
!      680.                 lheats                  1)                          
!     1.4142               sqroot2                  1)                          
!     2.0                    two                    2)-(i)                      
!    -3.0                  hm3pz                    2)-(a,b,d)                  
!    310.                  c31e2                    2)-(a,e,f,h)                
!   -0.7239e-9             hm723m1z                 2)-(a,c,d,e,f,g,h)          
!     0.0                   zero                    2)-(i)                      
!     0.1                   hp1                     2)-(a,b,h)                  
!     0.01                 h1m2                     2)-(a,e,f,h)                
!     30.                  thirty                   2)-(h,i)                    
!     0.5                  haf                      2)-(j)                      
!     9.0                  hnine                    2)-(j)                      
!                                                                               
!-------------------------------------------------------------------------------
   use hcon
!-------------------------------------------------------------------------------
!
!******the following are physical constants*****                                
!        arranged in alphabetical order                                         
!
   amolwt=28.9644                                                            
   csubp=1.00484e7                                                           
   diffctr=1.66                                                              
   g=980.665                                                                 
   ginv=1./g                                                                 
   gravdr=980.0                                                              
   o3difctr=1.90                                                             
   p0=1013250.                                                               
   p0inv=1./p0                                                               
   gp0inv=ginv*p0inv                                                         
   p0xzp2=202649.902                                                         
   p0xzp8=810600.098                                                         
   p0x2=2.*1013250.                                                          
   radcon=8.427                                                              
   radcon1=1./8.427                                                          
   ratco2mw=1.519449738                                                      
   rath2omw=.622                                                             
   rgas=8.3142e7                                                             
   rgassp=8.31432e7                                                          
   secpda=8.64e4                                                             
!                                                                               
!******the following are mathematical constants*******                          
!        arranged in decreasing order                                           
!
   hundred=100.                                                              
   hninety=90.                                                               
   sixty=60.                                                                 
   fifty=50.                                                                 
   ten=10.                                                                   
   eight=8.                                                                  
   five=5.                                                                   
   four=4.                                                                   
   three=3.                                                                  
   two=2.                                                                    
   one=1.                                                                    
   haf=0.5                                                                   
   quartr=0.25                                                               
   zero=0.                                                                   
!                                                                               
!******following are positive floating point constants(hs)                     
!       arranged in decreasing order                                            
!
   h83e26=8.3e26                                                             
   h71e26=7.1e26                                                             
   h1e15=1.e15                                                               
   h1e13=1.e13                                                               
   h1e11=1.e11                                                               
   h1e8=1.e8                                                                 
   h2e6=2.0e6                                                                
   h1e6=1.0e6                                                                
   h69766e5=6.97667e5                                                        
   h4e5=4.e5                                                                 
   h165e5=1.65e5                                                             
   h5725e4=57250.                                                            
   h488e4=48800.                                                             
   h1e4=1.e4                                                                 
   h24e3=2400.                                                               
   h20788e3=2078.8                                                           
   h2075e3=2075.                                                             
   h18e3=1800.                                                               
   h1224e3=1224.                                                             
   h67390e2=673.9057                                                         
   h5e2=500.                                                                 
   h3082e2=308.2                                                             
   h3e2=300.                                                                 
   h2945e2=294.5                                                             
   h29316e2=293.16                                                           
   h26e2=260.0                                                               
   h25e2=250.                                                                
   h23e2=230.                                                                
   h2e2=200.0                                                                
   h15e2=150.                                                                
   h1386e2=138.6                                                             
   h1036e2=103.6                                                             
   h8121e1=81.21                                                             
   h35e1=35.                                                                 
   h3116e1=31.16                                                             
   h28e1=28.                                                                 
   h181e1=18.1                                                               
   h18e1=18.                                                                 
   h161e1=16.1                                                               
   h16e1=16.                                                                 
   h1226e1=12.26                                                             
   h9p94=9.94                                                                
   h6p08108=6.081081081                                                      
   h3p6=3.6                                                                  
   h3p5=3.5                                                                  
   h2p9=2.9                                                                  
   h2p8=2.8                                                                  
   h2p5=2.5                                                                  
   h1p8=1.8                                                                  
   h1p4387=1.4387                                                            
   h1p41819=1.418191                                                         
   h1p4=1.4                                                                  
   h1p25892=1.258925411                                                      
   h1p082=1.082                                                              
   hp816=0.816                                                               
   hp805=0.805                                                               
   hp8=0.8                                                                   
   hp60241=0.60241                                                           
   hp602409=0.60240964                                                       
   hp6=0.6                                                                   
   hp526315=0.52631579                                                       
   hp518=0.518                                                               
   hp5048=0.5048                                                             
   hp3795=0.3795                                                             
   hp369=0.369                                                               
   hp26=0.26                                                                 
   hp228=0.228                                                               
   hp219=0.219                                                               
   hp166666=.166666                                                          
   hp144=0.144                                                               
   hp118666=0.118666192                                                      
   hp1=0.1                                                                   
!
!        (negative exponentials begin here)                                     
!
   h658m2=0.0658                                                             
   h625m2=0.0625                                                             
   h44871m2=4.4871e-2                                                        
   h44194m2=.044194                                                          
   h42m2=0.042                                                               
   h41666m2=0.0416666                                                        
   h28571m2=.02857142857                                                     
   h2118m2=0.02118                                                           
   h129m2=0.0129                                                             
   h1m2=.01                                                                  
   h559m3=5.59e-3                                                            
   h3m3=0.003                                                                
   h235m3=2.35e-3                                                            
   h1m3=1.0e-3                                                               
   h987m4=9.87e-4                                                            
   h323m4=0.000323                                                           
   h3m4=0.0003                                                               
   h285m4=2.85e-4                                                            
   h1m4=0.0001                                                               
   h75826m4=7.58265e-4                                                       
   h6938m5=6.938e-5                                                          
   h394m5=3.94e-5                                                            
   h37412m5=3.7412e-5                                                        
   h15m5=1.5e-5                                                              
   h1439m5=1.439e-5                                                          
   h128m5=1.28e-5                                                            
   h102m5=1.02e-5                                                            
   h1m5=1.0e-5                                                               
   h7m6=7.e-6                                                                
   h4999m6=4.999e-6                                                          
   h451m6=4.51e-6                                                            
   h25452m6=2.5452e-6                                                        
   h1m6=1.e-6                                                                
   h391m7=3.91e-7                                                            
   h1174m7=1.174e-7                                                          
   h8725m8=8.725e-8                                                          
   h327m8=3.27e-8                                                            
   h257m8=2.57e-8                                                            
   h1m8=1.0e-8                                                               
   h23m10=2.3e-10                                                            
   h14m10=1.4e-10                                                            
   h11m10=1.1e-10                                                            
   h1m10=1.e-10                                                              
   h83m11=8.3e-11                                                            
   h82m11=8.2e-11                                                            
   h8m11=8.e-11                                                              
   h77m11=7.7e-11                                                            
   h72m11=7.2e-11                                                            
   h53m11=5.3e-11                                                            
   h48m11=4.8e-11                                                            
   h44m11=4.4e-11                                                            
   h42m11=4.2e-11                                                            
   h37m11=3.7e-11                                                            
   h35m11=3.5e-11                                                            
   h32m11=3.2e-11                                                            
   h3m11=3.0e-11                                                             
   h28m11=2.8e-11                                                            
   h24m11=2.4e-11                                                            
   h23m11=2.3e-11                                                            
   h2m11=2.e-11                                                              
   h18m11=1.8e-11                                                            
   h15m11=1.5e-11                                                            
   h14m11=1.4e-11                                                            
   h114m11=1.14e-11                                                          
   h11m11=1.1e-11                                                            
   h1m11=1.e-11                                                              
   h96m12=9.6e-12                                                            
   h93m12=9.3e-12                                                            
   h77m12=7.7e-12                                                            
   h74m12=7.4e-12                                                            
   h65m12=6.5e-12                                                            
   h62m12=6.2e-12                                                            
   h6m12=6.e-12                                                              
   h45m12=4.5e-12                                                            
   h44m12=4.4e-12                                                            
   h4m12=4.e-12                                                              
   h38m12=3.8e-12                                                            
   h37m12=3.7e-12                                                            
   h3m12=3.e-12                                                              
   h29m12=2.9e-12                                                            
   h28m12=2.8e-12                                                            
   h24m12=2.4e-12                                                            
   h21m12=2.1e-12                                                            
   h16m12=1.6e-12                                                            
   h14m12=1.4e-12                                                            
   h12m12=1.2e-12                                                            
   h8m13=8.e-13                                                              
   h46m13=4.6e-13                                                            
   h36m13=3.6e-13                                                            
   h135m13=1.35e-13                                                          
   h12m13=1.2e-13                                                            
   h1m13=1.e-13                                                              
   h3m14=3.e-14                                                              
   h15m14=1.5e-14                                                            
   h14m14=1.4e-14                                                            
   h101m16=1.01e-16                                                          
   h1m16=1.0e-16                                                             
   h1m17=1.e-17                                                              
   h1m18=1.e-18                                                              
   h1m19=1.e-19                                                              
   h1m20=1.e-20                                                              
   h1m21=1.e-21                                                              
   h1m22=1.e-22                                                              
   h1m23=1.e-23                                                              
   h1m24=1.e-24                                                              
   h26m30=2.6e-30                                                            
   h14m30=1.4e-30                                                            
   h25m31=2.5e-31                                                            
   h21m31=2.1e-31                                                            
   h12m31=1.2e-31                                                            
   h9m32=9.e-32                                                              
   h55m32=5.5e-32                                                            
   h45m32=4.5e-32                                                            
   h4m33=4.e-33                                                              
   h62m34=6.2e-34                                                            
!                                                                               
!******following are negative floating point constants (hms)                   
!          arranged in descending order                                         
   hm2m2=-.02                                                                
   hm6666m2=-.066667                                                         
   hmp5=-0.5                                                                 
   hmp575=-0.575                                                             
   hmp66667=-.66667                                                          
   hmp805=-0.805                                                             
   hm1ez=-1.                                                                 
   hm13ez=-1.3                                                               
   hm19ez=-1.9                                                               
   hm1e1=-10.                                                                
   hm1597e1=-15.97469413                                                     
   hm161e1=-16.1                                                             
   hm1797e1=-17.97469413                                                     
   hm181e1=-18.1                                                             
   hm8e1=-80.                                                                
   hm1e2=-100.                                                               
!                                                                               
   return                                                                    
   end subroutine rad_constant_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_ozone_interp_step1(aa,bb,cc,dd,sigl)
!-------------------------------------------------------------------------------
!
! subroutine:    rad_ozone_interp_step1 
!
! abstract: this code written at gfdl...
!  - compute zonal mean ozone for sigma layers
!  - calculates seasonal zonal mean ozone,every 5 deg of latitude,
!    for current model vertical coordinate. output data in g/g * 1.e4
!    code is called only once.
!
! program history log:
!   1984-01-01  fels and schwarzkopf,gfdl.
!   1989-07-07  kenneth campana - adapted stand-alone code for in-line use.
!
! usage:    call rad_ozone_interp_step1(o3,sigl)
!   input argument list:
!     sigl     - layer sigma (k=1 is lowest model layer)
!   output argument list:
!     o3       - zonal mean ozone data in all model layers (g/g*1.e4)
!                dimensioned(l,n,is),where l(=37) is latitude between
!                n and s poles,n=num of vertical lyrs(k=1 is top lyr)
!                ,and is(=4) defines the season-win,spr,sum,fall.
!
!   output files:
!     output   - print file.
!
!-------------------------------------------------------------------------------
!
!....     program rad_ozone_interp_step1 from dan schwarzkopf-gets zonal mean o3
!..       code adapted for mrf use, in-line.....    k.a.c. june 1989
!....     launcher======subroutine rad_ozone_interp_step1(t41,o3o3)
!..    output o3 is winter,spring,summer,fall (northern hemisphere)
!-------------------------------------------------------------------------------
!       save data on permanent data set denoted by co222 ***!
!          ..... k.campana   october 1988
!     dimension t41(lp2,2),o3o3(37,l,4)
!-------------------------------------------------------------------------------
   use comio
#ifdef DFS
   use dfsvar, only : iope
#endif
   use rdparm
!-------------------------------------------------------------------------------
   real                 ::  sigl(l), o3o3(37,l,4)
   real                 ::  aa(37,l), bb(37,l), cc(37,l), dd(37,l)
!
   real                 ::  qi(82)
   real                 ::  dduo3n(19,l),ro31(10,41),ro32(10,41),duo3n(19,41)
   real                 ::  tempn(19)
   real                 ::  o3hi(10,25),o3lo1(10,16),o3lo2(10,16),             &
                            o3lo3(10,16),o3lo4(10,16)
   real                 ::  o3hi1(10,16),o3hi2(10,9),ph1(45),ph2(37),p1(48),   &
                            p2(33)
   real                 ::  o35deg(37,l)
   real                 ::  rstd(81),ro3(10,41),ro3m(10,40),rbar(l),rdata(81)
   real                 ::  phalf(lp1),pstd(lp2),p(81),ph(82)
   equivalence (o3hi1(1,1),o3hi(1,1)),(o3hi2(1,1),o3hi(1,17))
   equivalence (ph1(1),ph(1)),(ph2(1),ph(46))
   equivalence (p1(1),p(1)),(p2(1),p(49))
   data ph1/      0.,                                                          &
        0.1027246e-04, 0.1239831e-04, 0.1491845e-04, 0.1788053e-04,            &
        0.2135032e-04, 0.2540162e-04, 0.3011718e-04, 0.3558949e-04,            &
        0.4192172e-04, 0.4922875e-04, 0.5763817e-04, 0.6729146e-04,            &
        0.7834518e-04, 0.9097232e-04, 0.1053635e-03, 0.1217288e-03,            &
        0.1402989e-03, 0.1613270e-03, 0.1850904e-03, 0.2119495e-03,            &
        0.2423836e-03, 0.2768980e-03, 0.3160017e-03, 0.3602623e-03,            &
        0.4103126e-03, 0.4668569e-03, 0.5306792e-03, 0.6026516e-03,            &
        0.6839018e-03, 0.7759249e-03, 0.8803303e-03, 0.9987843e-03,            &
        0.1133178e-02, 0.1285955e-02, 0.1460360e-02, 0.1660001e-02,            &
        0.1888764e-02, 0.2151165e-02, 0.2452466e-02, 0.2798806e-02,            &
        0.3197345e-02, 0.3656456e-02, 0.4185934e-02, 0.4797257e-02/
   data ph2/                                                                   &
        0.5503893e-02, 0.6321654e-02, 0.7269144e-02, 0.8368272e-02,            &
        0.9644873e-02, 0.1112946e-01, 0.1285810e-01, 0.1487354e-01,            &
        0.1722643e-01, 0.1997696e-01, 0.2319670e-01, 0.2697093e-01,            &
        0.3140135e-01, 0.3660952e-01, 0.4274090e-01, 0.4996992e-01,            &
        0.5848471e-01, 0.6847525e-01, 0.8017242e-01, 0.9386772e-01,            &
        0.1099026e000, 0.1286765e000, 0.1506574e000, 0.1763932e000,            &
        0.2065253e000, 0.2415209e000, 0.2814823e000, 0.3266369e000,            &
        0.3774861e000, 0.4345638e000, 0.4984375e000, 0.5697097e000,            &
        0.6490189e000, 0.7370409e000, 0.8344896e000, 0.9421190e000,            &
        0.1000000e001/
   data p1/                                                                    &
        0.9300000e-05, 0.1129521e-04, 0.1360915e-04, 0.1635370e-04,            &
        0.1954990e-04, 0.2331653e-04, 0.2767314e-04, 0.3277707e-04,            &
        0.3864321e-04, 0.4547839e-04, 0.5328839e-04, 0.6234301e-04,            &
        0.7263268e-04, 0.8450696e-04, 0.9793231e-04, 0.1133587e-03,            &
        0.1307170e-03, 0.1505832e-03, 0.1728373e-03, 0.1982122e-03,            &
        0.2266389e-03, 0.2592220e-03, 0.2957792e-03, 0.3376068e-03,            &
        0.3844381e-03, 0.4379281e-03, 0.4976965e-03, 0.5658476e-03,            &
        0.6418494e-03, 0.7287094e-03, 0.8261995e-03, 0.9380076e-03,            &
        0.1063498e-02, 0.1207423e-02, 0.1369594e-02, 0.1557141e-02,            &
        0.1769657e-02, 0.2015887e-02, 0.2295520e-02, 0.2620143e-02,            &
        0.2989651e-02, 0.3419469e-02, 0.3909867e-02, 0.4481491e-02,            &
        0.5135272e-02, 0.5898971e-02, 0.6774619e-02, 0.7799763e-02/
   data p2/                                                                    &
        0.8978218e-02, 0.1036103e-01, 0.1195488e-01, 0.1382957e-01,            &
        0.1599631e-01, 0.1855114e-01, 0.2151235e-01, 0.2501293e-01,            &
        0.2908220e-01, 0.3390544e-01, 0.3952926e-01, 0.4621349e-01,            &
        0.5403168e-01, 0.6330472e-01, 0.7406807e-01, 0.8677983e-01,            &
        0.1015345e000, 0.1189603e000, 0.1391863e000, 0.1630739e000,            &
        0.1908004e000, 0.2235461e000, 0.2609410e000, 0.3036404e000,            &
        0.3513750e000, 0.4055375e000, 0.4656677e000, 0.5335132e000,            &
        0.6083618e000, 0.6923932e000, 0.7845676e000, 0.8875882e000,            &
        0.1000000e001/
   data o3hi1/                                                                 &
       .55,.50,.45,.45,.40,.35,.35,.30,.30,.30,                                &
       .55,.51,.46,.47,.42,.38,.37,.36,.35,.35,                                &
       .55,.53,.48,.49,.44,.42,.41,.40,.38,.38,                                &
       .60,.55,.52,.52,.50,.47,.46,.44,.42,.41,                                &
       .65,.60,.55,.56,.53,.52,.50,.48,.45,.45,                                &
       .75,.65,.60,.60,.55,.55,.55,.50,.48,.47,                                &
       .80,.75,.75,.75,.70,.70,.65,.63,.60,.60,                                &
       .90,.85,.85,.80,.80,.75,.75,.74,.72,.71,                                &
       1.10,1.05,1.00,.90,.90,.90,.85,.83,.80,.80,                             &
       1.40,1.30,1.25,1.25,1.25,1.20,1.15,1.10,1.05,1.00,                      &
       1.7,1.7,1.6,1.6,1.6,1.6,1.6,1.6,1.5,1.5,                                &
       2.1,2.0,1.9,1.9,1.9,1.8,1.8,1.8,1.7,1.7,                                &
       2.4,2.3,2.2,2.2,2.2,2.1,2.1,2.1,2.0,2.0,                                &
       2.7,2.5,2.5,2.5,2.5,2.5,2.4,2.4,2.3,2.3,                                &
       2.9,2.8,2.7,2.7,2.7,2.7,2.7,2.7,2.6,2.6,                                &
       3.1,3.1,3.0,3.0,3.0,3.0,3.0,3.0,2.9,2.8/
   data o3hi2/                                                                 &
       3.3,3.4,3.4,3.6,3.7,3.9,4.0,4.1,4.0,3.8,                                &
       3.6,3.8,3.9,4.2,4.7,5.3,5.6,5.7,5.5,5.2,                                &
       4.1,4.3,4.7,5.2,6.0,6.7,7.0,6.8,6.4,6.2,                                &
       5.4,5.7,6.0,6.6,7.3,8.0,8.4,7.7,7.1,6.7,                                &
       6.7,6.8,7.0,7.6,8.3,10.0,9.6,8.2,7.5,7.2,                               &
       9.2,9.3,9.4,9.6,10.3,10.6,10.0,8.5,7.7,7.3,                             &
       12.6,12.1,12.0,12.1,11.7,11.0,10.0,8.6,7.8,7.4,                         &
       14.2,13.5,13.1,12.8,11.9,10.9,9.8,8.5,7.8,7.5,                          &
       14.3,14.0,13.4,12.7,11.6,10.6,9.3,8.4,7.6,7.3/
   data o3lo1/                                                                 &
       14.9,14.2,13.3,12.5,11.2,10.3,9.5,8.6,7.5,7.4,                          &
       14.5,14.1,13.0,11.8,10.5,9.8,9.2,7.9,7.4,7.4,                           &
       11.8,11.5,10.9,10.5,9.9,9.6,8.9,7.5,7.2,7.2,                            &
       7.3,7.7,7.8,8.4,8.4,8.5,7.9,7.4,7.1,7.1,                                &
       4.1,4.4,5.3,6.6,6.9,7.5,7.4,7.2,7.0,6.9,                                &
       1.8,1.9,2.5,3.3,4.5,5.8,6.3,6.3,6.4,6.1,                                &
       0.4,0.5,0.8,1.2,2.7,3.6,4.6,4.7,5.0,5.2,                                &
       .10,.15,.20,.50,1.4,2.1,3.0,3.2,3.5,3.9,                                &
       .07,.10,.12,.30,1.0,1.4,1.8,1.9,2.3,2.5,                                &
       .06,.08,.10,.15,.60,.80,1.4,1.5,1.5,1.6,                                &
       .05,.05,.06,.09,.20,.40,.70,.80,.90,.90,                                &
       .05,.05,.06,.08,.10,.13,.20,.25,.30,.40,                                &
       .05,.05,.05,.06,.07,.07,.08,.09,.10,.13,                                &
       .05,.05,.05,.05,.06,.06,.06,.06,.07,.07,                                &
       .05,.05,.05,.05,.05,.05,.05,.06,.06,.06,                                &
       .04,.04,.04,.04,.04,.04,.04,.05,.05,.05/
   data o3lo2/                                                                 &
       14.8,14.2,13.8,12.2,11.0,9.8,8.5,7.8,7.4,6.9,                           &
       13.2,13.0,12.5,11.3,10.4,9.0,7.8,7.5,7.0,6.6,                           &
       10.6,10.6,10.7,10.1,9.4,8.6,7.5,7.0,6.5,6.1,                            &
       7.0,7.3,7.5,7.5,7.5,7.3,6.7,6.4,6.0,5.8,                                &
       3.8,4.0,4.7,5.0,5.2,5.9,5.8,5.6,5.5,5.5,                                &
       1.4,1.6,2.4,3.0,3.7,4.1,4.6,4.8,5.1,5.0,                                &
       .40,.50,.90,1.2,2.0,2.7,3.2,3.6,4.3,4.1,                                &
       .07,.10,.20,.30,.80,1.4,2.1,2.4,2.7,3.0,                                &
       .06,.07,.09,.15,.30,.70,1.2,1.4,1.6,2.0,                                &
       .05,.05,.06,.12,.15,.30,.60,.70,.80,.80,                                &
       .04,.05,.06,.08,.09,.15,.30,.40,.40,.40,                                &
       .04,.04,.05,.055,.06,.09,.12,.13,.15,.15,                               &
       .03,.03,.045,.052,.055,.06,.07,.07,.06,.07,                             &
       .03,.03,.04,.051,.052,.052,.06,.06,.05,.05,                             &
       .02,.02,.03,.05,.05,.05,.04,.04,.04,.04,                                &
       .02,.02,.02,.04,.04,.04,.03,.03,.03,.03/
   data o3lo3/                                                                 &
       14.5,14.0,13.5,11.3,11.0,10.0,9.0,8.3,7.5,7.3,                          &
       13.5,13.2,12.5,11.1,10.4,9.7,8.2,7.8,7.4,6.8,                           &
       10.8,10.9,11.0,10.4,10.0,9.6,7.9,7.5,7.0,6.7,                           &
       7.3,7.5,7.8,8.5,9.0,8.5,7.7,7.4,6.9,6.5,                                &
       4.1,4.5,5.3,6.2,7.3,7.7,7.3,7.0,6.6,6.4,                                &
       1.8,2.0,2.2,3.8,4.3,5.6,6.2,6.2,6.4,6.2,                                &
       .30,.50,.60,1.5,2.8,3.7,4.5,4.7,5.5,5.6,                                &
       .09,.10,.15,.60,1.2,2.1,3.0,3.5,4.0,4.3,                                &
       .06,.08,.10,.30,.60,1.1,1.9,2.2,2.9,3.0,                                &
       .04,.05,.06,.15,.45,.60,1.1,1.3,1.6,1.8,                                &
       .04,.04,.04,.08,.20,.30,.55,.60,.75,.90,                                &
       .04,.04,.04,.05,.06,.10,.12,.15,.20,.25,                                &
       .04,.04,.03,.04,.05,.06,.07,.07,.07,.08,                                &
       .03,.03,.04,.05,.05,.05,.05,.05,.05,.05,                                &
       .03,.03,.03,.04,.04,.04,.05,.05,.04,.04,                                &
       .02,.02,.02,.04,.04,.04,.04,.04,.03,.03/
   data o3lo4/                                                                 &
       14.2,13.8,13.2,12.5,11.7,10.5,8.6,7.8,7.5,6.6,                          &
       12.5,12.4,12.2,11.7,10.8,9.8,7.8,7.2,6.5,6.1,                           &
       10.6,10.5,10.4,10.1,9.6,9.0,7.1,6.8,6.1,5.9,                            &
       7.0,7.4,7.9,7.8,7.6,7.3,6.2,6.1,5.8,5.6,                                &
       4.2,4.6,5.1,5.6,5.9,5.9,5.9,5.8,5.6,5.3,                                &
       2.1,2.3,2.6,2.9,3.5,4.3,4.8,4.9,5.1,5.1,                                &
       0.7,0.8,1.0,1.5,2.0,2.8,3.5,3.6,3.7,4.0,                                &
       .15,.20,.40,.50,.60,1.4,2.1,2.2,2.3,2.5,                                &
       .08,.10,.15,.25,.30,.90,1.2,1.3,1.4,1.6,                                &
       .07,.08,.10,.14,.20,.50,.70,.90,.90,.80,                                &
       .05,.06,.08,.12,.14,.20,.35,.40,.60,.50,                                &
       .05,.05,.08,.09,.09,.09,.11,.12,.15,.18,                                &
       .04,.05,.06,.07,.07,.08,.08,.08,.08,.08,                                &
       .04,.04,.05,.07,.07,.07,.07,.07,.06,.05,                                &
       .02,.02,.04,.05,.05,.05,.05,.05,.04,.04,                                &
       .02,.02,.03,.04,.04,.04,.04,.04,.03,.03/
!
!***read in user-specified pressures,in mb. this can be output from ptz
!kac ..  file 41 from ptz(psfc=1013.25 mb)
606   format (5e16.9)
!cc      read (8,606) (pstd(k),k=1,lp2)
!cc      read (8,606) (phalf(k),k=1,lp1)
!o222  **************************************************!
!   
   pss=1013.250e0
   pstd(1) = 0.e0
   pstd(lp2) = pss
   do k = 2,lp1
     pstd(k) = sigl(lp2-k) * pss
   enddo
   phalf(1) = 0.e0
   phalf(lp1) = pss
   do k = 1,lm1
     phalf(k+1) = 0.5e0 * (pstd(k+1)+pstd(k+2))
   enddo
!
!kac  do 300 k=1,lp2
!kac    pstd(k) = t41(k,1)
!k300 continue
!kac  do 301 k=1,lp1
!kac    phalf(k) = t41(k,2)
!k301 continue
!c          rewind 66
!o222  **************************************************!
!
   nkk=41
   nk=81
   nkp=nk+1
   do k = 1,lp1
     phalf(k)=phalf(k)*1.0e03
     pstd(k)=pstd(k+1)*1.0e03
   enddo
   do k = 1,nk
     ph(k)=ph(k)*1013250.
     p(k)=p(k)*1013250.
   enddo
   ph(nkp)=ph(nkp)*1013250.
#ifndef RMP
#ifdef DBG
   if( iope ) write (6,3) (phalf(k),k=1,lp1)
   if( iope ) write (6,3) (pstd(k),k=1,lp1)
#endif
#endif
!
!***load arrays ro31,ro32,as in dicks pgm.
!
   do k = 1,25
     do n = 1,10
       ro31(n,k)=o3hi(n,k)
       ro32(n,k)=o3hi(n,k)
     enddo
   enddo
!
   do ncase = 1,4
     itape=ncase+50
     iplace=2
     if (ncase.eq.2) iplace=4
     if (ncase.eq.3) iplace=1
     if (ncase.eq.4) iplace=3
!
!***ncase=1: spring (in n.h.)
!***ncase=2: fall   (in n.h.)
!***ncase=3: winter (in n.h.)
!***ncase=4: summer (in n.h.)
!
     if (ncase.eq.1.or.ncase.eq.2) then
       do k = 26,41
         do n = 1,10
           ro31(n,k)=o3lo1(n,k-25)
           ro32(n,k)=o3lo2(n,k-25)
         enddo
       enddo
     endif
     if (ncase.eq.3.or.ncase.eq.4) then
       do k = 26,41
         do n = 1,10
           ro31(n,k)=o3lo3(n,k-25)
           ro32(n,k)=o3lo4(n,k-25)
         enddo
       enddo
     endif
     do kk = 1,nkk
       do n = 1,10
         duo3n(n,kk)=ro31(11-n,kk)
         duo3n(n+9,kk)=ro32(n,kk)
       enddo
       duo3n(10,kk)=.5*(ro31(1,kk)+ro32(1,kk))
     enddo
!
!***for ncase=2 or ncase=4,reverse latitude arrangement of corr. season
!
     if (ncase.eq.2.or.ncase.eq.4) then
       do kk = 1,nkk
         do n = 1,19
           tempn(n)=duo3n(20-n,kk)
         enddo
         do n = 1,19
           duo3n(n,kk)=tempn(n)
         enddo
       enddo
     endif
!
!***duo3n now is o3 profile for appropriate season,at std. pressure
!      levels
!***begin latitude (10 deg) loop
!
     do n = 1,19
       do kk = 1,nkk
         rstd(kk)=duo3n(n,kk)
       enddo
       nkm=nk-1
       nkmm=nk-3
!
!  bessels half-point interpolation formula
!
       do k = 4,nkmm,2
         ki=k/2
         rdata(k)=.5*(rstd(ki)+rstd(ki+1))-(rstd(ki+2)-rstd(ki+1)-         &
         rstd(ki)+rstd(ki-1))/16.
       enddo
       rdata(2)=.5*(rstd(2)+rstd(1))
       rdata(nkm)=.5*(rstd(nkk)+rstd(nkk-1))
!
!     put unchanged data into new array
!
       do k = 1,nk,2
         kq=(k+1)/2
         rdata(k)=rstd(kq)
       enddo
!
!---note to nmc: this write is commented out to reduce printout
!     calculate layer-mean ozone mixing ratio for each model level
!
       do 99 kk = 1,l
         rbar(kk)=0.
!     loop to calculate sums to get layer ozone mean
         do 98 k = 1,nk
           if(ph(k+1).lt.phalf(kk)) go to 98
           if(ph(k).gt.phalf(kk+1)) go to 98
           if(ph(k+1).lt.phalf(kk+1).and.ph(k).lt.phalf(kk)) rbar(kk)=         &
                rbar(kk)+rdata(k)*(ph(k+1)-phalf(kk))
           if(ph(k+1).lt.phalf(kk+1).and.ph(k).ge.phalf(kk)) rbar(kk)=         &
                rbar(kk)+rdata(k)*(ph(k+1)-ph(k))
           if(ph(k+1).gt.phalf(kk+1).and.ph(k).gt.phalf(kk)) rbar(kk)=         &
                rbar(kk)+rdata(k)*(phalf(kk+1)-ph(k))
98       continue
         rbar(kk)=rbar(kk)/(phalf(kk+1)-phalf(kk))
         if(rbar(kk).gt..0000) go to 99
!
!     code to cover case when model resolution is so fine that no value
!     of p(k) in the ozone data array falls between phalf(kk+1) and
!     phalf(kk).   procedure is to simply grab the nearest value from
!     rdata
!
         do k = 1,nk
           if(ph(k).lt.phalf(kk).and.ph(k+1).ge.phalf(kk+1)) rbar(kk)=         &
                rdata(k)
         enddo
99     continue
!
!     calculate total ozone
!
       o3rd=0.
       do kk = 1,80
         o3rd=o3rd+rdata(kk)*(ph(kk+1)-ph(kk))
       enddo
       o3rd=o3rd+rdata(81)*(p(81)-ph(81))
       o3rd=o3rd/980.
       o3tot=0.
       do kk = 1,l
         o3tot=o3tot+rbar(kk)*(phalf(kk+1)-phalf(kk))
       enddo
       o3tot=o3tot/980.
!
!     units are micrograms/cm**2
!
       o3du=o3tot/2.144
!
!     o3du units are dobson units (10**-3 atm-cm)
!--note to nmc: this is commented out to save printout
!     write (6,796) o3rd,o3tot,o3du
!
       do kk = 1,l
         dduo3n(n,kk)=rbar(kk)*.01
       enddo
     enddo
!
!***end of latitude loop
!
!***create 5 deg ozone quantities by linear interpolation of
!      10 deg values
!
     do kk = 1,l
       do n = 1,19
         o35deg(2*n-1,kk)=dduo3n(n,kk)
       enddo
       do n = 1,18
         o35deg(2*n,kk)=0.5*(dduo3n(n,kk)+dduo3n(n+1,kk))
       enddo
     enddo
!
!***output to unit (itape) the ozone values for later use
!c          write (66) o35deg
!
     do jj = 1,37
       do ken = 1,l
         o3o3(jj,ken,iplace) = o35deg(jj,ken)
       enddo
     enddo
   enddo
!
!***end of loop over cases
!
   aa(:,:)=o3o3(:,:,1)
   bb(:,:)=o3o3(:,:,2)
   cc(:,:)=o3o3(:,:,3)
   dd(:,:)=o3o3(:,:,4)
!
   return
   1  format(10f4.2)
   2  format(10x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x)
   3  format(10e12.5)
! 796 format(3e12.5)
  797 format(10f7.2)
! 798 format(20f6.2)
  799 format(19f6.4)
  800 format(19f6.2)
  101 format(5x,1h*,f6.5,1h,,f6.5,1h,,f6.5,1h,,f6.5,1h,,f6.5,1h,,f6.5,         &
   1h,,f6.5,1h,,f6.5,1h,,f6.5,1h,)
!
   end subroutine rad_ozone_interp_step1
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_ozone_interp_step2(aa,bb,cc,dd,pstd)
!-------------------------------------------------------------------------------
!
! subroutine:    rad_ozone_interp_step2
!
! abstract: this code written at gfdl...
!  - compute zonal mean ozone for sigma layers
!  - calculates seasonal zonal mean ozone,every 5 deg of latitude,
!    for current model vertical coordinate. output data in g/g * 1.e4
!    code is called only once.
!
! program history log:
!   84-01-01  fels and schwarzkopf,gfdl.
!   89-07-07  kenneth campana - adapted stand-alone code for in-line use.
!
! usage:    call rad_ozone_interp_step2(o3,sigl)
!   input argument list:
!     sigl     - layer sigma (k=1 is lowest model layer)
!   output argument list:
!     o3       - zonal mean ozone data in all model layers (g/g*1.e4)
!                dimensioned(l,n,is),where l(=37) is latitude between
!                n and s poles,n=num of vertical lyrs(k=1 is top lyr)
!                ,and is(=4) defines the season-win,spr,sum,fall.
!
!   output files:
!     output   - print file.
!
!-------------------------------------------------------------------------------
!
!....     program rad_ozone_interp_step2 from dan schwarzkopf-gets zonal mean o3
!..       code adapted for mrf use, in-line.....    k.a.c. june 1989
!....     launcher======subroutine rad_ozone_interp_step2(t41,o3o3)
!..     no input, just calculate complete 81 layer data for later
!..                         interpolation           k.a.c. dec 1994..
!..    output o3 is winter,spring,summer,fall (northern hemisphere)
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm
   use comfcst, only  : nl, nlp1
!-------------------------------------------------------------------------------
!
!     *********************************************************
!       save data on permanent data set denoted by co222 ****
!          ..... k.campana   october 1988
!     dimension t41(lp2,2),o3o3(37,l,4)
!o3   dimension sigl(l),o3o3(37,l,4)
!
   real                 ::  o3o3(37,nl,4)
   real                 ::  aa(37,nl),bb(37,nl),cc(37,nl),dd(37,nl)
!
   real                 ::  qi(82)
   real                 ::  dduo3n(19,nl),ro31(10,41),ro32(10,41),duo3n(19,41)
   real                 ::  tempn(19)
   real                 ::  o3hi(10,25),o3lo1(10,16),o3lo2(10,16),             &
                            o3lo3(10,16),o3lo4(10,16)
   real                 ::  o3hi1(10,16),o3hi2(10,9),ph1(45),                  &
                            ph2(37),p1(48),p2(33)
   real                 ::  o35deg(37,nl)
   real                 ::  rstd(81),ro3(10,41),ro3m(10,40),rbar(nl),          &
                            rdata(81),phalf(nl),pstd(nl),p(81),ph(82)
   equivalence (o3hi1(1,1),o3hi(1,1)),(o3hi2(1,1),o3hi(1,17))
   equivalence (ph1(1),ph(1)),(ph2(1),ph(46))
   equivalence (p1(1),p(1)),(p2(1),p(49))
   data ph1/      0.,                                                          &
        0.1027246e-04, 0.1239831e-04, 0.1491845e-04, 0.1788053e-04,            &
        0.2135032e-04, 0.2540162e-04, 0.3011718e-04, 0.3558949e-04,            &
        0.4192172e-04, 0.4922875e-04, 0.5763817e-04, 0.6729146e-04,            &
        0.7834518e-04, 0.9097232e-04, 0.1053635e-03, 0.1217288e-03,            &
        0.1402989e-03, 0.1613270e-03, 0.1850904e-03, 0.2119495e-03,            &
        0.2423836e-03, 0.2768980e-03, 0.3160017e-03, 0.3602623e-03,            &
        0.4103126e-03, 0.4668569e-03, 0.5306792e-03, 0.6026516e-03,            &
        0.6839018e-03, 0.7759249e-03, 0.8803303e-03, 0.9987843e-03,            &
        0.1133178e-02, 0.1285955e-02, 0.1460360e-02, 0.1660001e-02,            &
        0.1888764e-02, 0.2151165e-02, 0.2452466e-02, 0.2798806e-02,            &
        0.3197345e-02, 0.3656456e-02, 0.4185934e-02, 0.4797257e-02/
   data ph2/                                                                   &
        0.5503893e-02, 0.6321654e-02, 0.7269144e-02, 0.8368272e-02,            &
        0.9644873e-02, 0.1112946e-01, 0.1285810e-01, 0.1487354e-01,            &
        0.1722643e-01, 0.1997696e-01, 0.2319670e-01, 0.2697093e-01,            &
        0.3140135e-01, 0.3660952e-01, 0.4274090e-01, 0.4996992e-01,            &
        0.5848471e-01, 0.6847525e-01, 0.8017242e-01, 0.9386772e-01,            &
        0.1099026e000, 0.1286765e000, 0.1506574e000, 0.1763932e000,            &
        0.2065253e000, 0.2415209e000, 0.2814823e000, 0.3266369e000,            &
        0.3774861e000, 0.4345638e000, 0.4984375e000, 0.5697097e000,            &
        0.6490189e000, 0.7370409e000, 0.8344896e000, 0.9421190e000,            &
        0.1000000e001/
   data p1/                                                                    &
        0.9300000e-05, 0.1129521e-04, 0.1360915e-04, 0.1635370e-04,            &
        0.1954990e-04, 0.2331653e-04, 0.2767314e-04, 0.3277707e-04,            &
        0.3864321e-04, 0.4547839e-04, 0.5328839e-04, 0.6234301e-04,            &
        0.7263268e-04, 0.8450696e-04, 0.9793231e-04, 0.1133587e-03,            &
        0.1307170e-03, 0.1505832e-03, 0.1728373e-03, 0.1982122e-03,            &
        0.2266389e-03, 0.2592220e-03, 0.2957792e-03, 0.3376068e-03,            &
        0.3844381e-03, 0.4379281e-03, 0.4976965e-03, 0.5658476e-03,            &
        0.6418494e-03, 0.7287094e-03, 0.8261995e-03, 0.9380076e-03,            &
        0.1063498e-02, 0.1207423e-02, 0.1369594e-02, 0.1557141e-02,            &
        0.1769657e-02, 0.2015887e-02, 0.2295520e-02, 0.2620143e-02,            &
        0.2989651e-02, 0.3419469e-02, 0.3909867e-02, 0.4481491e-02,            &
        0.5135272e-02, 0.5898971e-02, 0.6774619e-02, 0.7799763e-02/
   data p2/                                                                    &
        0.8978218e-02, 0.1036103e-01, 0.1195488e-01, 0.1382957e-01,            &
        0.1599631e-01, 0.1855114e-01, 0.2151235e-01, 0.2501293e-01,            &
        0.2908220e-01, 0.3390544e-01, 0.3952926e-01, 0.4621349e-01,            &
        0.5403168e-01, 0.6330472e-01, 0.7406807e-01, 0.8677983e-01,            &
        0.1015345e000, 0.1189603e000, 0.1391863e000, 0.1630739e000,            &
        0.1908004e000, 0.2235461e000, 0.2609410e000, 0.3036404e000,            &
        0.3513750e000, 0.4055375e000, 0.4656677e000, 0.5335132e000,            &
        0.6083618e000, 0.6923932e000, 0.7845676e000, 0.8875882e000,            &
        0.1000000e001/
   data o3hi1/                                                                 &
    .55,.50,.45,.45,.40,.35,.35,.30,.30,.30,                                   &
    .55,.51,.46,.47,.42,.38,.37,.36,.35,.35,                                   &
    .55,.53,.48,.49,.44,.42,.41,.40,.38,.38,                                   &
    .60,.55,.52,.52,.50,.47,.46,.44,.42,.41,                                   &
    .65,.60,.55,.56,.53,.52,.50,.48,.45,.45,                                   &
    .75,.65,.60,.60,.55,.55,.55,.50,.48,.47,                                   &
    .80,.75,.75,.75,.70,.70,.65,.63,.60,.60,                                   &
    .90,.85,.85,.80,.80,.75,.75,.74,.72,.71,                                   &
    1.10,1.05,1.00,.90,.90,.90,.85,.83,.80,.80,                                &
    1.40,1.30,1.25,1.25,1.25,1.20,1.15,1.10,1.05,1.00,                         &
    1.7,1.7,1.6,1.6,1.6,1.6,1.6,1.6,1.5,1.5,                                   &
    2.1,2.0,1.9,1.9,1.9,1.8,1.8,1.8,1.7,1.7,                                   &
    2.4,2.3,2.2,2.2,2.2,2.1,2.1,2.1,2.0,2.0,                                   &
    2.7,2.5,2.5,2.5,2.5,2.5,2.4,2.4,2.3,2.3,                                   &
    2.9,2.8,2.7,2.7,2.7,2.7,2.7,2.7,2.6,2.6,                                   &
    3.1,3.1,3.0,3.0,3.0,3.0,3.0,3.0,2.9,2.8/
   data o3hi2/                                                                 &
    3.3,3.4,3.4,3.6,3.7,3.9,4.0,4.1,4.0,3.8,                                   &
    3.6,3.8,3.9,4.2,4.7,5.3,5.6,5.7,5.5,5.2,                                   &
    4.1,4.3,4.7,5.2,6.0,6.7,7.0,6.8,6.4,6.2,                                   &
    5.4,5.7,6.0,6.6,7.3,8.0,8.4,7.7,7.1,6.7,                                   &
    6.7,6.8,7.0,7.6,8.3,10.0,9.6,8.2,7.5,7.2,                                  &
    9.2,9.3,9.4,9.6,10.3,10.6,10.0,8.5,7.7,7.3,                                &
    12.6,12.1,12.0,12.1,11.7,11.0,10.0,8.6,7.8,7.4,                            &
    14.2,13.5,13.1,12.8,11.9,10.9,9.8,8.5,7.8,7.5,                             &
    14.3,14.0,13.4,12.7,11.6,10.6,9.3,8.4,7.6,7.3/
   data o3lo1/                                                                 &
    14.9,14.2,13.3,12.5,11.2,10.3,9.5,8.6,7.5,7.4,                             &
    14.5,14.1,13.0,11.8,10.5,9.8,9.2,7.9,7.4,7.4,                              &
    11.8,11.5,10.9,10.5,9.9,9.6,8.9,7.5,7.2,7.2,                               &
    7.3,7.7,7.8,8.4,8.4,8.5,7.9,7.4,7.1,7.1,                                   &
    4.1,4.4,5.3,6.6,6.9,7.5,7.4,7.2,7.0,6.9,                                   &
    1.8,1.9,2.5,3.3,4.5,5.8,6.3,6.3,6.4,6.1,                                   &
    0.4,0.5,0.8,1.2,2.7,3.6,4.6,4.7,5.0,5.2,                                   &
    .10,.15,.20,.50,1.4,2.1,3.0,3.2,3.5,3.9,                                   &
    .07,.10,.12,.30,1.0,1.4,1.8,1.9,2.3,2.5,                                   &
    .06,.08,.10,.15,.60,.80,1.4,1.5,1.5,1.6,                                   &
    .05,.05,.06,.09,.20,.40,.70,.80,.90,.90,                                   &
    .05,.05,.06,.08,.10,.13,.20,.25,.30,.40,                                   &
    .05,.05,.05,.06,.07,.07,.08,.09,.10,.13,                                   &
    .05,.05,.05,.05,.06,.06,.06,.06,.07,.07,                                   &
    .05,.05,.05,.05,.05,.05,.05,.06,.06,.06,                                   &
    .04,.04,.04,.04,.04,.04,.04,.05,.05,.05/
   data o3lo2/                                                                 &
    14.8,14.2,13.8,12.2,11.0,9.8,8.5,7.8,7.4,6.9,                              &
    13.2,13.0,12.5,11.3,10.4,9.0,7.8,7.5,7.0,6.6,                              &
    10.6,10.6,10.7,10.1,9.4,8.6,7.5,7.0,6.5,6.1,                               &
    7.0,7.3,7.5,7.5,7.5,7.3,6.7,6.4,6.0,5.8,                                   &
    3.8,4.0,4.7,5.0,5.2,5.9,5.8,5.6,5.5,5.5,                                   &
    1.4,1.6,2.4,3.0,3.7,4.1,4.6,4.8,5.1,5.0,                                   &
    .40,.50,.90,1.2,2.0,2.7,3.2,3.6,4.3,4.1,                                   &
    .07,.10,.20,.30,.80,1.4,2.1,2.4,2.7,3.0,                                   &
    .06,.07,.09,.15,.30,.70,1.2,1.4,1.6,2.0,                                   &
    .05,.05,.06,.12,.15,.30,.60,.70,.80,.80,                                   &
    .04,.05,.06,.08,.09,.15,.30,.40,.40,.40,                                   &
    .04,.04,.05,.055,.06,.09,.12,.13,.15,.15,                                  &
    .03,.03,.045,.052,.055,.06,.07,.07,.06,.07,                                &
    .03,.03,.04,.051,.052,.052,.06,.06,.05,.05,                                &
    .02,.02,.03,.05,.05,.05,.04,.04,.04,.04,                                   &
    .02,.02,.02,.04,.04,.04,.03,.03,.03,.03/
   data o3lo3/                                                                 &
    14.5,14.0,13.5,11.3,11.0,10.0,9.0,8.3,7.5,7.3,                             &
    13.5,13.2,12.5,11.1,10.4,9.7,8.2,7.8,7.4,6.8,                              &
    10.8,10.9,11.0,10.4,10.0,9.6,7.9,7.5,7.0,6.7,                              &
    7.3,7.5,7.8,8.5,9.0,8.5,7.7,7.4,6.9,6.5,                                   &
    4.1,4.5,5.3,6.2,7.3,7.7,7.3,7.0,6.6,6.4,                                   &
    1.8,2.0,2.2,3.8,4.3,5.6,6.2,6.2,6.4,6.2,                                   &
    .30,.50,.60,1.5,2.8,3.7,4.5,4.7,5.5,5.6,                                   &
    .09,.10,.15,.60,1.2,2.1,3.0,3.5,4.0,4.3,                                   &
    .06,.08,.10,.30,.60,1.1,1.9,2.2,2.9,3.0,                                   &
    .04,.05,.06,.15,.45,.60,1.1,1.3,1.6,1.8,                                   &
    .04,.04,.04,.08,.20,.30,.55,.60,.75,.90,                                   &
    .04,.04,.04,.05,.06,.10,.12,.15,.20,.25,                                   &
    .04,.04,.03,.04,.05,.06,.07,.07,.07,.08,                                   &
    .03,.03,.04,.05,.05,.05,.05,.05,.05,.05,                                   &
    .03,.03,.03,.04,.04,.04,.05,.05,.04,.04,                                   &
    .02,.02,.02,.04,.04,.04,.04,.04,.03,.03/
   data o3lo4/                                                                 &
    14.2,13.8,13.2,12.5,11.7,10.5,8.6,7.8,7.5,6.6,                             &
    12.5,12.4,12.2,11.7,10.8,9.8,7.8,7.2,6.5,6.1,                              &
    10.6,10.5,10.4,10.1,9.6,9.0,7.1,6.8,6.1,5.9,                               &
    7.0,7.4,7.9,7.8,7.6,7.3,6.2,6.1,5.8,5.6,                                   &
    4.2,4.6,5.1,5.6,5.9,5.9,5.9,5.8,5.6,5.3,                                   &
    2.1,2.3,2.6,2.9,3.5,4.3,4.8,4.9,5.1,5.1,                                   &
    0.7,0.8,1.0,1.5,2.0,2.8,3.5,3.6,3.7,4.0,                                   &
    .15,.20,.40,.50,.60,1.4,2.1,2.2,2.3,2.5,                                   &
    .08,.10,.15,.25,.30,.90,1.2,1.3,1.4,1.6,                                   &
    .07,.08,.10,.14,.20,.50,.70,.90,.90,.80,                                   &
    .05,.06,.08,.12,.14,.20,.35,.40,.60,.50,                                   &
    .05,.05,.08,.09,.09,.09,.11,.12,.15,.18,                                   &
    .04,.05,.06,.07,.07,.08,.08,.08,.08,.08,                                   &
    .04,.04,.05,.07,.07,.07,.07,.07,.06,.05,                                   &
    .02,.02,.04,.05,.05,.05,.05,.05,.04,.04,                                   &
    .02,.02,.03,.04,.04,.04,.04,.04,.03,.03/
!
   nkk=41
   nk=81
   nkp=nk+1
   do k = 1,nk
     ph(k)=ph(k)*1013250.
     p(k)=p(k)*1013250.
   enddo
   ph(nkp)=ph(nkp)*1013250.
   do k = 1,nl
     pstd(k)=p(k)
   enddo
!
!  write (6,3) ph
!  write (6,3) p
!  write (6,3) (pstd(k),k=1,lp1)
!***load arrays ro31,ro32,as in dicks pgm.
!
   do k = 1,25
     do n = 1,10
       ro31(n,k)=o3hi(n,k)
       ro32(n,k)=o3hi(n,k)
     enddo
   enddo
!                                                                               
   do ncase = 1,4
     itape=ncase+50
     iplace=2
     if (ncase.eq.2) iplace=4
     if (ncase.eq.3) iplace=1
     if (ncase.eq.4) iplace=3
!
!***ncase=1: spring (in n.h.)
!***ncase=2: fall   (in n.h.)
!***ncase=3: winter (in n.h.)
!***ncase=4: summer (in n.h.)
!
     if (ncase.eq.1.or.ncase.eq.2) then
       do k = 26,41
         do n = 1,10
           ro31(n,k)=o3lo1(n,k-25)
           ro32(n,k)=o3lo2(n,k-25)
         enddo
       enddo
     endif
     if (ncase.eq.3.or.ncase.eq.4) then
       do k = 26,41
         do n = 1,10
           ro31(n,k)=o3lo3(n,k-25)
           ro32(n,k)=o3lo4(n,k-25)
         enddo
       enddo
     endif                                                                     
     do kk = 1,nkk
       do n = 1,10
         duo3n(n,kk)=ro31(11-n,kk)
         duo3n(n+9,kk)=ro32(n,kk)
       enddo
       duo3n(10,kk)=.5*(ro31(1,kk)+ro32(1,kk))
     enddo
!
!***for ncase=2 or ncase=4,reverse latitude arrangement of corr. season
!
     if (ncase.eq.2.or.ncase.eq.4) then
       do kk = 1,nkk
         do n = 1,19
           tempn(n)=duo3n(20-n,kk)
         enddo
         do n = 1,19
           duo3n(n,kk)=tempn(n)
         enddo
       enddo
     endif
!
!***duo3n now is o3 profile for appropriate season,at std. pressure
!      levels
!  write (6,800) duo3n
!***begin latitude (10 deg) loop
!
     do n = 1,19
       do kk = 1,nkk
         rstd(kk)=duo3n(n,kk)
       enddo
       nkm=nk-1
       nkmm=nk-3
!
!     bessels half-point interpolation formula
!
       do k = 4,nkmm,2
         ki=k/2
         rdata(k)=.5*(rstd(ki)+rstd(ki+1))-(rstd(ki+2)-rstd(ki+1)-rstd(ki)+    &
                   rstd(ki-1))/16.
       enddo
       rdata(2)=.5*(rstd(2)+rstd(1))
       rdata(nkm)=.5*(rstd(nkk)+rstd(nkk-1))
!
!     put unchanged data into new array
!
       do k = 1,nk,2
         kq=(k+1)/2
         rdata(k)=rstd(kq)
       enddo
!
!---note to nmc: this write is commented out to reduce printout
!     write (6,798) rdata
!o3   do 23 kk=1,l
!
       do kk = 1,nl
         dduo3n(n,kk)=rdata(kk)*.01
       enddo
     enddo
!
!***end of latitude loop
!
!***create 5 deg ozone quantities by linear interpolation of
!      10 deg values
!o3   do 1060 kk=1,l
!
     do kk = 1,nl
       do n = 1,19
         o35deg(2*n-1,kk)=dduo3n(n,kk)
       enddo
       do n = 1,18
         o35deg(2*n,kk)=0.5*(dduo3n(n,kk)+dduo3n(n+1,kk))
       enddo
     enddo
!
!***output to unit (itape) the ozone values for later use
!o222  ***************************************************
!           write (66) o35deg
!
     do jj = 1,37
!o3    do 302 ken=1,l
       do ken = 1,nl
         o3o3(jj,ken,iplace) = o35deg(jj,ken)
       enddo
     enddo
!
!o222  ***************************************************
!     write (itape,101) o35deg
!     write (6,101) o35deg
   enddo
!
!***end of loop over cases
!
   aa(:,:)=o3o3(:,:,1)
   bb(:,:)=o3o3(:,:,2)
   cc(:,:)=o3o3(:,:,3)
   dd(:,:)=o3o3(:,:,4)
!
   return
   1  format(10f4.2)
   2  format(10x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x)
   3  format(10e12.5)
! 796 format(3e12.5)
  797 format(10f7.2)
! 798 format(20f6.2)
  799 format(19f6.4)
  800 format(19f6.2)
  101 format(5x,1h*,f6.5,1h,,f6.5,1h,,f6.5,1h,,f6.5,1h,,f6.5,1h,,f6.5,         &
   1h,,f6.5,1h,,f6.5,1h,,f6.5,1h,)
!
   end subroutine rad_ozone_interp_step2
!-------------------------------------------------------------------------------
#ifndef NIM
!
!-------------------------------------------------------------------------------
   subroutine rad_co2_read(nfile,rco2)
!-------------------------------------------------------------------------------
!
! subroutine:    rad_co2_read
!
! abstract:
!  -initializes arrays for -lw- radiation
!  -reads co2 transmission function data(from external file),
!   which has been pre-computed for current vertical coordinate on
!   the front-end machine. word conversion between front-end and c205
!   occurs here. also call tabl86 to set up tables for lw calculation
!   this code (rad_co2_read) is only called once...
!
! program history log:
!   84-01-01  fels and schwarzkopf,gfdl.
!   89-07-07  kenneth campana - removed unnecessary code and added
!                               reading and word conversion of co2 data.
!   89-11-29  kenneth campana - commented co2 reads because they
!                               are not yet ready for the new gfdl lw.
!
! usage:    call rad_co2_read(nfile)
!   input argument list:
!     nfile    - integer name of external co2 file.
!
!    *******************************************************************
!    *                           c o n r a d                           *
!    *    read co2 transmission data from unit(nfile)for new vertical  *
!    *      coordinate tests      ...                                  *
!    *    these arrays used to be in block data    ...k.campana-mar 90 *
!    *******************************************************************
!
!                 co2 data tables for user''s vertical coordinate
!
!   the following common blocks contain pretabulated co2 transmission
!       functions, evaluated using the methods of fels and
!       schwarzkopf (1981) and schwarzkopf and fels (1985),
!-----  the 2-dimensional arrays are
!                    co2 transmission functions and their derivatives
!        from 109-level line-by-line calculations made using the 1982
!        mcclatchy tape (12511 lines),consolidated,interpolated
!        to the nmc mrf vertical coordinatte,and re-consolidated to a
!        200 cm-1 bandwidth. the interpolation method is described in
!        schwarzkopf and fels (j.g.r.,1985).
!-----  the 1-dim arrays are
!                  co2 transmission functions and their derivatives
!          for tau(i,i+1),i=1,l,
!            where the values are not obtained by quadrature,but are the
!            actual transmissivities,etc,between a pair of pressures.
!          these used only for nearby layer calculations including qh2o.
!-----  the weighting function gtemp=p(k)**0.2*(1.+p(k)/30000.)**0.8/
!         1013250.,where p(k)=pressure,nmc mrf(new)  l18 data levels for
!         pstar=1013250.
!-----  stemp is us standard atmospheres,1976,at data pressure levels
!        using nmc mrf sigmas,where pstar=1013.25 mb (ptz program)
!====>   begin here to get constants for radiation package
!-------------------------------------------------------------------------------
#ifdef DFS
   use dfsvar, only : ib,jbw,iba,jbwa,iope
#endif
   use rdparm
   use co2dta
   use comio
#ifdef MP
   use commpi
#endif
!-------------------------------------------------------------------------------
   real                 ::  sgtmp(lp1,2),co21d(l,6),co22d(lp1,lp1,6)
   real                 ::  co21d3(lp1,6),co21d7(lp1,6)
!
   if( iope ) then
     rewind nfile
!
!       read in pre-computed co2 transmission data....
!
     do kk = 1,2
       read(nfile) (sgtmp(i,kk),i=1,lp1)
     enddo
     do kk = 1,6
       read(nfile) (co21d(i,kk),i=1,l)
     enddo
     do kk = 1,6                                                              
       read(nfile) ((co22d(i,j,kk),i=1,lp1),j=1,lp1)
     enddo
     do kk = 1,6
       read(nfile) (co21d3(i,kk),i=1,lp1)
     enddo
     do kk = 1,6                                                              
       read(nfile) (co21d7(i,kk),i=1,lp1)
     enddo
!
!  read co2 concentration in ppm (defaulted in gradfs if missing)               
!
     read(nfile,end=31) rco2
!
31   continue                                                                  
!
#ifndef NOPRINT
     write(6,*) ' co2 concentration is ',rco2
#endif
     rewind nfile
   endif
#ifdef MP
#ifdef RMP
   call rmpbcastr(sgtmp,lp1*2)
   call rmpbcastr(co21d,l*6)
   call rmpbcastr(co22d,lp1*lp1*6)
   call rmpbcastr(co21d3,lp1*6)
   call rmpbcastr(co21d7,lp1*6)
   call rmpbcastr(rco2,1)
#else
   call mpbcastr(sgtmp,lp1*2)
   call mpbcastr(co21d,l*6)
   call mpbcastr(co22d,lp1*lp1*6)
   call mpbcastr(co21d3,lp1*6)
   call mpbcastr(co21d7,lp1*6)
   call mpbcastr(rco2,1)
#endif
#endif
   do k = 1,lp1
     stemp(k) = sgtmp(k,1)
     gtemp(k) = sgtmp(k,2)
   enddo
   do k = 1,l
     cdtm51(k) = co21d(k,1)
     co2m51(k) = co21d(k,2)
     c2dm51(k) = co21d(k,3)
     cdtm58(k) = co21d(k,4)
     co2m58(k) = co21d(k,5)
     c2dm58(k) = co21d(k,6)
   enddo
   do j = 1,lp1
     do i = 1,lp1
       cdt51(i,j) = co22d(i,j,1)
       co251(i,j) = co22d(i,j,2)
       c2d51(i,j) = co22d(i,j,3)
       cdt58(i,j) = co22d(i,j,4)
       co258(i,j) = co22d(i,j,5)
       c2d58(i,j) = co22d(i,j,6)
     enddo
   enddo
   do k = 1,lp1
     cdt31(k) = co21d3(k,1)
     co231(k) = co21d3(k,2)
     c2d31(k) = co21d3(k,3)
     cdt38(k) = co21d3(k,4)
     co238(k) = co21d3(k,5)
     c2d38(k) = co21d3(k,6)
   enddo
   do k = 1,lp1
     cdt71(k) = co21d7(k,1)
     co271(k) = co21d7(k,2)
     c2d71(k) = co21d7(k,3)
     cdt78(k) = co21d7(k,4)
     co278(k) = co21d7(k,5)
     c2d78(k) = co21d7(k,6)
   enddo
#ifdef DBG
   if( iope ) then
     write(6,66)nfile
66   format(1h ,'----read co2 transmission functions from unit ',i2)
   endif
#endif
!
!......    define tables for lw radiation
!
   call rad_lw_table
!
   return
   end subroutine rad_co2_read
!-------------------------------------------------------------------------------
#endif /* no NIM */
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_table
!-------------------------------------------------------------------------------
!
!     subroutine rad_lw_table computes table entries used in the longwave radia
!     program. also calculated are indices used in strip-mining and for
!     some pre-computable functions.
!         inputs:
!         outputs:
!       em1,em1wde,table1,table2,table3         tabcom
!       em3,source,dsrce,ind,indx2,kmaxv        tabcom
!       kmaxvm,                                 tabcom
!       ao3rnd,bo3rnd,ab15                      bandta
!       ab15wd,skc1r,sko3r,sko2d                bdwide
!
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
   use rnddta
   use blockdata_rad_gfdl, only : delcm
!-------------------------------------------------------------------------------
#include <tabcom.h>
   real                 ::  sum(28,180),pertsm(28,180),sum3(28,180)
   real                 ::  sumwde(28,180),srcwd(28,nblx),src1nb(28,nblw)
   real                 ::  dbdtnb(28,nblw)
   real                 ::  zmass(181),zroot(181),sc(28),dsc(28),xtemv(28)
   real                 ::  tfour(28),fortcu(28),x(28),x1(28),x2(180),srcs(28)
   real                 ::  sum4(28),sum6(28),sum7(28),sum8(28),sum4wd(28)
   real                 ::  r1(28),r2(28),s2(28),t3(28),r1wd(28)
   real                 ::  expo(180),fac(180)
   real                 ::  cnusb(30),dnusb(30)
   real                 ::  alfanb(nblw),arotnb(nblw)
   real                 ::  anb(nblw),bnb(nblw),centnb(nblw),delnb(nblw)
   real                 ::  betanb(nblw)
!
!  common/tbltmp/ delcm(nbly)
!****************************************
!***compute local quantities and ao3,bo3,ab15
!
!....for narrow-bands...
!
   do n = 1,nblw
     anb(n)=arndm(n)
     bnb(n)=brndm(n)
     centnb(n)=haf*(bandlo(n)+bandhi(n))
     delnb(n)=bandhi(n)-bandlo(n)
     betanb(n)=betad(n)
   enddo
   ab15(1)=anb(57)*bnb(57)
   ab15(2)=anb(58)*bnb(58)
!
!....for wide bands...
!
   ab15wd=awide*bwide
!
!***compute indices: ind,indx2,kmaxv
!
   do i = 1,imax
     ind(i)=i
   enddo
!
   icnt=0
   do i1 = 1,l
     i2e=lp1-i1
     do i2 = 1,i2e
       icnt=icnt+1
       indx2(icnt)=lp1*(i2-1)+lp2*i1
     enddo
   enddo
!
   kmaxv(1)=1
   do i = 2,l
     kmaxv(i)=kmaxv(i-1)+(lp2-i)
   enddo
   kmaxvm=kmaxv(l)
!
!***compute ratios of cont. coeffs
!
   skc1r=betawd/betinw
   sko3r=betad(61)/betinw
   sko2d=one/betinw
!
!****begin table computations here***
!***compute temps, masses for table entries
!---note: the dimensioning and initialization of xtemv and other arrays
!   with dimension of 28 imply a restriction of model temperatures from
!   100k to 370k.
!---the dimensioning of zmass,zroot and other arrays with dimension of
!   180 imply a restriction of model h2o amounts such that optical paths
!   are between 10**-16 and 10**2, in cgs units.
!
   zmass(1)=h1m16
   do j = 1,180
     jp=j+1
     zroot(j)=sqrt(zmass(j))
     zmass(jp)=zmass(j)*h1p25892
   enddo
!
   do i = 1,28
     xtemv(i)=hninety+ten*i
     tfour(i)=xtemv(i)*xtemv(i)*xtemv(i)*xtemv(i)
     fortcu(i)=four*xtemv(i)*xtemv(i)*xtemv(i)
   enddo
!
!******the computation of source,dsrce is  needed only
!   for the combined wide-band case.to obtain them,the source
!   must be computed for each of the (nblx) wide bands(=srcwd)
!   then combined (using iband) into source.
!
   do n = 1,nbly
     do i = 1,28
       source(i,n)=zero
     enddo
   enddo
!
   do n = 1,nblx
     do i = 1,28
       srcwd(i,n)=zero
     enddo
   enddo
!
!---begin freq. loop (on n)
!
   do n = 1,nblx
     if (n.le.46) then
!
!***the 160-1200 band cases                                                     
!
       cent=centnb(n+16)
       del=delnb(n+16)
       bdlo=bandlo(n+16)
       bdhi=bandhi(n+16)
     endif
     if (n.eq.nblx) then
!
!***the 2270-2380 band case                                                     
!
       cent=centnb(nblw)
       del=delnb(nblw)
       bdlo=bandlo(nblw)
       bdhi=bandhi(nblw)
     endif
!
!***for purposes of accuracy, all evaluations of planck fctns are made
!  on 10 cm-1 intervals, then summed into the (nblx) wide bands.
!
     nsubds=(del-h1m3)/10+1
     do nsb = 1,nsubds
       if (nsb.ne.nsubds) then
         cnusb(nsb)=ten*(nsb-1)+bdlo+five
         dnusb(nsb)=ten
       else
         cnusb(nsb)=haf*(ten*(nsb-1)+bdlo+bdhi)
         dnusb(nsb)=bdhi-(ten*(nsb-1)+bdlo)
       endif
       c1=(h37412m5)*cnusb(nsb)**3
!
!---begin temp. loop (on i)
!
       do i = 1,28
         x(i)=h1p4387*cnusb(nsb)/xtemv(i)
         x1(i)=exp(x(i))
         srcs(i)=c1/(x1(i)-one)
         srcwd(i,n)=srcwd(i,n)+srcs(i)*dnusb(nsb)
       enddo
     enddo
   enddo
!
!***the following loops create the combined wide band quantities source
!   and dsrce
!
   do n = 1,40
     do i = 1,28
       source(i,iband(n))=source(i,iband(n))+srcwd(i,n)
     enddo
   enddo
!
   do n = 9,nbly
     do i = 1,28
       source(i,n)=srcwd(i,n+32)
     enddo
   enddo
!
   do n = 1,nbly
     do i = 1,27
       dsrce(i,n)=(source(i+1,n)-source(i,n))*hp1
     enddo
   enddo
!
   do n = 1,nblw
     alfanb(n)=bnb(n)*anb(n)
     arotnb(n)=sqrt(alfanb(n))
   enddo
!
!***first compute planck fctns (src1nb) and derivatives (dbdtnb) for
!   use in table evaluations. these are different from source,dsrce
!   because different frequency pts are used in evaluation, the freq.
!   ranges are different, and the derivative algorithm is different.
!
   do n = 1,nblw
     cent=centnb(n)
     del=delnb(n)
!
!---note: at present, the ia loop is only used for ia=2. the loop struct
!   is kept so that in the future, we may use a quadrature scheme for
!   the planck fctn evaluation, rather than use the mid-band frequency.
!
     do ia = 1,3
       anu=cent+haf*(ia-2)*del
       c1=(h37412m5)*anu*anu*anu+h1m20
!
!---temperature loop---
!
#ifdef SX6
!CDIR NOVECTOR
#endif
        do i = 1,28
          x(i)=h1p4387*anu/xtemv(i)
          x1(i)=exp(x(i))
          sc(i)=c1/((x1(i)-one)+h1m20)
          dsc(i)=sc(i)*sc(i)*x(i)*x1(i)/(xtemv(i)*c1)
        enddo
        if (ia.eq.2) then
          do i = 1,28
            src1nb(i,n)=del*sc(i)
            dbdtnb(i,n)=del*dsc(i)
          enddo
        endif
      enddo
   enddo
!
!***next compute r1,r2,s2,and t3- coefficients used for e3 function
!   when the optical path is less than 10-4. in this case, we assume a
!   different dependence on (zmass).
!---also obtain r1wd, which is r1 summed over the 160-560 cm-1 range
!
   do i = 1,28
     sum4(i)=zero
     sum6(i)=zero
     sum7(i)=zero
     sum8(i)=zero
     sum4wd(i)=zero
   enddo
   do n = 1,nblw
     cent=centnb(n)
!
!***perform summations for freq. ranges of 0-560,1200-2200 cm-1 for sum4
!   sum6,sum7,sum8
!
     if (cent.lt.560. .or. cent.gt.1200..and.cent.le.2200.) then
       do i = 1,28
         sum4(i)=sum4(i)+src1nb(i,n)
         sum6(i)=sum6(i)+dbdtnb(i,n)
         sum7(i)=sum7(i)+dbdtnb(i,n)*arotnb(n)
         sum8(i)=sum8(i)+dbdtnb(i,n)*alfanb(n)
       enddo
     endif
!
!***perform summations over 160-560 cm-1 freq range for e1 calcs (sum4wd
!
     if (cent.gt.160. .and. cent.lt.560.) then
       do i = 1,28
         sum4wd(i)=sum4wd(i)+src1nb(i,n)
       enddo
     endif
   enddo
!
   do i = 1,28
     r1(i)=sum4(i)/tfour(i)
     r2(i)=sum6(i)/fortcu(i)
     s2(i)=sum7(i)/fortcu(i)
     t3(i)=sum8(i)/fortcu(i)
     r1wd(i)=sum4wd(i)/tfour(i)
   enddo
!
   do j = 1,180
     do i = 1,28
       sum(i,j)=zero
       pertsm(i,j)=zero
       sum3(i,j)=zero
       sumwde(i,j)=zero
     enddo
   enddo
!
!---frequency loop begins---                                                    
!
   do n = 1,nblw
     cent=centnb(n)
!
!***perform calculations for freq. ranges of 0-560,1200-2200 cm-1               
!
     if (cent.lt.560. .or. cent.gt.1200..and.cent.le.2200.) then
       do j = 1,180
         x2(j)=arotnb(n)*zroot(j)
         expo(j)=exp(-x2(j))
       enddo
       do j = 1,180
         if (x2(j).ge.hundred) then
           expo(j)=zero
         endif
       enddo
       do j = 121,180
         fac(j)=zmass(j)*(one-(one+x2(j))*expo(j))/(x2(j)*x2(j))
       enddo
       do j = 1,180
         do i = 1,28
           sum(i,j)=sum(i,j)+src1nb(i,n)*expo(j)
           pertsm(i,j)=pertsm(i,j)+dbdtnb(i,n)*expo(j)
         enddo
       enddo
       do j = 121,180
         do i = 1,28
           sum3(i,j)=sum3(i,j)+dbdtnb(i,n)*fac(j)
         enddo
       enddo
     endif
!
!---compute sum over 160-560 cm-1 range for use in e1 calcs (sumwde)
!
     if (cent.gt.160. .and. cent.lt.560.) then
       do j = 1,180
          do i = 1,28
            sumwde(i,j)=sumwde(i,j)+src1nb(i,n)*expo(j)
          enddo
       enddo
     endif
   enddo
!
   do j = 1,180
     do i = 1,28
       em1(i,j)=sum(i,j)/tfour(i)
       table1(i,j)=pertsm(i,j)/fortcu(i)
     enddo
   enddo
!
   do j = 121,180
     do i = 1,28
       em3(i,j)=sum3(i,j)/fortcu(i)
     enddo
   enddo
!
   do j = 1,179
     do i = 1,28
       table2(i,j)=(table1(i,j+1)-table1(i,j))*ten
     enddo
   enddo
!
   do j = 1,180
     do i = 1,27
       table3(i,j)=(table1(i+1,j)-table1(i,j))*hp1
     enddo
   enddo
!
   do i = 1,28
     table2(i,180)=zero
   enddo
!
   do j = 1,180
     table3(28,j)=zero
   enddo
!
   do j = 1,2
     do i = 1,28
       em1(i,j)=r1(i)
     enddo
   enddo
!
   do j = 1,120
     do i = 1,28
       em3(i,j)=r2(i)/two-s2(i)*sqrt(zmass(j))/three+t3(i)*zmass(j)/eight
     enddo
   enddo
!
   do j = 121,180
     do i = 1,28
       em3(i,j)=em3(i,j)/zmass(j)
     enddo
   enddo
!
!***now compute e1 tables for 160-560 cm-1 bands only.
!   we use r1wd and sumwde obtained above.
!
   do j = 1,180
     do i = 1,28
       em1wde(i,j)=sumwde(i,j)/tfour(i)
     enddo
   enddo
!
   do j = 1,2
     do i = 1,28                                                             
       em1wde(i,j)=r1wd(i)
     enddo
   enddo
!
   return
   end subroutine rad_lw_table
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_aeros_init(si,sl)  
!-------------------------------------------------------------------------------
!
!  subprogram:  rad_aeros_init
!
!  abstract:    setup common block 'swaer' for aerosols and rayleigh scattering
!               optical properties in four uv+vis bands and four nir bands.
!             * band: 1. 0.225-0.285 (uv)       2. 0.175-0.225;0.285-0.300 (uv)
!                     3. 0.300-0.325 (uv0       4. 0.325-0.690 (par)
!                     5. 2.27 - 4.0  (nir)      6. 1.22 - 2.27 (nir)
!                     7. 0.70 - 1.22 (nir)      8. 0.70 - 4.0  (nir)
!
!  reference: wmo report wcp-112 (1986)
!                                                                               
!  the six typical aerosol profiles:
!     1.ubr; 2.cont-1; 3.mar-1; 4.cont-2; 5.mar-2; 6.conv
!
!     sigref  - ref. sigma level                   n/d   ndm*nae
!  arrays in the common block:
!     sig0    - ratio of sigl to domain boundary   n/d   l*nae
!     haer    - scale height of aerosols           km    ndm*nae
!     zaer    - ext. coef. of aerosols             1/km  ndm*nae
!     hh      - atmospheric scale height           km    l
!     dz      - layer thickness                    km    l
!     idm     - domain index                       n/d   l*nae
!     oaer    - single scattering albedo           n/d   ndm*nae*nbd
!     gaer    - asymmetry parameter                n/d   ndm*nae*nbd
!     raer    - ratio of band wavelength to the
!               refference wavelength (0.55 micron)n/d   nbd
!     taur    - rayleigh scattering optical depth  n/d   l*nbd
!
!-------------------------------------------------------------------------------
   use rdparm, only : l, lp1, nbd, nae
   use comswaer
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                 ::  si(lp1),sl(l),sigln(lp1)
   integer              ::  n,k,iaer,idom
   real                 ::  sbund
!
!===> ... compute layer distributions of rayleigh scattering
!                                                                               
   do n = 1,nbd
     do k = 1,l
       taur(k,n) = tauray(n) * (si(lp1-k)-si(lp1-k+1))
      enddo
   enddo
!
!===> ... setup log sigma array (set toa sigma=0.0001)
!         rem: si,sl k=1 is sfc; but in radiation k=1 is toa
!
   sigln(1) = alog(1.0e-4)
   do k = 1,l
     sigln(k+1) = alog(si(lp1-k))
   enddo
!
   do k = 1,l
     hh(k) = 6.05e0 + 2.5e0 * sl(lp1-k)
     dz(k) = hh(k) * (sigln(k+1)-sigln(k))
   enddo
!                                                                               
   do iaer = 1,nae
     sbund= si(1)
     idom = 1
     do k = l,1,-1
       if (si(lp1-k+1) .lt. sigref(idom,iaer)) then
         idom = idom + 1
         sbund= si(lp1-k)
       endif
       sig0(k,iaer) = sl(lp1-k) / sbund
       idm (k,iaer) = idom
     enddo
   enddo
!
   return
   end subroutine rad_aeros_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_aeros_init_ice(si,sl)
!-------------------------------------------------------------------------------
!
!  subprogmam:  rad_aeros_init_ice
!
!  abstract:    setup common block 'swaer' for aerosols and rayleigh scattering
!               optical properties in four uv+vis bands and four nir bands.
!           *  band: 1. 0.175-0.225 (uv-c)     2. 0.225-0.245;0.260-0.280 (uv-c)
!                    3. 0.245-0.260 (uv-c)     4. 0.280-0.295 (uv-b)
!                    5. 0.295-0.310 (uv-b)     6. 0.310-0.320 (uv-b)
!                    7. 0.320-0.400 (uv-a)     8. 0.400-0.700 (par)
!                    9. 2.27 - 4.0  (nir)     10. 1.22 - 2.27 (nir)
!                   11. 0.70 - 1.22 (nir)     12. 0.70 - 4.0  (nir)
!
!  reference: wmo report wcp-112 (1986)
!
!  the six typical aerosol profiles:
!     1.ubr; 2.cont-1; 3.mar-1; 4.cont-2; 5.mar-2; 6.conv
!
!     sigref  - ref. sigma level                   n/d   ndm*nae
!  arrays in the common block:
!     sig0    - ratio of sigl to domain boundary   n/d   l*nae
!     haer    - scale height of aerosols           km    ndm*nae
!     zaer    - ext. coef. of aerosols             1/km  ndm*nae
!     hh      - atmospheric scale height           km    l
!     dz      - layer thickness                    km    l
!     idm     - domain index                       n/d   l*nae
!     oaer    - single scattering albedo           n/d   ndm*nae*nbd
!     gaer    - asymmetry parameter                n/d   ndm*nae*nbd
!     raer    - ratio of band wavelength to the
!               refference wavelength (0.55 micron)n/d   nbd
!     taur    - rayleigh scattering optical depth  n/d   l*nbd
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm8
   use comswaer8
!-------------------------------------------------------------------------------
   real                 ::  si(lp1),sl(l),sigln(lp1)
!-------------------------------------------------------------------------------
!
!===> ... compute layer distributions of rayleigh scattering
!
   do n = 1,nbd
     do k = 1,l
       taur(k,n) = tauray(n) * (si(lp1-k)-si(lp1-k+1))
     enddo
   enddo
!
!===> ... setup log sigma array (set toa sigma=0.0001)
!         rem: si,sl k=1 is sfc; but in radiation k=1 is toa
!
   sigln(1) = alog(1.0e-4)
   do k = 1,l
     sigln(k+1) = alog(si(lp1-k))
   enddo
   do k = 1,l
     hh(k) = 6.05e0 + 2.5e0 * sl(lp1-k)
     dz(k) = hh(k) * (sigln(k+1)-sigln(k))
   enddo
!
   do iaer = 1,nae
     sbund= si(1)
     idom = 1
     do k = l,1,-1
       if (si(lp1-k+1) .lt. sigref(idom,iaer)) then
         idom = idom + 1
         sbund= si(lp1-k)
       endif
       sig0(k,iaer) = sl(lp1-k) / sbund
       idm (k,iaer) = idom
     enddo
   enddo
!
   return
   end subroutine rad_aeros_init_ice
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_aeros_init_nasa(si,sl)
!-------------------------------------------------------------------------------
!
!  subprogram:  rad_aeros_init_nasa
! 
!  abstract:    setup common block 'swaer' for aerosols and rayleigh scattering
!               optical properties in four uv+vis bands and four nir bands.
!           *  band: 1. 0.175-0.225 (uv-c)     2. 0.225-0.245;0.260-0.280 (uv-c)
!                    3. 0.245-0.260 (uv-c)     4. 0.280-0.295 (uv-b)
!                    5. 0.295-0.310 (uv-b)     6. 0.310-0.320 (uv-b)
!                    7. 0.320-0.400 (uv-a)     8. 0.400-0.700 (par)
!                    9. 2.27 - 4.0  (nir)     10. 1.22 - 2.27 (nir)
!                   11. 0.70 - 1.22 (nir)     12. 0.70 - 4.0  (nir)
!
!  reference: wmo report wcp-112 (1986)
!                                                                               
!  setup common block 'swaer' for aerosols and rayleigh scattering
!  optical properties in eight uv+vis bands and four nir bands.
!  ref: wmo report wcp-112 (1986)
!
!  there are seven typical vertical structures:
!     1.antarctic, 2.arctic, 3.continent, 4.maritime, 5.desert,
!     6.maritime with mineral overlay, 7.continent with mineral overlay
!
!     sigref  - ref. sigma level                   n/d   ndm*nae
!  arrays in the common block:
!     haer    - scale height of aerosols           km    ndm*nae
!     hh      - atmospheric scale height           km    l
!     hz      - level height                       km    l+1
!     dz      - layer thickness                    km    l
!     idm     - domain index                       n/d   l*nae
!     taur    - rayleigh scattering optical depth  n/d   l*nbd
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm99
   use comswaer99
!-------------------------------------------------------------------------------
   real                 ::  si(lp1),sl(l),sigln(lp1)
!
!===> ... compute layer distributions of rayleigh scattering
!
   do n = 1,nbd
     do k = 1,l
       taur(k,n) = tauray(n) * (si(lp1-k)-si(lp1-k+1))
     enddo
   enddo
!
!===> ... setup log sigma array (set toa sigma=0.0001)
!         rem: si,sl k=1 is sfc; but in radiation k=1 is toa
!
   sigln(1) = alog(1.0e-4)
   do k = 1,l
     sigln(k+1) = alog(si(lp1-k))
   enddo
!
   do k = 1,l
     hh(k) = 6.05e0 + 2.5e0 * sl(lp1-k)
     dz(k) = hh(k) * (sigln(k+1)-sigln(k))
   enddo
!
   hz(lp1) = 0.0
   do k = l,1,-1
     hz(k) = hz(k+1) + dz(k)
   enddo
!
   do iaer = 1,nae
     idom = 1
     do k = l,1,-1
       if (si(lp1-k+1) .lt. sigref(idom,iaer)) then
         idom = idom + 1
         if (idom.eq.2 .and.sigref(2,iaer).eq.sigref(3,iaer)) then
           idom = 3
         endif
       endif
       idm (k,iaer) = idom
     enddo
   enddo
!
   return
   end subroutine rad_aeros_init_nasa
!------------------------------------------------------------------------------
#ifndef RMP
!
!------------------------------------------------------------------------------
   subroutine rad_albedo_aeros_read(nflin,alvsf,alnsf,alvwf,alnwf,             &
                      facsf,facwf                                              &
#ifndef SWRMDC
                     ,paerf)
#else
                     )
#endif
!-------------------------------------------------------------------------------
!
!  the code reads in surface albedo and aerosol data. the albedo
!  data are derived from matthews vegetation index by using a
!  modified brieglebs scheme.  the aerosol distribution data is
!  based on matthews vegetation index.    --- y.hou  mar 7, 1995
!
!-------------------------------------------------------------------------------
#ifdef MP
   use paramodel, only    :  lonf2p_,latg2p_
#endif
#ifdef NIM
   use paramodel, only    :  LONF2S,LATG2S,nip_, ix=>nip_,iy=>latg2_
#else
   use paramodel, only    :  LONF2S,LATG2S, ix=>lonf2_,iy=>latg2_
   use module_trans, only :  dyn_trans2model_grid
#endif
#if defined(DFS)
   use dfsvar, only       :  iope
#else
#ifdef MP
   use commpi
#endif
#endif
   use comio
#ifdef MP
   real,allocatable      ::  grid(:,:,:)
#endif
#ifdef NIM
   real                  ::  alvsf(nip_,1,4),alnsf(nip_,1,4)
   real                  ::  alvwf(nip_,1,4),alnwf(nip_,1,4)
   real                  ::  facsf(nip_,1)  ,facwf(nip_,1)
#ifndef SWRMDC
   real                  ::  paerf(nip_,1,5)
#endif
#else
   real                  ::  alvsf(LONF2S,LATG2S,4),alnsf(LONF2S,LATG2S,4)
   real                  ::  alvwf(LONF2S,LATG2S,4),alnwf(LONF2S,LATG2S,4)
   real                  ::  facsf(LONF2S,LATG2S)  ,facwf(LONF2S,LATG2S)
#ifndef SWRMDC
   real                  ::  paerf(LONF2S,LATG2S,5)
#endif
#endif /* NIM */
!
#ifdef MP
   allocate(grid(ix,iy,5))
   if(iope) then
     rewind nflin
     do k = 1,4
       read(nflin) grid(:,:,k)
     enddo
     call dyn_trans2model_grid  (grid,4)
   endif
!
   call mpgf2p(grid,ix,iy,alvsf,lonf2p_,latg2p_,4)
   if (iope) then
     do k = 1,4
       read(nflin) grid(:,:,k)
     enddo
     call dyn_trans2model_grid  (grid,4)
   endif
!
   call mpgf2p(grid,ix,iy,alvwf,lonf2p_,latg2p_,4)
   if (iope) then
     do k = 1,4
       read(nflin) grid(:,:,k)
     enddo
     call dyn_trans2model_grid(grid,4)
   endif
!
   call mpgf2p(grid,ix,iy,alnsf,lonf2p_,latg2p_,4)
   if (iope) then
     do k = 1,4
       read(nflin) grid(:,:,k)
     enddo
     call dyn_trans2model_grid(grid,4)
   endif
!
   call mpgf2p(grid,ix,iy,alnwf,lonf2p_,latg2p_,4)
   if (iope) then
     read(nflin) grid(:,:,1)
     call dyn_trans2model_grid(grid,1)
   endif
!
   call mpgf2p(grid,ix,iy,facsf,lonf2p_,latg2p_,1)
   if (iope) then
     read(nflin) grid(:,:,1)
     call dyn_trans2model_grid(grid,1)
   endif
!
   call mpgf2p(grid,ix,iy,facwf,lonf2p_,latg2p_,1)
#ifndef SWRMDC
   if (iope) then
     do k = 1,5
       read(nflin) grid(:,:,k)
     enddo
     call dyn_trans2model_grid(grid,5)
   endif
!
   call mpgf2p(grid,ix,iy,paerf,lonf2p_,latg2p_,5)
#endif
   deallocate(grid)
!
#else
   rewind nflin
#ifdef NIM
#ifdef NIM_NOAER
   alvsf=0.001 ; alvwf=0.008 ; alnsf=0.002 ; alnwf=0.003
   facsf=0.008 ; facwf=0.008 ; paerf=0.02
#else
   do k = 1,4
     read(nflin) ((alvsf(i,j,k),i=1,ix),j=1,iy)
   enddo
   do k = 1,4
     read(nflin) ((alvwf(i,j,k),i=1,ix),j=1,iy)
   enddo
   do k = 1,4
     read(nflin) ((alnsf(i,j,k),i=1,ix),j=1,iy)
   enddo
   do k = 1,4
     read(nflin) ((alnwf(i,j,k),i=1,ix),j=1,iy)
   enddo
!
   read(nflin) facsf
   read(nflin) facwf
#endif
#else /* ~NIM */
   do k = 1,4
     read(nflin) ((alvsf(i,j,k),i=1,ix),j=1,iy)
   enddo
   call dyn_trans2model_grid  (alvsf,4)
   do k = 1,4
     read(nflin) ((alvwf(i,j,k),i=1,ix),j=1,iy)
   enddo
   call dyn_trans2model_grid  (alvwf,4)
   do k = 1,4
     read(nflin) ((alnsf(i,j,k),i=1,ix),j=1,iy)
   enddo
   call dyn_trans2model_grid  (alnsf,4)
   do k = 1,4
     read(nflin) ((alnwf(i,j,k),i=1,ix),j=1,iy)
   enddo
   call dyn_trans2model_grid  (alnwf,4)
!
   read(nflin) facsf
   read(nflin) facwf
   call dyn_trans2model_grid(facsf,1)
   call dyn_trans2model_grid(facwf,1)
#endif /* NIM end */
#ifndef SWRMDC
#ifndef NIM_NOAER
   do k = 1,5
     read(nflin) ((paerf(i,j,k),i=1,ix),j=1,iy)
   enddo
#endif
#ifndef NIM
   call dyn_trans2model_grid(paerf,5)
#endif
#endif /* no SWRMDC */
#endif /* MP */
!
#ifdef NIMAQUA
! homogeneous aersol distribution for aqua-planet
!
   paerf(:,:,2)=0.02 ; paerf(:,:,3)=0.2 ; paerf(:,:,4)=0.02
!
!-- If sea
!  alvsf(:,:,:)=0.001 ; alvwf(:,:,:)=0.008 ; alnsf(:,:,:)=0.002 ;
!  alnwf(:,:,:)=0.003 ; facsf(:,:)=0.008 ; facwf(:,:)=0.008 ; paerf(:,:,2)=0.02 ;
!  ; paerf(:,:,3)=0.2 ;  ; paerf(:,:,4)=0.02 ;
!-- If land
!  alvsf(:,:,:)=0.01 ; alvwf(:,:,:)=0.1 ; alnsf(:,:,:)=0.02 ; alnwf(:,:,:)=0.03 ;
!  facsf(:,:)=0.08 ; facwf(:,:)=0.08
!
#endif
   return
   end subroutine rad_albedo_aeros_read
!------------------------------------------------------------------------------
#endif  /* ~RMP end */
!
#ifndef RMP
!------------------------------------------------------------------------------
#ifndef NIM
   subroutine rad_prepare(idate,dorad)
#else
   subroutine rad_prepare(idate,dorad,lat,lon,igs,ige)
#endif
!------------------------------------------------------------------------------
   use  comfgrid        ! latdef 
   use  comfphys        ! coszen,coszer
   use  comgrad         ! runrad,kalb,albedr,cvr,cvtr,cvbr,slmskr,nalaer
   use  comio           ! iope
   use  comfver
   use  comfgsm
   use  radiag
#ifdef MP
#ifndef DFS
   use  commpi
#endif
   use  paramodel, only : lonf2p_,latg2p_,lonf2_,latg2_
#endif
#else   /* if RMP */
!------------------------------------------------------------------------------
   subroutine rad_prepare(dorad)
!------------------------------------------------------------------------------
#ifdef MP
   use  commpi
#endif
   use  paramodel, only : jgrd12_
   use  rscomf_rerun
   use  rscomltb
   use  rscommap
   use  rscomgrad
   use  rdparm
#endif
   use  paramodel, only : LATG2S,LONF2S,LONF2F,LATG2F,jcap_,levs_,ngases_
   use  varsfc, only    : lalbd_
   use  rnddta, only    : init_rnddta
   use  blockdata_rad_gfdl     ! subroutine : init_rad_gfdl_blockdata
                               !              init_coef_rad_gfdl_blockdata
                               ! season,lseason,lftype,jdnmc,ljdnmc
                               ! jtyme,ixxxx,lixxxx
                               ! fcstda,daz,fjdnmc,tslag,rlag,timin
                               ! tpi,hpi,year,day,dhr,sc
   use  comsfc
   use  comfcst, only   : loz, psnasa, o3nasa
   use  comreado3, only : o3out, pstr, jerr
#ifdef DFS
   use  dfsvar, only    : jgs,jls,iope,latdef
#endif
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
#include "abort.h"
#ifdef MP
   real                 ::  cvf  (LONF2F,LATG2F)
   real                 ::  cvtf (LONF2F,LATG2F)
   real                 ::  cvbf (LONF2F,LATG2F)
   real                 ::  cvrf (LONF2F,LATG2F)
   real                 ::  cvtrf(LONF2F,LATG2F)
   real                 ::  cvbrf(LONF2F,LATG2F)
#endif
!
#ifndef RMP
   integer              ::  idate(4)
#endif
#ifdef NIM
   integer              ::  igs,ige 
   real                 ::  lat(LONF2F),lon(LONF2F),cvmod
#endif
   logical              ::  dorad
   real                 ::  dtsmod,dtlmod,chour,xmin,temp
   integer              ::  i,j,n,iyr,imon,iday,iztim,kyear,munth
   integer              ::  im,id,iyear,ihr,jd
   integer              ::  its,ite,jts,jte
   integer              ::  latstart
   integer              ::  ifunit
#ifdef RMP
   integer              ::  k
#ifndef MP
   integer              ::  latdef(jgrd12_)
#endif
!
   do k = 1,jgrd12_
     latdef(k)=k
   enddo
#endif
#ifdef NIM
!
! sinlab, coslab : to cal zenith angle
!
   sinlab(:,1)=sin(lat(:))
   coslab(:,1)=cos(lat(:))
   xlat(:,1)=lat(:)
   xlon(:,1)=lon(:)
   latdef(:)=1   ! not used for NIM
!
! define the dimension constants for NIM
!
   its = 1
   ite = ige
   jts = 1
   jte = 1
#else
   its = 1
   ite = LONF2S
   jts = 1
   jte = LATG2S
#endif
#ifdef NIM
!
!   set switch for saving kuo data (for interactive clouds)..
!     ... from gmp_integrate.F90
!
   cvmod= mod (solhr+dthr,dtcvav)
   if(cvmod.lt.hdthr.or.cvmod.ge.dtcvav-hdthr) then
     clstp=min(dtcvav,(shour+deltim)/3600.)
   elseif(clstp.gt.0.) then
     clstp=0.
   else
     clstp=-10.
   endif
   call dyn_flux_zero_out(0)    ! write diagnostic files
#endif
!
   dthr = deltim / 3600.0
   hdthr = 0.5 * dthr
!
   itimsw = 0
   itimlw = 0
   dtsmod = amod(solhr,dtswav)
   if(inistp.ne.0 .or. dtsmod.lt.hdthr .or. dtsmod.ge.dtswav-hdthr)            &
     itimsw = 1
   dtlmod = amod(solhr,dtlwav)
   if(inistp.ne.0 .or. dtlmod.lt.hdthr .or. dtlmod.ge.dtlwav-hdthr)            &
     itimlw = 1
   dorad=.false.
   if(itimsw.eq.0 .and. itimlw.eq.0) return
   if(.not.runrad .and. inistp.ne.0) return
   dorad=.true.
!
!dg3     compute num sec between calls to cld code...for diagnostics
!
#ifdef DG
   dtacc = min(dtswav,dtlwav)*3600.
#endif
#ifdef DBG
#ifdef MP
   if(iope) then
#endif
     write(6,1001) jcap_, levs_
1001 format (1h0,'gfdl reduced rad ',i2,i2,'g,e typ, feb 20 1986')
#ifdef MP
   endif
#endif
#endif
!
!  paerf when SWRMDC =.false.
!  alvsf,alnsf,alvwf,alnwf from the rad_albedo_aerosol.snl will be discarded
!
   if(inistp.ne.0) then
     call init_rnddta()
     call init_coef_rad_gfdl_blockdata
     call init_rad_gfdl_blockdata
   endif
!
! setup coeff b0-b3 for restart
!
   if(fhour.ne.0.and.kdt.eq.1.and.inistp.eq.0) then
     if(iope) print*, 'call init_coef_rad for fhour=!0 and inistp=0'
     call init_coef_rad_gfdl_blockdata
   endif
#ifdef RMP
!
   if(kalb.eq.0) then
     do j = 1,LATG2S
       do i = 1,LONF2S
         albed (i,j) = albedo(i,j,1)
       enddo
     enddo
   endif
#endif
#ifndef RMP
!
   do j = jts,jte
     do i = its,ite
       if(kalb.eq.0) then
         albedr(i,j) = albedo(i,j,1)
       endif
       slmskr(i,j) =  slmsk(i,j)
     enddo
   enddo
#endif
!
!   astronomy calculations-once for each new radiation step
! 
!      get 4 digit year for julian day computation
!
   iyr = idate(4)
   imon = idate(2)
   iday = idate(3)
   iztim = idate(1)
#ifdef SOLAR_FIXED
   iyr = 2012
   imon = 3
   iday = 21
   iztim = 0
#endif
   if(iyr.lt.100) then
     write(6,*) 'iy .lt.100 in rad_prepare'
     call MPABORT
   endif
   kyear = iyr
#ifdef NIM
   chour=shour/3600.
#else
#ifndef RMP
   chour=fhour+shour/3600.
#else
   chour=fhour
   if(inistp.eq.0) chour=thour
#endif
#endif
#ifdef SOLAR_FIXED
   chour=amod(chour,24.)
#endif
!
! following broadcasts are to solve the bit reproducibility issue
!
   call rad_julian_day(kyear,imon,iday,iztim,0,jdnmc,fjdnmc)
   call rad_fcst_time(chour,imon,iday,iztim,jdnmc,fjdnmc,rlag,year,            &
                      rsin1,rcos1,rcos2,jd,fjd)
   if(itimsw.eq.1) then
     call rad_solar_time(jd,fjd,r1,dlt,alf,slag,sdec,cdec)
     solc=sc/(r1*r1)
     call rad_cos_zenith(LONF2S,LATG2S,                                        &
#ifndef RMP
                         dtswav,solhr,sinlab,coslab,sdec,cdec,slag,            &
                         xlon,coszer,.true.,coszdg)
#else
                         dtswav,solhr,sinlar,coslar,sdec,cdec,slag,            &
                         rlon,coszer,.true.,coszdg)
#endif
   endif
!
   call rad_calc_time(jd,fjd,munth,im,id,iyear,ihr,xmin)
   call rad_print_time(id,munth,iyear,ihr,xmin,jd,fjd,dlt,alf,r1,slag,solc)
!
!   nasa o3   calculations-once for each new radiation step      
!   get new climo from nasa 12-month sbuv data                  
!   if jo3=0 use old gfdl climo, jo3=1 use new nasa climo ozone
!
   if(thour.eq.0.or.ngases_.lt.1) then
     jo3=0
     o3nasa(:,:)=o3out(:,:)    ! from rad_ozone_time_interp in read_grims
     if (jerr.le.0) then
!
!  if nasa data file was available (jerr=0), not (jerr=1)
!
#ifdef AQUA_PLANET
       do n = 1,loz
         temp = 0.0
         do i = 1,37
           temp = temp + o3nasa(i,n)
         enddo
         do i = 1,37
           o3nasa(i,n) = temp/37.
         enddo
       enddo
#endif
       jo3=1
!
!   get nasa pressure in cb (flip vertical coordinate)
!
       do n = 1,loz
         psnasa(n) = pstr(loz+1-n)*1.e-1
       enddo
     endif
#ifdef DBG
     if(iope ) then
       write(6,167)
       if (jo3.eq.0) write(6,166)
       if (jo3.eq.1) write(6,168)
       write(6,167)
166 format('    using gfdl zonal seasonal ozone climo ')
167 format('       ------        ')
168 format('    using nasa zonal monthly ozone climo ')
     endif
#endif
   endif
#ifdef DBG
   if (iope) then
     call print_maxmin_six(tsea  ,LONF2S*LATG2S,1,1,1,                         &
                                                            'tsea rad_prepare')
     call print_maxmin_six(snoweq,LONF2S*LATG2S,1,1,1,                         &
                                                          'snoweq rad_prepare')
     call print_maxmin_six(albedo,LONF2S*LATG2S,lalbd_,1,lalbd_,               &
                                                          'albedo rad_prepare')
     call print_maxmin_six(slmsk ,LONF2S*LATG2S,1,1,1,                         &
                                                           'slmsk rad_prepare')
   endif
#endif
!
!  albedo and aerosol data interpolations
!
   if(itimsw.eq.1) then
#ifndef RMP
     do j = jts,jte
       do i = its,ite
         coszen(i,j)=coszer(i,j)
       enddo
     enddo
#endif
#ifdef DFS
     latstart=jls
#else
#ifdef MP
     latstart=latstr(mype)
#else
     latstart=1
#endif
#endif
!
     call rad_albedo_aerosol(LONF2S,LATG2S,                                    &
#ifndef RMP
                 im,slmsk,snoweq,z0cm,coszen,tsea,hprime,jsno,                 &
#else
                 im,slmsk,snoweq,z0cm,coszer,tsea,hprime,jsno,                 &
#endif
                 albedo(1,1,1),albedo(1,1,3),albedo(1,1,2),                    &
                 albedo(1,1,4),facalf(1,1,1),facalf(1,1,2),                    &
#ifndef SWRMDC
                 paer,                                                         &
#endif
#ifndef RMP
                 xlat,alvbr,alnbr,alvdr,alndr,latdef(latstart)                 &
#else
                 rlat,alvbr,alnbr,alvdr,alndr,latdef(latstart)                 &
#endif
#ifdef SWRMDC
                ,kprfi,idxci,cmixi,denni                                       &
                ,iswsrc(1),nfaer,kprfg,idxcg,cmixg,denng)
#else
                ,paerr)
#endif
   endif
#ifdef DBG
   if (iope) then
     call print_maxmin_six(tsea,LATG2S*LONF2S,1,1,1,'tsea in rad_prepare')
     call print_maxmin_six(snoweq,LATG2S*LONF2S,1,1,1,'snoweq in rad_prepare')
     call print_maxmin_six(slmsk,LATG2S*LONF2S,1,1,1,'slmsk in rad_prepare')
#ifdef RMP
     call print_maxmin_six(coszer,LATG2S*LONF2S,1,1,1,'coszer in rad_prepare')
#else
     call print_maxmin_six(albedr,LATG2S*LONF2S,1,1,1,'albedr in rad_prepare')
#endif
     call print_maxmin_six(alvbr,LATG2S*LONF2S,1,1,1,'alvbr in rad_prepare')
     call print_maxmin_six(alnbr,LATG2S*LONF2S,1,1,1,'alnbr in rad_prepare')
     call print_maxmin_six(alvdr,LATG2S*LONF2S,1,1,1,'alvdr in rad_prepare')
     call print_maxmin_six(alndr,LATG2S*LONF2S,1,1,1,'alndr in rad_prepare')
   endif
#endif
!
   tsmin=tsea(1,1)
   tsmax=tsea(1,1)
   shmin=snoweq(1,1)
   shmax=snoweq(1,1)
!
   do j = jts,jte
     do i = its,ite
       cvr  (i,j)=cv   (i,j)
       cvtr (i,j)=cvt  (i,j)
       cvbr (i,j)=cvb  (i,j)
     enddo
   enddo
!
   return
   end subroutine rad_prepare
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_julian_day(jyr,jmnth,jday,jhr,jmn,jd,fjd)
!-------------------------------------------------------------------------------
!
! subroutine: rad_julian_day
!
! abstract: this code written at gfdl ....
!   computes julian day and fraction
!   from year, month, day and time ut...accurate only between
!   march 1, 1900 and february 28, 2100.. based on julian calendar
!   corrected to correspond to gregorian calendar during this period.
!
! program history log:
!   1977-05-06  ray orzol,gfdl
!   1989-07-07  kenneth campana
! 
! usage:    call rad_julian_day(jyr,jmnth,jday,jhr,jmn,jd,fjd)
!   input argument list:
!     jyr      - year (4 digits)-intial fcst time.
!     jmnth    - month-initial fcst time.
!     jday     - day-initial fcst time.
!     jhr      - z-time of initial fcst time.
!     jmn      - minutes (zero passed from calling program).
!   output argument list:
!     jd       - julian day.
!     fjd      - fraction of the julian day.
!
!-------------------------------------------------------------------------------
!
!    *******************************************************************
!    *                           c o m p j d                           *
!    *    statement blocked by ray orzol                               *
!    *******************************************************************
!
   integer              ::  ndm(12)
!
   data                                                                        &
      jdor/2415019/,                                                           &
      jyr19/1900/
!                                                                               
   data  ndm/0,31,59,90,120,151,181,212,243,273,304,334/
!
!    *******************************************************************
!     computes julian day and fraction from year, month, day and time ut
!     accurate only between march 1, 1900 and february 28, 2100
!     based on julian calendar corrected to correspond to gregorian
!        calendar during this period
!    jdor=jd of december 30, 1899 at 12 hours ut
!    *******************************************************************
!
   jd=jdor
   jyrm9=jyr-jyr19
   lp=jyrm9/4
   if(lp.le.0) go to 4
   jd=jd+1461*lp
!
4  ny=jyrm9-4*lp
!
   ic=0
   if(ny.gt.0) go to 5
   if(jmnth.gt.2) ic=1
   go to 6
!
5  jd=jd+365*ny+1
!
6  jd=jd+ndm(jmnth)+jday+ic
!
   if(jhr.ge.12) go to 7
   jd=jd-1
   fjd=.5e0+.041666667e0*float(jhr)+.00069444444e0*float(jmn)
!
   return
!
7  fjd=.041666667e0*float(jhr-12)+.00069444444e0*float(jmn)
!
   return
   end subroutine rad_julian_day
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_fcst_time(fhour,imon,iday,iztim,jdnmc,fjdnmc,rlag,year,      &
                     rsin1,rcos1,rcos2,jd,fjd)
!-------------------------------------------------------------------------------
   use paramodel
   use comio
#ifdef DFS
   use dfsvar, only   : iope
#endif
   use constant, only : pi_
!-------------------------------------------------------------------------------
   real,parameter       ::  tpi=2.e0*pi_
   integer,save         ::  jmon(12)
   real,save            ::  two
   data jmon/31,28,31,30,31,30,31,31,30,31,30,31/
   data two/2.e0/
!
!....    first get number of days since beginning of year (no leap yrs)
!
   nnday =0
   imo = imon - 1
   if (imo.gt.0) then
     do i = 1,imo
       nnday = nnday + jmon(i)
     enddo
   endif
   nnday = nnday + iday
#ifndef RMP
#ifndef NOPRINT
   if(iope) write(6,1002)nnday
 1002 format(1h ,'*************** nnday of year = ',i4,'******')
#endif
#endif
!
!....    get number of days into fcst (dyfcst)
!....   following two cards changed on 10 apr 86 to fix slight error
!         in solar declination calc if initial hr not 00z or 12z.....
!
   dayini = nnday + float(iztim)/24.e0
   soltim = fhour + iztim
!
!...     reset to 24 hour clock
!
   fday = soltim / 24.e0
!>yh  soltim = soltim - int(fday) * 24.e0
   dyfcst = nnday + fday
#ifndef RMP
#ifdef DBG
   if(iope) write(6,1003)jdnmc,fjdnmc,fhour,dayini,dyfcst
 1003 format(2x,'jdnmc etc in rad_setup ',i9,2x,4(2x,f6.2))
#endif
#endif
   rang=tpi*(dyfcst-rlag)/year
   rsin1=sin(rang)
   rcos1=cos(rang)
   rcos2=cos(two*rang)
!
!....    update the julian date (initial in jdnmc,fjdnmc)
!
   dyinc = dyfcst - dayini
   idyin = dyinc
   fdyin = dyinc - idyin
   jd = jdnmc + idyin
   fjd = fjdnmc + fdyin
!
!.......need to reset if fraction (fjd) gt 1.
!
   ifjd = fjd
   if (ifjd.gt.0) then
     jd = jd + ifjd
     fjd = fjd - ifjd
   endif
!
   return
   end subroutine rad_fcst_time
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_solar_time(jd,fjd,r,dlt,alp,slag,sdec,cdec)
!-------------------------------------------------------------------------------
!
! subroutine:  rad_solar_time
!
! abstract: this code written at gfdl ....
!   - astronomical(solar) data - sw radiation
!   - computes radius vector,declination and right ascension of sun,
!     and equation of time. same as subroutine -solmrf-, but
!     with hour angle,fractional daylight,and
!     mean zenith angle calculations removed (subroutine -zenith-
!     calculates these for each point rather than each latitude).
!   --this code is to be used for other fcst models or for the mrf
!   --model if 'instantaneous' sw calculations desired.
!
! program history log:
!   77-07-21  robert white,gfdl.
!   89-07-07  kenneth campana-moved the hour angle calculations to
!                              subroutine -zenith-
!
! usage:    call solmrf(jd,fjd,r,dlt,alp,slag)
!   input argument list:
!     jd       - julian day for current fcst hour.
!     fjd      - fraction of the julian day.
!   output argument list:
!     r        - radius vector of the sun.
!     dlt      - declination of the sun (radians).
!     alp      - right ascension of the sun.
!     slag     - equation of time (radians).
!
!   output files:
!     output   - print file.
!
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_
!-------------------------------------------------------------------------------
   real,parameter       ::  tpi=2.0*pi,hpi=0.5*pi,rad=180.0/pi
!
!    *******************************************************************
!    *                            s o l a r                            *
!... *  patterned after original gfdl code---                          *
!... *     but no calculation of latitude mean cos solar zenith angle..*
!... *     zenith angle calculations done in subr zenith in this case..*
!... *  hr angle,mean cosz,and mean tauda calc removed--k.a.c. mar 89  *
!    *  updates by hualu pan to limit iterations in newton method and  *
!    *  also ccr reduced from(1.3e-7)--both to avoid nonconvergence in *
!    *  nmc s half precision version of gfdl s code   ----  fall 1988  *
!    *******************************************************************
!
!.....rad_solar_time computes radius vector, declination and right ascension of
!.....sun, equation of time
!
   data cyear/365.25/,      ccr/1.3e-6/
!
!.....tpp = days between epoch and perihelion passage of 1900
!.....svt6 = days between perihelion passage and march equinox of 1900
!.....jdor = jd of epoch which is january 0, 1900 at 12 hours ut
!
   data tpp/1.55/,          svt6/78.035/,       jdor/2415020/
!
!    *******************************************************************
!
   dat=float(jd-jdor)-tpp+fjd
!
!    computes time in julian centuries after epoch
!
   t=float(jd-jdor)/36525.e0
!
!    computes length of anomalistic and tropical years (minus 365 days)
!
   year=.25964134e0+.304e-5*t
   tyear=.24219879e0-.614e-5*t
!
!    computes orbit eccentricity and angle of earths inclination from t
!
   ec=.01675104e0-(.418e-4+.126e-6*t)*t
   angin=23.452294e0-(.0130125e0+.164e-5*t)*t
   ador=jdor
   jdoe=ador+(svt6*cyear)/(year-tyear)
!
!    deleqn=updated svt6 for current date
!
   deleqn=float(jdoe-jd)*(year-tyear)/cyear
   year=year+365.e0
   sni=sin(angin/rad)
   tini=1.e0/tan(angin/rad)
   er=sqrt((1.e0+ec)/(1.e0-ec))
   qq=deleqn*tpi/year
!
!    determine true anomaly at equinox
!
   e=1.e0
   iter = 0
32 ep=e-(e-ec*sin(e)-qq)/(1.e0-ec*cos(e))
   cd=abs(e-ep)
   e=ep
   iter = iter + 1
#ifndef NOPRINT
   if(iter.gt.10) then
     write(6,*) ' iteration count for loop 32 =', iter
     write(6,*) ' e, ep, cd =', e, ep, cd
   endif
#endif
   if(iter.gt.10) goto 1032
   if(cd.gt.ccr) go to 32
!
1032 continue
!
   he=.5e0*e
   eq=2.e0*atan(er*tan(he))
!
!    date=days since last perihelion passage
!
   date = mod(dat,year)
!
!    solve orbit equations by newtons method
!
   em=tpi*date/year
   e=1.e0
   iter = 0
!
31 ep=e-(e-ec*sin(e)-em)/(1.e0-ec*cos(e))
!
   cr=abs(e-ep)
   e=ep
   iter = iter + 1
#ifndef NOPRINT
   if(iter.gt.10) then
      write(6,*) ' iteration count for loop 31 =', iter
   endif
#endif
   if(iter.gt.10) goto 1031
   if(cr.gt.ccr) go to 31
!
1031 continue
!
   r=1.e0-ec*cos(e)
   he=.5e0*e
   w=2.e0*atan(er*tan(he))
! yh  sind=sni*sin(w-eq)                                                        
! yh  dlt=asin(sind)                                                            
   sdec=sni*sin(w-eq)
   cdec=sqrt(1.e0 - sdec*sdec)
   dlt=asin(sdec)
   alp=asin(tan(dlt)*tini)
   tst=cos(w-eq)
   if(tst.lt.0.e0) alp=pi-alp
   if(alp.lt.0.e0) alp=alp+tpi
   sun=tpi*(date-deleqn)/year
   if(sun.lt.0.e0) sun=sun+tpi
   slag=sun-alp-.03255e0
!
   return
   end subroutine rad_solar_time
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_cos_zenith(imx2,jmx2,                                        &
                             dtswav,solhr,sinlat,coslat,sdec,cdec,slag,        &
                             xlon,coszen,ldg,coszdg)                       
!-------------------------------------------------------------------------------
! 
! abstract :
!   compute mean cos solar zen angl over dtswav hrs
!   cosine of solar zen angl for both n. and s. hemispheres.
!   solhr=time(hrs) after 00z (greenwich time)..
!   xlon is east long(radians)..
!   sinlat, coslat are sin and cos of latitude (n. hemisphere)
!   sdec, cdec = the sine and cosine of the solar declination.
!   slag = equation of time
!
!-------------------------------------------------------------------------------
   use paramodel, only   :  ILOTS
!-------------------------------------------------------------------------------
   real                 ::  xlon(imx2,jmx2),coszen(imx2,jmx2)
   logical              ::  ldg
   real                 ::  coszdg(imx2,jmx2)
   real                 ::  sinlat(imx2,jmx2),coslat(imx2,jmx2)
!
! local array
!
   real                 ::  coszn(ILOTS)
   integer              ::  istsun(ILOTS)
!-------------------------------------------------------------------------------
   nstp = 6
   istp = nstp*dtswav
   pid12 = (2.e0 * asin(1.e0)) / 12.e0
!
   do j = 1,jmx2
     do i = 1,imx2
       coszen(i,j) = 0.e0
       istsun(i) = 0
     enddo
     do it = 1,istp
       cns = pid12 * (solhr-12.e0+(it-1)*1.e0/nstp) +slag
       do i = 1,imx2
         ss=sinlat(i,j)*sdec
         cc=coslat(i,j)*cdec
         coszn(i) = ss + cc * cos(cns + xlon(i,j))
         coszen(i,j) = coszen(i,j) + max (0.e0, coszn(i))
         if(coszn(i).gt.0.e0) istsun(i) = istsun(i) + 1
       enddo
     enddo
     do i = 1,imx2
       if(ldg) coszdg(i,j) = coszen(i,j) / istp
       if(istsun(i).gt.0) coszen(i,j) = coszen(i,j) / istsun(i)
     enddo
   enddo
!
   return
   end subroutine rad_cos_zenith
!-------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
   subroutine rad_calc_time(jd,fjd,munth,im,id,iyear,ihr,xmin)
!------------------------------------------------------------------------------
!
! subprogram:    rad_calc_time
!
! abstract: this code written at gfdl ....
!   computes month,day,year from julian day.
!   accurate only between march 1, 1900 and february 28, 2100....
!   based on julian calender corrected to correspond to gregorian
!   calender during this period.
!
! program history log:
!   1977-06-07  robert white,gfdl
!   1989-07-07  kenneth campana
!
! usage:    call rad_calc_time(jd,fjd,munth,im,id,iyear,ihr,xmin)
!   input argument list:
!     jd       - julian day for current fcst hour.
!     fjd      - fraction of the julian day.
!   output argument list:
!     munth    - month (character).
!     im       - month (integer).
!     id       - day of the month.
!     iyear    - year.
!     ihr      - hour of the day.
!     xmin     - minute of the hour.
!-------------------------------------------------------------------------------
   real                 ::  dy(13)
   integer              ::  month(12)
!
   data dy/                                                                    &
      0.,                 31.,                59.,                             &
      90.,                120.,               151.,                            &
      181.,               212.,               243.,                            &
      273.,               304.,               334.,                            &
      365.                /
!
   data month/                                                                 &
      4hjan.,             4hfeb.,             4hmar.,                          &
      4hapr.,             4hmay ,             4hjune,                          &
      4hjuly,             4haug.,             4hsep.,                          &
      4hoct.,             4hnov.,             4hdec.                           &
      /
!
!.....jdor = jd of december 30, 1899 at 12 hours ut
!
   integer              ::  jdor,iyr
   data                                                                        &
      jdor                /         2415019             /,                     &
      iyr                 /         1900                /
!                                                                               
   iyear=iyr
   nday=jd-jdor
   if(fjd.ge..5e0) nday=nday+1
!
61 if(nday.lt.1462) go to 62
!
   nday=nday-1461
   iyear=iyear+4
!
   go to 61
!
62 ndiy=365
!
   if(mod(iyear,4).eq.0) ndiy=366
   if(nday.le.ndiy) go to 65
   iyear=iyear+1
   nday=nday-ndiy
   go to 62
!
65 if(nday.gt.int(dy(2))) go to 66
!
   im=1
   id=nday
   go to 67
!
66 if(nday.ne.60) go to 68
!
   if(ndiy.eq.365) go to 68
   im=2
   id=29
   go to 67
!
68 if(nday.gt.(int(dy(3))+ndiy-365)) go to 69
!
   im=2
   id=nday-31
   go to 67
!
69 do 70 i = 3,12
!
     if(nday.gt.(int(dy(i+1))+ndiy-365)) go to 70
     im=i
     id=nday-int(dy(i))-ndiy+365
     go to 67
70 continue
67 munth=month(im)
   hr=24.e0*fjd
   ihr=hr
   xmin=60.e0*(hr-float(ihr))
   ihr=ihr+12
   if(ihr.ge.24) ihr=ihr-24
!
   return
   end subroutine rad_calc_time
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_print_time(id,munth,iyear,ihr,xmin,jd,fjd,                   &
                             dlt,alf,r1,slag,solc)
!-------------------------------------------------------------------------------
   use comio
   use constant, only : pi_
#ifdef DFS
   use dfsvar, only : iope
#endif
!-------------------------------------------------------------------------------
   real,parameter       ::  degrad=180.e0/pi_,hpi=0.5e0*pi_
   real,save            ::  sign,zero,six,sixty,q22855
   data       sign/1h-/,      sigb/1h /
   data zero,six,sixty,q22855/0.0,6.0,60.0,228.55735/
!
   dltd=degrad*dlt
   ltd=dltd
   dltm=sixty*(abs(dltd)-abs(float(ltd)))
   ltm=dltm
   dlts=sixty*(dltm-float(ltm))
   dsig=sigb
   if((dltd.lt.zero).and.(ltd.eq.0)) dsig=sign
   halp=six*alf/hpi
   ihalp=halp
   ymin=abs(halp-float(ihalp))*sixty
   iyy=ymin
   asec=(ymin-float(iyy))*sixty
   eqt=q22855*slag
   eqsec=sixty*eqt
!
#ifndef RMP
#ifndef NOPRINT
   if(iope) write(6,1004)                                                      &
           id,munth,iyear,ihr,xmin,jd,fjd,r1,halp,ihalp,                       &
          iyy,asec,dltd,dsig,ltd,ltm,dlts,eqt,eqsec,slag,solc
 1004 format('  forecast date',9x,i3,a5,i6,' at',i3,' hrs',f6.2,' mins'/       &
          '  julian day',12x,i8,2x,'plus',f11.6/                               &
          '  radius vector',9x,f10.7/                                          &
          '  right ascension of sun',f12.7,' hrs, or',i4,' hrs',i4,            &
                                    ' mins',f6.1,' secs'/                      &
          '  declination of the sun',f12.7,' degs, or',a2,i3,                  &
                                    ' degs',i4,' mins',f6.1,' secs'/           &
          '  equation of time',6x,f12.7,' mins, or',f10.2,' secs, or'          &
                              ,f9.6,' radians'/                                &
          '  solar constant',8x,f12.7)
#endif
#endif
!
   return
   end subroutine rad_print_time
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_albedo_aerosol(imx2,jmx2,                                    &
                     mon,slmsk,snowf,z0cmf,coszf,tseaf,hprif,jsno,             &
#ifndef SWRMDC
                     alvsf,alnsf,alvwf,alnwf,facsf,facwf,paerf,                &
                     xlat,alvbr,alnbr,alvdr,alndr,latdef,paerr)
#else
                     alvsf,alnsf,alvwf,alnwf,facsf,facwf,                      &
                     xlat,alvbr,alnbr,alvdr,alndr,latdef,                      &
                     kprfi,idxci,cmixi,denni,                                  &
                     iaer,nfaer,kprfg,idxcg,cmixg,denng)
#endif
!-------------------------------------------------------------------------------
!
! abstract:
!   this program computes four components of surface albedos (i.e.
!     vis-nir, direct-diffused) based on brieglebs scheme. and
!     bilinearly interpolates albedo and aerosol distribution to
!     radiation grid.
!
! program history log:
!   1999-09-01  yu-tai hou             updates to opac aerosol algorithm
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! input variables:
!     mon     - month of the year
!     slmsk   - sea(0),land(1),ice(2) mask on fcst model grid
!     snowf   - snow depth water equivalent in mm
!     z0cmf   - surface roughness in cm
!     coszf   - cosin of solar zenith angle
!     tseaf   - sea surface temperature in k
!1197 hprif   - topographic sdv in m
!1197 jsno    - lat at 70 deg - indicating extent of perm snow cover
!     alvsf   - mean vis albedo with strong cosz dependency
!     alnsf   - mean nir albedo with strong cosz dependency
!     alvwf   - mean vis albedo with weak cosz dependency
!     alnwf   - mean nir albedo with weak cosz dependency
#ifndef SWRMDC
!     paerf   - aerosol distribution factor on fcst grid
#else
!     nfaer   - file number for input aerosol data
#endif
!
! output variables:   (all on radiation grid)
!     alvbr   - vis beam surface albedo
!     alnbr   - nir beam surface albedo
!     alvdr   - vis diff surface albedo
!     alndr   - nir diff surface albedo
#ifndef SWRMDC
!     paerr   - aerosol distribution factor
#else
!     kprfg     - aerosol profile type index
!     idxcg     - aerosol components types indices
!     cmixg     - aerosol components mixing ratioes
!     denng     - first two layers aerosol number densities
#endif
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S
   use constant, only : t0c_,pi_
   use comio
#ifdef DFS
   use dfsvar, only : iope
#endif
#ifdef SWRMDC
   use aerparm
#endif
#include "abort.h"
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  slmsk(imx2,jmx2),snowf(imx2,jmx2),z0cmf(imx2,jmx2)
   real                 ::  tseaf(imx2,jmx2),coszf(imx2,jmx2),hprif(imx2,jmx2)
   real                 ::  alvsf(imx2,jmx2), alnsf(imx2,jmx2)
   real                 ::  alvwf(imx2,jmx2), alnwf(imx2,jmx2)
   real                 ::  facsf(imx2,jmx2),   facwf(imx2,jmx2)
#ifndef SWRMDC
   real                 ::  paerf(imx2,jmx2,5)
#endif
   real                 ::  xlat(imx2,jmx2)
!
! --- output
!
   real                 ::  alvbr(imx2,jmx2),alnbr(imx2,jmx2),alvdr(imx2,jmx2)
#ifndef SWRMDC
#ifdef RMPVECTORIZE
   real                 ::  alndr(imx2,jmx2),paerr(imx2,jmx2,5)
#else
   real                 ::  alndr(imx2,jmx2),paerr(imx2,5,jmx2)
#endif
#else
   real                 ::  alndr(imx2,jmx2)
   real                 ::  cmixi(imx2,jmx2,nxc)
   real                 ::  denni(imx2,jmx2,ndn)
   real                 ::  idxci(imx2,jmx2,nxc)
   real                 ::  kprfi(imx2,jmx2)
   real                 ::  cmixg(nxc,imx2,jmx2)
   real                 ::  denng(ndn,imx2,jmx2)
   integer              ::  idxcg(nxc,imx2,jmx2)
   integer              ::  kprfg(imx2,jmx2)
#endif
   integer              ::  latdef(jmx2)
!
! --- internal variables
!
   real                 ::  alvbf(LONF2S,LATG2S)
   real                 ::  alnbf(LONF2S,LATG2S)
   real                 ::  alvdf(LONF2S,LATG2S)
   real                 ::  alndf(LONF2S,LATG2S)
   real                 ::  asnvb(LONF2S),asnnb(LONF2S)
   real                 ::  asnvd(LONF2S),asnnd(LONF2S)
   real                 ::  asevb(LONF2S),asenb(LONF2S)
   real                 ::  asevd(LONF2S),asend(LONF2S)
   real                 ::  fsno (LONF2S),fsea (LONF2S)
   real                 ::  rfcs (LONF2S),rfcw (LONF2S)
#ifndef SWRMDC
   real                 ::  flnd (LONF2S)
   integer              ::  mm(4)
#else
   real                 ::  flnd (LONF2S)
   integer              ::  mm(4)
   integer              ::  idxc (nxc)
   real                 ::  cmix (nxc)
   character            ::  cline*80, ctyp*3, aerosol_file*40
   integer, save        ::  mon_save
   data                     mon_save/0/
#endif
!
   logical              ::  snochk
   real,parameter       ::  snodeg=70.* pi_ /180.
!
   ifcs=imx2
   jfcs=jmx2
   irad=imx2
   jrad=jmx2
!
#ifdef SWRMDC
! ==============================================
! ===  first section defines surface albedo  ===
! ==============================================
#endif
   do j = 1,jfcs
#if defined(SMP)
     nxlat = nint(xlat(1,j)*180./pi_)
     jj = 90 - nxlat
     if (jj.gt.90) jj = 180 - jj
#else
     jj=latdef(j)
#endif
!
     do i = 1,ifcs
       snochk = .false.
       if(jj.le.jsno.or.(jsno.eq.0.and.abs(xlat(i,j)).gt.snodeg))              &
                     snochk = .true.
!
! --- modified snow albedo scheme - units convert to m
!     (originally snowf in mm; z0cmf in cm)
!
       asnow = 0.02*snowf(i,j)
       argh  = min(1.0, max(.025, 0.01*z0cmf(i,j)))
       fsno0 = asnow / (argh + asnow)
       if (slmsk(i,j).eq.0.0 .and. tseaf(i,j).gt.271.2) fsno0 = 0.0
       fsno1 = 1.0 - fsno0
       flnd0 = min(1.0,facsf(i,j) + facwf(i,j))
       fsea0 = max(0.0, 1.0 - flnd0)
       fsno (i) = fsno0
       fsea (i) = fsea0 * fsno1
       flnd (i) = flnd0 * fsno1
!
! --- diffused sea surface albedo
!
       if (slmsk(i,j).eq.2.0 .or. tseaf(i,j) .lt. 271.2) then
         asevd(i) = 0.70
         asend(i) = 0.65
       else
         asevd(i) = 0.06
         asend(i) = 0.06
       endif
!
! --- diffused snow albedo
!
       if (slmsk(i,j).eq.1.0 .and. ((.not.snochk).or.tseaf(i,j).gt.271.2)) then
           hfac = max(0.10, min(1.0, 1.225-1.125e-3*hprif(i,j)))
       else
           hfac = max(0.75, min(1.0, 1.0625-0.3125e-3*hprif(i,j)))
       endif
       asnvd(i) = 0.90 * hfac
       asnnd(i) = 0.75 * hfac
     enddo
!
     do i = 1,ifcs
!
! --- direct snow albedo
!
       if (coszf(i,j) .lt. 0.5) then
         csnow = 0.5 * (3.0 / (1.0+4.0*coszf(i,j)) - 1.0)
         asnvb(i) = min( 0.98, asnvd(i)+(1.0-asnvd(i))*csnow )
         asnnb(i) = min( 0.98, asnnd(i)+(1.0-asnnd(i))*csnow )
       else
         asnvb(i) = asnvd(i)
         asnnb(i) = asnnd(i)
       endif
!
! --- direct sea surface albedo
!
       if (coszf(i,j) .gt.0.0) then
         rfcs(i) = 1.4 / (1.0 + 0.4*coszf(i,j))
         rfcw(i) = 1.1 / (1.0 + 0.2*coszf(i,j))
         if (tseaf(i,j) .ge. t0c_) then
           asevb(i) = max(0.055, 0.026/(coszf(i,j)**1.7+0.065)                 &
                    + 0.15 * (coszf(i,j)-0.1) * (coszf(i,j)-0.5)               &
                    * (coszf(i,j)-1.0))
           asenb(i) = asevb(i)
         else
           asevb(i) = asevd(i)
           asenb(i) = asend(i)
         endif
       else
         rfcs(i) = 1.0
         rfcw(i) = 1.0
         asevb(i) = asevd(i)
         asenb(i) = asend(i)
       endif
     enddo
     do i = 1,ifcs
       a1   = alvsf(i,j) * facsf(i,j)
       b1   = alvwf(i,j) * facwf(i,j)
       a2   = alnsf(i,j) * facsf(i,j)
       b2   = alnwf(i,j) * facwf(i,j)
       alvbf(i,j) = (a1*rfcs(i) + b1*rfcw(i))*flnd(i)                          &
                  + asevb(i)*fsea(i) + asnvb(i)*fsno(i)
       alvdf(i,j) = (a1         + b1        )*flnd(i)                          &
                  + asevd(i)*fsea(i) + asnvd(i)*fsno(i)
       alnbf(i,j) = (a2*rfcs(i) + b2*rfcw(i))*flnd(i)                          &
                  + asenb(i)*fsea(i) + asnnb(i)*fsno(i)
       alndf(i,j) = (a2         + b2        )*flnd(i)                          &
                  + asend(i)*fsea(i) + asnnd(i)*fsno(i)
     enddo
!
   enddo
!
! ... for one grid, no interpolation ...
!
   do j = 1,jfcs
     do i = 1,ifcs
       alvbr(i,j)=alvbf(i,j)
       alvdr(i,j)=alvdf(i,j)
       alnbr(i,j)=alnbf(i,j)
       alndr(i,j)=alndf(i,j)
     enddo
   enddo
#ifndef SWRMDC
!
   do k = 1,5
     do j = 1,jfcs
       do i = 1,ifcs
#ifndef RMPVECTORIZE
         paerr(i,k,j)=min(1.0,paerf(i,j,k))
         if (paerr(i,k,j).lt. 0.01) paerr(i,k,j) = 0.0
#else
         paerr(i,j,k)=min(1.0,paerf(i,j,k))
         if (paerr(i,j,k).lt. 0.01) paerr(i,j,k) = 0.0
#endif
       enddo
     enddo
   enddo
!
!....  final check to make total is one
!
   do j = 1,jrad
     do i = 1,irad
       if (slmsk(i,j).eq.0.0 .or. slmsk(i,j).eq.2.0) then
#ifndef RMPVECTORIZE
         psea = paerr(i,3,j) + paerr(i,5,j)
         paerr(i,1,j) = 0.0
         paerr(i,2,j) = 0.0
         paerr(i,4,j) = 0.0
         if (psea .gt. 1.0) then
           paerr(i,3,j) = min(1.0, paerr(i,3,j))
           paerr(i,5,j) = 1.0 - paerr(i,3,j)
         else if (psea .lt. 1.0) then
           paerr(i,3,j) = 1.0 - paerr(i,5,j)
         endif
       else
         plnd = paerr(i,1,j) + paerr(i,2,j) + paerr(i,4,j)
         paerr(i,5,j) = 0.0
         if (plnd .gt. 1.0) then
           paerr(i,3,j) = 0.0
           paerr(i,2,j) = min(1.0, paerr(i,2,j))
           if (paerr(i,1,j) .gt. 0.0) then
             paerr(i,4,j) = 0.0
             paerr(i,1,j) = 1.0 - paerr(i,2,j)
           else
             paerr(i,1,j) = 0.0
             paerr(i,4,j) = 1.0 - paerr(i,2,j)
           endif
         else if (plnd .lt. 1.0) then
!
! --- use mar-i as the background fill
!
           paerr(i,3,j) = 1.0 - plnd
#else
         psea = paerr(i,j,3) + paerr(i,j,5)
         paerr(i,j,1) = 0.0
         paerr(i,j,2) = 0.0
         paerr(i,j,4) = 0.0
         if (psea .gt. 1.0) then
           paerr(i,j,3) = min(1.0, paerr(i,j,3))
           paerr(i,j,5) = 1.0 - paerr(i,j,3)
         else if (psea .lt. 1.0) then
           paerr(i,j,3) = 1.0 - paerr(i,j,5)
         endif
       else
         plnd = paerr(i,j,1) + paerr(i,j,2) + paerr(i,j,4)
         paerr(i,j,5) = 0.0
         if (plnd .gt. 1.0) then
           paerr(i,j,3) = 0.0
           paerr(i,j,2) = min(1.0, paerr(i,j,2))
           if (paerr(i,j,1) .gt. 0.0) then
             paerr(i,j,4) = 0.0
             paerr(i,j,1) = 1.0 - paerr(i,j,2)
           else
             paerr(i,j,1) = 0.0
             paerr(i,j,4) = 1.0 - paerr(i,j,2)
           endif
         else if (plnd .lt. 1.0) then
!
! --- use mar-i as the background fill
!
           paerr(i,j,3) = 1.0 - plnd
#endif
         endif
       endif
     enddo
   enddo
#else /* SWRMDC */
! 
   if (iaer .eq. 1) then
     do j = 1,jrad
       do i = 1,irad
         kprfg(i,j) = nint(kprfi(i,j))
       enddo
     enddo
!
     do j = 1,jrad
       do i = 1,irad
         do k = 1,2
           denng(k,i,j) = denni(i,j,k)
         enddo
       enddo
     enddo
!
     do j = 1,jrad
       do i = 1,irad
         do k = 1,nxc
           idxcg(k,i,j) = nint(idxci(i,j,k))
           cmixg(k,i,j) = cmixi(i,j,k)
         enddo
       enddo
     enddo
   endif
#endif /* ~SWRMDC end */
!
   return
   end subroutine rad_albedo_aerosol
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_cloud_interp(cvin,cvtin,cvbin,iin,jtwidl,jin,                &
                     cvout,cvtout,cvbout,iout,jpout,jout,                      &
                     xx,wgt,tt,bb,sum,nn,                                      &
                     ltwidl,latrd1,latinb)
!-------------------------------------------------------------------------------
!
!     *  code bilinearly interpolates cld amt between gaussian grids--*
!     *  clone of ggintp for interpolation of convective cld amt (cv).*
!     *    special interp procedure for tops(cvt) and bots(cvb)...    *
!-    *  j = 1 is just belo n.pole, i = 1 is greenwich (then go east).*
!     * iin,jin are i,j dimensions of input grid--iout,jout for output*
!     * jin2,jout2=jin/2,jout/2                                       *
!     *                                     --k.campana - june 1988   *
!
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,jpout,jout,ltwidl,latrd1,latinb
   real                 ::  cvin(iin,jtwidl),cvtin(iin,jtwidl)
   real                 ::  cvbin(iin,jtwidl)
   real                 ::  cvout(iout,jpout)
   real                 ::  cvtout(iout,jpout),cvbout(iout,jpout)
   real                 ::  xx(iout,4),wgt(iout,4),tt(iout,4)
   real                 ::  bb(iout,4),sum(iout,4)
   integer              ::  nn(iout)
!
   iii = iin
   jbb = jtwidl
   jjj = jin
   iiiout = iout
   lbb = ltwidl
   lr1 = latrd1
!
   do latout = 1,jpout
     lat=latout+latinb-1
     if(lat.eq.1) then
       inslat=-1
       wgtlat=0.
     else
       inslat=lat-1
       wgtlat=1.
     endif
!
!===>    if output lat is poleward of input lat=1 ,then simpl average
!          (small region and cld amt wouldn t extrapolate well)
!
     call rad_cloud_interp_sub(iii,jbb,jjj,iiiout,inslat,wgtlat,               &
               cvin,cvtin,cvbin,cvout(1,latout),                               &
               cvtout(1,latout),cvbout(1,latout),                              &
               xx,wgt,tt,bb,sum,nn,lbb,lr1)
   enddo
!
! 100 format(1h ,' row =',i5,'  lat =',e15.5)
!
   return
   end subroutine rad_cloud_interp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_cloud_interp_sub(iin,jtwidl,jin,iout,                        &
                    inslat,wgtlat,                                             &
                    cv,cvt,cvb,camt,ctop,cbot,                                 &
                    xx,wgt,tt,bb,sum,nn,ltwidl,latrd1)
!-------------------------------------------------------------------------------
!
! subroutine: rad_cloud_interp_sub
!
!      simpl linear interpolation of cldamt, unless only 1,2 of the
!      surrounding pts has cv. then,if output gridpt not close enuf
!         do not interpolate to it(prevents spreading of cv clds)..
!           for 1 pt convection-intrp wgt ge (.7)**2 ...
!           for 2 pt convection-sum of intrp wgt ge .45...
!              .45 used rather than .5 to give better result for
!              diagonally opposed pts...
!===>    for tops(cvt) and bots(cvb) just take average of surrounding
!         non-zero cv points.....
!         nn will be number of surrounding pts with cld (gt zero)
!---     nhsh = 1,-1 for northern,southern hemisphere
!         here instead of an extrapolation,just do a simple mean....
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,ltwidl,latrd1,inslat
   real                 ::  wgtlat
   real                 ::  cv(iin,jtwidl),cvt(iin,jtwidl)
   real                 ::  cvb(iin,jtwidl)
   real                 ::  camt(iout),ctop(iout),cbot(iout)
   real                 ::  xx(iout,4),wgt(iout,4),tt(iout,4)
   real                 ::  bb(iout,4),sum(iout,4)
   integer              ::  nn(iout)
!
! local
!
   integer              ::  i,j,kpt,lonf,inth,inth1,ileft,irght,iout2
   integer              ::  ileft2,irght2,ja,ltop,lbot
   real                 ::  wgtlon
!
   lonf=iin/2
   if (inslat.lt.0) go to 600
   inth = mod(ltwidl + inslat + jtwidl - latrd1 - 1,jtwidl) + 1
   inth1 = mod(inth,jtwidl) + 1
   if (inslat.eq.jin) go to 105
!
   do i = 1,iout
     ileft=i
     irght=i+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     wgtlon=0.0
     wgtlat=1.0
!
!----   normalized distance from upper lat to gaussian lat
!
     xx(i,1) = cv(ileft,inth)
     xx(i,2) = cv(ileft,inth1)
     xx(i,3) = cv(irght,inth)
     xx(i,4) = cv(irght,inth1)
     wgt(i,1) = (1.e0-wgtlon)*(1.e0-wgtlat)
     wgt(i,2) = (1.e0-wgtlon)*wgtlat
     wgt(i,3) = wgtlon*(1.e0-wgtlat)
     wgt(i,4) = wgtlon*wgtlat
     tt(i,1) = cvt(ileft,inth)
     tt(i,2) = cvt(ileft,inth1)
     tt(i,3) = cvt(irght,inth)
     tt(i,4) = cvt(irght,inth1)
     bb(i,1) = cvb(ileft,inth)
     bb(i,2) = cvb(ileft,inth1)
     bb(i,3) = cvb(irght,inth)
     bb(i,4) = cvb(irght,inth1)
   enddo
!
   go to 130
!
105 continue
!
   do i = 1,iout
     ileft=i
     irght=i+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     wgtlon=0.0
     wgtlat=1.0
!
!----   normalized distance from upper lat to gaussian lat
!
     xx(i,1) = cv(ileft,inth)
     xx(i,3) = cv(irght,inth)
     wgt(i,1) = (1.e0-wgtlon)*(1.e0-wgtlat)
     wgt(i,2) = (1.e0-wgtlon)*wgtlat
     wgt(i,3) = wgtlon*(1.e0-wgtlat)
     wgt(i,4) = wgtlon*wgtlat
     tt(i,1) = cvt(ileft,inth)
     tt(i,3) = cvt(irght,inth)
     bb(i,1) = cvb(ileft,inth)
     bb(i,3) = cvb(irght,inth)
   enddo
!
   iout2 = iout / 2
   do i = 1,iout2
     ileft=i
     irght=ileft+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     ileft2=i+iout2
     irght2=ileft2+1
     if(mod(ileft2,lonf).eq.0) irght2=irght2-lonf
     xx(i      ,2) = cv(ileft2,inth)
     xx(i+iout2,2) = cv(ileft ,inth)
     xx(i      ,4) = cv(irght2,inth)
     xx(i+iout2,4) = cv(irght ,inth)
     bb(i      ,2) = cvb(ileft2,inth)
     bb(i+iout2,2) = cvb(ileft ,inth)
     bb(i      ,4) = cvb(irght2,inth)
     bb(i+iout2,4) = cvb(irght ,inth)
     tt(i      ,2) = cvt(ileft2,inth)
     tt(i+iout2,2) = cvt(ileft ,inth)
     tt(i      ,4) = cvt(irght2,inth)
     tt(i+iout2,4) = cvt(irght ,inth)
   enddo
!
!---      nn will be number of surrounding pts with cld (gt zero)
!
130 continue
!
   do i = 1,iout
     nn(i) = 0                                                               
   enddo
!
   do j = 1,4
     do i = 1,iout
       sum(i,j) = 0.e0
     enddo
   enddo
!
   do kpt = 1,4
     do i = 1,iout
       if (xx(i,kpt).gt.0.e0) then
         nn(i) = nn(i) + 1
         sum(i,1) = sum(i,1) + wgt(i,kpt)
         sum(i,2) = sum(i,2) + tt(i,kpt)
         sum(i,3) = sum(i,3) + bb(i,kpt)
       endif                                                               
     enddo
     do i = 1,iout
       sum(i,4) = sum(i,4) + wgt(i,kpt) * xx(i,kpt)
     enddo
   enddo
!
   do i = 1,iout
     if (nn(i).eq.1.and.sum(i,1).gt.0.49e0) go to 17
     if (nn(i).eq.2.and.sum(i,1).ge.0.45e0) go to 17
     if (nn(i).ge.3) go to 17
!
     ctop(i) = 0.e0
     cbot(i) = 100.e0
     camt(i) = 0.e0
!
     cycle
!
17   continue
!
     ltop = sum(i,2)/nn(i) + 0.5e0
     lbot = sum(i,3)/nn(i) + 0.5e0
     ctop(i) = ltop
     cbot(i) = lbot
     camt(i) = sum(i,4)
   enddo
!
   return                                                                    
!
!--- polar region-no extrapolation
!
600 continue
!
   ja = iabs(inslat)
   do i = 1,iout
     ileft=i
     irght=ileft+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     wgtlon=0.0
!
!----    get left point on nearest latitude
!
     xx(i,1) = cv(ileft,ja)
     xx(i,2) = cv(irght,ja)
     wgt(i,1) = 1.e0-wgtlon
     wgt(i,2) = wgtlon
     tt(i,1) = cvt(ileft,ja)
     tt(i,2) = cvt(irght,ja)
     bb(i,1) = cvb(ileft,ja)
     bb(i,2) = cvb(irght,ja)
   enddo
!
!---      nn will be number of surrounding pts with cld (gt zero)
!
   do i = 1,iout
     nn(i) = 0
   enddo
!
   do j = 1,4
     do i = 1,iout
       sum(i,j) = 0.e0
     enddo
   enddo
!
   do kpt = 1,2
     do i = 1,iout
       if (xx(i,kpt).gt.0.e0) then
         nn(i) = nn(i) + 1
         sum(i,1) = sum(i,1) + wgt(i,kpt)
         sum(i,2) = sum(i,2) + tt(i,kpt)
         sum(i,3) = sum(i,3) + bb(i,kpt)
       endif
     enddo
     do i = 1,iout
       sum(i,4) = sum(i,4) + wgt(i,kpt) * xx(i,kpt)
     enddo
   enddo
!
   do i = 1,iout
     if (nn(i).eq.1.and.sum(i,1).gt.0.7e0) go to 27
     if (nn(i).eq.2) go to 27
!
     ctop(i) = 0.e0
     cbot(i) = 100.e0
     camt(i) = 0.e0
     cycle
!
27   continue
!
     ltop = sum(i,2)/nn(i) + 0.5e0
     lbot = sum(i,3)/nn(i) + 0.5e0
     ctop(i) = ltop
     cbot(i) = lbot
     camt(i) = sum(i,4)
   enddo
!
   return
   end subroutine rad_cloud_interp_sub
!-------------------------------------------------------------------------------
#ifdef NIM
!
!-------------------------------------------------------------------------------
   subroutine rad_co2_read1(nfile,rco2,sgtmp,co21d,co22d,co21d3,co21d7)
!-------------------------------------------------------------------------------
!
! subroutine:    rad_co2_read1
!
! abstract:
!  -initializes arrays for -lw- radiation
!  -reads co2 transmission function data(from external file),
!   which has been pre-computed for current vertical coordinate on
!   the front-end machine. word conversion between front-end and c205
!   occurs here. also call tabl86 to set up tables for lw calculation
!   this code (rad_co2_read) is only called once...
!
! program history log:
!   84-01-01  fels and schwarzkopf,gfdl.
!   89-07-07  kenneth campana - removed unnecessary code and added
!                               reading and word conversion of co2 data.
!   89-11-29  kenneth campana - commented co2 reads because they
!                               are not yet ready for the new gfdl lw.
!
!   02-09-09  yifeng cui      - RMP mp    org:sdsc
!
! usage:    call rad_co2_read1(nfile,rco2,sgtmp,co21d,co22d,co21d3,co21d7)
!   input argument list:
!     nfile    - integer name of external co2 file.
!
!    *******************************************************************
!    *                           c o n r a d                           *
!    *    read co2 transmission data from unit(nfile)for new vertical  *
!    *      coordinate tests      ...                                  *!    *
!    these arrays used to be in block data    ...k.campana-mar 90 *
!    *******************************************************************
!
!                 co2 data tables for user''s vertical coordinate
!
!   the following common blocks contain pretabulated co2 transmission
!       functions, evaluated using the methods of fels and
!       schwarzkopf (1981) and schwarzkopf and fels (1985),
!-----  the 2-dimensional arrays are
!                    co2 transmission functions and their derivatives
!        from 109-level line-by-line calculations made using the 1982
!        mcclatchy tape (12511 lines),consolidated,interpolated
!        to the nmc mrf vertical coordinatte,and re-consolidated to a
!        200 cm-1 bandwidth. the interpolation method is described in
!        schwarzkopf and fels (j.g.r.,1985).
!-----  the 1-dim arrays are
!                  co2 transmission functions and their derivatives
!          for tau(i,i+1),i=1,l,                                               !
!          where the values are not obtained by quadrature,but are the
!            actual transmissivities,etc,between a pair of pressures.
!          these used only for nearby layer calculations including qh2o.
!-----  the weighting function gtemp=p(k)**0.2*(1.+p(k)/30000.)**0.8/
!         1013250.,where p(k)=pressure,nmc mrf(new)  l18 data levels for
!         pstar=1013250.
!-----  stemp is us standard atmospheres,1976,at data pressure levels
!        using nmc mrf sigmas,where pstar=1013.25 mb (ptz program)
!====>   begin here to get constants for radiation package
!
!-------------------------------------------------------------------------------
   use rdparm
   use co2dta
   use comio
!-------------------------------------------------------------------------------
   real                 ::  sgtmp(lp1,2),co21d(l,6),co22d(lp1,lp1,6)
   real                 ::  co21d3(lp1,6),co21d7(lp1,6)
!
   rewind nfile
!
!  read in pre-computed co2 transmission data....
!
   do kk = 1,2
     read(nfile) (sgtmp(i,kk),i=1,lp1)
   enddo
   do kk = 1,6
     read(nfile) (co21d(i,kk),i=1,l)
   enddo
   do kk = 1,6
     read(nfile) ((co22d(i,j,kk),i=1,lp1),j=1,lp1)
   enddo
   do kk = 1,6
     read(nfile) (co21d3(i,kk),i=1,lp1)
   enddo
   do kk = 1,6
     read(nfile) (co21d7(i,kk),i=1,lp1)
   enddo
#ifndef NIMAQUA
   do kk = 1,26
     read(nfile)
   enddo
#endif
!
!  read co2 concentration in ppm (defaulted in gradfs if missing)
!
   read(nfile,end=31) rco2
!
31 continue
!
#ifndef NOPRINT
   write(6,*) 'co2 concentration is ',rco2
#endif
   rewind nfile
!
   return
   end subroutine rad_co2_read1
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_co2_read2(sgtmp,co21d,co22d,co21d3,co21d7)
!-------------------------------------------------------------------------------
   use co2dta
   use comio,   only : iope
!-------------------------------------------------------------------------------
   real                 ::  sgtmp(lp1,2),co21d(l,6),co22d(lp1,lp1,6)
   real                 ::  co21d3(lp1,6),co21d7(lp1,6)
!
   do k = 1,lp1
     stemp(k) = sgtmp(k,1)
     gtemp(k) = sgtmp(k,2)
   enddo
   do k = 1,l
     cdtm51(k) = co21d(k,1)
     co2m51(k) = co21d(k,2)
     c2dm51(k) = co21d(k,3)
     cdtm58(k) = co21d(k,4)
     co2m58(k) = co21d(k,5)
     c2dm58(k) = co21d(k,6)
   enddo
   do j = 1,lp1
     do i = 1,lp1
       cdt51(i,j) = co22d(i,j,1)
       co251(i,j) = co22d(i,j,2)
       c2d51(i,j) = co22d(i,j,3)
       cdt58(i,j) = co22d(i,j,4)
       co258(i,j) = co22d(i,j,5)
       c2d58(i,j) = co22d(i,j,6)
     enddo
   enddo
   do k = 1,lp1
     cdt31(k) = co21d3(k,1)
     co231(k) = co21d3(k,2)
     c2d31(k) = co21d3(k,3)
     cdt38(k) = co21d3(k,4)
     co238(k) = co21d3(k,5)
     c2d38(k) = co21d3(k,6)
   enddo
   do k = 1,lp1
     cdt71(k) = co21d7(k,1)
     co271(k) = co21d7(k,2)
     c2d71(k) = co21d7(k,3)
     cdt78(k) = co21d7(k,4)
     co278(k) = co21d7(k,5)
     c2d78(k) = co21d7(k,6)
   enddo
!
!......    define tables for lw radiation
!
   call rad_lw_table
!
   return
   end subroutine rad_co2_read2
!-------------------------------------------------------------------------------
#endif
