#include <define.h>
#include "abort.h"
   subroutine file_write_flux(nn)
!-------------------------------------------------------------------------------
!
#ifdef NOAHYDRO
   use constant, only : hvap=>hvap_
#endif
   use varsfc, only   : lsoil_,nsoil_,msub_
   use comsfc
#ifdef DFS
   use dfsvar, only   : iope
#endif
   use paramodel, only : ILEN2S=>LONF2S,JLEN2S=>LATG2S,                        &
                         ILEN=>LONFD,JLEN=>LATGD,                              &
                         ILEN2=>LONF2F,JLEN2=>LATG2F,lonf_,latg_
#ifdef RIVER
   use paramodel, only : io2_,jo2_
   use comfrivh
#endif
#ifndef RMP
   use comio
#endif
#ifdef MP
#ifndef DFS
   use commpi
#else
   use commpi, only : mype, master
#endif
#endif
#ifdef RMP
   use rscomf_rerun
   use rscommap
#else		/* not RMP */
   use comfphys
   use comfspec_vr
   use comfver
   use comfgrid
   use radiag
#endif			/* RMP */
   use module_trans
!-------------------------------------------------------------------------------
#ifdef RIVER
   real                   ::  gdriv2(io2_,jo2_)
   integer                ::  lenr
   logical, allocatable   ::  lbmr(:)
   character, allocatable ::  gr(:)
#endif
!
   integer,parameter  ::  iprs=1,itemp=11,iznlw=33,imerw=34,isphum=51,ipwat=54,&
                          ipcpr=59,isnowd=65,icldf=71,iccldf=72,               &
                          islmsk=81,iz0cm=83,ialbdo=84,isoilm=144,icemsk=91,   &
                          ilhflx=121,ishflx=122,izws=124,imws=125,ighflx=155,  &
                          iuswfc=160,idswfc=161,iulwfc=162,idlwfc=163,         &
                          inswfc=164,inlwfc=165,                               &
                          idswvb=166,idswvd=167,idswnb=168,idswnd=169,         &
                          itmx=15,itmn=16,irnof=90,iep=145,                    &
                          isnwfll=64,isnwevp=230,isnwmlt=229,                  &
                          icldwk=146,izgw=147,imgw=148,ihpbl=221,              &
                          idswf=204,idlwf=205,iuswf=211,iulwf=212,icpcpr=214,  &
                          icanopy=223,iqull=202,iqvll=203,                     &
                          ivegtyp=225,ivegcov=87,isoltyp=226,                  &
                          iz0cmt=118,iomld=190,                                &
#ifdef RIVER
                          igdriv=21,irflow=22,iimap=23,iroff=24,               &
#endif
#ifdef GWDC
                          izgwc=170,imgwc=171,                                 &
#endif
#ifdef VICLSM1
                          igdh=188,                                            &
#endif
#ifndef NOAHYDRO
                          ialhtfl=236,ievcnp=180,ibgrun=234
#else
                          ialhtfl=236,ievcnp=180,ibgrun=234,iraintot2=255
#endif
   integer,parameter  ::  isfc=1,itoa=8,ielev=105,                             &
                          isglev=107,idbls=111,i2dbls=112,icolmn=200,          &
                          ilcbl=212,ilctl=213,ilclyr=214,                      &
                          imcbl=222,imctl=223,imclyr=224,                      &
                          ihcbl=232,ihctl=233,ihclyr=234
   integer,parameter  ::  inst=10,itr=3,iacc=4
   integer,parameter  ::  ifhour=1,ifday=2
!
#ifdef MP
#ifdef RMP
#define MPGP2F rmpgp2f
#else
#define MPGP2F mpgp2f
#endif
#endif
!
! full arrays
!
   integer                ::  len
   logical,allocatable    ::  lbm(:)
   character,allocatable  ::  g(:)
   !
   integer,parameter    ::  nfld=16
   integer              ::  ipur(nfld),itlr(nfld)
   data ipur/iulwf , iuswf , iuswf , idswf ,  icldf,   iprs,                   &
             iprs, itemp ,  icldf,   iprs,   iprs, itemp ,                     &
             icldf,   iprs,   iprs, itemp /
   data itlr/itoa  , itoa  , isfc  , isfc  , ihclyr, ihctl ,                   &
             ihcbl , ihctl , imclyr, imctl , imcbl , imctl ,                   &
             ilclyr, ilctl , ilcbl , ilctl /
   real rtimer(nfld)
!
   integer              ::  lens
!
!  full working arrays
!
   real, allocatable  ::  work(:),slmsep(:),fluxf(:,:)
   real               ::  fluxw(ILEN2,JLEN2,26)
!
   integer            ::  lsoimlev(lsoil_+1)
#ifdef OSULSM1
   data lsoimlev/0,10,200/
#endif
#ifdef OSULSM2
   data lsoimlev/0,10,200/
#endif
#ifdef NOALSM1
   data lsoimlev/0,10,40,100,200/
#endif
#ifdef VIC
   data lsoimlev/0,1,2,3/
#endif
#ifdef VIC
   integer            ::  lsoitlev(nsoil_+1)
   data lsoitlev/0,1,2,3,4,5/
#endif
!
#if defined(RMP) && defined(MP)
   real, allocatable  ::  fullat(:),fullon(:)
#endif
!
   integer            ::  ids(255),                                            &
                          iens(5)
#ifndef RMP
#undef MRG_POST
#endif
!
   character(len=128)  ::  fno
#ifdef ASSIGN
   character(len=128)  ::   asgnstr
#endif
!
   lenr=io2_*(jo2_)
   len=ILEN*JLEN
   lens=ILEN2S*JLEN2S
#ifdef RIVER
   allocate(lbmr(lenr))
   allocate(gr(200+lenr*(16+1)/8))
#endif
   allocate(lbm(len))
   allocate(g(300+2*len*(16+1)/8))
   allocate(work(len),slmsep(len),fluxf(len,4))
!
   if(nn.le.0) return
!
   if(iope) then
#ifdef RMP
#ifndef MRG_POST
     call file_name('r_flx',5,thour,fno,ncho)
#endif
#else			/* nor RMP */
     call file_name('flx',3,thour,fno,ncho)
#endif			/* RMP */
#ifndef MRG_POST
#ifdef ASSIGN
     write(asgnstr,'(23hassign -s unblocked  u:,I2,)') nn
     call assign('assign -R')
     call assign(asgnstr)
#endif
     open(unit=nn,file=fno(1:ncho),form='unformatted',err=900)
     go to 901
 900 continue
     write(6,*) ' error in opening file ',fno(1:ncho)
     call MPABORT
 901 continue
#ifndef NOPRINT
     write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
#endif
     rewind nn
#endif /* ~MRG_POST end */
   endif
!
   call idsdef(1,ids)
   ilpds=28
   if(icen2.eq.2) ilpds=45
!
   ienst=0
   iensi=0
   iens(1)=1
   iens(2)=ienst
   iens(3)=iensi
   iens(4)=1
   iens(5)=255
   iyr=idate(4)
   imo=idate(2)
   ida=idate(3)
   ihr=idate(1)
   ifhr=nint(fhour)
   ithr=nint(thour)
!
! -- Never show 0-0hr averages, always round up to 1 hour
!
   if(ithr.eq.0) then
     iavg=1
   else
     iavg=itr
   endif
!
   dhour=dtpost/3600.
   if(dtpost.gt.0) then
     rtime=1./dtpost
   else
     rtime=0.
   endif
!
   secswr=max(dhour,dtswav) * 3600.
   seclwr=max(dhour,dtlwav) * 3600.
   if(secswr.gt.0.) then
     rtimsw=1./secswr
   else
     rtimsw=1.
   endif
   if(seclwr.gt.0.) then
     rtimlw=1./seclwr
   else
     rtimlw=1.
   endif
   do n = 1,nfld
     rtimer(n)=rtimsw
   enddo
   rtimer(1)=rtimlw
!
#if defined(DFS)
   idrt=0
   rlat1=90.-(180./(2.*latg_))
   rlon1=0.
   rlat2=-rlat1
   rlon2=0.
   delx=360./lonf_
   dely=180./latg_
   ortru=0.
   proj=0.
#elif defined(RMP)
   proj=rproj
   delx=rdelx
   dely=rdely
#ifdef MP
   allocate (fullat(ILEN*JLEN),fullon(ILEN*JLEN))
   call rmpgp2f(flat,ILEN2S,JLEN2S,fullat,ILEN2,JLEN2,1)
   call rmpgp2f(flon,ILEN2S,JLEN2S,fullon,ILEN2,JLEN2,1)
#define FLAT fullat
#define FLON fullon
#else
#define FLAT flat
#define FLON flon
#endif
   call rmp_trans2output_grid(FLAT,1)
   call rmp_trans2output_grid(FLON,1)
   rlat1=FLAT(1)
   rlat2=FLAT(ILEN*JLEN)
   rlon1=FLON(1*1)
   rlon2=FLON(ILEN*JLEN)
#ifndef MP
   call rmp_trans2model_grid(FLAT,1)
   call rmp_trans2model_grid(FLON,1)
#else
   deallocate (fullat,fullon)
#endif
   nproj=rproj
   if( nproj.eq.0 ) then   ! mercater
     idrt=1
     ortru=rtruth
     proj=0.0
   elseif( abs(nproj).eq.1 ) then  ! polar projection
     idrt=5
     ortru=rorient
     proj=rproj
   elseif( abs(nproj).eq.2 ) then  ! lambert conformal
     idrt=3
     ortru=rorient
     proj=rproj
   else
     print *,' error in rwrtsfc for projection type....'
   endif
#else /* ~RMP */
#ifdef SMP
   idrt=0
#else
   idrt=4
#endif
   rlat1=0.
   rlon1=0.
   rlat2=0.
   rlon2=0.
   delx=0.
   dely=0.
   ortru=0.
   proj=0.
#endif /* RMP end */
!
#ifdef MP
   call MPGP2F(slmsk,ILEN2S,JLEN2S,slmsep,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       slmsep(n)=slmsk(i,j)
     enddo
   enddo
!
#endif
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(slmsep,1)
#else
     call dyn_trans2output_grid(slmsep,1)
#ifdef SMP
     cl1=90.
#else
     cl1=colrad(1)
#endif
#endif /* RMP end */
   endif
!
#ifndef RMP
#undef MRG_POST
#endif
#ifndef MRG_POST
#define GRIBIT  write(nn)
#else
#define GRIBIT  call file_make_grib(
#endif
! 1  surface u-stress
!
#ifdef MP
   call MPGP2F(dusfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dusfc(i,j)*rtime
#endif
       enddo
     enddo
     maxbit=16
#ifndef NOAHYDRO
!    maxbit2=16
     maxbit2=32
#else
     maxbit2=32
#endif
#if defined(DFS)
     colat=rlat1
#elif defined (RMP)
     colat=0.
#else
     colat=colrad(1)
#endif
     iptv=2
     ibm0=0
     ibm1=1
     il1k=0
     il2k=0
     ip1=0
     ip2=0
     ina=0
     inm=0
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'uflx')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,izws   ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(izws   ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'uflx grib created    '
#endif
     else
       print *,'uflx grib make failed'
     endif
#endif
   endif
!
! 2  surface v-stress
!
#ifdef MP
   call MPGP2F(dvsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dvsfc(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'vflx')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,imws   ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(imws   ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'vflx grib created    '
#endif
     else
       print *,'vflx grib make failed'
     endif
#endif
   endif
!
! 3  sensible heat flux
!
#ifdef MP
   call MPGP2F(dtsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dtsfc(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'tflx')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ishflx ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ishflx ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'tflx grib created    '
#endif
     else
       print *,'tflx grib make failed'
     endif
#endif
   endif
!
! 4  latent heat flux
!
#ifdef MP
   call MPGP2F(dqsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dqsfc(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'qflx')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ilhflx ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ilhflx ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'qflx grib created    '
#endif
     else
       print *,'qflx grib make failed'
     endif
#endif
   endif
!
! 5  surface temperature
!
#ifdef MP
   call MPGP2F(tsea,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       work(n)=tsea(i,j)
     enddo
   enddo
#endif
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'tsfc')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,itemp  ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(itemp  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'tsfc grib created    '
#endif
     else
       print *,'tsfc grib make failed'
     endif
#endif
   endif
!
! 6  soil moisture (multi-level)
!
   do k = 1,lsoil_
#ifdef MP
     call MPGP2F(smc(1,1,k),ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
         work(n)=smc(i,j,k)
       enddo
     enddo
#endif
     if(iope) then
#ifdef RMP
       call rmp_trans2output_grid(work,1)
#else
       call dyn_trans2output_grid(work,1)
#endif
       il1k=lsoimlev(k)
       il2k=lsoimlev(k+1)
#ifdef DBG
       write(6,*)'soil level=',k
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'soilw')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isoilm,i2dbls,il1k,il2k,iyr,imo,ida,ihr,                   &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(isoilm ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'soilm level',k,' grib created'
#endif
       else
         print *,'soilm level',k,' grib make failed'
       endif
#endif
     endif
   enddo
!
   il1k=0
   il2k=0
!
! 8  soil temperature (multi-level)
!
#ifndef VIC
   do k = 1,lsoil_
#else
   do k = 1,nsoil_
#endif
#ifdef MP
     call MPGP2F(stc(1,1,k),ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
         work(n)=stc(i,j,k)
       enddo
     enddo
#endif
     if(iope) then
#ifdef RMP
       call rmp_trans2output_grid(work,1)
#else
       call dyn_trans2output_grid(work,1)
#endif
#ifdef VIC
       il1k=lsoitlev(k)
       il2k=lsoitlev(k+1)
#else
       il1k=lsoimlev(k)
       il2k=lsoimlev(k+1)
#endif
#ifdef DBG
       print *,'soil level=',k
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'soilt')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,itemp,i2dbls,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(itemp  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'soil temperature level 1 grib created    '
#endif
       else
         print *,'soil temperature level 1 grib make failed'
       endif
#endif
     endif
   enddo
!
   il1k=0
   il2k=0
!
! 10  soil water equivalent snow depth
!
#ifdef MP
   call MPGP2F(snoweq,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       work(n)=snoweq(i,j)
     enddo
   enddo
#endif
!
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'snow')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isnowd ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(isnowd ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
!
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'snow depth grib created    '
#endif
     else
       print *,'snow depth grib make failed'
     endif
#endif
   endif
!
! 11 downward longwave radiation flux at the surface
!
#ifdef MP
   call MPGP2F(dlwsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dlwsfc(i,j)*rtime
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
!
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'dlwsfc')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,idlwf  ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(idlwf),iens,            &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'dlwflx sfc grib created    '
#endif
     else
       print *,'dlwflx sfc grib make failed'
     endif
#endif
   endif
!
! 12 upword longwave radiation flux at the surface
!
#ifdef MP
   call MPGP2F(ulwsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=ulwsfc(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
!
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'ulwsfc')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iulwf  ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iulwf  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'ulwflx sfc grib created    '
#endif
     else
       print *,'ulwflx sfc grib make failed'
     endif
#endif
   endif
!
#ifdef MP
   call MPGP2F(fluxr,ILEN2S,JLEN2S,fluxw,ILEN2,JLEN2,26)
#else
   do k = 1,26
     do j = 1,JLEN2
       do i = 1,ILEN2
         fluxw(i,j,k)=fluxr(i,j,k)
       enddo
     enddo
   enddo
#endif
!
   if(iope) then
!
!   13 upward long wave radiation flux at the top of the atmosphere
!   14 upward short wave radiation flux at the top of the atmosphere
!   15 upward short wave radiation flux at the surface
!   16 downward short wave radiation flux at the surface
!
     do k = 1,4
       ii=0
       do j = 1,JLEN2
         do i = 1,ILEN2
           ii=ii+1
           work(ii)=fluxw(i,j,k)
         enddo
       enddo
       do n = 1,len
         work(n)=work(n)*rtimer(k)
       enddo
#ifdef SMP_NUDGING_RAD 
!
! SWUP
!
       if (k.eq.3) then
#ifdef MP
         call MPGP2F(scmsup,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
         n=0
         do j = 1,JLEN2
           do i = 1,ILEN2
             n=n+1
#ifdef MP
             work(n)=work(n)*rtime
#else
             work(n)=scmsup(i,j)*rtime
#endif
           enddo
         enddo
       endif
!
! SWDN
!
       if (k.eq.4) then
#ifdef MP
         call MPGP2F(scmsdn,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
         n=0
         do j = 1,JLEN2
           do i = 1,ILEN2
             n=n+1
#ifdef MP
             work(n)=work(n)*rtime
#else
             work(n)=scmsdn(i,j)*rtime
#endif
           enddo
         enddo
       endif
#endif /* SMP_NUDGING_RAD end */
#ifdef RMP
       call rmp_trans2output_grid(work,1)
#else
       call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'radflxes')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipur(k),itlr(k),il1k,il2k,iyr,imo,ida,ihr,                 &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(k)),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'radiation fluxes grib created    '
#endif
       else
         print *,'radiation fluxes grib make failed'
       endif
#endif
     enddo  ! k=1,4
!
     iw1=1
     iw2=iw1+4*ILEN2
     iw3=iw2+4*ILEN2
     iw4=iw3+4*ILEN2
     iw5=iw4+4*ILEN2
     do k = 5,7
       ij=0
       do j = 1,JLEN2
         do i = 1,ILEN2
           ij=ij+1
           if(fluxw(i,j,k).gt.0.) then
             fluxf(ij,2) = fluxw(i,j,k+3) / fluxw(i,j,k)
             fluxf(ij,3) = fluxw(i,j,k+6) / fluxw(i,j,k)
             fluxf(ij,4) = fluxw(i,j,k+9) / fluxw(i,j,k)
             fluxf(ij,1) = fluxw(i,j,k) * rtimsw
           else
!
!  zero cld top temp if no clds--safety, cause use zero in ggintt
!
             fluxw(i,j,k+9) = 0.
             fluxf(ij,2) = fluxw(i,j,k+3)
             fluxf(ij,3) = fluxw(i,j,k+6)
             fluxf(ij,4) = fluxw(i,j,k+9)
             fluxf(ij,1) = fluxw(i,j,k) * rtimsw
           endif
         enddo
       enddo
!
!   17 high level cloud amount
!   21 middle level cloud amount
!   25 low level cloud amount
!
       k4=4+(k-5)*4
#ifdef RMP
       call rmp_trans2output_grid(fluxf(1,1),1)
#else
       call dyn_trans2output_grid(fluxf(1,1),1)
#endif
       do n = 1,len
         work(n)=fluxf(n,1)*1.e2
       enddo
       l=k4+1
#ifdef DBG
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'cloud amount')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,                 &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'cloud amount grib created    '
#endif
       else
         print *,'cloud amount grib make failed'
       endif
#endif
!
!   18 high level cloud top pressure
!   22 middle level cloud top pressure
!   26 low level cloud top pressure
!
#ifdef RMP
       call rmp_trans2output_grid(fluxf(1,2),1)
#else
       call dyn_trans2output_grid(fluxf(1,2),1)
#endif
       do n = 1,len
         work(n)=fluxf(n,2)*1.e3
       enddo
       l=k4+2
#ifdef DBG
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'cloud top pres')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,                 &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'cloud top pressure grib created    '
#endif
       else
         print *,'cloud top pressure grib make failed'
       endif
#endif
!
!   19 high level cloud bottom pressure
!   23 middle level cloud bottom pressure
!   27 low level cloud bottom pressure
!
#ifdef RMP
       call rmp_trans2output_grid(fluxf(1,3),1)
#else
       call dyn_trans2output_grid(fluxf(1,3),1)
#endif
       do n = 1,len
         work(n)=fluxf(n,3)*1.e3
       enddo
       l=k4+3
#ifdef DBG
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'cloud bot pres')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,                 &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'cloud bottom pressure grib created    '
#endif
       else
         print *,'cloud bottom pressure grib make failed'
       endif
#endif
!
!   20 high level cloud top tempeature
!   24 middle level cloud top tempeature
!   28 low level cloud top tempeature
!
#ifdef RMP
       call rmp_trans2output_grid(fluxf(1,4),1)
#else
       call dyn_trans2output_grid(fluxf(1,4),1)
#endif
       do n = 1,len
         work(n)=fluxf(n,4)
       enddo
       l=k4+4
#ifdef DBG
       call print_maxmin_six(work,ILEN*JLEN,1,1,1,'cloud top tmp')
#endif
       GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                            &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,                 &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
       if(ierr.eq.0) then
         call file_write_byte(nn,lg,g)
#ifdef DBG
         print *,'cloud top temp grib created    '
#endif
       else
         print *,'cloud top temp grib make failed'
       endif
#endif
     enddo  ! k=5,7
   endif     ! iope
!
! 29  precipitation rate
!
#ifdef MP
   call MPGP2F(raintot,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
       n=n+1
#ifdef MP
       work(n)=work(n)*1.e3*rtime
#else
       work(n)=raintot(i,j)*1.e3*rtime
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
!
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'prate')
#endif
#ifndef ES
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,                             &
#else
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
#endif
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipcpr,isfc,il1k,il2k,iyr,imo,ida,ihr,                      &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipcpr),iens,            &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'precip grib created    '
#endif
     else
       print *,'precip grib make failed'
     endif
#endif
   endif
!
! 30  convective precpipitation rate
!
#ifdef MP
   call MPGP2F(raincps,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*1.e3*rtime
#else
         work(n)=raincps(i,j)*1.e3*rtime
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
!
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'cprate')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,icpcpr ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icpcpr ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'convective precip grib created    '
#endif
     else
       print *,'convective precip grib make failed'
     endif
#endif
   endif
!
! 31  ground flux
!
#ifdef MP
   call MPGP2F(gflux,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=gflux(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
!
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'gflux')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ighflx ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ighflx ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'ground heat flux grib created    '
#endif
     else
       print *,'ground heat flux grib make failed'
     endif
#endif
   endif
!
! 32  land sea mask
!
   if(iope) then
     do n = 1,len
       work(n)=mod(slmsep(n),2.)
     enddo
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'lsmask')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,islmsk ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(islmsk ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'land sea mask grib created    '
#endif
     else
       print *,'land sea mask grib make failed'
     endif
#endif
   endif
!
! 33  ice mask
!
   if(iope) then
     do n = 1,len
       work(n)=max(slmsep(n)-1.,0.)
     enddo
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'ice mask')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,icemsk ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(icemsk ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'ice mask grib created    '
#endif
     else
       print *,'ice mask grib make failed'
     endif
#endif
   endif
!
! 34  10m wind u
!
#ifdef MP
   call MPGP2F(u10m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=u10m(i,j)
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
!
     il2k=10
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'u10m')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iznlw,ielev,il1k,il2k,iyr,imo,ida,ihr,                     &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(iznlw  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'10m u grib created    '
#endif
     else
       print *,'10m u grib make failed'
     endif
#endif
     il2k=0
   endif
!
! 35  10m wind v
!
#ifdef MP
   call MPGP2F(v10m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=v10m(i,j)
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
     il2k=10
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'v10m')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,imerw,ielev ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(imerw  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'10m v grib created    '
#endif
     else
       print *,'10m v grib make failed'
     endif
#endif
     il2k=0
   endif
!
! 36  temperature 2m
!
#ifdef MP
   call MPGP2F(t2m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=t2m(i,j)
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
     il2k=2
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'t2m')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,itemp,ielev ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(itemp ),iens,           &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'2m temp grib created    '
#endif
     else
       print *,'2m temp grib make failed'
     endif
#endif
     il2k=0
   endif
!
! 37  specific humidity 2m
!
#ifdef MP
   call MPGP2F(q2m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
       n=n+1
#ifdef MP
       work(n)=work(n)
#else
       work(n)=q2m(i,j)
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
     il2k=2
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'q2m')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isphum,ielev,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(isphum ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'2m specific humidity grib created    '
#endif
     else
       print *,'2m specific humidity grib make failed'
     endif
#endif
     il2k=0
   endif
!
! 39  max temperature
!
#ifdef MP
   call MPGP2F(tmpmax,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       work(n)=tmpmax(i,j)
     enddo
   enddo
#endif
!
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
     il2k=2
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'max temp')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,itmx,ielev  ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(itmx   ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'max temp grib created    '
#endif
     else
       print *,'max temp grib make failed'
     endif
#endif
     il2k=0
   endif
!
!  reset tmpmax
!
   do j = 1,JLEN2S
     do i = 1,ILEN2S
       tmpmax(i,j) = 0.
     enddo
   enddo
!
! 40  min temperature
!
#ifdef MP
   call MPGP2F(tmpmin,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2S
     do i = 1,ILEN2S
       n=n+1
       work(n)=tmpmin(i,j)
     enddo
   enddo
#endif
!
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
     il2k=2
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'min temp')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,itmn,ielev  ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(itmn   ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'min temp grib created    '
#endif
     else
       print *,'min temp grib make failed'
     endif
#endif
     il2k=0
   endif
!
!  reset tmpmin
!
   do j = 1,JLEN2S
     do i = 1,ILEN2S
       tmpmin(i,j) = 1.e10
     enddo
   enddo
!
! 41  runoff
!
#ifdef MP
   call MPGP2F(runoff,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*1.e3*rtime
#else
         work(n)=runoff(i,j) * 1.e3 * rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'runoff')
#endif
#ifndef ES
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,                             &
#else
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
#endif
               ilpds,iptv,icen,igen,                                           &
               ibm0,irnof,isfc  ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irnof  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'runoff grib created    '
#endif
     else
       print *,'runoff grib make failed'
     endif
#endif
   endif
!
! 42  potential evaporation
!
#ifdef MP
   call MPGP2F(ep,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=ep(i,j) * rtime
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'potevap')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iep,isfc    ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iep    ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'potential evap grib created    '
#endif
     else
       print *,'potential evap grib make failed'
     endif
#endif
   endif
#if defined SAS || defined KSAS
!
! 43  cloud work function
!
#ifdef MP
   call MPGP2F(cldwrk,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=cldwrk(i,j) * rtime
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'cloud wk func.')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,icldwk,isfc ,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icldwk ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'cloud work func grib created    '
#endif
     else
       print *,'cloud work func grib make failed'
     endif
#endif
   endif
#endif         /* FOR SAS */
!
! 44  gravity wave drag (u)
!
#ifdef MP
   call MPGP2F(dugwd,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dugwd(i,j)*rtime
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'ugwd')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,izgw   ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(izgw   ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'gravity wave drag u grib created    '
#endif
     else
       print *,'gravity wave drag u grib make failed'
     endif
#endif
   endif
!
! 45  gravity wave drag (v)
!
#ifdef MP
   call MPGP2F(dvgwd,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dvgwd(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'vgwd')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,imgw   ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(imgw   ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'gravity wave drag v grib created    '
#endif
     else
       print *,'gravity wave drag v grib make failed'
     endif
#endif
   endif
!
! 46  pbl height
!
#ifdef MP
   call MPGP2F(hpbl,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       work(n)=hpbl(i,j)
     enddo
   enddo
!
#endif
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'pblh')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ihpbl,isfc,il1k,il2k,iyr,imo,ida,ihr,                      &
               ifhour,ithr,0,inst,ina,inm,icen2,ids(ihpbl),iens,               &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'pbl height grib created    '
#endif
     else
        print *,'pbl height grib make failed'
     endif
#endif
   endif
!
! 47  precipitable water
!
#ifdef MP
   call MPGP2F(pwat,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       work(n)=pwat(i,j)
     enddo
   enddo
#endif
!
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'pwater')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ipwat,icolmn,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0,inst,ina,inm,icen2,ids(ipwat),iens,               &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'precipitable water grib created    '
#endif
     else
       print *,'precipitable water grib make failed'
     endif
#endif
   endif
!
! 48  snow fall amount
!
#ifdef MP
   call MPGP2F(snowfall,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*1.e3*rtime
#else
         work(n)=snowfall(i,j)*1.e3*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'snowfall')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isnwfll,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isnwfll),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'snow fall grib created    '
#endif
     else
       print *,'snow fall grib make failed'
     endif
#endif
   endif
!
! 49  snow evaporation
!
#ifdef MP
   call MPGP2F(snowevap,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=snowevap(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'snowevap')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isnwevp,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isnwevp),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'snow evap grib created    '
#endif
     else
       print *,'snow evap grib make failed'
     endif
#endif
   endif
!
! 50  snow melt
!
#ifdef MP
   call MPGP2F(snowmelt,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=snowmelt(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'snow melt')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isnwmlt,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isnwmlt),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'snow melt grib created    '
#endif
     else
       print *,'snow melt grib make failed'
     endif
#endif
   endif
!
! 51  moisture flux (u)
!
#ifdef MP
   call MPGP2F(qull,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=qull(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'qfulx-u')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iqull,icolmn,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iqull),iens,            &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'moisture flux u grib created    '
#endif
     else
       print *,'moisture flux u grib make failed'
     endif
#endif
   endif
!
! 52  moisture flux (v)
!
#ifdef MP
   call MPGP2F(qvll,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=qvll(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'qflux-v')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iqvll,icolmn,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iqvll),iens,            &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'moisture flux v grib created    '
#endif
     else
       print *,'moisture flux v grib make failed'
     endif
#endif
   endif
!
! 53  canopy water content
!
#ifdef MP
   call MPGP2F(canopy,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=canopy(i,j)
#endif
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'canopy water')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,icanopy,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0,inst,ina,inm,icen2,ids(icanopy),iens,             &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'canopy water content grib created    '
#endif
     else
       print *,'canopy water content grib make failed'
     endif
#endif
   endif
!
! 54 downward short wave flux at the top of the atmosphere
!
   if(iope) then
     ii=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         ii=ii+1
         work(ii)=fluxw(i,j,18)
       enddo
     enddo
     do n = 1,len
       work(n)=work(n)*rtime
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'dswtoa')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,idswf,itoa,il1k,il2k,iyr,imo,ida,ihr,                      &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(idswf),iens,            &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'dswflx top grib created    '
#endif
     else
       print *,'dswflx top grib make failed'
     endif
#endif
   endif
!
! 55  total cloud cover
!
   if(iope) then
     ij=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         ij=ij+1
         fluxf(ij,1) = fluxw(i,j,26) * rtimsw
       enddo
     enddo
!
#ifdef RMP
     call rmp_trans2output_grid(fluxf(1,1),1)
#else
     call dyn_trans2output_grid(fluxf(1,1),1)
#endif
     do n = 1,len
       work(n)=fluxf(n,1)*1.e2
     enddo
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'total cloud')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,icldf,icolmn,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icldf),iens,            &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'total cloud cover grib created    '
#endif
     else
       print *,'total cloud cover grib make failed'
     endif
#endif
   endif
!
! 56  albedo
!
   if(iope) then
     ii=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         ii=ii+1
         work(ii)=fluxw(i,j,17)
       enddo
     enddo
     do n = 1,len
       work(n)=work(n)*rtimsw * 100.
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'albedo')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ialbdo,isfc,il1k,il2k,iyr,imo,ida,ihr,                     &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ialbdo),iens,           &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'albedo grib created    '
#endif
     else
       print *,'albedo grib make failed'
     endif
#endif
   endif
!
! 57  surface roughness
!
#ifdef MP
   call MPGP2F(z0cm,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=z0cm(i,j)
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'z0')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iz0cm,isfc,il1k,il2k,iyr,imo,ida,ihr,                      &
               ifhour,ithr,0,inst,ina,inm,icen2,ids(iz0cm),iens,               &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'albedo grib created    '
#endif
     else
       print *,'surface roughness grib make failed'
     endif
#endif
   endif
!
! 58  vegetation type
!
#ifndef OSULSM1
#ifdef MP
   call MPGP2F(vtype,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=vtype(i,j)
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'vegtype')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ivegtyp,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ivegtyp),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'vegitation type grib created    '
#endif
     else
       print *,'vegitation type grib make failed'
     endif
#endif
   endif
#endif /* ~OSULSM1 end */
!
! 59  plantr/vegetation cover
!
#ifdef OSULSM1
#define VARIAB plantr
#else
#define VARIAB vfrac
#endif
!
#ifdef MP
   call MPGP2F(VARIAB,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=VARIAB(i,j)
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'vegcover')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ivegcov,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ivegcov),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'vegetation cover grib created    '
#endif
     else
       print *,'vegetation cover grib make failed'
     endif
#endif
   endif
#ifndef VIC
#ifndef OSULSM1
!
! 60  soil type
!
#ifdef MP
   call MPGP2F(stype,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=stype(i,j)
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'soiltype')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,isoltyp,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isoltyp),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'soil type grib created    '
#endif
     else
       print *,'soil type grib make failed'
     endif
#endif
   endif
#endif /* ~OSULSM1 end */
#endif /* ~VIC end */
!
! 61  Another latent heat flux (from phys_lsm_osu2)
!
#ifdef MP
   call MPGP2F(alhtfl,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=alhtfl(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'qflx phys_lsm_osu1')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ialhtfl,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ialhtfl),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'another latent heat flux grib created    '
#endif
     else
       print *,'another latent heat flux grib make failed'
     endif
#endif
   endif
!
! 62  Canopy evaporation
!
#ifdef HYDRO
#define EVCNP evcnp
#else
#ifndef NOAHYDRO
#define EVCNP ep
#else
#define EVCNP evcnp
#endif
#endif
#ifdef MP
   call MPGP2F(EVCNP,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
#ifndef NOAHYDRO
         work(n)=EVCNP(i,j)*rtime
#else
         work(n)=EVCNP(i,j)*rtime*1.e3*hvap
#endif
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'canopy evap')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ievcnp,isfc,il1k,il2k,iyr,imo,ida,ihr,                     &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ievcnp),iens,           &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'canopy evap grib created    '
#endif
     else
       print *,'canopy evap grib make failed'
     endif
#endif
   endif
!
! 62  Baseflow Groundwater Runoff
!
#ifdef MP
   call MPGP2F(bgrun,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*1.e03*rtime
#else
         work(n)=bgrun(i,j)*1.e03*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'base flow')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,ibgrun,isfc,il1k,il2k,iyr,imo,ida,ihr,                     &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ibgrun),iens,           &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'bgrun grib created    '
#endif
     else
       print *,'bgrun grib make failed'
     endif
#endif
   endif
!
! 63  Surface Pressure
!
#ifdef MP
   call MPGP2F(psurf,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*1.e3
#else
         work(n)=psurf(i,j)*1.e3
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'psurf')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iprs,isfc,il1k,il2k,iyr,imo,ida,ihr,                       &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(iprs),iens,             &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'surface pressure grib created    '
#endif
     else
       print *,'surface pressure grib make failed'
     endif
#endif
   endif
#ifdef NOAHYDRO
!
! 64  Second Precip
!
#ifdef MP
   call MPGP2F(raintot2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*1.e03*rtime
#else
         work(n)=raintot2(i,j)*1.e03*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'base flow')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iraintot2,isfc,il1k,il2k,iyr,imo,ida,ihr,                  &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iraintot2),iens,        &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'raintot2 grib created    '
#endif
     else
       print *,'raintot2 grib make failed'
     endif
#endif
   endif
#endif
!
! 65  surface thermal roughness
!
#ifdef MP
   call MPGP2F(z0cmt,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)
#else
         work(n)=z0cmt(i,j)
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'z0t')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iz0cmt ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(iz0cmt ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'albedo grib created    '
#endif
     else
       print *,'surface roughness(thermal) grib make failed'
     endif
#endif
   endif
!
! 66  convection induced gravity wave drag (u)
!
#ifdef GWDC
#ifdef MP
   call MPGP2F(dugwdc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dugwdc(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'ugwdc')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,izgwc  ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(izgwc  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'gravity wave drag u grib created    '
#endif
     else
       print *,'gravity wave drag u grib make failed'
     endif
#endif
   endif
!
! 67  convection induced gravity wave drag (v)
!
#ifdef MP
   call MPGP2F(dvgwdc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(iope) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=dvgwdc(i,j)*rtime
#endif
       enddo
     enddo
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'vgwdc')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,imgwc  ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(imgwc  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'gravity wave drag v grib created    '
#endif
     else
       print *,'gravity wave drag v grib make failed'
     endif
#endif
   endif
#endif
#ifdef VICLSM1
!
! 68  ground heat storage
!
#ifdef MP
   call MPGP2F(gflux,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#endif
   if(IOPE) then
     n=0
     do j = 1,JLEN2
       do i = 1,ILEN2
         n=n+1
#ifdef MP
         work(n)=work(n)*rtime
#else
         work(n)=gheat(i,j)*rtime
#endif
       enddo
     enddo
#ifndef DFS
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'gheat')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,igdh ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(igdh ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'ground heat storage grib created    '
#endif
     else
       print *,'ground heat storage grib make failed'
     endif
#endif
   endif
#endif 
!
! 69  ocean mixed layer depth (m)
!
#ifdef MP
   call MPGP2F(hml,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
#else
   n=0
   do j = 1,JLEN2
     do i = 1,ILEN2
       n=n+1
       work(n)=hml(i,j)
     enddo
   enddo
#endif
!
   if(iope) then
#ifdef RMP
     call rmp_trans2output_grid(work,1)
#else
     call dyn_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,ILEN*JLEN,1,1,1,'omld')
#endif
     GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,                              &
               ilpds,iptv,icen,igen,                                           &
               ibm0,iomld  ,isfc,il1k,il2k,iyr,imo,ida,ihr,                    &
               ifhour,ithr,0   ,inst,ina,inm,icen2,ids(iomld  ),iens,          &
#ifndef MRG_POST
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#else
               rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,       &
               g,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,g)
#ifdef DBG
       print *,'omld grib created    '
#endif
     else
       print *,'omld grib make failed'
     endif
#endif
   endif
! 
#ifndef NOPRINT
   if(iope) write(6,*)'grib flux file created.'
!
#endif
   if(iope) close(nn)
!
#ifdef RIVER
#define GRIBIT2  call file_make_grib_river(
!
! output for river files
!
#ifdef MP
   if(mype.eq.master) then
#endif
#ifndef RMP
     call file_name('rivgb',5,thour,fno,ncho)
#else
     call file_name('r_rivgb',7,thour,fno,ncho)
#endif
     open(unit=nn,file=fno(1:ncho),form='unformatted',err=902)
     go to 903
 902 continue
     write(6,*) ' error in opening file ',fno(1:ncho)
#ifdef MP
     call MPABORT
#else
     call abort
#endif
 903 continue
#ifndef NOPRINT
     write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
#endif
     rewind nn
#ifdef MP
   endif
#endif
!
! 1 riv storage
!
#ifdef MP
   if(mype.eq.master) then
#endif
     do j = 1,jo2_
       do i = 1,io2_
         if (gdriv(i,j).gt.1.E-10) then
           gdriv2(i,j)=log10(gdriv(i,j))
         else
           gdriv2(i,j)=99
         endif
       enddo
     enddo
     GRIBIT2 gdriv2,lbmr,0,io2_,jo2_,maxbit,90,                                &
                  ilpds,iptv,icen,igen,                                        &
                  ibm0,igdriv,isfc,il1k,il2k,iyr,imo,ida,ihr,                  &
                  ifhour,ithr,0,inst,ina,inm,icen2,ids(igdriv),iens,           &
                  rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,    &
                  gr,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,gr)
     else
       print *,'gdriv grib make failed'
     endif
#ifdef MP
   endif
#endif
!
! 2 riv discharge
!
#ifdef MP
   if(mype.eq.master) then
#endif
     do j = 1,jo2_
       do i = 1,io2_
         rflow(i,j)=rflow(i,j)*rtime
         if (rflow(i,j).gt.1.E-10) then
           rflow(i,j)=log10(rflow(i,j))
         else
           rflow(i,j)=99
         endif
       enddo
     enddo
     GRIBIT2 rflow,lbmr,0,io2_,jo2_,maxbit,90,                                 &
                  ilpds,iptv,icen,igen,                                        &
                  ibm0,irflow,isfc,il1k,il2k,iyr,imo,ida,ihr,                  &
                  ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irflow),iens,        &
                  rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,    &
                  gr,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,gr)
     else
       print *,'rflow grib make failed'
     endif
#ifdef MP
   endif
#endif
!
! 3 riv direction
!
#ifdef MP
   if(mype.eq.master) then
#endif
     GRIBIT2 imap,lbmr,0,io2_,jo2_,maxbit,90,                                  &
                  ilpds,iptv,icen,igen,                                        &
                  ibm0,iimap,isfc,il1k,il2k,iyr,imo,ida,ihr,                   &
                  ifhour,ithr,0,inst,ina,inm,icen2,ids(iimap),iens,            &
                  rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,    &
                  gr,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,gr)
     else
       print *,'imap grib make failed'
     endif
#ifdef MP
   endif
#endif
!
! 4 total runoff
!
#ifdef MP
   if(mype.eq.master) then
#endif
     do j = 1,jo2_
        do i = 1,io2_
          roff(i,j)=roff(i,j)*rtime
        enddo
     enddo
     GRIBIT2 roff,lbmr,0,io2_,jo2_,maxbit,90,                                  &
                  ilpds,iptv,icen,igen,                                        &
                  ibm0,iroff,isfc,il1k,il2k,iyr,imo,ida,ihr,                   &
                  ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iroff),iens,         &
                  rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,    &
                  gr,lg,ierr)
     if(ierr.eq.0) then
       call file_write_byte(nn,lg,gr)
     else
       print *,'roff grib make failed'
     endif
#ifdef MP
   endif
#endif
!
#ifdef MP
   if(mype.eq.master) then
#endif
     print *,'grib riv file created.'
#ifdef MP
   endif
#endif
!
#ifdef MP
   if(mype.eq.master) then
#endif
     do j = 1,jo2_
       do i = 1,io2_
         rflow(i,j)=0.e0
         roff(i,j)=0.e0
       enddo
     enddo
     close(nn)
#ifdef MP
   endif
#endif
#endif /* for RIVER */
!
#ifndef NOPRINT
   if(iope) write(6,*)'river file created.'
#endif
!
#ifndef RMP
#undef MRG_POST
#endif
!
#ifdef RIVER
   deallocate( lbmr, gr)
#endif
   deallocate( lbm, g, work,slmsep,fluxf)
!
   return
   end subroutine file_write_flux
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   subroutine rad_cloud_smooth(cvin,iin,jtwidl,jin,cvout,iout,jpout,jout,      &
                               xx,wgt,sum,nn,ltwidl,latrd1,latinb)                      
!-------------------------------------------------------------------------------
!   *  code bilinearly interpolates cld amt between gaussian grids--      *         
!   *  clone of rad_cloud_interp without the cloud top/base interpolation *
!-  *  j = 1 is just belo n.pole, i = 1 is greenwich (then go east).      *         
!   * iin,jin are i,j dimensions of input grid--iout,jout for output      *         
!   * jin2,jout2=jin/2,jout/2                                             *         
!   *                           campana+katz+campana(again) nov94         *         
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,jpout,jout,ltwidl,latrd1,latinb
   real                 ::  cvin(iin,jtwidl)                                                
   real                 ::  cvout(iout,jpout)                                               
   real                 ::  xx(iout,4),wgt(iout,4),sum(iout,4)        
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
       wgtlat=1.0
     endif
!
!     print 100,lat,xlat                                                      
!===>    if output lat is poleward of input lat=1 ,then simpl average           
!          (small region and cld amt wouldn t extrapolate well)                 
!
     call rad_cloud_smooth_sub(iii,jbb,jjj,iiiout,inslat,wgtlat,               &
               cvin,cvout(1,latout),xx,wgt,sum,nn,lbb,lr1)    
   enddo
! 100 format(1h ,' row =',i5,'  lat =',e15.5)                                   
   return                                                                    
   end subroutine rad_cloud_smooth
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   subroutine rad_cloud_smooth_sub(iin,jtwidl,jin,iout,inslat,wgtlat,cv,camt,  &
                      xx,wgt,sum,nn,ltwidl,latrd1)
!-------------------------------------------------------------------------------
!..        clone of rad_cloud_interp_sub without cldtop/base interpolation.....   
!        simpl linear interpolation of cldamt, unless only 1,2 of the           
!         surrounding pts has cv. then,if output gridpt not close enuf          
!         do not interpolate to it(prevents spreading of clds)..                
!           for 1 pt convection-intrp wgt ge (.7)**2 ...                        
!           for 2 pt convection-sum of intrp wgt ge .45...                      
!              .45 used rather than .5 to give better result for                
!              diagonally opposed pts...                                        
!         nn will be number of surrounding pts with cld (gt zero)               
!---     nhsh = 1,-1 for northern,southern hemisphere                           
!         here instead of an extrapolation,just do a simple mean....            
!                                                                               
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,inslat,ltwidl,latrd1
   real                 ::  cv(iin,jtwidl)                                                  
   real                 ::  camt(iout)                                                      
   real                 ::  xx(iout,4),wgt(iout,4),sum(iout,4)        
   integer              ::  nn(iout)                                                        
!
   lonf=iin/2
!
   if (inslat.lt.0) go to 600                                                
!
   inth = mod(ltwidl + inslat + jtwidl - latrd1 - 1,jtwidl) + 1              
   inth1 = mod(inth,jtwidl) + 1                                              
!
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
     camt(i) = 0.e0                                                     
     cycle
!
17   continue
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
     camt(i) = 0.e0                                                     
     cycle
!
27   continue                                                              
     camt(i) = sum(i,4)                                                  
!
   enddo
!
   return                                                                    
   end subroutine rad_cloud_smooth_sub
