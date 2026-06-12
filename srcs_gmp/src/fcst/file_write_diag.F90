#include "define.h"
   subroutine file_write_diag
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [file_write]
!      |
!      |--- [file_write_cld_grib] *
!      |--- [file_write_3d] *
!      |--- [file_write_slq] *
!
!-------------------------------------------------------------------------------
   end subroutine file_write_diag
!-------------------------------------------------------------------------------
!
#ifdef EXPLICIT_CLOUDINESS
!-------------------------------------------------------------------------------
   subroutine file_write_cld_grib(fhour,thour,idate,sl,colrad,                 &
        fluxr,cvavg,qcicps,qrscps,taucld,cldwp,cldip,ndg,lastep)
!-------------------------------------------------------------------------------
   use paramodel, only      : levs_,jcap_,LONF2S,LATG2S,LONF22S,lonf2_,latg2_
#ifdef DFS
   use dfsvar, only         : iope
#else
   use comio , only         : iope
#endif
#ifdef MP
   use commpi
#endif
   use comgda
   use module_trans, only   : dyn_trans2output_grid
   use diag_3d_module, only : diag_3d_get
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!     to write diagnostic data concerning with cloud & radiation
!     y.-h. byun               18 Aug 2004
!
!-------------------------------------------------------------------------------
!
#include "abort.h"
!
#ifndef RMP	
#define MPGP2F mpgp2f
#else 
#define MPGP2F rmpgp2f
#endif
!
! passing variables
!
   integer                               ::  idate(4)
   real                                  ::  fhour,thour
   real                                  ::  sl(levs_),colrad(latg2_)
   real   , dimension(LONF2S,LATG2S,26)  ::  fluxr
   real   , dimension(LONF2S,LATG2S)     ::  cvavg
   real                                  ::  qcicps(LONF2S,levs_,LATG2S),      &
                                             qrscps(LONF2S,levs_,LATG2S)
   real                                  ::  cldwp(LONF2S,levs_,LATG2S),       &
                                             cldip(LONF2S,levs_,LATG2S)
   real                                  ::  taucld(LONF2S,levs_,LATG2S)
   logical                               ::  lastep
!
   integer,parameter    ::  iprs=1,itemp=11,iznlw=33,imerw=34,isphum=51
   integer,parameter    ::  ipcpr=59,isnowd=65,icldf=71,iccldf=72
   integer,parameter    ::  islmsk=81,iz0cm=83,ialbdo=84,isoilm=144,icemsk=91
   integer,parameter    ::  ilhflx=121,ishflx=122,izws=124,imws=125,ighflx=155
   integer,parameter    ::  iuswfc=160,idswfc=161,iulwfc=162,idlwfc=163
   integer,parameter    ::  inswfc=164,inlwfc=165
   integer,parameter    ::  idswvb=166,idswvd=167,idswnb=168,idswnd=169
   integer,parameter    ::  isglyr=175,icnpy=145
   integer,parameter    ::  idswf=204,idlwf=205,iuswf=211,iulwf=212,icpcpr=214
   integer,parameter    ::  isfc=1,itoa=8,ielev=105
   integer,parameter    ::  isglev=107,idbls=111,i2dbls=112,icolmn=200
   integer,parameter    ::  ilcbl=212,ilctl=213,ilclyr=214
   integer,parameter    ::  imcbl=222,imctl=223,imclyr=224
   integer,parameter    ::  ihcbl=232,ihctl=233,ihclyr=234
   integer,parameter    ::  inst=10,iavg=3,iacc=4
   integer,parameter    ::  ifhour=1,ifday=2
   integer              ::  lonb,latb
   logical              ::  lbm(lonf2_,latg2_)
   character            ::  g(200+lonf2_*latg2_*(16+1)/8)
   integer              ::  ids(255)
   integer              ::  iens(5)
   integer              ::  iclyr(3),ictl(3),icbl(3),itlcf(3)
   data iclyr/ihclyr,imclyr,ilclyr/
   data ictl /ihctl ,imctl ,ilctl /
   data icbl /ihcbl ,imcbl ,ilcbl /
   data itlcf/itoa,isfc,icolmn/
   real,save            ::  phour
   data phour/0.0/
#ifdef ASSIGN
   character(len=120)   ::  asgnstr
#endif
!
   integer              ::  idrt,ncho,ndg,iyr,imo,ida,ihr,iftime,ifhr,ithr,isl,iensi,ienst
   real                 ::  dhour
   integer              ::  nvar,k,j,i,ipunv,idsnv,ierr,lg
   real                 ::  work(lonf2_,latg2_),slmsep(lonf2_,latg2_)
   real                 ::  work2(lonf2_,latg2_,levs_)
   real                 ::  fluxw(lonf2_,latg2_,26)
#ifdef MP
   real                 ::  work2p(LONF2S,LATG2S,levs_)
#endif
   integer,parameter    ::  icen=7,icen2=0,igen=99
   real                 ::  cl1
!
   character(len=80)    ::  fno
!
! initialize local variables
!
#ifdef MP
   work2p=0.
#endif
!
   lonb=lonf2_/2
   latb=latg2_*2
#ifdef DFS
   idrt=0
#else
   idrt=4
#endif
   cl1=colrad(1)
!
   if(iope) then
     print *,' start file_write_3d '
     call file_name('cpsgb',5,thour,fno,ncho)
#ifdef ASSIGN
     write(asgnstr,'(23hassign -s unblocked  u:,I2,)') ndg
     call assign('assign -R')
     call assign(asgnstr)
#endif
     open(unit=ndg,file=fno(1:ncho),form='unformatted',err=900)
!
     go to 901
!
900  continue
!
     write(6,*) ' error in opening file ',fno(1:ncho)
     call MPABORT
!
901  continue
!
#ifndef NOPRINT
     write(6,*) ' file ',fno(1:ncho),' opened. unit=',ndg
#endif
     rewind ndg
   endif  ! iope
!
   if( phour.eq.0.0 ) phour=fhour
!
   call idsdef(1,ids)
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
   iftime=ifhour
   ifhr=nint(phour)
   ithr=nint(thour)
   dhour=thour-phour
!
   if(dhour.gt.0.) then
     if(iope) then
       if(lastep) then
         call file_name('cpscld',6,thour,fno,ncho)
         close(44)
         open(44,file=fno(1:ncho),form='unformatted')
       endif
     endif
!
     do nvar = 1,5
       do k = 1,levs_
         do j = 1,LATG2S
           do i = 1,LONF2S
#ifdef MP
             if (nvar.eq.1) work2p(i,j,k)=qcicps(i,k,j)
             if (nvar.eq.2) work2p(i,j,k)=qrscps(i,k,j)
             if (nvar.eq.3) work2p(i,j,k)=taucld(i,k,j)
             if (nvar.eq.4) work2p(i,j,k)=cldwp (i,k,j)
             if (nvar.eq.5) work2p(i,j,k)=cldip (i,k,j)
#else
             if (nvar.eq.1) work2(i,j,k)=qcicps(i,k,j)
             if (nvar.eq.2) work2(i,j,k)=qrscps(i,k,j)
             if (nvar.eq.3) work2(i,j,k)=taucld(i,k,j)
             if (nvar.eq.4) work2(i,j,k)=cldwp (i,k,j)
             if (nvar.eq.5) work2(i,j,k)=cldip (i,k,j)
#endif
           enddo
         enddo
       enddo
#ifdef MP
       call MPGP2F(work2p,LONF2S,LATG2S,work2,lonf2_,latg2_,levs_)
#endif
!
       if(iope) then
!
         if(lastep) then
           write(44) (((work2(i,j,k),i=1,lonf2_),j=1,latg2_),k=1,levs_)
         endif
!
         ipunv = 185 + nvar - 1
         idsnv = 3
         if (nvar.le.2) then
           idsnv = 4
         elseif (nvar.ge.4) then
           idsnv = 1
         endif
!
         call dyn_trans2output_grid(work2,levs_)
!
         do k = 1,levs_
           isl=nint(sl(k)*1.e4)
           call file_make_grib(work2(1,1,k),lbm,idrt,lonb,latb,16,cl1,28,2,    &
                icen,igen,                                                     &
                0,ipunv,isglev,0,isl,iyr,imo,ida,ihr,                          &
                iftime,ifhr,ithr,iavg,0,0,icen2,idsnv,iens,                    &
                0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)
           if(ierr.eq.0) call file_write_byte(ndg,lg,g)
         enddo
!
       endif   ! iope
     enddo     ! nvar
!
! yhb
!
   endif ! dhour.gt.0.
!
   phour=thour
!
   if( iope ) then
     close(ndg)
   endif
!
#undef MPGP2F
!
   return
   end subroutine file_write_cld_grib
#endif /* EXPLICIT_CLOUDINESS end */
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine file_write_3d(slmask,ndg)
#ifdef DG3
!-------------------------------------------------------------------------------
   use paramodel, only      : LONF2S,LATG2S,LONF2F,LATG2F,latg_,latg2_,levs_,lonf_,lonf2_
#ifdef DFS
   use dfsvar, only         : iope,sl=>sigma
#else
   use comio , only         : iope
#endif
#ifdef RMP
   use paramodel, only      : igrd1_,jgrd1_
#ifdef MP
   use paramodel, only      : igrd12_,jgrd12_, igrd12p_,jgrd12p_
   use commpi
#endif
   use rscomf_rerun
   use rscommap
   use module_trans, only   : rmp_trans2model_grid, rmp_trans2output_grid
#else
   use comfphys
   use comfspec_vr
   use comfver
   use comfgrid
   use radiag
#if defined(MP) && !defined(DFS)
   use commpi
#endif
   use module_trans, only   : dyn_trans2model_grid, dyn_trans2output_grid
#endif
   use comsfc
   use comgda
   use comfcst, only        : dfvbr, dfnbr, dfvdr, dfndr
   use diag_3d_module, only : diag_3d_get
!
#ifndef RMP
#define ILEN lonf_
#define JLEN latg_
#ifdef MP
#define MPGP2F mpgp2f
#define MPABORT mpabort
#else         /* not MP */
#define MPABORT abort
#endif         /* MP */
#else         /* RMP */
#define ILEN igrd1_
#define JLEN jgrd1_
#ifdef MP
#define MPGP2F rmpgp2f
#define MPABORT rmpabort
#else         /* MP */
#define MPABORT abort
#endif         /* MP */
#endif         /* RMP */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   integer                                          ::  ndg
   real   , dimension(LONF2S,LATG2S)                ::  slmask
!
! local variables
!
#ifdef CLR
   real   , dimension(LONF2F,LATG2F,3)              ::  cfsw,cflw
#endif
   real   , dimension(LONF2F,LATG2F)                ::  work,slmsep
   real   , dimension(LONF2F,LATG2F,26)             ::  fluxw
!   real   , dimension(LONF2F,LATG2F,levs_+2/levs_)  ::  work2
   real   , dimension(LONF2F,LATG2F,levs_)          ::  work2
#ifdef MP
!   real   , dimension(LONF2S,LATG2S,levs_+2/levs_)  ::  work2p
   real   , dimension(LONF2S,LATG2S,levs_)           ::  work2p
#endif
   real                                             ::  gda(nwgda)
   integer,parameter    ::  iprs=1,itemp=11,iznlw=33,imerw=34,isphum=51,       &
                            ipcpr=59,isnowd=65,icldf=71,iccldf=72,             &
                            islmsk=81,iz0cm=83,ialbdo=84,isoilm=144,icemsk=91, &
                            ilhflx=121,ishflx=122,izws=124,imws=125,ighflx=155,&
                            iuswfc=160,idswfc=161,iulwfc=162,idlwfc=163,       &
                            inswfc=164,inlwfc=165,                             &
                            idswvb=166,idswvd=167,idswnb=168,idswnd=169,       &
                            isglyr=175,icnpy=145,                              &
                            idswf=204,idlwf=205,iuswf=211,iulwf=212,icpcpr=214
   integer,parameter    ::  isfc=1,itoa=8,ielev=105,                           &
                            isglev=107,idbls=111,i2dbls=112,icolmn=200,        &
                            ilcbl=212,ilctl=213,ilclyr=214,                    &
                            imcbl=222,imctl=223,imclyr=224,                    &
                            ihcbl=232,ihctl=233,ihclyr=234
   integer,parameter    ::  inst=10,iavg=3,iacc=4
   integer,parameter    ::  ifhour=1,ifday=2
   integer              ::  lonb, latb
!   logical,allocatable    ::  lbm(:,:)
!   character,allocatable  ::  g(:)
   logical              ::  lbm(LONF2F,LATG2F)
   character            ::  g(200+LONF2F*LATG2F*(16+1)/8)
   integer              ::  ids(255)
   integer              ::  iens(5)
   integer              ::  iclyr(3),ictl(3),icbl(3),itlcf(3)
   data iclyr/ihclyr,imclyr,ilclyr/
   data ictl /ihctl ,imctl ,ilctl /
   data icbl /ihcbl ,imcbl ,ilcbl /
   data itlcf/itoa,isfc,icolmn/
#ifdef MP
   real, allocatable    :: fullat(:),fullon(:)
#endif
!
#ifndef RMP
   integer,parameter    ::  icen=7,icen2=0,igen=99
#endif
!
#ifdef ASSIGN
   character(len=120)   ::  asgnstr
#endif
!
   integer              ::  ncho,ienst,iensi,iyr,imo,ida,ihr,iftime,ifhr,ithr,&
                            idrt,i,j,k,kd,kgda,igda,ipu,ibm,isl,ibitmap,lg,ierr
   real                 ::  cl1,dhour,rtime,rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,&
                            truth,cotru
!-------------------------------------------------------------------------------
   character(len=128)   ::   fno
!
! initialize local variables
!
#ifdef CLR
   cfsw=0. ; cflw=0.
#endif
   gda=0.  ; ids=0 ; iens=0
   work=0. ; slmsep=0. ; work2=0. ; fluxw=0.
#ifdef MP
   work2p=0.
#endif
!
   lonb=lonf2_/2
   latb=latg2_*2
!
#ifdef SMP
   colrad(:)=90.
#else
   cl1=colrad(1)
#endif
!
#ifdef MP
   if(iope) then
     print *,' start file_write_3d '
#endif
#ifndef RMP
     call file_name('dg3gb',5,thour,fno,ncho)
#else
     call file_name('r_dg3gb',7,thour,fno,ncho)
#endif
#ifdef ASSIGN
     write(asgnstr,'(23hassign -s unblocked  u:,I2,)') ndg
     call assign('assign -R')
     call assign(asgnstr)
#endif
     open(unit=ndg,file=fno(1:ncho),form='unformatted',err=900)
!
     go to 901
!
900  continue
     write(6,*) ' error in opening file ',fno(1:ncho)
     call MPABORT
!
901  continue
#ifndef NOPRINT
     write(6,*) ' file ',fno(1:ncho),' opened. unit=',ndg
#endif
     rewind ndg
#ifdef MP
   endif   ! iope
#endif
!
   call idsdef(1,ids)
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
   iftime=ifhour
   ifhr=nint(fhour)
   ithr=nint(thour)
!
! -- Never show 0-0hr averages in wgrib label, always round up to 1 hour
!
   if(ithr.eq.0) ithr=1
!
   dhour=dtpost/3600.
   if(dtpost.gt.0.) then
     rtime=1./dtpost
   else
     rtime=0.
   endif
!
#ifdef RMP
   colrad=0.
   proj=rproj
   delx=rdelx
   dely=rdely
#ifdef MP
   allocate (fullat(ILEN*JLEN),fullon(ILEN*JLEN))
   call MPGP2F(flat,igrd12p_,jgrd12p_,fullat,igrd12_,jgrd12_,1)
   call MPGP2F(flon,igrd12p_,jgrd12p_,fullon,igrd12_,jgrd12_,1)
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
!
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
#else
#ifdef SMP
   idrt=0
#else
#ifdef DFS
   idrt=0
#else
   idrt=4
#endif
#endif
   rlat1=0.
   rlon1=0.
   rlat2=0.
   rlon2=0.
   delx=0.
   dely=0.
   ortru=0.
   proj=0.
#endif /* end RMP */
!
#ifdef MP
   call MPGP2F(slmask,LONF2S,LATG2S,slmsep,LONF2F,LATG2F,1)
!
   if(iope) then
#else
   do j = 1,LATG2F
     do i = 1,LONF2F
       slmsep(i,j)=slmask(i,j)
     enddo
   enddo
!
#endif
#ifndef RMP
   call dyn_trans2output_grid(slmsep,1)
#else
   call rmp_trans2output_grid(slmsep,1)
#endif
#ifdef MP
   endif
#endif
!
#ifdef MP
   call MPGP2F(fluxr,LONF2S,LATG2S,fluxw,LONF2F,LATG2F,26)
!
   if(iope) then
#else
     do k = 1,26
       do j = 1,LATG2F
         do i = 1,LONF2F
           fluxw(i,j,k)=fluxr(i,j,k)
         enddo
       enddo
     enddo
!
#endif
#ifndef RMP
     call dyn_trans2output_grid(fluxw,26)
#else
     call rmp_trans2output_grid(fluxw,26)
#endif
#ifdef CLR
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         work(i,j)=fluxw(i,j,21)*rtime
       enddo
     enddo
!
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'fluxw-21')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,iulwfc,itoa,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(iulwfc),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         work(i,j)=fluxw(i,j,22)*rtime
       enddo
     enddo
!
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'fluxw-22')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,iuswfc,itoa,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(iuswfc),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         work(i,j)=fluxw(i,j,25)*rtime
       enddo
     enddo
!
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'fluxw-25')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,idlwfc,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(idlwfc),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         work(i,j)=fluxw(i,j,23)*rtime
       enddo
     enddo
!
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'fluxw-23')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,idswfc,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswfc),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         work(i,j)=fluxw(i,j,24)*rtime
       enddo
     enddo
!
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'fluxw-24')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,iuswfc,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(iuswfc),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)                                                    
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
!   compute sw cloud forcing at toa (cld-clear),
!   flip sign so positive means cld is warming relative to clear
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         cfsw(i,j,1) = - (fluxw(i,j,2)-fluxw(i,j,22))
       enddo
     enddo 
!
!   compute cloud forcing at sfc (cld-clear)
!   again flip sign so positive means cld is warming relative to clear
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         cfsw(i,j,2) =-(fluxw(i,j,3)-fluxw(i,j,4)                              &
                         -(fluxw(i,j,24)-fluxw(i,j,23)))
       enddo 
     enddo
!
!   flip sign so positive means cld is warming relative to clear
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         cfsw(i,j,3) = - (cfsw(i,j,2) - cfsw(i,j,1))
       enddo 
     enddo
!
     do k = 1,3
       do j = 1,LATG2F
         do i = 1,LONF2F
           work(i,j)=cfsw(i,j,k)*rtime
         enddo
       enddo
!
#ifdef DBG
       call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'cfsw')
#endif
       call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,   &
                   0,inswfc,itlcf(k),0,0,iyr,imo,ida,ihr,                      &
                   iftime,ifhr,ithr,iavg,0,0,icen2,ids(inswfc),iens,           &
                   rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,               &
                   truth,cotru,                                                &
                   g,lg,ierr)
       if(ierr.eq.0) call file_write_byte(ndg,lg,g)
     enddo 
!
!   compute lw cloud forcing at toa (cld-clear),
!   flip sign so positive means cld is warming relative to clear
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         cflw(i,j,1) = - (fluxw(i,j,1)-fluxw(i,j,21))
       enddo 
     enddo
!
!    compute cloud forcing at sfc (cld-clear)
!    again flip sign so positive means cld is warming relative to clear
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         cflw(i,j,2) = - (fluxw(i,j,25)-fluxw(i,j,19))
       enddo
     enddo
!
!   flip sign so positive means cld is warming relative to clear
!
     do j = 1,LATG2F
       do i = 1,LONF2F
         cflw(i,j,3) = - (cflw(i,j,2) - cflw(i,j,1))
       enddo
     enddo 
!
     do k = 1,3
       do j = 1,LATG2F 
         do i = 1,LONF2F
           work(i,j)=cflw(i,j,k)
           work(i,j)=work(i,j)*rtime
         enddo
       enddo
#ifdef DBG
       call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'cflw')
#endif
       call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,   &
                   0,inlwfc,itlcf(k),0,0,iyr,imo,ida,ihr,                      &
                   iftime,ifhr,ithr,iavg,0,0,icen2,ids(inlwfc),iens,           &
                   rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,               &
                   truth,cotru,                                                &
                   g,lg,ierr)
       if(ierr.eq.0) call file_write_byte(ndg,lg,g)
     enddo 
#endif         /* CLR */
!
#ifdef MP
   endif ! iope
!
   call MPGP2F(dfvbr,LONF2S,LATG2S,work,LONF2F,LATG2F,1)
!
   if(iope) then
#endif
!
     do j = 1,LATG2F
       do i = 1,LONF2F
#ifdef MP
         work(i,j)= work(i,j)*rtime
#else
         work(i,j)=dfvbr(i,j)*rtime
#endif
       enddo
     enddo
!
#ifndef RMP
     call dyn_trans2output_grid(work,1)
#else
     call rmp_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'dfvbr')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,idswvb,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswvb),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
#ifdef MP
   endif ! mype.eq.master
!
   call MPGP2F(dfvdr,LONF2S,LATG2S,work,LONF2F,LATG2F,1)
!
   if(iope) then
#endif
!
     do j = 1,LATG2F
       do i = 1,LONF2F
#ifdef MP
         work(i,j)= work(i,j)*rtime
#else
         work(i,j)=dfvdr(i,j)*rtime
#endif
       enddo
     enddo
!
#ifndef RMP
     call dyn_trans2output_grid(work,1)
#else
     call rmp_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'dfvdr')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,idswvd,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswvd),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
#ifdef MP
   endif ! mype.eq.master
!
   call MPGP2F(dfnbr,LONF2S,LATG2S,work,LONF2F,LATG2F,1)
!
   if(iope) then
#endif
!
     do j = 1,LATG2F
       do i = 1,LONF2F
#ifdef MP
         work(i,j)= work(i,j)*rtime
#else
         work(i,j)=dfnbr(i,j)*rtime
#endif
       enddo
     enddo
!
#ifndef RMP
     call dyn_trans2output_grid(work,1)
#else
     call rmp_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'dfnbr')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,idswnb,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswnb),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
#ifdef MP
   endif ! mype.eq.master
!
   call MPGP2F(dfndr,LONF2S,LATG2S,work,LONF2F,LATG2F,1)
!
   if(iope) then
#endif
!
     do j = 1,LATG2F
       do i = 1,LONF2F
#ifdef MP
         work(i,j)= work(i,j)*rtime
#else
         work(i,j)=dfndr(i,j)*rtime
#endif
       enddo
     enddo
!
#ifndef RMP
     call dyn_trans2output_grid(work,1)
#else
     call rmp_trans2output_grid(work,1)
#endif
#ifdef DBG
     call print_maxmin_six(work,LONF2S*LATG2S,1,1,1,'dfndr')
#endif
     call file_make_grib(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,     &
                 0,idswnd,isfc,0,0,iyr,imo,ida,ihr,                            &
                 iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswnd),iens,             &
                 rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,                 &
                 truth,cotru,                                                  &
                 g,lg,ierr)
     if(ierr.eq.0) call file_write_byte(ndg,lg,g)
!
#ifdef MP
   endif ! mype.eq.master
!
#endif
!
   do kd = 1,kdgda
!
     do j = 1,nrgda
       kgda=(j-1)*kdgda+kd
       call diag_3d_get(kgda,nwgda,gda)
       igda=0
       do k = 1,levs_
         do i = 1,LONF2S
           igda=igda+1
#ifdef MP
           work2p(i,j,k)=gda(igda)*rtime
#else
           work2(i,j,k)=gda(igda)*rtime
#endif
         enddo
       enddo
     enddo 
!  
#ifdef MP
     call MPGP2F(work2p,LONF2S,LATG2S,work2,LONF2F,LATG2F,levs_)
!
     if(iope) then
#endif
#ifndef RMP
       call dyn_trans2output_grid(work2,levs_)
#else
       call rmp_trans2output_grid(work2,levs_)
#endif
       ipu=ipugda(kd)
       ibm=ibmgda(kd)
!
       do k = 1,levs_
         isl=nint(sl(k)*1.e4)
#ifdef DBG
         print *,'kd=',kd,' k=',k
         call print_maxmin_six(work2(1,1,k),LONF2F*LATG2F,1,1,1,'diags')
#endif
         ibitmap=0
!
         if(ibm.ne.0) then
           ibitmap=1
           do j = 1,LATG2F
             do i = 1,LONF2F
               lbm(i,j)=work2(i,j,k).ne.0.
             enddo
           enddo
         endif
!
         call file_make_grib(work2(1,1,k),lbm,idrt,lonb,latb,16,cl1,28,2,   &
                   icen,igen,                                                  &
!    &            ibitmap,ipu,isglev,0,isl,iyr,imo,ida,ihr,                    &
                   0,ipu,isglev,0,isl,iyr,imo,ida,ihr,                         &
                   iftime,ifhr,ithr,iavg,0,0,icen2,ids(ipu),iens,              &
                   rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,               &
                   truth,cotru,                                                &
                   g,lg,ierr)
         if(ierr.eq.0) call file_write_byte(ndg,lg,g)
       enddo
!
#ifdef MP
     endif ! mype.eq.master
!
#endif
   enddo   
!
#ifdef MP
   if(iope) then
#endif
     close(ndg)
#ifdef MP
   endif
#endif
#undef ILEN
#undef JLEN
#undef MPGP2F
#undef MPABORT
#undef FLAT
#undef FLON
!
   return
#endif /* DG3 end */
   end subroutine file_write_3d
!-------------------------------------------------------------------------------
#ifdef NISLQ_GRIB
!
!-------------------------------------------------------------------------------
   subroutine file_write_slq(idate,itpdt,thour)
!-------------------------------------------------------------------------------
   use paramodel, only      : levs_,LONF2S,LATG2S,lonf2_,latg2_,ntotal_       ,&
                              nwmass_,nwsize_,ngases_
#ifdef DFS
   use dfsvar, only         : iope,sl=>sigma
#else
   use comio , only         : iope
   use comfver, only        : sl
#endif
   use comfgrid, only       : colrad
#ifdef MP
   use commpi
#endif
   use comgda
   use module_trans, only   : dyn_trans2output_grid
   use nislq, only : slq_q1,slq_q2
#ifdef DCMIP
   use dcmip_grims, only : dcmip_grims_init,icase,h
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!     to write moisture in case of nislq
!       Myung-Seo Koo           2 Feb 2015
!
!-------------------------------------------------------------------------------
#include "abort.h"
#ifndef RMP	
#define MPGP2F mpgp2f
#else 
#define MPGP2F rmpgp2f
#endif
!
! passing variables
!
   integer                               ::  idate(4),itpdt
   real                                  ::  fhour,thour
!
! local variables
!
   integer, dimension(ntotal_)           ::  ipu
   integer, parameter                    ::  ndg=11
   integer, parameter                    ::  iavg=3,ifhour=1                  ,&
                                             icen=7,icen2=0,igen=99           ,&
                                             ipuq=51                          ,&
                                             ipuqc=153                        ,&
                                             ipuqr=152                        ,&
                                             ipuqi=151                        ,&
                                             ipuqs=150                        ,&
                                             ipuqg=255                        ,&
                                             ipunccn=181                      ,&
                                             ipunc=182                        ,&
                                             ipunr=183                        ,&
                                             ipuo3=154                        ,&
                                             iputke=158
   logical                               ::  lpqc,lpqr,lpqi,lpqs,lpqg,lpnccn  ,&
                                             lpnc,lpnr,lpo3,lptke
#ifdef DFS
   integer, parameter                    ::  idrt=0   ! uniform grid
#else
   integer, parameter                    ::  idrt=4   ! gaussian grid
#endif
   logical                               ::  lastep
!
   integer                               ::  isglev,lonb,latb,kk
   logical                               ::  lbm(lonf2_,latg2_)
   character                             ::  g(200+lonf2_*latg2_*(16+1)/8)
   integer                               ::  ids(255)
   integer                               ::  iens(5)
#ifdef ASSIGN
   character(len=120)   ::  asgnstr
#endif
!
   integer              ::  ksq,ksqc,ksqr,ksqi,ksqs,ksqg,ksnccn,ksnc,ksnr,kso3,kstke
   integer              ::  ncho,iftime,ifhr,ithr,isl,iensi,ienst,ki
   integer              ::  nvar,k,j,i,ipunv,ierr,lg
   real                 ::  dhour
   real   , dimension(lonf2_,latg2_)        ::  work,slmsep
   real   , dimension(lonf2_,latg2_,levs_)  ::  work2
#ifdef MP
   real   , dimension(LONF2S,LATG2S,levs_)  ::  work2p
#endif
   real                 ::  cl1
   character(len=80)    ::  fno
#ifdef DCMIP
   real                ::  dum1(levs_+1),dum2(levs_+1)
#endif
!
! model layer index (eta or sig)
!
  isglev=107
!
! initialize local variables
!
#ifdef MP
   work2p=0.
#endif
   lonb=lonf2_/2
   latb=latg2_*2
   cl1=colrad(1)
!
#ifdef DCMIP
! dcmip setting
!
   call dcmip_grims_init(dum1,dum2)
   if(icase.eq.11 .or. icase.eq.12) isglev=105 ! height
   ipu(1)=51
   ipu(2)=154
   ipu(3)=158
   ipu(4)=196
   ipu(5)=197
   ipu(6)=198
   ipu(7)=199
#else
! logical
!
   lpqc   = nwmass_.ge.2
   lpqr   = nwmass_.ge.3
   lpqi   = nwmass_.ge.4
   lpqs   = nwmass_.ge.5
   lpqg   = nwmass_.ge.6
   lpnccn = nwsize_.gt.0
   lpnc   = nwsize_.gt.0
   lpnr   = nwsize_.gt.0
   lpo3   = ngases_.ge.1
   lptke  = ngases_.ge.2
!
   ki=1
   ipu(ki)=ipuq    ; ki=ki+1
   if(lpqc) then
     ipu(ki)=ipuqc ; ki=ki+1
   endif
   ! wsm3
   if(lpqr) then
     ipu(ki)=ipuqr ; ki=ki+1
   endif
   ! wsm5
   if(lpqi) then
     ipu(ki)=ipuqi ; ki=ki+1
     ipu(ki)=ipuqs ; ki=ki+1
   endif
   ! wsm6
   if(lpqg) then
     ipu(ki)=ipuqg ; ki=ki+1
   endif
!
! --  wdm 5 and 6
!
   if(lpnccn) then
     ipu(ki)=ipunccn ; ki=ki+1
     ipu(ki)=ipunc   ; ki=ki+1
     ipu(ki)=ipunr   ; ki=ki+1
   endif
!
! --  o3
!
   if(lpo3) then
     ipu(ki)=ipuo3   ; ki=ki+1
   endif
!
! --  tke
!
   if(lptke) then
     ipu(ki)=iputke  ; ki=ki+1
   endif
#endif
!
! start
!
   if(iope) then
     print *,' start file_write_slqgb'
     call file_name('slqgb',5,thour,fno,ncho)
#ifdef ASSIGN
     write(asgnstr,'(23hassign -s unblocked  u:,I2,)') ndg
     call assign('assign -R')
     call assign(asgnstr)
#endif
     open(unit=ndg,file=fno(1:ncho),form='unformatted',err=900)
!
     go to 901
!
900  continue
!
     write(6,*) ' error in opening file ',fno(1:ncho)
     call MPABORT
!
901  continue
!
#ifndef NOPRINT
     write(6,*) ' file ',fno(1:ncho),' opened. unit=',ndg
#endif
     rewind ndg
   endif  ! iope
!
   call idsdef(1,ids)
#ifdef DCMIP
   ids(196)=6
   ids(197)=6
   ids(198)=6
   ids(199)=6
#endif
   ienst=0
   iensi=0
   iens(1)=1
   iens(2)=ienst
   iens(3)=iensi
   iens(4)=1
   iens(5)=255
   iftime=ifhour
   fhour=0.
   ifhr=nint(fhour)
   ithr=nint(thour)
!
#ifdef MP
#define WORK2 WORK2P
#endif
   do nvar = 1,ntotal_
     do k = 1,levs_
       kk=(nvar-1)*levs_+k
       do j = 1,LATG2S
#ifdef DFS
         do i = 1,LONF2S/2
           WORK2(         i,j,k)=slq_q1(i,         j,kk)
           WORK2(LONF2S/2+i,j,k)=slq_q1(i,LATG2S+1-j,kk)
         enddo
#else
         do i = 1,LONF2S
           WORK2(i,j,k)=slq_q1(i,kk,j)
         enddo
#endif /* DFS end */
#undef WORK2
       enddo
     enddo
#ifdef MP
     call MPGP2F(work2p,LONF2S,LATG2S,work2,lonf2_,latg2_,levs_)
#endif
!
     if(iope) then
!
       call dyn_trans2output_grid(work2,levs_)
       ipunv = ipu(nvar)
       do k = 1,levs_
         isl=nint(sl(k)*1.e4)
#ifdef DCMIP
         if(icase.eq.11.or.icase.eq.12) isl=H*log(1./sl(k))+1.
#endif
         call file_make_grib(work2(1,1,k),lbm,idrt,lonb,latb,16,cl1,28,2,         &
              icen,igen,0,ipunv,isglev,0,isl,idate(4),idate(2),idate(3),idate(1), &
              iftime,ifhr,ithr,iavg,0,0,icen2,ids(ipunv),iens,                    &
              0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)
         if(ierr.eq.0) call file_write_byte(ndg,lg,g)
       enddo
!
     endif   ! iope
   enddo     ! nvar
!
   if( iope ) then
     close(ndg)
   endif
!
#undef MPGP2F
!
   return
   end subroutine file_write_slq
!-------------------------------------------------------------------------------
#endif /* NISLQ_GRIB end */
