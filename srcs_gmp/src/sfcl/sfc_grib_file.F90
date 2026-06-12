#include "define.h"
   subroutine sfc_grib_file
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [sfc_grib_file]
!           |
!           |-- [sfc_grib_driver] *
!           |-- [sfc_grib_solver] *
!           |-- [sfc_grib_numchar] *
!           |-- [sfc_read_grib] *
!           |-- [sfc_read_mask] *
!           |-- [sfc_read_mtn] *
!           |-- [sfc_field2mask] *
!           |-- [sfc_get_area_gmp] *
!           |-- [sfc_get_area_rmp] *
!           |-- [sfc_interp_date] *
!
!-------------------------------------------------------------------------------
   end subroutine sfc_grib_file
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_grib_driver(lugb,idim,jdim,                                  &
                       numsfcs,numsfcv,ksfc,                                   &
                       is2g,                                                   &
                       iy,im,id,ih,fh,                                         &
                       fnmask,fnorog,fnmskg,fn,numgrbs,                        &
                       orog,slmask,                                            &
                       grbfld,lsf)
!-------------------------------------------------------------------------------
!
! subroutine: sfc_grib_driver
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
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,latg_,npes_
#if defined(MP)
   use commpi, only : mype,master
#endif
!
!  read surface grib files
!
!    lugb:  integer unit number 
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                 ::  lugb,idim,jdim,numsfcs,numsfcv
   integer                 ::  iy,im,id,ih
   integer                 ::  numgrbs
   real                    ::  fh
   real                    ::  orog(idim*jdim),slmask(idim*jdim)
   real                    ::  grbfld(idim*jdim,numsfcs)
   character(len=128)      ::  fnmask,fnorog,fnmskg
   character(len=128)      ::  fn(numgrbs)
   integer                 ::  is2g(numsfcs)
!
   logical                 ::  lsf(numsfcs)
   integer                 ::  ksfc(numsfcs)
!
   integer                 ::  nsfc,nsfcv,ngrb
   integer                 ::  nch
   integer                 ::  ifp
   data                        ifp/0/
   save                        ifp
!
!  read fixed fiedls once at first pass and at fh=0.
!
   if(ifp.eq.0) then
!
!  get land-sea mask on model from mtn program binary output.
!
#ifdef MP
     if(mype.eq.master) then
#endif
       call sfc_read_mask(lugb,fnmask,idim*jdim,slmask)
#ifdef MP
     endif
#endif
!
!  get orography on model grid from mtn program binary output.
!
#ifdef MP
     if(mype.eq.master) then
#endif
       call sfc_read_mtn(lugb,fnorog,idim*jdim,slmask,orog)
#ifdef MP
     endif
#endif
     ifp=1
   endif
!
!   read all other grib files
!   note that grbfld is in sfc rec order
!
   nsfc=1
   do nsfcv = 1,numsfcv
     ngrb=is2g(nsfc)
#if defined(DBG) && !defined(MP)
     if(ngrb.ne.9999) then
       call sfc_grib_numchar(fn(ngrb),nch)
       write(6,100)nsfcv,nsfc,fn(ngrb)(1:nch),ngrb,ksfc(nsfc),numsfcs
       call flush(6)
100    format('sfc_grib_driver:nsfcv=',i3,' nsfc=',i3,                         &
                ' fn=',a50,' ngrb=',i4,' ksfc(nsfc)=',i3,' num sfc=',i3)
     endif
#endif
     if(ngrb.ne.9999) then
       if(fn(ngrb)(1:4).ne.'    ') then
         lsf(nsfc)=.true.
         call sfc_grib_solver(lugb,fn(ngrb),idim,jdim,slmask,                  &
                          ngrb,iy,im,id,ih,fh,                                 &
                          fnmskg,grbfld(1,nsfc))
       else
         lsf(nsfc)=.false.
       endif
     endif
     nsfc=nsfc+ksfc(nsfc)
   enddo
!
   return
   end subroutine sfc_grib_driver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_grib_solver(lugb,fn,idim,jdim,slmask,                        &
                       indx,iy,im,id,ih,fh,                                    &
                       fnmskg,out)
!-------------------------------------------------------------------------------
!
!  read grib file, interpolates in time and to model grid.
!  complex interpolation by using masks are applied
!  control of interpolation is specified in inlude/vargrb.h
!
!  input:
!
!    lugb    int                   .. unit number for reading grib
!    fn      char                  .. grib file name
!    idim    int                   .. output x-dimension (model x-dimension)
!    jdim    int                   .. output y-dimension (model y-dimension)
!    slmask  real                  .. land sea mask on model grid
!    iy      int                   .. 4-digit year
!    im      int                   .. month
!    id      int                   .. day
!    ih      int                   .. hour
!    fh      real                  .. forecast hour
!    indx    int                   .. pointer of grb fields (see vargrb.h)
!    fnmskg character             .. hi-res lat/lon land sea global mask 
!
!  output:
!   
!    out    real array (idim*jdim,kdim) .. output field
!    iclmtyp integer               .. =1 yearly mean, 4=seasonal mean,
!                                     =12 monthly mean, -1=not a climatology
!
!  the map related parameters appear in this program are for grib
!  files and not those of the forecast program
!
!-------------------------------------------------------------------------------
   use vargrb
#if defined(MP)
   use commpi, only            : mype,master
#endif
   use module_file_write, only : file_write_bin
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include "abort.h"
!
   character(len=128)   ::  fn
   character(len=128)   ::  fnmskg
   real, parameter      ::  undef=1.e30
   integer              ::  idim,jdim
   real                 ::  slmask(idim*jdim)
   real                 ::  fh
   real                 ::  out(idim*jdim,*)
   real                 ::  proj,orient,truth,cotru
   real                 ::  delx,dely,rlat1,rlat2,rlon1,rlon2
   real                 ::  dlon,dlat,wlon,rnlat
   real                 ::  dummy
   real, allocatable    ::  data(:),rslmsk(:)
   integer              ::  lugb
   integer              ::  iy,im,id,ih
   integer              ::  indx
   integer              ::  k,kgau,ij,kpds5,kgds1,ijmdim,imaxgrb,jmaxgrb
   integer              ::  iret
   logical              ::  lmask
!
#ifdef DBG
   print *,'lugb,nlev,indx,fn=',lugb,kgrb(indx),indx,fn
   print *,'idim,jdim=',idim,jdim,kgrb(indx)
   print *,'iy,im,id,ih,fh=',iy,im,id,ih,fh
   print *,'indx=',indx
#endif
   if(fn(1:4).eq.'    ') then
     print *,gvar(indx),' file name empty.'
     return
   endif
!
   ijmdim=imdata(indx)*jmdata(indx)
   allocate (data(ijmdim))
!
   lmask=mask(indx).gt.0
!
!  read grib file
!
!  convert land mask to ocean mask if mask=2
!
   if(mask(indx).eq.2) then
     do ij = 1,idim*jdim
       slmask(ij)=1.-slmask(ij)
     enddo
   endif
!
   do k = 1,kgrb(indx)
#ifdef MP
     if(mype.eq.master) then
#endif
       write(6,*)'========== start reading ',                                  &
              gvar(indx),kgrb(indx),' ========='
       call flush(6)

       call sfc_read_grib(lugb,fn,kpd(1,k,indx),                               &
                 iy,im,id,ih,fh,                                               &
                 imdata(indx),jmdata(indx),                                    &
                 imaxgrb,jmaxgrb,                                              &
                 kpds5,kgds1,                                                  &
                 proj,orient,truth,cotru,                                      &
                 delx,dely,rlat1,rlat2,rlon1,rlon2,                            &
                 dlon,dlat,wlon,rnlat,                                         &
                 data,iret)
       if(iret.gt.0) then
         print *,gvar(indx),' read failed at k=',k
         call MPABORT
       endif
#ifdef DBG
       print *,'imaxgrb,jmaxgrb=',imaxgrb,jmaxgrb
       print *,'kpds5,kgds1=',kpds5,kgds1
       print *,'proj,orient,truth,cotru=',proj,orient,truth,cotru
       print *,'delx,dely,rlat1,rlat2,rlon1,rlon2=',                           &
              delx,dely,rlat1,rlat2,rlon1,rlon2
       print *,'dlon,dlat,wlon,rnlat=',dlon,dlat,wlon,rnlat
       print *,'data=',(data(ij),ij=1,10)
       print *,'lmask,limask(indx),mask(indx)=',                               &
              lmask,limask(indx),mask(indx)
#endif
       close(lugb)
#ifdef MP
     endif
#endif
!
!  define mask on input grid 
!
!  create slmask over input grid using ready made hires lat/lon mask
!  when no mask info is available from grib file itself
!
     allocate(rslmsk(imaxgrb*jmaxgrb))
     !
     if(lmask) then
       if(.not.limask(indx)) then
         if(k.eq.1) then
           call sfc_interp_ll2mask(lugb,fnmskg,                                &
                           2,rslmsk,imaxgrb,jmaxgrb,                           &
                           proj,orient,truth,cotru,                            &
                           delx,dely,rlat1,rlat2,rlon1,rlon2)
!
!  convert land mask to ocean mask if mask=2
!
           if(mask(indx).eq.2) then
             do ij = 1,imaxgrb*jmaxgrb
               rslmsk(ij)=1.-rslmsk(ij)
             enddo
           endif
         endif
       else
         call sfc_field2mask(data,imaxgrb*jmaxgrb,                             &
                         cvalin(indx),condin(indx),rslmsk)
!  convert land mask to ocean mask if mask=2
         if(mask(indx).eq.2) then
           do ij = 1,imaxgrb*jmaxgrb
             rslmsk(ij)=1.-rslmsk(ij)
           enddo
         endif
       endif
     else
       do ij = 1,imaxgrb*jmaxgrb
         rslmsk(ij)=1.
       enddo
       if(mask(indx).eq.2) then
         do ij = 1,imaxgrb*jmaxgrb
           rslmsk(ij)=1.-rslmsk(ij)
         enddo
       endif
     endif
!
!  spacial interpolation 
!
#ifdef RMP
     call sfc_interp_ll2rmp(data,imaxgrb,jmaxgrb,                              &
                  inttyp(indx),out(1,k),idim,jdim,                             &
                  lmask,rslmsk,slmask,                                         &
                  proj,orient,truth,cotru,                                     &
                  delx,dely,rlat1,rlat2,rlon1,rlon2)
#else
     kgau=0
     if(kgds1.eq.4) kgau=1
     call sfc_interp_file_gmp(data,imaxgrb,jmaxgrb,                            &
                 abs(dlon),abs(dlat),wlon,rnlat,                               &
                 out(1,k),idim,jdim,lmask,rslmsk,slmask,kgau,                  &
                 inttyp(indx))
#endif
#ifdef BIN_DBG
     call file_write_bin(119,out(1,k),idim,jdim,1,0)
#endif
     deallocate (rslmsk)
   enddo
!
   deallocate (data)
!
!  get original mask value if mask=2
!
   if(mask(indx).eq.2) then
     do ij = 1,idim*jdim
       slmask(ij)=1.-slmask(ij)
     enddo
   endif
!
!  scale the value if necessary
!
   if(scale(indx).ne.1.) then
     do k = 1,kgrb(indx)
       do ij = 1,idim*jdim
         out(ij,k)=out(ij,k)*scale(indx)
       enddo
     enddo
   endif
!
!  process mask output (0 or 1 field)
! 
   if(lomask(indx)) then
     call sfc_field2mask(out,idim*jdim,cvalout(indx),condout(indx),out)
   endif
#ifdef DBG 
   do k = 1,kgrb(indx)
     print *,gvar(indx),' for k=',k
     call sfc_quick_print(out(1,k),idim,jdim)
   enddo
#endif 
!
   return
   end subroutine sfc_grib_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_grib_numchar(char,num)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   character(len=128)   ::  char
   integer              ::  num,n
!
   num = 1
!
   do while(char(num:num).ne.' ')
     num = num+1
   enddo
!
   num = num-1
   do n = num+1,128
     char(n:n) = ' '
   enddo
!
   return
   end subroutine sfc_grib_numchar
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_read_grib(lugb,fngrib,kkpds,                                 &
                    iy,im,id,ih,fh,imdata,jmdata,                              &
                    imax,jmax,kkpds5,kkgds1,                                   &
                    proji,orienti,truthi,cotrui,                               &
                    delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i,                   &
                    dlon,dlat,wlon,rnlat,                                      &
                    data,iret)
#include "abort.h"
!-------------------------------------------------------------------------------
!
! read in grib climatology files.
!
! interpolate climatology to the dates
!
! grib file should allow all the necessary parameters to be extracted fr
! the description records.
!
!
! nvalid:  analysis later than (current date - nvalid) is regarded as
!          valid for current analysis
!
!-------------------------------------------------------------------------------
   integer, parameter      ::  nvalid=5
   integer, parameter      ::  mbuf=1024*128*128
   character(len=128)      ::  fngrib
   character(len=1)        ::  cbuf(mbuf)
   character(len=80)       ::  asgnstr
!   integer, parameter      ::  mbuf=1024*128*64
!
! iret=0  ... successfully read
!     =-1 ... record with unmatched dates exist
!
   real                    ::  data(imdata*jmdata)
   real, allocatable       ::  data_next(:)
   integer                 ::  kpds(25),kgds(22),kens(5)
   integer                 ::  jpds(25),jgds(22),jens(5)
   integer                 ::  kkpds(25),kpds0(25)
   logical                 ::  lbms(imdata*jmdata)
   data msk1/3840000/,msk2/2400000/
!
#ifdef REAL4_W3LIB
   integer*4                   lugb4,msk14,msk24,mnum4,mbuf4
   integer*4                   nlen4,nnum4,iret4
   integer*4                   ndata4
   real*4, allocatable     ::  data4(:)
   integer*4                   lskip4,lgrib4,lret4
   integer*4                   n4,jpds4(25),jgds4(22),jens4(5)
   integer*4                   k4,kpds4(25),kgds4(22),kens4(5)
#endif
#ifdef FASTBAREAD
   character, allocatable  ::  bbuf(:)
#endif
   logical                 ::  lclim
!
! julian day of the middle of each month
!
   real                    ::  dayhf(13)
   data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,                             &
              196.5,227.5,258.0,288.5,319.0,349.5,380.5/
   save dayhf
!
#ifdef REAL4_W3LIB
   allocate (data4(imdata*jmdata))
#endif
!
   call sfc_interp_date(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
   close(lugb)
   call sfc_grib_numchar(fngrib,nch)
#ifdef ASSIGN
!
   if(lugb.lt.10) then
     write(asgnstr,'(22hassign -s unblocked u:,i1)') lugb
   else
     write(asgnstr,'(22hassign -s unblocked u:,i2)') lugb
   endif
!
   print*,' assign = ',asgnstr
   call assign('assign -R')
   call assign(asgnstr)
   open(unit=lugb,file=fngrib(1:nch),status='old',form='unformatted',          &
        err=910)
#else
#ifdef FASTBAREAD
   call nainit(lugb,fngrib(1:nch),isize,iret)
   allocate (bbuf(isize))
   call naopen(lugb,fngrib(1:nch),bbuf,isize,iret)
#else
#ifdef REAL4_W3LIB
   lugb4=lugb
   call baopen(lugb4,fngrib(1:nch),iret4)
   iret=iret4
#else
   call baopen(lugb,fngrib(1:nch),iret)
#endif
#endif
   if(iret.ne.0) go to 910
#endif
   !
   go to 911
   !
910 continue
   print *,'error in opening file ',fngrib(1:nch)
   call MPABORT
!
911 continue
   write(6,*) ' file ',fngrib(1:nch),' opened. unit=',lugb
!
! get grib index buffer
!
   mnum=0
#ifdef REAL4_W3LIB
   lugb4=lugb
   msk14=msk1
   msk24=msk2
   mnum4=mnum
   mbuf4=mbuf
#ifdef FASTBAREAD
   call ngetgir(bbuf,isize,msk14,msk24,mnum4,mbuf4,                            &
                cbuf,nlen4,nnum4,iret4)
#else
   call getgir(lugb4,msk14,msk24,mnum4,mbuf4,                                  &
               cbuf,nlen4,nnum4,iret4)
#endif
   nlen=nlen4
   nnum=nnum4
   iret=iret4
#else
#ifdef FASTBAREAD
   call ngetgir(bbuf,isize,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#else
   call getgir(lugb,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#endif
#endif
#ifdef DBG
   print *,'nlen=',nlen,' nnum=',nnum
#endif
!
   if(iret.ne.0) then
     print *,'error.  cbuf length too short in ngetgir'
     call MPABORT
   endif
!
   if(nnum.eq.0) then
     print *,'error.  not a grib file. detected in ngetgir'
     call MPABORT
   endif
!
   if(nlen.eq.0) then
     print *,'error.  nlen=0.  detected in ngetgir'
     call MPABORT
   endif
!
   if(nlen.gt.imdata*jmdata) then
     print *,'error.  nlen .gt. imdata*jmdata  detected in ngetgir'
     call MPABORT
   endif
!
! find file type climatology or analysis
!
   do i = 1,25
     jpds(i)=-1
   enddo
!
   do i = 1,22
     jgds(i)=-1
   enddo
!
   do i = 1,5
     jens(i)=-1
   enddo
!
   do i = 1,25
     jpds(i)=kkpds(i)
   enddo
!
! fix for ecmwf grib
!
   if ( (jpds(5).eq.11.or.jpds(5).eq.91) .and.                                 &
        (jpds(7).eq.-1.or.jpds(7).eq.0) ) then
     jpds(5)=-1
     jpds(7)=0
     write(6,*) 'Changing jpds5 and jpds7 to allow ',                          &
                 'for ECMWF and NCEP ice and sst'
   endif
!
#ifdef DBG
   write(6,*) ' Searching the following fields with kpds=.'
   write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
   write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
   write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
!
   n=0
#ifdef REAL4_W3LIB
   nlen4=nlen
   nnum4=nnum
   n4=n
!
   do i = 1,25
     jpds4(i)=jpds(i)
   enddo
!
   do i = 1,22
     jgds4(i)=jgds(i)
   enddo
!
   do i = 1,5
     jens4(i)=jens(i)
   enddo
!
   call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,                         &
                k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
!
   k=k4
   do i = 1,25
     kpds(i)=kpds4(i)
   enddo
!
   do i = 1,22
     kgds(i)=kgds4(i)
   enddo
!
   do i = 1,5
     kens(i)=kens4(i)
   enddo
!
   lskip=lskip4
   lgrib=lgrib4
   iret=iret4
#else
   call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,                               &
                k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
!
#ifdef DBG
   write(6,*) 'first grib record in the  grib file',fngrib(1:nch)
   write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
   write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
   write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
!
   if(lgrib.eq.0) then
     write(6,*) ' error in getgbss.  No matching records.'
     call MPABORT
   endif
!
   do i = 1,25
     kpds0(i)=kpds(i)
   enddo
!
   kpds0(4)=-1
   kpds0(18)=-1
!
! lu_rev: manipulate kpds(13) and kpds(15) for fixed fields
!
! neet to clean up  mk      
!
   if(kpds(5).eq.236 .or. kpds(5).eq.255) then
     kpds(13) = 4
     kpds(15) = 1
     kpds(16) = 51
     print *,' quick fix for noah fixed fields'
   endif
!
   if(kpds(16).eq.51) then
     write(6,*) ' climatology file.'
     lclim=.true.
   else
     write(6,*) ' analysis file.'
     lclim=.false.
   endif
!
! handling climatology file
!
   if(lclim) then
!
! find average type
! weekly,biweekly,monthly,seasonal,annual
!
!  kpds(13)=4 & kpds(15)=1 .. annual mean
!  kpds(13)=3 & kpds(15)=1 .. monthly mean
!  kpds(13)=2 & kpds(15)=7 .. weekly mean
!  kpds(13)=2 & kpds(15)=14.. bi-weekly mean
!
     if(kpds(13).eq.2.and.kpds(15).eq.7) then
       write(6,*) ' this is weekly mean climatology'
       write(6,*) ' cannot process.'
       call MPABORT
     elseif(kpds(13).eq.2.and.kpds(15).eq.14) then
       write(6,*) ' this is bi-weekly mean climatology'
       write(6,*) ' cannot process.'
       call MPABORT
     elseif(kpds(13).eq.3.and.kpds(15).le.1) then
       write(6,*) ' this is monthly mean climatology'
       monend=12
       do mm = 1,monend
         mmm=mm
         mmp=mm+1
         if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
           mon1=mmm
           mon2=mmp
           !
           go to 20
           !
         endif
       enddo
       print *,'wrong rjday',rjday
       call MPABORT
20     continue
       ijmax=0
       do nn = 1,2
         if(nn.eq.2) then
           allocate (data_next(ijmax))
         endif
         lskip = -1
         n=0
         do i = 1,25
           jpds(i)=kpds0(i)
         enddo
         jpds(24)=-1
         jpds(25)=-1
         if(nn.eq.1) jpds( 9)=mon1
         if(nn.eq.2) jpds( 9)=mon2
         if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
         write(6,*) ' Searching the following fields with kpds=.'
         write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
         write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
         write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
         nlen4=nlen
         nnum4=nnum
         n4=n
         do i = 1,25
           jpds4(i)=jpds(i)
         enddo
         do i = 1,22
           jgds4(i)=jgds(i)
         enddo
         do i = 1,5
           jens4(i)=jens(i)
         enddo
         call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,                   &
                         k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
         k=k4
         do i = 1,25
           kpds(i)=kpds4(i)
         enddo
         do i = 1,22
           kgds(i)=kgds4(i)
         enddo
         do i = 1,5
           kens(i)=kens4(i)
         enddo
         lskip=lskip4
         lgrib=lgrib4
         iret=iret4
#else
         call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,                         &
                         k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
!
         if(lgrib.eq.0) then
#ifdef DBG
           write(6,*) ' No matching record found.'
           write(6,*) ' The last header records read:'
           write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
           write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
           write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
           call MPABORT
         endif
!
         write(6,*) ' Matching record found for mon=',jpds(9)
#ifdef DBG
         write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
         write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
         write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
#endif
         do i = 1,imdata*jmdata
           lbms(i)=.true.
         enddo
#ifdef REAL4_W3LIB
         lugb4=lugb
#ifdef FASTBAREAD
         call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,          &
                       data4,6)
#else
         call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,                &
                      data4,6)
#endif
         do i = 1,25
           kpds(i)=kpds4(i)
         enddo
         do i = 1,22
           kgds(i)=kgds4(i)
         enddo
         ndata=ndata4
         if(nn.eq.1) then
           do i = 1,ndata
             data(i)=data4(i)
           enddo
         else
           do i = 1,ndata
             data_next(i)=data4(i)
           enddo
         endif
#else
         if(nn.eq.1) then
#ifdef FASTBAREAD
           call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,             &
                          data,6)
#else
           call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,                    &
                         data,6)
#endif
         else
#ifdef FASTBAREAD
           call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,             &
                          data_next,6)
#else
           call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,                    &
                         data_next,6)
#endif
         endif
#endif
         if(ndata.eq.0) then
           write(6,*) ' error in nrdgb'
           call MPABORT
         endif
         imax=kgds(2)
         jmax=kgds(3)
         ijmax=imax*jmax
#ifdef DBG
         write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
       enddo
!
       wei1=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
       wei2=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
       if(mon2.eq.13) mon2=1
#ifdef DBG
       print *,'rjday,mon1,mon2,wei1,wei2=',                                   &
                  rjday,mon1,mon2,wei1,wei2
#endif
       do i = 1,ijmax
         data(i)=wei1*data(i)+wei2*data_next(i)
       enddo
       deallocate (data_next)
     elseif(kpds(13).eq.5.and.kpds(15).le.1) then
       write(6,*) ' this is monthly mean climatology (no-weighting)'
       monend=12
       do mm = 1,monend
         mmm=mm
         mmp=mm+1
         if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
           mon1=mmm
           mon2=mmp
           !
           go to 70
           !
         endif
       enddo
       print *,'wrong rjday',rjday
       call MPABORT
70     continue
       ijmax=0
       do nn = 1,2
         if(nn.eq.2) then
           allocate (data_next(ijmax))
         endif
         lskip = -1
         n=0
         do i = 1,25
           jpds(i)=kpds0(i)
         enddo
         jpds(24)=-1
         jpds(25)=-1
         if(nn.eq.1) jpds( 9)=mon1
         if(nn.eq.2) jpds( 9)=mon2
         if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
         write(6,*) ' Searching the following fields with kpds=.'
         write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
         write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
         write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
         nlen4=nlen
         nnum4=nnum
         n4=n
         do i = 1,25
           jpds4(i)=jpds(i)
         enddo
         do i = 1,22
           jgds4(i)=jgds(i)
         enddo
         do i = 1,5
           jens4(i)=jens(i)
         enddo
         call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,                   &
                         k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
         k=k4
         do i = 1,25
           kpds(i)=kpds4(i)
         enddo
         do i = 1,22
           kgds(i)=kgds4(i)
         enddo
         do i = 1,5
           kens(i)=kens4(i)
         enddo
         lskip=lskip4
         lgrib=lgrib4
         iret=iret4
#else
         call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,                         &
                         k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
!
         if(lgrib.eq.0) then
#ifdef DBG
           write(6,*) ' No matching record found.'
           write(6,*) ' The last header records read:'
           write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
           write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
           write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
           call MPABORT
         endif
!
         write(6,*) ' Matching record found for mon=',jpds(9)
#ifdef DBG
         write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
         write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
         write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
#endif
         do i = 1,imdata*jmdata
           lbms(i)=.true.
         enddo
#ifdef REAL4_W3LIB
         lugb4=lugb
#ifdef FASTBAREAD
         call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,          &
                       data4,6)
#else
         call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,                &
                      data4,6)
#endif
         do i = 1,25
           kpds(i)=kpds4(i)
         enddo
         do i = 1,22
           kgds(i)=kgds4(i)
         enddo
         ndata=ndata4
         if(nn.eq.1) then
           do i = 1,ndata
             data(i)=data4(i)
           enddo
         else
           do i = 1,ndata
             data_next(i)=data4(i)
           enddo
         endif
#else
         if(nn.eq.1) then
#ifdef FASTBAREAD
           call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,             &
                          data,6)
#else
           call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,                    &
                         data,6)
#endif
         else
#ifdef FASTBAREAD
           call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,             &
                          data_next,6)
#else
           call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,                    &
                         data_next,6)
#endif
         endif
#endif
         if(ndata.eq.0) then
           write(6,*) ' error in nrdgb'
           call MPABORT
         endif
         imax=kgds(2)
         jmax=kgds(3)
         ijmax=imax*jmax
#ifdef DBG
         write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
       enddo
!
       if(mon2.eq.13) mon2=1
#ifdef DBG
       print *,'rjday,mon1,mon2=',                                             &
                  rjday,mon1,mon2
#endif
       do i = 1,ijmax
         data(i)=data_next(i)
       enddo
       deallocate (data_next)
     elseif(kpds(13).eq.4.and.kpds(15).eq.3) then
       write(6,*) ' this is seasonal mean climatology'
       monend=4
       do mm = 1,monend
         mmm=mm*3-2
         mmp=(mm+1)*3-2
         if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
           mon1=mmm
           mon2=mmp
           !
           go to 30
           !
         endif
       enddo
       print *,'wrong rjday',rjday
       call MPABORT
30     continue
       is=im/3+1
       if(is.eq.5) is=1
       is1=mon1/3+1
       if(is1.eq.5) is1=1
       is2=mon2/3+1
       if(is2.eq.5) is2=1
       !
       ijmax=0
       do nn = 1,2
         if(nn.eq.2) then
           allocate (data_next(ijmax))
         endif
         do i = 1,25
           jpds(i)=kpds0(i)
         enddo
         jpds(24)=-1
         jpds(25)=-1
         n=0
         if(nn.eq.1) then
           isx=is1
         else
           isx=is2
         endif
         if(isx.eq.1) jpds(9)=12
         if(isx.eq.2) jpds(9)=3
         if(isx.eq.3) jpds(9)=6
         if(isx.eq.4) jpds(9)=9
         if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
         write(6,*) ' Searching the following fields with kpds=.'
         write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
         write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
         write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
         nlen4=nlen
         nnum4=nnum
         n4=n
         do i = 1,25
           jpds4(i)=jpds(i)
         enddo
         do i = 1,22
           jgds4(i)=jgds(i)
         enddo
         do i = 1,5
           jens4(i)=jens(i)
         enddo
         call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,                   &
                       k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
         k=k4
         do i = 1,25
           kpds(i)=kpds4(i)
         enddo
         do i = 1,22
           kgds(i)=kgds4(i)
         enddo
         do i = 1,5
           kens(i)=kens4(i)
         enddo
         lskip=lskip4
         lgrib=lgrib4
         iret=iret4
#else
         call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,                         &
                         k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
!
         if(lgrib.eq.0) then
#ifdef DBG
           write(6,*) ' No matching record found.'
           write(6,*) ' The last header records read:'
           write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
           write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
           write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
           call MPABORT
         endif
         do i = 1,imdata*jmdata
           lbms(i)=.true.
         enddo
#ifdef REAL4_W3LIB
         lugb4=lugb
#ifdef FASTBAREAD
         call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,          &
                     data4,6)
#else
         call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,                &
                    data4,6)
#endif
         do i = 1,25
           kpds(i)=kpds4(i)
         enddo
         do i = 1,22
           kgds(i)=kgds4(i)
         enddo
         ndata=ndata4
         if(nn.eq.1) then
           do i = 1,ndata
             data(i)=data4(i)
           enddo
         else
           do i = 1,ndata
             data_next(i)=data4(i)
           enddo
         endif
#else
         if(nn.eq.1) then
#ifdef FASTBAREAD
           call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,             &
                        data,6)
#else
           call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,                    &
                       data,6)
#endif
         else
#ifdef FASTBAREAD
           call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,             &
                        data_next,6)
#else
           call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,                    &
                       data_next,6)
#endif
         endif
#endif
         if(ndata.eq.0) then
           write(6,*) ' error in nrdgb'
           call MPABORT
         endif
         imax=kgds(2)
         jmax=kgds(3)
         ijmax=imax*jmax
       enddo
       wei1=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
       wei2=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
       if(mon2.eq.13) mon2=1
#ifdef DBG
       print *,'rjday=',rjday
       print *,'mon1 =',mon1 ,' mon2=',mon2
       print *,'wei1 =',wei1 ,' wei2=',wei2
       print *,'ses1 =', is1 ,' ses2=', is2
#endif
       do i = 1,ijmax
         data(i)=wei1*data(i)+wei2*data_next(i)
       enddo
       deallocate (data_next)
     elseif(kpds(13).eq.4.and.kpds(15).eq.1) then
       write(6,*) ' this is annual mean climatology'
       monend=-1
#ifdef REAL4_W3LIB
       lugb4=lugb
#ifdef FASTBAREAD
       call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,            &
                  data4,6)
#else
       call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,                  &
                 data4,6)
#endif
       do i = 1,25
         kpds(i)=kpds4(i)
       enddo
       do i = 1,22
         kgds(i)=kgds4(i)
       enddo
       ndata=ndata4
       do i = 1,ndata
         data(i)=data4(i)
       enddo
#else
#ifdef FASTBAREAD
       call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#else
       call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#endif
#endif
       if(ndata.eq.0) then
         write(6,*) ' error in nrdgb'
         write(6,*) ' kpds=',kpds
         write(6,*) ' kgds=',kgds
         call MPABORT
       endif
       imax=kgds(2)
       jmax=kgds(3)
       ijmax=imax*jmax
     else
       write(6,*) ' climatology file average period unknown.'
       write(6,*) ' kpds(13)=',kpds(13),' kpds(15)=',kpds(15)
       call MPABORT
     endif
   else
!-------------------------------------------------------------------
!
!  handling analysis file
!
!  find record for the given hour/day/month/year
!
     monend=0
     nrept=0
     do i = 1,25
       kpds(i)=kpds0(i)
     enddo
     iyr=jy
     imo=jm
     idy=jd
     ihr=jh
50   continue
     jpds( 8)=iyr
     jpds( 9)=imo
     jpds(10)=idy
     if(ihr.eq.-1) then
! t.r. jpds(11)=0
       jpds(11)=-1
     else
! t.r. jpds(11)=ihr
       jpds(11)=-1
     endif
#ifdef DBG
     write(6,*) ' Searching the following fields with kpds=.'
     write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
     write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
     write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
     n=0
#ifdef REAL4_W3LIB
     nlen4=nlen
     nnum4=nnum
     n4=n
     do i = 1,25
       jpds4(i)=jpds(i)
     enddo
     do i = 1,22
       jgds4(i)=jgds(i)
     enddo
     do i = 1,5
       jens4(i)=jens(i)
     enddo
     call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,                       &
                 k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
     k=k4
     do i = 1,25
       kpds(i)=kpds4(i)
     enddo
     do i = 1,22
       kgds(i)=kgds4(i)
     enddo
     do i = 1,5
       kens(i)=kens4(i)
     enddo
     lskip=lskip4
     lgrib=lgrib4
     iret=iret4
#else
     call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,                             &
                 k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
!
     if(lgrib.ne.0) then
       if(nrept.le.iabs(nvalid).and.nrept.ne.0) then
         write(6,*) '<warning:cycl> grib record date does not match',          &
                     ' but within ',iabs(nvalid),' days.'
         write(6,*) '<warning:cycl> assume it is valid.'
       endif
       write(6,*) ' Matching record found'
#ifdef DBG
       write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
       write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
       write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
       write(6,*) ' Now reading grib record'
#endif
#ifdef REAL4_W3LIB
       lugb4=lugb
#ifdef FASTBAREAD
       call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,            &
                  data4,6)
#else
       call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,data4,6)
#endif
       do i = 1,25
         kpds(i)=kpds4(i)
       enddo
       do i = 1,22
         kgds(i)=kgds4(i)
       enddo
       ndata=ndata4
       do i = 1,ndata
         data(i)=data4(i)
       enddo
#else
#ifdef FASTBAREAD
       call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#else
       call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#endif
#endif
       if(ndata.eq.0) then
         write(6,*) ' error in nrdgb'
         write(6,*) ' kpds=',kpds
         write(6,*) ' kgds=',kgds
         write(6,*) ' lgrib,lskip=',lgrib,lskip
         call MPABORT
       endif
#ifdef DBG
       write(6,*) ' matching grib record successfully read.'
#endif
       imax=kgds(2)
       jmax=kgds(3)
       ijmax=imax*jmax
     else
       if(nrept.eq.0) then
         write(6,*) ' no matching dates found.  start searching',              &
                     ' the nearest matching past date.'
       endif
!
!  no matching ih found. search nearest hour
!
       if(ihr.eq.6) then
         ihr=0
         write(6,*) ' hour set to 0'
         !
         go to 50
         !
       elseif(ihr.eq.12) then
         ihr=0
         write(6,*) ' hour set to 0'
         !
         go to 50
         !
       elseif(ihr.eq.18) then
         ihr=12
         write(6,*) ' hour set to 12'
         !
         go to 50
         !
       elseif(ihr.eq.0.or.ihr.eq.-1) then
         idy=idy-1
         if(idy.eq.0) then
           imo=imo-1
           if(imo.eq.0) then
             iyr=iyr-1
             if(iyr.lt.0) iyr=99
             imo=12
           endif
           idy=31
           if(imo.eq.4.or.imo.eq.6.or.imo.eq.9.or.imo.eq.11) idy=30
           if(imo.eq.2) then
             if(mod(iyr,4).eq.0) then
               idy=29
             else
               idy=28
             endif
           endif
         endif
         ihr=-1
         write(6,*) ' dates decremented to:',iyr,imo,idy,ihr
         nrept=nrept+1
         if(nrept.gt.nvalid) then
           iret=-1
           write(6,*) ' <warning:surface> search range exceeded.'
           write(6,*) ' <warning:surface> terminating search.'
           call MPABORT
         endif
         !
         go to 50
         !
       else
#ifdef DBG
         write(6,*) ' Matching date could not be fond.'
         write(6,*) ' The last header records read:'
         write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
         write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
         write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
         call MPABORT
       endif
     endif
   endif
!
!  finished reading grib file. 
!
!  substitution to the output array
!
80 continue
   call print_maxmin_six(data,imax*jmax,1,1,1,'input data')
!
#ifdef DBG
   if(kkpds(5).eq.225) then
     print *,'input uninterpolated vegtype'
     call sfc_quick_print(data,imax,jmax)
   endif
   if(kkpds(5).eq.87) then
     print *,'input uninterpolated vegfrac'
     call sfc_quick_print(data,imax,jmax)
   endif
   if(kkpds(5).eq.230) then
     print *,'input uninterpolated soiltype'
     call sfc_quick_print(data,imax,jmax)
   endif
#endif
!
#ifdef DBG
   write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
   proji=0.
   orienti=0.
   truthi=0.
   cotrui=0.
   call sfc_get_area_rmp(kgds,proji,orienti,truthi,cotrui,                     &
                 delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i)
#ifdef DBG
   write(6,*) 'proji,orienti,truthi,cotrui='
   write(6,*)  proji,orienti,truthi,cotrui
   write(6,*) 'delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i='
   write(6,*)  delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i
#endif
#ifndef RMP
   call sfc_get_area_gmp(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr)
#ifdef DBG
   write(6,*) 'dlat,dlon,rslat,rnlat,wlon,elon,ijordr ='
   write(6,*)  dlat,dlon,rslat,rnlat,wlon,elon,ijordr
#endif
   call sfc_file_subst(data,imax,jmax,dlon,dlat,ijordr)
#endif
!
   write(6,*) 'imax,jmax of grib=',imax,jmax
!
   kkpds5=kkpds(5)
   kkgds1=kgds(1)
!
#ifdef REAL4_W3LIB
   deallocate (data4)
#endif
#ifdef FASTBAREAD
   deallocate (bbuf)
#endif
   close(lugb)
!
   return
   end subroutine sfc_read_grib
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_read_mask(lugb,fnmask,ijdim,slmask)
!-------------------------------------------------------------------------------
!
!  land/sea mask
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   character(len=132)      ::  fnmask
   integer                 ::  lugb,ijdim,nch
   real                    ::  slmask(ijdim)
   integer                 ::  ij
!
   if(fnmask(1:6).eq.'      ') then
     print *,'land-sea mask file name empty'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
   endif
!
   close(lugb)
#ifdef ASSIGN
   call assign('assign -R')
#endif
   call sfc_grib_numchar(fnmask,nch)
   open(unit=lugb,file=fnmask(1:nch),status='old',form='unformatted',          &
        err=920)
   !
   go to 921
   !
920 continue
   print *,'error in opening file ',fnmask(1:nch)
#ifdef MP
#ifdef RMP
   call rmpabort
#else
   call mpabort
#endif
#else
   call abort
#endif
921 continue
   write(6,*) ' file ',fnmask(1:nch),' opened. unit=',lugb
   read(lugb) slmask
!
   do ij = 1,ijdim
     slmask(ij)=nint(slmask(ij))
   enddo
!
   close(lugb)
!
   return
   end subroutine sfc_read_mask
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_read_mtn(lugb,fnorog,ijdim,slmask,orog)
!-------------------------------------------------------------------------------
!
!  read model orography
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include <abort.h>
!
   character(len=132)   ::  fnorog
   integer              ::  lugb,ijdim
   real                 ::  slmask(ijdim)
   real                 ::  orog(ijdim)
!
   integer, allocatable ::  loismsk(:)
   integer              ::  ij,nch
!
!  max/min over ocean, no snow land, no-snow sea ice, snow covered land
!  and snow covered sea ice.  at this stage, slmask contain only 0 or 1.
!
   real                 ::  vmaxmin(2,5)
   data vmaxmin/9000.,-5000.,                                                  &
                9000.,-5000.,                                                  &
                   0.,    0.,                                                  &
                   0.,    0.,                                                  &
                   0.,    0.                                                   &
                /
!
   if(fnorog(1:6).eq.'      ') then
     print *,'orography file name empty'
     call MPABORT
   endif
!
   close(lugb)
#ifdef ASSIGN
   call assign('assign -R')
#endif
   call sfc_grib_numchar(fnorog,nch)
   open(unit=lugb,file=fnorog(1:nch),status='old',form='unformatted',          &
        err=910)
   !
   go to 911
   !
910 continue
   print *,'error in opening file ',fnorog(1:nch)
   call MPABORT
911 continue
   write(6,*) ' file ',fnorog(1:nch),' opened. unit=',lugb
   read(lugb) orog
!
   allocate (loismsk(ijdim))
!
   do ij = 1,ijdim
     loismsk(ij)=nint(slmask(ij))
   enddo
!
!  quality control by checking max/min over various surfaces
!
!hoon call sfc_check_maxmin(orog,ijdim,loismsk,vmaxmin,'orog    ')
!
   deallocate (loismsk)
   close(lugb)
!
   return
   end subroutine sfc_read_mtn
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_field2mask(data,ijdim,cval,cond,out)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ijdim,ij
   real                 ::  cval
   real                 ::  data(ijdim),out(ijdim)
   character(len=2)     ::  cond
!
   if(cond.eq.'ge') then
     do ij = 1,ijdim
       if(data(ij).ge.cval) then
         out(ij)=1.
       else
         out(ij)=0.
       endif
     enddo
   elseif(cond.eq.'gt') then
     do ij = 1,ijdim
       if(data(ij).gt.cval) then
         out(ij)=1.
       else
         out(ij)=0.
       endif
     enddo
   elseif(cond.eq.'le') then
     do ij = 1,ijdim
       if(data(ij).le.cval) then
         out(ij)=1.
       else
         out(ij)=0.
       endif
     enddo
   elseif(cond.eq.'lt') then
     do ij = 1,ijdim
       if(data(ij).lt.cval) then
         out(ij)=1.
       else
         out(ij)=0.
       endif
     enddo
   elseif(cond.eq.'eq') then
     do ij = 1,ijdim
       if(data(ij).eq.cval) then
         out(ij)=1.
       else
         out(ij)=0.
       endif
     enddo
   elseif(cond.eq.'ne') then
     do ij = 1,ijdim
       if(data(ij).ne.cval) then
         out(ij)=1.
       else
         out(ij)=0.
       endif
     enddo
   else
     write(6,*) ' illegal cond in sfc_field2mask.  cond=',cond
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
   endif
!
   return
   end subroutine sfc_field2mask
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_get_area_gmp(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr)
!-------------------------------------------------------------------------------
!
!  get area of the grib record
!
!-------------------------------------------------------------------------------
   integer, dimension(22)  ::  kgds
   logical                 ::  ijordr
!
#ifdef DBG
   write(6,*) ' kgds( 1-10)=',(kgds(j),j= 1,10)
   write(6,*) ' kgds(11-20)=',(kgds(j),j=11,20)
   write(6,*) ' kgds(21-  )=',(kgds(j),j=21,22)
#endif
!
   if(kgds(1).eq.0) then
!
!  lat/lon grid
!
     write(6,*) 'lat/lon grid'
     dlat=float(kgds(10))/1000.0
     dlon=float(kgds( 9))/1000.0
     f0lon=float(kgds(5))/1000.0
     f0lat=float(kgds(4))/1000.0
     kgds11=kgds(11)
!
!  increase accuracy of dlon if possible (11-24-00)
!
     dlonx=360./float(kgds(2))
     if(abs(dlonx-dlon).lt.0.001) then
       if(abs(dlon).eq.abs(dlat)) then
         dlat=dlat/abs(dlat)*dlonx
       endif
       dlon=dlon/abs(dlon)*dlonx
     endif
     if(kgds11.ge.128) then
       wlon=f0lon-dlon*(kgds(2)-1)
       elon=f0lon
       if(dlon*kgds(2).gt.359.99) then
         wlon=f0lon-dlon*kgds(2)
       endif
       dlon=-dlon
       kgds11=kgds11-128
     else
       wlon=f0lon
       elon=f0lon+dlon*(kgds(2)-1)
       if(dlon*kgds(2).gt.359.99) then
         elon=f0lon+dlon*kgds(2)
       endif
     endif
     if(kgds11.ge.64) then
       rnlat=f0lat+dlat*(kgds(3)-1)
       rslat=f0lat
       kgds11=kgds11-64
     else
       rnlat=f0lat
       rslat=f0lat-dlat*(kgds(3)-1)
       dlat=-dlat
     endif
     if(kgds11.ge.32) then
       ijordr=.false.
     else
       ijordr=.true.
     endif
!
     if(wlon.gt.180.) wlon=wlon-360.
     if(elon.gt.180.) elon=elon-360.
     wlon=nint(wlon*1000.)/1000.
     elon=nint(elon*1000.)/1000.
     rslat=nint(rslat*1000.)/1000.
     rnlat=nint(rnlat*1000.)/1000.
     return
!
!  mercator projection
!
   elseif(kgds(1).eq.1) then
     write(6,*) 'mercator grid'
     write(6,*) 'cannot process'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
!
!  gnomonic projection
!
   elseif(kgds(1).eq.2) then
     write(6,*) 'gnomonic grid'
     write(6,*) 'error!! gnomonic projection not coded'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
!
!  lambert conformal
!
   elseif(kgds(1).eq.3) then
     write(6,*) 'lambert conformal'
     write(6,*) 'cannot process'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
   elseif(kgds(1).eq.4) then
!
!  gaussian grid
!
     write(6,*) 'gaussian grid'
     dlat=99.
     dlon=float(kgds( 9))/1000.0
     f0lon=float(kgds(5))/1000.0
     f0lat=99.
     kgds11=kgds(11)
     if(kgds11.ge.128) then
       wlon=f0lon
       elon=f0lon
       if(dlon*kgds(2).gt.359.99) then
         wlon=f0lon-dlon*kgds(2)
       endif
       dlon=-dlon
       kgds11=kgds11-128
     else
       wlon=f0lon
       elon=f0lon+dlon*(kgds(2)-1)
       if(dlon*kgds(2).gt.359.99) then
         elon=f0lon+dlon*kgds(2)
       endif
     endif
     if(kgds11.ge.64) then
       rnlat=99.
       rslat=99.
       kgds11=kgds11-64
     else
       rnlat=99.
       rslat=99.
       dlat=-99.
     endif
     if(kgds11.ge.32) then
       ijordr=.false.
     else
       ijordr=.true.
     endif
     return
!
!  polar strereographic
!
   elseif(kgds(1).eq.5) then
     write(6,*) 'polar stereographic grid'
     write(6,*) 'cannot process'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
     return
!
!  oblique lambert conformal
!
   elseif(kgds(1).eq.13) then
     write(6,*) 'oblique lambert conformal grid'
     write(6,*) 'cannot process'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
!
!  spherical coefficient
!
   elseif(kgds(1).eq.50) then
     write(6,*) 'spherical coefficient'
     write(6,*) 'cannot process'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
     return
!
!  space view perspective (orthographic grid)
!
   elseif(kgds(1).eq.90) then
     write(6,*) 'space view perspective grid'
     write(6,*) 'cannot process'
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
     return
!
!  unknown projection.  abort.
!
   else
     write(6,*) 'error!! unknown map projection'
     write(6,*) 'kgds(1)=',kgds(1)
     print *,'error!! unknown map projection'
     print *,'kgds(1)=',kgds(1)
#ifdef MP
#ifdef RMP
     call rmpabort
#else
     call mpabort
#endif
#else
     call abort
#endif
   endif
!
   return
   end subroutine sfc_get_area_gmp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_get_area_rmp(kgds,proj,orient,truth,cotru,                   &
                       delx,dely,rlat1,rlat2,rlon1,rlon2)
!-------------------------------------------------------------------------------
!
!  get area of the grib record
!
!-------------------------------------------------------------------------------
!fpp$ noconcur r
   integer, dimension(25)    ::  kgds
   logical                   ::  ijordr
!-------------------------------------------------------------------------------
#ifdef DBG
   write(6,*) ' kgds( 1-10)=',(kgds(j),j= 1,10)
   write(6,*) ' kgds(11-20)=',(kgds(j),j=11,20)
   write(6,*) ' kgds(21-  )=',(kgds(j),j=21,25)
!
#endif
   idrt=kgds(1)
   rlat1=kgds(4)*1.e-3
   rlon1=kgds(5)*1.e-3
!
   if(idrt.eq.0) then
     proj=3
     rlat2=kgds(07)*1.e-3  ! latitude of last point
     rlon2=kgds(08)*1.e-3  ! longitude of last point
     delx=kgds(09)  ! dx (milledegree) on truth latitude
     dely=kgds(10)  ! dy (milledegree) on truth latitude
!    orient=kgds(11)  ! scanning mode ( 0 : n->s, 64 : s->n)
     if(kgds(11).eq.0.) proj = -3.
   endif
!
   if( idrt.eq.4 ) then      ! gaussian projection
     proj=4
     rlat2=kgds(07)*1.e-3  ! latitude of last point
     rlon2=kgds(08)*1.e-3  ! longitude of last point
     delx=kgds(09)  ! dx (milledegree) on truth latitude
     dely=kgds(10)  ! dy (milledegree) on truth latitude
!    orient=kgds(11)  ! scanning mode ( 0 : n->s, 64 : s->n)
     if(kgds(11).eq.0.) proj = -4.
   endif
!
   if( idrt.eq.1 ) then      ! mercater projection
     proj=0
     rlat2=kgds(07)*1.e-3  ! latitude of last point
     rlon2=kgds(08)*1.e-3  ! longitude of last point
     delx=kgds(12)  ! dx (meter) on truth latitude
     dely=kgds(13)  ! dy (meter) on truth latitude
     orient=rlon1
     truth=kgds(09)*1.e-3 ! truth latitude
     cotru=truth          ! co-truth latitude
   endif
!
   if( idrt.eq.5 ) then      ! polar projection
     truth=60.0              ! truth latitude
     cotru=60.0              ! co-truth latitude
     orient=kgds(07)*1.e-3   ! orientation
     delx=kgds(08)  ! dx (meter) on 60 deg
     dely=kgds(09)  ! dy (meter) on 60 deg
     iproj=kgds(10) ! polar projection (first bit 0:north;1=south)
     if( iproj.eq.0 ) then
       proj=1.0
     else
       proj=-1.0
     endif
   endif
!
   if( idrt.eq.3 ) then      ! lambert projection
     orient=kgds(07)*1.e-3   ! orientation
     delx=kgds(08)  ! dx (meter) on 60 deg
     dely=kgds(09)  ! dy (meter) on 60 deg
     iproj=kgds(10) !  projection (first bit 0:north;1=south)
     if( iproj.eq.0 ) then
       proj=2.0
     else
       proj=-2.0
     endif
     truth=kgds(12)*1.e-3 ! the 1st lat from pole to cut
     cotru=kgds(13)*1.e-3 ! the 2nd lat from pole to cut
   endif
!
   return
   end subroutine sfc_get_area_rmp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_date(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
#include "abort.h"
!-------------------------------------------------------------------------------
!
! julian day of the middle of each month
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  iy,im,id,ih,jy,jm,jd,jh
   real                 ::  fh,rjday
   real, dimension(13)  ::  dayhf
   data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,                             &
              196.5,227.5,258.0,288.5,319.0,349.5,380.5/
   save dayhf
!
! number of days in a month
!
   integer,dimension(12)    ::  mjday
   data mjday/31,28,31,30,31,30,31,31,30,31,30,31/
   save mjday
!
! julian day of the first day of the month
!
   real, dimension(12)  ::  fjday
   integer              ::  iret,mon,monend,imm,mfmon,incdy,incd,mondy
!
   logical ijordr
!
   iret=0
!
   monend=9999
!
   fjday(1)=1.
   do mon = 2,12
     fjday(mon)=fjday(mon-1)+float(mjday(mon-1))
   enddo
!
!  get julian day of the iy/im/id/ih provided.
!
   rjday=0.
   imm=im-1
   if(imm.ge.1) then
     do mon = 1,imm
       rjday=rjday+mjday(mon)
     enddo
   endif
!
   rjday=rjday+id+float(ih)/24.0
   rjday=rjday+fh/24.0
   rjday = mod(rjday,365.)
   if(rjday.eq.0.) rjday = 365.
!
   if(rjday.le.0..or.rjday.gt.365.) then
     print *,'wrong rjday',rjday
     call MPABORT
   endif
!
   do mon = 1,11
     if(rjday.lt.fjday(mon+1)) then
       mfmon=mon
       !
       go to 10
       !
     endif
   enddo
!
   mfmon=12
10 continue
   if(rjday.lt.dayhf(1)) rjday=rjday+365.
!
!  compute jy,jm,jd,jh of forecast
!
   jy=iy
   jm=im
   jd=id
   incdy=int(fh/24.)
   jh=ih+mod(fh,24.)
   incdy=incdy+jh/24
   jh=mod(jh,24)
   if(incdy.ge.1) then
     do incd = 1,incdy
       jd=jd+1
       if(jm.eq.4.or.jm.eq.6.or.jm.eq.9.or.jm.eq.11) then
         mondy=30
       elseif(jm.eq.2) then
         if((mod(jy,4).eq.0.and.mod(jy,100).ne.0).or.mod(jy,400).eq.0) then
           mondy=29
         else
           mondy=28
         endif
       else
         mondy=31
       endif
       if(jd.gt.mondy) then
         jm=jm+1
         jd=1
         if(jm.gt.12) then
           jy=jy+1
           jm=1
         endif
       endif
     enddo
   endif
!
!  write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
!
   return
   end subroutine sfc_interp_date
!-------------------------------------------------------------------------------
