#include "define.h"
   subroutine chgrp_read_pressuregrib(lugb,fn,iy,im,id,ih,fh,                  &
                                      stdprs,kdimprs,kdimq,nfld,               &
                                      iupd5,ivpd5,itpd5,iqpd5,izpd5,           &
#ifdef RMP
                                      flat,flon,fm2,rcsln,rsnln,               &
#endif
                                      out,idim,jdim)
!-------------------------------------------------------------------------------
!
! subprogram: chgrp_read_pressuregrib
! 
! abstract:   read standard pressure level grib file, interpolates to model grid.
!
!  input:
!
!    lugb    int                   .. unit number for reading grib
!    fn      char                  .. grib file name
!    idim    int                   .. output x-dimension (model x-dimension)
!    jdim    int                   .. output y-dimension (model y-dimension)
!    iy      int                   .. 4-digit year
!    im      int                   .. month
!    id      int                   .. day
!    ih      int                   .. hour
!    fh      real                  .. forecast hour
!    stdprs  real                  .. standard pressure levels to extract
!    kdimprs int                   .. number of standard pressure levels
!    nfld    int                   .. number of dependent variables (normally 5)
!                                                                    1=relative humidity
!
!  output:
!   
!    out    real array (idim*jdim,kdimprs) .. output field
!
!  the map related parameters appear in this program are for grib
!  files and not those of the forecast program
!
!  note on map projection:
!
!   abs(proj)=0  : mercator
!            =1  : polar stereopgraphic
!            =2  : lambert
!            =3  : lat/lon
!            =4  : gaussian
!    negative value indicates origin in southern hemisphere
!
!-------------------------------------------------------------------------------
   use constant, only     : pi_,qmin_
#ifdef RMP
   use rscomloc
   use module_trans, only : rmp_trans2output_grid
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  lugb
   integer              ::  iupd5,ivpd5,itpd5,iqpd5,izpd5
   integer              ::  idim,jdim
   integer              ::  iy,im,id,ih
   integer              ::  kdimprs,kdimq
   integer              ::  nfld
   character(len=128)   ::  fn
   real                 ::  fh
   real                 ::  stdprs(kdimprs)
   real                 ::  out(idim*jdim,kdimprs,*)
!
   real, allocatable    ::  data(:)
#ifdef RMP
   real, allocatable    ::  rrcsln(:),rrsnln(:)
   real                 ::  flat(idim*jdim),flon(idim*jdim),fm2(idim*jdim)
   real                 ::  rcslsn(idim*jdim),rsnln(idim*jdim)
#endif
!
   integer,allocatable  ::  kpd(:,:,:)
   integer              ::  k,kgau,i,ij,j,kpds5,kgds1,ijmdim,imaxgrb,jmaxgrb,ijdim
   integer              ::  iret
   real                 ::  proj,orient,truth,cotru
   real                 ::  delx,dely
   real                 ::  dlon,dlat,wlon,rnlat
   real                 ::  rlat1,rlat2,rlon1,rlon2
   real,allocatable     ::  slmsk1(:,:)
   real                 ::  slmsko(idim,jdim)
!
!  max possible number of grid points on input grib file
!
   integer,parameter    ::  imdata=1440,jmdata=721
!
   integer              ::  n
!-------------------------------------------------------------------------------
#ifdef DBG
   print *,'lugb,fn=',lugb,fn
   print *,'idim,jdim=',idim,jdim
   print *,'iy,im,id,ih,fh=',iy,im,id,ih,fh
#endif
!
   ijdim=idim*jdim
!
   if(fn(1:4).eq.'    ') then
     print *,'pressure grib file name empty.'
     call abort
   endif
!
   allocate (kpd(25,kdimprs,nfld))
#ifdef RMP
   allocate(rrcsln(ijdim),rrsnl(ijdim))
#endif
!
   ijmdim=imdata*jmdata
   allocate (data(ijmdim))
!
!  set kpd array
!
!  n=1; u  (kpds5=33)
!  n=2: v  (kpds5=34)
!  n=3: t  (kpds5=11)
!  n=4: q  (kpds5=51)
!  n=4: rh (kpds5=52)
!  n=5: z  (kpds5= 7)
!
!   kpds6=100 for isobaric level
!   kpds7=isobaric level pressure in hPa
!
   do n = 1,nfld
     do k = 1,kdimprs
       do i = 1,25
         kpd(i,k,n)=-1
       enddo
     enddo
   enddo
!
   do k = 1,kdimprs
     kpd(5,k,1)=iupd5
     kpd(5,k,2)=ivpd5
     kpd(5,k,3)=itpd5
     kpd(5,k,4)=iqpd5
     kpd(5,k,5)=izpd5
     do n = 1,nfld
       kpd(6,k,n)=100
       kpd(7,k,n)=stdprs(k)
     enddo
   enddo
!      kpd(5,2,5)=1
!      kpd(6,2,5)=1
!      kpd(7,2,5)=0
!      kpd(5,3,5)=7
!      kpd(6,3,5)=1
!      kpd(7,3,5)=0
!
!  read grib file
!
   do n = 1,nfld
     do k = 1,kdimprs
       if (k.le.kdimq.or.iqpd5.ne.kpd(5,k,n)) then
         print *,'== start reading ',kpd(5,k,n),kpd(6,k,n),kpd(7,k,n)
         call sfc_read_grib(lugb,fn,kpd(1,k,n),                                &
                  iy,im,id,ih,fh,                                              &
                  imdata,jmdata,                                               &
                  imaxgrb,jmaxgrb,                                             &
                  kpds5,kgds1,                                                 &
                  proj,orient,truth,cotru,                                     &
                  delx,dely,rlat1,rlat2,rlon1,rlon2,                           &
                  dlon,dlat,wlon,rnlat,                                        &
                  data,iret)
         if(iret.gt.0) then
           print *,kpd(5,k,n),kpd(6,k,n),kpd(7,k,n),' read failed'
           call abort
         endif
#ifdef DBG
         print *,'imaxgrb,jmaxgrb=',imaxgrb,jmaxgrb
         print *,'kpds5,kgds1=',kpds5,kgds1
         print *,'proj,orient,truth,cotru=',proj,orient,truth,cotru
         print *,'delx,dely,rlat1,rlat2,rlon1,rlon2=',                         &
                  delx,dely,rlat1,rlat2,rlon1,rlon2
         print *,'dlon,dlat,wlon,rnlat=',dlon,dlat,wlon,rnlat
         print *,'data=',(data(ij),ij=1,10)
#endif
!
!  spatial interpolation 
!
#ifdef DBG
         print *,'imaxgrb,jmaxgrb,dlon,dlat,wlon,rnlat,idim,jdim=',            &
                  imaxgrb,jmaxgrb,dlon,dlat,wlon,rnlat,idim,jdim
#endif
         allocate (slmsk1(imaxgrb,jmaxgrb))
         do j = 1,jmaxgrb
           do i = 1,imaxgrb
             slmsk1(i,j)=1.
           enddo
         enddo
         do j = 1,jdim
           do i = 1,idim
             slmsko(i,j)=1.
           enddo
         enddo
#ifdef RMP
         call sfc_interp_ll2rmp(data,imaxgrb,jmaxgrb,                          &
                   0,out(1,k,n),idim,jdim,                                     &
                   .false.,slmsk1,slmsko,                                      &
                   proj,orient,truth,cotru,                                    &
                   delx,dely,rlat1,rlat2,rlon1,rlon2)
#else
!
         kgau=0
         if(kgds1.eq.4) kgau=1
         call sfc_interp_file_gmp(data,imaxgrb,jmaxgrb,                        &
                  abs(dlon),abs(dlat),wlon,rnlat,                              &
                  out(1,k,n),idim,jdim,.false.,slmsk1,slmsko,kgau,             &
                  0)
#endif
         close(lugb)
         deallocate (slmsk1)
       else
         print *,'== set zero ',kpd(5,k,n),kpd(6,k,n),kpd(7,k,n)
         out(:,k,n)=0.0
       endif
     enddo
   enddo
#ifndef RMP
!
!  input non lat/lon, non gaussian and output global not allowed.
!
   if(abs(proj).ne.3..and.abs(proj).ne.4.) then
     print *,'input non lat/lon non gaussian and output global',' not allowed.'
     call abort
   endif
#endif
#ifdef RMP
!
!  rotate u and v for non lat/lon, non gaussian grid
!
    rad = pi_ / 180.0
    hfpi= 0.5 * pi_
    truth = rtruth
    nproj = rproj
    dlamda0 = orient * rad + hfpi
    ncproj = cproj
    cdlamda0 = corient * rad + hfpi
    call setgrd2(blat,blon,delx,dely,dlamda0)
    call rmp_trans2output_grid(blat,1)
    call rmp_trans2output_grid(blon,1)
    do ij = 1,ijdim
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
        cmap(n)=(1.+sin( nproj * blat(n) ))/a2  ! 1/(map factor)
        a2 = ( 1.0 + sin( nproj * truth * rad ) )
        !  need to define blon blat
        bmap=a2/(1.+sin( nproj * blat(ij) ))
        rgln = blon(ij) - dlamda0
        rgln = nproj * rgln
        csmm = nproj * cos( rgln )
        snmm = nproj * sin( rgln )
      else if( nproj.eq.0 ) then
        cmap(n)=cos( blat(n) )/a2               ! 1/(map factor)
        a2 = cos( truth * rad )
        bmap=a2/cos( blat(ij) )
        csmm =  0.0
        snmm = -1.0
      endif
      if( ncproj.eq.1 .or. ncproj.eq.-1 ) then
        call ll2xy2(blat(ij),blon(ij),x,y)
        call xy2ll2(x,y,clat,clon)
        rgln = clon - cdlamda0
        rgln = ncproj * rgln
        csll = ncproj * cos( rgln )
        snll = ncproj * sin( rgln )
      else if( ncproj.eq.0 ) then
        csll =  0.0
        snll = -1.0
      else if( nproj.eq.2 .or. nproj.eq.-2 ) then
        is = 1
        if( nproj.lt.0 ) is = -1
        truth  = truth * rad
        cotru  = cotru * rad
        orient  = rorient * rad
        dlamda0 = orient
        if( truth.eq.cotru ) then
           cone= cos (pi_*0.5-is*truth)
        else
           cone=(log(cos(truth))-log(cos(cotru)))/                             &
                (log(tan(qtpi-is*truth/2))                                     &
                -log(tan(qtpi-is*cotru/2)))
        endif
        a2 =  rerth_/cone*cos(truth)/(tan(qtpi-is*truth/2))**cone
        bmap = a2/(cos(( is * blat(n) )/cone/(tan(qtpi-is*blat(n)/2)))**cone)
        cmap(n)=reath_/cone/a2*cos(blat(n)/(tan(qtpi-is*blat(n)/2))**cone
        rgln = blon(n) - dlamda0
        rgln = is * rgln
        csmm = is * cos( rgln )
        snmm = is * sin( rgln )
      endif
      if( ncproj.eq.1 .or. ncproj.eq.-1 ) then
        call ll2xy2(blat(ij),blon(ij),x,y)
        call xy2ll2(x,y,clat,clon)
        rgln = clon - cdlamda0
        rgln = ncproj * rgln
        csll = ncproj * cos( rgln )
        snll = ncproj * sin( rgln )
      else if( ncproj.eq.0 ) then
        csll =  0.0
        snll = -1.0
      endif
!
!  need to define cmap
!
      cmorm = sqrt( cmap(ij) )/ bmap
      bcsln(ij) = ( csll * csmm + snll * snmm ) * cmorm
      bsnln(ij) = ( csll * snmm - snll * csmm ) * cmorm
    enddo
!
!  rotate
!
   do k = 1,kdimprs
     do ij = 1,ijdim
       hold = out(ij,k,1)
       out(ij,k,1)=hold*bcsln(ij)-out(ij,k,2)*bsnln(ij)
       out(ij,k,2)=hold*bsnln(ij)+out(ij,k,2)*bcsln(ij)
     enddo
   enddo
#endif
!
   deallocate (data)
#ifdef RMP
   deallocate(rrcsln,rrsnl)
#endif
   deallocate (kpd)
!
#ifdef DBG 
   do n = 1,nfld
     do k = 1,kdimprs
       print *,'k=',k,' n=',n
!       call post_quick_print(out(1,k,n),idim,jdim)
     enddo
   enddo
#endif 
!
   return
   end
