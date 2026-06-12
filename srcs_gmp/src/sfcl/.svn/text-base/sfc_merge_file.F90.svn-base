#include "define.h"
   subroutine sfc_merge_file
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [sfc_merge_file]
!           |
!           |-- [sfc_merge_field_driver] *
!           |-- [sfc_merge_field] *
!           |-- [sfc_merge_snow] *
!
!-------------------------------------------------------------------------------
   end subroutine sfc_merge_file
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_merge_field_driver(sfcftypin,grbfld,idim,jdim,               &
                     iy,im,id,ih,fh,orog,slmask,glacier,amxice,lsf,is2g)
!#define BIN_DBG
!-------------------------------------------------------------------------------
!
! subroutine: sfc_merge_field_driver
!
! abstract:  this program takes grib climatology and analysis fields,
!    merge them with the surface varible array in the forecast program.
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          cvs version, physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only         : mtnvar_,ntotal_,nwater_,ngases_
   use constant, only          : pi_
   use vargrb
   use varsfc
   use comsfc                  ! sfcfcs
#ifndef RMP
   use comfgrid
   use module_trans, only      : dyn_trans2model_grid, dyn_trans2output_grid
#ifdef MP
   use comfphys
#endif
#else
   use module_trans, only      : rmp_trans2output_grid
#ifdef MP
   use rscommap
#endif
#endif
#ifdef BIN_DBG
   use module_file_write, only : file_write_bin
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include <abort.h>
#ifdef MP
#undef BIN_DBG
#undef DBG
#endif
   character(len=4)     ::  sfcftypin
   character(len=8)     ::  svar
   integer              ::  idim,jdim,iy,im,id,ih,ijdim
   real                 ::  fh
   real                 ::  grbfld(idim*jdim,numsfcs)
   real                 ::  orog(idim*jdim),slmask(idim*jdim)
   real                 ::  glacier(idim*jdim),amxice(idim*jdim)
   logical              ::  lsf(numsfcs)
   integer              ::  is2g(*)
#ifdef USGS
   integer              ::  jy,jm,jd,jh
#endif
   integer, allocatable :: loismsk(:,:)
#if defined(DBG) || defined(BIN_DBG)
   character(len=20)    :: fmt
   real, allocatable    :: rloismsk(:,:)
#endif
!
!  character(len=4)     ::  sfcftyp
   real, parameter      ::  rlapse=0.65e-2
   real                 ::  rjday
   integer              ::  ij,nsfc,nsfcv,ngrb,k
#ifndef MP
   real                 ::  rlat(idim,jdim),rlon(idim,jdim)
#endif
   real                 ::  rdelx,rdely,dlamda0
   real                 ::  smcfracf, smcfracg
   integer              ::  i,j,isep
!
#ifdef SM_SPINUP
   real                 ::  sotmod(idim*jdim)
#ifdef USGS_SFC
   real                 ::  smref(16),smwlt(16)
   data smref/0.236, 0.283, 0.312, 0.360, 0.360, 0.329, 0.314,                 &
              0.387, 0.382, 0.338, 0.404, 0.412, 0.329, 0.0, 0.108, 0.283/
   data smwlt/0.010, 0.028, 0.047, 0.084, 0.084, 0.066, 0.067,                 &
              0.120, 0.103, 0.100, 0.126, 0.138, 0.066, 0.0, 0.006, 0.028/
#else
   real                 ::  smref(9),smwlt(9)
   data smref/.283,.387,.412,.312,.338,.382,.315,.329,.283/
   data smwlt/.029,.119,.139,.047,.010,.103,.069,.066,.029/
#endif
#endif
   logical              ::  lmerge(numsfcs)
!
   allocate (loismsk(idim,jdim))
   ijdim = idim*jdim
!
   isep=1
!
#ifdef BIN_DBG
   call file_write_bin(129,slmask,idim,jdim,1,isep)
#ifndef OSULSM1
   call file_write_bin(129,grbfld(1,jtsf),idim,jdim,1,isep)
   call file_write_bin(129,grbfld(1,jveg),idim,jdim,1,isep)
   call file_write_bin(129,sfcfcs(1,jsli),idim,jdim,1,isep)
#endif
#endif
!
!  temporary, set loismsk to land sea mask
!
   ij=0
   do j = 1,jdim
     do i = 1,idim
       ij=ij+1
       loismsk(i,j)=nint(slmask(ij))
     enddo
   enddo
!
!  quality control of sea ice
!
   call sfc_check_seaice(grbfld(1,jsli),glacier,amxice,slmask,idim,jdim)
!
!  merge sea ice
!
!    always use slmask as land sea mask, discarding the mask
!    on original sfc file
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).lt.2.) sfcfcs(ij,jsli)=slmask(ij)
   enddo
#ifdef BIN_DBG
   call file_write_bin(129,slmask,idim,jdim,1,isep)
   call file_write_bin(129,sfcfcs(1,jsli),idim,jdim,1,isep)
#endif
!
!    then add new sea ice
!
   if (lflgrb(3,iais)) then
     do ij = 1,ijdim
       if(grbfld(ij,jsli).eq.1.) sfcfcs(ij,jsli)=2.
       if(grbfld(ij,jsli).eq.0 ) sfcfcs(ij,jsli)=slmask(ij)
     enddo
   endif
#ifdef BIN_DBG
   call file_write_bin(129,slmask,idim,jdim,1,isep)
   call file_write_bin(129,sfcfcs(1,jsli),idim,jdim,1,isep)
#endif
!
!  temporary, set loismsk to land sea sea-ice mask
!
   ij=0
   do j = 1,jdim
     do i = 1,idim
       ij=ij+1
       if(sfcfcs(ij,jsli).eq.2.) then
         loismsk(i,j)=2
       endif
     enddo
   enddo
!
!  process snow
!   if snow depth available (snow cover analysis not used)
!
   ngrb=is2g(jsno)
!
   if(lsf(jsno).and.ngrb.ne.9999) then
     ij=1
     do while (grbfld(ij,jsno).ne.-999..and.ij.lt.ijdim)
       ij=ij+1
     enddo
     if(ij.eq.ijdim) then
       call sfc_merge_field(sfcfcs(1,jsno),grbfld(1,jsno),                     &
                     ijdim,ksfc(jsno),                                         &
                     lsf(jsno),lflgrb(1,ngrb),loismsk)
#ifdef NOALSM1
       do ij = 1,ijdim
         sfcfcs(ij,jsnd)=sfcfcs(ij,jsno)*5.
       enddo
#endif
     else
!
!  merge snow
!
       call sfc_merge_snow(grbfld(1,jsno),sfcfcs(1,jsno),sfcfcs(1,jtsf),       &
#ifdef NOALSM1
                      sfcfcs(1,jsnd),                                          &
#endif
                      ijdim)
     endif
   endif
!
!  generate loismsk
!        0=sea
!        1=snow free land
!        2=snow free sea-ice
!        3=land with snow
!        4=sea ice with snow mask
!
   ij=0
   do j = 1,jdim
     do i = 1,idim
       ij=ij+1
       if(sfcfcs(ij,jsno).gt.0.) then
         if(loismsk(i,j).eq.1) then
           loismsk(i,j)=3
         elseif(loismsk(i,j).eq.2) then
           loismsk(i,j)=4
         endif
       endif
     enddo
   enddo
!
#ifdef DBG
   allocate (rloismsk(idim,jdim))
   do j = 1,jdim
     do i = 1,idim
       rloismsk(i,j)=loismsk(i,j)
     enddo
   enddo
#ifndef RMP
   call dyn_trans2output_grid(rloismsk,1)
#else
   call rmp_trans2output_grid(rloismsk,1)
#endif
   call sfc_quick_print(rloismsk,idim,jdim)
#endif
#ifdef BIN_DBG
   if (.not.allocated(rloismsk)) allocate (rloismsk(idim,jdim))
   do j = 1,jdim
     do i = 1,idim
       rloismsk(i,j)=loismsk(i,j)
     enddo
   enddo
   call file_write_bin(129,rloismsk,idim,jdim,1,isep)
#endif
!
! make sst and tg3 correction before merge
!
   if(lsf(jtg3)) then
     call sfc_tsfc_correc(grbfld(1,jtg3),orog,slmask,1.,ijdim)
   endif
   if(lsf(jtsf)) then
     call sfc_tsfc_correc(grbfld(1,jtsf),orog,slmask,0.,ijdim)
   endif
!
!  now ready to merge 
!
   nsfc=1
   do nsfcv = 1,numsfcv
     lmerge(nsfc)=.true.
     if(nsfc.eq.jsli) lmerge(nsfc)=.false.
     if(nsfc.eq.jsno) lmerge(nsfc)=.false.
#ifdef VICLSM1
     if(nsfc.eq.jzor) lmerge(nsfc)=.false.
#endif
#ifdef ANL_SOIL
     if(nsfc.eq.jsmc) lmerge(nsfc)=.false.
#endif
     nsfc=nsfc+ksfc(nsfc)
   enddo
!
!  soil moisture  
!
#ifdef SM_SPINUP
#ifndef OSULSM1
   if (fh.eq.0..and.sfcftypin.eq.'osu1') then
#ifdef DBG
     print*, 'soil moisture spinup at fh  = ',fh
#endif
     call sot2one(grbfld(1,jsot),sotmod,ijdim)
     do k = jsmc,jsmc+lsoil_-1
       do ij = 1,ijdim
         if(sfcfcs(ij,k).gt.0.47) sfcfcs(ij,k)=0.47
         if(sfcfcs(ij,k).lt.0.1) sfcfcs(ij,k)=0.1
         if(grbfld(ij,k).gt.0.47) grbfld(ij,k)=0.47
         if(grbfld(ij,k).lt.0.1) grbfld(ij,k)=0.1
         if(sotmod(ij).ne.0.) then
#ifdef DBG
           print*,'read sotmod(ij)=',sotmod(ij),sfcfcs(ij,k),                  &
                       grbfld(ij,k)
#endif
           smcfracf= (sfcfcs(ij,k) - 0.1)/(0.47 - 0.1)
           smcfracg= (grbfld(ij,k) - 0.1)/(0.47 - 0.1)
#ifdef WET_SOIL
           smcfracf=1.
           smcfracg=1.
#endif
#ifdef DRY_SOIL
           smcfracf=0.
           smcfracg=0.
#endif
           sfcfcs(ij,k)=smcfracf*(smref(sotmod(ij))-smwlt(sotmod(ij)))         &
                           +smwlt(sotmod(ij))
           grbfld(ij,k)=smcfracg*(smref(sotmod(ij))-smwlt(sotmod(ij)))         &
                           +smwlt(sotmod(ij))
#ifdef DBG
           print*,'adjusted sotmod(ij)=',sotmod(ij),sfcfcs(ij,k),              &
                  grbfld(ij,k),smcfracf,smref(sotmod(ij)),smwlt(sotmod(ij))
#endif
#ifdef CLM_SOIL
           if(slmask(ij).eq.1) sfcfcs(ij,k) = grbfld(ij,k)
#endif
         endif
       enddo
     enddo
   endif
#else
   if (fh.eq.0.) then
     do k = jsmc,jsmc+lsoil_-1
       do ij = 1,ijdim
#ifdef WET_SOIL
         if(slmask(ij).eq.1) sfcfcs(ij,k) = 0.47
#endif
#ifdef DRY_SOIL
         if(slmask(ij).eq.1) sfcfcs(ij,k) = 0.1
#endif
#ifdef CLM_SOIL
         if(slmask(ij).eq.1) sfcfcs(ij,k) = grbfld(ij,k)
#endif
       enddo
     enddo
   endif
#endif
!
#endif
!
   nsfc=1
   do nsfcv = 1,numsfcv
     if(lmerge(nsfc)) then ! skp merge done earlier.
       ngrb=is2g(nsfc)
       if(lsf(nsfc).and.ngrb.ne.9999) then
#ifdef DBG
#ifndef OSULSM1
         if(nsfc.eq.jvet) then
           call print_maxmin_six(grbfld(1,jvet),idim,jdim,1,jdim,              &
                                    'before grb veg')
           call print_maxmin_six(sfcfcs(1,jvet),idim,jdim,1,jdim,              &
                                    'before sfc veg')
         endif
#endif
#endif
         call sfc_merge_field(sfcfcs(1,nsfc),grbfld(1,nsfc),                   &
                        ijdim,ksfc(nsfc),                                      &
                        lsf(nsfc),lflgrb(1,ngrb),loismsk)
#ifdef DBG
#ifndef OSULSM1
         if(nsfc.eq.jvet) then
           call print_maxmin_six(grbfld(1,jvet),idim,jdim,1,jdim,              &
                                    'after vegt')
           call print_maxmin_six(sfcfcs(1,jvet),idim,jdim,1,jdim,              &
                                    'after vegt2')
         endif
#endif
#endif
       endif
     endif
     nsfc=nsfc+ksfc(nsfc)
   enddo
!
#ifdef ANL_SOIL
   do k = jsmc,jsmc+lsoil_-1
     do ij = 1,ijdim
       if(grbfld(ij,k).gt.0.47) grbfld(ij,k)=0.47
       if(grbfld(ij,k).lt.0.1) grbfld(ij,k)=0.1
       if(slmask(ij).eq.1) sfcfcs(ij,jsmc+1) = grbfld(ij,jsmc+1)
     enddo
   enddo
#endif
#ifdef BIN_DBG
   call file_write_bin(129,sfcfcs(1,1),idim,jdim,1,isep)
#ifndef OSULSM1
   call file_write_bin(129,sfcfcs(1,jvet),idim,jdim,1,isep)
   call file_write_bin(129,sfcfcs(1,jalb),idim,jdim,1,isep)
#endif
#endif
#if defined NOA || defined OSULSM2
!
!  quality control of vegataion and soil type
!
   call sfc_land_type(sfcfcs(1,jvet),sfcfcs(1,jsot),sfcfcs(1,jsli),            &
                      idim,jdim)
#else
!
!  check vic land surface parameters
!
   call sfc_check_vicland(idim,jdim)
#endif
#ifdef USGS
!
!  in case of USGS, zor and albedo are computed from vegetation type
!
   call sfc_usgs_rough(sfcfcs(1,jvet),ijdim,sfcfcs(1,jzor))
   call sfc_interp_date(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
#ifdef DBG
   print*, 'before dyn_trans2model_grid'
   call print_maxmin_six(sfcfcs(1,jalb),idim,jdim,1,jdim,'before albedo')
   call print_maxmin_six(sfcfcs(1,jvet),idim,jdim,1,jdim,'jvet bef albedo')
#endif
#ifndef MP
#ifdef RMP
   call rmp_setup_rsm_grid(rlat,rlon,rdelx,rdely,dlamda0)
#else
   call dyn_lon_lat(rlon,rlat,colrad,lonf_,latg_)
   call dyn_trans2model_grid(rlon,1)
   call dyn_trans2model_grid(rlat,1)
   do j = 1,latg2_
     do i = 1,lonf2_
       rlon(i,j)=rlon(i,j)*180./pi_
       rlat(i,j)=rlat(i,j)*180./pi_
     enddo
   enddo
#ifdef DBG
   print*, 'after dyn_trans2model_grid'
   call print_maxmin_six(sfcfcs(1,jalb),idim,jdim,1,jdim,'before albedo')
   call print_maxmin_six(sfcfcs(1,jvet),idim,jdim,1,jdim,'jvet bef albedo')
#endif
#endif
#endif
!
#ifdef MP
#ifdef RMP
#define RLAT rlat
#else
#define RLAT xlat
#endif
#else
#define RLAT rlat
#endif
   call sfc_usgs_albedo(sfcfcs(1,jvet),RLAT,ijdim,jm,                          &
                  sfcfcs(1,jalb),sfcfcs(1,jalf))
#ifdef DBG
   call print_maxmin_six(sfcfcs(1,jvet),idim,jdim,1,jdim,'jvet aft albedo')
   call print_maxmin_six(sfcfcs(1,jalb),idim,jdim,1,jdim,'after albedo')
#ifdef BIN_DBG
   call file_write_bin(129,sfcfcs(1,jvet),idim,jdim,1,isep)
   call file_write_bin(129,sfcfcs(1,jalb),idim,jdim,1,isep)
#endif
#endif
#endif
!
!  quality control by checking max/min over various surfaces
!
   nsfc=1
   do nsfcv = 1,numsfcv
     write(svar,100) nsfcv
100  format('nv=',i3)
     do k = 1,ksfc(nsfc)
       call sfc_check_maxmin(sfcfcs(1,nsfc+k-1),ijdim,loismsk,                 &
                     vmaxmin(1,1,nsfc),svar)
     enddo
     nsfc=nsfc+ksfc(nsfc)
   enddo
#ifdef BIN_DBG
   call file_write_bin(129,sfcfcs(1,1),idim,jdim,1,isep)
#endif
#ifdef NOALSM1
!
!  check snodpth to make sure that it is nonzero if snow is nonzero
!
   call sfc_check_snowdepth(sfcfcs(1,jsno),sfcfcs(1,jsnd),ijdim)
#endif
   deallocate (loismsk)
#ifdef DBG
   deallocate (rloismsk)
#endif
!
   return
   end subroutine sfc_merge_field_driver
!-------------------------------------------------------------------------------
#ifdef SM_SPINUP
!
!-------------------------------------------------------------------------------
   subroutine sot2one(sotclmx,sotmod,ijdim)
!-------------------------------------------------------------------------------
   real                 ::   sotclmx(ijdim),sotmod(ijdim)
!
   sotmod=sotclmx
!
   return
   end subroutine sot2one
!-------------------------------------------------------------------------------
#endif
!
!-------------------------------------------------------------------------------
   subroutine sfc_merge_field(sfcfcs,grbfld,                                   &
                     ijdim,kdim,                                               &
                     lsf,lflgrb,loismsk)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
#include "define.h"
   integer        ::  ijdim,kdim
   real           ::  sfcfcs(ijdim,kdim)
   real           ::  grbfld(ijdim,kdim)
   integer        ::  loismsk(ijdim)
   logical        ::  lsf,lflgrb(5)
   logical        ::  lall
   integer        ::  i,ij,k
#ifdef DBG
   integer        ::  nmrg
!
   nmrg=0
#endif
!
   if(.not.lsf) return
   lall=.true.
   do i = 1,5
     if(.not.lflgrb(i)) lall=.false.
   enddo
!
   if(lall) then
     do k = 1,kdim
       do ij = 1,ijdim
         sfcfcs(ij,k)=grbfld(ij,k)
       enddo
     enddo
#ifdef DBG
     print*,'in sfc_merge_field all ',ijdim,'x',kdim,'merged'
#endif
     return
   endif
!
   do i = 1,5
     if(lflgrb(i)) then
       do k = 1,kdim
         do ij = 1,ijdim
           if(loismsk(ij).eq.i-1) then
#ifdef DBG
             nmrg=nmrg+1
#endif
             sfcfcs(ij,k)=grbfld(ij,k)
           endif
         enddo
       enddo
     endif
   enddo
!
#ifdef DBG
   print*,'in sfc_merge_field:nmrg=',nmrg
#endif
!
   return
   end subroutine sfc_merge_field
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_merge_snow(grbsnow,sfcsnow,sfctemp,                          &
#ifdef NOALSM1
                sfcsnod,                                                       &
#endif
                ijdim)
!-------------------------------------------------------------------------------
!
!   if snow cover available but no snow depth available
!     1. set fcst snow depth to zero where snow cover mask is 0
!     2. use fcst snow depth where snow cover mask is 1
!     3. compute snow depth where fcst snow cover mask is 0
!
!-------------------------------------------------------------------------------
   implicit none  
!-------------------------------------------------------------------------------
   real, parameter      ::  snwmin=25.0,snwmax=100.
   integer              ::  ijdim
   real                 ::  grbsnow(ijdim),sfcsnow(ijdim),sfctemp(ijdim)
#ifdef NOALSM1
   real                 ::  sfcsnod(ijdim)
#endif
   integer              ::  ij
!
   do ij = 1,ijdim
!
! 1.
!
     if(grbsnow(ij).eq.0.) then
       sfcsnow(ij) = 0.
#ifdef NOALSM1
       sfcsnod(ij) = 0.
#endif
     endif
!
! 2.
!    if(grbsnow(ij).eq.1..and.sfcsnow(ij).gt.0.) then
!      sfcsnow(ij)=sfcsnow(ij)
!    endif
! 3.
!
     if(grbsnow(ij).eq.-999..and.sfcsnow(ij).eq.0.) then
       if(sfctemp(ij).lt.243.0) then
         sfcsnow(ij) = snwmax
       elseif(sfctemp(ij).lt.273.0) then
         sfcsnow(ij) = snwmin+                                                 &
                       (snwmax-snwmin)*(273.0-sfctemp(ij))/30.0
       else
         sfcsnow(ij) = snwmin
       endif
#ifdef NOALSM1
       sfcsnod(ij) = sfcsnow(ij)*5.
#endif
     endif
   enddo
!
   return
   end subroutine sfc_merge_snow
!-------------------------------------------------------------------------------
