#include "define.h"
   subroutine sfc_cycle_driver(sfcftypin,lugb,iy,im,id,ih,fh)
!-------------------------------------------------------------------------------
!
! subroutine: sfc_cycle_driver
!
!  surface program version 3 
!    to be called inside forecast program
!    complete clean-up.
!    supposed to be easier to add additional variables for new 
!    land models
!
! abstract: this program reads climatology and analysis from grib files
!      interpolates in time and space, merge the two, and then
!      merges with predicted sfc field array in common block
!      the program minimizes the grib file read
!
!  lugb  integer unit number
!  iy .. integer year of initial time
!  im .. integer month of initial time
!  id .. integer day of initial time
!  ih .. integer hour  of initial time
!  fh .. real forecast hour
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
#ifdef SMP_RA2SFC
   use paramodel, only           :  LONF2S=>lonf2scm_,LATG2S=>latg2scm_,       &
                                    LONF2F=>lonf2scm_,LATG2F=>latg2scm_
#else
   use paramodel, only           :  LONF2S,LATG2S,LONF2F,LATG2F
#endif
   use varsfc
   use vargrb
#ifdef MP
   use commpi, only              :  mype,master
#endif
#ifndef RMP
   use module_trans, only        :  dyn_trans2model_grid
   use module_trans, only        :  dyn_trans2output_grid
#else
   use module_trans, only        :  rmp_trans2model_grid
   use module_trans, only        :  rmp_trans2output_grid
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include "abort.h"
!
   integer                      ::  idim,jdim,ijdim
!
   real,allocatable,save        ::  orog(:), slmask(:)
   real,allocatable,save        ::  glacier(:),amxice(:)
   real,allocatable             ::  grbfld(:,:)
#ifdef MP
   real                         ::  grbfldp(LONF2S*LATG2S,numsfcs)
   real                         ::  orogp(LONF2S*LATG2S),slmaskp(LONF2S*LATG2S)
   real                         ::  glacierp(LONF2S*LATG2S),amxicep(LONF2S*LATG2S)
#endif
!
!  namelist variables
!
!  fnmskg .. character file namr of grib high resolution global 
!            lat/lon land ocean mask
!  fnorog .. character file name of binary orography (mtn output) 
!  fnmask .. character file name of binary land sea mask (mtn output) 
!  fngt0  .. character array, file name of grib sfc file for ft=0.
!  fngfc  .. character array, file name of grib sfc file for ft>0.
!  
!
! default file name for initial and forecast
!
   character(len=128)   ::  fnmskg,fnorog,fnmask                           
   data fnmskg/'clim.maskh.grib          '/
   data fnorog/'orogrd.smth              '/
   data fnmask/'slmsk                    '/
!
!   initial sfc fields.  blank indicates not available or simply use
!   what is available in initial sfc file.
!
   character(len=4)                 ::  sfcftypin
   character(len=128)   ::  fnt0(numgrbs)
   data fnt0(igla)/'clim.glacier.grib    '/
   data fnt0(imxi)/'clim.maxice.grib     '/
   data fnt0(imsk)/'                     '/ ! fnmask used
   data fnt0(iais)/'clim.ice.grib        '/
   data fnt0(iscv)/'                     '/ ! not used
   data fnt0(isno)/'clim.snow.grib       '/ ! to be overwritten by anl
   data fnt0(itsf)/'clim.sst.grib        '/ ! to be overwritten by anl
   data fnt0(iab4)/'clim.yhalbedo.grib   '/
   data fnt0(iab1)/'clim.sibalbedo.grib  '/
   data fnt0(ialf)/'clim.yhalbedo.grib   '/
   data fnt0(iso2)/'                     '/ ! not available
   data fnt0(isn1)/'                     '/ ! not available
   data fnt0(izor)/'clim.sibrough.grib   '/
   data fnt0(iplr)/'clim.sibresis.grib   '/
   data fnt0(itg3)/'clim.tg3.grib        '/
   data fnt0(ito2)/'                     '/ ! not available
   data fnt0(itn1)/'                     '/ ! not available
   data fnt0(iveg)/'clim.vegfrac40m.grib '/
   data fnt0(ivet)/'clim.vegtype40m.grib '/
   data fnt0(isot)/'clim.soiltype40m.grib'/
   data fnt0(islo)/'clim.slptyp.grib     '/
   data fnt0(isna)/'clim.snoalb.grib     '/
   data fnt0(ismx)/'clim.shdmax.grib     '/
   data fnt0(ismn)/'clim.shdmin.grib     '/
#ifdef VIC
   data fnt0(ivgv)/'glob.veg1typ.vegcov.grib'/
   data fnt0(ivtv)/'glob.veg1typ.vegtyp.grib'/
   data fnt0(ivrt)/'glob.veg1typ.vegrt.grib'/
   data fnt0(ibif)/'glob.soil.binf.grib  '/
   data fnt0(ids )/'glob.soil.ds.grib    '/
   data fnt0(idsm)/'glob.soil.dsmax.grib '/
   data fnt0(iws )/'glob.soil.ws.grib    '/
   data fnt0(icef)/'glob.soil.cef.grib   '/
   data fnt0(iexp)/'glob.soil.expt.grib  '/
   data fnt0(ikst)/'glob.soil.ksat.grib  '/
   data fnt0(idph)/'glob.soil.dph.grib   '/
   data fnt0(ibub)/'glob.soil.bble.grib  '/
   data fnt0(iqrt)/'glob.soil.qurtz.grib '/
   data fnt0(ibkd)/'glob.soil.blkd.grib  '/
   data fnt0(isld)/'glob.soil.sden.grib  '/
   data fnt0(iwcr)/'glob.soil.wcr.grib   '/
   data fnt0(iwpw)/'glob.soil.wpwp.grib  '/
   data fnt0(ismr)/'glob.soil.smr.grib   '/
   data fnt0(islz)/'glob.soil.silz.grib  '/
   data fnt0(isnz)/'glob.soil.snwz.grib  '/
   data fnt0(ilai)/'glob.veg.lai.grib    '/
   data fnt0(ivg2)/'glob.veg.vegcov.grib '/
   data fnt0(irt1)/'glob.veg.rt1.grib    '/
   data fnt0(irt2)/'glob.veg.rt2.grib    '/
   data fnt0(irt3)/'glob.veg.rt3.grib    '/
#else
   data fnt0(ivgv)/'                     '/
   data fnt0(ivtv)/'                     '/
   data fnt0(ivrt)/'                     '/
   data fnt0(ibif)/'                     '/
   data fnt0(ids )/'                     '/
   data fnt0(idsm)/'                     '/
   data fnt0(iws )/'                     '/
   data fnt0(icef)/'                     '/
   data fnt0(iexp)/'                     '/
   data fnt0(ikst)/'                     '/
   data fnt0(idph)/'                     '/
   data fnt0(ibub)/'                     '/
   data fnt0(iqrt)/'                     '/
   data fnt0(ibkd)/'                     '/
   data fnt0(isld)/'                     '/
   data fnt0(iwcr)/'                     '/
   data fnt0(iwpw)/'                     '/
   data fnt0(ismr)/'                     '/
   data fnt0(islz)/'                     '/
   data fnt0(isnz)/'                     '/
   data fnt0(ilai)/'                     '/
   data fnt0(ivg2)/'                     '/
   data fnt0(irt1)/'                     '/
   data fnt0(irt2)/'                     '/
   data fnt0(irt3)/'                     '/
#endif
   data fnt0(isv1)/'                     '/
   data fnt0(itv1)/'                     '/
   data fnt0(isv2)/'                     '/
   data fnt0(itv2)/'                     '/
!
   data fnt0(ioml)/'clim.omld.grib'/          ! ocean mixed layer depth
   data fnt0(iaer)/'glob.aerosol.paer.grib'/          ! aerosol distribution
   data fnt0(ikpr)/'glob.aerosol.kprf.grib'/          ! aerosol distribution
   data fnt0(iden)/'glob.aerosol.denn.grib'/          ! aerosol distribution
   data fnt0(idxc)/'glob.aerosol.idxc.grib'/          ! aerosol distribution
   data fnt0(imix)/'glob.aerosol.cmix.grib'/          ! aerosol distribution

!
!  forecast sfc fields.  blank indicate time invariant or not available.
!  for time invariant climatology, initial sfc fields needs to be persisted,
!  thus blank to avoid unnecessary io.
!
   character(len=128)   ::  fnfc(numgrbs)
   data fnfc(igla)/'                     '/ ! time invariant
   data fnfc(imxi)/'                     '/ ! time invariant
   data fnfc(imsk)/'                     '/ ! fnmask used 
   data fnfc(iais)/'                     '/ ! monthly varying
   data fnfc(iscv)/'                     '/ ! predicted
   data fnfc(isno)/'                     '/ ! predicted
   data fnfc(itsf)/'                     '/ ! monthly varying
   data fnfc(iab4)/'clim.yhalbedo.grib   '/ ! monthly varying
   data fnfc(iab1)/'clim.sibalbedo.grib  '/ ! monthly varying
   data fnfc(ialf)/'clim.yhalbedo.grib   '/ ! monthly varying
   data fnfc(iso2)/'                     '/ ! predicted
   data fnfc(isn1)/'                     '/ ! predicted
   data fnfc(izor)/'clim.sibrough.grib   '/ ! monthly varying
   data fnfc(iplr)/'clim.sibresis.grib   '/ ! monthly varying
   data fnfc(itg3)/'                     '/ ! time invariant
   data fnfc(ito2)/'                     '/ ! predicted
   data fnfc(itn1)/'                     '/ ! predicted
   data fnfc(iveg)/'clim.vegfrac40m.grib '/ ! monthly varying
   data fnfc(ivet)/'                     '/ ! time invariant
   data fnfc(isot)/'                     '/ ! time invariant
   data fnfc(islo)/'                     '/ ! time invariant
   data fnfc(isna)/'                     '/ ! time invariant
   data fnfc(ismx)/'                     '/ ! time invariant
   data fnfc(ismn)/'                     '/ ! time invariant
   data fnfc(ivgv)/'                     '/ ! time invariant
   data fnfc(ivtv)/'                     '/ ! time invariant
   data fnfc(ivrt)/'                     '/ ! time invariant
   data fnfc(ibif)/'                     '/ ! time invariant
   data fnfc(ids )/'                     '/ ! time invariant
   data fnfc(idsm)/'                     '/ ! time invariant
   data fnfc(iws )/'                     '/ ! time invariant
   data fnfc(icef)/'                     '/ ! time invariant
   data fnfc(iexp)/'                     '/ ! time invariant
   data fnfc(ikst)/'                     '/ ! time invariant
   data fnfc(idph)/'                     '/ ! time invariant
   data fnfc(ibub)/'                     '/ ! time invariant
   data fnfc(iqrt)/'                     '/ ! time invariant
   data fnfc(ibkd)/'                     '/ ! time invariant
   data fnfc(isld)/'                     '/ ! time invariant
   data fnfc(iwcr)/'                     '/ ! time invariant
   data fnfc(iwpw)/'                     '/ ! time invariant
   data fnfc(ismr)/'                     '/ ! time invariant
   data fnfc(islz)/'                     '/ ! time invariant
   data fnfc(isnz)/'                     '/ ! time invariant
#ifdef VIC
   data fnfc(ilai)/'glob.veg.lai.grib    '/ ! monthly varying
#else
   data fnfc(ilai)/'                     '/ ! monthly varying
#endif
   data fnfc(ivg2)/'                     '/ ! time invariant
   data fnfc(irt1)/'                     '/ ! time invariant
   data fnfc(irt2)/'                     '/ ! time invariant
   data fnfc(irt3)/'                     '/ ! time invariant
   data fnfc(isv1)/'                     '/ ! predicted 
   data fnfc(itv1)/'                     '/ ! predicted
   data fnfc(isv2)/'                     '/ ! predicted
   data fnfc(itv2)/'                     '/ ! predicted
!
   data fnfc(ioml)/'clim.omld.grib'/        ! ocean mixed layer depth
   data fnfc(iaer)/'glob.aerosol.paer.grib'/        ! aerosol distribution
   data fnfc(ikpr)/'glob.aerosol.kprf.grib'/        ! aerosol distribution
   data fnfc(iden)/'glob.aerosol.denn.grib'/        ! aerosol distribution
   data fnfc(idxc)/'glob.aerosol.idxc.grib'/        ! aerosol distribution
   data fnfc(imix)/'glob.aerosol.cmix.grib'/        ! aerosol distribution
!
   character(len=128)   ::  fngrb(numgrbs)
   character(len=128)   ::  condir,bindir
   real                 ::  fh
   integer              ::  lugb
   integer              ::  iy,im,id,ih
   integer              ::  is2g(numsfcs)
   integer              ::  ij,k,n
   logical              ::  lsf(numsfcs)
   logical              ::  lfgrb(5,numgrbs)
!
   data condir/'    '/,bindir/'    '/
!
   namelist/namsfc/fnmskg,fnorog,fnmask,                                       &
                   fngrb,condir,bindir,                                        &
                   lfgrb
   save  fnmskg,fnorog,fnmask,fngrb,condir,bindir
!
   integer              ::  ifp,jfp
   data ifp/0/,jfp/0/
   save                     ifp,jfp
!
   idim=LONF2F/2
   jdim=LATG2F*2
   ijdim=idim*jdim
!
   allocate(grbfld(ijdim,numsfcs))
   if(.not.allocated(orog))    allocate(orog(ijdim))
   if(.not.allocated(slmask))  allocate(slmask(ijdim))
   if(.not.allocated(glacier)) allocate(glacier(ijdim))
   if(.not.allocated(amxice))  allocate(amxice(ijdim))
!
#ifdef MP
   if( mype.eq.master ) then
#endif
!
     if(ifp.eq.0) then
       if(fh.eq.0.) then
         do n = 1,numgrbs
           fngrb(n)=fnt0(n)
         enddo
       else
         do n = 1,numgrbs
           fngrb(n)=fnfc(n)
         enddo
       endif
       lfgrb=lflgrb
       read(95,namsfc)
       lflgrb=lfgrb
#ifndef MP
       ifp=1
#endif
     endif
#ifdef MP
   endif
#ifdef RMP
#define MPBCASTL rmpbcastl
#else
#define MPBCASTL mpbcastl
#endif
!
   if (ifp.eq.0) then
     call MPBCASTL(lflgrb,5*numgrbs)
     ifp=1
   endif
#endif
!
!  check input sfc file names for consistency
!  if condir and bindir are given, process them
!
!  note that if filenames starts from '/', it is assumed that
!  full directory is already given and condir, bindir are not added.
!
   call file_check_name(condir,bindir,fngrb,numgrbs,                           &
                fnmskg,fnorog,fnmask)
!
!  translate unpacked grib climatology/analysis file field record 
!  index and sfc file record index.
!     is2g(index of sfc file)=index of grb file
!
   call file_check_index(is2g)
#ifdef DBG
#ifdef MP
   if (mype.eq.master) then
#endif
     print *,'is2g'
     do n = 1,numsfcs
        print *,n,is2g(n)
     enddo
     print *,'ksfc'
     do n = 1,numsfcs
        print *,n,ksfc(n)
     enddo
#ifdef MP
   endif
#endif
#endif
!
!  initialize file exist flag
!
   do n = 1,numsfcs
     lsf(n)=.false.
   enddo
!
   lsf(jtg3)=.true.
   lsf(jtsf)=.true.
!
!  initialize grib read in array
!
   do k = 1,numsfcs
     do ij = 1,idim*jdim
       grbfld(ij,k)=0.
     enddo
   enddo
#ifdef MP
!
   if (mype.eq.master) then
#endif
!
!  read sfc fields from grib (climatology as well as analysis)
!  in the order of forecast sfc file record
!
     call sfc_grib_driver(lugb,idim,jdim,                                      &
                 numsfcs,numsfcv,ksfc,                                         &
                 is2g,                                                         &
                 iy,im,id,ih,fh,                                               &
                 fnmask,fnorog,fnmskg,fngrb,numgrbs,                           &
                 orog,slmask,                                                  &
                 grbfld,lsf)
!
!  process multiple grib files that correspond to a single
!  forecast sfc record.  currently, snow/snow-cover belongs
!  to this category.
!
!  read snow cover field (snow depth is read by default)
!  if snow depth is missing and snow cover is given,
!  do special processing in sfc_merge_field_driver.  set snow depth to -999.,
!  where the snow cover is observed.
!
     if(.not.lsf(jsno).and.fngrb(iscv)(1:4).ne.'    ') then
#ifdef DBG
       print *,'sfc_grib_solver, ngrb,nsfc,fn=',iscv,jsno,trim(fngrb(iscv))
#endif
       call sfc_grib_solver(lugb,fngrb(iscv),idim,jdim,slmask,                 &
                    iscv,iy,im,id,ih,fh,                                       &
                    fnmskg,grbfld(1,jsno))                                     
       lsf(jsno)=.true.
!
!  in order to distinguish from actual snow depth, set snow mask
!  as -999.
!
       do ij = 1,idim*jdim
         grbfld(ij,jsno)=grbfld(ij,jsno)*(-999.)
       enddo
     endif
!
!  obtain glacier points and maximum possible ice extent from
!  grib file
!
!  get glacier points from grib input file
!
     if(jfp.eq.0) then
#ifdef DBG
       print *,'sfc_grib_solver, ngrb,nsfc,fn=',igla,jsno,trim(fngrb(igla))
       print *,'igla sfc_grib_solver, igla=',igla
#endif
       if(fngrb(igla)(1:4).eq.'    ') then
         print *,'glacier grib file name empty'
         call MPABORT
       endif
       call sfc_grib_solver(lugb,fngrb(igla),idim,jdim,slmask,                 &
                    igla,iy,im,id,ih,fh,                                       &
                    fnmskg,glacier)
!
!  get maximum possible sea ice extent from grib file
!
#ifdef DBG
       print *,'sfc_grib_solver, ngrb,nsfc,fn=',imxi,trim(fngrb(imxi))
       print *,'imxi sfc_grib_solver, imxi=',imxi
#endif
       if(fngrb(imxi)(1:4).eq.'    ') then
         print *,'max ice grib file name empty'
         call MPABORT
       endif
       call sfc_grib_solver(lugb,fngrb(imxi),idim,jdim,slmask,                 &
                    imxi,iy,im,id,ih,fh,                                       &
                    fnmskg,amxice)
       jfp=1
     endif
!
!  rearrange all the arrays to transform format before merge
!
#ifndef RMP
     call dyn_trans2model_grid(orog,1)
     call dyn_trans2model_grid(slmask,1)
     call dyn_trans2model_grid(glacier,1)
     call dyn_trans2model_grid(amxice,1)
     call dyn_trans2model_grid(grbfld,numsfcs)
#else
     call rmp_trans2model_grid(orog,1)
     call rmp_trans2model_grid(slmask,1)
     call rmp_trans2model_grid(glacier,1)
     call rmp_trans2model_grid(amxice,1)
     call rmp_trans2model_grid(grbfld,numsfcs)
#endif
!
#ifdef MP
   endif
#ifndef RMP
#define MPBCASTL mpbcastl
#define MPBCASTI mpbcasti
#define MPBCASTR mpbcastr
#define MPGF2P   mpgf2p
#define MPSYNALL mpsynall
#else
#define MPBCASTL rmpbcastl
#define MPBCASTI rmpbcasti
#define MPBCASTR rmpbcastr
#define MPGF2P   rmpgf2p
#define MPSYNALL rmpsynall
#endif
!
   call MPBCASTL(lsf    ,numsfcs)
   call MPBCASTI(is2g   ,numsfcs)
   !call MPBCASTR(orog   ,  ijdim)
   !call MPBCASTR(slmask ,  ijdim)
   !call MPBCASTR(glacier,  ijdim)
   !call MPBCASTR(amxice ,  ijdim)
   !call MPBCASTR(grbfld ,ijdim*numsfcs)
!
   call MPGF2P(orog  ,LONF2F,LATG2F,  orogp,LONF2S,LATG2S,1)
   call MPGF2P(slmask,LONF2F,LATG2F,slmaskp,LONF2S,LATG2S,1)
   call MPGF2P(glacier,LONF2F,LATG2F,glacierp,LONF2S,LATG2S,1)
   call MPGF2P(amxice,LONF2F,LATG2F,amxicep,LONF2S,LATG2S,1)
   call MPGF2P(grbfld,LONF2F,LATG2F,grbfldp,LONF2S,LATG2S,numsfcs)
#endif
!
!  merge with forecast sfc array
!
#ifdef MP
#define GRBFLD grbfldp
#define OROG   orogp
#define SLMASK slmaskp
#define GLACIER glacierp
#define AMAXICE amxicep
#else
#define GRBFLD grbfld
#define OROG   orog
#define SLMASK slmask
#define GLACIER glacier
#define AMAXICE amxice
#endif
#ifndef SKIPSFCMRG
   call sfc_merge_field_driver(sfcftypin,GRBFLD,LONF2S,LATG2S,iy,im,id,ih,fh,  &
               OROG,SLMASK,GLACIER,AMAXICE,lsf,is2g)
#else
#ifdef SMP
#ifdef SMP_RA2SFC
   call sfc_merge_field_driver(sfcftypin,GRBFLD,LONF2S,LATG2S,iy,im,id,ih,fh,  &
               OROG,SLMASK,GLACIER,AMAXICE,lsf,is2g)
#endif
#endif
#endif
!
!  rearrange back to input format 
!
#ifndef RMP
#ifdef MP
   if(mype.eq.master) then
#endif
     call dyn_trans2output_grid(orog,1)
     call dyn_trans2output_grid(slmask,1)
     call dyn_trans2output_grid(glacier,1)
     call dyn_trans2output_grid(amxice,1)
#ifdef MP
   endif
#endif
#else
#ifdef MP
   if(mype.eq.master) then
#endif
     call rmp_trans2output_grid(orog,1)
     call rmp_trans2output_grid(slmask,1)
     call rmp_trans2output_grid(glacier,1)
     call rmp_trans2output_grid(amxice,1)
#ifdef MP
   endif
#endif
#endif            /* RMP */
!
   deallocate(grbfld)
!
   return
   end subroutine sfc_cycle_driver
!-------------------------------------------------------------------------------
