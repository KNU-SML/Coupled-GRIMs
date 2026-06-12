#include "define.h"
   program change_resolution_pressure
!#define BIN_DBG
#undef SFC_CORRECTION
!-------------------------------------------------------------------------------
!
! subprogram:    change_resolution_pressure      
!                           vertically interpolate standard pressure level
!                           data to sigma level and convert to either
!                           spectral coefficients or RMP grid suitable
!                           for G-RMP input.
!
!                           This program assumes that all the input standard
!                           pressure level fields are given in one grib format
!                           file (file name fn given as a namelist input).
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
! abstract: this routine spatially interpolates standard pressure level data
!           to forecast model gaussian grid, or RMP grid, interpolates surface
!           pressure to the new orography, vertically interpolates wind,
!           temperature and moisture, fixes moisture in the stratosphere if
!           necessary, then transform the data either to spectral coeffs
!           or RMP grid suitable for G-RMP input.
!           if this program is used to create lateral forcing for RMP
!           on RMP grid, you need to run this program using larger domain
!           RMP as a model configuration.
!
!-------------------------------------------------------------------------------
   use constant, only     : g_,pi_,rd_,rv_,qmin_
#ifdef HYBRID
   use constant, only     : akapa_
#endif
#ifdef DFS
   use dfsvar, only       : get_dfs_dim,mt,jl,jla,jlg,coslat,jge
#else
   use module_sph_legendre, only : sph_gaussian_lat, sph_poly_epsilon1
#endif
   use paramodel, only    : para_init, latg_,lonf_,ngases_,nwater_,            &
                            ntotal_,jcap_,levs_
   use comchgr, only      : comchgr_init,                                      &
                            idim,jdim,idimt,jdimhf,ijdim,                      &
                            kdim,kdimp,kdimm,mdimv,indxmm,                     &
                            dio,zeo,teo,rqo,gzo,qo
#ifdef DFS
!
! for dfs
!
   use comchgr, only      : mdim=>mdimd
#ifdef ALIASED
   use comchgr, only      : mw
#endif
#else
!
! for sph
!
   use comchgr, only      : mdim,eps,colrad,wgt,wgtcs,rcs2,                    &
                            zss,pss,tts,qqs,uus,vvs,                           &
                            zsa,psa,tta,qqa,uua,vva
#endif
   use module_trans, only : dyn_trans2model_grid, dyn_trans2output_grid
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                 ::  cslat,dlat
   real, allocatable    ::  xlat(:)
   real, allocatable    ::  ps(:),zs(:),expps(:)
#ifdef HYBRID
   real, allocatable    ::  exppso(:)
#endif
   real, allocatable    ::  spino(:,:),splno(:,:)
   real, allocatable    ::  fsig(:,:)
   real, allocatable    ::  qe(:,:)
!
   character(len=32)    ::  lab
   real                 ::  waves,xlayers,trun,order,realform,gencode
   real                 ::  rlond,rlatd,rlonp,rlatp,rlonr,rlatr,gases
   real                 ::  water,pdryini,subcen,ppid,slid,vcid,vmid,vtid
   real                 ::  runid,usrid
   integer,parameter    ::  kdum2=21,kens=2,nfld=5
   integer              ::  kdum
   real, allocatable    ::  dummy(:)
   real                 ::  dummy2(kdum2),ensemble(kens)
   integer              ::  idate(4)
   real                 ::  fhour
!
!  arrays for global spectral conversion
!
#ifdef DFS
   real, allocatable    ::  uua(:,:,:),vva(:,:,:)
   real, allocatable    ::  psa(:,:),tta(:,:,:),qqa(:,:,:)
   real, allocatable    ::  zsa(:,:)
   real, allocatable    ::  grid(:,:)
   real, allocatable    ::  tave(:)
   real                 ::  rlat
#else
   real, allocatable    ::  qlnt(:),qlnv(:),qdert(:),qlnwct(:)
   real, allocatable    ::  b(:)
#endif
!
!  model orography
!
   real, allocatable    ::  gzso(:)
   real, allocatable    ::  zso(:,:)
#ifdef SFC_CORRECTION
   real, allocatable    ::  sfcti(:),sfcto(:,:)
   real                 ::  wei
#endif
!
!  model sigma level
!
   real, allocatable    ::  si(:),sl(:),del(:)
   real, allocatable    ::  ci(:),cl(:),rpi(:)
#ifdef HYBRID
   real, allocatable    ::  ak5(:),bk5(:)
   real, allocatable    ::  ak5x(:),bk5x(:)
   real                 ::  rk1,rkinv
   real                 ::  dif
#endif
!
#if defined(DBG) || defined(BIN_DBG)
   real, allocatable    ::  zout(:,:),pout(:,:)
   real, allocatable    ::  tout(:,:,:),qout(:,:,:)
   real, allocatable    ::  uout(:,:,:),vout(:,:,:)
#endif
!
   real,parameter       ::  tensn=10.0
   real,parameter       ::  gamma=6.5e-3
   real                 ::  ps1,rh1,psk
!
   integer              ::  iy,im,id,ih,ifmoisrh,ifheight,ifmtngrb
   integer              ::  iupd5,ivpd5,itpd5,iqpd5,izpd5 
   real                 ::  fh
   integer,parameter    ::  maxstdprs=100
   real                 ::  stdprs(maxstdprs)
   character(len=128)   ::  fn,fnmtngrb
   data ifmoisrh/0/
   data ifheight/1/
   data iupd5,ivpd5,itpd5,iqpd5,izpd5 /33,34,11,52,7/
!
   integer              ::  kdimprs,ij,iret
!
   real es,fpvs0,qs,xeps,xepsm1
   parameter(xeps=rd_/rv_,xepsm1=rd_/rv_-1.)
!
   integer i,j,in,in2,k,n,m,lat,kdimi
!
   real, allocatable :: out(:,:,:,:),out2(:,:,:)
   real, allocatable :: slin(:)
   real, allocatable :: siin(:)
!
!   iugrb : unit number for grib read file
!   fn   : input standard pressure level grid file name
!   iy,im,id,ih,fh   : input year,month,day,hour,forecast hour
!   stdprs : standard pressure level to be read in mb unit
!   kdimprs : number of standard pressure levels
!   ifmoisrh : =0 specific humidity ; =1 relative humidity
!   ifheight : =1 ncep  ; =0 ecmwf grib file
!   iorsmgsm : choice of output. =0 global spectral ; =1 regional
!
   integer kdimq,ifcubic
   integer iugrb,iumtn,iusfc,iuatm,iusig              ! input & output file unit
   data iugrb,iumtn,iusfc,iuatm,iusig/10,11,14,15,51/
   data ifmtngrb,kdimq/0,0/,ifcubic/1/
   namelist/namp2s/fn,iy,im,id,ih,fh,stdprs,ifmoisrh,ifheight,                 &
                   iupd5,ivpd5,itpd5,iqpd5,izpd5,ifcubic,                      &
                   iugrb,iumtn,iusfc,iuatm,iusig,ifmtngrb,fnmtngrb,kdimq
!-------------------------------------------------------------------------------
#ifdef RMP
   print *,'This program does not work under RMP model configuration.'
   call abort
#endif
   call para_init
   call comchgr_init
!
#ifndef HYBRID
   kdum=201-kdim-1-kdim
#else
   kdum=201-kdim-2-kdim
   rk1 = akapa_ + 1.
   rkinv=1./akapa_
#endif
!
   allocate(xlat(idimt))
   allocate(ps(idimt),zs(idimt),expps(idimt))
#ifdef HYBRID
   allocate(exppso(idimt))
#endif
   allocate(spino(idimt,kdimp),splno(idimt,kdim))
   allocate(fsig(idimt,kdim*nfld+2))
   allocate(qe(idimt,kdim))
   allocate(dummy(kdum))
#ifdef DFS
   allocate(uua(idim,jdim,kdim),vva(idim,jdim,kdim))
   allocate(psa(idim,jdim),tta(idim,jdim,kdim),qqa(idim,jdim,kdim))
   allocate(zsa(idim,jdim))
   allocate(grid(idim*jdim,kdim))
   allocate(tave(kdim))
#else
   allocate(qlnt(mdim),qlnv(mdimv),qdert(mdim),qlnwct(mdim))
   allocate(b(mdim))
#endif
!
   allocate(gzso(mdim))
   allocate(zso(idimt,jdimhf))
#ifdef SFC_CORRECTION
   allocate(sfcti(ijdim),sfcto(idimt,jdimhf))
#endif
   allocate(si(kdimp),sl(kdim),del(kdim))
   allocate(ci(kdimp),cl(kdim),rpi(kdim-1))
#ifdef HYBRID
   allocate(ak5(kdimp),bk5(kdimp))
   allocate(ak5x(kdimp),bk5x(kdimp))
#endif
#if defined(DBG) || defined(BIN_DBG)
   allocate(zout(idimt,jdimhf),pout(idimt,jdimhf))
   allocate(tout(idimt,jdimhf,kdim),qout(idimt,jdimhf,kdim))
   allocate(uout(idimt,jdimhf,kdim),vout(idimt,jdimhf,kdim))
#endif
!
   do k = 1,maxstdprs
     stdprs(k)=0.
   enddo
!
   read (5,namp2s)
   write(6,namp2s)
!
   kdimprs=1
   do while(stdprs(kdimprs).ne.0.)
     kdimprs=kdimprs+1
   enddo
   kdimprs=kdimprs-1
   if (kdimq.eq.0) kdimq=kdimprs
!
!  define constants for global spherical transform
!
#ifdef DFS
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
#ifdef ALIASED
#ifdef ALIASED2
   if ((jdim-2+1).ne.jl) then
#else
   if ((mw-mod(mw,2)+1).ne.jl) then
#endif
     print*,'jl.ne.nw',jl,mw-mod(mw,2)+1
     call exit(1)
   endif
#endif
!
! calculate coslat
!
   call dfs_sincos_lat
#else
   call sph_gaussian_lat (jdimhf, colrad, wgt, wgtcs, rcs2)
   call sph_poly_epsilon1(eps,jcap_)
   call sph_comp_index
#endif
!
!  define model sigma levels
!
#ifndef HYBRID
   call chgr_new_sigma(ci,si,del,sl,cl,rpi)
#else
   call chgr_hybrid_setup(ak5,bk5,ci,si,del,sl,cl,rpi)
   do k = 1,kdimp
     print 1002, k,ak5(k),bk5(k)
   enddo
   1002 format (1x,'ak bk',i4,2f15.7)
#endif
!
! define constants for saturation vapor pressure
!
   call funct_svp_init
!
!  read model surface topography (grid space).
!
   read(iumtn) zso
!
#ifdef SFC_CORRECTION
!
!  read sfc temperature  (grid space).
!
   read(iusfc) 
   read(iusfc) 
   read(iusfc) sfcto
!
   call dyn_trans2model_grid(sfcto,1)
#endif
!
#ifdef DFS
#define ISO_THERMAL
#ifndef ISO_THERMAL
   call dfs_standard_temp(iuatm,sl,tave,kdim,1000.)
#else
   tave(:)=300.
#endif
   print*,'tave=',tave
#endif
!
   in=0
   in2=0
!
! interpolate/extrapolate from slin to sl
!
   print *,' this program converts from pres to sl '
   print *,' pres=',(stdprs(k),k=1,kdimprs)
!
!  zero out spectral array for accumulation in lat loop.
!
   do i = 1,mdim
     gzo(i)=0.0
     qo(i)=0.0
   enddo
!
   do k = 1,kdim
     do i = 1,mdim
       teo(i,k)=0.0
       dio(i,k)=0.0
       zeo(i,k)=0.0
       rqo(i,k)=0.0
     enddo
   enddo
!
   allocate (out(idimt,jdimhf,kdimprs,nfld))
   allocate (out2(idimt,kdimprs,nfld))
   allocate (slin(kdimprs))
   allocate (siin(kdimprs+1))
!
   do k = 1,kdimprs
     slin(k)=stdprs(k)/stdprs(1)
   enddo
!
   siin(1)=1.D0
   do k = 2,kdimprs+1
     siin(k)=(stdprs(k-1)+stdprs(k))/2.D0/stdprs(1)
   enddo
!
   print *,' slin=',slin
   print *,' siin=',siin
   print *,' sl=  ',sl
!
!  fsig(1,1)=zg
!  fsig(1,2)=ps
!  fsig(1,3)=t
!  fsig(1,3+kdim)=q
!  fsig(1,3+kdim*2)=u
!  fsig(1,3+kdim*3)=v
!
   call chgrp_read_pressuregrib(iugrb,fn,iy,im,id,ih,fh,                       &
                                stdprs,kdimprs,kdimq,nfld,                     &
                                iupd5,ivpd5,itpd5,iqpd5,izpd5,                 &
                                out,idim,jdim)
   if(ifmtngrb.ne.0) then
     call chgrp_read_pressuregrib_2d(iugrb,fnmtngrb,iy,im,id,ih,fh,1,1,izpd5,  &
                                     zso,idim,jdim,iret)
     if(iret.eq.0.and.ifheight.ne.1) then
       do ij = 1,idim*jdim
         zso(ij,1)=zso(ij,1)/g_
       enddo
#ifdef DFS
       call dfs_fft_driver( 1,zso,idim,jdim,1,gzo,mt,jl,1,jl,                  &
                        1,1,0,coslat,1)
       call dfs_fft_driver(-1,zso,idim,jdim,1,gzo,mt,jl,1,jl,                  &
                        1,1,0,coslat,1)
#else
       !
       ! grid to wave
       !
       call dyn_trans2model_grid(zso,1)
       do lat = 1,jdimhf
         call chgr_sph_funct(qlnt,qlnv,colrad,lat)
         do i = 1,mdim
           qdert(i) = qlnt(i) * wgt(lat)
         enddo
         call sph_fft2grid(zso(1,lat),zss,2,-1)
         call chgr_fourier_diff(zss,zsa,1)
         call sph_coeff_lat(zss,zsa,gzo,qdert,1)
         !
         ! wave to grid
         !
         call chgr_sum_coeff(gzo,zso(1,lat),qlnt,     1)
         call sph_fft2grid(zso(1,lat),zso(1,lat),2,1)
       enddo
       call dyn_trans2output_grid(zso,1)
#endif
     endif
   endif
!
#ifdef BIN_DBG
   call file_write_bin(101,zso,idim,jdim,1,0)
   call file_write_bin(101,out,idim,jdim,kdimprs*nfld,0)
#endif
!
   if(ifheight.ne.1) then
     do k = 1,kdimprs
       do ij = 1,idim*jdim
         out(ij,1,k,5) = out(ij,1,k,5)/g_   ! ECMWD data z*g_
       enddo
     enddo
   endif
!
   call dyn_trans2model_grid(out,nfld*kdimprs)
   call dyn_trans2model_grid(zso,1)
!
!  at this point grid, values are available
!      out(1,1,1,1) ... wind x-direction
!      out(1,1,1,2) ... wind y-direction
!      out(1,1,1,3) ... temperature
!      out(1,1,1,4) ... moisture
!      out(1,1,1,5) ... geopotential height
!
!  note that idimt=idim and jdimhf=jdim for RMP
!
   do i = 1,mdim
     gzo(i)=0.0
   enddo
!
   dlat=pi_/jdim
   do lat = 1,jdimhf !!! latitude loop
!
#ifndef DFS
     do i = 1,idim
       xlat(i)=0.5*pi_-colrad(lat)
       xlat(i+idim)=-0.5*pi_+colrad(lat)
     enddo
     call chgr_sph_funct(qlnt,qlnv,colrad,lat)
     cslat=cos(xlat(1))
#else
     do i = 1,idim
       xlat(i)=0.5*pi_-(lat-0.5)*dlat
       xlat(i+idim)=-0.5*pi_+(lat-0.5)*dlat
     enddo
     cslat=coslat(lat,1)
#endif
!
!  compute surface pressure from model topography
!
!  ps=1000mb=100cb (stdprs(1))
!
     do i = 1,idimt
       ps(i)=log(stdprs(1)/10.)
       zs(i)=out(i,lat,1,5)
       fsig(i,1)=zso(i,lat)
     enddo
!
!  u,v multiplied by coslat
!
     do i = 1,idimt
       do k = 1,kdimprs
         out(i,lat,k,1)=out(i,lat,k,1)*cslat
         out(i,lat,k,2)=out(i,lat,k,2)*cslat
       enddo
     enddo
!
     do n = 1,nfld
       do k = 1,kdimprs
         do i = 1,idimt
           out2(i,k,n)=out(i,lat,k,n)
         enddo
       enddo
     enddo
!
!  convert rh to q for newps
!
     if(ifmoisrh.eq.1) then
       do k = 1,kdimprs
         do i = 1,idimt
           es=fpvs0(out2(i,k,3)) !! cb
           qs=xeps*es/(stdprs(k)*0.1+xepsm1*es)
           out2(i,k,4)=max(qmin_,out2(i,k,4)*0.01*qs)
         enddo
       enddo
     endif
!
!  compute surface pressure for new topography
!  chgrp_sfc_pressure: uses all Z to estimate sfc pressure. (better)
!                                !! unit of slin : [ND]
!    
     call chgrp_sfc_pressure(out2(1,1,3),out2(1,1,4),                          &
          kdimprs,ps,out2(1,1,5),fsig(1,1),fsig(1,2),siin,slin)
#ifdef DBG
     print *,'old new zs ps ',out2(1,1,5),fsig(1,1),exp(fsig(1,2))*10.
     print *,'t :(out2(1,k,3),k=1,kdimprs)=',                                  &
          (out2(1,k,3),k=1,kdimprs)
     print *,'r :(out2(1,k,4),k=1,kdimprs)=',                                  &
          (out2(1,k,4),k=1,kdimprs)
#endif
!
!! unit of ps and fsig(i,2) : ln(cb)
!! unit of slin : [ND]
!
! re-calculate sigma values of output using new sfc pressure (fsig(i,2))
!
#ifndef HYBRID
     do k = 1,kdimp
       do i = 1,idimt
         spino(i,k)=si(k)
       enddo
     enddo
     do k = 1,kdim
       do i = 1,idimt
         splno(i,k)=sl(k)
       enddo
     enddo
#else
     do i = 1, idimt
!      exppso(i)=exp(pso(i))
       exppso(i)=exp(fsig(i,2))
     enddo
!
     do k = 1,kdimp
       do i = 1,idimt
         spino(i,k) = ak5(k)/exppso(i) + bk5(k)
       enddo
     enddo
     do k = 1,kdim
       do i = 1,idimt
         dif = spino(i,k)**rk1 - spino(i,k+1)**rk1
         dif = dif / (rk1*(spino(i,k)-spino(i,k+1)))
         splno(i,k) = dif**rkinv
!        splno(i,k) = (spino(i,k)+spino(i,k+1))*0.5
       enddo
     enddo
#endif
!
!   extrapolate q above kdimqi when kdimqi < kdimprs
!
     if (kdimq.ne.kdimprs) then
       call chgrp_extrap_moisture(stdprs,out2(1,1,4),out2(1,1,3),              &
                    out2(1,1,4),idimt,kdimprs,kdimq,xlat)
     endif
#ifdef DBG
     print *,'after chgrp_extrap_moistures '
     print *,'t :(out2(1,k,3),k=1,kdimprs)=',                                  &
          (out2(1,k,3),k=1,kdimprs)
     print *,'r :(out2(1,k,4),k=1,kdimprs)=',                                  &
          (out2(1,k,4),k=1,kdimprs)
#endif
!
!  convert q to rh before vertical interpolation
!
     do k = 1,kdimprs
       do i = 1,idimt
         es=fpvs0(out2(i,k,3)) !! cb
         qs=xeps*es/(stdprs(k)/10.+xepsm1*es)
         out2(i,k,4)=min(max(out2(i,k,4)/qs,0.),1.)*100.
       enddo
     enddo
#ifdef DBG
     print *,'after rh  '
     print *,'r :(out2(1,k,4),k=1,kdimprs)=',                                  &
          (out2(1,k,4),k=1,kdimprs)
#endif
!
!  vertical interpolation in the order of u, v, t, q
!  chgrp_sigma2sigma sets values to constant outside of input domain.
!  fix lapse rate and relative humidity below input surface.
!
     if(ifcubic.ne.1) then 
       call chgrp_pressure2sigma(fsig(1,2),stdprs,splno,out2(1,1,1),           &
                  fsig(1,2*kdim+3),                                            &
                  idimt,kdimprs,kdim,2)
       call chgrp_pressure2sigma(fsig(1,2),stdprs,splno,out2(1,1,2),           &
                  fsig(1,3*kdim+3),                                            &
                  idimt,kdimprs,kdim,3)
       call chgrp_pressure2sigma(fsig(1,2),stdprs,splno,out2(1,1,3),           &
                  fsig(1,       3),                                            &
                  idimt,kdimprs,kdim,4)
       call chgrp_pressure2sigma(fsig(1,2),stdprs,splno,out2(1,1,4),           &
                  fsig(1,  kdim+3),                                            &
                  idimt,kdimprs,kdim,5)
     else
       call chgrp_sigma2sigma(ps,slin,out2(1,1,1),fsig(1,2),splno,             &
                  fsig(1,2*kdim+3),                                            &
                  idimt,1,kdimprs,kdim,in,tensn,1)
       call chgrp_sigma2sigma(ps,slin,out2(1,1,2),fsig(1,2),splno,             &
                  fsig(1,3*kdim+3),                                            &
                  idimt,1,kdimprs,kdim,in,tensn,1)
       call chgrp_sigma2sigma(ps,slin,out2(1,1,3),fsig(1,2),splno,             &
                  fsig(1,       3),                                            &
                  idimt,1,kdimprs,kdim,in,tensn,2)
       call chgrp_sigma2sigma(ps,slin,out2(1,1,4),fsig(1,2),splno,             &
                  fsig(1,  kdim+3),                                            &
                  idimt,1,kdimprs,kdim,in,tensn,3)
     endif
!
#ifdef DBG
!
!        print *,'t :(out2(i,10,3),i=1,10)=',(out2(i,10,3),i=1,10)
!        print *,'rh:(out2(i,10,4),i=1,10)=',(out2(i,10,4),i=1,10)
!        print *,'u :(out2(i,10,1),i=1,10)=',(out2(i,10,1),i=1,10)
!        print *,'v :(out2(i,10,2),i=1,10)=',(out2(i,10,2),i=1,10)
!
     print *,'t :(out2(1,k,3),k=1,kdimprs)=',                                  &
          (out2(1,k,3),k=1,kdimprs)
     print *,'t :(fsig(1,2+k),k=1,kdim)=',                                     &
          (fsig(1,2+k),k=1,kdim)
     print *,'rh:(out2(1,k,4),k=1,kdimprs)=',                                  &
          (out2(1,k,4),k=1,kdimprs)
     print *,'rh:(fsig(1,2+kdim+k),k=1,kdim)=',                                &
          (fsig(1,2+kdim+k),k=1,kdim)
     print *,'u :(out2(1,k,1),k=1,kdimprs)=',                                  &
          (out2(1,k,1),k=1,kdimprs)
     print *,'u :(fsig(1,2+2*kdim+k),k=1,kdim)=',                              &
          (fsig(1,2+2*kdim+k),k=1,kdim)
     print *,'v :(out2(1,k,2),k=1,kdimprs)=',                                  &
          (out2(1,k,2),k=1,kdimprs)
     print *,'v :(fsig(1,2+3*kdim+k),k=1,kdim)=',                              &
          (fsig(1,2+3*kdim+k),k=1,kdim)
#endif
!
! --- Extrapolate the T and rh below 1000 hPa
!
     do i = 1,idimt
       ps1=log(slin(1))+ps(i)
       do k = 1,kdim
         psk=log(splno(i,k))+fsig(i,2)
         if(psk.ge.ps1) then !! below 1000hPa level
           fsig(i,k+2)=out2(i,1,3)*exp(gamma*rd_/g_*(psk-ps1)) !! t
           fsig(i,kdim+k+2)=out2(i,1,4) !! rh
#ifdef DBG
           if(i.eq.1) then
             print*,' osk ps1 ',psk,ps1
             print *,' t-correction  :',k,out2(i,k,3),fsig(i,2+k)
           endif
#endif
         endif
       enddo
     enddo
#ifdef SFC_CORRECTION
     do i = 1,idimt
       ps1=log(slin(1))+ps(i)
       do k = 1,kdim
         psk=log(splno(i,k))+fsig(i,2)
         if(psk.ge.ps1) then !! below 1000hPa level
#ifdef DBG
           if(i.eq.1) then
             print*,' before : weighting ',wei
             print *,' t-correction  :',k,out2(i,1,3),fsig(i,2+k)
             print *,' u-correction  :',k,out2(i,1,1),fsig(i,kdim*2+2+k)
             print *,' v-correction  :',k,out2(i,1,2),fsig(i,kdim*3+2+k)
           endif
#endif
           wei = (psk-ps1)/(fsig(i,2)-ps1)
           fsig(i,k+2)=out2(i,1,3)+wei*(sfcto(i,lat)-out2(i,1,3))    !! T
           fsig(i,kdim+k+2)=out2(i,1,4) !! rh
           fsig(i,kdim*2+k+2)=out2(i,1,1)+wei*(-out2(i,1,1)) !! u
           fsig(i,kdim*3+k+2)=out2(i,1,2)+wei*(-out2(i,1,2)) !! v
#ifdef DBG
           if(i.eq.1) then
             print*,' after : weighting ',wei
             print *,' t-correction  :',k,fsig(i,2+k)
             print *,' u-correction  :',k,fsig(i,kdim*2+2+k)
             print *,' v-correction  :',k,fsig(i,kdim*3+2+k)
           endif
#endif
         endif
       enddo
     enddo
#endif
!
!  convert rh to q and temp to virtual temp
!
     do k = 1,kdim
       do i = 1,idimt
         es=fpvs0(fsig(i,k+2)) !! es is in cb
         qs=xeps*es/((exp(fsig(i,2))*sl(k))+xepsm1*es)
         fsig(i,kdim+k+2)=max(qmin_,fsig(i,kdim+k+2)*0.01*qs)
         fsig(i,k+2)=fsig(i,k+2)*(1+0.61*fsig(i,kdim+k+2))
       enddo
     enddo
!
     do k = 1,kdimprs
       do i = 1,idimt
         es=fpvs0(out2(i,k,3)) !! es is in cb
         qs=xeps*es/(stdprs(k)/10.+xepsm1*es)
         out2(i,k,4)=max(qmin_,out2(i,k,4)*0.01*qs)
         out2(i,k,3)=out2(i,k,3)*(1+0.61*out2(i,k,4))
       enddo
     enddo
#ifdef DBG
     print *,'Tv:(out2(1,k,3),k=1,kdimprs)=',                                  &
          (out2(1,k,3),k=1,kdimprs)
     print *,'Tv:(fsig(1,2+k),k=1,kdim)=',                                     &
          (fsig(1,2+k),k=1,kdim)
     print *,'q :(out2(1,k,4),k=1,kdimprs)=',                                  &
          (out2(1,k,4),k=1,kdimprs)
     print *,'q :(fsig(1,2+kdim+k),k=1,kdim)=',                                &
          (fsig(1,2+kdim+k),k=1,kdim)
#endif
!
#if defined(BIN_DBG) || defined(DBG)
     do i = 1,idimt
       zout(i,lat)=fsig(i,1)
       pout(i,lat)=fsig(i,2)
     enddo
!
     do k = 1,kdim
       do i = 1,idimt
         tout(i,lat,k)=fsig(i,2+kdim*0+k)
         qout(i,lat,k)=fsig(i,2+kdim*1+k)
         uout(i,lat,k)=fsig(i,2+kdim*2+k)
         vout(i,lat,k)=fsig(i,2+kdim*3+k)
       enddo
     enddo
#endif
#ifdef DFS
     !
     ! save grid value
     ! reverse for dfs grid : south to north
     !
     !  fsig(1,1)=zg
     !  fsig(1,2)=ps
     !  fsig(1,3)=t
     !  fsig(1,3+kdim)=q
     !  fsig(1,3+kdim*2)=u
     !  fsig(1,3+kdim*3)=v
     zsa(1:idim,lat)=fsig(1:idim,1)
     zsa(1:idim,jdim-lat+1)=fsig(idim+1:2*idim,1)
     psa(1:idim,lat)=fsig(1:idim,2)
     psa(1:idim,jdim-lat+1)=fsig(idim+1:2*idim,2)            ! kPa
     do k = 1,kdim
       uua(1:idim,lat,k)=fsig(1:idim,k+2+kdim*2)
       uua(1:idim,jdim-lat+1,k)=fsig(idim+1:2*idim,k+2+kdim*2)
       !
       vva(1:idim,lat,k)=fsig(1:idim,k+2+kdim*3)
       vva(1:idim,jdim-lat+1,k)=fsig(idim+1:2*idim,k+2+kdim*3)
       !
       tta(1:idim,lat,k)=fsig(1:idim,k+2)-tave(k)
       tta(1:idim,jdim-lat+1,k)=fsig(idim+1:2*idim,k+2)-tave(k)
       !
       qqa(1:idim,lat,k)=fsig(1:idim,k+2+kdim)
       qqa(1:idim,jdim-lat+1,k)=fsig(idim+1:2*idim,k+2+kdim)
     enddo
   enddo                !!! end of lat loop
!
#ifdef BIN_DBG
   call file_write_bin(112,zsa,idim,jdim,1,0)
   call file_write_bin(112,psa,idim,jdim,1,0)
   call file_write_bin(112,uua,idim,jdim,kdim,0)
   call file_write_bin(112,vva,idim,jdim,kdim,0)
   call file_write_bin(112,tta,idim,jdim,kdim,0)
   call file_write_bin(112,qqa,idim,jdim,kdim,0)
#endif
!
!  dfs transform
!
   call dfs_fft_driver( 1,zsa,idim,jdim,1,gzo,mt,jl,1,jl,                      &
                           1,1,0,coslat,1)
   call dfs_fft_driver(1,psa,idim,jdim,1,qo ,mt,jl,1,jl,                       &
                        1,1,0,coslat,1)
   call dfs_fft_driver(1,tta,idim,jdim,kdim,teo,mt,jl,kdim,jl,                 &
                        kdim,kdim,1,coslat,1)
   call dfs_fft_driver(1,qqa,idim,jdim,kdim,rqo,mt,jl,kdim,jl,                 &
                        kdim,kdim,1,coslat,1)
   call dfs_wind2divor(uua,vva,idim,jdim,kdim,zeo,dio,mt,jl,kdim,coslat,1)
#ifdef ALIASED
   call dfs_cut_alias(qo,mt,0,jla,1   )
   call dfs_cut_alias(dio,mt,0,jla,kdim)
   call dfs_cut_alias(teo,mt,0,jla,kdim)
   call dfs_cut_alias(zeo,mt,0,jla,kdim)
   call dfs_cut_alias(rqo,mt,0,jla,kdim)
   call dfs_cut_alias(gzo,mt,0,jla,1   )
#endif
#ifdef BIN_DBG
   call dfs_fft_driver(-1,uua,idim,jdim,kdim,dio,mt,jl,kdim,jl,                &
                        kdim,kdim,1,coslat,1)
   call dfs_fft_driver(-1,vva,idim,jdim,kdim,zeo,mt,jl,kdim,jl,                &
                        kdim,kdim,1,coslat,1)
!
!  div
!
   call file_write_bin(122,dio,mt,jl,kdim,0)
!
!  vor
!
   call file_write_bin(122,zeo,mt,jl,kdim,0)
   call file_write_bin(122,teo,mt,jl,kdim,0)
   call file_write_bin(122,qo,mt,jl,kdim,0)
#endif
#else /* else of DFS */
!
!  gsm output
!  spherical transform
!
     call sph_fft2grid(fsig(1,1),zss,2,-1)
     call sph_fft2grid(fsig(1,2),pss,2,-1)
     call sph_fft2grid(fsig(1,3),tts,2*kdim,-1)
     call sph_fft2grid(fsig(1,3+kdim),qqs,2*kdim,-1)
     call sph_fft2grid(fsig(1,3+2*kdim),uus,2*kdim,-1)
     call sph_fft2grid(fsig(1,3+3*kdim),vvs,2*kdim,-1)
!
     do i = 1,mdim
       qdert(i) = qlnt(i) * wgt(lat)
     enddo
!
     call chgr_fourier_diff(zss,zsa,1)
     call sph_coeff_lat(zss,zsa,gzo,qdert,1)
!
     call chgr_fourier_diff(pss,psa,1)
     call sph_coeff_lat(pss,psa,qo ,qdert,1)
!
     call chgr_fourier_diff(tts,tta,kdim)
     call sph_coeff_lat(tts,tta,teo,qdert,kdim)
!
     call chgr_fourier_diff(qqs,qqa,kdim)
     call sph_coeff_lat(qqs,qqa,rqo,qdert,kdim)
!
     call chgr_sph_derivative(qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgt(lat))
!
     call chgr_fourier_diff(uus,uua,kdim)
     call chgr_fourier_diff(vvs,vva,kdim)
     call chgr_sph_operator(uua,uus,vva,vvs,zeo,qlnwct,qdert,kdim)
     call chgr_div_coeff_lat(uua,uus,vva,vvs,dio,qlnwct,qdert,kdim)
!
   enddo !                !!! end of lat loop
!
   do k = 1,kdim
     do m = 1,mdim
       dio(m,k)=-dio(m,k)
     enddo
   enddo
#endif /* DFS end */
!
#ifdef DBG
   call print_maxmin_six(zout,idim*jdim,1,1,1,'zout')
   call print_maxmin_six(pout,idim*jdim,1,1,1,'pout')
   call print_maxmin_six(tout,idim*jdim,kdim,1,kdim,'tout')
   call print_maxmin_six(qout,idim*jdim,kdim,1,kdim,'qout')
   call print_maxmin_six(uout,idim*jdim,kdim,1,kdim,'uout')
   call print_maxmin_six(vout,idim*jdim,kdim,1,kdim,'vout')
#endif
#ifdef BIN_DBG
   call file_write_bin(112,zout,idim,jdim,1,1)
   call file_write_bin(112,pout,idim,jdim,1,1)
   call file_write_bin(112,uout,idim,jdim,kdim,1)
   call file_write_bin(112,vout,idim,jdim,kdim,1)
   call file_write_bin(112,tout,idim,jdim,kdim,1)
   call file_write_bin(112,qout,idim,jdim,kdim,1)
#endif
!
   idate(4)=iy
   idate(2)=im
   idate(3)=id
   idate(1)=ih
   fhour=fh
   lab='        '
!
   do k = 1,kdum
     dummy(k)=0.
   enddo
!
   waves=jcap_
   xlayers=levs_
   trun=1.
   order=2.
   realform=1.
   gencode=7.
   rlond=lonf_
   rlatd=latg_
   rlonp=lonf_
   rlatp=latg_
   rlonr=lonf_
   rlatr=latg_
   gases=0.
   water=1.
   pdryini=1. !pdryini=0 in chgr_write_sigma
   subcen=0.
   ensemble(1)=0.
   ensemble(2)=0.
   ppid=0.
   slid=0.
   vcid=0.
   vmid=0.
   vtid=0.
   runid=0.
   usrid=0.
!
   do k = 1,kdum2
     dummy2(k)=0.
   enddo
!
   write(iusig) lab
#ifndef HYBRID
   write(iusig) fhour,idate,si,sl                                              &
#else
   do k = 1,levs_+1
     ak5x(k)=ak5(k)*1000. ! cb -> Pa
     bk5x(k)=bk5(k)
   enddo
!
   write(iusig)fhour,idate,ak5x,bk5x                                           &
#endif
          ,dummy,waves,xlayers,trun,order,realform,gencode                     &
          ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                           &
          ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                &
          ,pdryini,dummy2,gases  !same as chgr
!          ,dummy2,gases
#ifdef DFS
   write(iusig) gzo(:)
   write(iusig) qo(:),tave(:)
   call print_maxmin_six(teo,mdim,kdim,1,kdim,'temp in chgr_write_sigma')
!
   do k = 1,kdim
     write(iusig) teo(:,k)
   enddo
!
   do k = 1,kdim
     write(iusig) dio(:,k)
     write(iusig) zeo(:,k)
   enddo
!
   call print_maxmin_six(rqo,mdim,kdim,1,kdim,'rq in chgr_write_sigma')
   do k = 1,kdim
     write(iusig) rqo(:,k)
   enddo
#else
!
   do i = 1,mdim
     b(indxmm(i)) = gzo(i)
   enddo
   write(iusig) b
   do i = 1,mdim
     b(indxmm(i)) = qo(i)
   enddo
   write(iusig) b
   do k = 1,kdim
     do i = 1,mdim
       b(indxmm(i)) = teo(i,k)
     enddo
     write(iusig) b
   enddo
   do k = 1,kdim
     do i = 1,mdim
       b(indxmm(i)) = dio(i,k)
     enddo
     write(iusig) b
     do i = 1,mdim
       b(indxmm(i)) = zeo(i,k)
     enddo
     write(iusig) b
   enddo
   do k = 1,kdim
     do i = 1,mdim
       b(indxmm(i)) = rqo(i,k)
     enddo
     write(iusig) b
   enddo
#endif
!
   deallocate(xlat,ps,zs,expps)
   deallocate(fsig, qe, dummy)
#ifdef DFS
   deallocate(psa,tta,qqa,zsa,grid,tave)
#else
   deallocate(qlnt,qlnv,qdert,qlnwct,b)
#endif
   deallocate(gzso,zso)
#ifdef SFC_CORRECTION
   deallocate(sfcti,sfcto)
#endif
   deallocate(si,sl,del,ci,cl,rpi)
#ifdef HYBRID
   deallocate(ak5,bk5,ak5x,bk5x)
#endif
#if defined(DBG) || defined(BIN_DBG)
   deallocate(zout,pout,tout,qout,uout,vout)
#endif
   deallocate (out,out2,slin,siin)
#ifdef HYBRID
   deallocate(exppso)
#endif
   deallocate(spino,splno)
!
   stop
   end program change_resolution_pressure
