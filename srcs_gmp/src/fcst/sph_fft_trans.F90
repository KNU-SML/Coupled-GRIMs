#include <define.h>
   subroutine sph_fft_trans
!-------------------------------------------------------------------------------
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : latg_,lonf_
#ifdef REDUCE_GRID
   use comreduce
#define LATS latg_/2
#else
#define LATS 1
#endif
!
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifndef RMP
   use comfcst, only : crscale,rcscale
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
   use comfcst, only : scale,trig
#endif
#ifdef FFTW
#undef DEFAULT
   use comfcst, only : scale, iplan_c_to_r,iplan_r_to_c
#endif
#ifdef DEFAULT
   use comfcst, only : ifax, trigs
#endif
!-------------------------------------------------------------------------------
#include <abort.h>
!-------------------------------------------------------------------------------
!                                                                               
   do lat = 1,LATS
#ifdef REDUCE_GRID
     lonff=lonfd(lat)
#else
     lonff=lonf_
#endif
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
     crscale(lat)=1.0
     rcscale(lat)=1./float(lonff)
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
     scale(lat)=1./float(lonff)
     call dzfftm1dui (lonff,trigs(1,lat))                                      
#endif
#ifdef FFTW
#undef DEFAULT
     ifftw_real_to_complex=-1
     ifftw_complex_to_real=1
     ifftw_estimate=0
     scale(lat)=1./float(lonff)
     call rfftw_f77_create_plan(iplan_c_to_r(lat),lonff,                       &
                        ifftw_complex_to_real,ifftw_estimate)
     call rfftw_f77_create_plan(iplan_r_to_c(lat),lonff,                       &
                        ifftw_real_to_complex,ifftw_estimate)
!
! ---- plan can be destryed as following, but may not be necessary ----
! ---  call rfftw_f77_destroy_plan(int_plan)
!
#endif
#ifdef RFFTMLT
#undef DEFAULT
     call fftfax (lonff,ifax(1,lat),trigs(1,1,lat))                                      
#endif
#ifdef ASLES
#undef DEFAULT
     call ldfrmfb(lonff,0,0,0,0,0,ifax(1,lat),                                 &
                                   trigs(1,1,lat),0,ierr)
#endif
#ifdef DEFAULT
     call    fax (ifax(1,lat), lonff,3)                                         
     call fftrig (trigs(1,1,lat),lonff,3)                                         
#endif
#ifdef ASLES
     if (ierr .ge. 3000)  then
       write(6,*)' error in asles fft initialization'
       call MPABORT
     endif
#else
#ifdef DEFAULT
     if (ifax(1,lat) .eq. -99)  then
       write(6,'(A,I2,A)')'error in sph_fft_trans. lonff=',lonff,' not factorable.'
       call MPABORT
     endif
#endif
#endif
   enddo
!
   return                                                                    
   end subroutine sph_fft_trans
#ifndef RMP
!-------------------------------------------------------------------------------
   subroutine sph_fft_driver(isgn,grid,wave,ntotal)
!-------------------------------------------------------------------------------
! abstract: Spherical harmonics transform on the sphere
!
!     isgn   =  1 : grid to wave
!     isgn   = -1 : wave to grid
!-------------------------------------------------------------------------------
   use paramodel, only : LCAP22S,LNT22S,JCAP1S,levs_,latg2_,jcap_,jcap1_      ,&
                         LONF2S,LATG2S
   use comio, only     : iope
#ifdef MP
   use paramodel, only : levsp_,lln22p_,lnt22p_,lcap22p_              ,&
                         lonf2_,lonf22_,lonf22p_,lonf2p_,latg2p_
   use commpi, only    : latdef, lwvdef, lwvstr, lwvlen, mype         ,&
                         latstr, latlen
#else
   use comfgrid, only : latdef,lwvdef
#endif
#ifdef REDUCE_GRID
   use comreduce
#else
   use paramodel, only : lonf_
#endif
   use comgpln , only : qtt,qww
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                 ::  isgn,ntotal
   integer                                 ::  nlevs,nlevsp
   integer                                 ::  lat, lat1, lat2, latdon        ,&
                                               lcapf, lonff,llstr,llens,llensd,&
                                               lan,i,k,jstr,jend
   real, dimension(LONF2S,levs_*ntotal,LATG2S) ::  grid
   real, dimension(LNT22S,levs_*ntotal)        ::  wave
   real, allocatable, dimension(:,:,:)     ::  anl
   real, allocatable, dimension(:,:,:,:)   ::  flp,flm
#ifdef MP
   real, allocatable, dimension(:,:)       ::  pwave
   real, allocatable, dimension(:,:,:)     ::  anf,grs
#endif
!
! setting variables
!
   nlevs =levs_*ntotal
#ifdef MP
   nlevsp=levsp_*ntotal
#else
   nlevsp=nlevs
#endif
!
! allocation
!
   allocate(anl(LCAP22S ,nlevsp,latg2_)   ,& 
            flp(2,JCAP1S,nlevsp,latg2_)   ,&
            flm(2,JCAP1S,nlevsp,latg2_)     )
#ifdef MP
   allocate(pwave(lln22p_,nlevsp)           )
   allocate(anf(lonf22_ ,nlevsp ,latg2p_) ,&
            grs(lonf22p_,nlevs  ,latg2p_)   )
#endif
!
! initialize
!
   anl=0.   ;  flp=0. ;  flm=0.
#ifdef MP
   pwave=0. ;  anf=0. ;  grs=0.
#endif
!
! Latitude bands
!
   lat1=1
   lat2=latg2_
   latdon=0
!
! index for wave and grid
!
#ifdef MP
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
#define ANLS anf
#define WAVE pwave
#else
   llstr=0
   llens=jcap1_
#define ANLS anl
#define WAVE wave
#endif
!-------------------------------------------------------------------------------
   if(isgn.eq.1) then   ! grid to wave
!-------------------------------------------------------------------------------
     wave=0.
!
     do lat = 1,LATG2S
       do k = 1,nlevs
         do i = 1,LONF2S
#ifdef MP
           grs(i,k,lat)=grid(i,k,lat)
#else
           anl(i,k,lat)=grid(i,k,lat)
#endif
         enddo
       enddo
     enddo
#ifdef MP
!
! mpnk2nx
!
     call mpnx2nk(grs,lonf22p_,nlevs,anf,lonf22_,nlevsp,latg2p_,               &
                  levs_,levsp_,1,1,ntotal)
!
! Latitude bands for MPI
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#endif
!
! grid to wave
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(ANLS(1,1,lan),ANLS(1,1,lan),2*nlevsp,lcapf,          &
                          lonff,latdef(lat),-1)
     enddo
#ifdef MP
!
! mpnl2ny
!
     call mpny2nl(ANLS,lonf22_,latg2p_,anl,lcap22p_,latg2_,nlevsp,1,nlevsp)
!
! Latitude bands
!
     lat1=1
     lat2=latg2_
     latdon=0
#endif
!
! grid to wave (2)
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_grid2wave(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),llensd,   &
                          nlevsp)
     enddo
!
! sph_comp_coeff2
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_comp_coeff2(flp(1,1,1,lan),flm(1,1,1,lan),WAVE,qww(1,lat),     &
                              llstr,llensd,lwvdef,nlevsp)
     enddo
#ifdef MP
!
! 2D decomposition
!
     call mpnk2nn(pwave,lln22p_,levsp_,wave,lnt22p_,levs_,ntotal)
#endif
!-------------------------------------------------------------------------------
   else if(isgn.eq.-1) then   ! wave to grid
!-------------------------------------------------------------------------------
     grid=0.
#ifdef MP
!
! 2D decomposition
!
     call mpnn2nk(wave,lnt22p_,levs_,pwave,lln22p_,levsp_,ntotal)
#endif
!
! sph_sum_coeff
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_sum_coeff(WAVE,anl(1,1,lan),qtt(1,lat),llstr,llensd,lwvdef,nlevsp)
     enddo
#ifdef MP
!
! mpnl2ny
!
     call mpnl2ny(anl,lcap22p_,latg2_,anf,lonf22_,latg2p_,nlevsp,1,nlevsp)
!
! Latitude bands for MPI
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#endif
!
! wave to grid
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(ANLS(1,1,lan),ANLS(1,1,lan),2*nlevsp,lcapf,lonff,latdef(lat),1)
     enddo
#ifdef MP
!
! mpnk2nx
!
     call mpnk2nx(anf,lonf22_,nlevsp,grs,lonf22p_,nlevs,latg2p_,levsp_,levs_,1,1,ntotal)
#endif
     do lat = 1,LATG2S
       do k=1,nlevs
         do i = 1,LONF2S
#ifdef MP
           grid(i,k,lat)=grs(i,k,lat)
#else
           grid(i,k,lat)=anl(i,k,lat)
#endif
         enddo
       enddo
     enddo
!-------------------------------------------------------------------------------
   endif
!
#undef ANLS
#undef WAVE
   deallocate(anl,flp,flm)
#ifdef MP
   deallocate(pwave,anf,grs)
#endif
   end subroutine sph_fft_driver
!
!-------------------------------------------------------------------------------
   subroutine sph_fft_driver_2d(isgn,grid,wave,ntotal)
!-------------------------------------------------------------------------------
! abstract: Spherical harmonics transform on the sphere
!
!     isgn   =  1 : grid to wave
!     isgn   = -1 : wave to grid
!-------------------------------------------------------------------------------
   use paramodel, only : LCAP22S,LNT22S,JCAP1S,levs_,latg2_,jcap_,jcap1_      ,&
                         LONF2S,LATG2S
   use comio, only     : iope
#ifdef MP
   use paramodel, only : levsp_,lln22p_,lnt22p_,lcap22p_              ,&
                         lonf2_,lonf22_,lonf22p_,lonf2p_,latg2p_
   use commpi, only    : latdef, lwvdef, lwvstr, lwvlen, mype         ,&
                         latstr, latlen
#else
   use comfgrid, only : latdef,lwvdef
#endif
#ifdef REDUCE_GRID
   use comreduce
#else
   use paramodel, only : lonf_
#endif
   use comgpln , only : qtt,qww
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                 ::  isgn,ntotal
   integer                                 ::  lat, lat1, lat2, latdon        ,&
                                               lcapf, lonff,llstr,llens,llensd,&
                                               lan,i,k,jstr,jend
   real, dimension(LONF2S,LATG2S)          ::  grid
   real, dimension(LNT22S)                 ::  wave
   real, allocatable, dimension(:,:)       ::  anl
   real, allocatable, dimension(:,:,:)     ::  flp,flm
#ifdef MP
   real, allocatable, dimension(:)         ::  pwave
   real, allocatable, dimension(:,:)       ::  anf,grs
#endif
!
! allocation
!
   allocate(anl(LCAP22S ,latg2_)   ,& 
            flp(2,JCAP1S,latg2_)   ,&
            flm(2,JCAP1S,latg2_)     )
#ifdef MP
   allocate(pwave(lln22p_)           )
   allocate(anf(lonf22_ ,latg2p_) ,&
            grs(lonf22p_,latg2p_)   )
#endif
!
! initialize
!
   anl=0.   ;  flp=0. ;  flm=0.
#ifdef MP
   pwave=0. ;  anf=0. ;  grs=0.
#endif
!
! Latitude bands
!
   lat1=1
   lat2=latg2_
   latdon=0
!
! index for wave and grid
!
#ifdef MP
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
#define ANLS anf
#define WAVE pwave
#else
   llstr=0
   llens=jcap1_
#define ANLS anl
#define WAVE wave
#endif
!-------------------------------------------------------------------------------
   if(isgn.eq.1) then   ! grid to wave
!-------------------------------------------------------------------------------
     wave=0.
!
     do lat = 1,LATG2S
       do i = 1,LONF2S
#ifdef MP
         grs(i,lat)=grid(i,lat)
#else
         anl(i,lat)=grid(i,lat)
#endif
       enddo
     enddo
#ifdef MP
!
! mpnk2x
!
     call mpnx2x(grs,lonf22p_,1,anf,lonf22_,1,latg2p_,1,1,1)
!
! Latitude bands for MPI
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#endif
!
! grid to wave
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(ANLS(1,lan),ANLS(1,lan),2,lcapf,lonff,latdef(lat),-1)
     enddo
#ifdef MP
!
! mpnl2ny
!
     call mpny2nl(ANLS,lonf22_,latg2p_,anl,lcap22p_,latg2_,1,1,1)
!
! Latitude bands
!
     lat1=1
     lat2=latg2_
     latdon=0
#endif
!
! grid to wave (2)
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_grid2wave(flp(1,1,lan),flm(1,1,lan),anl(1,lan),llensd,1)
     enddo
!
! sph_comp_coeff2
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_comp_coeff2(flp(1,1,lan),flm(1,1,lan),WAVE,qww(1,lat), &
                              llstr,llensd,lwvdef,1)
     enddo
#ifdef MP
!
! 2D decomposition
!
     call mpn2nn(pwave,lln22p_,wave,lnt22p_,1)
#endif
!-------------------------------------------------------------------------------
   else if(isgn.eq.-1) then   ! wave to grid
!-------------------------------------------------------------------------------
     grid=0.
#ifdef MP
!
! 2D decomposition
!
     call mpnn2n(wave,lnt22p_,pwave,lln22p_,1)
#endif
!
! sph_sum_coeff
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_sum_coeff(WAVE,anl(1,lan),qtt(1,lat),llstr,llensd,lwvdef,1)
     enddo
#ifdef MP
!
! mpnl2ny
!
     call mpnl2ny(anl,lcap22p_,latg2_,anf,lonf22_,latg2p_,1,1,1)
!
! Latitude bands for MPI
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#endif
!
! wave to grid
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(ANLS(1,lan),ANLS(1,lan),2,lcapf,lonff,latdef(lat),1)
     enddo
#ifdef MP
!
! mpnk2nx
!
     call mpx2nx(anf,lonf22_,1,grs,lonf22p_,1,latg2p_,1,1,1)
#endif
     do lat = 1,LATG2S
       do i = 1,LONF2S
#ifdef MP
         grid(i,lat)=grs(i,lat)
#else
         grid(i,lat)=anl(i,lat)
#endif
       enddo
     enddo
!-------------------------------------------------------------------------------
   endif
!
#undef ANLS
#undef WAVE
   deallocate(anl,flp,flm)
#ifdef MP
   deallocate(pwave,anf,grs)
#endif
   end subroutine sph_fft_driver_2d
#endif /* ~RMP end */
