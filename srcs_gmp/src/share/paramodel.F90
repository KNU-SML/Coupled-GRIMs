#include <define.h>
module paramodel
!-------------------------------------------------------------------------------
   use constant
#ifdef NIM
   use ReadNamelist
#endif
!-------------------------------------------------------------------------------
   integer      jcap_     ,&    ! number of spectral wave for gmp
                levs_     ,&    ! number of full vertical layer
                lonf_     ,&    ! number of i-grid  for gmp
                latg_     ,&    ! number of j-grid  for gmp
                nvect_    ,&    ! number of vector length (=64
                lonf2_    ,&    ! lonf * 2
                lonf22_   ,&    ! lonf * 2 + 2 
                latg2_    ,&    ! latg / 2 
!
                ncpus_    ,&    ! number of cpu when thread model is on
                jcap1_    ,&    ! jcap + 1
                jcap2_    ,&    ! jcap + 1
                lnt_      ,&
                lnt2_     ,&
                lnt22_    ,&
                lnut22_   ,&
                lnut_     ,&
                lnut2_    ,&
                lnuv_     ,&
                lnu2_     ,&
                twoj1_    ,&
!
                nwmass_   ,&    ! number of water specicies mass including water vapor
                ncloud_   ,&    ! number of hydrometeors
                icloud_   ,&    ! starting index for hydrometeors
                kcloud_   ,&    ! starting k- index for hydrometeors
                nwsize_   ,&    ! number of water specicies size distribution 
                iwsize_   ,&    ! starting index for water specicies size distribution 
                kwsize_   ,&    ! starting k- index for water specicies size distribution 
                ngases_   ,&    ! number of gases
                igases_   ,&    ! starting index for gases 
                kgases_   ,&    ! starting k- index for gases 
                ibltke_   ,&    ! starting index for tke 
                kbltke_   ,&    ! starting k- index for tke pbl  
                nwater_   ,&    ! number of water species 
                ntotal_   ,&    ! nwater_ + ngases_
!
                levh_     ,&    ! levs * ntotal_
                levp1_    ,&    ! levs_ + 1
                levp2_    ,&    ! levs_ + 2
                levm1_    ,&    ! levs_ - 1
                lotanu_     
   real         mtnres_   ,&    ! resolution of mountain data
                rmtnres_        ! resolution of mountain data for rmp
   integer      mtnvar_   ,&    ! number of mountain statistics file
                rmtnvar_  ,&    ! number of mountain statistics file for rmp
                imn_      ,&    ! i-grid number of mountain data
                jmn_      ,&    ! j-grid number of mountain data
                imnr_     ,&    ! i-grid number of mountain data for rmp
                jmnr_     ,&    ! j-grid number of mountain data for rmp
                slvark_   ,&    ! number of surface variable for time series output
                mlvark_   ,&    ! number of upper level data for time series output
                lpnt_     ,&    ! number of point output : currently not used
                ltstp_    ,&    ! number of point output : currently not used
                igen_     ,&    ! index for source code generation (not used)
!
                ilonf_    ,&    ! input grid resolution : i-grid number
                ilatg_    ,&    ! input grid resolution : j-grid number
                ijcap_    ,&    ! input grid resolution : wave number
                ilevs_    ,&    ! input grid resolution : k-grid number
!
                io_       ,&    ! output-pressure grid : i-grid number
                jo_       ,&    ! output-pressure grid : j-grid number
                iod_      ,&    ! io * 2 
                iodd_     ,&    ! io * 2 + 6
                siodd_    ,&    ! lonf_ * 2 + 6
                johf_     ,&    ! (jo+1) / 2
                ko_       ,&    ! output-pressure grid : k-grid number 
                kt_       ,&    ! output-pressure grid: tropopause information
                klot_     ,&    ! dimension of k-direction
                nstype_   ,&    ! number of soil type
                nvtype_   ,&    ! number of vegetation type
                ILOTS           ! dimension of i-direction
#if defined(SMP_RA2SFC)
   integer      jcapscm_  ,&    ! wave number for scm
                lonfscm_  ,&    ! i-grid number for scm
                latgscm_  ,&    ! j-grid number for scm
                lonf2scm_ ,&    ! lonfscm * 2
                latg2scm_       ! latgscm / 2
#endif
#ifdef RIVER
   integer      io2_      ,&    ! i-grid number for river model
                jo2_            ! j-grid number for river model
   real         delxo2_   ,&    ! river model resolution in degree
                delyo2_         ! river model resolution in degree
#endif
#ifdef MUL_CLDTOP
   integer      ncldtop_        ! number of cloud top in sas multi clouds (=3)
#endif
#ifdef GMPDAMP
   integer      gdamp_crlat_    ! gmp damping critical latitude  
#endif
#ifdef SMP
   real         gxlon_    ,&    ! location of i-grid
                gylat_          ! location of j-grid
   integer      lsmsk_    ,&    ! land-sea mask
                lobsl_          ! observation data
#endif
#ifdef RMP
   integer      igrd_     ,&    ! number of i-grid volume for rmp
                igrd1_    ,&
                igrd2_    ,&
                igrd12_   ,&
                igrd1i_   ,&
                igrd1o_   ,&
                igrdcut_  ,&
                jgrd_     ,&    ! number of j-grid volume for rmp
                jgrd1_    ,&
                jgrd12_   ,&
                jgrd1i_   ,&
                jgrd1o_   ,&
                jgrdcut_  ,&
                igrdp_    ,&    ! model grid for rmp : mpi partial grid
                jgrdp_    ,&    ! model grid for rmp : mpi partial grid
                cigrd1_   ,&    ! number of i-grid for coarse rmp
                cjgrd1_   ,&    ! number of j-grid for coarse rmp 
                irelx_    ,&    ! relaxation factor in implicit diffusion
                border_   ,&    ! number of buffer zone in the base grid
                bgf_      ,&    ! ratio of rmp grid resol over the base grid
                iwav_     ,&    ! number of spectral wave
                iwav1_    ,&    ! iwav_ + 1
                jwav_     ,&
                jwav1_    ,&
                lnwav_    ,&
                levr_     ,&    ! number of k-grid for rmp
                jumpr_    ,&    ! igrd_*2 + 3
                rgen_     ,&    ! generation model id
                lngrd_    ,&    ! number of physics grid
                clngrd_   ,&    ! number of physics grid for coarse grid
                rlpnt_    ,&      
                rltstp_   ,&
                rslvark_  ,&
                rmlvark_  ,&
                difuh_    ,&   ! factor of horizontal diffusions for T and q
                difum_    ,&   ! factor of horizontal diffusions for momentum 
                bsmooth_  ,&   ! smoothing factor for base fields
                nst_      ,&   ! parameter for time step (=1)
                rlevmax_       ! maximum number of model levels for output files 
   real         kvdif_         ! second-order vertical diffusion coefficient (m2/s)
#endif
#ifdef SSI
   integer      nt_       ,&
                nq_       ,&
                nu_       ,&
                nv_       ,&
                np_       ,&
                lonssi_   ,&
                latssi_   ,&
                maxlev_   ,&
                maxrep_   ,&
                maxcyl_   ,&
                maxpak_   ,&
                nqc_      ,&
                nmx_      ,&
                n_        ,&
                ntdata_   ,&
                nsdata_   ,&
                nwdata_   ,&
                npdata_   ,&
                nqdata_   ,&
                npwdat_   ,&
                nsprof_   ,&
                nsigdivt_ ,&
                jcapdivt_ ,&
                jcapstat_ ,&
                nmdszh_   ,&
                nsigsat_
#endif
#ifdef SMP_RA2SFC
   integer      ilonfscm_ ,&
                ilatgscm_ ,&
                ijcapscm_
#endif
#ifdef MP
   integer      npes_     ,&
                ncol_     ,&
                nrow_     ,&
                levsp_    ,&
                levsr_    ,&
                levhp_    ,&
                levp1p_   ,&
                levm1p_   ,&
                lonfp_    ,&
                latgp_    ,&
                lonf2p_   ,&
                lonf22p_  ,&
                latg2p_   ,&
                jcapp_    ,&
                jcap1p_   ,&
                lcapp_    ,&
                lcap22p_  ,&
                lntp_     ,&
                lnt2p_    ,&
                lnt22p_   ,&
                llnp_     ,&
                lln2p_    ,&
                lln22p_
#ifdef RMP
   integer      igrd1p_   ,&
                igrd12p_  ,&
                jgrd1p_   ,&
                jgrd12p_  ,&
                lngrdp_   ,&
                iwav1p_   ,&
                levrp_    ,&
                llwavp_   ,&
                lnwavp_
#endif
#else
   integer      npes_
#endif
!
#ifndef RMP
   integer      LONFD     ,&
                LATGD     ,&
                LONF2F    ,&
                LATG2F    ,&
                LONF22F   ,&
                LATGS     ,&
                LONF2S    ,&
                LONF22S   ,&
                LATG2S    ,&
                LNT22S    ,&
                LNT2S     ,&
                TWOJ1S    ,&
                JCAPS     ,&
                JCAP1S    ,&
                LEVSS     ,&
                LEVHS     ,&
                LCAPS     ,&
                LCAP22S   ,&
                LLN2S     ,&
                LLN22S    ,&
                NCPUSS
#ifdef NIM
   integer      :: rc_
   integer      :: nip_
#endif
#else   /* RMP */ 
   integer      LONFD     ,&
!               LATGS     ,&
                LATGD     ,&
                LONF2F    ,&
                LATG2F    ,&
                LONF2S    ,&
                LONF22S   ,&
                LATG2S    ,&
                LNT22S    ,&
                LEVSS     ,&
                LEVHS     ,&
                LCAP22S   ,&
                JCAPS     ,&
                JCAP1S    ,&
                LNGRDS    ,&
                LNWAVS    ,&
                LLWAVS    ,&
                NCPUSS    ,&
                IWAV1S    ,&
                TWOJ1S
#ifndef MP
   integer      LCAPS     ,&
                LLN2S     ,&
                LLN22S
#endif
#endif

#ifdef RMP
   integer      IGRD12S   ,&
                JGRD12S   ,&
                IGRD1S
#endif
!
   contains
!-------------------------------------------------------------------------------
#ifdef NIM
   subroutine para_init(igme)
#else
   subroutine para_init
#endif
!-------------------------------------------------------------------------------
#ifndef NIM
   namelist /paralist/                                                         &
   ncpus_   ,                                                                  &
   npes_    ,                                                                  &
   ncol_    ,                                                                  &
#ifdef RMP
   igrd_    ,                                                                  &
   jgrd_    ,                                                                  &
#ifdef C2R
   cigrd1_  ,                                                                  &
   cjgrd1_  ,                                                                  &
#endif
#endif
#ifdef SSI
   nt_      ,                                                                  &
   nq_      ,                                                                  &
   nu_      ,                                                                  &
   nv_      ,                                                                  &
   np_      ,                                                                  &
   lonssi_  ,                                                                  &
   maxlev_  ,                                                                  &
   maxrep_  ,                                                                  &
   maxcyl_  ,                                                                  &
   maxpak_  ,                                                                  &
   nqc_     ,                                                                  &
   nmx_     ,                                                                  &
   n_       ,                                                                  &
   ntdata_  ,                                                                  &
   nsdata_  ,                                                                  &
   nwdata_  ,                                                                  &
   npdata_  ,                                                                  &
   nqdata_  ,                                                                  &
   npwdat_  ,                                                                  &
   nsprof_  ,                                                                  &
   nsigdivt_,                                                                  &
   jcapdivt_,                                                                  &
   jcapstat_,                                                                  &
   nmdszh_  ,                                                                  &
   nsigsat_ ,                                                                  &
#endif
   jcap_    ,                                                                  &
   levs_    ,                                                                  &
   ilonf_   ,                                                                  &
   ilatg_   ,                                                                  &
   ijcap_   ,                                                                  &
   ilevs_   ,                                                                  &
   irelx_   ,                                                                  &
   io_      ,                                                                  &
   jo_      ,                                                                  &
   ko_
#endif /* ~NIM */
!
!--- pre-defined value
!
   igen_    = 80
!
   kt_      = 6
!
   slvark_  = 80
   mlvark_  = 8
   lpnt_    = 30
!
   mtnres_  = _mtnres_
   mtnvar_  = _mtnvar_
   rmtnres_ = _rmtnres_
   rmtnvar_ = _rmtnvar_
!
   nstype_  = _nstype_
   nvtype_  = _nvtype_
!
   nwmass_   = _nwmass_
   nwsize_   = _nwsize_
   ngases_   = _ngases_
#ifdef GMPDAMP 
   gdamp_crlat_  = _gdamp_crlat_
#endif
#ifdef MUL_CLDTOP
   ncldtop_ = _ncldtop_
#endif
#ifdef SMP
   gxlon_   = _gxlon_
   gylat_   = _gylat_
   lsmsk_   = _lsmsk_
   lobsl_   = _lobsl_
#endif
#ifdef SMP_RA2SFC
   jcapscm_ = 62
   lonfscm_ = 192
   latgscm_ = 94
   lonf2scm_= lonfscm_*2
   latg2scm_= latgscm_/2
#endif
#ifdef RMP
   rgen_    = 80
   rlpnt_   = 50
   rltstp_  = 50
   rslvark_ = 50
   rmlvark_ = 50
   bsmooth_ = 1
   nst_     = 1
!
   irelx_   = 10
   border_  = 3
   bgf_     = 2
   difuh_   = 6.
   difum_   = 4.
   kvdif_   = _kvdif_
#endif
#ifdef RMP
   ltstp_   = 42
#else
   ltstp_   = 48
#endif
#ifdef SSI
   nt_      = 15
   nq_      = 10
   nu_      = 11
   nv_      = 11
   np_      = 6
#endif
#ifdef RIVER
   io2_     = _io2_
   jo2_     = _jo2_
   delxo2_  = _delxo2_
   delyo2_  = _delxo2_
#endif
!
#ifndef NIM
   read(1,paralist)     ! read basic parameters
#else
   call readnl(rc_)
!
   ncpus_ = 1
   jcap_  = 12
   levs_  = nz
   nip_   = 10*(2**glvl)**2+2
   ilonf_ = 1
   ilatg_ = 1
   ijcap_ = 1
   ilevs_ = nz
   ko_    = 17
#endif
!
! calculate other parameters
!
#ifndef SMP
   nvect_=64
   lonf_=3*jcap_+1
!
   do while(mod(lonf_,8).ne.0)
     lonf_=lonf_+1
   enddo
!
   latg_=lonf_/2
#ifndef DFS
!
   if(jcap_.eq.62) then
     latg_=94
   endif
!
#endif
   klot_=levs_
#ifdef IOSGB
   io_=lonf_
   jo_=latg_
#endif
!
!   if(jcap_<=62)then
!     io_= 144; jo_=  73
!   else if(jcap_>62 .and. jcap_<=170)then
!     io_= 360; jo_= 181
!   else if(jcap_>170 .and. jcap_<=426)then
!     io_= 720; jo_= 361
!   else if(jcap_>426 .and. jcap_<=1000)then
!     io_=1440; jo_= 721
!   else
!     io_=jcap_*2
!     jo_=io_/2+1
!   end if
!
#else /* SMP */
   nvect_=2
   lonf_=1
   latg_=2
   io_=1
   jo_=2
#endif /* SMP end */
   iod_=2*io_
   iodd_=2*io_+6
   siodd_=2*lonf_+6
   johf_=(jo_+1)/2
!
#ifndef NIM
   lonf2_=lonf_*2
#else
   lonf2_=igme
#endif
#if defined(SMP)
   lonf22_=lonf_*2+1
#elif defined(DFS)
   lonf22_=lonf_*2
#else
   lonf22_=lonf_*2+2
#endif
#ifdef NIM
   lonf22_=igme
#endif
#ifdef LATG2SSI
   latg2_=latg_/2+1
#else
   latg2_=latg_/2
#endif
#ifdef NIM
   latg2_=1
#endif
   jcap1_=jcap_+1
   jcap2_=jcap_+2
   lnt_=(jcap_+1)*(jcap_+2)/2
   lnt2_=(jcap_+1)*(jcap_+2)
   lnt22_=(jcap_+1)*(jcap_+2)+1
   lnut22_=((jcap_+2)*(jcap_+1)/2+(jcap_+1))*2+1
   lnut_=(jcap_+1)*(jcap_+2)/2+jcap_+1
   lnut2_=2*((jcap_+1)*(jcap_+2)/2+jcap_+1)
   lnuv_=(jcap_+1)*(jcap_+2)
   lnu2_=(jcap_+1)*(jcap_+2)*2
   twoj1_=2*(jcap_+1)
#ifdef RMP
   igrd1_=igrd_+1
   igrd2_=igrd_*2
   igrd12_=(igrd_+1)*2
   igrd1i_=igrd_+1
   igrd1o_=igrd_+1
   igrdcut_=igrd_-12
   jgrd1_=jgrd_+1
   jgrd12_=(jgrd_+1)/2
   jgrd1i_=jgrd_+1
   jgrd1o_=jgrd_+1
   jgrdcut_=jgrd_-12
   igrdp_=1
   jgrdp_=1
   iwav_=(igrd_-12)/3*2
   iwav1_=(igrd_-12)/3*2+1
#ifdef IDEAL
   jwav_=iwav_*jgrd_/(igrd_*2)*2
   jwav1_=jwav_+1
   lnwav_=iwav1_*jwav1_
#else
   jwav_=(jgrd_-12)/3*2
   jwav1_=(jgrd_-12)/3*2+1
   lnwav_=((igrd_-12)/3*2+1)*((jgrd_-12)/3*2+1)
#endif
   jumpr_=igrd_*2+3
   lngrd_=(igrd_+1)*(jgrd_+1)
!   cigrd1_=lonf_
!   cjgrd1_=latg_
   clngrd_=cigrd1_*cjgrd1_
   levr_=levs_
   if(levs_.le.100) then
     rlevmax_=100
   else
     rlevmax_=200
   endif
#endif
!
   nwater_   = nwmass_+nwsize_
   ntotal_   = nwater_+ngases_
   ncloud_   = nwmass_-1
   icloud_   = 2
   kcloud_   = levs_+1
   iwsize_   = nwmass_+1
   kwsize_   = nwmass_*levs_+1
   igases_   = nwater_+1
   kgases_   = nwater_*levs_+1
   ibltke_   = igases_+1
   kbltke_   = igases_*levs_+1
!
   levh_=levs_*ntotal_
   levp1_=levs_+1
   levp2_=levs_+2
   levm1_=levs_-1
   lotanu_=4*levs_+1
   imn_=360*60/mtnres_
   jmn_=180*60/mtnres_
   imnr_=360*60/rmtnres_
   jmnr_=180*60/rmtnres_
!
#ifdef SSI
   latssi_=lonssi_/2+1
#endif
#ifdef SMP_RA2SFC
   ilonfscm_=lonfscm_
   ilatgscm_=latgscm_
   ijcapscm_=jcapscm_
#endif
#if defined (MP) ||  defined (OPENMP)
#ifdef MP
   nrow_=npes_/ncol_
   call mpiset_sub(jcap_,levs_,lonf_,latg_,npes_,ncol_,nrow_,                  &
               levsp_,levsr_,lonfp_,latgp_,jcapp_,lntp_,llnp_) 
   levhp_=levsp_*ntotal_
   levp1p_=levsp_+1
   levm1p_=levsp_-1
   lonf2p_=lonfp_*2
#ifdef DFS
   lonf22p_=lonfp_*2
#else
   lonf22p_=lonfp_*2+2
#endif
   latg2p_=latgp_/2
   jcap1p_=jcapp_+1
   lcapp_=(jcapp_+1)*2
   lcap22p_=(jcapp_+1)*4+2
   lnt2p_=lntp_*2
   lnt22p_=lntp_*2+1
   lln2p_=llnp_*2
   lln22p_=llnp_*2+1
#ifdef RMP
   if(npes_.eq.2) ncol_=1
   call rmpiset_sub(levr_,igrd_,jgrd_,npes_,ncol_,nrow_,                       &
                levrp_,igrd1p_,jgrd1p_,iwav1p_,lnwavp_,llwavp_)
   igrd12p_=igrd1p_*2
   jgrd12p_=jgrd1p_/2
   lngrdp_=igrd1p_*jgrd1p_
#endif
#endif /* MP end */
#else /* ~MP */
   ncpus_=1
   npes_=1
   ncol_=1
#endif /* MP or OPENMP end */
#ifdef SMP
   ncpus_=1
   npes_=1
   ncol_=1
#endif
!
#ifndef RMP
   LONFD=lonf_
   LATGD=latg_
   LONF2F=lonf2_
   LATG2F=latg2_
   LONF22F=LONF22_
#ifdef MP
   LATGS=latgp_
   LONF2S=lonf2p_
   LONF22S=lonf22p_
   LATG2S=latg2p_
   LNT22S=lnt22p_
   LNT2S=lnt2p_
   TWOJ1S=lcapp_
   JCAPS=jcapp_
   JCAP1S=jcap1p_
   LEVSS=levsp_
   LEVHS=levhp_
   LCAPS=lcapp_
   LCAP22S=lcap22p_
   LLN2S=lln2p_
   LLN22S=lln22p_
   NCPUSS=latg2_
#else
   LATGS=latg_
   LONF2S=lonf2_
   LONF22S=lonf22_
   LATG2S=latg2_
   LNT22S=lnt22_
   LNT2S=lnt2_
   TWOJ1S=twoj1_
   JCAPS=jcap_
   JCAP1S=jcap1_
   LEVSS=levs_
   LEVHS=levh_
   LCAPS=lonf_
   LCAP22S=lonf22_
   LLN2S=lnt2_
   LLN22S=lnt22_
   NCPUSS=ncpus_
#endif
#else   /* RMP */ 
!  LATGS=jgrd1_
   LONFD=igrd1_
   LATGD=jgrd1_
   LONF2F=igrd12_
   LATG2F=jgrd12_
#ifdef MP
   LONF2S=igrd12p_
#ifdef RMPVECTORIZE
   LONF22S=igrd12p_*jgrd12p_
#else
   LONF22S=igrd12p_
#endif
   LATG2S=jgrd12p_
   LNT22S=lnt22p_
   LEVSS=levsp_
   LEVHS=levhp_
   LCAP22S=lcap22p_
   JCAPS=jcapp_
   JCAP1S=jcap1p_
   LNGRDS=lngrdp_
   LNWAVS=lnwavp_
   LLWAVS=llwavp_
   NCPUSS=latg2_
   IWAV1S=iwav1p_
   TWOJ1S=lcapp_
#else
   LONF2S=igrd12_
   LONF22S=igrd12_
   LATG2S=jgrd12_
   LNT22S=lnt22_
   LEVSS=levs_
   LEVHS=levh_
   LCAP22S=lonf22_
   JCAPS=jcap_
   JCAP1S=jcap1_
   LCAPS=lonf_
   LLN2S=lnt2_
   LLN22S=lnt22_
   LNGRDS=lngrd_
   LNWAVS=lnwav_
   LLWAVS=lnwav_
   NCPUSS=ncpus_
   IWAV1S=iwav1_
   TWOJ1S=twoj1_
#endif
#endif

#ifdef RMP
#ifdef MP
   IGRD12S=igrd12p_
   JGRD12S=jgrd12p_
   IGRD1S=igrd1p_
#else
   IGRD12S=igrd12_
   JGRD12S=jgrd12_
   IGRD1S=igrd1_
#endif
#endif
#if defined(RMP) && defined(RMPVECTORIZE) && defined(MP)
   ILOTS=LONF2S*LATG2S
#else
   ILOTS=LONF2S
#endif
!
   end subroutine para_init
!-------------------------------------------------------------------------------
!
end module paramodel
