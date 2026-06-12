#include <define.h>
#ifndef NIM
   subroutine phys_main_driver
#else
   subroutine phys_main_driver(                                                &
                     uugrsp,vvgrsp,plgrp,tvgrsp,rqgrsp,wwgrsp,                 &
                     ggt0k,ggq0k,ggu0k,ggv0k,                                  &
                     rn2d,rc2d,sn2d,                                           &
#ifdef NIM_DIAG
                     evap_nim,shflx_nim,lhflx_nim,dlwsfc_nim,ulwsfc_nim,       &
                     rnet_nim,swhr_nim,xmu_nim,lwhr_nim,ttenpbl_nim,           &
                     o3_nim,olr_nim,toar_usw,sfcr_usw,sfcr_dsw,                &
                     toar_dsw,tnet_dsw,                                        &
#endif
                     si1d,sl1d,prsi_nim,prsl_nim,nim_z,nim_zm,                 &
                     nimlat,nimlon,its,dorad,                                  &
                     ips,ipe,igms,igme,igs,ige)
#endif
!-------------------------------------------------------------------------------
!
! subroutine: phys_main_driver         
!
! abstract:
! - computes model physics tendency terms
! - program starts with spectral coefficients temp. of vorticity, 
!   divergence, specific humidity, and ln((psfc).  
!   converts them to the gaussian grid at each latitude and 
!   calls solver_physcis for the northern and southern hemispheres 
!   at the same time.  after return from fidisr.  
!   completes calculation of tendencies of temp. div. and lnps.
!   specific humidity, and vorticity are predicted by sr. sigvor
!   all input/output  is via commons.
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatble with namelist input
!
! references : 
!   hong et al. (2013, apjas): global/regional integrated model system (grims)
!   park et al. (2013, mwr): dfs dynamical core
!   byun and hong (2007, j. climate): single-column model (SMP)
!   kanamitsu et al. (2002, bams): ncep dynamical seasonal forecast system 2000
!
!-------------------------------------------------------------------------------
   use varsfc, only    : msub_,lalbd_,lsoil_,nsoil_
   use comkindex, only : ky,kys
!soojin_couple
#ifdef AOMG
   use couple_oasis, only : cpl_send
   use coupling
   use mod_oasis
   use paramodel, only : lonf_, latg_
#ifdef MP
   use paramodel, only : lonfp_, latgp_, lonf2p_, latg2p_
#endif
#endif /* AOMG */
!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef NIM
   use varsfc, only    : numsfcs
   use comsfc
#endif
   use constant, only  : g_,cp_,rd_,hvap_,hfus_,pi_,rerth_
   use comfibm
#ifdef DFS
   use dfsvar, only    : mt,jl,jla,jlha,jlg,lnt2,ngs,nge                      ,&
                         ib,jb,jbw,levsp,levs,levhp,levh,ntotal,jls,jle       ,&
                         V2,D2,T2,Q2,psl2                                     ,&
                         V3,D3,T3,Q3,PSL3                                     ,&
#ifdef CLM_CWF
                         cgs                                                  ,&
#endif
                         VJ,DJ,TJ,QJ,PSJ                                      ,&
                         VOR,DIO,VXC,VYC,TAI,QAI,PRS,XPSL,YPSL,PLAP           ,&
                         coslat,tavexy                                        ,&
                         AMATm,DMATm,AMATm_,DMATm_,aMSQUAR,jcol2js            ,&
                         si=>sigmafull,del=>delsig,sl=>sigma,iope
#else
   use module_sph_semi_implicit, only : sph_semi_impl_integrate
#endif
   use paramodel, only : JCAP1S,LEVSS,LEVHS,LCAPS,LCAP22S,LONF2S,LATG2S       ,&
                         LATG2F,LONF22F,LONF22S,ncpus_,lnt2_,ngases_          ,&
#ifdef MUL_CLDTOP
                         ncldtop_                                             ,&
#endif
#ifndef DYNAMIC_ALLOC
                         ncpus=>ncpus_                                        ,&
#endif
                         levs_,levh_,jcap1_,lonf2_,latg2_
#ifdef NIM
   use paramodel, only : nip_,levp1_
#endif
!
#ifdef MP
#ifndef DFS
   use paramodel, only : lln22p_,lnt22p_,ntotal_,lcap22p_,                     &
                         lonf22_,lonf22p_,latg2p_,levsp_, levhp_
#endif
   use commpi
   use compspec
#endif
#ifdef DG
   use comgda
#endif
#ifdef DG3
   use diag_3d_module, only : diag_3d_get, diag_3d_arrange
   use comgda
#endif
#ifdef REDUCE_GRID
   use comreduce
#else
   use paramodel, only : lonf_
#endif
#ifdef PSPLIT
   use compsplit
#endif
   use comgrad, only : runrad
!
#ifdef SMP
   use paramodel, only : ntotal_
#ifndef CLM_CWF
   use comfcst, only : dtbdy,curtime,vvel,hour1,hour2,hours,houre,ftim1,ftim2
#else
   use comfcst, only : wdiv, hadq
#endif
#ifdef SMP_NUDGING_FLX
   use comfcst, only : scmsh, scmlh
#endif
#ifdef SMP_NUDGING_RAD
   use comfcst, only : psfc,tsfc,qsfc,prec,scmlwup,scmlwdn,scmswup,scmswdn
#endif
#endif /* SMP end */
#ifdef SAS_DIAG
   use comfcst, only : dcu,dcv,dct,dcq,dch,fcu,fcd,                            &
                       deltb,delqb,delhb,cbmf,dlt,dlq,dlh
#else
#ifdef GWDC
   use comfcst, only : dct,fcu,fcd
#endif
#endif /* SAS_DIAG end */
#ifdef EXPLICIT_CLOUDINESS
   use comfcst, only : qcicps,qrscps,taucld,cldwp,cldip
#endif
#if defined(CLM_CWF) && !defined(DFS)
   use comfcst, only : cgs
#endif
#ifdef RAS
   use phys_ras_module, only  : ras_setup
#endif
#ifdef RASV2
#undef RAS
   use phys_ras2_module, only : ras2_setup
   use comfcst, only : mct, krmin, krmax, kfmax, ncrnd, kctop,                 &
                       sig, sgb, prj, rasal,dsfc
#endif
#ifdef NISLQ
   use nislq    , only : slq_q3
#endif
#ifdef STOCH
   use comfver, only : kdt
   use stochastic
   use module_file_write, only : file_write_bin
#endif
   use module_file_write, only : file_write_bin2
#ifdef CHEM
   use flexaod, only : calc_aod
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#ifdef NIM
!
! ---- NIM passing argument declaration
!
   integer,intent(in )   ::  ips,ipe,igms,igme,igs,ige
   real*8,intent(in  )   ::  plgrp(igms:igme)
   real*8,intent(in  )   ::  uugrsp(igms:igme,levs_)
   real*8,intent(in  )   ::  vvgrsp(igms:igme,levs_)
   real*8,intent(in  )   ::  tvgrsp(igms:igme,levs_)
   real*8,intent(in  )   ::  rqgrsp(igms:igme,levh_)
   real*8,intent(in  )   ::  wwgrsp(igms:igme,levp1_)
   real*8,intent(in  )   ::  si1d(levp1_),sl1d(levs_)
   real*8,intent(in  )   ::  prsi_nim(igms:igme,levp1_),prsl_nim(igms:igme,levs_)
   real*8,intent(in  )   ::  nim_z(igms:igme,levp1_), nim_zm(igms:igme,levs_)
   real*8,intent(in  )   ::  nimlat(igms:igme), nimlon(igms:igme)
   real*8,intent(out )   ::  ggt0k(igms:igme,levs_)
   real*8,intent(out )   ::  ggq0k(igms:igme,levh_)
   real*8,intent(out )   ::  ggu0k(igms:igme,levs_)
   real*8,intent(out )   ::  ggv0k(igms:igme,levs_)
   real*8,intent(out )   ::  rn2d(igms:igme), rc2d(igms:igme), sn2d(igms:igme)
   integer,intent(in  )  ::  its
!
!  local
!
   real*8               ::  gt0k(igms:igme,levs_), gq0k(igms:igme,levh_)
   real*8               ::  gu0k(igms:igme,levs_), gv0k(igms:igme,levs_)
   real*8               ::  rainl(igms:igme), rainc(igms:igme)
   real*8               ::  snow1(igms:igme)
#ifdef NIM_DIAG
   real*8,dimension(igms:igme,levs_),intent(out) :: swhr_nim,lwhr_nim,ttenpbl_nim
   real*8,dimension(igms:igme,levs_),intent(out) :: o3_nim
   real*8,dimension(igms:igme),intent(out)       :: evap_nim,shflx_nim,lhflx_nim
   real*8,dimension(igms:igme),intent(out)       :: dlwsfc_nim,ulwsfc_nim
   real*8,dimension(igms:igme),intent(out)       :: rnet_nim,xmu_nim
   real*8,dimension(igms:igme),intent(out)       :: olr_nim,toar_usw
   real*8,dimension(igms:igme),intent(out)       :: sfcr_usw,sfcr_dsw
   real*8,dimension(igms:igme),intent(out)       :: toar_dsw,tnet_dsw
   real*8,dimension(igms:igme,1:levs_)           :: swhr,lwhr,ttenpbl
   real*8,dimension(igms:igme,1:6)               :: rflux
#endif
   integer              ::  ign
#endif /* NIM ---- */
#include <abort.h>
!
! syn(1, 0*levs_+1, lan)  y
! syn(1, 1*levs_+1, lan)  rt
! syn(1, 2*levs_+1, lan)  x
! syn(1, 3*levs_+1, lan)  uln
! syn(1, 4*levs_+1, lan)  vln
! syn(1, 5*levs_+1, lan)  dpdphi
! syn(1, 5*levs_+2, lan)  dpdlam
! syn(1, 5*levs_+3, lan)  q
! syn(1, 5*levs_+4, lan)  qlap
!
! anl(1, 0*levs_+1, lan)  uu    dvdt
! anl(1, 1*levs_+1, lan)  vv    dtdt
! anl(1, 2*levs_+1, lan)  te    dtdt
! anl(1, 3*levs_+1, lan)  rq    drdt
!
   integer,parameter    ::  no3p=28 ,&
                            no3l=29
   integer              ::  lots    ,&
                            lotst   ,&
                            kst     ,&
                            ksr     ,&
                            ksd     ,&
                            ksu     ,&
                            ksv     ,&
                            kspphi  ,&
                            ksplam  ,&
                            ksp     ,&
                            ksplap  
   integer              ::  lota    ,&
                            lotat   ,&
                            kau     ,&
                            kav     ,&
                            kat     ,&
                            kar     
   integer              ::  klev
#ifndef DFS
   integer              ::  lotss   ,&
                            lotsts  ,&
                            ksts    ,&
                            ksrs    ,&
                            ksds    ,&
                            ksus    ,&
                            ksvs    ,&
                            kspphis ,&
                            ksplams ,&
                            ksps    ,&
                            ksplaps 
   integer              ::  lotas   ,&
                            lotats  ,&
                            kaus    ,&
                            kavs    ,&
                            kats    ,&
                            kars    
#ifdef MP
   real, allocatable    ::  syf(:,:,:)
   real, allocatable    ::  grs(:,:,:)
   real, allocatable    ::  anf(:,:,:)
   real, allocatable    ::  gra(:,:,:)
#endif
   real, allocatable    ::  syn(:,:,:)
   real, allocatable    ::  syntop(:,:,:)
   real, allocatable    ::  anl(:,:,:)
   real, allocatable    ::  anltop(:,:,:)
   real, allocatable    ::  flp(:,:,:,:)
   real, allocatable    ::  flm(:,:,:,:)
#endif /* ~DFS end */
#ifdef NISLQ
   real, dimension(LONF22S,levh_,LATG2S)  ::  slq_in,slq_out
#endif
#ifdef STOCH
   real, dimension(LONF2S,levs_)  ::  temp
#ifdef DFS
   real, dimension(LONF2S,levs_,LATG2S) ::  grid_u1,grid_v1,grid_t1,grid_q1
#endif
#endif /* STOCH end */
#ifdef STOCH_PRINT
   real, dimension(LONF2S,LATG2S,levs_) ::  stem
   real, dimension(LONF2S,levs_,LATG2S) ::  stc_utend,stc_vtend,stc_ttend,stc_qtend
#endif
!
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
#ifdef DFS
#define NCPUSS LATG2S
#endif
#ifdef DG
   real                 ::  tgmxl(NCPUSS),tgmx
   real                 ::  tgmnl(NCPUSS),tgmn
   integer              ::  igmxl(NCPUSS),kgmxl(NCPUSS),igmx,kgmx
   integer              ::  igmnl(NCPUSS),kgmnl(NCPUSS),igmn,kgmn
   integer              ::  jgmx,jgmn
   real                 ::  cldt(LONF2S,levs_,NCPUSS)
   real                 ::  clcv(LONF2S,levs_,NCPUSS)
#endif
#ifdef DG3
   real                 ::  gda(nwgda,kdgda,NCPUSS)
#endif
#undef NCPUSS
#if defined(SMP) || defined(DFS)
   logical,parameter    ::  ladj=.false.
#else
   logical,parameter    ::  ladj=.true.
#endif
!
#ifdef DFS
   real                 ::  sanm1(mt,jlg)
   real, allocatable    ::  syn(:,:,:),anl(:,:,:)
#endif
!
! local 
!
   integer              ::  jstr,jend,lons2,lats2,k,j,m,kk,lonsd2
#ifdef DFS
   integer              ::  nvar,lotl,lotg
#endif
   integer              ::  last,nggs,inclat,lat2,lat1,latdon,ngg,lat,lan,latrue
   integer              ::  llstr,llens,lnts2,lnoffset,n,l,lcapf
   integer              ::  i,lonff,kc,kaa,llensd
   real                 ::  rdelt2
!soojin_couple
#ifdef AOMG
   integer              :: intsec, jj
   real, allocatable    :: tevp(:,:), qnso(:,:), qsr(:,:), taux(:,:),            &
                           tauy(:,:), qnsi(:,:), ievp(:,:), train(:,:), tsnow(:,:)
#endif /* AOMG */
!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CLM_CWF
#ifndef DFS
   real                 :: cgsx(LONF2S,levs_)
#else
   real                 :: cgsx(ib*2,levs_)
#endif
#endif
   logical              ::  dorad
!
#ifdef MUL_CLDTOP
   real                 ::  xkt2(LONF2S*ncldtop_,LATG2S)
   real                 ::  me,xkt2m
   real,save            ::  seed0
   integer,allocatable,save  ::  NRND(:)
   integer              ::  iseed
   integer              ::  krsize
   logical,save         ::  first
!
   data first/.true./
   data krsize/20/
!
   if (first) then
     first = .false.
!
! machine dependency 
!  krsize ????
!
     CALL RANDOM_SEED(SIZE=krsize)
     allocate (nrnd(krsize))
     seed0 = idate(1) + idate(2) + idate(3) + idate(4)
   endif
   me = 0.0
   iseed = mod(100.0*sqrt(fhour*3600+100.0*me),1.0E9) + 1 + seed0
   nrnd(1:krsize) = iseed
#ifdef IBM
   call random_seed(generator=2)
#else
   CALL RANDOM_SEED(PUT=NRND)
#endif
   CALL RANDOM_NUMBER(xkt2)
#ifdef SMP
   xkt2m = 0.0
   do jj = 1, LATG2S
     do ii = 1, LONF2S*ncldtop_
       xkt2m = xkt2m + xkt2(ii,jj)
     enddo
   enddo
   xkt2m = xkt2m/real(LATG2S*LONF2S*ncldtop_)
   do jj = 1, LATG2S
     do ii = 1, LONF2S*ncldtop_
       xkt2(ii,jj) = xkt2m
     enddo
   enddo
#endif                  /* SMP */
#endif                  /* MUL_CLDTOP */
!
#ifdef RAS
#ifdef NIM
   real*8               ::  ci(levp1_), del(levs_), dtp, frain
#endif
   logical,parameter    ::  ras=.true.
   integer,parameter    ::  nsphys=1
   real,parameter       ::  cp=cp_, alhl=hvap_, grav=g_, rgas=rd_
   real                 ::  sig(levs_+1), prj(levs_+1), prh(levs_)
   real                 ::  fpk(levs_),   hpk(levs_),   sgb(levs_)
   real                 ::  ods(levs_), rasal(levs_),  prns(levs_/2)
   real                 ::  rannum(200),afac,ufac
   integer              ::  lm, krmin, krmax, nstrp, ncrnd
!
   sig=0. ; prj=0.   ; prh=0. ; fpk=0.    ; hpk=0. ; sgb=0. 
   ods=0. ; rasal=0. ; prns=0.; rannum=0.
#ifdef NIM
   do k = 1,levp1_
     ci(k)=1. - si1d(k)
   enddo
   do k = 1,levs_
     del(k)=ci(k+1)-ci(k)
   enddo
#endif
   call ras_setup(levs_, si, sl, del, cp, rgas, deltim                         &
               ,nsphys, thour                                                  &
               ,sig, sgb, prh, prj, hpk, fpk, ods, prns                        &
               ,rasal, lm, krmin, krmax, nstrp                                 &
               ,ncrnd, rannum, afac, ufac)
#endif
#ifdef RASV2
   logical,parameter    ::  ras=.true.
   integer,parameter    ::  nsphys=1
   real,parameter       ::  CP=Cp_,alhl=hvap_,grav=g_,rgas=rd_,hfus=hfus_
   real                 ::  rannum(200,5)
!
!  pdd is the lowest pressure above which downdraft is allowed
!  mct is a number of cloud types
!  
   real,parameter       ::  pdd=600.0
   real                 ::  frain, dtp
!
   frain = 0.5
   if(inistp.eq.1) frain = 1.
   dtp= deltim / frain
!
!  check fhour and thour
!
   call ras2_setup(levs_,  si,  sl, cp_, rd_, dtp, nsphys, fhour               &
                ,sig,   sgb, prj                                               &
                ,rasal, krmin, krmax, kfmax, ncrnd, rannum                     &
                ,mct,   kctop, deltim, dsfc, LONF2S*LATG2S)
#endif
!
#ifdef NISLQ
   lots    =4*levs_+4
#else
   lots    =4*levs_+levh_+4
#endif
   lotst   =2*levs_+1
   kst     =1
#ifndef NISLQ
   ksr     =kst+levs_
   ksd     =ksr+levh_
#else
   ksd     =kst+levs_
#endif
   ksu     =ksd+levs_
   ksv     =ksu+levs_
   kspphi  =ksv+levs_
   ksplam  =kspphi+1
   ksp     =ksplam+1
   ksplap  =ksp   +1
!
#ifndef NISLQ
   lota    =3*levs_+levh_
#else
   lota    =3*levs_
#endif
   lotat   =2*levs_
   kau     =1
   kav     =kau+levs_
   kat     =kav+levs_
#ifndef NISLQ
   kar     =kat+levs_
#endif
!
#ifndef DFS
#ifndef NISLQ
   lotss   =4*LEVSS+LEVHS+4
#else
   lotss   =4*LEVSS+4
#endif
   lotsts  =2*LEVSS+1
   ksts    =1
#ifndef NISLQ
   ksrs    =ksts+LEVSS
   ksds    =ksrs+LEVHS
#else
   ksds    =ksts+LEVSS
#endif
   ksus    =ksds+LEVSS
   ksvs    =ksus+LEVSS
   kspphis =ksvs+LEVSS
   ksplams =kspphis+1
   ksps    =ksplams+1
   ksplaps =ksps   +1
!
#ifndef NISLQ
   lotas   =3*LEVSS+LEVHS
#else
   lotas   =3*LEVSS
#endif
   lotats  =2*LEVSS
   kaus    =1
   kavs    =kaus+LEVSS
   kats    =kavs+LEVSS
#ifndef NISLQ
   kars    =kats+LEVSS
#endif
!
#ifdef MP
   allocate(syf(LONF22F,lotss,LATG2S))
   allocate(grs(LONF22S,lots,LATG2S))
   allocate(anf(LONF22F,lotas,LATG2S))
   allocate(gra(LONF22S,lota,LATG2S))
#endif
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus
#endif
#ifdef AOMG
!soojin_couple
!   allocate(syn(LCAP22S,lotss,NCPUSS))
#else
   allocate(anf(LONF22F,lotas,LATG2S))
#endif /*AOMG*/
   allocate(syn(LCAP22S,lotss,NCPUSS))
   allocate(syntop(2,JCAP1S,lotsts))
   allocate(anl(LCAP22S,lotas,NCPUSS))
   allocate(anltop(2,JCAP1S,lotats))
   allocate(flp(2,JCAP1S,lotas,NCPUSS))
   allocate(flm(2,JCAP1S,lotas,NCPUSS))
#endif/* ~DFS end */
#undef NCPUSS
!
#ifdef DFS
#ifdef MP
   allocate(syn(LONF22S,lots,LATG2S),anl(LONF22S,lota,LATG2S))
#else
   allocate(syn(LCAP22S,lots,ncpus),anl(LCAP22S,lota,ncpus))
#endif
#endif /* DFS end */
#ifdef STOCH
!
! stochastic forcing
!
#if defined(DYN_U) || defined(PHY_U)
   call stc_forcing(stc_u)
#endif
#if defined(DYN_V) || defined(PHY_V)
   call stc_forcing(stc_v)
#endif
#if defined(DYN_T) || defined(PHY_T)
   call stc_forcing(stc_t)
#endif
#if defined(DYN_Q) || defined(PHY_Q)
   call stc_forcing(stc_q)
#endif
#endif /* STOCH end */
!
! initialize local variables
!
   syn=0. ; anl=0.
#ifndef DFS
#ifdef MP
   syf=0.    ; grs=0.    ; anf=0. ; gra=0.
#endif
   syntop=0. ; anltop=0. ; flp=0. ; flm=0.
#endif /* ~DFS end */
#ifdef DG
   tgmxl=0. ;  tgmnl=0. ; igmxl=0. ;  kgmxl=0. ; igmnl=0. ; kgmnl=0.
   cldt=0. ; clcv=0.
#endif
#ifdef DG3
   gda  =0.
#endif
#ifdef DFS
   sanm1=0.
#endif
#ifdef CLM_CWF
   cgsx=0.
#endif
#ifdef NIM
#ifdef NIM_DBG
   if(iope) then
     print*, '0igs,ige iope',iope,igs,ige,levs_,lbound(prsl_nim),ubound(prsl_nim)
     print*, '0prsl_nim 1-Phys_main_driver', (prsl_nim(1,k),k=1,levs_)
     print*, '0prsl_nim 2-Phys_main_driver', (prsl_nim(2,k),k=1,levs_)
   else
     print*, '0igs,ige iope',iope,igs,ige,levs_,lbound(prsl_nim),ubound(prsl_nim)
     print*, '0prsl_nim 5122-Phys_main_driver', (prsl_nim(1,k),k=1,levs_)
     print*, '0prsl_nim 5123-Phys_main_driver', (prsl_nim(2,k),k=1,levs_)
   endif
#endif
#ifdef NIM_DIAG
   shflx_nim=0.  ; lhflx_nim=0.  ; dlwsfc_nim=0.; ulwsfc_nim=0. ; rnet_nim=0.
   swhr_nim=0.   ; lwhr_nim=0.   ; ttenpbl=0.   ; ttenpbl_nim=0.; olr_nim=0.
   toar_usw=0.   ; sfcr_usw=0.   ; sfcr_dsw=0.  ; toar_dsw=0.   ; tnet_dsw=0.
#endif
   gt0k=0.       ; gq0k=0.       ; gu0k=0.       ; gv0k=0.
   rainl=0.      ; rainc=0.      ; snow1=0.
#endif  /* NIM */
#ifdef NISLQ
   slq_in=0.
   slq_out=0.
#endif
!
#ifdef NIM
   call rad_prepare(idate,dorad,nimlat,nimlon,igs,ige)
#else
   call rad_prepare(idate,dorad)
#endif
#ifndef NIM
   if(ngases_.ge.1) call rad_ozone_setup(idate,fhour,no3p,no3l)
#endif
!
#ifdef CHEM
   if (dorad) call calc_aod
#endif
!
#ifndef NIM
#ifdef MP
#ifdef DFS
   jstr=jls
   jend=jle
   lons2=ib*2
   lats2=jb
#else
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   lons2=lonlen(mype)*2
   lats2=latlen(mype)
   lnts2=lntlen(mype)*2
   lnoffset=lntstr(mype)*2
#endif
#else
   llstr=0
   llens=jcap1_
   lons2=lonf2_
   lats2=latg2_
   lnts2=lnt2_
   lnoffset=0
#endif
#ifdef PSPLIT
!
   do k = 1,levs_
     do i = 1,LNT22S
         y (i,k)=tem(i,k)
         rt(i,k)= rm(i,k)
         x (i,k)=dim(i,k)
         w (i,k)=zem(i,k)
     enddo
   enddo
   do i = 1,LNT22S
     q(i)=qm(i)
   enddo
#endif
#ifdef CLM_CWF
   call sas_bh2007_avg(cgs,colrad)
#endif					/* CLM_CWF */
#else                   /* for NIM */
#define LATX 1
#endif					/* ~NIM end */
!=================================================================
#if !defined(SMP) && !defined(NIM)
#ifndef DFS
#ifdef MP
   call mpnn2n (q ,lnt22p_,       qa ,lln22p_,       1)
#ifndef NISLQ
!
!  y,rt,x,w
!
   call mpnn2nk(z(1,ky),lnt22p_,levs_,za(1,kys),lln22p_,levsp_,3+ntotal_)
#else
!
!  y,x,w
!
   call mpnn2nk(z(1,ky),lnt22p_,levs_,za(1,kys),lln22p_,levsp_,3)
#endif
#define QS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define XS xa
#define WS wa
#define ULNS ulna
#define VLNS vlna
#else /* ~MP */
#define QS q
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define XS x
#define WS w
#define ULNS uln
#define VLNS vln
#endif /* MP end */
   call sph_del_log_ps(QS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,               &
               llstr,llens,lwvdef)
   call sph_solve_laplacian(QS,QLAPS,llstr,llens,lwvdef)
   call sph_divor2wind(XS,WS,ULNS,VLNS,syntop(1,1,1),                          &
               syntop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef XS
#undef WS
#undef ULNS
#undef VLNS
#ifndef PSPLIT
   do n=1,2
     do j = 1,JCAP1S
       do l=1,lotats
         anltop(n,j,l)=0.0
       enddo
     enddo
   enddo
#endif
#else /* DFS */
!
! wave to grid transform (vxc,vyc,dio,tai,qai,prs,xpsl,ypsl,plap)
!
   call dfs_wave2grid_physics
!
#ifdef BIN_DBG
   call file_write_bin2(334,dio,ib,jbw,4*levs+1,0)
   call file_write_bin2(334,vxc,ib,jbw,2*levs  ,0)
   call MPABORT
#endif
#endif /* ~DFS end */
#endif /* ~SMP ~NIM end */
#ifdef DG
   tgmx=-1.e20
   tgmn= 1.e20
#endif
#ifndef NIM
#ifdef DFS
   do k = 1,levsp
     do j = 1,jlg
       do m = 1,mt
         v2(m,j,k)=0.0
         d2(m,j,k)=0.0
         t2(m,j,k)=0.0
       enddo
     enddo
   enddo
#ifndef NISLQ
   do k = 1,levhp
     do j = 1,jlg
       do m = 1,mt
         q2(m,j,k)=0.0
       enddo
     enddo
   enddo
#endif
#else /* SPH */
!
#ifdef MP
#define LNT2X lln22p_
#define ZES zea
#define DIS dia
#define TES tea
#define RQS rqa
#else
#define LNT2X lnt2_
#define ZES ze
#define DIS di
#define TES te
#define RQS rq
#endif /* MP end */
#ifdef OPENMP

!$omp parallel do private(j,k)
#endif
   do k = 1,LEVSS
     do j = 1,LNT2X
       ZES(j,k)=0.0
       DIS(j,k)=0.0
       TES(j,k)=0.0
     enddo
   enddo
#ifndef NISLQ
   do k = 1,LEVHS
     do j = 1,LNT2X
       RQS(j,k)=0.0
     enddo
   enddo
#endif
#undef LNT2X
#undef ZES
#undef DIS
#undef TES
#undef RQS
#endif /* DFS end */
!
#ifndef MP
! compute latitude band limits
!
   last=mod(latg2_,ncpus)
   nggs=(latg2_-last)/ncpus
   if(last.ne.0)nggs=nggs+1
   inclat=ncpus
   lat1=1-ncpus
   lat2=0
   latdon=0
!
   do ngg = 1,nggs
     if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
     lat1=lat1+ncpus
     lat2=lat2+inclat
#endif
!
#ifdef DG3
#ifdef MP
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define LATX lan
#else
#define LATX lat
#endif
     do lat=lat1,lat2
       lan=lat-latdon
       call diag_3d_get(LATX,nwgda*kdgda,gda(1,1,lan))
     enddo
#undef LATX
#endif /* DG3 end */
#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif /* MP end */
#ifdef DFS
#define NCPUSS LATG2S
#endif
     do k = 1,NCPUSS
       tgmxl(k)=tgmx
       tgmnl(k)=tgmn
     enddo
#undef NCPUSS
#endif /* DG end */
!
#ifndef DFS
#ifndef SMP
#ifdef MP
     lat1=1
     lat2=latg2_
     latdon=0
#define ZS za
#define ZES zea
#else
#define ZS z
#define ZES ze
#endif /* MP end */
! first lat loop
#ifdef OPENMP

!$omp parallel do private(lat,lan,llensd)
#endif
     do lat=lat1,lat2
       lan=lat-latdon
       latrue=latdef(lat)
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
#ifndef NISLQ
!
!      y,rt,x
!
       call sph_sum_coeff(ZS(1,2*LEVSS+2)       ,syn(1,ksts,lan),qtt(1,lat)   ,&
                                       llstr,llensd,lwvdef,2*LEVSS+LEVHS)
!
!      uln,vln,pphi,plam,q,qlap
!
       call sph_sum_coeff(ZES(1,3*LEVSS+LEVHS+1),syn(1,ksus,lan),qtt(1,lat)   ,&
                                       llstr,llensd,lwvdef,2*LEVSS+4)
#else
!
!      y,x
!
       call sph_sum_coeff(ZS(1,2*LEVSS+2)       ,syn(1,ksts,lan),qtt(1,lat)   ,&
                                       llstr,llensd,lwvdef,2*LEVSS)
!
!      uln,vln,pphi,plam,q,qlap
!
       call sph_sum_coeff(ZES(1,3*LEVSS+1)      ,syn(1,ksus,lan),qtt(1,lat)   ,&
                                       llstr,llensd,lwvdef,2*LEVSS+4)
#endif
!
!      uln,vln,pphi top
!
       call sph_sum_coeff_top(syn(1,ksus,lan)   ,syntop(1,1,1)                ,&
                                       qvv(1,lat),llstr,llensd,lwvdef,2*LEVSS+1)
     enddo
#undef ZS
#undef ZES
!
#ifdef MP
     call mpnl2ny(syn,lcap22p_,latg2_,syf,lonf22_,latg2p_,lotss,1,lotss)
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define SYNS syf
#else
#define SYNS syn
#endif /* MP end */
!
#ifdef OPENMP

!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
     do lat=lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,                 &
                   lcapf,lonff,latdef(lat),1)
     enddo
!
#undef SYNS
#undef DYNS
#endif                          /* SMP */
!
#ifdef MP
#ifndef NISLQ
     klev=4+ntotal_
#else
     klev=4
#endif
     call mpnk2nx(syf,lonf22_,lotss,                                           &
                 grs,lonf22p_,lots,latg2p_,levsp_,levs_,ksts,kst,klev)
     call mpx2nx (syf,lonf22_,lotss,                                           &
                 grs,lonf22p_,lots,latg2p_,kspphis,kspphi,4)
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define SYNS grs
#define ANLS gra
#define LATX lan
#else
#define SYNS syn
#define ANLS anl
#define LATX lat
#endif
!
#else /* DFS */
!
#define SYNS syn
#define ANLS anl
#ifdef MP
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define LATX lan
#else
#define LATX lat
#endif
!
#endif /* ~DFS end */
!
#ifdef OPENMP

!$omp parallel do private(lat,lan,lonsd2)
#endif
     do lat=lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       lonsd2=lonfdp(lan,mype)*2
#else
       lonsd2=lonfd(latdef(lat))*2
#endif
#else
       lonsd2=lons2
#endif
#ifdef DFS
       do i = 1,ib
         SYNS(i   ,ksp   ,lan) = prs (i,LATX      )              ! in kPa
         SYNS(ib+i,ksp   ,lan) = prs (i,jbw-LATX+1)              ! in kPa
         SYNS(i   ,ksplam,lan) = xpsl(i,LATX      )
         SYNS(ib+i,ksplam,lan) = xpsl(i,jbw-LATX+1)
         SYNS(i   ,kspphi,lan) = ypsl(i,LATX      )
         SYNS(ib+i,kspphi,lan) = ypsl(i,jbw-LATX+1)
         SYNS(i   ,ksplap,lan) = plap(i,LATX      )
         SYNS(ib+i,ksplap,lan) = plap(i,jbw-LATX+1)
       enddo
       do k = 1,levs
         do i = 1,ib
           SYNS(i   ,ksd-1+k,lan) = dio(i,LATX      ,k)
           SYNS(ib+i,ksd-1+k,lan) = dio(i,jbw-LATX+1,k)
           SYNS(i   ,kst-1+k,lan) = tai(i,LATX      ,k)+tavexy(k)
           SYNS(ib+i,kst-1+k,lan) = tai(i,jbw-LATX+1,k)+tavexy(k)
           SYNS(i   ,ksu-1+k,lan) = vxc(i,LATX      ,k)
           SYNS(ib+i,ksu-1+k,lan) = vxc(i,jbw-LATX+1,k) ! ucos(lat)
           SYNS(i   ,ksv-1+k,lan) = vyc(i,LATX      ,k)
           SYNS(ib+i,ksv-1+k,lan) = vyc(i,jbw-LATX+1,k)
#ifdef STOCH
           grid_u1(i   ,k,lan) = dfsg_u1(i,LATX      ,k)
           grid_u1(ib+i,k,lan) = dfsg_u1(i,jbw-LATX+1,k)
           grid_v1(i   ,k,lan) = dfsg_v1(i,LATX      ,k)
           grid_v1(ib+i,k,lan) = dfsg_v1(i,jbw-LATX+1,k)
           grid_t1(i   ,k,lan) = dfsg_t1(i,LATX      ,k)+tavexy(k)
           grid_t1(ib+i,k,lan) = dfsg_t1(i,jbw-LATX+1,k)+tavexy(k)
           grid_q1(i   ,k,lan) = dfsg_q1(i,LATX      ,k)
           grid_q1(ib+i,k,lan) = dfsg_q1(i,jbw-LATX+1,k)
#endif
         enddo
       enddo
       do k = 1,levh
         do i=1,ib
#ifdef NISLQ
!
!          nislq moisture at n+1 time
!
           slq_in(i   ,k,lan) = slq_q3(i,LATX      ,k)
           slq_in(ib+i,k,lan) = slq_q3(i,jbw-LATX+1,k)
#else
           SYNS(i   ,ksr-1+k,lan) = qai(i,LATX      ,k)
           SYNS(ib+i,ksr-1+k,lan) = qai(i,jbw-LATX+1,k)
#endif
         enddo
       enddo
#else /* SPH */
#ifdef NISLQ
!
!      nislq moisture at n+1 time
!
!
       do k = 1,levh_
         do i=1,lonsd2
           slq_in(i,k,lan) = slq_q3(i,k,LATX)
         enddo
       enddo
#endif /* NISLQ end */
#endif /* DFS end*/
!
#ifdef CLM_CWF
#ifdef DFS
       jj=2*(LATX-1)+1
       do kk = 1,levs
         cgsx(1:ib,kk)=cgs(1:ib,jj,kk)
         cgsx(ib+1:2*ib,kk)=cgs(1:ib,jj+1,kk)
       enddo
#else
       do kk=1,levs_
         cgsx(1:lonsd2,kk)=cgs(1:lonsd2,LATX,kk)
       enddo
#endif
#endif /* CLM_CWF end */
#endif /* NIM end */
!
#undef NO_PHYSICS
#ifndef NO_PHYSICS
#ifdef MP
       if( lonsd2.gt.0 ) then
#endif
#ifdef NIM
         lonsd2=ige
#endif
         if(dorad) then
           call phys_rad_solver(lonsd2,lats2,                                  &
#ifdef NIM
#ifndef NONHYD
              uugrsp,vvgrsp,plgrp,tvgrsp,rqgrsp,xxgrsp,                        &
#else
              plgrp,tvgrsp,rqgrsp,wwgrsp,                                      &
#ifdef NIM_DIAG
              o3_nim,rflux,                                                    &
#endif
#endif
              ips,ipe,1,1,nimlat,nimlon,si1d,sl1d,prsi_nim,prsl_nim)
#else /* else of NIM */
#ifdef SMP
              x,w,z,y,rt,vvel,                                                 &
#else /* else of SMP */
#ifndef NONHYD
              SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),                 &
#ifdef NISLQ
              SYNS(1,kst,lan),slq_in(1,1,lan),SYNS(1,ksd,lan),                 &
#else
              SYNS(1,kst,lan),SYNS(1,ksr,lan),SYNS(1,ksd,lan),                 &
#endif
#ifndef SMP
              SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),                           &
#endif
#else /* else of NONHYD */
              SYNS(1,ksp,lan),SYNS(1,kst,lan),                                 &
#ifdef NISLQ
              slq_in(1,1,lan),SYNS(1,ksd,lan),                                 &
#else
              SYNS(1,ksr,lan),SYNS(1,ksd,lan),                                 &
#endif
#endif /* endif NONHYD */
#endif /* endif SMP */
#ifdef DG3
              gda(1,1,lan),                                                    &
#endif
#ifdef DG
              cldt(1,1,lan),clcv(1,1,lan),                                     &              
#endif
#ifdef EXPLICIT_CLOUDINESS
              qcicps(1,1,LATX), qrscps(1,1,LATX),                              &
              taucld(1,1,LATX), cldwp(1,1,LATX), cldip(1,1,LATX),              &
#endif
              LATX,latrue)
#endif /* endif NIM */
#ifdef NIM_DIAG
           do ign=igs,ige
             olr_nim(ign)  = rflux(ign,1)
             toar_usw(ign) = rflux(ign,2)
             sfcr_usw(ign) = rflux(ign,3)
             sfcr_dsw(ign) = rflux(ign,4)
             toar_dsw(ign) = rflux(ign,5)
             tnet_dsw(ign) = rflux(ign,6)
           enddo
#endif
#ifdef DG3
           call diag_3d_arrange(LATX,nwgda*kdgda,gda(1,1,lan))
#endif
         endif
#ifndef NO_PHYSICS
         if(inistp.ne.0)then
           runrad=.false.
         else
           runrad=.true.
         endif
#endif
!
#ifdef STOCH
#if defined(DYN_U) || defined(DYN_V) || defined(DYN_T) || defined(DYN_Q)
! stochastic forcing in total dynamical tendency
!
   if(kdt.gt.2) then
     do k=1,levs_
       do i=1,lonsd2
#ifdef DYN_U
         temp(i,k)=(SYNS(i,ksu+k-1,lan)-grid_u1(i,k,lan))*stc_u(i,k,lan)
         SYNS(i,ksu+k-1,lan)=grid_u1(i,k,lan)+temp(i,k)
#endif
#ifdef DYN_V
         temp(i,k)=(SYNS(i,ksv+k-1,lan)-grid_v1(i,k,lan))*stc_v(i,k,lan)
         SYNS(i,ksv+k-1,lan)=grid_v1(i,k,lan)+temp(i,k)
#endif
#ifdef DYN_T
         temp(i,k)=(SYNS(i,kst+k-1,lan)-grid_t1(i,k,lan))*stc_t(i,k,lan)
         SYNS(i,kst+k-1,lan)=grid_t1(i,k,lan)+temp(i,k)
#endif
#ifdef DYN_Q
#ifdef NISLQ
         temp(i,k)=(slq_in(i,k,lan)-grid_q1(i,k,lan))*stc_q(i,k,lan)
         slq_in(i,k,lan)=grid_q1(i,k,lan)+temp(i,k)
#else
         temp(i,k)=(SYNS(i,ksr+k-1,lan)-grid_q1(i,k,lan))*stc_q(i,k,lan)
         SYNS(i,ksr+k-1,lan)=grid_q1(i,k,lan)+temp(i,k)
#endif
#endif /* DYN_Q end */
       enddo
     enddo
   endif
#endif /* DYN end */
#endif /* STOCH end */
!
         call phys_main_solver(lonsd2,                                         &
#ifdef NIM
              uugrsp,vvgrsp,plgrp,tvgrsp,rqgrsp,wwgrsp,                        &
              gt0k,gq0k,gu0k,gv0k,                                             &
#else
#ifdef SMP
              x,w,z,y,rt,vvel,te,rq,di,ze,                                     &
#ifdef SMP_NUDGING_FLX
              scmsh,scmlh,                                                     &
#endif
#ifdef SMP_NUDGING_RAD
              psfc,tsfc,qsfc,prec,                                             &
              scmlwup,scmlwdn,scmswup,scmswdn,                                 &
#endif
#else
              SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),                 &
#ifdef NISLQ
              SYNS(1,kst,lan),slq_in(1,1,lan),SYNS(1,ksd,lan),                 &
#else
              SYNS(1,kst,lan),SYNS(1,ksr,lan),SYNS(1,ksd,lan),                 &
#endif
              SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),SYNS(1,ksplap,lan),        &
#ifdef PSPLIT
              dtdtm(1,1,LATX),drdtm(1,1,LATX),                                 &
              dudtm(1,1,LATX),dvdtm(1,1,LATX),                                 &
#else
#ifdef NISLQ
              ANLS(1,kat,lan),slq_out(1,1,lan),                                &
#else
              ANLS(1,kat,lan),ANLS(1,kar,lan),                                 &
#endif
              ANLS(1,kau,lan),ANLS(1,kav,lan),                                 &
#endif				/* PSPLIT */
#endif				/* SMP */
#endif				/* NIM */
#ifdef DG
              tgmxl(lan),igmxl(lan),kgmxl(lan),                                &
              tgmnl(lan),igmnl(lan),kgmnl(lan),                                &
#endif
#ifdef DG3
              gda(1,1,lan),                                                    &
#endif
#ifdef RAS
              ras,levs_,cp,alhl,grav,rgas,                                     &
              sig, sgb, prh, prj, hpk, fpk, ods, prns,                         &
              rasal, lm, krmin, krmax, nstrp,                                  &
              ncrnd, rannum, afac, ufac,                                       &
#endif
#ifdef RASV2
              ras,rgas, cp, grav, alhl,                                        &
              sig,prj,sgb,rasal, rannum, dsfc(1,LATX),                         &
              pdd,krmin, krmax, kfmax, ncrnd, mct,kctop,                       &
#endif
#ifdef MUL_CLDTOP
              xkt2(1,LATX),ncldtop_,                                           &
#endif
#ifdef CLM_CWF
              cgsx,                                                            &
#endif
#ifdef SAS_DIAG
              dcu(1,1,LATX), dcv(1,1,LATX), dct(1,1,LATX),                     &
              dcq(1,1,LATX), dch(1,1,LATX),                                    &
              fcu(1,1,LATX), fcd(1,1,LATX),                                    &
              deltb(1,LATX), delqb(1,LATX), delhb(1,LATX), cbmf(1,LATX),       &
              dlt(1,1,LATX), dlq(1,1,LATX), dlh(1,1,LATX),                     &
#else
#ifdef GWDC
              dct(1,1,LATX), fcu(1,1,LATX), fcd(1,1,LATX),                     &
#endif
#endif
#ifdef NIM
              rainl,rainc,snow1,                                               &
#ifdef NIM_DIAG
              evap_nim,shflx_nim,lhflx_nim,dlwsfc_nim,ulwsfc_nim,              &
              rnet_nim,swhr_nim,xmu_nim,lwhr_nim,ttenpbl_nim,                  &
#endif
              si1d,sl1d,prsi_nim,prsl_nim,nim_z,nim_zm,its,ips,ipe,            &
#endif
#ifdef RMP
              LATX,1.0)
#else
#ifndef VIC
              LATX,1.0)
#else
              LATX,1.0,idate)
#endif
#endif			/* RMP */
!
#ifdef STOCH
#if defined(PHY_U) || defined(PHY_V) || defined(PHY_T) || defined(PHY_Q)
! stochastic forcing in total physical tendency
!
     do k = 1,levs_
       do i = 1,lonsd2
#ifdef PHY_U
         temp(i,k)=(ANLS(i,kau+k-1,lan)-SYNS(i,ksu+k-1,lan))*stc_u(i,k,lan)
         ANLS(i,kau+k-1,lan)=SYNS(i,ksu+k-1,lan)+temp(i,k)
#endif
#ifdef PHY_V
         temp(i,k)=(ANLS(i,kav+k-1,lan)-SYNS(i,ksv+k-1,lan))*stc_v(i,k,lan)
         ANLS(i,kav+k-1,lan)=SYNS(i,ksv+k-1,lan)+temp(i,k)
#endif
#ifdef PHY_T
         temp(i,k)=(ANLS(i,kat+k-1,lan)-SYNS(i,kst+k-1,lan))*stc_t(i,k,lan)
         ANLS(i,kat+k-1,lan)=SYNS(i,kst+k-1,lan)+temp(i,k)
#endif
#ifdef PHY_Q
#ifdef NISLQ
         temp(i,k)=(slq_out(i,k,lan)-slq_in(i,k,lan))*stc_q(i,k,lan)
         slq_out(i,k,lan)=slq_in(i,k,lan)+temp(i,k)
#else
         temp(i,k)=(ANLS(i,kar+k-1,lan)-SYNS(i,ksr+k-1,lan))*stc_q(i,k,lan)
         ANLS(i,kar+k-1,lan)=SYNS(i,ksr+k-1,lan)+temp(i,k)
#endif
#endif /* PHY_Q end */
       enddo
     enddo
#endif /* PHY end */
#ifdef STOCH_PRINT
     do k = 1,levs_
       do i = 1,lonsd2
         stc_utend(i,k,lan)=(ANLS(i,kau+k-1,lan)-grid_u1(i,k,lan))*0.5/deltim
         stc_vtend(i,k,lan)=(ANLS(i,kav+k-1,lan)-grid_v1(i,k,lan))*0.5/deltim
         stc_ttend(i,k,lan)=(ANLS(i,kat+k-1,lan)-grid_t1(i,k,lan))*0.5/deltim
         stc_qtend(i,k,lan)=(ANLS(i,kar+k-1,lan)-grid_q1(i,k,lan))*0.5/deltim
       enddo
     enddo
#endif
#endif /* STOCH end */
!
#ifdef DG
         if(tgmxl(lan).gt.tgmx) then
           tgmx=tgmxl(lan)
           igmx=igmxl(lan)
           kgmx=kgmxl(lan)
           jgmx=lat
         else if(tgmnl(lan).lt.tgmn) then
           tgmn=tgmnl(lan)
           igmn=igmnl(lan)
           kgmn=kgmnl(lan)
           jgmn=lat
         endif
#endif
#ifdef MP
       endif
#endif
#endif /* ~NO_PHYSICS end */
!
#ifndef NIM
#ifdef DFS
       do k = 1,levs
           do i = 1,ib
             tai(i,LATX      ,k)=ANL(i   ,kat-1+k,lan)-tavexy(k)
             tai(i,jbw-LATX+1,k)=ANL(ib+i,kat-1+k,lan)-tavexy(k)
             vxc(i,LATX      ,k)=ANL(i   ,kau-1+k,lan)
             vxc(i,jbw-LATX+1,k)=ANL(ib+i,kau-1+k,lan)
             vyc(i,LATX      ,k)=ANL(i   ,kav-1+k,lan)
             vyc(i,jbw-LATX+1,k)=ANL(ib+i,kav-1+k,lan)
           enddo
       enddo
       do k = 1,levh
           do i = 1,ib
#ifdef NISLQ
             slq_q3(i,LATX      ,k)=slq_out(i   ,k,lan)
             slq_q3(i,jbw-LATX+1,k)=slq_out(ib+i,k,lan)
#else
             qai(i,LATX      ,k)=ANL(i   ,kar-1+k,lan)
             qai(i,jbw-LATX+1,k)=ANL(ib+i,kar-1+k,lan)
#endif
           enddo
       enddo
#else /* SPH */
#ifdef NISLQ
!
!      nislq moisture update at n+1
!
       do k = 1,levh_
         do i=1,lonsd2
           slq_q3(i,k,LATX)=slq_out(i,k,lan)
         enddo
       enddo
#endif
#endif /* DFS end */
#undef LATX
     enddo
!

!soojin_couple
! do k = 1, levs_
!    do i = 1, LONF2S
!      write(*,124) "soojin grs after phy solver",k,lan,i,grs(i,ksu+k-1,lan),&
!                                              grs(i,ksv+k-1,lan),&
!                                              grs(i,kst+k-1,lan),&
!                                              gra(i,kau+k-1,lan),&
!                                              gra(i,kav+k-1,lan),&
!                                              gra(i,kat+k-1,lan)
!    enddo
!  enddo
!124 format (a,x,3(i3,x),6(f20.15,x))
!write(*,*) "soojin check var", maxval(cpl_dusfc(:,:)),minval(cpl_dusfc(:,:)) &
!                             , "dvsfc",maxval(cpl_dvsfc(:,:)),minval(cpl_dvsfc(:,:)) &
!                             , "evp",maxval(cpl_evp(:,:)),minval(cpl_evp(:,:)) &
!                             , "rain",maxval(cpl_raint(:,:)),minval(cpl_raint(:,:)) 

#ifdef AOMG
!soojin_coupled
!    total integration time (sec)
!! intsec last part -300 is -(dt)
     intsec = ( int(fhour)*3600 ) + int(shour) -300.
     print*, "intsec =", intsec

! first step intsec -> skip exchange
    if ( intsec == -150 ) then
     cpl_evp(:,:)    = 0.0    ; cpl_raint(:,:)  = 0.0
     cpl_dlwsfc(:,:) = 0.0    ; cpl_ulwsfc(:,:) = 0.0
     cpl_dtsfc(:,:)  = 0.0    ; cpl_dqsfc(:,:)  = 0.0
     cpl_sswdn(:,:)  = 0.0    ; cpl_sswup(:,:)  = 0.0
     cpl_dusfc(:,:)  = 0.0    ; cpl_dvsfc(:,:)  = 0.0
     cpl_ievp(:,:)   = 0.0    ; cpl_snow(:,:)   = 0.0
    else

    if ( mod(intsec, 3600) .eq. 0 ) then
     print*, "=== start AOMG exchange part ==="
     allocate ( tevp(LONF2S, LATG2S),  qnso(LONF2S, LATG2S) )
     allocate ( taux(LONF2S, LATG2S), tauy(LONF2S, LATG2S) )
     allocate ( qsr(LONF2S, LATG2S),  qnsi(LONF2S, LATG2S) )
     allocate ( train(LONF2S, LATG2S),tsnow(LONF2S, LATG2S))
     allocate ( ievp(LONF2S, LATG2S) )
     print*, "=== end of alloacte ==="

!!send variables
     if (intsec == 0) then
      cpl_evp(:,:)   = cpl_evp(:,:)/150. !dt
      cpl_raint(:,:) = cpl_raint(:,:)/150.
      cpl_dlwsfc(:,:)= cpl_dlwsfc(:,:)/150.
      cpl_ulwsfc(:,:)= cpl_ulwsfc(:,:)/150.
      cpl_dtsfc(:,:) = cpl_dtsfc(:,:)/150.
      cpl_dqsfc(:,:) = cpl_dqsfc(:,:)/150.
      cpl_sswdn(:,:) = cpl_sswdn(:,:)/3600.
      cpl_sswup(:,:) = cpl_sswup(:,:)/3600.
      cpl_dusfc(:,:) = cpl_dusfc(:,:)/300.
      cpl_dvsfc(:,:) = cpl_dvsfc(:,:)/300.
      cpl_ievp(:,:)  = cpl_ievp(:,:)/150.
      cpl_snow(:,:)  = cpl_snow(:,:)/150.
     else
      cpl_evp(:,:)   = cpl_evp(:,:)/3600. !exchange frequancy
      cpl_raint(:,:) = cpl_raint(:,:)/3600.
      cpl_dlwsfc(:,:)= cpl_dlwsfc(:,:)/3600.
      cpl_ulwsfc(:,:)= cpl_ulwsfc(:,:)/3600.
      cpl_dtsfc(:,:) = cpl_dtsfc(:,:)/3600.
      cpl_dqsfc(:,:) = cpl_dqsfc(:,:)/3600.
      cpl_sswdn(:,:) = cpl_sswdn(:,:)/3600.
      cpl_sswup(:,:) = cpl_sswup(:,:)/3600.
      cpl_dusfc(:,:) = cpl_dusfc(:,:)/3600.
      cpl_dvsfc(:,:) = cpl_dvsfc(:,:)/3600.
      cpl_ievp(:,:)  = cpl_ievp(:,:)/3600.
      cpl_snow(:,:)  = cpl_snow(:,:)/3600.
     endif
     tevp(:,:)=cpl_evp(:,:)
     qnso(:,:)=(cpl_dlwsfc(:,:)-cpl_ulwsfc(:,:))-cpl_dtsfc(:,:)-cpl_dqsfc(:,:)
     qnsi(:,:)=qnso(:,:)
     qsr(:,:)=cpl_sswdn(:,:)-cpl_sswup(:,:)
     taux(:,:)=cpl_dusfc(:,:)
     tauy(:,:)=cpl_dvsfc(:,:)
     ievp(:,:)=cpl_ievp(:,:)
     train(:,:)=cpl_raint(:,:)*1.e3
     tsnow(:,:)=cpl_snow(:,:)*1.e3

     cpl_evp(:,:)    = 0.0    ; cpl_raint(:,:)  = 0.0
     cpl_dlwsfc(:,:) = 0.0    ; cpl_ulwsfc(:,:) = 0.0
     cpl_dtsfc(:,:)  = 0.0    ; cpl_dqsfc(:,:)  = 0.0
     cpl_sswdn(:,:)  = 0.0    ; cpl_sswup(:,:)  = 0.0
     cpl_dusfc(:,:)  = 0.0    ; cpl_dvsfc(:,:)  = 0.0
     cpl_ievp(:,:)   = 0.0    ; cpl_snow(:,:)   = 0.0

!cpl_send(input,intsec,var_id)
!if (.false.) then
!     call cpl_send(taux, intsec, 4)
!     call cpl_send(tauy, intsec, 5)
!     call cpl_send(taux, intsec, 6)
!     call cpl_send(tauy, intsec, 7)
!     call cpl_send(qsr, intsec, 8)
!     call cpl_send(qnso, intsec, 9)
!     call cpl_send(qnsi, intsec, 10)
!     call cpl_send(tevp, intsec, 11)
!     call cpl_send(ievp, intsec, 12)
!     call cpl_send(train, intsec, 13)
!     call cpl_send(tsnow, intsec, 14)
!endif
!if (.true.) then
     call cpl_send(taux, intsec, 2)
     call cpl_send(tauy, intsec, 3)
     call cpl_send(qsr, intsec, 4)
     call cpl_send(qnso, intsec, 5)
     call cpl_send(tevp, intsec, 6)
!endif
     write(*,*) "grims_cpl_snd 1", minval(taux), maxval(taux)
     write(*,*) "grims_cpl_snd 2", minval(tauy), maxval(tauy)
     write(*,*) "grims_cpl_snd 3", minval(qsr ), maxval(qsr )
     write(*,*) "grims_cpl_snd 4", minval(qnso), maxval(qnso)
     write(*,*) "grims_cpl_snd 5", minval(qnsi), maxval(qnsi)
     write(*,*) "grims_cpl_snd 6", minval(tevp), maxval(tevp)
     write(*,*) "grims_cpl_snd 7", minval(ievp), maxval(ievp)
     write(*,*) "grims_cpl_snd 8", minval(train), maxval(train)
     write(*,*) "grims_cpl_snd 9", minval(tsnow), maxval(tsnow)

!!!!!!!!!!!!!!!!end send variables

     deallocate( tevp, qnso, qnsi, qsr, taux, tauy, ievp, train, tsnow )
     print*, "end of aomg part"
    endif !---mod
    endif !--end for skip exchange forward 
#endif /* AOMG */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef STOCH_PRINT
! print random & total tendency
!
#if defined(NLT_U) || defined(DYN_U) || defined(PHY_U)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_u(i,k,j)
#elif defined(NLT_V) || defined(DYN_V) || defined(PHY_V)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_v(i,k,j)
#elif defined(NLT_T) || defined(DYN_T) || defined(PHY_T)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_t(i,k,j)
#elif defined(NLT_Q) || defined(DYN_Q) || defined(PHY_Q)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_q(i,k,j)
#else  /* control */
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=1.
#endif
       call file_write_bin(333,stem,LONF2S,LATG2S,levs_,1)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_utend(i,k,j)
       call file_write_bin(333,stem,LONF2S,LATG2S,levs_,1)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_vtend(i,k,j)
       call file_write_bin(333,stem,LONF2S,LATG2S,levs_,1)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_ttend(i,k,j)
       call file_write_bin(333,stem,LONF2S,LATG2S,levs_,1)
       forall(i=1:LONF2S,j=1:LATG2S,k=1:levs_) stem(i,j,k)=stc_qtend(i,k,j)
       call file_write_bin(333,stem,LONF2S,LATG2S,levs_,1)
#endif /* STOCH_PRINT end */
!
#ifndef PSPLIT
#ifdef DG3
#ifdef MP
#define LATX lan
#else
#define LATX lat
#endif
     do lat=lat1,lat2
       lan=lat-latdon
       call diag_3d_arrange(LATX,nwgda*kdgda,gda(1,1,lan))
     enddo
#undef LATX
#endif
!
#undef SYNS
#undef DYNS
#undef ANLS

#ifndef SMP
#ifndef DFS
#ifdef MP
#ifndef NISLQ
     klev=3+ntotal_
#else
     klev=3
#endif
     call mpnx2nk(gra,lonf22p_,lota,                                           &
                anf,lonf22_ ,lotas,latg2p_,levs_,levsp_,kau,kaus,klev)
     if( levlen(mype).lt.levsp_ ) then
       do lan=1,lats2
         do i=1,lonf22_
           anf(i,kaus+levsp_-1,lan)=0.0
           anf(i,kavs+levsp_-1,lan)=0.0
           anf(i,kats+levsp_-1,lan)=0.0
         enddo
       enddo
#ifndef NISLQ
       do lan=1,lats2
         do i=1,lonf22_
           do kc = 1,ntotal_
             anf(i,kars+levsp_*kc-1,lan)=0.0
           enddo
         enddo
       enddo
#endif
     endif
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define ANLS anf
#else
#define ANLS anl
#endif
#ifdef OPENMP

!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       call sph_wave2grid(ANLS(1,1,lan),ANLS(1,1,lan),2*lotas,                 &
                  lcapf,lonff,latdef(lat),-1)
     enddo
#undef ANLS
!
#ifdef MP
     call mpny2nl(anf,lonf22_ ,latg2p_,                                      &
                anl,lcap22p_,latg2_ ,lotas,kaus,lotas)
!
     lat1=1
     lat2=latg2_
     latdon=0
#define ZES zea
#else
#define ZES ze
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,kaa,llensd)
#endif
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
!
       do k = 1,LEVSS*2
         kaa=kaus-1+k
         do j = 1,llensd*2
           anl(      j,kaa,lan)=anl(      j,kaa,lan)*rcs2(lat)
           anl(LCAPS+j,kaa,lan)=anl(LCAPS+j,kaa,lan)*rcs2(lat)
         enddo
       enddo
       call sph_grid2wave(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),          &
                             llensd,lotas)
     enddo
!
!    no multi-threads should be given in following loop
!    otherwise results will be non-reproduceble.
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
!
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
       call sph_grid2wave(flp(1,1,kaus,lan),flm(1,1,kaus,lan),ZES,qww(1,lat),  &
                                                      llstr,llensd,lwvdef,lotas)
#endif
#ifdef DEFAULT
       call sph_comp_coeff2(flp(1,1,kaus,lan),flm(1,1,kaus,lan),ZES,qww(1,lat),&
                                                      llstr,llensd,lwvdef,lotas)
#endif
!
       call sph_sum_wind_coeff(flp(1,1,kaus,lan),flm(1,1,kaus,lan),            &
                 flp(1,1,kavs,lan),flm(1,1,kavs,lan),                          &
                 anltop(1,1,1),anltop(1,1,LEVSS+1),                            &
                 qvv(1,lat),wgt(lat),                                          &
                 llstr,llensd,lwvdef,LEVSS)
     enddo
#undef ZES
#endif			/* DFS */
!
#endif			/* not SMP */
!
#endif  /* endif for ifndef psplit */
#if !defined(MP)
     latdon=latdon+(lat2-lat1+1)
   enddo
!
#endif
!
#ifdef PSPLIT
#ifndef DFS
#ifdef MP
   deallocate(syf, grs, anf, gra)
#endif
   deallocate(syn, syntop, anl, anltop, flp, flm)
#else /* SPH */
   deallocate(syn, anl)
#endif /* ~DFS end */
#else /* ~PSPLIT */
#ifndef SMP
!
!     input : ze=d(u)/d(t) di=d(v)/d(t)
!     output: uln=d(di)/d(t) vln=d(ze)/d(t)
!
#ifdef DFS
!
! grid to wave transform (Vj,Dj,Tj,Qj)
!
   call dfs_grid2wave_physics
!
   if(ladj) then
     do j = 1,jlg
       do m = 1,mt
         PSJ(m,j)=0.0
        enddo
     enddo
   endif
   rdelt2=1./(2*deltim)
   if(ladj) then
#ifdef OPENMP
!$omp parallel do private(m,j,k)
#endif
     do k = 1,levsp
       do j = 1,jlg
         do m = 1,mt
           DJ(m,j,k)=(Dj(m,j,k)-d3(m,j,k))*rdelt2
           TJ(m,j,k)=(tj(m,j,k)-t3(m,j,k))*rdelt2
           v3(m,j,k)=Vj(m,j,k)
         enddo
       enddo
     enddo
   else
     do k = 1,levsp
       do j = 1,jlg
         do m = 1,mt
           d3(m,j,k)=Dj(m,j,k)
           t3(m,j,k)=tj(m,j,k)
           v3(m,j,k)=Vj(m,j,k)
         enddo
       enddo
     enddo
   endif
!
#ifndef NISLQ
   do k = 1,levhp
     do j = 1,jlg
       do m = 1,mt
         q3(m,j,k)=qj(m,j,k)
       enddo
     enddo
   enddo
#endif
   if(ladj) then
     call dfs_semi_impl_integrate(D3,T3,PSL2,DJ,TJ,PSJ,deltim)
   endif
!
#else /* SPH */
#ifdef MP
#define UUS zea
#define VVS dia
#define XS ulna
#define WS vlna
#else
#define UUS ze
#define VVS di
#define XS uln
#define WS vln
#endif
   call sph_wind2divor(UUS,VVS,XS,WS,anltop(1,1,1),                            &
               anltop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef UUS
#undef VVS
#undef XS
#undef WS
!
#ifdef MP
#ifndef NISLQ
   klev=3+ntotal_   ! te,rq,uln,vln
#else
   klev=3           ! te,uln,vln
#endif
   call mpnk2nn(zea(1,2*levsp_+1),lln22p_,levsp_,ze(1,2*levs_+1),lnt22p_,levs_,klev)
#endif /* MP end */
#endif /* DFS end */
#endif /* SMP end */
!
#ifndef DFS
   if(ladj) then
     do j = 1,lnts2
       z(j,1)=0.0
     enddo
   endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
   do k = 1,levs_
     do j = 1,lnts2
       if(ladj) then
         di(j,k)=uln(j,k)-x(j,k)
         te(j,k)=te(j,k)-y(j,k)
       else
#ifdef SMP
         x(j,k)=di(j,k)
         y(j,k)=te(j,k)
       endif
       w(j,k)=ze(j,k)
#else
         x(j,k)=uln(j,k)
         y(j,k)=te(j,k)
       endif
       w(j,k)=vln(j,k)
#endif
     enddo
   enddo
#ifndef NISLQ
   do k = 1,levh_
      do j = 1,lnts2
         rt(j,k)=rq(j,k)
      enddo
   enddo
#endif
!
#ifdef SMP
   if(inistp.eq.1.and.ntotal_.ge.2) then
     do k = levs_+1,levh_
       do j = 1,lnts2
         rm(j,k)=rq(j,k)
       enddo
     enddo
   endif
#endif
   if(ladj) then
     call sph_semi_impl_integrate(x,y,q,di,te,z,uln,vln,lnts2,lnoffset)
   endif
#endif /* ~DFS end */
!
#ifdef DG
#ifdef MP
   if(iope) then
#endif
      write(6,'("phys_main_driver t range ",                                   &
                        2(4x,f6.1," @i,k,lat ",3i4))')                         &
           tgmx,igmx,kgmx,jgmx,tgmn,igmn,kgmn,jgmn
#ifdef MP
   endif
#endif
#endif					/* DG */
!
#ifndef DFS
#ifdef MP
   deallocate(syf, grs, anf, gra)
#endif
   deallocate(syn, syntop, anl, anltop, flp, flm)
#else /* DFS */
   deallocate(syn,anl)
#endif /* ~DFS end */
#endif/* PSPLIT end */
!
#else
!
   do ign = igs,ige
     do k = 1,levs_
       ggt0k(ign,k)=gt0k(ign,k)
       ggu0k(ign,k)=gu0k(ign,k)
       ggv0k(ign,k)=gv0k(ign,k)
     enddo
!    if(1.eq.1) stop
     do k = 1,levh_
       ggq0k(ign,k)=gq0k(ign,k)
     enddo
     rn2d(ign)=rainl(ign)*1.e3
     rc2d(ign)=rainc(ign)*1.e3
     sn2d(ign)=snow1(ign)*0.5*1.e3
   enddo
#endif /* ~NIM end */
!
   return
!
   end subroutine phys_main_driver
