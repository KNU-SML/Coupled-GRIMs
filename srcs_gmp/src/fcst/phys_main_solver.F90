#include <define.h>
#ifndef RMP
   subroutine phys_main_solver(lons2,                                          &
#ifdef NONHYD
                     uugrs,vvgrs,plgr,tvgrs,rqgrs,wwgrs,                       &
#else
                     uugrs,vvgrs,plgr,tvgrs,rqgrs,xxgrs,                       &
#endif
#if !defined (SMP) && !defined (NIM) 
                     pllamgr,plphigr,pslap,                                    &
#endif
                     ggt0,ggq0,ggu0,ggv0,                                      &
#ifdef SMP_NUDGING_FLX
                     scmsh,scmlh,                                              &
#endif
#ifdef SMP_NUDGING_RAD
                     psfc,tsfc,qsfc,prec,                                      &
                     scmlwup,scmlwdn,scmswup,scmswdn,                          &
#endif
#ifdef DG
                     tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,                      &
#endif
#ifdef DG3
                     gda,                                                      &
#endif
#else /* else of ~RMP */
   subroutine phys_main_solver(                                                &
                     ugrs,vgrs,pgr,tgrs,qgrs,                                  &
#ifndef NONHYD
                     xgrs,                                                     &
#ifndef SMP
                     plamgr,pphigr,                                            &
#endif /* ~SMP */
#endif /* ~NONHYD */ 
#ifdef PSPLIT
                     ggt0,ggq0,ggu0,ggv0,                                      &
#else
                     gt0,gq0,gu0,gv0,                                          &
#endif
#ifdef NONHYD_HYD
                     pgrs,wgrs,gp0,gw0,tgrh,                                   &
#endif /* NONHYD_HYD */
#ifdef DG
                     tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,                      &
#endif
#ifdef DG3
                     gda,                                                      &
#endif
#endif /* ~RMP */
#ifdef RAS
                     ras,lmx,cp,alhl,grav,rgas,                                &
                     sig, sgb, prh, prj, hpk, fpk, ods, prns,                  &
                     rasal, lm, krmin, krmax, nstrp,                           &
                     ncrnd, rannum, afac, ufac,                                &
#endif
#ifdef RASV2
                     ras,rgas, cp, grav, alhl,                                 &
                     sig,prj,sgb,rasal, rannum, dsfc,                          &
                     pdd,krmin, krmax, kfmax, ncrnd, mct,kctop,                &
#endif
#ifdef MUL_CLDTOP
                     xkt2,NCLDTOP,                                             &
#endif
#ifdef CLM_CWF
                     cgs2,                                                     &
#endif
#ifdef SAS_DIAG
                     dcu,dcv,dct,dcq,dch,fcu,fcd,                              &
                     deltb,delqb,delhb,cbmf,dlt,dlq,dlh,                       &
#else
#ifdef GWDC
                     dct,fcu,fcd,                                              &
#endif
#endif
#ifdef NIM
                     rainl,rainc,snow1,                                        &
#ifdef NIM_DIAG
                     evap_nim,shflx_nim,lhflx_nim,dlwsfc_nim,ulwsfc_nim,       &
                     rnet_nim,swhr_nim,xmu_nim,lwhr_nim,ttend_nim,             &
#endif
                     si1d,sl1d,prsi_nim,prsl_nim,nim_z,nim_zm,nits,ips,ipe,    &
#endif
#ifdef VIC
                     lat,count,idate)
#else
                     lat,count)
#endif
!
!-------------------------------------------------------------------------------
!
! subroutine: phys_main_solver
!
! abstract: main solver calling physics algorithms
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2008-10-01  kei yoshimura          add river model
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!   2012-02-24  suryun ham             river model implementation for rmp
!
! references : 
!   hong et al. (2013, apjas): global/regional integrated model system (grims)
!   park et al. (2013, mwr): dfs dynamical core
!   byun and hong (2007, j. climate): single-column model (SMP)
!   kanamitsu et al. (2002, bams): ncep dynamical seasonal forecast system 2000
!
!-------------------------------------------------------------------------------
   use varsfc, only    : lsoil_,nsoil_,msub_
#ifdef NIM
   use varsfc, only    : numsfcs
   use paramodel, only : glvl
   use comnim
#endif
   use constant, only  : g_,pi_,sbc_,rd_,rv_,qmin_,qmin8_,cal_,cp_,rdog_,rdorv_,      &
                         karman_,rhoh2o_,hfus_,hvap_,fv_,akapa_,rvordm1_,t0c_, &
                         omega_,cvap_,hsub_,rhoair0_,cliq_,cice_,psat_
#ifdef KSAS
   use constant, only  : cv_
#endif
   use paramodel, only : ILOTS,LATG2S,LONF22S,levs_,levp1_,levh_,              &
#if defined (RIVER) || defined (SMP)
                         LONF2S,                                               &
#endif
                         nwater_,ngases_,ntotal_,                              &
                         nwmass_,icloud_,kcloud_,ncloud_,igases_,kgases_,      &
                         iwsize_,kwsize_,nwsize_,ibltke_,kbltke_,              &
#ifndef RMP
                         jcap_
#else
                         igrd12_
#endif
#ifdef WDM
   use comccn,only     : iccn_
#endif
#ifndef MP
   use paramodel, only : latg2_
#endif
#ifdef DFS
   use dfsvar, only    : si=>sigmafull,del=>delsig,sl=>sigma,iope,ib,jbw
#endif
#ifndef RMP
   use comfgrid
   use comfphys
   use comfver
   use comznl
   use comgpd
   use comio
#else
   use rscomf_rerun
   use rscomltb
   use rscommap
   use rscomgpd
#ifdef MP
   use commpi
#endif
#endif
   use comsfc
#ifdef EXPLICIT_CLOUDINESS
   use comfcst, only                 : qcicps, qrscps
#endif
#ifdef DG3
   use comgda
   use diag_3d_module, only          : diag_3d_archive
#endif
#ifdef RIVER
   use comfrivh
#endif
#ifdef OMLWRF
   use module_ocean_slab_wrf
#endif
#ifdef WSM3
   use module_mp_wsm3
#endif
#ifdef WSM5
   use module_mp_wsm5
#endif
#ifdef WSM6
   use module_mp_wsm6
#endif
#ifdef WDM5
   use module_mp_wdm5
#endif
#ifdef WDM6
   use module_mp_wdm6
#endif
   use comconsts, only : levshcin,kpblmaxin
#ifdef DGP
   use diag_point_module, only : diag_point_arrange
#endif
#ifdef CHEM
   use dao_mod, only : c_eflux=>eflux, c_hflux=>hflux,   c_lwi=>lwi
   use dao_mod, only : c_pbl=>pbl,     c_preacc=>preacc, c_precon=>precon
   use dao_mod, only : c_ps1=>ps1,     c_snomas=>snomas, c_ts=>ts
   use dao_mod, only : c_u10m=>u10m,   c_ustar=>ustar,   c_v10m=>v10m
   use dao_mod, only : c_z0=>z0
   use dao_mod, only : c_cmfmc=>cmfmc, c_moistq=>moistq, c_sphu=>sphu
   use dao_mod, only : c_t=>t,         c_uwnd=>uwnd,     c_vwnd=>vwnd
#endif
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   real,parameter       ::  cnwatt=-cal_*1.e4/60.,epsq=2.e-12
   real,parameter       ::  critsnow=273.16
#ifdef OMLWRF
   real,parameter       ::  oml_gamma=0.14
#endif
#ifdef DGP
   integer              ::  iptlat(nptken)
   integer              ::  jptlat(nptken)
#endif
   integer              ::  lat,lons2
#ifndef RMP
   integer              ::  idate(4)
!
#if !defined (SMP) && !defined (NIM)
   real                 ::  pllamgr(LONF22S)
   real                 ::  plphigr(LONF22S)
   real                 ::  pslap(LONF22S)
#endif /* ~SMP and ~NIM end */
   real                 ::  plgr(LONF22S)
   real                 ::  uugrs(LONF22S,levs_)
   real                 ::  vvgrs(LONF22S,levs_)
   real                 ::  tvgrs(LONF22S,levs_)
   real                 ::  rqgrs(LONF22S,levh_)
#ifndef NONHYD
   real                 ::  xxgrs(LONF22S,levs_)
#else
   real                 ::  wwgrs(LONF22S,levp1_)
#endif
   real                 ::  ggt0(LONF22S,levs_)
   real                 ::  ggq0(LONF22S,levh_)
   real                 ::  ggu0(LONF22S,levs_)
   real                 ::  ggv0(LONF22S,levs_)
#ifdef NIM
#ifdef NIM_DIAG
   real*8               ::  evap_nim(LONF22S)
   real*8               ::  shflx_nim(LONF22S)
   real*8               ::  lhflx_nim(LONF22S)
   real*8               ::  dlwsfc_nim(LONF22S)
   real*8               ::  ulwsfc_nim(LONF22S)
   real*8               ::  rnet_nim(LONF22S)
   real*8               ::  swhr_nim(LONF22S,levs_)
   real*8               ::  lwhr_nim(LONF22S,levs_)
   real*8               ::  ttend_nim(LONF22S,levs_)
   real*8               ::  tots_pbl(LONF22S)
   real*8               ::  toth_pbl(LONF22S)
#endif
#endif
#else /* RMP start */
#ifdef PSPLIT
   real                 ::  ggt0(ILOTS,levs_)
   real                 ::  ggq0(ILOTS,levh_)
   real                 ::  ggu0(ILOTS,levs_)
   real                 ::  ggv0(ILOTS,levs_)
#endif
#endif /* ~RMP end */
#ifdef NIM_DIAG
   real                 ::  xmu(ILOTS)
   real                 ::  xmu_nim(ILOTS)
#endif
#ifdef YSUTKE
#ifndef RMP
   real                 ::  corf(ILOTS)
#else
   real                 ::  rcorf(ILOTS)
#endif
#endif
#ifndef NONHYD
#ifndef SMP
   real                 ::  plamgr(ILOTS)
   real                 ::  pphigr(ILOTS)
#endif
#endif /* ~NONHYD end */
!
! warning: gt0,gq0,gu0,gv0 may overlay tgrs,qgrs,ugrs,vgrs respectively.
! use local arrays to save the inputs
!
#ifdef CPS_KSAS_DIUR
   real :: hpbl_hold(ILOTS)
#endif
   real                 ::  ugrs(ILOTS,levs_)
   real                 ::  vgrs(ILOTS,levs_)
   real                 ::  pgr(ILOTS)
   real                 ::  tgrs(ILOTS,levs_)
   real                 ::  qgrs(ILOTS,levh_)
#ifdef YSUTKE
#ifdef RMP
   real                 ::  qgrs_tke(ILOTS,levs_)
#endif
#endif
#ifndef NONHYD
   real                 ::  xgrs(ILOTS,levs_)
#endif
   real                 ::  gt0(ILOTS,levs_)
   real                 ::  gq0(ILOTS,levh_)
   real                 ::  gu0(ILOTS,levs_)
   real                 ::  gv0(ILOTS,levs_)
#ifdef CPS_ENS
   real                 ::  gt0i(ILOTS,levs_)
   real                 ::  gq0i(ILOTS,levh_)
   real                 ::  gu0i(ILOTS,levs_)
   real                 ::  gv0i(ILOTS,levs_)
!
   real                 ::  gtsas(ILOTS,levs_),gtras(ILOTS,levs_)
   real                 ::  gtccm(ILOTS,levs_),gtkf2(ILOTS,levs_)
   real                 ::  gtkuo(ILOTS,levs_)
   real                 ::  gqsas(ILOTS,levh_),gqras(ILOTS,levh_)
   real                 ::  gqccm(ILOTS,levh_),gqkf2(ILOTS,levh_)
   real                 ::  gqkuo(ILOTS,levh_)
   real                 ::  gusas(ILOTS,levs_),guras(ILOTS,levs_)
   real                 ::  guccm(ILOTS,levs_),gukf2(ILOTS,levs_)
   real                 ::  gukuo(ILOTS,levs_)
   real                 ::  gvsas(ILOTS,levs_),gvras(ILOTS,levs_)
   real                 ::  gvccm(ILOTS,levs_),gvkf2(ILOTS,levs_)
   real                 ::  gvkuo(ILOTS,levs_)
   real                 ::  cld1dsas(ILOTS),cld1dras(ILOTS)
   real                 ::  cld1dccm(ILOTS),cld1dkuo(ILOTS)
   real                 ::  cld1dkf2(ILOTS)
   real                 ::  rain1sas(ILOTS),rain1ras(ILOTS)
   real                 ::  rain1ccm(ILOTS),rain1kuo(ILOTS)
   real                 ::  rain1kf2(ILOTS)
   integer              ::  icpssas(ILOTS),icpsras(ILOTS)
   integer              ::  icpsccm(ILOTS),icpskuo(ILOTS)
   integer              ::  icpskf2(ILOTS)
   integer              ::  kbotsas(ILOTS),kbotras(ILOTS)
   integer              ::  kbotccm(ILOTS),kbotkuo(ILOTS)
   integer              ::  kbotkf2(ILOTS)
   integer              ::  ktopsas(ILOTS),ktopras(ILOTS)
   integer              ::  ktopccm(ILOTS),ktopkuo(ILOTS)
   integer              ::  ktopkf2(ILOTS)
#endif
#ifdef NONHYD_HYD
   real                 ::  pgrs(ILOTS,levs_)
   real                 ::  wgrs(ILOTS,levp1_)
   real                 ::   gp0(ILOTS,levs_)
   real                 ::   gw0(ILOTS,levp1_)
   real                 ::  tgrh(ILOTS,levs_)
#endif /* NONHYD_HYD */
#ifdef NIM
   integer              ::  nits,ips,ipe
   real                 ::  wgrs(ILOTS,levp1_)
   real                 ::  si1d(levp1_)
   real                 ::  sl1d(levs_)
   real                 ::  prsl_nim(ILOTS,levs_)
   real                 ::  prsi_nim(ILOTS,levp1_)
   real                 ::  nim_z(ILOTS,levp1_)
   real                 ::  nim_zm(ILOTS,levs_)
#endif
   real                 ::  count
   real                 ::  slk(levs_)
   real                 ::  sik(levp1_)
   integer              ::  kbot(ILOTS)
   integer              ::  ktop(ILOTS)
   integer              ::  icps(ILOTS)
   real                 ::  dudt(ILOTS,levs_)
   real                 ::  dvdt(ILOTS,levs_)
   real                 ::  dtdt(ILOTS,levs_)
   real                 ::  dqdt(ILOTS,levh_)
#ifdef NONHYD_HYD
   real                 ::  dwdt(ILOTS,levp1_)
#endif /* NONHYD_HYD */
   real                 ::  hsw(ILOTS,levs_)
   real                 ::  gflx(ILOTS)
   real                 ::  rain(ILOTS)
   real                 ::  rainc(ILOTS)
   real                 ::  rainl(ILOTS)
   real                 ::  rain1(ILOTS)
   real                 ::  rain2(ILOTS)
   real                 ::  snow1(ILOTS)
   real                 ::  graupel1(ILOTS)
   real                 ::  snowncv(ILOTS)
   real                 ::  graupelncv(ILOTS)
   real                 ::  rainncv(ILOTS)
   real                 ::  sr(ILOTS)
   real                 ::  evapc(ILOTS)
   real                 ::  work1(ILOTS)
   real                 ::  temp
   real                 ::  temp_i(ILOTS)
   real                 ::  temp_ik(ILOTS,levs_)
   real                 ::  wscalek(ILOTS,levs_)
   real                 ::  wstar(ILOTS)
   real                 ::  delta(ILOTS)
   real                 ::  wind(ILOTS)
   real                 ::  w2(ILOTS,levs_)
   real                 ::  vvel(ILOTS,levs_)
#ifdef KSAS
   real                 ::  vvel_pa(ILOTS,levs_)
#endif
   real                 ::  prsi(ILOTS,levs_+1)
   real                 ::  prsi_pa(ILOTS,levs_+1)
   real                 ::  prsik(ILOTS,levs_+1)
   real                 ::  prsl(ILOTS,levs_)
   real                 ::  prsl_pa(ILOTS,levs_)
   real                 ::  prslk(ILOTS,levs_)
   real                 ::  phii(ILOTS,levs_+1)
   real                 ::  phil(ILOTS,levs_)
   real                 ::  zl(ILOTS,levs_)
   real                 ::  zi(ILOTS,levs_+1)
#ifdef YSUTKE
   real                 ::  el_pbl(ILOTS,levs_)
#endif
   real                 ::  delz(ILOTS,levs_)
   real                 ::  delprsi(ILOTS,levs_)
   real                 ::  delprsi_pa(ILOTS,levs_)
   real                 ::  ww(ILOTS,levs_)
   real                 ::  dtaux(ILOTS,levs_)
   real                 ::  dtauy(ILOTS,levs_)
!
#ifdef HYBRID
   real                 ::  sihyb(ILOTS,levp1_),cihyb(ILOTS,levp1_)
   real                 ::  slhyb(ILOTS,levs_) ,clhyb(ILOTS,levs_)
   real                 ::  delhyb(ILOTS,levs_),slkhyb(ILOTS,levs_)
#endif
   real                 ::  tstar(ILOTS) ,psexp(ILOTS) ,psfcpa(ILOTS)
   real                 ::  snowmt(ILOTS),snowev(ILOTS),snowfl(ILOTS)
   real                 ::  fm(ILOTS)    ,fh(ILOTS)
   real                 ::  cd(ILOTS)    ,cdq(ILOTS)
   real                 ::  qss(ILOTS)   ,radsl(ILOTS)
   real                 ::  dusfcg(ILOTS),dvsfcg(ILOTS)
#ifdef GWDC
   real                 ::  tauctx(ILOTS),taucty(ILOTS)
#endif
   real                 ::  dusfc1(ILOTS),dvsfc1(ILOTS)
   real                 ::  dtsfc1(ILOTS),dqsfc1(ILOTS)
   real                 ::  dlwsf1(ILOTS),ulwsf1(ILOTS)
#ifdef OMLWRF
   real                 ::  nswsf1(ILOTS)
#endif
#ifdef VIC
   real                 ::  dswsf1(ILOTS)
#endif
!
   real                 ::  smsoil(ILOTS,lsoil_)
#ifndef VIC
   real                 ::  stsoil(ILOTS,lsoil_)
#else
   real                 ::  stsoil(ILOTS,nsoil_)
#endif
   integer              ::  soiltyp(ILOTS)
   integer              ::  vegtype(ILOTS)
   integer              ::  kpbl(ILOTS)
   real                 ::  sigmaf(ILOTS)
   real                 ::  rb(ILOTS)
   real                 ::  rhscnpy(ILOTS)
   real                 ::  ai(ILOTS,lsoil_)
   real                 ::  bi(ILOTS,lsoil_)
   real                 ::  drain(ILOTS)
   real                 ::  cci(ILOTS,lsoil_)
   real                 ::  rhsmc(ILOTS,lsoil_)
   real                 ::  runof(ILOTS)
   real                 ::  zsoil(ILOTS,lsoil_)
   real                 ::  cld1d(ILOTS)
   real                 ::  evap(ILOTS)
   real                 ::  evapol(ILOTS)
   real                 ::  hflx(ILOTS)
   real                 ::  rnet(ILOTS)
   real                 ::  t850(ILOTS)
   real                 ::  pk(ILOTS)
   real                 ::  pkp(ILOTS)
   real                 ::  ep1d(ILOTS)
   real                 ::  gamt(ILOTS)
   real                 ::  gamq(ILOTS)
   real                 ::  var(ILOTS)
   real                 ::  oc(ILOTS)
   real                 ::  oa4(ILOTS,4)
   real                 ::  ol4(ILOTS,4)
#ifdef FLOW_BLOCKING
   real                 ::  omax(ILOTS)
#endif
   real                 ::  ustar(ILOTS)
#ifdef HYDRO
   real                 ::  hydrow(ILOTS,levs_)
   real                 ::  evcn(ILOTS)
#else
#ifdef NOAHYDRO
   real                 ::  evcn(ILOTS)
#endif
#endif
   real                 ::  qu00(ILOTS),qv00(ILOTS)
#ifdef VIC
   real                 ::  evp(ILOTS),heat(ILOTS)
   real                 ::  rho(ILOTS)
#ifdef VICLSM1
   real                 ::  ght(ILOTS)
#endif
#endif
#ifdef NOALSM1
   integer              ::  slptyp(ILOTS)
   real                 ::  slsoil(ILOTS,lsoil_)
#endif
   real                 ::  fm10   (ILOTS),fh2      (ILOTS)
#ifndef OSU 
   real                 ::  oldprcp(ILOTS),oldsrflag(ILOTS)
#endif
#if defined(VICLSM1) || defined(VICLSM2)
   real                 ::  expld(ILOTS,lsoil_)
   real                 ::  dphld(ILOTS,lsoil_),bubld(ILOTS,lsoil_)
   real                 ::  qrtld(ILOTS,lsoil_),bkdld(ILOTS,lsoil_)
   real                 ::  sldld(ILOTS,lsoil_),wcrld(ILOTS,lsoil_)
   real                 ::  wpwld(ILOTS,lsoil_),smrld(ILOTS,lsoil_)
   real                 ::  smxld(ILOTS,lsoil_)
   real                 ::  dpnld(ILOTS,nsoil_),sxnld(ILOTS,nsoil_)
   real                 ::  epnld(ILOTS,nsoil_),bbnld(ILOTS,nsoil_)
   real                 ::  apnld(ILOTS,nsoil_),btnld(ILOTS,nsoil_)
   real                 ::  gmnld(ILOTS,nsoil_)
   real                 ::  kstld(ILOTS,lsoil_)
#ifdef VICLSM1
   real                 ::  sicld(ILOTS,lsoil_)
   real                 ::  vrtld(ILOTS,lsoil_)
   integer              ::  lstsnld(ILOTS)
#endif
#ifdef VICLSM2
   real                 ::  flaild(ILOTS,msub_)
   real                 ::  vfrld(ILOTS,msub_), vtyld(ILOTS,msub_)
   real                 ::  cnpld(ILOTS,msub_), snold(ILOTS,msub_)
   real                 ::  csnld(ILOTS,msub_), rsnld(ILOTS,msub_)
   real                 ::  tsfld(ILOTS,msub_), tpkld(ILOTS,msub_)
   real                 ::  sfwld(ILOTS,msub_), pkwld(ILOTS,msub_)
!
   real                 ::  sicld(ILOTS,lsoil_,msub_)
   real                 ::  smcld(ILOTS,lsoil_,msub_)
   real                 ::  stcld(ILOTS,nsoil_,msub_)
!
   real                 ::  vrtld(ILOTS,lsoil_,msub_)
   integer              ::  lstsnld(ILOTS,msub_), nvegld(ILOTS)
#endif
#endif /* VICLSM1 or VICLMS2 end */
#ifdef WSM2
   real                 ::  tp(ILOTS,levs_,LATG2S),tp1(ILOTS,levs_,LATG2S),    &
                            qp(ILOTS,levs_,LATG2S),qp1(ILOTS,levs_,LATG2S),    &
                            psp(ILOTS,LATG2S),psp1(ILOTS,LATG2S)
#endif
#ifdef SMP
   integer              :: jcapset
#endif
#ifdef DG
   real                 ::  tgmxl,tgmnl
   integer              ::  ismax, ismin
   integer              ::  igmxl,kgmxl,igmnl,kgmnl
#endif
#ifdef DG3
   real                 ::  gda(nwgda,kdgda)
#endif
#ifdef CLM_CWF
#ifdef DFS
   real                 ::  cgs2(lons2,levs_)
#else
   real                 ::  cgs2(ILOTS,levs_)
#endif
#endif
#ifdef MUL_CLDTOP
   real                 ::  xkt2(ILOTS*NCLDTOP)
#endif
#ifdef SAS_DIAG
   real                 ::  dcu(ILOTS,levs_),dcv(ILOTS,levs_)
   real                 ::  dct(ILOTS,levs_)
   real                 ::  dcq(ILOTS,levs_),dch(ILOTS,levs_)
   real                 ::  fcu(ILOTS,levs_),fcd(ILOTS,levs_)
   real                 ::  deltb(ILOTS),delqb(ILOTS),delhb(ILOTS)
   real                 ::  cbmf(ILOTS)
   real                 ::  dlt(ILOTS,levs_),dlq(ILOTS,levs_)
   real                 ::  dlh(ILOTS,levs_)
#else
#if defined GWDC || defined KSAS
   real                 ::  dct(ILOTS,levs_)
   real                 ::  fcu(ILOTS,levs_),fcd(ILOTS,levs_)
#endif
#endif
#ifdef RAS
   logical              ::  ras
   integer              ::  lmx,lm,krmin,krmax,nstrp,ncrnd,il1r,il2r,il3r
   real                 ::  cp,alhl,grav,rgas,afac,ufac,dtras
   real                 ::  sig(lmx+1),prj(lmx+1),prh(lmx),fpk(lmx)
   real                 ::  hpk(lmx),sgb(lmx),ods(lmx),rasal(lmx),prns(lmx/2)
   real                 ::  rannum(200)
#endif
#ifdef RASV2
   logical              ::  ras
   real                 ::  rgas,cp,grav,alhl,dtras
   real                 ::  sig(levs_+1),prj(levs_+1),sgb(levs_)
   real                 ::  rasal(levs_)
   real                 ::  rannum(ncrnd),dsfc(ILOTS),pdd
   integer              ::  krmin, krmax,kfmax,ncrnd,mct,kctop(mct+1)
   real                 ::  gt00(levs_),gq00(levs_),gu00(levs_),gv00(levs_)
   real                 ::  clw0(levs_),clt(levs_),clw(ILOTS,levs_)
#endif
#ifdef EXPLICIT_CLOUDINESS
   real                 ::  qci(ILOTS,levs_),qrs(ILOTS,levs_)
#endif
#ifdef CLW_DET
   real                 ::  clw(ILOTS,levs_)
#endif
#ifdef SMP_NUDGING_FLX
   real                 ::  scmsh(LONF2S),scmlh(LONF2S)
#endif
#ifdef SMP_NUDGING_RAD
   real                 ::  psfc(LONF2S),tsfc(LONF2S),qsfc(LONF2S)
   real                 ::  prec(LONF2S)
   real                 ::  scmlwup(LONF2S),scmlwdn(LONF2S)
   real                 ::  scmswup(LONF2S),scmswdn(LONF2S)
   real                 ::  scmalb(LONF2S)
#endif
!
   integer              ::  ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
   integer              ::  ktotal,nvdiff,levshc,i,ic,k,kc,jcapr,n,kpblmax
#ifdef DBG
   integer              ::  latd,lond,level
#endif
   real                 ::  coriolis,frain,dtf,xscale,dxmeter,rcs
   real                 ::  sigshc,p850
#ifndef RMP
   real                 ::  dt2,rcl
#endif
#ifdef RMP
   real                 ::  rbs2
#endif
   logical              ::  ifvmix(ntotal_)
!-------------------------------------------------------------------------------
!
! define the LONS2 i-grid vector number
!
#ifndef RMP
#define LONS2 lons2
#else
#ifdef MP
#ifdef RMPVECTORIZE
#define LONS2 ILOTS
#else
#define LONS2 lonlen(mype)*2
#endif
#else
#define LONS2 igrd12_
#endif
#endif
!
! initialize local variables
!
#ifndef RMP
   do i = 1,ILOTS
#if !defined(SMP) && !defined(NIM)
     plamgr(i)=0.
     pphigr(i)=0.
#endif
     pgr(i)=0.
   enddo
!
#endif   /* ~RMP */
   do k = 1,levs_
     do i = 1,ILOTS
#ifndef RMP 
#if !defined(SMP) && !defined(NONHYD)
       xgrs(i,k)=0.
#endif
       ugrs(i,k)=0.   ;  gu0(i,k)=0.
       vgrs(i,k)=0.   ;  gv0(i,k)=0.
       tgrs(i,k)=0.   ;  gt0(i,k)=0.
#endif /* ~RMP end */
       dudt(i,k)=0.   ;  dvdt(i,k)=0.     ;  dtdt(i,k)=0.
       hsw(i,k)=0.
       temp_ik(i,k)=0.
       w2(i,k)=0.
       vvel(i,k)=0.
       prsl(i,k)=0.   ;  prsl_pa(i,k)=0.  ;  prslk(i,k)=0.
       phil(i,k)=0.   ;  zl(i,k)=0.       ;  delz(i,k)=0.
       ww(i,k)=0.
       delprsi(i,k)=0.;  delprsi_pa(i,k)=0.
       dtaux(i,k)=0.  ;  dtauy(i,k)=0.
#ifdef HYBRID
       slhyb(i,k)=0.  ;  clhyb(i,k)=0.
       delhyb(i,k)=0. ;  slkhyb(i,k)=0.
#endif
#ifdef HYDRO
       hydrow(i,k)=0.
#endif
#if defined(CLW_DET) || defined(RASV2)
       clw(i,k)=0.
#endif
#ifdef CPS_ENS
       gt0i(i,k)=0. ; gu0i(i,k)=0. ; gv0i(i,k)=0.
       gtsas(i,k)=0.; gtras(i,k)=0.; gtccm(i,k)=0.; gtkf2(i,k)=0.; gtkuo(i,k)=0.
       gusas(i,k)=0.; guras(i,k)=0.; guccm(i,k)=0.; gukf2(i,k)=0.; gukuo(i,k)=0.
       gvsas(i,k)=0.; gvras(i,k)=0.; gvccm(i,k)=0.; gvkf2(i,k)=0.; gvkuo(i,k)=0.
#endif
#ifdef YSUTKE
       el_pbl(i,k)=0.
#endif
     enddo
     slk(k)=0. 
#ifdef RASV2
     gt00(k)=0.       ;  gq00(k)=0.      ;  gu00(k)=0     ;  gv00(k)=0.
     clw0(k)=0.       ;  clt(k)=0.
#endif 
   enddo
!
   do k = 1,levh_
     do i = 1,ILOTS
#ifndef RMP
       qgrs(i,k)=0.   ;  gq0(i,k)=0.
#endif
       dqdt(i,k)=0.
#ifdef CPS_ENS
       gq0i(i,k)=0.
       gqsas(i,k)=0.; gqras(i,k)=0.; gqccm(i,k)=0.; gqkf2(i,k)=0.; gqkuo(i,k)=0.
#endif
     enddo
   enddo
!
! just i
!
   temp=0.
   do i = 1,ILOTS
     kbot(i)=0        ;  ktop(i)=0        ;  icps(i)=0
     gflx(i)=0.
     rain(i)=0.       ;  rainc(i)=0.      ;  rainl(i)=0.
     rain1(i)=0.      ;  rain2(i)=0.      ;  snow1(i)=0.
     graupel1(i)=0.   ;  snowncv(i)=0.
     graupelncv(i)=0. ;  rainncv(i)=0.
     sr(i)=0.         ;  evapc(i)=0.      ;  work1(i)=0.
     temp_i(i)=0.
     wind(i)=0.
     tstar(i)=0.      ;  psexp(i)=0.      ;  psfcpa(i)=0.
     snowmt(i)=0.     ;  snowev(i)=0.     ;  snowfl(i)=0.  
     fm(i)=0.         ;  fh(i)=0.         ;  cd(i)=0.      ;  cdq(i)=0.
     qss(i)=0.        ;  radsl(i)=0.      ;  dusfcg(i)=0.  ;  dvsfcg(i)=0.
#ifdef GWDC
     tauctx(i)=0.     ;  taucty(i)=0.
#endif
     dusfc1(i)=0.     ;  dvsfc1(i)=0.     ;  dtsfc1(i)=0.
     dqsfc1(i)=0.     ;  dlwsf1(i)=0.     ;  ulwsf1(i)=0.
#ifdef OMLWRF
     nswsf1(i)=0.
#endif
#ifdef VIC
     dswsf1(i)=0.
#endif
     soiltyp(i)=0     ;  vegtype(i)=0 
     kpbl(i)=0        ;  sigmaf(i)=0.     ;  rb(i)=0.      ; rhscnpy(i)=0.
     drain(i)=0.
     runof(i)=0.
     cld1d(i)=0.      ;  evap(i)=0.       ;  evapol(i)=0.
     hflx(i)=0.       ;  rnet(i)=0.       ;  t850(i)=0.
     pk(i)=0.         ;  pkp(i)=0.
     ep1d(i)=0.       ;  gamt(i)=0.       ;  gamq(i)=0.
     var(i)=0.        ;  oc(i)=0.         ;  ustar(i)=0.
#ifdef FLOW_BLOCKING
     omax(i)=0.
#endif
#ifdef HYDRO
     evcn(i)=0.
#else
#ifdef NOAHYDRO
     evcn(i)=0.
#endif
#endif
     qu00(i)=0.       ;  qv00(i)=0.
#ifdef VIC
     evp(i)=0.        ;  heat(i)=0.       ;  rho(i)=0.
#ifdef VICLSM1
     ght(i)=0.
#endif
#endif
#ifdef NOALSM1
     slptyp(i)=0
#endif
     fm10(i)=0.       ;  fh2(i)=0.
#ifndef OSU 
     oldprcp(i)=0.    ;  oldsrflag(i)=0.
#endif
#ifdef VICLSM1
     lstsnld(i)=0 
#endif
#ifdef VICLSM2
     nvegld(i)=0
#endif
#ifdef CPS_ENS
     cld1dsas(i)=0.;cld1dras(i)=0.;cld1dccm(i)=0.;cld1dkuo(i)=0.;cld1dkf2(i)=0.
     rain1sas(i)=0.;rain1ras(i)=0.;rain1ccm(i)=0.;rain1kuo(i)=0.;rain1kf2(i)=0.
     icpssas(i)=0  ;icpsras(i)=0  ;icpsccm(i)=0  ;icpskuo(i)=0  ;icpskf2(i)=0
     kbotsas(i)=0  ;kbotras(i)=0  ;kbotccm(i)=0  ;kbotkuo(i)=0  ;kbotkf2(i)=0
     ktopsas(i)=0  ;ktopras(i)=0  ;ktopccm(i)=0  ;ktopkuo(i)=0  ;ktopkf2(i)=0
#endif
#ifdef GRIMSCV
     wstar(i)=0.
     delta(i)=0.
#endif
#ifdef NIM_DIAG
     xmu(i)=0.        ;  xmu_nim(i)=0.
#endif
#ifdef YSUTKE
#ifndef RMP
     corf(i)=0.
#else
     rcorf(i)=0.
#endif
#endif
   enddo
!
! levs+1
!
   do k = 1,levp1_
     do i = 1,ILOTS
#ifdef NONHYD_HYD
       dwdt(i,k)=0.
#endif
       prsi(i,k)=0.   ;  prsi_pa(i,k)=0.  ;  prsik(i,k)=0.
       phii(i,k)=0.   ;  zi(i,k)=0.
#ifdef HYBRID
       sihyb(i,k)=0.  ;  cihyb(i,k)=0.
#endif
     enddo
     sik(k)=0.
   enddo
!
! lsoil
!
   do k = 1,lsoil_
     do i = 1,ILOTS
       smsoil(i,k)=0. ;  ai(i,k)=0.       ;  bi(i,k)=0.
#ifndef VIC
       stsoil(i,k)=0.
#endif
       cci(i,k)=0.    ;  rhsmc(i,k)=0.    ;  zsoil(i,k)=0.
#ifdef NOALSM1
       slsoil(i,k)=0.
#endif
#if defined(VICLSM1) || defined(VICLSM2)
       expld(i,k)=0.  ;  dphld(i,k)=0.    ;  bubld(i,k)=0.
       qrtld(i,k)=0.  ;  bkdld(i,k)=0.
       sldld(i,k)=0.  ;  wcrld(i,k)=0.    ;  wpwld(i,k)=0.
       smrld(i,k)=0.  ;  smxld(i,k)=0.    :  kstld(i,k)=0.
#ifdef VICLSM1
       sicld(i,k)=0.  ;  vrtld(i,k)=0.
#endif
#endif
     enddo
   enddo
!
! nsoil
!
#if defined(VICLSM1) || defined(VICLSM2)
   do k = 1,nsoil_
     do i = 1,ILOTS
       dpnld(i,k)=0.  ;  sxnld(i,k)=0.    ;  epnld(i,k)=0.
       bbnld(i,k)=0.  ;  apnld(i,k)=0.
       btnld(i,k)=0.  ;  gmnld(i,k)=0.
#ifdef VIC
       stsoil(i,k)=0.
#endif
     enddo
   enddo
!
#endif
   do k = 1,4
     do i = 1,ILOTS
       oa4(i,k)=0.    ;  ol4(i,k)=0.
     enddo
   enddo
!
#ifdef VICLSM2
! msub
!
   do k = 1,msub_
     do i = 1,ILOTS
       flaild(i,k)=0. ;  vfrld(i,k)=0.    ;  vtyld(i,k)=0.
       cnpld(i,k)=0.  ;  snold(i,k)=0.
       csnld(i,k)=0.  ;  rsnld(i,k)=0.    ;  tsfld(i,k)=0.
       tpkld(i,k)=0.  ;  sfwld(i,k)=0.    ;  pkwld(i,k)=0.
       lstsnld(i,k)=0
     enddo
     do j = 1,lsoil_
       do i = 1,ILOTS
         sicld(i,j,k)=0.  ;  smcld(i,j,k)=0.
         vrtld(i,j,k)=0.
       enddo
     enddo
     do j = 1,nsoil_
       do i = 1,ILOTS
         stcld(i,j,k)=0.
       enddo
     enddo
   enddo
!
#endif
#ifdef WSM2
! latg2s
!
   do k = 1,LATG2S
     do j = 1,levs_
       do i = 1,ILOTS
         tp(i,j,k)=0.   ;  tp1(i,j,k)=0.   ;  qp(i,j,k)=0.
         qp1(i,j,k)=0.
       enddo
     enddo
     do i = 1,ILOTS
       psp(i,k)=0.
       psp1(i,k)=0.
     enddo
   enddo
!
#endif
#ifdef NIM_DIAG
   do i = 1,LONF22S
     evap_nim(i)=0.      ; shflx_nim(i)=0.     ; lhflx_nim(i)=0.
     dlwsfc_nim(i)=0.    ; ulwsfc_nim(i)=0.    ; rnet_nim(i)=0.
     tots_pbl(i)=0.      ; toth_pbl(i)=0.
   enddo
!
   do k = 1,levs_
     do i = 1,LONF22S
       swhr_nim(i,k)=0.  ; lwhr_nim(i,k)=0.
       ttend_nim(i,k)=0.
     enddo
   enddo
#endif
!
! define the dimension constants
!
   ids = 1
   ide = ILOTS
   jds = 1
   jde = 1
   kds = 1
   kde = levs_
!
   ims = 1
   ime = ILOTS 
   jms = 1
   jme = 1
   kms = 1
   kme = levs_
!
   its = 1
   ite = LONS2
   jts = 1
   jte = 1
   kts = 1
   kte = levs_
!
! define the water and gas substance dimension constants
!
   ktotal = kte*ntotal_
!
   nvdiff = ntotal_
   ifvmix(:) = .true.
#ifdef YSUTKE
   ifvmix(ibltke_) = .false.
#endif
   if(nwater_.ge.3) then
     do ic = 3, nwater_
       if(ic.ne.4) ifvmix(ic) = .false.
     enddo
   endif
!
!  time step constants
!
   frain = .5
#ifndef NIM
!soojin_couple
   if(inistp.eq.1) frain = 0.25
!   if(inistp.eq.1) frain = 1.
#endif
   dt2    = deltim*2.e0
   dtf    = frain*dt2
#ifndef RMP
   dtf    = dtf*count
#endif
!
!  define physics related constants
!
#ifndef SMP
#ifndef NIM
#ifdef RMP
   rcl     = 1.
   rbs2    = 1.0
   jcapr   = 126.*106000./ rdelx
   dxmeter = rdelx
#else
   rcl     = rbs2(lat)
   jcapr   = jcap_
   dxmeter = 12600000./jcap_
#endif
#else
   rcl     = 1.
   xscale  = 2.
   jcapr   = exp(log(xscale)*(glvl-5))*62.      ! scale for SAS
   dxmeter = 7071000./exp(log(xscale)*(glvl))   ! dx
#endif
#else
   rcl     = rbs2(lat)
   jcapset = 126
   dxmeter = 12600000./jcapset
#endif
   rcs = sqrt(rcl)
!
#ifdef YSUTKE
   do i = its,ite
#ifndef RMP
     corf(i) = 2.*omega_*sinlab(i,lat)
#else
     rcorf(i) = 2.*omega_*sinlar(i,lat)
#endif
   enddo
#endif
#ifdef FEB2013_GRIMSCV
   sigshc = 0.6
#else
   sigshc = 0.7
#endif
   levshc  = levshcin
   kpblmax = kpblmaxin
   p850    = 850. * .1
!
#ifdef NIM
   kdt=nits
!
   do k = 1,levs_
     si(k)=si1d(k)
     sl(k)=sl1d(k)
   enddo
!
   si(levp1_)=si1d(levp1_)
#endif
   do k = kts,kte
     slk(k)    = sl(k)**(akapa_)
   enddo
!
   do k = kts,kte+1
     sik(k)    = si(k)**(akapa_)
   enddo
#ifdef DBG
!
!  define the  i j points for printing
!
   latd = -1
   lond = -1
#ifdef OPENMP
!$omp master
#endif
   if (iope ) then
     latd = 1
     lond = 10
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifndef RMP
!
!   assign the input variables without changes
!
   do i = its,ite
#if !defined(SMP) && !defined(NIM)
     plamgr(i) = pllamgr(i)
     pphigr(i) = plphigr(i)
#endif
     pgr(i) = plgr(i)
   enddo
!
   do k = kts,kte
     do i = its,ite
       ugrs(i,k) = uugrs(i,k)
       vgrs(i,k) = vvgrs(i,k)
       tgrs(i,k) = tvgrs(i,k)
       qgrs(i,k) = rqgrs(i,k)
#ifndef NONHYD
       xgrs(i,k) = xxgrs(i,k)
#endif
     enddo
   enddo
#ifdef NISLQ_CHECK
!
! check negative before phys_main_solver
!
   do k = 1,levh_
     do i = its,ite
       if(qgrs(i,k).lt.0.) print *,'negative before phys_main_solver',k,qgrs(i,k)
     enddo
   enddo
#endif /* NISLQ_CHECK */
#ifdef NONHYD
   do k = kts,kte+1
     do i = its,ite
          wgrs(i,k) = wwgrs(i,k)
     enddo
   enddo
#endif
!
   if(ntotal_.gt.1) then
     do k = kcloud_,levh_
       do i = its,ite
         qgrs(i,k) = max(rqgrs(i,k),qmin_)
       enddo
     enddo
   endif
!
    if(ngases_.ge.1.and.kdt.le.2) then
#ifdef DBG
!
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       if(lat.eq.latd) then
         call print_maxmin_seven(qgrs(1,kgases_),ite,ime,kte,kts,kte,          &
                                'gases before')
       endif
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
     do k = kts,kte
       do i = its,ite
#ifndef AQUA_PLANET
         qgrs(i,k+kgases_-1) = ozon(i,k,lat)
#else
         qgrs(i,k+kgases_-1) = 0.0
#endif
       enddo
     enddo
!
#ifdef WDM
     do k = kts,kte
       do i = its,ite
         qgrs(i,k+kwsize_-1) = iccn_
       enddo
     enddo
!
#endif
#ifdef YSUPBL
#ifdef YSUTKE
     do k = kts,kte
       do i = its,ite
         qgrs(i,k+kbltke_-1) = 0.01
       enddo
     enddo
#endif
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       if(lat.eq.latd)                                                         &
         call print_maxmin_seven(qgrs(1,kgases_),ite,ime,kte,kts,kte,          &
                                'gases after')
     endif
#ifdef OPENMP
!$omp end master
#endif
!
#endif
   endif
#else    /* RMP */
!
!   reserve the input without change
!
   if(ntotal_.gt.1) then
     do k = kcloud_,levh_
       do i = its,ite
         qgrs(i,k)=max(qgrs(i,k),qmin_)
       enddo
     enddo
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       if(lat.eq.latd)                                                         &
         call print_maxmin_seven(qgrs(1,kgases_),ite,ime,kte,kts,kte,          &
                                'gases before')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
   endif
!
   if(inistp.eq.1.and.ngases_.ge.1) then
     do k = kts,kte
       do i = its,ite
         qgrs(i,k+kgases_-1) = ozon(i,k,lat)
       enddo
     enddo
!
#ifdef WDM
     do k = kts,kte
       do i = its,ite
         qgrs(i,k+kwsize_-1) = iccn_
       enddo
     enddo
!
#endif
#ifdef YSUPBL
#ifdef YSUTKE
     do k = kts,kte
       do i = its,ite
         qgrs(i,k+kbltke_-1) = 0.01
       enddo
     enddo
!
#endif
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       if(lat.eq.latd)                                                         &
         call print_maxmin_seven(qgrs(1,kgases_),ite,ime,kte,kts,kte,          &
                               'gases after')
     endif
#ifdef OPENMP
!$omp end master
#endif
!
#endif
   endif
#endif
!
!  transfer vegetation fraction from global to local location
!  as well as vegetation type and soil type
!
   do i = its,ite
#ifndef R2_PHYSICS
#ifdef OSULSM2
     soiltyp(i)= int(stype(i,lat)+.5)
#ifdef NOALSM1
     sigmaf(i) = max(vfrac(i,lat),0.)
#else
#ifdef NCAR_EDIR
     sigmaf(i) = max(vfrac(i,lat),0.)
#else
     sigmaf(i) = max(vfrac(i,lat),.3)
#endif
#endif  /* ifndef NOALSM1 */
#endif  /* ifdef OSULSM2 */
#ifdef NOA
     soiltyp(i)= stype(i,lat)+.5
#ifdef NOALSM1
     sigmaf(i) = max(vfrac(i,lat),0.)
#else
#ifdef NCAR_EDIR
     sigmaf(i) = max(vfrac(i,lat),0.)
#else
     sigmaf(i) = max(vfrac(i,lat),.3)
#endif
#endif  /* ifndef NOALSM1 */
#endif  /* ifdef NOA */
#ifdef VIC
     soiltyp(i) = 0.
     sigmaf(i) = max(vfrac(i,lat),0.)
#endif
#ifdef OSULSM1
     soiltyp(i) = 7
     sigmaf(i) = .7
#else
     vegtype(i)= int(vtype(i,lat)+.5)
     fm(i) = ffmm(i,lat)
     fh(i) = ffhh(i,lat)
     ustar(i) = uustar(i,lat)
#endif
#ifdef OSULSM2
!
     if(slmsk(i,lat).eq.2.) then
#ifdef STATSGO_SOIL
       vegtype(i) = 12
       soiltyp(i) = 16
#else
       vegtype(i) = 13
       soiltyp(i) = 9
#endif
     endif
!
#endif		/* OSULSM2 end */
#endif   	/* ifndef R2_PHYSICS  */
#ifdef R2_PHYSICS
     soiltyp(i) = 7
     sigmaf(i) = .7
#ifndef OSULSM1
     vegtype(i)= int(vtype(i,lat)+.5)
     fm(i) = ffmm(i,lat)
     fh(i) = ffhh(i,lat)
     ustar(i) = uustar(i,lat)
!
     if(slmsk(i,lat).eq.2.) then
       vegtype(i) = 13
       soiltyp(i) = 9
     endif
!
#endif		/* OSULSM1 */
#endif		/* R2_PHYSICS */
#ifdef NOALSM1
     slptyp(i) = slope(i,lat)+.5
     if(slmsk(i,lat).eq.2.) then
       slptyp(i) = 0
     endif
#endif
   enddo
!
!  transfer soil moisture and temperature from global to local variables
!
   do k = 1,lsoil_
     do i = its,ite
       smsoil(i,k) = smc(i,lat,k)
#ifndef VIC
       stsoil(i,k) = stc(i,lat,k)
#endif
#ifdef NOALSM1
       slsoil(i,k) = slc(i,lat,k)
#endif
#ifdef VIC
       expld(i,k) = expt(i,lat,k)
       kstld(i,k) = kest(i,lat,k)
       dphld(i,k) = dph(i,lat,k)
       bubld(i,k) = bub(i,lat,k)
       qrtld(i,k) = qrt(i,lat,k)
       bkdld(i,k) = bkd(i,lat,k)
       sldld(i,k) = sld(i,lat,k)
       wcrld(i,k) = wcr(i,lat,k)
       wpwld(i,k) = wpw(i,lat,k)
       smrld(i,k) = smr(i,lat,k)
       smxld(i,k) = smx(i,lat,k)
       vrtld(i,k) = vroot(i,lat,k)
#endif
#ifdef VICLSM1
       sicld(i,k) = sic(i,lat,k)
#endif
     enddo
   enddo
!
#ifdef VICLSM2
   do i = its,ite
     nvegld(i) = nveg(i,lat)
     do nv = 1,nvegld(i)
       flaild(i,nv) = flai(i,lat,nv)
       vfrld(i,nv)  = mvfr(i,lat,nv)      ! tile coverage
       vtyld(i,nv)  = mvty(i,lat,nv)
       cnpld(i,nv)  = mcnp(i,lat,nv)
       snold(i,nv)  = msno(i,lat,nv)      ! tile snow depth
       csnld(i,nv)  = csno(i,lat,nv)      ! canopy snow
       rsnld(i,nv)  = rsno(i,lat,nv)
       tsfld(i,nv)  = tsf(i,lat,nv)
       tpkld(i,nv)  = tpk(i,lat,nv)
       sfwld(i,nv)  = sfw(i,lat,nv)
       pkwld(i,nv)  = pkw(i,lat,nv)
       lstsnld(i,nv)= lstsn(i,lat,nv)
     enddo
   enddo
!
   do k = 1,lsoil_
     do i = its,ite
       do nv = 1,nvegld(i)
         vrtld(i,k,nv) = vroot(i,lat,(nv-1)*lsoil_+k)
         smcld(i,k,nv) = msmc(i,lat,(nv-1)*lsoil_+k)
         sicld(i,k,nv) = msic(i,lat,(nv-1)*lsoil_+k)
       enddo
     enddo
   enddo
!
#endif
#ifdef VIC
   do k = 1,nsoil_
     do i = its,ite
       stsoil(i,k)= stc(i,lat,k)
       dpnld(i,k) = dphn(i,lat,k)
       sxnld(i,k) = smxn(i,lat,k)
       epnld(i,k) = expn(i,lat,k)
       bbnld(i,k) = bubn(i,lat,k)
       apnld(i,k) = alpn(i,lat,k)
       btnld(i,k) = betn(i,lat,k)
       gmnld(i,k) = gamn(i,lat,k)
     enddo
   enddo
!
#endif
#ifdef VICLSM2
   stcld = 0.
   do k = 1,nsoil_
     do i = its,ite
       do nv = 1,nvegld(i)
         stcld(i,k,nv) = mstc(i,lat,(nv-1)*nsoil_+k)
       enddo
     enddo
   enddo
!
#endif
#ifdef VICLSM1
   lstsnld = 0.
   do i = its,ite
     lstsnld(i) = -1
   enddo
!
   do i = its,ite
     if(slmsk(i,lat).eq.1) then
       if(lstsn(i,lat).gt.1.e10) then
         lstsnld(i) = -1
       else
         lstsnld(i) = lstsn(i,lat)
       endif
     endif
   enddo
#endif
!
!  get dry temperature from virtual temperature
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
       call print_maxmin_seven(tgrs(1,1),ite,ime,kte,1,kte,'tgrs virtual')
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
   do k = kts,kte
     do i = its,ite
       w2(i,k) = 1.+fv_*qgrs(i,k)
       work1(i) = 0.
       if(ncloud_.ge.1) then
         do ic = icloud_,nwmass_
           kc = (ic-1)*kte + k
           work1(i) = max(qgrs(i,kc),qmin_) + work1(i)
         enddo
         w2(i,k) = w2(i,k) - work1(i)
       endif
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       tgrs(i,k) = tgrs(i,k)/w2(i,k)
#ifdef NONHYD_HYD
       tgrh(i,k) = tgrh(i,k)/w2(i,k)
#endif /* NONHYD_HYD */
     enddo
   enddo
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(tgrs(1,1),ite,ime,kte,1,kte,'tgrs dry')
   endif
#ifdef OPENMP
!$omp end master
#endif
!
#endif
#ifdef DG
   do k = kts,kte
     i = ismax(ite,tgrs(1,k),1)
     if(tgrs(i,k).gt.tgmxl) then
       tgmxl = tgrs(i,k)
       igmxl = i
       kgmxl = k
     endif
     i = ismin(ite,tgrs(1,k),1)
     if(tgrs(i,k).lt.tgmnl) then
       tgmnl = tgrs(i,k)
       igmnl = i
       kgmnl = k
     endif
   enddo
#endif
!
!  get surface pressure (cb : centi bar)
!
   do i = its,ite
#ifndef NIM
     psexp(i)      = exp(pgr(i))
#else
     psexp(i)      = pgr(i)
#endif
     psurf(i,lat)  = psexp(i)
     psfcpa(i)     = psexp(i) * 1.e3
     psmean(i,lat) = psmean(i,lat)+psexp(i)*dtf
   enddo
#ifndef NIM
#ifdef HYBRID
!
! ak5 : Pa -> cb in dyn_hybrid_setup.F90
!
   do i = its,ite
     temp = (psexp(i)/100.)**(akapa_)
     prsik(i,1)     = temp
     do k = kts,kte+1
       sihyb(i,k)   = ak5(levp1_-k+1)/psexp(i)+bk5(levp1_-k+1)
       cihyb(i,k)   = 1.-sihyb(i,k)
     enddo
#ifdef LINUX_PGI
!pgi$l novector
#endif
     do k = kts,kte
       slhyb(i,k)   = 0.5*(sihyb(i,k)+sihyb(i,k+1))
       slkhyb(i,k)  = slhyb(i,k)**(akapa_)
       delhyb(i,k)  = sihyb(i,k)-sihyb(i,k+1)
       clhyb(i,k)   = 1.-slhyb(i,k)
       prsi(i,k)    = sihyb(i,k)*psexp(i)
       prsik(i,k+1) = temp * (sihyb(i,k+1)**(akapa_))
       prsl(i,k)    = slhyb(i,k)*psexp(i)
       prslk(i,k)   = temp * slkhyb(i,k)
       prsl_pa(i,k) = prsl(i,k)*1.e3
       prsi_pa(i,k) = prsi(i,k)*1.e3
     enddo
     prsi(i,kte+1) = psexp(i) * sihyb(i,kte+1)
     prsi_pa(i,kte+1) = prsi(i,kte+1)*1.e3
   enddo
#else /* else of HYBRID */
!
   do i = its,ite
!
!  get 3-d pressure  (cb) : _pa (pascal)
!
     temp         = (psexp(i)/100.)**(akapa_)
     prsik(i,1)    = temp
     do k = kts,kte
       prsl(i,k)   = psexp(i)*sl(k)
       prsi(i,k)   = psexp(i)*si(k)
       prslk(i,k)  = temp   *slk(k)
       prsik(i,k+1)= temp   *sik(k+1)
       prsl_pa(i,k) = prsl(i,k)*1.e3
       prsi_pa(i,k) = prsi(i,k)*1.e3
     enddo
     prsi(i,kte+1) = psexp(i)*si(kte+1)
     prsi_pa(i,kte+1) = prsi(i,kte+1)*1.e3
   enddo
#endif /* HYBRID */
#else /* else of ~NIM */
!
   do i = its,ite
     temp         = exp(log(psexp(i)/100.)*(akapa_))
     do k = kts,kte
       prsl(i,k)    = prsl_nim(i,k)
       prsi(i,k)    = prsi_nim(i,k)
       prsl_pa(i,k) = prsl(i,k)*1.e3
       prsi_pa(i,k) = prsi(i,k)*1.e3
     enddo
     prsi(i,kte+1) = prsi_nim(i,kte+1)
     prsi_pa(i,kte+1) = prsi(i,kte+1)*1.e3
   enddo
#endif /* ~NIM */
!
#ifndef NONHYD_HYD
   call dyn_get_pressure(ite,ime,kte,akapa_,cp_,fv_,tgrs,qgrs,psexp,si,        &
                prsi,prsik,prsl,prslk,phii,phil,delprsi)
#else
   call dyn_get_pressure(ite,ime,kte,akapa_,cp_,fv_,tgrh,qgrs,psexp,si,        &
                prsi,prsik,prsl,prslk,phii,phil,delprsi)
#endif /* ~NONHYD_HYD */
!
#ifndef NIM
   do k = kts,kte
     do i = its,ite
       zl(i,k)=phil(i,k)/g_
       zi(i,k)=phii(i,k)/g_
       delprsi_pa(i,k) = delprsi(i,k)*1.e3
     enddo
   enddo
!
   do i = its,ite
     zi(i,kte+1)=phii(i,kte+1)/g_
   enddo
#else
   do k = kts,kte
     do i = its,ite
       zl(i,k) =nim_zm(i,k)
       zi(i,k) =nim_z(i,k)
       delprsi_pa(i,k) = delprsi(i,k)*1.e3
     enddo
   enddo
!
   do i = its,ite
       zi(i,kte+1) =nim_z(i,kte+1)
   enddo
#endif
!
   do k = kts,kte
     do i = its,ite
       delz(i,k) = zi(i,k+1)-zi(i,k)
     enddo
   enddo
#ifdef DBG
!
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(prsi,ite,ime,kme+1,kts,kte+1,                     &
                            'prsi in phys_main_solver')
     call print_maxmin_seven(delprsi,ite,ime,kme,kts,kte,                      &
                            'delprsi in phys_main_solver')
     call print_maxmin_seven(delprsi_pa,ite,ime,kme,kts,kte,                   &
                            'delprsi_pa in phys_main_solver')
     call print_maxmin_seven(prsl,ite,ime,kme,kts,kte,                         &
                            'prsl in phys_main_solver')
     call print_maxmin_seven(prsik,ite,ime,kme+1,kts,kte+1,                    &
                            'prsik in phys_main_solver')
     call print_maxmin_seven(prslk,ite,ime,kme,kts,kte,                        &
                            'prslk in phys_main_solver')
     call print_maxmin_seven(slk,kte,kme,1,1,1,                                &
                            'slk in phys_main_solver')
     call print_maxmin_seven(phii,ite,ime,kme+1,kts,kte+1,                     &
                            'phii in phys_main_solver')
     call print_maxmin_seven(phil,ite,ime,kme,kts,kte,                         &
                            'phil in phys_main_solver')
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
!  initialize dtdt with heating rate from rad_diurnal_cycle and 
!  get radsl for phys_lsm_osu1
!
#ifdef LFC
   if(1.ne.1) then
#endif
#ifdef SMP_NUDGING_RAD
     do i = 1,LONS2
!
! prescribe the surface albedo
!
#ifdef OSULSM1
       scmalb(i)=0.1957
#endif
#ifdef OSULSM2
       scmalb(i)=0.2506
#endif
#ifdef NOALSM1
       scmalb(i)=0.20
#endif
#ifdef VICLSM1
       scmalb(i)=0.20
#endif
!
! calculate SWUP
!
       scmswup(i)=scmalb(i)*scmswdn(i)
       sfcdlw(i,lat)=scmlwdn(i)/cnwatt
       sfcnsw(i,lat)=(scmswdn(i)-scmswup(i))/cnwatt
#ifdef VIC
       sfcdsw(i,lat)=scmswdn(i)/cnwatt
#endif
     enddo
!
#endif /* SMP_NUDGING_RAD end*/
     call rad_diurnal_cycle(ite,ime,kte,                                       &
#ifndef RMP
                solhr,slag,sinlab(1,lat),coslab(1,lat),                        &
                sdec,cdec,                                                     &
                xlon(1,lat),coszen(1,lat),                                     &
#else
                solhr,slag,sinlar(1,lat),coslar(1,lat),                        &
                sdec,cdec,                                                     &
                rlon(1,lat),coszer(1,lat),                                     &
#endif
                sfcdlw(1,lat),sfcnsw(1,lat),                                   &
#ifdef VIC
                sfcdsw(1,lat),dswsf1,                                          &
#endif
                tgrs(1,1),                                                     &
                tsea(1,lat),tsflw(1,lat),swh(1,1,lat),hlw(1,1,lat),            &
#ifdef NIM_DIAG
                xmu,                                                           &
#endif
                dlwsf1,ulwsf1,radsl,dtdt)
#ifdef OMLWRF
!
! net downward SW flux (W^m^2)
!
     do i = its,ite
       nswsf1(i)=radsl(i)*cnwatt-dlwsf1(i)
     enddo
!
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if(iope) then
#ifndef RMP
       call print_maxmin_seven(sinlab(1,lat),ite,ime,1,1,1,'sinlab')
#endif
       call print_maxmin_seven(tsflw(1,lat),ite,ime,1,1,1,'tsflw')
       call print_maxmin_seven(dlwsf1,ite,ime,1,1,1,'dlwsf1')
       call print_maxmin_seven(ulwsf1,ite,ime,1,1,1,'ulwsf1')
       call print_maxmin_seven(radsl,ite,ime,1,1,1,'radsl')
#ifdef OMLWRF
       call print_maxmin_seven(nswsf1,ite,ime,1,1,1,'nswsf1')
#endif
       call print_maxmin_seven(swh(1,1,lat),ite,ime,kte,kts,kte,'swh')
       call print_maxmin_seven(hlw(1,1,lat),ite,ime,kte,kts,kte,'hlw')
       call print_maxmin_seven(dtdt(1,1),ite,ime,kte,kts,kte,'dtdt')
     endif
#ifndef MP
     if(lat.eq.latg2_) then
       call print_maxmin_seven(tgrs(1,1),ite,ime,1,1,1,'rad_diurnal_cycle')
     endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, '::: NIM_DBGM : point(ndbg1,kdbg1)=',ndbg1,kdbg1
   print*, 'after rad_diurnal_cycle dlwsf1',ips,ipe,dlwsf1(ndbg1)
   print*, 'after rad_diurnal_cycle ulwsf1',ips,ipe,ulwsf1(ndbg1)
   print*, 'after rad_diurnal_cycle radsl ',ips,ipe,radsl(ndbg1)
   print*, 'after rad_diurnal_cycle swh   ',ips,ipe,swh(ndbg1,kdbg1,lat)
   print*, 'after rad_diurnal_cycle hlw   ',ips,ipe,hlw(ndbg1,kdbg1,lat)
   print*, 'after rad_diurnal_cycle dtdt  ',ips,ipe,dtdt(ndbg1,kdbg1)
#endif
#ifdef LFC
   else
     do i = its,ite
       radsl(i)  = sfcnsw(i,lat)+sfcdlw(i,lat)
       ulwsf1(i) = sbc_*tsea(i,lat)**4
     enddo
!
     do k = kts,kte
       do i = its,ite
         dtdt(i,k) = dtdt(i,k)+swh(i,k,lat)+hlw(i,k,lat)
       enddo
     enddo
   endif
#endif
!
   do i = its,ite
     dlwsfc(i,lat) = dlwsfc(i,lat)+dlwsf1(i)*dtf
#ifndef SMP_NUDGING_RAD
!
! change the step of updating LWUP by skh
!
     ulwsfc(i,lat) = ulwsfc(i,lat)+ulwsf1(i)*dtf
#endif
#ifdef SMP_NUDGING_RAD
     scmsdn(i,lat)=scmsdn(i,lat)+scmswdn(i)*dtf
     scmsup(i,lat)=scmsup(i,lat)+scmswup(i)*dtf
#endif
   enddo
!soojin_couple
#ifdef AOMG
   cpl_dlwsfc(:,lat) = cpl_dlwsfc(:,lat)+dlwsf1(:)*dtf
   cpl_ulwsfc(:,lat) = cpl_ulwsfc(:,lat)+ulwsf1(:)*dtf
#endif /*AOMG*/
!
   do k = kts,kte
     do i = its,ite
       hsw(i,k) = dtdt(i,k)-hlw(i,k,lat)
     enddo
   enddo
#ifdef DG3
   call diag_3d_archive(ite,ime,hsw, dtf,kdthsw,gda)
   call diag_3d_archive(ite,ime,hlw(1,1,lat), dtf,kdthlw,gda)
#endif
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
      write(6,*)'----- phys_main_solver -----  lat=',lat
   endif
   level=1
!
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(ugrs(1,1),ite,ime,kte,1,kte,                      &
                  'ugrs in phys_main_solver')
     call print_maxmin_seven(vgrs(1,1),ite,ime,kte,1,kte,                      &
                  'vgrs in phys_main_solver')
     call print_maxmin_seven(tgrs(1,1),ite,ime,kte,1,kte,                      &
                  'tgrs in phys_main_solver')
     call print_maxmin_seven(qgrs(1,1),ite,ime,ktotal,1,ktotal,                &
                  'qgrs in phys_main_solver')
     call print_maxmin_seven(hsw(1,1),ite,ime,kte,1,kte,                       &
                  'hsw in phys_main_solver')
     call print_maxmin_seven(hlw(1,1,lat),ite,ime,kte,1,kte,                   &
                  'hlw in phys_main_solver')
     call print_maxmin_seven(snoweq(1,lat),ite,ime,1,1,1,                      &
                  'snoweq in phys_main_solver')
     call print_maxmin_seven(tsea(1,lat),ite ,ime,1,1,1,                       &
                  'tsea in phys_main_solver')
     call print_maxmin_seven(smsoil(1,1),ite,ime,1,1,1,                        &
                  'smsoil in phys_main_solver')
     call print_maxmin_seven(stsoil(1,1),ite,ime,1,1,1,                        &
                  'stsoil in phys_main_solver')
     call print_maxmin_seven(canopy(1,lat),ite,ime,1,1,1,                      &
                  'canopy in phys_main_solver')
     call print_maxmin_seven(z0cm(1,lat),ite,ime,1,1,1,                        &
                  'z0cm in phys_main_solver')
     call print_maxmin_seven(z0cmt(1,lat),ite,ime,1,1,1,                       &
                  'z0cmt in phys_main_solver')
#ifdef OSULSM1
     call print_maxmin_seven(plantr(1,lat),ite,ime,1,1,1,                      &
                  'plantr in phys_main_solver')
#endif
     call print_maxmin_seven(tg3(1,lat),ite,ime,1,1,1,                         &
                  'tg3 in phys_main_solver')
     call print_maxmin_seven(slmsk(1,lat),ite,ime,1,1,1,                       &
                  'slmsk in phys_main_solver')
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
! Run merged LSM: lsm_exch_coeff -> lsm_driver -> lsm_diagnostic
!
#ifndef OSU
! assign precipitation and srflag from global to local array
!
   do i = its,ite
     oldprcp(i)   = prcp(i,lat)
     oldsrflag(i) = srflag(i,lat)
   enddo
#endif
#ifdef VIC
!
   call sfc_interp_date(idate(4),idate(2),idate(3),idate(1),                   &
                         thour,jyy,jmm,jdd,jhh,dummy)
#endif
#ifdef DBG
!
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     write(6,*)' before lsm_driver lat = ',lat
     write(6,*)ite,ime,                                                        &
        lsoil_,psexp(1),                                                       &
        ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),                               &
        snoweq(1,lat),tsea(1,lat),                                             &
        smsoil(1,1),stsoil(1,1),soiltyp(1),sigmaf(1),                          &
        vegtype(1),canopy(1,lat),                                              &
        dlwsf1(1),radsl(1),                                                    &
        deltim,z0cm(1,lat),z0cmt(1,lat),tg3(1,lat),                            &
        slk(1),slmsk(1,lat),                                                   &
        wind(1)
#ifdef VIC
     write(6,*)'for viclsm variables'
     write(6,*)nsoil_, jmm,                                                    &
        binf(1,lat),ds(1,lat),dsm(1,lat),ws(1,lat),cef(1,lat),                 &
        expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,wpwld,smrld,           &
        smxld,sicld,dpnld,sxnld,epnld,bbnld,apnld,btnld,gmnld,                 &
        flai(1,lat),vrtld,lstsnld,                                             &
        silz(1,lat),snwz(1,lat),csno(1,lat),rsno(1,lat),                       &
        tsf(1,lat),tpk(1,lat),sfw(1,lat),pkw(1,lat)
#endif
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
! compute exchange coeff (CD, CDQ, FFMM, FFHH)
!
   call lsm_exch_coeff(ite,ime,lsoil_,                                         &
               ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),                        &
               tsea(1,lat), z0cm(1,lat),z0cmt(1,lat),                          &
               cd,cdq,rb,rcl,                                                  &
               zl(1,1),prsi(1,1),prsl(1,1),prsik(1,1),prslk(1,1),              &
               slmsk(1,lat),inistp,lat,                                        &
#ifdef VSGD
               hpbl(1,lat),dxmeter,rainc,                                      &
#endif
               fm,fh,ustar,wind,                                               &
#ifdef VIC
               fm10,fh2,rho)
#else
               fm10,fh2)
#endif
!
! compute sfc energy/water balance
!
   call lsm_driver(ite,ime,lsoil_,                                             &
               tgrs(1,1),qgrs(1,1),                                            &
               snoweq(1,lat),tsea(1,lat),qss,                                  &
               smsoil,stsoil,evapc,soiltyp,sigmaf,                             &
#ifndef OSULSM1
               vegtype,canopy(1,lat),                                          &
#else
               canopy(1,lat),                                                  &
#endif
               dlwsf1,radsl,snowmt,snowev,                                     &
               deltim,z0cm(1,lat),tg3(1,lat),                                  &
               gflx,zsoil,cd,cdq,rcl,                                          &
               prsi(1,1),prsl(1,1),prsik(1,1),prslk(1,1),                      &
               zl(1,1),slmsk(1,lat),inistp,lat,                                &
               ims,ime,its,ite,                                                &
#ifndef NOAHYDRO
#ifdef OSU
               rhscnpy,rhsmc,ai,bi,cci,                                        &
#ifdef OSULSM1
               plantr(1,lat),cowave,                                           &
#endif
               drain,evap,hflx,ep1d,wind)
#else
               drain,evap,hflx,ep1d,wind,                                      &
#endif
#else
               drain,evap,hflx,ep1d,wind,evcn,                                 &
#endif
#ifdef NOALSM1
               snwdph(1,lat),slsoil,snoalb(1,lat),                             &
               slptyp,shdmin(1,lat),shdmax(1,lat),                             &
#endif
#ifdef VICLSM1
               nsoil_, jmm,dswsf1,                                             &
               binf(1,lat),ds(1,lat),dsm(1,lat),ws(1,lat),cef(1,lat),          &
               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,wpwld,smrld,    &
               smxld,sicld,dpnld,sxnld,epnld,bbnld,apnld,btnld,gmnld,          &
               flai(1,lat),vrtld,lstsnld,                                      &
               silz(1,lat),snwz(1,lat),csno(1,lat),rsno(1,lat),                &
               tsf(1,lat),tpk(1,lat),sfw(1,lat),pkw(1,lat),ght,                &
#endif
#ifdef VICLSM2
               nsoil_,msub_,jmm,dswsf1,                                        &
               binf(1,lat),ds(1,lat),dsm(1,lat),ws(1,lat),cef(1,lat),          &
               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,wpwld,          &
               smrld,smxld,dpnld,sxnld,epnld,bbnld,apnld,btnld,gmnld,          &
               nvegld,flaild,vfrld,vtypld,cnpld,snold,csnld,rsnld,             &
               tsfld,tpkld,sfwld,pkwld,lstsnld,vrtld,smcld,sicld,stcld,        &
               silz(1,lat),snwz(1,lat),                                        &
#endif
#ifndef OSU
               snowfl, runof,                                                  &
               oldprcp,oldsrflag)
#endif
!
!  update surface layer properties
!
   call lsm_diagnostic(ite,ime,lsoil_,                                         &
            ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),tsea(1,lat),               &
            radsl,rnet,                                                        &
            f10m(1,lat),u10m(1,lat),v10m(1,lat),t2m(1,lat),q2m(1,lat),         &
            rcl,prsi(1,1),prsik(1,1),prslk(1,1),slmsk(1,lat),inistp,lat,       &
            qss,evap,fm,fh,fm10,fh2)
#ifdef DBG
!
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     write(6,*)' after lsm_driver lat = ',lat
     write(6,*)ite,ime,lsoil_,                                                 &
        psexp(1),ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),                      &
        snoweq(1,lat),tsea(1,lat),qss(1),                                      &
        smsoil(1,1),stsoil(1,1),evapc(1),soiltyp(1),sigmaf(1),                 &
        vegtype(1),canopy(1,lat),                                              &
        dlwsf1(1),radsl(1),snowmt(1),snowev(1),                                &
        deltim,z0cm(1,lat),tg3(1,lat),                                         &
        gflx(1),f10m(1,lat),u10m(1,lat),v10m(1,lat),t2m(1,lat),                &
        q2m(1,lat),                                                            &
        (zsoil(1,k),k=1,lsoil_),cd(1),cdq(1),                                  &
        rcl,slk(1),slmsk(1,lat),lat,                                           &
        drain(1),evap(1),hflx(1),rnet(1),ep1d(1),fm(1),fh(1),ustar(1),         &
        wind(1)
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef OMLWRF
!
!  do ocean mixed layer interaction
!
   do i = its,ite
     if(slmsk(i,lat).eq.0.) then
#ifndef RMP
     coriolis = 2.*omega_*sinlab(i,lat)
#else
     coriolis = 2.*omega_*sinlar(i,lat)
#endif
     call ocean_mixed_layer_wrf(i,lat,tml(i,lat),t0ml(i,lat),                  &
         hml(i,lat),h0ml(i,lat),HUML(i,lat),HVML(i,lat),tsea(i,lat),hflx(i),   &
         evap(i),nswsf1(i),dlwsf1(i),ulwsf1(i),TMOML(i,lat),                   &
         ugrs(i,1),vgrs(i,1),ustar(i),                                         &
         coriolis,temp,temp,g_,deltim,oml_gamma,                               &
           ids,ide, jds,jde, kds,kde,                                          &
           ims,ime, jms,jme, kms,kme,                                          &
           its,ite, jts,jte, kts,kte                                          )
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       if(i.eq.10) print*,' ocean_mixed_layer_after ',i,lat,tml(i,lat),        &
         t0ml(i,lat),hml(i,lat),h0ml(i,lat),HUML(i,lat),HVML(i,lat),           &
         tsea(i,lat),hflx(i),evap(i),nswsf1(i),dlwsf1(i),ulwsf1(i),            &
         ugrs(i,1),vgrs(i,1),ustar(i),coriolis,temp,temp,deltim,oml_gamma                     
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
     call sst_skin_update(dlwsf1(i),ulwsf1(i),nswsf1(i),                       &
         hflx(i),evap(i),tsea(i,lat),ustar(i),temp,dtw1(i,lat),temp,           &
                deltim,temp,                                                   &
                ids, ide, jds, jde, kds, kde,                                  &
                ims, ime, jms, jme, kms, kme,                                  &
                its, ite, jts, jte, kts, kte                                  )
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       if(i.eq.10) print*,' sst_skin_update ',dlwsf1(i),ulwsf1(i),nswsf1(i),   &
         hflx(i),evap(i),tsea(i,lat),ustar(i),temp,dtw1(i,lat),temp,deltim                                                    
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
     endif
   enddo
#endif
!
!  Run 2-step LSM (delete phys_lsm_osu1 or phys_lsm_osu2)
!
#ifdef RASV2
   do i = its,ite
     dsfc(i)=0.
   enddo
!
#endif
#ifdef MRGLSM
   do i = its,ite
     snowfall(i,lat) = snowfall(i,lat) + snowfl(i)
   enddo
!
#endif
   do i = its,ite
#ifdef HYDRO
     snowevap(i,lat) = snowevap(i,lat) + snowev(i)*dtf
#else
#ifdef NOAHYDRO
     if (slmsk(i,lat).eq.1)                                                    &
       snowevap(i,lat) = snowevap(i,lat) + snowev(i)*deltim*(hvap_+hfus_)
     if (slmsk(i,lat).eq.2)                                                    &
       snowevap(i,lat) = snowevap(i,lat) + snowev(i)*deltim*(hvap_+hfus_)/hvap_
#else
     snowevap(i,lat) = snowevap(i,lat) + snowev(i)*deltim
#endif
#endif /* HYDRO end */
     gflux(i,lat) = gflux(i,lat)+gflx(i)*dtf
     tmpmax(i,lat) = max(tmpmax(i,lat),t2m(i,lat))
     tmpmin(i,lat) = min(tmpmin(i,lat),t2m(i,lat))
     alhtfl(i,lat) = alhtfl(i,lat) + evap(i) * dtf
#ifdef HYDRO
     evcnp(i,lat) = evcnp(i,lat) + evcn(i) * dtf
#else
#ifdef NOAHYDRO
     if (slmsk(i,lat).eq.1)                                                    &
       evcnp(i,lat) = evcnp(i,lat) + evcn(i) * deltim
#endif
     ep(i,lat) = ep(i,lat) + ep1d(i) * dtf
#endif
#ifdef VICLSM1
     gheat(i,lat) =gheat(i,lat)+ght(i)*dtf
#endif
#ifndef OSULSM1
     uustar(i,lat) = ustar(i)
     ffmm(i,lat) = fm(i)
     ffhh(i,lat) = fh(i)
#endif
   enddo
#ifdef SMP_NUDGING_RAD
!
! update the LWUP using updated Ts
!
   do i = its,ite
     ulwsf1(i)=sbc_*tsea(i,lat)**4
     ulwsfc(i,lat)=ulwsfc(i,lat)+ulwsf1(i)*dtf
   enddo
!
! end of update LWUP
#endif
!
!     compute coefficient of evaporation in evapc
!
   do i = its,ite
     tstar(i) = tsea(i,lat)
     if (evapc(i).gt.1.e0) evapc(i) = 1.0e0
     evapol(i) = evap(i)/hvap_
#ifdef SMP_NUDGING_FLX
     hflx(i) = scmsh(1)
     evapol(i) = scmlh(1)/hvap_
#endif
   enddo
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(tstar,ite,ime,1,1,1,                              &
                            'tstar in phys_main_solver')
     call print_maxmin_seven(qss,ite,ime,1,1,1,                                &
                            'qss in phys_main_solver')
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
!     over snow cover or ice or sea, coef of evap =1.0. e 0
!
   do i = its,ite
     if ((snoweq(i,lat).gt.0.e0) .or. (slmsk(i,lat).ne.1.e0))                  &
       evapc(i)=1.e0
   enddo
!
!  do vertical diffusion
!
#ifdef CPS_KSAS_DIUR
   do i = its,ite
     hpbl_hold(i) = hpbl(i,lat)     ! for sas
   enddo
#endif

#ifdef DG3
   call diag_3d_archive(ite,ime,dtdt,-dtf,kdtvrdf,gda)
#endif
#ifdef VIC
   do i = 1,ime
     if (tsea(i,lat).lt.273.15) then
       ls = (677.-0.07*(tsea(i,lat)-273.15))*4.1868*1000.0
       evp(i) = evap(i)/(rho(i)*ls)
     else
       le = (2.501-0.002361*(tsea(i,lat)-273.15))*1.e6
       evp(i) = evap(i)/(rho(i)*le)
     endif
   enddo
!
#endif
#ifdef NIM_DIAG
   do k = kts,kte
     do i = its,ite
       ttend_nim(i,k) = dtdt(i,k)       ! pbl tmp tendency
     enddo
   enddo
!
#endif
#ifdef YSUPBL
#ifdef YSUTKE
#ifdef RMP
     do k = kts,kte
       do i = its,ite
         qgrs_tke(i,k) = qgrs(i,k+kbltke_-1)
       enddo
     enddo
!
#endif
#endif
   call phys_pbl_ysu(lat,ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),              &
                  prsl_pa(1,1),prsi_pa(1,1),prslk(1,1),                        &
                  dudt,dvdt,dtdt,dqdt,nvdiff,nwmass_,ifvmix,                   &
                  cp_,g_,akapa_,rd_,rdog_,fv_,rdorv_,karman_,hvap_,rv_,        & 
                  delz(1,1),psfcpa(1),                                         &
                  z0cm(1,lat),ustar,hpbl(1,lat),fm,fh,                         &
                  slmsk(1,lat),hflx,evapol,wind,rb,                            &
                  dusfc1,dvsfc1,dtsfc1,dqsfc1,                                 &
                  deltim,rcl,kpbl,                                             &
                  gamt,gamq,                                                   &
                  temp_ik,                                                     &
                  wstar,delta,                                                 &
#ifdef YSUTKE     
#ifndef RMP
                  qgrs(1,kbltke_),el_pbl,corf,                                 &
#else
                  qgrs(1,kbltke_),el_pbl,rcorf,                                &
#endif
#endif
                  u10m(1,lat),v10m(1,lat),                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte)
#endif
#ifdef MRFPBL
   call phys_pbl_mrf(lat,ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),              &
                  prsl_pa(1,1),prsi_pa(1,1),prslk(1,1),                        &
                  dudt,dvdt,dtdt,dqdt,nvdiff,                                  &
                  cp_,g_,akapa_,rd_,rdog_,fv_,rdorv_,karman_,hvap_,rv_,        & 
                  delz(1,1),psfcpa(1),                                         &
                  z0cm(1,lat),ustar,hpbl(1,lat),fm,fh,                         &
                  slmsk(1,lat),hflx,evapol,wind,rb,                            &
                  dusfc1,dvsfc1,dtsfc1,dqsfc1,                                 &
                  deltim,rcl,kpbl,                                             &
                  gamt,gamq,                                                   &
                  temp_ik,                                                     &
                  u10m(1,lat),v10m(1,lat),                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte)
#endif
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(dudt(1,1),ite,ime,kte,1,kte,                      &
                            'dudt after phys_pbl')
     call print_maxmin_seven(dvdt(1,1),ite,ime,kte,1,kte,                      &
                            'dvdt after phys_pbl')
     call print_maxmin_seven(dtdt(1,1),ite,ime,kte,1,kte,                      &
                            'dtdt after phys_pbl')
     call print_maxmin_seven(dqdt(1,1),ite,ime,ktotal,1,ktotal,                &
                            'dqdt after phys_pbl')
   endif
!
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, 'after pbl ugrs,vgrs',ips,ipe,ugrs(ndbg1,kdbg1),vgrs(ndbg1,kdbg1)
   print*, 'after pbl tgrs,qgrs',ips,ipe,tgrs(ndbg1,kdbg1),qgrs(ndbg1,kdbg1)
   print*, 'after pbl gu0,dqdt ',ips,ipe,gu0(ndbg1,kdbg1),dudt(ndbg1,kdbg1)
   print*, 'after pbl gv0,dqdt ',ips,ipe,gv0(ndbg1,kdbg1),dvdt(ndbg1,kdbg1)
   print*, 'after pbl gt0,dtdt ',ips,ipe,gt0(ndbg1,kdbg1),dtdt(ndbg1,kdbg1)
   print*, 'after pbl gq0,dqdt ',ips,ipe,gq0(ndbg1,kdbg1),dqdt(ndbg1,kdbg1)
#endif


#ifdef NIM_DIAG
   do k = kts,kte
     do i = its,ite
       ttend_nim(i,k) = dtdt(i,k) - ttend_nim(i,k)
     enddo
   enddo
!
#endif
#ifdef SMP
#ifdef SMP_NUDGING_RAD
   do i = its,ite
     print*,'in phys_main_solver replace the heat fluxes'
     dtsfc1(i)=hflx(i)
     dqsfc1(i)=evap(i)
   enddo
!
#endif
#endif
!soojin_couple
#ifdef AOMG
   cpl_dusfc(:,lat)=
   cpl_dvsfc(:,lat)=
   cpl_dtsfc(:,lat)=
   cpl_dqsfc(:,lat)=
#endif /*AOMG*/
!
!  do gravity wave drag induced by orography
!
#ifndef NOGWD
   do i = its,ite
     var(i) = hprime(i,lat,1)
   enddo
!
#ifndef KAGWD
   call phys_gwd_alpert(dudt,dvdt,                                             &
              ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),                         &
              delprsi,prsi,prsl,prslk,zl,rcl,                                  &
              deltim,lat,kdt,var,dusfcg,dvsfcg,                                &
              g_,cp_,rd_,rv_,fv_,                                              &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte)
#else
   do i = its,ite
     oc(i) = hprime(i,lat,2)
   enddo
!
   do k = 1, 4
     do i = its,ite
       oa4(i,k) = hprime(i,lat,k+2)
       ol4(i,k) = hprime(i,lat,k+6)
     enddo
   enddo
#ifdef FLOW_BLOCKING
   do i = its,ite
     omax(i) = hprime(i,lat,11)
   enddo
#endif
!
   call phys_gwd_kimarakawa(dudt,dvdt,                                         &
             dtaux(1,1),dtauy(1,1),                                            &
             ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),                          &
             prsi_pa(1,1),prsl_pa(1,1),prslk(1,1),zl(1,1),rcl,kpblmax,         &
#ifdef FLOW_BLOCKING
             omax(1),                                                          &
#endif
             var(1),oc(1),oa4(1,1),ol4(1,1),dusfcg,dvsfcg,                     &
             g_,cp_,rd_,rv_,fv_,pi_,dxmeter,deltim,kpbl,kdt,lat,               &
#if defined(RMP) && defined(FLOW_BLOCKING)
             rdelx,rdely,                                                      &
#endif
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte)
#endif
#ifdef DBG
!
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(var,ite,ime,1,1,1,                                &
                            'hprime after phys_gwdo')
     call print_maxmin_seven(dudt,ite,ime,kte,1,kte,                           &
                            'dudt after phys_gwdo')
     call print_maxmin_seven(dvdt,ite,ime,kte,1,kte,                           &
                            'dvdt after phys_gwdo')
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, 'after gwdo dudt',ips,ipe,dudt(ndbg1,kdbg1)
   print*, 'after gwdo dvdt',ips,ipe,dvdt(ndbg1,kdbg1)
#endif
!
   do i = its,ite
     dugwd(i,lat) = dugwd(i,lat)+dusfcg(i)*dtf
     dvgwd(i,lat) = dvgwd(i,lat)+dvsfcg(i)*dtf
   enddo
#endif /* ~NOGWD */
!
#ifdef NONHYD_HYD
!
! ustr = u / m    vstr = v /m
! dwdt=m**2(d ustr/dt*dzdx + d vstr/dt*dzdy)
! dwdt=m(dudt*dgzdx+dvdt*dgzdy)/g
!
   do i = its,ite
     dwdt(i,1)=dudt(i,1)*gzdx(i,lat)+dvdt(i,1)*gzdy(i,lat)
     dwdt(i,1)=dwdt(i,1)*xm(i,lat)/g_
   enddo
#endif /* NONHYD_HYD */
!
#ifdef DG3
   call diag_3d_archive(ite,ime,dudt, dtf,kduvrdf,gda)
   call diag_3d_archive(ite,ime,dvdt, dtf,kdvvrdf,gda)
   call diag_3d_archive(ite,ime,dtdt, dtf,kdtvrdf,gda)
   call diag_3d_archive(ite,ime,dqdt, dtf,kdqvrdf,gda)
#endif
#ifndef SMP
   if(inistp.le.1) then
!
!   get vertical motion (cb/sec) in vvel for phys_cps_sas
!
#ifdef NONHYD
     do k = kts,kte
       do i = its,ite
         ww(i,k) = 0.5*(wgrs(i,k)+wgrs(i,k+1))
         vvel(i,k) = -ww(i,k)*delprsi(i,k)/delz(i,k)
#ifdef KSAS
         vvel_pa(i,k) = vvel(i,k) * 1.e3
#endif
       enddo
     enddo
!
#else /* else of NONHYD */
     call dyn_get_omega(ite,ime,kte,                                           &
               pphigr(1),plamgr(1),w2,                                         &
               ugrs(1,1),vgrs(1,1),                                            &
#ifndef RMP
               xgrs(1,1),delprsi,rbs2(lat),vvel,                               &
#else
               xgrs(1,1),delprsi,rbs2,vvel,                                    &
#endif
               psexp,prsl)
!
     do k = kts,kte
       do i = its,ite
         ww(i,k) = -vvel(i,k)*delz(i,k)/delprsi(i,k)
#ifdef KSAS
         vvel_pa(i,k) = vvel(i,k) * 1.e3
#endif
       enddo
     enddo
!
#endif /* endif NONHYD */
   endif
#else /* else of ndef SMP */
!
!  OMEGA : vvel -> xgrs in SMP 
!
   do k = kts,kte
     do i = its,ite
       ww(i,k) = -xgrs(i,k)*delz(i,k)/delprsi(i,k)
     enddo
   enddo
#endif /* endif SMP */
!
!   w2 is used to store wind speed
!
   do k = kts,kte
      do i = its,ite
        w2(i,k) = rcs*sqrt(ugrs(i,k)**2+vgrs(i,k)**2)
      enddo
   enddo
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifndef MP
   if(lat.eq.latd.and.iope) then
      call print_maxmin_seven(qgrs(1,1),ite,ime,ktotal,1,ktotal,               &
                             'gq0 before update')
   endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
!  update u v t q due to vertical diffusion and gwdo
!
   do k = kts,kte
     do i = its,ite
       gt0(i,k) = tgrs(i,k)+dtdt(i,k)*dt2
       gu0(i,k) = ugrs(i,k)+dudt(i,k)*dt2
       gv0(i,k) = vgrs(i,k)+dvdt(i,k)*dt2
     enddo
   enddo
!
#ifdef NONHYD_HYD
   do k = kts,kte+1
     do i = its,ite
       gw0(i,k)=wgrs(i,k)+dwdt(i,k)*dt2
     enddo
   enddo
!
#endif /* NONHYD_HYD */
   do k = 1,levh_
     do i = its,ite
       gq0(i,k) = qgrs(i,k)+dqdt(i,k)*dt2
     enddo
   enddo
#ifdef YSUTKE
#ifdef RMP
!
   do k = kts,kte
     do i = its,ite
       qgrs(i,k+kbltke_-1) = qgrs_tke(i,k)
     enddo
   enddo
#endif
#endif
!
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifndef MP
   if(lat.eq.latd.and.iope) then
     call print_maxmin_seven(gq0(1,1),ite,ime,ktotal,1,ktotal,                 &
                             'gq0 after update')
   endif
!
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
!  do rad_ozone physcis
!
#ifndef AQUA_PLANET
   if(ngases_ .ge. 1) then
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if(lat.eq.latd.and.iope) then
       call print_maxmin_seven(gq0(1,kgases_),ite,ime,kte,1,kte,               &
                              'o3 before rad_ozone_physics')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
     call rad_ozone_physics(ite,ime,kte,deltim,                                &
                 gq0(1,kgases_),gq0(1,kgases_),prsi,prsl,                      &
#ifndef RMP
                 xlat(1,lat),                                                  &
#else
                 rlat(1,lat),                                                  &
#endif
#ifndef O3CHEM
                 pgr,lat)
#else
                 tgrs(1,1),qgrs(1,1),psexp,coszen(1,lat),lat,inistp)
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if(lat.eq.latd.and.iope) then
       call print_maxmin_seven(gq0(1,kgases_),ite,ime,kte,1,kte,               &
                                'o3 after rad_ozone_physics')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
   endif
#endif /* ~AQUA_PLANET */
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if (iope) then
     call print_maxmin_seven(gt0,ite,ime,kte,1,kte,'gt0 update')
     call print_maxmin_seven(gq0,ite,ime,ktotal,1,ktotal,'gq0 update')
     call print_maxmin_seven(gu0,ite,ime,kte,1,kte,'gu0 update')
     call print_maxmin_seven(gv0,ite,ime,kte,1,kte,'gv0 update')
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
#ifndef NONHYD
! update zl, zi, delz in case of hydrostatic dynamics
!
   call dyn_get_height(ite,ime,kte,cp_,fv_,gt0(1,1),gq0(1,1),                  &
                       prsik,prslk,phii,phil)
!
   do k = kts,kte
     do i = its,ite
       zl(i,k)=phil(i,k)/g_
       zi(i,k)=phii(i,k)/g_
     enddo
   enddo
!
   do i = its,ite
     zi(i,kte+1)=phii(i,kte+1)/g_
   enddo
!
   do k = kts,kte
     do i = its,ite
       delz(i,k) = zi(i,k+1)-zi(i,k)
     enddo
   enddo
#endif /* ~NONHYD */
!
#ifdef CHEM
   do k = kts,kte
     do i = its,ite
       c_moistq(k,i,lat) = gq0(i,k)
     enddo
   enddo
#endif
!
!  do cumulus parameterization physics
!
   if(inistp.le.1) then
#ifdef DG3
     call diag_3d_archive(ite,ime,gt0(1,1),-frain,kdtconv,gda)
     call diag_3d_archive(ite,ime,gq0(1,1),-frain,kdqconv,gda)
#endif
!
!  for cps ensemble (sas ras kf2 ccmcnv or kuo)
!
#ifdef CPS_ENS
     gq0i = gq0
     gt0i = gt0
     gu0i = gu0
     gv0i = gv0
#endif
!
#ifdef SAS
#ifdef SMP
#define JCAPS jcapset
#define VVELS xgrs
#else
#define JCAPS jcapr
#define VVELS vvel
#endif
     call phys_cps_sas(deltim,dxmeter,delprsi,prsl,prsi,prslk,prsik,zl,zi,     &
#ifdef DET_HYDRO  
               gq0(1,kcloud_),                                                 &
#endif 
#ifdef SAS_CCN
               gq0(1,kwsize_),                                                 &
#endif
               gq0(1,1),gt0(1,1),gu0(1,1),gv0(1,1),rcs,slmsk(1,lat),VVELS,     &
               cld1d,rain1,                                                    &
#ifdef MUL_CLDTOP
               xkt2,NCLDTOP,                                                   &
#endif
#ifdef CLM_CWF
               cgs2,                                                           &
#endif
#ifdef SAS_DIAG
               dcu,dcv,dct,dcq,dch,fcu,fcd,                                    &
               deltb,delqb,delhb,cbmf,                                         &
#else
#ifdef GWDC
               dct,fcu,fcd,                                                    &
#endif
#endif
               hpbl(1,lat),kpbl,gamt,gamq,                                     &
#ifdef EXPLICIT_CLOUDINESS
               qci,qrs,                                                        &
#endif
               JCAPS,ncloud_,lat,kbot,ktop,icps,                               &
               ids,ide, jds,jde, kds,kde,                                      &
               ims,ime, jms,jme, kms,kme,                                      &
               its,ite, jts,jte, kts,kte)
!
#ifdef EXPLICIT_CLOUDINESS
     do i = its,ite
       do k = kts,kte
         qcicps(i,k,lat) = qci(i,k)
         qrscps(i,k,lat) = qrs(i,k)
       enddo
     enddo
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       call print_maxmin_seven(gt0,ite,ime,kte,1,kte,                          &
                              'gt0 after phys_cps')
       call print_maxmin_seven(gq0,ite,ime,ktotal,1,ktotal,                    &
                              'gq0 after phys_cps')
       call print_maxmin_seven(gu0,ite,ime,kte,1,kte,                          &
                              'gu0 after phys_cps')
       call print_maxmin_seven(gv0,ite,ime,kte,1,kte,                          &
                              'gv0 after phys_cps')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, 'after cps gt0',ips,ipe,gt0(ndbg1,kdbg1)
   print*, 'after cps gq0',ips,ipe,gq0(ndbg1,kdbg1)
   print*, 'after cps gu0',ips,ipe,gu0(ndbg1,kdbg1)
   print*, 'after cps gv0',ips,ipe,gv0(ndbg1,kdbg1)
#endif
#ifdef CPS_ENS
!
! for cps ens (sas)
!
     gqsas = gq0
     gtsas = gt0
     gusas = gu0
     gvsas = gv0
     cld1dsas = cld1d
     rain1sas = rain1
     icpssas  = icps
     kbotsas  = kbot
     ktopsas  = ktop
     gq0 = gq0i
     gt0 = gt0i
     gu0 = gu0i
     gv0 = gv0i
!
#endif
#endif /* SAS end */
!
!============================= yng ADD KSAS ===========================
#ifdef KSAS
     call phys_cps_ksas(deltim,dxmeter,delprsi_pa,prsl_pa,prsi_pa,zl,zi,     &
#ifdef DET_HYDRO  
               gq0(1,kcloud_),                                                 &
#endif 
#ifdef SAS_CCN
               gq0(1,kwsize_),                                                 &
#endif
               gq0(1,1),gt0(1,1),gu0(1,1),gv0(1,1),rcs,slmsk(1,lat),vvel_pa,      &
               cld1d,rain1,                                                    &
               dct,fcu,fcd,                                                    &
               hpbl(1,lat), &
#ifdef CPS_KSAS_DIUR
               hpbl_hold,                                         &
#endif
#ifdef EXPLICIT_CLOUDINESS
               qci,qrs,                                                        &
#endif
               ncloud_,lat,kbot,ktop,icps,                               &
               fv_,cp_,cv_,cvap_,cliq_,hvap_,                          &
               g_,pi_,rd_,rv_,t0c_,qmin_,  &
               ids,ide, jds,jde, kds,kde,                                      &
               ims,ime, jms,jme, kms,kme,                                      &
               its,ite, jts,jte, kts,kte)
!
#ifdef EXPLICIT_CLOUDINESS
     do i = its,ite
       do k = kts,kte
         qcicps(i,k,lat) = qci(i,k)
         qrscps(i,k,lat) = qrs(i,k)
       enddo
     enddo
#endif
#endif /* KSAS end */
!============================= yng ADD KSAS ===========================
!
#ifdef RAS
     dtras = deltim / frain
     il1r = 1 + 8 * ime
     il2r = il1r + 10 * ime * lm
     il3r = il2r +  ime * (5*lm+16)
     call phys_cps_ras(ite, ime,                                               &
               kte, lm, nstrp, dtras,                                          &
               prsi, prsl,                                                     &
               krmin, krmax, ncrnd, afac, rannum, ufac,                        &
               rgas, cp, grav, alhl,                                           &
               gt0(1,1), gq0(1,1),                                             &
               gu0(1,1), gv0(1,1),                                             &
               rain1, kbot, ktop, icps, lat, cd,                               &
               sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns)
!
     do i = its,ite
       cld1d(i) = 0.
     enddo
#ifdef CPS_ENS
!
! for cps ens (ras)
!
     gqras = gq0
     gtras = gt0
     guras = gu0
     gvras = gv0
     cld1dras = cld1d
     rain1ras = rain1
     icpsras  = icps
     kbotras  = kbot
     ktopras  = ktop
     gq0 = gq0i
     gt0 = gt0i
     gu0 = gu0i
     gv0 = gv0i
#endif
#endif /* RAS end */
!
#ifdef RASV2
     dtras = deltim / frain
     do i = its,ite
       do k = kts,kte
         gt00(k) = gt0(i,k)
         gq00(k) = gq0(i,k)
         gu00(k) = gu0(i,k)
         gv00(k) = gv0(i,k)
       enddo
!
       call phys_cps_ras2(kte, dtras                                           &
                  ,krmin, krmax, kfmax, ncrnd, rannum, pdd                     &
                  ,mct, kctop                                                  &
                  ,rgas, cp, grav, alhl, hfus_                                 &
                  ,gt00, gq00, gu00, gv00, psexp(i), kpbl(i)                   &
                  ,rain1(i), kbot(i), ktop(i), icps(i), lat, cd(i)             &
                  ,sig, prj, sgb, rasal, clw0, clt, dsfc(i))
!
       do k = kts,kte
         gt0(i,k) = gt00(k)
         gq0(i,k) = gq00(k)
         gu0(i,k) = gu00(k)
         gv0(i,k) = gv00(k)
         clw(i,k) = clw0(kte-k+1)
       enddo
     enddo
!
     do i = its,ite
       cld1d(i) = 0.
     enddo
#ifdef CPS_ENS
!
! for cps ens (rasv2)
!
     gqras = gq0
     gtras = gt0
     guras = gu0
     gvras = gv0
     cld1dras = cld1d
     rain1ras = rain1
     icpsras  = icps
     kbotras  = kbot
     ktopras  = ktop
     gq0 = gq0i
     gt0 = gt0i
     gu0 = gu0i
     gv0 = gv0i
#endif
#endif /* RASV2 end */
!
#ifdef CCMCNV
     call phys_cps_ccm3(ite,ime,kte,                                           &
#ifdef SMP
               jcapset,deltim,kdt,prsi,prsl,                                   &
               gq0(1,1),gt0(1,1),cld1d,                                        &
               rain1,kbot,ktop,icps,w2,lat,slmsk(1,lat),xgrs,                  &
#else
               jcapr,deltim,kdt,prsi,prsl,                                     &
               gq0(1,1),gt0(1,1),cld1d,                                        &
               rain1,kbot,ktop,icps,w2,lat,slmsk(1,lat),vvel,                  &
#endif
               hpbl(1,lat),gamt,gamq,dxmeter)
#ifdef CPS_ENS
!
! for cps ens (ccm3)
!
     gqccm = gq0
     gtccm = gt0
     guccm = gu0
     gvccm = gv0
     cld1dccm = cld1d
     rain1ccm = rain1
     icpsccm  = icps
     kbotccm  = kbot
     ktopccm  = ktop
     gq0 = gq0i
     gt0 = gt0i
     gu0 = gu0i
     gv0 = gv0i
!
#endif
#endif /* CCMCNV end */
!
#ifdef KF2
     call phys_cps_kf2(ite,ime,kte,gu0(1,1),gv0(1,1),                          &
              gt0(1,1),gq0(1,1),delprsi,prsi,prsl,                             &
              vvel,deltim,dxmeter,rain1,kbot,ktop,icps,lat)
     !
     do i = its,ite
       cld1d(i) = 0.
     enddo
#ifdef CPS_ENS
!
! for cps ens (kf2)
!
     gqkf2 = gq0
     gtkf2 = gt0
     gukf2 = gu0
     gvkf2 = gv0
     cld1dkf2 = cld1d
     rain1kf2 = rain1
     icpskf2  = icps
     kbotkf2  = kbot
     ktopkf2  = ktop
     gq0 = gq0i
     gt0 = gt0i
     gu0 = gu0i
     gv0 = gv0i
#endif
#endif /* KF2 end */
!
#ifdef KUO
     call phys_cps_kuo(deltim,delprsi,prsl,prsi,prslk,prsik,                   &
                       qgrs(1,1),gq0(1,1),gt0(1,1),rain1,                      &
                       jcap_,kbot,ktop,icps,rd_,rv_,cp_,t0c_,g_,hvap_,         &
                             ids,ide, jds,jde, kds,kde,                        &
                             ims,ime, jms,jme, kms,kme,                        &
                             its,ite, jts,jte, kts,kte)
#ifdef CPS_ENS
!
! for cps ens (kuo)
!
     gqkuo = gq0
     gtkuo = gt0
     gukuo = gu0
     gvkuo = gv0
     cld1dkuo = cld1d
     rain1kuo = rain1
     icpskuo  = icps
     kbotkuo  = kbot
     ktopkuo  = ktop
#endif
#endif /* KUO end */
!
#ifdef CPS_ENS
!
! for CPS ensemble
!
#ifdef KUO
     gq0 = (gqsas+gqras+gqccm+gqkf2+gqkuo)*0.2
     gt0 = (gtsas+gtras+gtccm+gtkf2+gtkuo)*0.2
     gu0 = (gusas+guras+guccm+gukf2+gukuo)*0.2
     gv0 = (gvsas+gvras+gvccm+gvkf2+gvkuo)*0.2
     cld1d = (cld1dsas+cld1dras+cld1dccm+cld1dkf2+cld1dkuo)*0.2
     rain1 = (rain1sas+rain1ras+rain1ccm+rain1kf2+rain1kuo)*0.2
     icps  = nint((icpssas+icpsras+icpsccm+icpskf2+icpskuo)*0.2)
     kbot  = nint((kbotsas+kbotras+kbotccm+kbotkf2+kbotkuo)*0.2)
     ktop  = nint((ktopsas+ktopras+ktopccm+ktopkf2+ktopkuo)*0.2)
#else
     gq0 = (gqsas+gqras+gqccm+gqkf2)*0.25
     gt0 = (gtsas+gtras+gtccm+gtkf2)*0.25
     gu0 = (gusas+guras+guccm+gukf2)*0.25
     gv0 = (gvsas+gvras+gvccm+gvkf2)*0.25
     cld1d = (cld1dsas+cld1dras+cld1dccm+cld1dkf2)*0.25
     rain1 = (rain1sas+rain1ras+rain1ccm+rain1kf2)*0.25
     icps  = nint((icpssas+icpsras+icpsccm+icpskf2)*0.25)
     kbot  = nint((kbotsas+kbotras+kbotccm+kbotkf2)*0.25)
     ktop  = nint((ktopsas+ktopras+ktopccm+ktopkf2)*0.25)
#endif
#endif
!
!  do gwdc - gravity wave drag induced by convection
!
#ifdef GWDC
     call phys_gwd_chunbaik(ite,ime,kte,                                       &
                  gu0(1,1),gv0(1,1),gt0(1,1),                                  &
                  prsi,prsl,                                                   &
#if defined SAS || defined KSAS
                  dct, fcu, fcd,                                               &
#else
#ifdef CCMCNV
                  hrate, netflux,                                              &
#endif
#endif
                  rcl,deltim,ktop,kbot,                                        &
                  tauctx,taucty,                                               &
#ifdef DG3
                  dudt,dvdt,                                                   &
#endif
#ifdef RMP
                  rdelx,rdely,                                                 &
#endif
                  lat)
!
     do i = its,ite
       dugwdc(i,lat) = dugwdc(i,lat)+tauctx(i)*dtf
       dvgwdc(i,lat) = dvgwdc(i,lat)+taucty(i)*dtf
     enddo
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       call print_maxmin_seven(gu0,ite,ime,kte,kts,kte,                        &
                              'gu0 aft phys_gwdc')
       call print_maxmin_seven(gv0,ite,ime,kte,kts,kte,                        &
                              'gv0 aft phys_gwdc')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, 'after gwdc gu0',ips,ipe,gu0(ndbg1,kdbg1)
   print*, 'after gwdc gv0',ips,ipe,gv0(ndbg1,kdbg1)
#endif
#endif
!
     do  i = its,ite
       cldwrk(i,lat) = cldwrk(i,lat) + cld1d(i) * dtf
     enddo
#ifdef DG3
     call diag_3d_archive(ite,ime,gt0(1,1), frain,kdtconv,gda)
     call diag_3d_archive(ite,ime,gq0(1,1), frain,kdqconv,gda)
#endif
#ifdef DG3
     call diag_3d_archive(ite,ime,gt0(1,1),-frain,kdtshal,gda)
     call diag_3d_archive(ite,ime,gq0(1,1),-frain,kdqshal,gda)
#endif
!
!  do shallow convection
!
#ifdef CCMSCV
     call phys_scv_ccm2(ite,ime,kte,                                           &
#ifdef SMP
               jcapset,deltim,kdt,prsi,                                        &
               gq0(1,1),gt0(1,1),                                              &
               rain2,kbot,ktop,icps,w2,lat,slmsk(1,lat),xgrs,                  &
#else
               jcap_,deltim,kdt,prsi,                                          &
               gq0(1,1),gt0(1,1),                                              &
               rain2,kbot,ktop,icps,w2,lat,slmsk(1,lat),vvel,                  &
#endif
               hpbl(1,lat),gamt,gamq,dxmeter)
!
     do i = its,ite
       rainc(i)      = frain*(rain1(i)+rain2(i))
       raincps(i,lat) = raincps(i,lat)+rainc(i)
     enddo
#endif
#ifdef MRFSCV
     call phys_scv_mrf(levshc,deltim,                                          &
                    delprsi,prsi,prsik,prsl,prslk,zl,icps,gq0(1,1),gt0(1,1),   &
                    cp_,g_,hvap_,rd_,rv_,                                      &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte)
!
     do i = its,ite
       rainc(i)      = frain*rain1(i)
#ifdef SMP_NUDGING_RAD
       rainc(i)=prec(i)
#endif
       raincps(i,lat) = raincps(i,lat)+rainc(i)
     enddo
#endif
#ifdef GRIMSCV
     call phys_scv_grims(gq0(1,1),gt0(1,1),prsi_pa(1,1),prsik(1,1),            &
                    delprsi_pa(1,1),prsl_pa(1,1),prslk(1,1),zl(1,1),           & 
                    wstar(1),hpbl(1,lat),delta(1),                             &
                    deltim,cp_,g_,hvap_,rd_,rv_,                               &
                    icps,kpbl,levshc,                                          &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte)
!
     do i = its,ite
       rainc(i)      = frain*rain1(i)
#ifdef SMP_NUDGING_RAD
       rainc(i)=prec(i)
#endif
       raincps(i,lat) = raincps(i,lat)+rainc(i)
     enddo
#endif
#ifdef HANPANSCV
!
     call phys_scv_hanpan(jcapr,deltim,delprsi,prsl,prsi,prslk,prsik,zl,zi,    &
               phil,ncloud_,                                                   &
#ifdef DET_HYDRO
               gq0(1,kcloud_),                                                 &
#endif
               gq0(1,1),gt0(1,1),rain2,                                        &
               kbot,ktop,icps,slmsk(1,lat),vvel,hpbl(1,lat),                   &
               gu0(1,1),gv0(1,1),                                              &
               rcs,                                                            &
               hflx,evapol,                                                    &
               ids,ide, jds,jde, kds,kde,                                      &
               ims,ime, jms,jme, kms,kme,                                      &
               its,ite, jts,jte, kts,kte)
     do i = its,ite
       rainc(i)      = frain*(rain1(i)+rain2(i))
       raincps(i,lat) = raincps(i,lat)+rainc(i)
     enddo
#endif
#if !defined(CCMSCV) && !defined(MRFSCV) && !defined(HANPANSCV) && !defined(GRIMSCV)
     do i = its,ite
       rainc(i)      = frain*rain1(i)
       raincps(i,lat) = raincps(i,lat)+rainc(i)
     enddo
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       call print_maxmin_seven(gt0,ite,ime,kte,kts,kte,                        &
                              'gt0 aft phys_scv')
       call print_maxmin_seven(gq0,ite,ime,ktotal,1,ktotal,                    &
                              'gq0 aft phys_scv')
       call print_maxmin_seven(gu0,ite,ime,kte,kts,kte,                        &
                              'gu0 aft phys_scv')
       call print_maxmin_seven(gv0,ite,ime,kte,kts,kte,                        &
                              'gv0 aft phys_scv')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, 'after scv gu0',ips,ipe,gu0(ndbg1,kdbg1)
   print*, 'after scv gv0',ips,ipe,gv0(ndbg1,kdbg1)
   print*, 'after scv gt0',ips,ipe,gt0(ndbg1,kdbg1)
   print*, 'after scv gq0',ips,ipe,gq0(ndbg1,kdbg1)
#endif
#ifdef DG3
     call diag_3d_archive(ite,ime,gt0(1,1),frain,kdtshal,gda)
     call diag_3d_archive(ite,ime,gq0(1,1),frain,kdqshal,gda)
#endif
!
!  do convective cloud cover
!
     call cps_cloudiness_slingo(ite,ime,                                       &
               clstp,rainc,kbot,ktop,                                          &
#ifndef RMP
               cv(1,lat),cvb(1,lat),cvt(1,lat),count)
#else
               cv(1,lat),cvb(1,lat),cvt(1,lat),1.)
#endif
#ifdef DG3
     call diag_3d_archive(ite,ime,gt0(1,1),-frain,kdtlarg,gda)
#endif
!
!  do large-scale precipitation physics
!
#ifdef WSM1
     call phys_mps_wsm1(deltim,gt0(1,1),gq0(1,1),                              &
                   prsl,delprsi,rain1,lat,                                     &
                   cp_,g_,rd_,rv_,hvap_,                                       &
#ifdef SAS_DIAG
               dlt,dlq,dlh,                                                    &
#endif
#ifdef RASV2
               clw,                                                            &
#endif
               ids,ide, jds,jde, kds,kde,                                      &
               ims,ime, jms,jme, kms,kme,                                      &
               its,ite, jts,jte, kts,kte)

#endif
#ifdef WSM2
#ifdef ZHAO_CLD2
!
!  check for first time step of the segment run for wsm2
!
     if (shour.eq.0.) then
       do k = 1,kte
         do i = 1,ite
           tp(i,k,lat)  = gt0(i,k)
           qp(i,k,lat)  = max(gq0(i,k),epsq)
           tp1(i,k,lat) = gt0(i,k)
           qp1(i,k,lat) = max(gq0(i,k),epsq)
         enddo
       enddo
!
       do i = 1,ite
         psp(i,lat)  = psexp(i)
         psp1(i,lat) = psexp(i)
       enddo
     endif
!
     call phys_mps_cld2_cloud(ite,ime,kte,deltim,prsi,prsl,                    &
               gq0(1,1),gq0(1,kcloud_),gt0(1,1),lat,                           &
               tp(1,1,lat),qp(1,1,lat),psp(1,lat),                             &
               tp1(1,1,lat),qp1(1,1,lat),psp1(1,lat),rcs)
!
     call phys_mps_cld2_precip(ite,ime,kte,deltim,delprsi,prsi,prsl,           &
               gq0(1,1),gq0(1,kcloud_),gt0(1,1),rain1,lat,rcs)
#else
     call phys_mps_wsm2(ite,ime,kte,deltim,                                    &
               gt0(1,1),gq0(1,1),gq0(1,kcloud_),ncloud_,                       &
               delprsi,prsi,prsl,                                              &
#ifdef SMP
               rain1,lat,xgrs,kdt,fhour)
#else
               rain1,lat,vvel,kdt,fhour)
#endif
#endif /* ZHAO_CLD2 end */
#endif /* WSM2 */
#ifdef WSM3
     call phys_mps_wsm3(gt0(1,1),gq0(1,1),gq0(1,kcloud_),ncloud_,              &
                  ww(1,1),prsl_pa(1,1),                                        &
                  delz(1,1),deltim,g_,cp_,cvap_,rd_,rv_,t0c_,                  &
                  fv_,rdorv_,qmin_,hsub_,hvap_,hfus_,rhoair0_,rhoh2o_,         &
                  cliq_,cice_,psat_,                                           &
                  lat, rain1,rainncv,snow1,snowncv,                            &
                  sr,                                                          &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte)
#endif
!
#ifdef WSM5
     call phys_mps_wsm5(gt0(1,1),gq0(1,1),gq0(1,kcloud_),ncloud_,prsl_pa(1,1), &
                  delz(1,1),deltim,g_,cp_,cvap_,rd_,rv_,t0c_,                  &
                  fv_,rdorv_,qmin_,hsub_,hvap_,hfus_,rhoair0_,rhoh2o_,         &
                  cliq_,cice_,psat_,                                           &
                  lat, rain1,rainncv,sr,                                       &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                  snow1,snowncv)
#endif
#ifdef WSM6
     call phys_mps_wsm6(gt0(1,1),gq0(1,1),gq0(1,kcloud_),ncloud_,prsl_pa(1,1), &
                  delz(1,1),deltim,g_,cp_,cvap_,rd_,rv_,t0c_,                  &
                  fv_,rdorv_,qmin_,hsub_,hvap_,hfus_,rhoair0_,rhoh2o_,         &
                  cliq_,cice_,psat_,                                           &
                  lat, rain1,rainncv,sr,                                       &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                  snow1, snowncv, graupel1, graupelncv)
#endif
#ifdef WDM5
     call phys_mps_wdm5(gt0,gq0(1,1),gq0(1,kcloud_),ncloud_,gq0(1,kwsize_),    &
                  prsl_pa(1,1),delz(1,1),deltim,g_,cp_,cvap_,rd_,rv_,t0c_,     &
                  fv_,rdorv_,qmin_,hsub_,hvap_,hfus_,rhoair0_,rhoh2o_,         &
                  cliq_,cice_,psat_,                                           &
                  lat, rain1,rainncv,sr,                                       &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                  snow1, snowncv)
#endif
#ifdef WDM6
     call phys_mps_wdm6(gt0,gq0(1,1),gq0(1,kcloud_),ncloud_,gq0(1,kwsize_),    &
                  prsl_pa(1,1),delz(1,1),deltim,g_,cp_,cvap_,rd_,rv_,t0c_,     &
                  fv_,rdorv_,qmin_,hsub_,hvap_,hfus_,rhoair0_,rhoh2o_,         &
                  cliq_,cice_,psat_,                                           &
                  lat, rain1,rainncv,sr,                                       &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                  snow1, snowncv, graupel1, graupelncv)
#endif
#ifdef DG3
     call diag_3d_archive(ite,ime,gt0(1,1), frain,kdtlarg,gda)
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       call print_maxmin_seven(gt0,ite,ime,kte,kts,kte,                        &
                              'gt0 aft phys_mps')
       call print_maxmin_seven(gq0,ite,ime,ktotal,1,ktotal,                    &
                              'gq0 aft phys_mps')
       call print_maxmin_seven(gu0,ite,ime,kte,kts,kte,                        &
                              'gu0 aft phys_mps')
       call print_maxmin_seven(gv0,ite,ime,kte,kts,kte,                        &
                              'gv0 aft phys_mps')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
#ifdef NIM_DBGM
   print*, 'after mps gu0',ips,ipe,gu0(ndbg1,kdbg1)
   print*, 'after mps gv0',ips,ipe,gv0(ndbg1,kdbg1)
   print*, 'after mps gt0',ips,ipe,gt0(ndbg1,kdbg1)
   print*, 'after mps gq0',ips,ipe,gq0(ndbg1,kdbg1)
#endif
!
     do i = its,ite
#if defined(WSM1) || defined(WSM2) || defined(WSM3) ||                         \
    defined(WSM5) || defined(WSM6) || defined(WDM5) || defined(WDM6)
       rainl(i)      = frain*rain1(i)
#endif
       rain(i)        = rainc(i)+rainl(i)
       raintot(i,lat) = raintot(i,lat)+rain(i)
!soojin_couple

#ifdef NOAHYDRO
       raintot2(i,lat)= raintot2(i,lat)+prcp(i,lat)
#endif
     enddo
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       call print_maxmin_seven(rainl,ite,ime,1,1,1,                            &
                                'rainl aft phys_mps')
       call print_maxmin_seven(rainc,ite,ime,1,1,1,                            &
                                'rainc aft phys_mps')
       call print_maxmin_seven(raincps(1,lat),ite,ime,1,1,1,                   &
                                'raincps aft phys_mps')
       call print_maxmin_seven(raintot(1,lat),ite,ime,1,1,1,                   &
                                'raintot aft phys_mps')
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
!
!  estimate t850 for rain-snow decision
!
     do i = its,ite
       pk(i)   = prsl(i,1)
       t850(i) = gt0(i,1)
     enddo
!
     do k = kts,kte-1
       do i = its,ite
         pkp(i) = prsl(i,k+1)
         if(pk(i) .gt. p850 .and. pkp(i) .le. p850) then
           t850(i) = gt0(i,k) - (pk(i) - p850) / (pk(i) - pkp(i)) *            &
                     (gt0(i,k) - gt0(i,k+1))
         endif
         pk(i) = pkp(i)
       enddo
     enddo
!
!  factor=weighted mean temperature
!
     do i = its,ite
#ifdef MRGLSM
       prcp(i,lat) = rain(i)
       srflag(i,lat) = 0.
       if(t850(i).le.critsnow) srflag(i,lat)=1.
#else
       if(t850(i).le.critsnow) then
         snowfl(i) = rain(i)
         if(slmsk(i,lat).ne.0.) then
           snoweq(i,lat)   = snoweq(i,lat)+1.e3*rain(i)
           snowfall(i,lat) = snowfall(i,lat)+rain(i)
         endif
         rain(i)=0.
       else
         snowfl(i)=0.
       endif
#endif
     enddo

#ifdef OSU
!
!  update soil moisture and canopy water after precipitation has been computed
!
     call lsm_osu_precip(ite,ime,lsoil_,                                       &
         rhscnpy,rhsmc,ai,bi,cci,smsoil,                                       &
         slmsk(1,lat),canopy(1,lat),rain,runof,snowmt,                         &
#ifdef HYDRO
         zsoil,soiltyp,sigmaf,deltim,lat,hydrow)
#else
         zsoil,soiltyp,sigmaf,deltim,lat)
#endif
#endif
!
!  total runoff is composed of drainage into water table and
!  runoff at the surface and is accumulated in unit of meters
!
     do i = its,ite                          
       runoff(i,lat) = (drain(i) + runof(i)) * dtf / 1000.+ runoff(i,lat)
       bgrun(i,lat) = drain(i) * dtf / 1000.+ bgrun(i,lat)
       snowmelt(i,lat) = snowmelt(i,lat) +                                     &
#ifdef HYDRO
                       snowmt(i)*dtf*hfus_*rhoh2o_
#else
#ifdef NOAHYDRO
                       snowmt(i)*hfus_*rhoh2o_
#else
                       snowmt(i)*deltim*hfus_*rhoh2o_
#endif
#endif
     enddo                                          
#ifdef RIVER
!
!  for river discharge calculation, snapshots of total runoff are stored
!
     do i = its,ite                          
       trunof(i,lat)=runof(i)+drain(i)
     enddo                                          
#endif
!
!  return updated smsoil and stsoil to global arrays
!
     do k = 1,lsoil_
       do i = its,ite
         smc(i,lat,k) = smsoil(i,k)
#ifndef VIC
         stc(i,lat,k) = stsoil(i,k)
#endif
#ifdef NOALSM1
         slc(i,lat,k) = slsoil(i,k)
#endif
#ifdef VICLSM1
         sic(i,lat,k) = sicld(i,k)
#endif
       enddo
     enddo
#ifdef VICLSM1
!
     do k = 1,nsoil_
       do i = its,ite
         stc(i,lat,k) = stsoil(i,k)
       enddo
     enddo
!
     do i = its,ite
       lstsn(i,lat) = lstsnld(i)
     enddo
#endif
!
   endif   !!! end of physics computations for inistp .ne. 1
#ifdef DGP
!
!  grid point diagnostics
!
   call diag_point_arrange(lat,rcl,slmsk(1,lat),psexp,                         &
            tg3(1,lat),snoweq(1,lat),radsl,dlwsf1,                             &
#ifndef OSULSM1
            tsea(1,lat),qss,vegtype,gflx,z0cm(1,lat),cd,cdq,                   &
#else
            tsea(1,lat),qss,plantr(1,lat),gflx,z0cm(1,lat),cd,cdq,             &
#endif
            rnet,hflx,stsoil,                                                  &
            canopy(1,lat),drain,smsoil,runof,cld1d,                            &
            u10m(1,lat),v10m(1,lat),t2m(1,lat),q2m(1,lat),                     &
            hpbl(1,lat),gamt,gamq,                                             &
            dqsfc1,dtsfc1,dusfc1,dvsfc1,                                       &
            dusfcg,dvsfcg,                                                     &
            rainc,rainl,                                                       &
            gu0(1,1),gv0(1,1),                                                 &
            gt0(1,1),gq0(1,1),ntotal_,                                         &
            hsw,hlw(1,1,lat),vvel,                                             &
#ifdef EXPLICIT_CLOUDINESS
            qcicps,qrscps,                                                     &
#endif
            snowmt,snowev,snowfl)
#endif
!
#ifndef NIM
   do i = its,ite
     pwat(i,lat)=0.
   enddo
#endif
!
!  restore dry temp back to virtual temperature
!
#undef VTEMP
#ifdef RMP
#define VTEMP
#endif
#ifdef PSPLIT
#define VTEMP
#endif
#ifdef VTEMP
   do k = kts,kte
     do i = its,ite
       w2(i,k)  = 1.+fv_*qgrs(i,k)
       work1(i) = 0.
       if(ncloud_.ge.1) then
         do ic = icloud_,nwmass_
           kc = (ic-1)*kte + k
           work1(i) = max(qgrs(i,kc),qmin_) + work1(i)
         enddo
         w2(i,k) = w2(i,k) - work1(i)
       endif
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       tgrs(i,k) = tgrs(i,k)*w2(i,k)
     enddo
   enddo
#endif
   do k = kts,kte
     do i = its,ite
       w2(i,k)  = 1.+fv_*gq0(i,k)
       work1(i) = 0.
       if(ncloud_.ge.1) then
         do ic = icloud_,nwmass_
           kc = (ic-1)*kte + k
           work1(i) = max(gq0(i,kc),qmin_) + work1(i)
         enddo
         w2(i,k) = w2(i,k) - work1(i)
       endif
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       gt0(i,k) = gt0(i,k)*w2(i,k)
#ifdef NONHYD_HYD
       gp0(i,k) = pgrs(i,k)
#endif /* NONHYD_HYD */
     enddo
   enddo
!
#ifndef NIM
   do k = kts,kte
     do i = its,ite
       pwat(i,lat) = pwat(i,lat)+delprsi(i,k)/prsi(i,1)*gq0(i,k)
       qu00(i)     = qu00(i)+del(k)*gq0(i,k)*gu0(i,k)*rcs
       qv00(i)     = qv00(i)+del(k)*gq0(i,k)*gv0(i,k)*rcs
     enddo
   enddo
!
   do i = its,ite
     pwat(i,lat)=pwat(i,lat)*prsi(i,1)*(1.e2/g_)
#ifdef HYDRO
     qull(i,lat)=qull(i,lat)+qu00(i)*prsi(i,1)*(1.e2/g_)*dtf
     qvll(i,lat)=qvll(i,lat)+qv00(i)*prsi(i,1)*(1.e2/g_)*dtf
#else
     qull(i,lat)=qull(i,lat)+qu00(i)*prsi(i,1)*(1.e2/g_)*deltim
     qvll(i,lat)=qvll(i,lat)+qv00(i)*prsi(i,1)*(1.e2/g_)*deltim
#endif
   enddo
#ifndef RMP
#if !defined(SMP) && !defined(NIM)
!
!  quasi-pressure surface horizontal diffusion correction for t and q
!
#ifndef NISLQ
   call phys_grid_diffusion(deltim,prsl,pslap,gt0,gq0,ntotal_,                 &
                              ids,ide, jds,jde, kds,kde,                       &
                              ims,ime, jms,jme, kms,kme,                       &
                              its,ite, jts,jte, kts,kte)
#endif /* ~NISLQ end */
#endif
#endif /* ~RMP end */
#endif /* ~NIM end */
#ifdef NISLQ_CHECK
!
! check negative after phys_main_solver
!
   do k = 1,levh_
     do i = its,ite
       if(gq0(i,k).lt.0.) print *,'negative after phys_main_solver',k,gq0(i,k)
     enddo
   enddo
#endif /* NISLQ_CHECK end */
#ifdef NISLQ_PHYS
!
! replace negative with zero
!
   do k = 1,levh_
     do i = its,ite
       if(gq0(i,k).lt.0.) gq0(i,k) = 0.
     enddo
   enddo
#endif
#ifndef RMP
!
   do k = kts,kte
     do i = its,ite
#ifndef PSPLIT
       ggu0(i,k) = gu0(i,k)
       ggv0(i,k) = gv0(i,k)
       ggt0(i,k) = gt0(i,k)
#else
       ggu0(i,k) = (gu0(i,k)-ugrs(i,k))/dt2
       ggv0(i,k) = (gv0(i,k)-vgrs(i,k))/dt2
       ggt0(i,k) = (gt0(i,k)-tgrs(i,k))/dt2
#endif
     enddo
   enddo
!
   do k = kts,levh_
     do i = its,ite
#ifndef PSPLIT
       ggq0(i,k)=gq0(i,k)
#else
       ggq0(i,k)=(gq0(i,k)-qgrs(i,k))/dt2
#endif
     enddo
   enddo
#else
#ifdef PSPLIT
!
   do k = kts,kte
     do i = its,ite
       ggu0(i,k)=(gu0(i,k)-ugrs(i,k))/dt2
       ggv0(i,k)=(gv0(i,k)-vgrs(i,k))/dt2
       ggt0(i,k)=(gt0(i,k)-tgrs(i,k))/dt2
     enddo
   enddo
!
   do k = 1,levh_
     do i = its,ite
       ggq0(i,k)=(gq0(i,k)-qgrs(i,k))/dt2
     enddo
   enddo
#endif
#endif /* RMP end */
#ifdef NIM
#ifdef NIM_BIN
   wor1d=-99999
   wor2d=-99999
   do i = its,ite
     wor1d(i,1)=tsea(i,1)
     wor1d(i,2)=smc(i,1,1)
     wor1d(i,3)=snoweq(i,1)
     wor1d(i,4)=stc(i,1,1)
     wor1d(i,5)=tg3(i,1)
     wor1d(i,6)=z0cm(i,1)
     wor1d(i,7)=albedo(i,1,1)
     wor1d(i,8)=albedo(i,1,2)
     wor1d(i,9)=albedo(i,1,3)
     wor1d(i,10)=albedo(i,1,4)
     wor1d(i,11)=slmsk(i,1)
     wor1d(i,12)=vfrac(i,1)
     wor1d(i,13)=canopy(i,1)
     wor1d(i,14)=f10m(i,1)
     wor1d(i,15)=vtype(i,1)
     wor1d(i,16)=stype(i,1)
     wor1d(i,17)=u10m(i,1)
     wor1d(i,18)=v10m(i,1)
     wor1d(i,19)=evap(i)
     wor1d(i,20)=dtsfc1(i)
     wor1d(i,21)=dqsfc1(i)
     wor1d(i,22)=dlwsf1(i)
     wor1d(i,23)=ulwsf1(i)
     wor1d(i,24)=rnet(i)
   enddo
   do k=1,levs_
     do i = its,ite
       wor2d(i,k,1)=ggu0(i,k)
       wor2d(i,k,2)=ggv0(i,k)
       wor2d(i,k,3)=ggt0(i,k)
     enddo
   enddo
#endif
#ifdef NIM_DIAG
!
   toth_pbl=0. ; tots_pbl=0.
   do i = its,ite
     evap_nim(i)   = evap(i)      ! evaporation from sfc
     shflx_nim(i)  = dtsfc1(i)    ! sensible heat flux at every time step
     lhflx_nim(i)  = dqsfc1(i)    ! latent heat flux at every time step
     dlwsfc_nim(i) = dlwsf1(i)    ! downward longwave flux at the sfc
     ulwsfc_nim(i) = ulwsf1(i)    ! upward longwave flux at the sfc
     rnet_nim(i)   = rnet(i)      ! net sfc radiation
     tots_pbl(i)   = shflx_nim(i) ! for energy budget related to pbl
   enddo
   do i = its,ite
     do k = kts,kte
       swhr_nim(i,k)  = hsw(i,k)        ! sw heating rate * xmu
       lwhr_nim(i,k)  = hlw(i,k,lat)    ! lw heating rate
       xmu_nim(i)     = xmu(i)
       toth_pbl(i)    = ttend_nim(i,k)*cp_/g_*(prsi_pa(i,k)-prsi_pa(i,k+1))/prslk(i,k)
     enddo
   enddo
!
#ifdef NIM_DBGM
   if(ncpu.eq.1) then  ! for 1cpu
     print*, 'SPBL 10', tots_pbl(10),toth_pbl(10),(toth_pbl(10)-tots_pbl(10))/tots_pbl(10)
     print*, 'SPBL 100', tots_pbl(100),toth_pbl(100),(toth_pbl(100)-tots_pbl(100))/tots_pbl(100)
     print*, 'SPBL 1000', tots_pbl(1000),toth_pbl(1000),&
                       (toth_pbl(1000)-tots_pbl(1000))/tots_pbl(1000)
   endif
#endif  /* NIM_DBGM */
#endif  /* NIM_DIAG */
#endif  /* NIM */
!
#ifdef CHEM
   ! save 2-D variables for chemistry
   do i = its,ite
     c_eflux(i,lat) = dqsfc1(i)
     c_hflux(i,lat) = dtsfc1(i)
     c_lwi(i,lat) = slmsk(i,lat)
     c_pbl(i,lat) = hpbl(i,lat)
     c_preacc(i,lat) = rain(i) / dtf * 1d3 * 86400d0
     c_precon(i,lat) = rainc(i) / dtf * 1d3 * 86400d0
     c_ps1(i,lat) = psurf(i,lat) * 1d1
     c_snomas(i,lat) = snoweq(i,lat)
     c_ts(i,lat) = tsea(i,lat)
     c_u10m(i,lat) = u10m(i,lat)
     c_ustar(i,lat) = uustar(i,lat)
     c_v10m(i,lat) = v10m(i,lat)
     c_z0(i,lat) = z0cm(i,lat) * 1d-2
   enddo
   ! save 3-D variables for chemistry
   do k = kts,kte
     do i = its,ite
       c_cmfmc(i,lat,k+1) = fcu(i,k) + fcd(i,k)
       c_moistq(k,i,lat) = ( c_moistq(k,i,lat) - ggq0(i,k) ) / dt2
       c_sphu(i,lat,k) = ggq0(i,k) * 1d3
       c_t(i,lat,k) = ggt0(i,k)
       c_uwnd(i,lat,k) = ggu0(i,k)
       c_vwnd(i,lat,k) = ggv0(i,k)
     enddo
   enddo
#endif
!
   return
   end subroutine phys_main_solver
