#include <define.h>
   subroutine lsm_driver(ims2,imx2,kmx,                                        &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,soiltyp,sigmaf,                                   &
#ifndef OSULSM1
                  vegtype,canopy,                                              &
#else
                  canopy,                                                      &
#endif
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,cm,ch,rcl,                                       &
                  prsi1,prsl1,prsik1,prslk1,                                   &
                  zl1,slimsk,inistp,lat,                                       &
                  ims,ime,its,ite,                                             &
#ifdef NOAHYDRO
                  drain,evap,hflx,ep,wind,ec,                                  &
#else
#ifdef OSU
                  rhscnpy,rhsmc,aim,bim,cim,                                   &
#ifdef OSULSM1
                  plantr,cowave,                                               &
#endif
                  drain,evap,hflx,ep,wind)
#else
                  drain,evap,hflx,ep,wind,                                     &
#endif
#endif
#ifdef NOALSM1
                  snwdph,slc,snoalb,slope,shdmin,shdmax,                       &
#endif
#ifdef VICLSM1
                  nsl,month,dswflx,binf,ds,dsm,ws,cef,                         &
                  expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,             &
                  wpwld,smrld,smxld,sicld,dpnld,sxnld,epnld,bbnld,             &
                  apnld,btnld,gmnld,flaild,vrtld,lstsnld,                      &
                  silz,snwz,csno,rsno,tsf,tpk,sfw,pkw,gheat,                   &
#endif
#ifdef VICLSM2
                  nsl,msub,month,dswflx,binf,ds,dsm,ws,cef,                    &
                  expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,             &
                  wpwld,smrld,smxld,dpnld,sxnld,epnld,bbnld,apnld,             &
                  btnld,gmnld,                                                 &
                  nvegld,flaild,vfrld,vtypld,cnpld,snold,csnld,                &
                  rsnld,tsfld,tpkld,sfwld,pkwld,lstsnld,vrtld,                 &
                  smcld, sicld,stcld,silz,snwz,                                &
#endif
#ifndef OSU
                  snowfl, runoff, precip, srflag)
#endif
!-------------------------------------------------------------------------------
!
! subroutine: lsm_driver
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
   use constant, only  : cal_
   use paramodel, only : ILOTS,levs_
#ifdef DFS
   use dfsvar, only    : iope
#else
   use comio, only     : iope        ! RMP, GSM share
#endif
#ifdef CRAY_THREAD
!fpp$ noconcur r
!fpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   integer              ::  ims2,imx2,kmx
! passing array
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  zl1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2)
   real                 ::  dlwflx(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),tg3(imx2),canopy(imx2)
#ifndef VIC
   real                 ::  stc(imx2,kmx)
#else
   real                 ::  stc(imx2,nsl),dswflx(imx2)
#endif
!
   integer              ::  soiltyp(imx2),vegtype(imx2)
   real                 ::  z0cm(imx2),gflux(imx2),slimsk(imx2)
   real                 ::  drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),ep(imx2),wind(imx2)
#ifdef NOAHYDRO
   real                 ::  ec(imx2)
#endif
   real                 ::  runoff(imx2)
   real                 ::  precip(imx2), srflag(imx2),snowfl(imx2)
#ifdef OSU
   real                 ::  rhscnpy(imx2),rhsmc(imx2,kmx)
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
#ifdef OSULSM1
   real                 ::  plantr(imx2),cowave
#endif
#endif
#ifdef NOALSM1
   integer              ::  slope(imx2)
   real                 ::  snoalb(imx2), shdmin(imx2), shdmax(imx2)
   real                 ::  slc(imx2,kmx), snwdph(imx2)
#endif
#ifdef VICLSM1
   real                 ::  binf(imx2),ds(imx2),dsm(imx2),ws(imx2),cef(imx2)
   real                 ::  expld(imx2,kmx)
   real                 ::  kstld(imx2,kmx) 
   real                 ::  dphld(imx2,kmx),bubld(imx2,kmx) 
   real                 ::  qrtld(imx2,kmx),bkdld(imx2,kmx) 
   real                 ::  sldld(imx2,kmx),wcrld(imx2,kmx) 
   real                 ::  wpwld(imx2,kmx),smrld(imx2,kmx) 
   real                 ::  smxld(imx2,kmx),sicld(imx2,kmx) 
   real                 ::  dpnld(imx2,nsl),sxnld(imx2,nsl) 
   real                 ::  epnld(imx2,nsl),bbnld(imx2,nsl) 
   real                 ::  apnld(imx2,nsl),btnld(imx2,nsl) 
   real                 ::  gmnld(imx2,nsl)
   real                 ::  flaild(imx2), vrtld(imx2,kmx)
   integer              ::  lstsnld(imx2)
   real                 ::  silz(imx2),snwz(imx2),csno(imx2),rsno(imx2)
   real                 ::  tsf(imx2),tpk(imx2),sfw(imx2),pkw(imx2),gheat(imx2)
#endif
#ifdef VICLSM2
   real                 ::  binf(imx2),ds(imx2),dsm(imx2),ws(imx2),cef(imx2)
   real                 ::  expld(imx2,kmx)
   real                 ::  kstld(imx2,kmx),
   real                 ::  dphld(imx2,kmx),bubld(imx2,kmx),
   real                 ::  qrtld(imx2,kmx),bkdld(imx2,kmx),
   real                 ::  sldld(imx2,kmx),wcrld(imx2,kmx),
   real                 ::  wpwld(imx2,kmx),smrld(imx2,kmx),
   real                 ::  smxld(imx2,kmx),
   real                 ::  dpnld(imx2,nsl),sxnld(imx2,nsl),
   real                 ::  epnld(imx2,nsl),bbnld(imx2,nsl),
   real                 ::  apnld(imx2,nsl),btnld(imx2,nsl),
   real                 ::  gmnld(imx2,nsl)
!
   integer              ::  nvegld(imx2)
   real                 ::  flaild(imx2,msub)
   real                 ::  vfrld(imx2,msub)
   integer              ::  vtypld(imx2,msub)
   real                 ::  cnpld(imx2,msub),snold(imx2,msub),csnld(imx2,msub)
   real                 ::  rsnld(imx2,msub),tsfld(imx2,msub),tpkld(imx2,msub)
   real                 ::  sfwld(imx2,msub),pkwld(imx2,msub)
   integer              ::  lstsnld(imx2,msub)
   real                 ::  vrtld(imx2,kmx,msub)
   real                 ::  smcld(imx2,kmx,msub)
   real                 ::  sicld(imx2,kmx,msub)
   real                 ::  stcld(imx2,nsl,msub)
   real                 ::  silz(imx2),snwz(imx2)
!
#endif
   real                 ::  delt,rcl
   integer              ::  inistp,lat
   integer              ::  ims,ime,its,ite
!
! local array
!
   integer              ::  im,km
   integer              ::  i, j, k
#ifdef NOALSM1
!
!  declare zsoil(4)
!
   real                 ::  zsoil_noah(4)
   data zsoil_noah/-0.1, -0.4, -1.0, -2.0/
#endif
!
#ifdef DBG
   character(len=4)     ::  sfcty
#ifdef OSULSM1
      sfcty = 'osu1'
#endif
#ifdef OSULSM2
      sfcty = 'osu2'
#endif
#ifdef NOALSM1
      sfcty = 'noa1'
#endif
#ifdef VICLSM1
      sfcty = 'vic1'
#endif
#ifdef VICLSM2
      sfcty = 'vic2'
#endif
#endif
!
   im = ims2
   km = kmx
!
!-------------------------------------------------------------------------------
!
! initialize
!
#ifdef DBG
   if(iope) print *,'--- enter sfcdrv ---'
#endif      
   runoff=0.
   drain=0.
   dm=0.
   evap=0.
   hflx=0.
   ep=0.
   snowfl=0.
   snowmt=0.
   snowev=0.
   gflux=0.
!
! set up soil layer configuration
!
   do i = 1,im
     if(slimsk(i).eq.0.) then
       zsoil(i,1) = 0.
     elseif(slimsk(i).eq.1.) then
       zsoil(i,1) = -.10
     else
       zsoil(i,1) = -3. / km
     endif
   enddo
!
#ifdef DBG
   if (iope) write(6,*)'in sfcdrv before compute zsoil'
!
#endif
   do k = 2,km
     do i = 1,im
       if(slimsk(i).eq.0.) then
         zsoil(i,k) = 0.
       elseif(slimsk(i).eq.1.) then
#ifdef OSU
!
!  use default soil layer configuration
!
         zsoil(i,k) = zsoil(i,k-1)+ (-2. - zsoil(i,1)) / (km - 1)
#endif
#ifdef NOALSM1
!
!  override soil layer configuration using noah_zsoil
!
         zsoil(i,k) = zsoil_noah(k)
#endif
#ifdef VIC
         zsoil(i,k) = -dphld(i,k)
#endif
       else
         zsoil(i,k) = - 3. * float(k) / float(km)
       endif
     enddo
   enddo
!
!  surface energy/water balance over land and seaice
!
#ifdef DBG
   if (iope) write(6,*) ' in sfcdrv surface type: '//sfcty
!
#endif
#ifdef OSU
   call lsm_osu_main(ims2,imx2,kmx,                                            &
                  t1,q1,snoweq,tskin,qsurf,                                    &
#ifdef OSULSM1
                  smc,stc,dm,soiltyp,sigmaf,canopy,                            &
#else
                  smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,                    &
#endif
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,cm,ch,                                           &
                  rhscnpy,rhsmc,aim,bim,cim,                                   &
#ifdef OSULSM1
                  plantr,cowave,                                               &
#endif
                  prsi1,prsl1,prsik1,prslk1,                                   &
                  rcl,zl1,slimsk,inistp,lat,                                   &
                  drain,evap,hflx,ep,wind,                                     &
                  ims,ime,its,ite)
#endif /* OSU end */
!
#ifdef NOALSM1
#ifdef DBG
   if (iope) write(6,*)' in sfcdrv before lsm_noah_main'
!
#endif
!
!  call Noah LSM over land (use sfc_noah_ver2.F)
!
   call lsm_noah_main(ims2,imx2,kmx,                                           &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,                    &
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,cm,ch,rcl,                                       &
                  prsi1,prsl1,prsik1,prslk1,                                   &
                  zl1,slimsk,inistp,lat,                                       &
                  drain,evap,hflx, ep, wind,                                   &
                  snwdph, slc, snoalb, slope, shdmin, shdmax,                  &
#ifndef NOAHYDRO
                  snowfl, runoff, precip, srflag)
#else
                  snowfl, runoff, precip, srflag, ec)
#endif
#endif /* NOALSM1 end */
!
#ifdef VIC
   do i = 1,imx2
     dswflx(i)=-dswflx(i)*cal_*1.e04/60.
   enddo
#endif
!
#ifdef VICLSM1
!
!  call UW VIC over land (use fortran-version VIC)
!
#ifdef DBG
   if (iope) write(6,*)' in sfcdrv before sfc_vic1'
!
#endif
   call lsm_vic1_main(ims2,imx2,kmx,lat,                                       &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,sigmaf,vegtype,canopy,                            &
                  dlwflx,dswflx,snowmt,snowev,delt,z0cm,tg3,                   &
                  gflux,zsoil,cm,ch,                                           &
                  prsi1,prsl1,prsik1,prslk1,zl1,slimsk,                        &
                  drain,evap,hflx, ep, wind,                                   &
                  nsl,month,binf,ds,dsm,ws,cef,                                &
                  expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,             &
                  wpwld,smrld,smxld,sicld,dpnld,sxnld,epnld,bbnld,             &
                  apnld,btnld,gmnld,flaild,vrtld,lstsnld,                      &
                  silz,snwz,csno,rsno,tsf,tpk,sfw,pkw,gheat,                   &
#ifdef SMP_RA2SFC
                  slrad, snowfl, runoff, precip)
#else
                  snowfl, runoff, precip)
#endif
#endif /* VICLSM1 end */
! 
#ifdef VICLSM2
!
!  call UW VIC over land (use fortran-version VIC)
!
#ifdef DBG
   if (iope) write(6,*)' in sfcdrv before sfc_vic2'
#endif
!
   call lsm_vic2_main(ims2,imx2,kmx,lat,                                       &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,sigmaf,vegtype,canopy,                            &
                  dlwflx,dswflx,snowmt,snowev,delt,z0cm,tg3,                   &
                  gflux,zsoil,cm,ch,                                           &
                  prsi1,prsl1,prsik1,prslk1,slimsk,                            &
                  drain,evap,hflx, ep, wind,                                   &
                  nsl,msub,month,binf,ds,dsm,ws,cef,                           &
                  expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,             &
                  wpwld,smrld,smxld,dpnld,sxnld,epnld,bbnld,apnld,             &
                  btnld,gmnld,                                                 &
                  nvegld,flaild,vfrld,vtypld,cnpld,snold,csnld,                &
                  rsnld,tsfld,tpkld,sfwld,pkwld,lstsnld,vrtld,                 &
                  smcld, sicld,stcld,silz,snwz,                                &
                  snowfl, runoff, precip)
!
#endif
!
!  call OSU LSM over seaice
!
#ifndef OSU
!
#ifdef DBG
   if (iope) write(6,*)' in sfcdrv before noah_vic_seaice'
!
#endif
   call noah_vic_seaice(ims2,imx2,kmx,                                         &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,                    &
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,cm,ch,rcl,                                       &
                  prsi1,prsl1,prsik1,prslk1,                                   &
                  zl1,slimsk,inistp,lat,                                       &
                  drain,evap,hflx,ep,wind,                                     &
                  snowfl, runoff, precip, srflag)
!
!  surface energy/water balance over ocean
!
#ifdef DBG
   if (iope) write(6,*)' in sfcdrv before noah_vic_ocean'
!
#endif
    call noah_vic_ocean(ims2,imx2,kmx,                                         &
                   t1,q1,tskin,qsurf,                                          &
                   dm, cm, ch,                                                 &
                   prsi1,prsl1,prsik1,prslk1,                                  &
                   slimsk,inistp,lat,                                          &
                   evap,hflx,wind)
!
#endif
   return
   end subroutine lsm_driver
