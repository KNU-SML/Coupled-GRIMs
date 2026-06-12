#include <define.h>
   subroutine lsm_common_subprogram
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!      [phys_main_solver] --- [lsm_common_subprogram]
!                                  |
!                                  |--- [lsm_exch_coeff] *
!                                  |--- [lsm_diagnostic] *
!
!         [lsm_driver]    --- [lsm_common_subprogram]
!                                  |
!                                  |--- [noah_vic_ocean] *
!                                  |--- [noah_vic_seaice] *
!                                            |--- [noah_vic_seaicetm] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!   2011-10-01  song-you hong          revisions z0,wind,ustar kim hong (2010)
!
!-------------------------------------------------------------------------------
   end subroutine lsm_common_subprogram
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine lsm_exch_coeff(ims2,imx2,kmx,                                    &
                     u1,v1,t1,q1,tskin,z0cm,z0cmt,                             &
                     cm, ch, rb, rcl,                                          &
                     zl1,prsi1,prsl1,prsik1,prslk1,                            &
                     slimsk,inistp,lat,                                        &
#ifdef VSGD
                     hpbl,rdelx,rainc,                                         &
#endif
                     fm,fh,ustar,wind,                                         &
#ifdef VIC
                     fm10,fh2,rho)
#else
                     fm10,fh2)
#endif
!-------------------------------------------------------------------------------
   use constant, only  : cal_,cp=>cp_,g=>g_,hfus=>hfus_,hvap=>hvap_,           &
                         rd=>rd_,rv=>rv_,sigma=>sbc_,                          &
                         akapa_,qmin8_,cb2pa=>cb2pa_,rhoh2o=>rhoh2o_,          &
                         convrad=>convrad_,elocp=>elocp_,                      &
                         rdorv=>rdorv_,rdorvm1=>rdorvm1_,rvordm1=>rvordm1_
   use paramodel, only : ILOTS,levs_,nstype_,nvtype_
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
#ifdef VSGD
#define VCON
#define VGRD
#undef VCPR
#else
#undef VCON
#undef VGRD
#undef VCPR
#endif
!
#include "abort.h"
!
   real,parameter       ::  charnock=.014,ca=.4
   real,parameter       ::  alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
   real,parameter       ::  a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
   real,parameter       ::  aa1=-1.076,bb1=.7045,cc1=-.05808
   real,parameter       ::  bb2=-.1954,cc2=.009999
   real,parameter       ::  dfsnow=.31,ch2o=4.2e6,csoil=1.26e6
   real,parameter       ::  scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2
   real,parameter       ::  cice=1880.*917.,topt=298.
   real,parameter       ::  oz0min = 1.59e-5, oz0max = 2.85e-3
   real,parameter       ::  ctfil1=.5,ctfil2=1.-ctfil1
   real,parameter       ::  rnu=1.51e-5,arnu=.135*rnu
!
! passing array
!
   integer              ::  ims2,imx2,kmx
   real                 ::  rcl
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  zl1(imx2),u1(imx2),v1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),z0cmt(imx2),plantr(imx2)
   real                 ::  gflux(imx2)
   real                 ::  u10m(imx2),v10m(imx2),t2m(imx2),q2m(imx2)
   real                 ::  slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
   real                 ::  f10m(imx2),drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),rnet(imx2),ep(imx2)
   real                 ::  fm(imx2),fh(imx2),ustar(imx2),wind(imx2)
   real                 ::  dlwflx(imx2)
#ifdef VSGD
   real                 ::  hpbl(imx2)
   real                 ::  rdelx
   real                 ::  rainc(imx2)
#endif
   integer              ::  inistp,lat
   integer              ::  soiltyp(imx2)
   integer              ::  vegtype(imx2)
!
! local array
!
   real                 ::  xrcl,dthv,vconv,vcpr,vgrd
   real                 ::  hl0inf,hltinf,aa,aa0,bb,bb0,fms,fhs,hl0,hlt
   real                 ::  ktsoil
   real                 ::  theta1(ILOTS)
   real                 ::  tv1(ILOTS),tvs(ILOTS)
   real                 ::  z1(ILOTS),thv1(ILOTS)
   real                 ::  rho(ILOTS),qs1(ILOTS)
   real                 ::  qss(ILOTS)
   real                 ::  tsurf(ILOTS)
   real                 ::  q0(ILOTS),cq(ILOTS)
   real                 ::  z0max(ILOTS),ztmax(ILOTS)
   real                 ::  dtv(ILOTS),adtv(ILOTS)
   real                 ::  fm10(ILOTS),fh2(ILOTS)
   real                 ::  hlinf(ILOTS)
   real                 ::  hl1(ILOTS),pm(ILOTS)
   real                 ::  ph(ILOTS)
   real                 ::  hl110(ILOTS),hl12(ILOTS)
   real                 ::  pm10(ILOTS),ph2(ILOTS)
   real                 ::  olinf(ILOTS)
   real                 ::  z0(ILOTS)
   real                 ::  restar(ILOTS),rat(ILOTS)
#ifdef ICE
   real                 ::  fpvs
   external fpvs
#else
   real                 ::  fpvs0
   external fpvs0
#endif
   logical              ::  flag(ILOTS)
   logical              ::  flagsnw(ILOTS)
   integer              ::  im,km,i
!-------------------------------------------------------------------------------
!
! initialize local variables
!
   theta1=0.   ;  tv1=0.   ;  tvs=0.   ;  z1=0.    ;  thv1=0.
   rho=0.      ;  qs1=0.   ;  qss=0.   ;  tsurf=0. ;   q0=0.
   cq=0.       ;  z0max=0. ;  ztmax=0. ;  dtv=0.   ;  adtv=0.
   fm10=0.     ;  fh2=0.   ;  hlinf=0. ;  hl1=0.   ;  pm=0.
   ph=0.       ;  hl110=0. ;  hl12=0.  ;  pm10=0.  ;  ph2=0.
   olinf=0.    ;  z0=0.    ;  restar=0.;  rat=0.
   hl0inf=0.   ;  hltinf=0.;  aa=0.    ;  aa0=0.   ;  bb=0.
   bb0=0.      ;  fms=0.   ;  fhs=0.   ;  hl0=0.   ;  hlt=0.
   ktsoil=0.
!
   im = ims2
   km = kmx
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  prsi(i,1)*1000. is in pascals.
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
!  surface. convert slrad to the civilized unit from langley minute-1 k-4
!  surface roughness length is converted to m from cm
!
   xrcl = sqrt(rcl)
#ifndef VGRD
   vgrd = 0.0
#else
   vgrd = max(0.32 * (max(rdelx/5000.-1.,0.))**.33,0.1)
#endif
   do i = 1,im
     q0(i) = max(q1(i),qmin8_)
     tsurf(i) = tskin(i)
     theta1(i) = t1(i) / prslk1(i) * prsik1(i) 
#ifndef VSGD
     wind(i) = xrcl * sqrt(u1(i) * u1(i) + v1(i) * v1(i))
     wind(i) = max(wind(i),1.)
#else                   /* VSGD */
!
! Convective velocity scale Vc and subgrid-scale velocity Vsg
! following Beljjars (1995, QJRMS) and Mahrt and Sun (1995, MWR)
!
     dthv = ( tsurf(i) - theta1(i) )
#ifdef V31_201203
     if(dthv.gt.0.)  then
#else
     if(dthv.gt.0..and.slimsk(i).ne.0)  then
#endif
       vconv = max(2.0*sqrt(dthv),.1)
     else
       vconv = 0.0
     endif
     vcpr = 0.0
#ifdef VCPR
     xcon = rainc(i)*8640.               ! mm/s => cm/day
     vcpr = alog( 1. + 67.*xcon - 0.48*(xcon**2) )
#endif
     wind(i) = xrcl * sqrt(u1(i) * u1(i) + v1(i) * v1(i))
     wind(i)= wind(i) + sqrt(vconv*vconv+vgrd*vgrd+vcpr*vcpr)
#endif                   /* VSGD */
     tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
     thv1(i) = theta1(i) * (1. + rvordm1 * q0(i))
     tvs(i) = tsurf(i) * (1. + rvordm1 * q0(i))
     rho(i) = (cb2pa * prsl1(i)) / (rd * tv1(i))
#ifdef ICE
     qs1(i) = cb2pa * fpvs(t1(i))
#else
     qs1(i) = cb2pa * fpvs0(t1(i))
#endif
     qs1(i) = rdorv * qs1(i) / (cb2pa * prsl1(i) + rdorvm1 * qs1(i))
     qs1(i) = max(qs1(i), qmin8_)
     q0(i) = min(qs1(i),q0(i))
#ifdef ICE
     qss(i) = cb2pa * fpvs(tsurf(i))
#else
     qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
     qss(i) = rdorv * qss(i) / (cb2pa * prsi1(i) + rdorvm1 * qss(i))
     z0(i) = .01 * z0cm(i)
     z1(i) = zl1(i)
   enddo
!
!  compute stability dependent exchange coefficients
!
   do i = 1,im
     ustar(i) = .1 * wind(i)
     if(slimsk(i).eq.0.) then
#ifndef KIMHONG2010
       ustar(i) = sqrt(g * z0(i) / charnock)
#else
       z0(i) = max(min(z0(i),oz0max),oz0min)
       ustar(i) = max(sqrt(g * max(z0(i)-oz0min,1.e-6) / charnock),0.01)
#endif
     endif
   enddo
!
!  compute stability indices (rb and hlinf)
!
   do i = 1,im
     z0max(i) = min(z0(i),1. * z1(i))
#ifndef Z0T
     ztmax(i) = z0max(i)
     if(slimsk(i).eq.0.) then
       restar(i) = ustar(i) * z0max(i) / vis
       restar(i) = max(restar(i),.000001)
       restar(i) =  log(restar(i))
       restar(i) = min(restar(i),5.)
       restar(i) = max(restar(i),-5.)
       rat(i) = aa1 + bb1 * restar(i) + cc1 * restar(i) ** 2
       rat(i) = rat(i) / (1. + bb2 * restar(i)                                 &
                         + cc2 * restar(i) ** 2)
       ztmax(i) = z0max(i) * exp(-rat(i))
     endif
#else
     restar(i) = ustar(i) * z0max(i) / vis
     restar(i) = max(restar(i),.000001)
     if(slimsk(i).eq.0.) then
       rat(i) = 2.67 * restar(i) ** .25 - 2.57
     else
       restar(i) = min(restar(i),1000.)
       rat(i) = 0.13 * restar(i) ** 0.45
     endif
     rat(i) = min(rat(i),5.)
     ztmax(i) = z0max(i) * exp(-rat(i))
#endif
     z0cmt(i) = ztmax(i) * 100.
   enddo
!
   do i = 1,im
     dtv(i) = thv1(i) - tvs(i)
     adtv(i) = abs(dtv(i))
     adtv(i) = max(adtv(i),.001)
     dtv(i) = sign(1.,dtv(i)) * adtv(i)
     rb(i) = g * dtv(i) * z1(i) / (.5 * (thv1(i) + tvs(i))* wind(i) * wind(i))
     rb(i) = max(rb(i),-5000.)
     if(z0max(i).eq.0.) then
       write(6,*)'z0max=0. at i=',i
       call MPABORT
     elseif(ztmax(i).eq.0.) then
       write(6,*)'ztmax=0. at i=',i
       call MPABORT
     endif
     fm(i) = log((z0max(i)+z1(i)) / z0max(i))
     fh(i) = log((ztmax(i)+z1(i)) / ztmax(i))
     fm10(i) = log((z0max(i)+10.) / z0max(i))
     fh2(i) = log((ztmax(i)+2.) / ztmax(i))
     hlinf(i) = rb(i) * fm(i) * fm(i) / fh(i)
   enddo
!
!  stable case
!
   do i = 1,im
     if(dtv(i).ge.0.) then
       hl1(i) = hlinf(i)
     endif
     if(dtv(i).ge.0..and.hlinf(i).gt..25) then
       hl0inf = z0max(i) * hlinf(i) / z1(i)
       hltinf = ztmax(i) * hlinf(i) / z1(i)
       aa = sqrt(1. + 4. * alpha * hlinf(i))
       aa0 = sqrt(1. + 4. * alpha * hl0inf)
       bb = aa
       bb0 = sqrt(1. + 4. * alpha * hltinf)
       pm(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
       ph(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
       fms = fm(i) - pm(i)
       fhs = fh(i) - ph(i)
       hl1(i) = fms * fms * rb(i) / fhs
     endif
   enddo
!
!  second iteration
!
   do i = 1,im
     if(dtv(i).ge.0.) then
       hl0 = z0max(i) * hl1(i) / z1(i)
       hlt = ztmax(i) * hl1(i) / z1(i)
       aa = sqrt(1. + 4. * alpha * hl1(i))
       aa0 = sqrt(1. + 4. * alpha * hl0)
       bb = aa
       bb0 = sqrt(1. + 4. * alpha * hlt)
       pm(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
       ph(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
       hl110(i) = hl1(i) * 10. / z1(i)
       aa = sqrt(1. + 4. * alpha * hl110(i))
       pm10(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
       hl12(i) = hl1(i) * 2. / z1(i)
       bb = sqrt(1. + 4. * alpha * hl12(i))
       ph2(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
     endif
   enddo
!
!  unstable case
!
!
!  check for unphysical obukhov length
!
   do i = 1,im
     if(dtv(i).lt.0.) then
       olinf(i) = z1(i) / hlinf(i)
       if(abs(olinf(i)).le.50. * z0max(i)) then
         hlinf(i) = -z1(i) / (50. * z0max(i))
       endif
     endif
   enddo
!
!  get pm and ph
!
   do i = 1,im
     if(dtv(i).lt.0..and.hlinf(i).ge.-.5) then
       hl1(i) = hlinf(i)
       pm(i) = (a0 + a1 * hl1(i)) * hl1(i)                                     &
             / (1. + b1 * hl1(i) + b2 * hl1(i) * hl1(i))
       ph(i) = (a0p + a1p * hl1(i)) * hl1(i)                                   &
             / (1. + b1p * hl1(i) + b2p * hl1(i) * hl1(i))
       hl110(i) = hl1(i) * 10. / z1(i)
       pm10(i) = (a0 + a1 * hl110(i)) * hl110(i)                               &
             / (1. + b1 * hl110(i) + b2 * hl110(i) * hl110(i))
       hl12(i) = hl1(i) * 2. / z1(i)
       ph2(i) = (a0p + a1p * hl12(i)) * hl12(i)                                &
             / (1. + b1p * hl12(i) + b2p * hl12(i) * hl12(i))
     endif
     if(dtv(i).lt.0.and.hlinf(i).lt.-.5) then
       hl1(i) = -hlinf(i)
       pm(i) = log(hl1(i)) + 2. * hl1(i) ** (-.25) - .8776
       ph(i) = log(hl1(i)) + .5 * hl1(i) ** (-.5) + 1.386
       hl110(i) = hl1(i) * 10. / z1(i)
       pm10(i) = log(hl110(i)) + 2. * hl110(i) ** (-.25) - .8776
       hl12(i) = hl1(i) * 2. / z1(i)
       ph2(i) = log(hl12(i)) + .5 * hl12(i) ** (-.5) + 1.386
     endif
   enddo
!
!  finish the exchange coefficient computation to provide fm and fh
!
   do i = 1,im
     fm(i) = fm(i) - pm(i)
     fh(i) = fh(i) - ph(i)
     fm10(i) = fm10(i) - pm10(i)
     fh2(i) = fh2(i) - ph2(i)
     cm(i) = ca * ca / (fm(i) * fm(i))
     ch(i) = ca * ca / (fm(i) * fh(i))
     cq(i) = ch(i)
#ifndef KIMHONG2010
     ustar(i) = max(sqrt(cm(i) * wind(i) * wind(i)),0.1)
#else
     ustar(i) = sqrt(cm(i) * wind(i) * wind(i))
#endif
   enddo
!
!  update z0 over ocean
!
   do i = 1,im
     if(slimsk(i).eq.0.) then
#ifndef KIMHONG2010
       z0(i) = (charnock / g) * ustar(i) ** 2
       z0(i) = min(z0(i),.1)
       z0(i) = max(z0(i),1.e-7)
#else
       z0(i) = (charnock / g) * ustar(i) ** 2 + oz0min
#endif
       z0cm(i) = 100. * z0(i)
     endif
   enddo
!      
   return
   end subroutine lsm_exch_coeff
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine lsm_diagnostic(ims2,imx2,kmx,                                    &
                      u1,v1,t1,q1,tskin,                                       &
                      slrad,rnet,                                              &
                      f10m,u10m,v10m,t2m,q2m,                                  &
                      rcl,prsi1,prsik1,prslk1,slimsk,inistp,lat,               &
                      qsurf,evap,fm,fh,fm10,fh2)
!-------------------------------------------------------------------------------
   use constant,  only : cal_,cp=>cp_,g=>g_,hfus=>hfus_,hvap=>hvap_,           &
                         rd=>rd_,rv=>rv_,sigma=>sbc_,rhoh2o=>rhoh2o_,          &
                         convrad=>convrad_,elocp=>elocp_,cb2pa=>cb2pa_,        &
                         rdorv=>rdorv_,rdorvm1=>rdorvm1_,rvordm1=>rvordm1_
   use constant,  only : qmin8_
   use paramodel, only : LONF2S,LATG2S,levs_,nstype_,nvtype_
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   real,parameter       ::  charnock=.014,ca=.4
   real,parameter       ::  alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
   real,parameter       ::  a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
   real,parameter       ::  aa1=-1.076,bb1=.7045,cc1=-.05808
   real,parameter       ::  bb2=-.1954,cc2=.009999
   real,parameter       ::  dfsnow=.31,ch2o=4.2e6,csoil=1.26e6
   real,parameter       ::  scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2
   real,parameter       ::  cice=1880.*917.,topt=298.
   real,parameter       ::  ctfil1=.5,ctfil2=1.-ctfil1
   real,parameter       ::  rnu=1.51e-5,arnu=.135*rnu
! passing array
   integer              ::  ims2,imx2,kmx
   real                 ::  prsi1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  u1(imx2),v1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),plantr(imx2),gflux(imx2)
   real                 ::  u10m(imx2),v10m(imx2),t2m(imx2),q2m(imx2)
   real                 ::  slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
   real                 ::  f10m(imx2),drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),rnet(imx2),ep(imx2)
   real                 ::  fm(imx2),fh(imx2),ustar(imx2),wind(imx2)
   real                 ::  dlwflx(imx2)
   integer              ::  vegtype(imx2)
   integer              ::  soiltyp(imx2)
!
   integer              ::  inistp,lat
   real                 ::  rcl
   real                 ::  fm10(imx2),fh2(imx2)
!
! local array
!
#if defined(RMPVECRORIZE) && defined(MP)
#define ILOTS LONF2S*LATG2S
#else
#define ILOTS LONF2S
#endif
   real                 ::  sig2k
   real                 ::  rs(ILOTS)
   real                 ::  theta1(ILOTS)
   real                 ::  tv1(ILOTS),tvs(ILOTS)
   real                 ::  z1(ILOTS),thv1(ILOTS)
   real                 ::  rho(ILOTS),qs1(ILOTS)
   real                 ::  qss(ILOTS),snowd(ILOTS)
   real                 ::  etpfac(ILOTS),tsurf(ILOTS)
   real                 ::  q0(ILOTS),cq(ILOTS)
   real                 ::  stsoil(ILOTS,levs_)
   real                 ::  dew(ILOTS)
   real                 ::  edir(ILOTS),et(ILOTS,levs_)
   real                 ::  ec(ILOTS)
   real                 ::  z0max(ILOTS),ztmax(ILOTS)
   real                 ::  dtv(ILOTS),adtv(ILOTS)
   real                 ::  hlinf(ILOTS)
   real                 ::  hl1(ILOTS),pm(ILOTS)
   real                 ::  ph(ILOTS)
   real                 ::  hl110(ILOTS),hl12(ILOTS)
   real                 ::  rcap(ILOTS),rsmall(ILOTS)
   real                 ::  pm10(ILOTS),ph2(ILOTS)
   real                 ::  olinf(ILOTS),rch(ILOTS)
   real                 ::  dft0(ILOTS)
   real                 ::  t12(ILOTS),t14(ILOTS)
   real                 ::  delta(ILOTS)
   logical              ::  flag(ILOTS)
   real                 ::  tref(ILOTS)
   real                 ::  twilt(ILOTS),df1(ILOTS)
   real                 ::  kt1(ILOTS),fx(ILOTS)
   real                 ::  gx(ILOTS),canfac(ILOTS)
   real                 ::  smcz(ILOTS),dmdz(ILOTS)
   real                 ::  ddz(ILOTS),dmdz2(ILOTS)
   real                 ::  ddz2(ILOTS),df2(ILOTS)
   real                 ::  kt2(ILOTS)
   real                 ::  xx(ILOTS),yy(ILOTS)
   real                 ::  zz(ILOTS)
   real                 ::  dtdz2(ILOTS),dft2(ILOTS)
   real                 ::  dtdz1(ILOTS),dft1(ILOTS)
   real                 ::  hcpct(ILOTS)
   real                 ::  ai(ILOTS,levs_)
   real                 ::  bi(ILOTS,levs_)
   real                 ::  ci(ILOTS,levs_)
   real                 ::  rhstc(ILOTS,levs_)
   real                 ::  factsnw(ILOTS),z0(ILOTS)
   real                 ::  slwd(ILOTS)
   logical              ::  flagsnw(ILOTS)
   real                 ::  term1(ILOTS),term2(ILOTS)
   real                 ::  partlnd(ILOTS)
   real                 ::  restar(ILOTS),rat(ILOTS)
!   real                 ::  fm10(ILOTS),fh2(ILOTS)
!
   real                 ::  snet(ILOTS),smcdry(ILOTS)
   real                 ::  rsmax(nvtype_),rgl(nvtype_)
   real                 ::  rsmin(nvtype_),hs(nvtype_)
   real                 ::  smdry(nstype_)
   real                 ::  smref(nstype_),smwlt(nstype_)
   real                 ::  ktsoil
!
   real                 ::  stcx(ILOTS,levs_)
   integer              ::  im,km,i
   real                 ::  xrcl
#ifdef ICE
   real                 ::  fpvs
   external fpvs
#else
   real                 ::  fpvs0
   external fpvs0
#endif
!-------------------------------------------------------------------------------
   im = ims2
   km = kmx
!
!     estimate sigma ** k at 2 m
!
   sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables.
!
   q0=0.
!
   xrcl = sqrt(rcl)
   do i = 1,im
     slwd(i) = slrad(i) * convrad
     theta1(i) = t1(i) / prslk1(i) * prsik1(i)
   enddo
!
! update surface layer properties
!
   do i = 1,im
     f10m(i) = fm10(i) / fm(i)
     u10m(i) = f10m(i) * xrcl * u1(i)
     v10m(i) = f10m(i) * xrcl * v1(i)
     t2m(i) = tskin(i) * (1. - fh2(i) / fh(i)) + theta1(i) * fh2(i) / fh(i)
     t2m(i) = t2m(i) * sig2k
!
!  in case of evaporation, use the inferred qsurf to deduce q2m
!
     if(evap(i).ge.0.) then
       q0(i) = max(q1(i),qmin8_)
       q2m(i) = qsurf(i) * (1. - fh2(i) / fh(i)) + q0(i) * fh2(i) / fh(i)
!
!  for dew formation situation, use saturated q at tskin
!
     else
#ifdef ICE
       qss(i) = cb2pa * fpvs(tskin(i))
#else
       qss(i) = cb2pa * fpvs0(tskin(i))
#endif
       qss(i) = rdorv * qss(i) / (cb2pa * prsi1(i) + rvordm1 * qss(i))
       q2m(i) = qss(i) * (1. - fh2(i) / fh(i)) + q0(i) * fh2(i) / fh(i)
     endif
#ifdef ICE
     qss(i) = cb2pa * fpvs(t2m(i))
#else
     qss(i) = cb2pa * fpvs0(t2m(i))
#endif
     qss(i) = rdorv * qss(i) / (cb2pa * prsi1(i) + rvordm1 * qss(i))
     q2m(i) = min(q2m(i),qss(i))
   enddo
!
!     compute net radiation
!
   do i = 1,im
     rnet(i) = -slwd(i) - sigma * tskin(i) **4
   enddo
!
   return
   end subroutine lsm_diagnostic
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine noah_vic_ocean(ims2,imx2,kmx,                                    &
                       t1,q1,tskin,qsurf,                                      &
                       dm, cm, ch,                                             &
                       prsi1,prsl1,prsik1,prslk1,                              &
                       slimsk,inistp,lat,                                      &
                       evap,hflx,wind)
!-------------------------------------------------------------------------------
   use constant, only  : cp=>cp_,g=>g_,rd=>rd_,rv=>rv_,                        &
                         hfus=>hfus_,hvap=>hvap_,sigma=>sbc_,                  &
                         cal_,rhoh2o=>rhoh2o_,cb2pa=>cb2pa_,                   &
                         convrad=>convrad_,elocp=>elocp_,                      &
                         rdorv=>rdorv_,rdorvm1=>rdorvm1_,rvordm1=>rvordm1_
   use constant,  only : qmin8_
   use paramodel, only : ILOTS,nstype_,nvtype_,levs_
!fpp$ noconcur r
!fpp$ expand(fpvs,fpvs0)
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   real,parameter       ::  charnock=.014,ca=.4
   real,parameter       ::  alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
   real,parameter       ::  a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
   real,parameter       ::  aa1=-1.076,bb1=.7045,cc1=-.05808
   real,parameter       ::  bb2=-.1954,cc2=.009999
   real,parameter       ::  dfsnow=.31,ch2o=4.2e6,csoil=1.26e6
   real,parameter       ::  scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2
   real,parameter       ::  cice=1880.*917.,topt=298.
   real,parameter       ::  ctfil1=.5,ctfil2=1.-ctfil1
   real,parameter       ::  rnu=1.51e-5,arnu=.135*rnu
!
!   passing array
!
   integer              ::  ims2,imx2,kmx
   integer              ::  inistp,lat
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  u1(imx2),v1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),plantr(imx2),gflux(imx2)
   real                 ::  u10m(imx2),v10m(imx2),t2m(imx2),q2m(imx2)
   real                 ::  slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
   real                 ::  f10m(imx2),drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),rnet(imx2),ep(imx2)
   real                 ::  fm(imx2),fh(imx2),ustar(imx2),wind(imx2)
   real                 ::  dlwflx(imx2)
   integer              ::  vegtype(imx2)
   integer              ::  soiltyp(imx2)
!
! local array
!
   real                 ::  theta1(ILOTS)
   real                 ::  tv1(ILOTS)
   real                 ::  rho(ILOTS),qs1(ILOTS)
   real                 ::  qss(ILOTS)
   real                 ::  tsurf(ILOTS)
   real                 ::  q0(ILOTS)
   real                 ::  rch(ILOTS)
   integer              ::  i,j,k,im,km
#ifdef ICE
   real                 ::  fpvs
#else
   real                 ::  fpvs0
#endif
!-------------------------------------------------------------------------------
   im = ims2
   km = kmx
!
!  initialize local variables
!
   theta1=0.            ;   tv1=0.           ;   rho=0.           ;   qs1=0.
   qss=0.               ;   tsurf=0.         ;   q0=0.            ;   rch=0. 
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  1000.*prsi1(i) is in pascals
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
!  surface
!  convert slrad to the civilized unit from langley minute-1 k-4
!  surface roughness length is converted to m from cm
!
   do i = 1,im
     if(slimsk(i).eq.0.) then
       q0(i) = max(q1(i),qmin8_)
       tsurf(i) = tskin(i)
       theta1(i) = t1(i) / prslk1(i) * prsik1(i)
       tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
       rho(i) = (cb2pa*prsl1(i)) / (rd * tv1(i))
#ifdef ICE
       qs1(i) = cb2pa * fpvs(t1(i))
#else
       qs1(i) = cb2pa * fpvs0(t1(i))
#endif
       qs1(i) = rdorv * qs1(i) / (cb2pa*prsl1(i) + rdorvm1 * qs1(i))
       qs1(i) = max(qs1(i),qmin8_)
       q0(i) = min(qs1(i),q0(i))
#ifdef ICE
       qss(i) = cb2pa * fpvs(tsurf(i))
#else
       qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
       qss(i) = rdorv * qss(i) / (cb2pa*prsi1(i) + rdorvm1 * qss(i))
     endif
   enddo
!
!  rcp = rho cp ch v
!
   do i = 1,im
     if(slimsk(i).eq.0.) then
       rch(i) = rho(i) * cp * ch(i) * wind(i)
     endif
   enddo
!
!  sensible and latent heat flux over open water
!
   do i = 1,im
     if(slimsk(i).eq.0.) then
       evap(i) = elocp * rch(i) * (qss(i) - q1(i))
       hflx(i) = rch(i) * (tskin(i) - theta1(i))
       qsurf(i) = qss(i)
       dm(i) = 1.
     endif
   enddo
!
   return
   end subroutine noah_vic_ocean
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine noah_vic_seaice(ims2,imx2,kmx,                                   &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,                    &
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,cm,ch,rcl,                                       &
                  prsi1,prsl1,prsik1,prslk1,                                   &
                  zl1,slimsk,inistp,lat,                                       &
                  drain,evap,hflx,ep,wind,                                     &
                  snowfl, runoff,precip,srflag)
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
!
! ca is the von karman constant
!
   integer              ::  ims2,imx2,kmx,lat,inistp
   real                 ::  delt,rcl
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  zl1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),plantr(imx2),gflux(imx2)
   real                 ::  slimsk(imx2)
   real                 ::  drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),ep(imx2)
   real                 ::  wind(imx2)
   real                 ::  dlwflx(imx2)
   integer              ::  vegtype(imx2)
   integer              ::  soiltyp(imx2)
   real                 ::  runoff(imx2)
   real                 ::  precip(imx2),   srflag(imx2)
   real                 ::  snowfl(imx2)
!
! local array
!
   integer              ::  i, j, k, im, km
   real                 ::  rhscnpy(imx2),rhsmc(imx2,kmx)
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
!  real                 ::  osu_funct_ktsoil
!-------------------------------------------------------------------------------
!
   im = ims2
   km = kmx
!
! initialize
!
   rhscnpy=0.           ;   rhsmc=0.
   aim=0.               ;   bim=0.           ;   cim=0.
!
! surface energy/water balance over land and seaice
!
! perform snow-rain detection
!
   do j = 1,im
     if(slimsk(j) .eq. 2.) then
!
! perform snow-rain detection using precip and srflag
!
       if(srflag(j) .eq. 1.) then
         snowfl(j)=precip(j)
         if(slimsk(j) .ne. 0.) then
           snoweq(j) = snoweq(j) + 1.e3 * precip(j)
         endif
         precip(j) = 0.
       else
         snowfl(j)=0.
       endif
     endif
   enddo
!
! compute sfc energy balance
!
   call noah_vic_seaicetm(ims2,   imx2,    kmx,                                &
                     t1,     q1,                                               &
                 snoweq,  tskin,  qsurf,                                       &
                    smc,    stc,     dm, soiltyp, sigmaf, vegtype, canopy,     &
                 dlwflx,  slrad, snowmt,  snowev,                              &
                   delt,   z0cm,    tg3,                                       &
                  gflux,                                                       &
                  zsoil, cm, ch,rhscnpy,   rhsmc,    aim, bim, cim,            &
                    rcl,                                                       &
                  prsi1,  prsl1, prsik1,  prslk1,                              &
                    zl1, slimsk, inistp,     lat,                              &
                  drain,   evap,   hflx,      ep,   wind)
!
   return
   end subroutine noah_vic_seaice
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine noah_vic_seaicetm(ims2,imx2,kmx,                                 &
                      t1,     q1,                                              &
                  snoweq,  tskin,  qsurf,                                      &
                     smc,    stc,     dm, soiltyp, sigmaf, vegtype, canopy,    &
                  dlwflx,  slrad, snowmt,  snowev,                             &
                    delt,   z0cm,    tg3,                                      &
                   gflux,                                                      &
                   zsoil, cm, ch,rhscnpy,   rhsmc,    aim, bim, cim,           &
                     rcl,                                                      &
                   prsi1,  prsl1, prsik1,  prslk1,                             &
                     zl1, slimsk, inistp,     lat,                             &
                   drain,   evap,   hflx,      ep,   wind)
!-------------------------------------------------------------------------------
   use constant, only  : cal_,cp=>cp_,g=>g_,hfus=>hfus_,hvap=>hvap_,           &
                         rd=>rd_,rv=>rv_,sig=>sbc_,t0c=>t0c_,                  &
                         rhoh2o=>rhoh2o_,convrad=>convrad_,cb2pa=>cb2pa_,      &
                         rdorv=>rdorv_,rdorvm1=>rdorvm1_,rvordm1=>rvordm1_,    &
                         elocp=>elocp_
   use constant,  only : qmin8_
#ifdef DFS
   use dfsvar, only    : iope
#else
   use comio, only     : iope
#endif
   use paramodel, only : ILOTS,levs_,nstype_,nvtype_
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   real,parameter       ::  charnock=.014,ca=.4
   real,parameter       ::  alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
   real,parameter       ::  a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
   real,parameter       ::  aa1=-1.076,bb1=.7045,cc1=-.05808
   real,parameter       ::  bb2=-.1954,cc2=.009999
   real,parameter       ::  dfsnow=.31,ch2o=4.2e6,csoil=1.26e6
!soojin
   real,parameter       ::  scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2
!   real,parameter       ::  scanop=.5,cfactr=.5,zbot=-3.,tgice=291.2
   real,parameter       ::  cice=1880.*917.,topt=298.
   real,parameter       ::  ctfil1=.5,ctfil2=1.-ctfil1
   real,parameter       ::  rnu=1.51e-5,arnu=.135*rnu
!
!  passing array
!
   integer              ::  ims2,imx2,kmx,inistp,lat
   real                 ::  delt,rcl
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  zl1(imx2),u1(imx2),v1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),plantr(imx2),gflux(imx2)
   real                 ::  u10m(imx2),v10m(imx2),t2m(imx2),q2m(imx2)
   real                 ::  slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
   real                 ::  f10m(imx2),drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),rnet(imx2),ep(imx2)
   real                 ::  fm(imx2),fh(imx2),ustar(imx2),wind(imx2)
   real                 ::  dlwflx(imx2)
   integer              ::  vegtype(imx2)
   integer              ::  soiltyp(imx2)
!
! local array
!
   integer              ::  minvegtype,maxvegtype,minsoiltyp,maxsoiltyp
   integer              ::  i,k,im,km
   real                 ::  xrcl,sig2k,tflx,bfact,cc
#ifdef ICE
   real                 ::  fpvs
   external fpvs
#else
   real                 ::  fpvs0
   external fpvs0
#endif
   real                 ::  rs(ILOTS)
   real                 ::  theta1(ILOTS)
   real                 ::  tv1(ILOTS),tvs(ILOTS)
   real                 ::  z1(ILOTS),thv1(ILOTS)
   real                 ::  rho(ILOTS),qs1(ILOTS)
   real                 ::  qss(ILOTS),snowd(ILOTS)
   real                 ::  etpfac(ILOTS),tsurf(ILOTS)
   real                 ::  q0(ILOTS),cq(ILOTS)
   real                 ::  stsoil(ILOTS,levs_),dew(ILOTS)
   real                 ::  edir(ILOTS),et(ILOTS,levs_)
   real                 ::  ec(ILOTS)
   real                 ::  z0max(ILOTS),ztmax(ILOTS)
   real                 ::  dtv(ILOTS),adtv(ILOTS)
   real                 ::  fm10(ILOTS),fh2(ILOTS)
   real                 ::  hlinf(ILOTS)
   real                 ::  hl1(ILOTS)
   real                 ::  pm(ILOTS),ph(ILOTS)
   real                 ::  hl110(ILOTS),hl12(ILOTS)
   real                 ::  rcap(ILOTS),rsmall(ILOTS)
   real                 ::  pm10(ILOTS),ph2(ILOTS)
   real                 ::  olinf(ILOTS),rch(ILOTS)
   real                 ::  dft0(ILOTS)
   real                 ::  t12(ILOTS),t14(ILOTS)
   real                 ::  delta(ILOTS)
   logical              ::  flag(ILOTS)
   real                 ::  tref(ILOTS)
   real                 ::  twilt(ILOTS),df1(ILOTS)
   real                 ::  kt1(ILOTS),fx(ILOTS)
   real                 ::  gx(ILOTS),canfac(ILOTS)
   real                 ::  smcz(ILOTS),dmdz(ILOTS)
   real                 ::  ddz(ILOTS),dmdz2(ILOTS)
   real                 ::  ddz2(ILOTS)
   real                 ::  df2(ILOTS),kt2(ILOTS)
   real                 ::  xx(ILOTS),yy(ILOTS)
   real                 ::  zz(ILOTS)
   real                 ::  dtdz2(ILOTS),dft2(ILOTS)
   real                 ::  dtdz1(ILOTS),dft1(ILOTS)
   real                 ::  hcpct(ILOTS)
   real                 ::  ai(ILOTS,levs_)
   real                 ::  bi(ILOTS,levs_)
   real                 ::  ci(ILOTS,levs_)
   real                 ::  rhstc(ILOTS,levs_)
   real                 ::  factsnw(ILOTS),z0(ILOTS)
   real                 ::  slwd(ILOTS)
   logical              ::  flagsnw(ILOTS)
   real                 ::  term1(ILOTS),term2(ILOTS)
   real                 ::  partlnd(ILOTS)
   real                 ::  restar(ILOTS),rat(ILOTS)
!
   real                 ::  snet(ILOTS),smcdry(ILOTS)
   real                 ::  rsmax(nvtype_),rgl(nvtype_)
   real                 ::  rsmin(nvtype_),hs(nvtype_)
   real                 ::  smdry(nstype_)
   real                 ::  smref(nstype_),smwlt(nstype_)
!
   real                 ::  stcx(ILOTS,levs_)
   real                 ::  osu_funct_ktsoil
!
#ifndef USGS
   do i = 1,13
     rsmax(i)=5000.
   enddo
!
   rsmin=(/150.,100.,125.,150.,100.,70.,40.,300.,400.,150.,999.,40.,999./)
   rgl=(/30.,30.,30.,30.,30.,65.,100.,100.,100.,100.,999.,100.,999./)
   hs=(/41.69,54.53,51.93,47.35,47.35,54.53,36.35,42.00,42.00,42.00,999.,36.35,999./)
#else
   do i = 1,12
     rsmax(i)=5000.
   enddo
!
   rsmin=(/40.,300.,40.,150.,999.,999.,150.,150.,400.,100.,125.,999./)
   rgl=(/100.,100.,100.,100.,999.,999.,30.,30.,100.,30.,30.,999./)
   hs=(/36.35,42.00,36.35,42.00,999.,999.,41.69,47.35,42.00,54.53,51.93,999./)
#endif
!
#ifdef STATSGO_SOIL
   smref=(/0.236, 0.283, 0.312, 0.360, 0.360, 0.329, 0.314,               &
           0.387, 0.382, 0.338, 0.404, 0.412, 0.329, 0.0, 0.108, 0.283/)
   smwlt=(/0.010, 0.028, 0.047, 0.084, 0.084, 0.066, 0.067,               &
           0.120, 0.103, 0.100, 0.126, 0.138, 0.066, 0.0, 0.006, 0.028/)
#else
!    smmax=(/.421,.464,.468,.434,.406,.465,.404,.439,.421/)
   smdry=(/.07,.14,.22,.08,.18,.16,.12,.10,.07/)
   smref=(/.283,.387,.412,.312,.338,.382,.315,.329,.283/)
   smwlt=(/.029,.119,.139,.047,.010,.103,.069,.066,.029/)
#endif
!-------------------------------------------------------------------------------
   rho = 0.  ; rch = 0. ; tsurf = 0.  
!
! initialize local variables
!
   rs=0.      ;  theta1=0.   ;  tv1=0.   ;  tvs=0.   ;  z1=0.     ;  thv1=0.
   rho=0.     ;  qs1=0.      ;  qss=0.   ;  snowd=0. ;  etpfac=0. ;  tsurf=0.
   q0=0.      ;  cq=0.       ;  stsoil=0.;  dew=0.   ;  edir=0.   ;  et=0.
   ec=0.      ;  z0max=0.    ;  ztmax=0. ;  dtv=0.   ;  adtv=0.   ;  fm10=0.
   fh2=0.     ;  hlinf=0.    ;  hl1=0.   ;  pm=0.    ;  ph=0.     ;  hl110=0.
   hl12=0.    ;  rcap=0.     ;  rsmall=0.;  pm10=0.  ;  ph2=0.    ;  olinf=0.
   rch=0.     ;  dft0=0.     ;  t12=0.   ;  t14=0.   ;  delta=0.  ;  tref=0.
   twilt=0.   ;  df1=0.      ;  kt1=0.   ;  fx=0.    ;  gx=0.     ;  canfac=0.
   smcz=0.    ;  dmdz=0.     ;  ddz=0.   ;  dmdz2=0. ;  ddz2=0.   ;  df2=0.
   kt2=0.     ;  xx=0.       ;  yy=0.    ;  zz=0.    ;  dtdz2=0.  ;  dft2=0.
   dtdz1=0.   ;  dft1=0.     ;  hcpct=0. ;  ai=0.    ;  bi=0.     ;  ci=0.
   rhstc=0.   ;  factsnw=0.  ;  z0=0.    ;  slwd=0.  ;  term1=0.  ;  term2=0.
   partlnd=0. ;  restar=0.   ;  rat=0.   ;  snet=0.  ;  smcdry=0. ;  rsmax=0.
   rgl=0.     ;  rsmin=0.    ;  hs=0.    ;  smdry=0. ;  smref=0.  ;  smwlt=0.
   stcx=0.
!
!  minvegtype and minsoiltype are defined over land only 
!  (not used over ocean)
!
   minvegtype=1
   maxvegtype=nvtype_
!
   minsoiltyp=1
   maxsoiltyp=nstype_
!
   im = ims2
   km = kmx
   xrcl = sqrt(rcl)
!
!     estimate sig ** k at 2 m
!
   sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  1000.*prsi1(i) is in pascals
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
!  surface
!  convert slrad to the civilized unit from langley minute-1 k-4
!  surface roughness length is converted to m from cm
!
   flag(1:ILOTS) = .false.
   flagsnw(1:ILOTS) = .false.
!
   do i = 1,im
     if(slimsk(i) .eq. 2.) then
       flag(i) = .true.
     endif
   enddo
!
   do i = 1,im
     if(flag(i)) then
       slwd(i) = slrad(i) * convrad
!
!  dlwflx has been given a negative sign for downward longwave
!  snet is the net shortwave flux
!
       snet(i) = -slwd(i) - dlwflx(i)
!
!   wind is now computed inside lsm_exch_coeff routine
!        wind(i) = xrcl * sqrt(u1(i) * u1(i) + v1(i) * v1(i))
!     1            +max(0.0,min(tem*dsfc(i)/(1000.*prsi1(i)),dsfcmax))
!         wind(i) = max(wind(i),1.)
!

       q0(i) = max(q1(i),qmin8_)
       tsurf(i) = tskin(i)
       theta1(i) = t1(i) / prslk1(i) * prsik1(i)
       tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
       thv1(i) = theta1(i) * (1. + rvordm1 * q0(i))
       tvs(i) = tsurf(i) * (1. + rvordm1 * q0(i))
       rho(i) = (cb2pa*prsl1(i)) / (rd * tv1(i))
#ifdef ICE
       qs1(i) = cb2pa * fpvs(t1(i))
#else
       qs1(i) = cb2pa * fpvs0(t1(i))
#endif
       qs1(i) = rdorv * qs1(i) / (cb2pa*prsl1(i) + rdorvm1 * qs1(i))
       qs1(i) = max(qs1(i), qmin8_)
       q0(i) = min(qs1(i),q0(i))
#ifdef ICE
       qss(i) = cb2pa * fpvs(tsurf(i))
#else
       qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
       qss(i) = rdorv * qss(i) / (cb2pa*prsi1(i) + rdorvm1 * qss(i))
       rs(i) = 0.
       if(vegtype(i).gt.0) rs(i) = rsmin(vegtype(i))
       z0(i) = .01 * z0cm(i)
       canopy(i)= max(canopy(i),0.)
       dm(i) = 1.
       factsnw(i) = 3.
!
!  snow depth in water equivalent is converted from mm to m unit
!
       snowd(i) = snoweq(i) / 1000.
!
!  when snow depth is less than 1 mm, a patchy snow is assumed and
!  soil is allowed to interact with the atmosphere.
!  we should eventually move to a linear combination of soil and
!  snow under the condition of patchy snow.
!
       if(snowd(i).gt..001.or.slimsk(i).eq.2) rs(i) = 0.
       if(snowd(i).gt..001) flagsnw(i) = .true.
     endif
   enddo
!
! initialization is done inside phys_noah_vic_driver routine
!
   do i = 1,im
     if(flag(i)) then
       z1(i) = zl1(i)
     endif
   enddo
!
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
          et(i,k) = 0.
          rhsmc(i,k) = 0.
          aim(i,k) = 0.
          bim(i,k) = 1.
          cim(i,k) = 0.
          stsoil(i,k) = stc(i,k)
          stcx(i,k) = stc(i,k)
       endif
     enddo
   enddo
!
   do i = 1,im
     if(flag(i)) then
       edir(i) = 0.
       ec(i) = 0.
       rhscnpy(i) = 0.
       fx(i) = 0.
       etpfac(i) = 0.
       canfac(i) = 0.
       gflux(i) = 0.
       dft0(i) = 2.2
       dew(i) = 0.
     endif
   enddo
!
!  rcp = rho cp ch v
!
   do i = 1,im
     if(flag(i)) then
       rch(i) = rho(i) * cp * ch(i) * wind(i)
     endif
   enddo
!
!  compute soil/snow/ice heat flux in preparation for surface energy
!  balance calculation
!  df for ice is taken from maykut and untersteiner
!  df is in si unit of w k-1 m-1
!  when snow covered, ground heat flux comes from snow
!
   do i = 1,im
     if(flag(i)) then
       if(flagsnw(i)) then
         tflx = min(t1(i), tsurf(i))
         gflux(i) = -dfsnow * (tflx - stsoil(i,1))                             &
                   / (factsnw(i) * max(snowd(i),.001))
       else
         gflux(i) = dft0(i) * (stsoil(i,1) - tsurf(i))                         &
                 / (-.5 * zsoil(i,1))
       endif
       gflux(i) = max(gflux(i),-200.)
       gflux(i) = min(gflux(i),+200.)
     endif
   enddo
!
   do i = 1,im
     if(flag(i)) then
       partlnd(i) = 1.
       if(snowd(i).gt.0..and.snowd(i).le..001) then
         partlnd(i) = 1. - snowd(i) / .001
       endif
     endif
   enddo
!
   do i = 1,im
     if(flag(i)) then
       snowev(i) = 0.
       if(snowd(i).gt..001) partlnd(i) = 0.
     endif
   enddo
!
!  compute potential evaporation for land and sea ice
!
   do i = 1,im
     if(flag(i)) then
       t12(i) = t1(i) * t1(i)
       t14(i) = t12(i) * t12(i)
!
!  rcap = fnet - sig t**4 + gflx - rho cp ch v (t1-theta1)
!
       rcap(i) = -slwd(i) - sig * t14(i) + gflux(i)                            &
                 - rch(i) * (t1(i) - theta1(i))
!
!  rsmall = 4 sig t**3 / rch + 1
!
       rsmall(i) = 4. * sig * t1(i) * t12(i) / rch(i) + 1.
!
!  delta = l / cp * dqs/dt
!
       delta(i) = elocp * rdorv * hvap * qs1(i) / (rd * t12(i))
!
!  potential evapotranspiration ( watts / m**2 ) and
!  potential evaporation
!
       term1(i) = elocp * rsmall(i) * rch(i)*(qs1(i)-q0(i))
       term2(i) = rcap(i) * delta(i)
       ep(i) = (elocp * rsmall(i) * rch(i) * (qs1(i) - q0(i))                  &
               + rcap(i) * delta(i))
       ep(i) = ep(i) / (rsmall(i) + delta(i))
     endif
   enddo
!
!  evaporation over bare sea ice
!
   do i = 1,im
     if(flag(i)) then
       evap(i) = partlnd(i) * ep(i)
     endif
   enddo
!
!  treat downward moisture flux situation
!  (evap was preset to zero so no update needed)
!  dew is converted from kg m-2 to m to conform to precip unit
!
   do i = 1,im
     if(flag(i).and.ep(i).le.0.) then
       dew(i) = -ep(i) * delt / (hvap * rhoh2o)
       evap(i) = ep(i)
       dew(i) = dew(i) * partlnd(i)
       evap(i) = evap(i) * partlnd(i)
       dm(i) = 1.
     endif
   enddo
!
!  snow covered land and sea ice
!  change of snow depth due to evaporation or sublimation
!  convert evap from kg m-2 s-1 to m s-1 to determine the reduction of s
!
   do i = 1,im
     if(flag(i).and.snowd(i).gt.0.) then
       bfact = snowd(i) / (delt * ep(i) / ((hvap+hfus) * rhoh2o))
       bfact = min(bfact,1.)
!
!  the evaporation of snow
!
       if(ep(i).le.0.) bfact = 1.
       if(snowd(i).le..001) then
         snowev(i) = bfact * ep(i) * (1. - partlnd(i))
         evap(i) = evap(i) + snowev(i)
       else
         snowev(i) = bfact * ep(i)
         evap(i) = snowev(i)
       endif
       tsurf(i) = t1(i)+(rcap(i)-evap(i))/(rsmall(i)*rch(i)+                   &
                 dfsnow/(factsnw(i)*max(snowd(i),.001)))
       snowd(i) = snowd(i) - snowev(i)*delt/(rhoh2o * (hvap+hfus))
       snowd(i) = max(snowd(i),0.)
     endif
   enddo
!
!  snow melt (m s-1)
!
   do i = 1,im
     if(flag(i).and.snowd(i).gt..0.and.tsurf(i).gt.t0c) then
       snowmt(i) = rch(i) * rsmall(i)                                          &
               * (tsurf(i) - t0c) / (rhoh2o * hfus)
       snowmt(i) = min(snowmt(i),snowd(i)/delt)
       snowd(i) = snowd(i) - snowmt(i) * delt
       snowd(i) = max(snowd(i),0.)
       tsurf(i) = max(t0c,tsurf(i)                                             &
              -hfus*snowmt(i)*rhoh2o/(rch(i)*rsmall(i)))
!
!  we need to re-evaluate evaporation because of snow melt
!    the skin temperature is now bounded to 0 deg c
!
#ifdef ICE
       qss(i) = cb2pa * fpvs(tsurf(i))
#else
       qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
       qss(i) = rdorv * qss(i) / (cb2pa*prsi1(i) + rdorvm1 * qss(i))
       evap(i) = elocp * rch(i) * (qss(i) - q0(i))
     endif
   enddo
!
!  update soil temperature and sea ice temperature
!  surface temperature is part of the update when snow is absent
!  heat flux from snow is explicit in time
!
   do i = 1,im
     if(flag(i)) then
       if((.not.flagsnw(i))) then
         yy(i) = t1(i) +                                                       &
                (rcap(i)-gflux(i) - evap(i)) / (rsmall(i) * rch(i))
         zz(i) = 1. + dft0(i) / (-.5 * zsoil(i,1) * rch(i) * rsmall(i))
         xx(i) = dft0(i) * (stsoil(i,1) - yy(i)) /(.5 * zsoil(i,1) * zz(i))
       else
         yy(i) = stsoil(i,1)
         zz(i) = 1.
         xx(i) = dfsnow * (stsoil(i,1) - tsurf(i))                             &
               / (-factsnw(i) * max(snowd(i),.001))
       endif
     endif
   enddo
!
!  compute the forcing and the implicit matrix elements for update
!  ch2o is the heat capacity of water and csoil is the heat capacity of
!
   do i = 1,im
     if(flag(i)) then
       smcz(i) = max(smc(i,1), smc(i,2))
       dtdz1(i) = (stsoil(i,1) - stsoil(i,2)) / (-.5 * zsoil(i,2))
       dft1(i) = dft0(i)
       hcpct(i) = cice
       dft2(i) = dft1(i)
       ddz(i) = 1. / (-.5 * zsoil(i,2))
!
!  ai, bi, and ci are the elements of the tridiagonal matrix for the
!  implicit update of the soil temperature
!
       ai(i,1) = 0.
       bi(i,1) = dft1(i) * ddz(i) / (-zsoil(i,1) * hcpct(i))
       ci(i,1) = -bi(i,1)
       bi(i,1) = bi(i,1)                                                       &
               + dft0(i) / (.5 * zsoil(i,1) **2 * hcpct(i) * zz(i))
       rhstc(i,1) = (dft1(i) * dtdz1(i) - xx(i))/ (zsoil(i,1) * hcpct(i))
     endif
   enddo
!
   do k = 2,km
     do i = 1,im
       if(flag(i)) then
         hcpct(i) = cice
       endif
     enddo
!
     if(k.lt.km) then
       do i = 1,im
         if(flag(i)) then
           dtdz2(i) = (stsoil(i,k) - stsoil(i,k+1))                            &
                  / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
           smcz(i) = max(smc(i,k), smc(i,k+1))
           ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
           ci(i,k) = -dft2(i) * ddz2(i)                                        &
               / ((zsoil(i,k-1) - zsoil(i,k)) * hcpct(i))
         endif
       enddo
     else
!
!  at the bottom, climatology is assumed at 2m depth for land and
!  freezing temperature is assumed for sea ice at z(i,km)
!
       do i = 1,im
         if(flag(i)) then
           dtdz2(i) = (stsoil(i,k) - tgice)                                    &
                  / (.5 * zsoil(i,k-1) - .5 * zsoil(i,k))
           dft2(i) = dft1(i)
           ci(i,k) = 0.
         endif
       enddo
     endif
!
     do i = 1,im
       if(flag(i)) then
         rhstc(i,k) = (dft2(i) * dtdz2(i) - dft1(i) * dtdz1(i))                &
                 / ((zsoil(i,k) - zsoil(i,k-1)) * hcpct(i))
         ai(i,k) = -dft1(i) * ddz(i)                                           &
                / ((zsoil(i,k-1) - zsoil(i,k)) * hcpct(i))
         bi(i,k) = -(ai(i,k) + ci(i,k))
         dft1(i) = dft2(i)
         dtdz1(i) = dtdz2(i)
         ddz(i) = ddz2(i)
       endif
     enddo
   enddo ! k
!
!  solve the tri-diagonal matrix
!
   do k = 1,km
     do i = 1,im
       if(flag(i))  then
         rhstc(i,k) = rhstc(i,k) * delt
         ai(i,k) = ai(i,k) * delt
         bi(i,k) = 1. + bi(i,k) * delt
         ci(i,k) = ci(i,k) * delt
       endif
     enddo
   enddo
!
!  forward elimination
!
   do i = 1,im
     if(flag(i)) then
       ci(i,1) = -ci(i,1) / bi(i,1)
       rhstc(i,1) = rhstc(i,1) / bi(i,1)
     endif
   enddo
!
   do k = 2,km
     do i = 1,im
       if(flag(i)) then
         cc = 1. / (bi(i,k) + ai(i,k) * ci(i,k-1))
         ci(i,k) = -ci(i,k) * cc
         rhstc(i,k) = (rhstc(i,k) - ai(i,k) * rhstc(i,k-1)) * cc
       endif
     enddo
   enddo
!
!  backward substituttion
!
   do i = 1,im
     if(flag(i)) then
       ci(i,km) = rhstc(i,km)
     endif
   enddo
!
   do k = km-1,1
     do i = 1,im
       if(flag(i)) then
         ci(i,k) = ci(i,k) * ci(i,k+1) + rhstc(i,k)
       endif
     enddo
   enddo
!
!  update soil and ice temperature
!
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
         stsoil(i,k) = stsoil(i,k) + ci(i,k)
       endif
     enddo
   enddo
!
!  update surface temperature for snow free surfaces
!
   do i = 1,im
     if((flag(i)).and..not.flagsnw(i)) then
       tsurf(i) = t1(i)+(rcap(i)-evap(i))/                                     &
               (rsmall(i)*rch(i)+dft0(i)/(-.5 * zsoil(i,1)))
       tsurf(i) = min(tsurf(i),t0c)
     endif
   enddo
!
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
         stsoil(i,k) = min(stsoil(i,k),t0c)
       endif
     enddo
   enddo
!
!  time filter for soil and skin temperature
!
#ifdef DBG
   if(iope) then
     call print_maxmin_seven(tsurf,im,im,1,1,1,'tsurf in seaicetm')
   endif
#endif
   if(inistp.eq.0) then
     do i = 1,im
       if(flag(i)) then
!!!!     tskin(i) = ctfil1 * tsurf(i) + ctfil2 * tskin(i)
         tskin(i) = tsurf(i)
       endif
     enddo
   endif
!
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
         stc(i,k) = stsoil(i,k)
       endif
     enddo
   enddo
!
!  gflux calculation
!
   do i = 1,im
     if(flag(i)) then
       if(flagsnw(i)) then
         gflux(i) = dfsnow * (stcx(i,1) - tskin(i))                            &
                  /(factsnw(i) * max(snowd(i),.001))
       else
         gflux(i) = dft0(i) * (stcx(i,1) - tskin(i)) / (-.5 * zsoil(i,1))
       endif
     endif
   enddo
!
!  calculate sensible heat flux
!
   do i = 1,im
     if(flag(i)) then
       hflx(i) = rch(i) * (tskin(i) - theta1(i))
     endif
   enddo
!
!  the rest of the output
!
   do i = 1,im
     if(flag(i)) then
       qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
       dm(i) = 1.
!
!  convert snow depth back to mm of water equivalent
!
       snoweq(i) = snowd(i) * 1000.
     endif
   enddo
#ifdef DBG
   if(iope) then
     call print_maxmin_seven(qsurf,im,im,1,1,1,'qsurf in seaicetm')
     call print_maxmin_seven(evap,im,im,1,1,1,'evap in seaicetm')
     call print_maxmin_seven(rch,im,im,1,1,1,'rch in seaicetm')
   endif
#endif
!
   return
   end subroutine noah_vic_seaicetm
!-------------------------------------------------------------------------------
