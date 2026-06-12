#include <define.h>
#ifndef OSU
   subroutine lsm_osu_main
   end subroutine lsm_osu_main
#else
   subroutine lsm_osu_main(ims2,imx2,kmx,                                      &
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
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [lsm_driver] --- [lsm_osu_main] * ----- [phys_lsm_osu1] *
!                                        |-- [phys_lsm_osu2] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2008-10-01  kyeong hee seol        implementation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!   2012-05-01  song-you hong          clean up
!
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS
!-------------------------------------------------------------------------------
! ca is the von karman constant
   integer              ::  i, j, k
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  zl1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),plantr(imx2),gflux(imx2)
   real                 ::  slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx)        
   real                 ::  aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
   real                 ::  drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),ep(imx2)
   real                 ::  wind(imx2)
   real                 ::  dlwflx(imx2)
#ifdef OSULSM2
   integer              ::  vegtype(imx2)
#endif
   integer              ::  soiltyp(imx2)
   real                 ::  runoff(imx2)
   real                 ::  precip(imx2),   srflag(imx2)
   real                 ::  snowfl(imx2)
#ifdef HYDRO
   real                 ::  hydrw(imx2,kmx)
#endif
   integer              ::  ims,ime,its,ite
!-------------------------------------------------------------------------------
!
! local array
!
   im = ims2
   km = kmx
!
! ---------------------------------------------------
! compute sfc energy balance
!
#ifndef OSULSM1
   call phys_lsm_osu2(kmx,t1,q1,                                               &
                  snoweq,tskin,qsurf,inistp,lat,                               &
                  smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,                    &
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,                                                 &
                  cm, ch, rhscnpy,rhsmc,aim,bim,cim,                           &
                  rcl,prsi1,prsl1,prsik1,prslk1,                               &
                  zl1,slimsk,                                                  &
                  drain,evap,hflx,ep,wind,                                     &
                  ims,ime,its,ite)
#else
   call phys_lsm_osu1(kmx,t1,q1,                                               &
                  snoweq,tskin,qsurf,                                          &
                  smc,stc,dm,soiltyp,sigmaf,canopy,                            &
                  slrad,snowmt,snowev,delt,z0cm,plantr,tg3,                    &
                  gflux,zsoil,                                                 &
                  cm, ch, rhscnpy,rhsmc,aim,bim,cim,                           &
                  rcl,prsi1,prsl1,prsik1,prslk1,                               &
                  zl1,slimsk,inistp,lat,                                       &
                  drain,evap,hflx,ep,cowave,wind,                              &
                  ims,ime,its,ite)
#endif
! ---------------------------------------------------
   return
   end subroutine lsm_osu_main
!
!-------------------------------------------------------------------------------
#include <define.h>
   subroutine phys_lsm_osu1(lsoil,t1,q1,                                       &
                         snoweq,tskin,qsurf,                                   &
                         smc,stc,dm,soiltyp,sigmaf,canopy,                     &
                         slrad,snowmt,snowev,delt,z0cm,plantr,tg3,             &
                         gflux,zsoil,                                          &
                         cm, ch, rhscnpy,rhsmc,aim,bim,cim,                    &
                         rcl,prsi1,prsl1,prsik1,prslk1,                        &
                         zl1,slimsk,inistp,lat,                                &
                         drain,evap,hflx,ep,cowave,wind,                       &
                         ims,ime,its,ite)
!-------------------------------------------------------------------------------
   use paramodel, only : nvtype_,nstype_
   use constant, only : cp=>cp_,rd=>rd_,rv=>rv_,g=>g_,rhoh2o=>rhoh2o_,         &
                        rvordm1=>rvordm1_,cal_,hvap=>hvap_,hfus=>hfus_,        &
                        sbc=>sbc_,t0c=>t0c_,akapa_,elocp=>elocp_,eps=>rdorv_,  &
                        convrad=>convrad_,epsm1=>rdorvm1_,cb2pa=>cb2pa_
#ifdef DFS
   use dfsvar, only : iope
#else
   use comio, only : iope
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include "abort.h"
!
! passing array
!
   integer              ::  lsoil,inistp,lat,                                  &
                            ims,ime,its,ite
   integer              ::  soiltyp(ims:ime)
!
   real                 ::  rcl,rdelx,delt
   real                 ::  prsi1(ims:ime),prsl1(ims:ime)
   real                 ::  prsik1(ims:ime),prslk1(ims:ime),zl1(ims:ime)
   real                 ::  u1(ims:ime),v1(ims:ime),t1(ims:ime),q1(ims:ime)
   real                 ::  snoweq(ims:ime),snowmt(ims:ime),snowev(ims:ime)
   real                 ::  cm(ims:ime),ch(ims:ime)
   real                 ::  tskin(ims:ime),qsurf(ims:ime)
   real                 ::  dm(ims:ime),slrad(ims:ime)
   real                 ::  smc(ims:ime,1:lsoil),stc(ims:ime,1:lsoil)
   real                 ::  tg3(ims:ime),canopy(ims:ime)
   real                 ::  z0cm(ims:ime),z0cmt(ims:ime)
   real                 ::  plantr(ims:ime),gflux(ims:ime)
   real                 ::  u10m(ims:ime),v10m(ims:ime)
   real                 ::  t2m(ims:ime),q2m(ims:ime)
   real                 ::  slimsk(ims:ime),rhscnpy(ims:ime)
   real                 ::  rhsmc(ims:ime,1:lsoil)
   real                 ::  aim(ims:ime,1:lsoil),bim(ims:ime,1:lsoil)
   real                 ::  cim(ims:ime,1:lsoil)
   real                 ::  f10m(ims:ime),drain(ims:ime)
   real                 ::  zsoil(ims:ime,1:lsoil),sigmaf(ims:ime)
   real                 ::  evap(ims:ime),hflx(ims:ime)
   real                 ::  ep(ims:ime)
   real                 ::  wind(ims:ime)
!
! parameters
!
   real, parameter      ::  charnock=.014
   real, parameter      ::  ca=.4         ! ca is the von karman constant
   real, parameter      ::  dfsnow=.31,ch2o=4.2e6,csoil=1.26e6
   real, parameter      ::  alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
   real, parameter      ::  a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
   real, parameter      ::  aa1=-1.076,bb1=.7045,cc1=-.05808
   real, parameter      ::  bb2=-.1954,cc2=.009999
   real, parameter      ::  scanop=2.0,cfactr=.5,zbot=-3.,tgice=271.2
   real, parameter      ::  cice=1880.*917.           
   real, parameter      ::  ctfil1=.5,ctfil2=1.-ctfil1
   real, parameter      ::  rnu=1.51e-5,arnu=.135*rnu
!
! local array
!
   integer              ::  i,k,latd,lond,delt2
   integer              ::  jx1
   real                 ::  dthv,vconv,vsgd,vcpr
   real                 ::  aa,aa0,bb,bb0,cc,ff,fms,fhs
   real                 ::  tflx,hlt,hl0,hltinf,hl0inf
   real                 ::  rsi,eth,bfact,xj1,sig2k,xrcl
   real                 ::  fpvs,fpvs0,fpvs01
   real                 ::  osu_funct_ktsoil,osu_funct_df,osu_funct_kt
   real                 ::  rs(its:ite),theta1(its:ite)
   real                 ::  tv1(its:ite),tvs(its:ite)
   real                 ::  z1(its:ite),thv1(its:ite)
   real                 ::  rho(its:ite),qs1(its:ite)
   real                 ::  qss(its:ite),snowd(its:ite)
   real                 ::  etpfac(its:ite),tsurf(its:ite)
   real                 ::  q0(its:ite)
   real                 ::  stsoil(its:ite,1:lsoil),dew(its:ite)
   real                 ::  edir(its:ite),et(its:ite,1:lsoil),ec(its:ite)
   real                 ::  rcap(its:ite),rsmall(its:ite)
   real                 ::  rch(its:ite)
   real                 ::  dft0(its:ite),t12(its:ite),t14(its:ite)
   real                 ::  delta(its:ite),tref(its:ite)
   real                 ::  twilt(its:ite),df1(its:ite),etp(its:ite)
   real                 ::  kt1(its:ite),fx(its:ite)
   real                 ::  gx(its:ite),canfac(its:ite)
   real                 ::  smcz(its:ite),dmdz(its:ite)
   real                 ::  ddz(its:ite),dmdz2(its:ite)
   real                 ::  ddz2(its:ite),df2(its:ite),kt2(its:ite)
   real                 ::  xx(its:ite),yy(its:ite),zz(its:ite)
   real                 ::  dtdz2(its:ite),dft2(its:ite)
   real                 ::  dtdz1(its:ite),dft1(its:ite),hcpct(its:ite)
   real                 ::  ai(its:ite,1:lsoil),bi(its:ite,1:lsoil)
   real                 ::  ci(its:ite,1:lsoil)
   real                 ::  rhstc(its:ite,1:lsoil)
   real                 ::  factsnw(its:ite),z0(its:ite)
   real                 ::  slwd(its:ite)
   logical              ::  flagsnw(its:ite)
   logical              ::  flag(its:ite)
   real                 ::  term1(its:ite),term2(its:ite),partlnd(its:ite)
   real                 ::  p0,ths,thvs,pa,tha,thva,tva,rhoa
   real                 ::  fluxc,cowave
   real                 ::  osu_funct_thsat,osu_funct_twlt
   real                 ::  stcx(its:ite,1:lsoil)
!
! initialize local variables
!
   rs=0.     ;  theta1=0.     ;  tv1=0.    ;  tvs=0.     ;  z1=0.
   thv1=0.   ;  rho=0.        ;  qs1=0.    ;  qss=0.     ;  snowd=0.
   etpfac=0. ;  tsurf=0.      ;  q0=0.     ;  stsoil=0.  ;  dew=0.
   edir=0.   ;  et=0.         ;  ec=0.     ;  rcap=0.    ;  rsmall=0.
   rch=0.    ;  dft0=0.       ;  t12=0.    ;  t14=0.     ;  delta=0.
   tref=0.   ;  twilt=0.      ;  df1=0.    ;  etp=0.     ;  kt1=0.
   fx=0.     ;  gx=0.         ;  canfac=0. ;  smcz=0.    ;  dmdz=0.
   ddz=0.    ;  dmdz2=0.      ;  ddz2=0.   ;  df2=0.     ;  kt2=0.
   xx=0.     ;  yy=0.         ;  zz=0.     ;  dtdz2=0.   ;  dft2=0.
   dtdz1=0.  ;  dft1=0.       ;  hcpct=0.  ;  ai=0.      ;  bi=0.
   ci=0.     ;  rhstc=0.      ;  factsnw=0.;  z0=0.      ;  slwd=0.
   term1=0.  ;  term2=0.      ;  partlnd=0.;  stcx=0.
!
!  constants
!
   latd = 23
   lond = 100
   delt2 = delt * 2.
!
!     estimate sigma ** k at 2 m
!
   sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  1000*prsi1(i) is in pascals
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
!  surface
!  convert slrad to the civilized unit from langley minute-1 k-4
!  surface roughness length is converted to m from cm
!
   xrcl = sqrt(rcl)
   do i = its,ite
     slwd(i) = slrad(i) * convrad
     q0(i) = max(q1(i),1.e-9)
     tsurf(i) = tskin(i)
     theta1(i) = t1(i) / prslk1(i) * prsik1(i)
     tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
     thv1(i) = theta1(i) * (1. + rvordm1 * q0(i))
     tvs(i) = tsurf(i) * (1. + rvordm1 * q0(i))
     rho(i) = cb2pa * prsl1(i) / (rd * tv1(i))
#ifdef ICE
     qs1(i) = cb2pa * fpvs(t1(i))
#else
     qs1(i) = cb2pa * fpvs0(t1(i))
#endif
     qs1(i) = eps * qs1(i) / (cb2pa * prsl1(i) + epsm1 * qs1(i))
     qs1(i) = max(qs1(i), 1.e-8)
#ifdef ICE
     qss(i) = cb2pa * fpvs(tsurf(i))
#else
     qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
     qss(i) = eps * qss(i) / (cb2pa * prsi1(i) + epsm1 * qss(i))
     rs(i) = plantr(i)
     z0(i) = .01 * z0cm(i)
     canopy(i)= max(canopy(i),0.)
     dm(i) = 1.
     factsnw(i) = 10.
     if(slimsk(i).eq.2.) factsnw(i) = 3.
!
!  snow depth in water equivalent is converted from mm to m unit
!
     snowd(i) = snoweq(i) / 1000.
     flagsnw(i) = .false.
!
!  when snow depth is less than 1 mm, a patchy snow is assumed and
!  soil is allowed to interact with the atmosphere.
!  we should eventually move to a linear combination of soil and
!  snow under the condition of patchy snow.
!
     if(snowd(i).gt..001.or.slimsk(i).eq.2) rs(i) = 0.
     if(snowd(i).gt..001) flagsnw(i) = .true.
   enddo
!
   do i = its,ite
     if(slimsk(i).eq.0.) then
       zsoil(i,1) = 0.
     elseif(slimsk(i).eq.1.) then
       zsoil(i,1) = -.10
     else
       zsoil(i,1) = -3. / lsoil 
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
       if(slimsk(i).eq.0.) then
         zsoil(i,k) = 0.
       elseif(slimsk(i).eq.1.) then
         zsoil(i,k) = zsoil(i,k-1) + (-2.-zsoil(i,1)) / (lsoil-1)
       else
         zsoil(i,k) = - 3. * float(k) / float(lsoil)
       endif
      enddo
   enddo
!
   do i = its,ite
     z1(i) = zl1(i)
     drain(i) = 0.
   enddo
!
   do k = 1,lsoil
     do i = its,ite
       rhsmc(i,k) = 0.
        aim(i,k) = 0.
        bim(i,k) = 1.
        cim(i,k) = 0.
        stsoil(i,k) = stc(i,k)
       stcx(i,k) = stc(i,k)
     enddo
   enddo
!
   do i = its,ite
     evap(i) = 0.
     ep(i) = 0.
     snowmt(i) = 0.
     snowev(i) = 0.
     gflux(i) = 0.
     rhscnpy(i) = 0.
   enddo
!
!  rcp = rho cp ch v
!
   do i = its,ite
     rch(i) = rho(i) * cp * ch(i) * wind(i)
   enddo
!
!  sensible and latent heat flux over open water
!
   do i = its,ite
     if(slimsk(i).eq.0.) then
       evap(i) = elocp * rch(i) * (qss(i) - q0(i))
       dm(i) = 1.
       qsurf(i) = qss(i)
     endif
   enddo
!
!  compute soil/snow/ice heat flux in preparation for surface energy
!  balance calculation
!
   do i = its,ite
     gflux(i) = 0.
     if(slimsk(i).eq.1.) then
       smcz(i) = .5 * (smc(i,1) + .20)
       dft0(i) = osu_funct_ktsoil(smcz(i),soiltyp(i))
     elseif(slimsk(i).eq.2.) then
!
!  df for ice is taken from maykut and untersteiner
!  df is in si unit of w k-1 m-1
!
       dft0(i) = 2.2
     endif
   enddo
!
   do i = its,ite
     if(slimsk(i).ne.0.) then
       if(flagsnw(i)) then
!
!  when snow covered, ground heat flux comes from snow
!
         tflx = min(t1(i), tsurf(i))
         gflux(i) = -dfsnow * (t1(i) - stsoil(i,1))                            &
                   / (factsnw(i) * max(snowd(i),.001))
       else
         gflux(i) = dft0(i) * (stsoil(i,1) - t1(i))                            &
                    / (-.5 * zsoil(i,1))
       endif
       gflux(i) = max(gflux(i),-200.)
       gflux(i) = min(gflux(i),+200.)
     endif
   enddo
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0.
     partlnd(i) = 1.
     if(snowd(i).gt.0..and.snowd(i).le..001) then
       partlnd(i) = 1. - snowd(i) / .001
     endif
   enddo
!
!  compute potential evaporation for land and sea ice
!
   do i = its,ite
     if(flag(i)) then
       t12(i) = t1(i) * t1(i)
       t14(i) = t12(i) * t12(i)
!
!  rcap = fnet - sig t**4 + gflx - rho cp ch v (t1-theta1)
!
       rcap(i) = -slwd(i) - sbc * t14(i) + gflux(i)                            &
                 - rch(i) * (t1(i) - theta1(i))
!
!  rsmall = 4 sig t**3 / rch + 1
!
       rsmall(i) = 4. * sbc * t1(i) * t12(i) / rch(i) + 1.
!
!  delta = l / cp * dqs/dt
!
       delta(i) = elocp * eps * hvap * qs1(i) / (rd * t12(i))
!
!  potential evapotranspiration ( watts / m**2 ) and
!  potential evaporation
!
       term1(i) = elocp * rsmall(i) * rch(i)*(qs1(i)-q0(i))
       term2(i) = rcap(i) * delta(i)
       ep(i) = (elocp * rsmall(i) * rch(i) * (qs1(i) - q0(i))                  &
               + rcap(i) * delta(i))
       etp(i) = ep(i) /                                                        &
                (rsmall(i) * (1. + rs(i) * wind(i) * ch(i)) + delta(i))
       ep(i) = ep(i) / (rsmall(i) + delta(i))
     endif
   enddo
!
!  actual evaporation over land in three parts : edir, et, and ec
!  direct evaporation from soil, the unit goes from m s-1 to kg m-2 s-1
!
   do i = its,ite
     flag(i) = slimsk(i).eq.1..and.ep(i).gt.0.
   enddo
!
   do i = its,ite
     if(flag(i)) then
       tref(i) = .75 * osu_funct_thsat(soiltyp(i))
       twilt(i) = osu_funct_twlt(soiltyp(i))
       df1(i) = osu_funct_df(smc(i,1),soiltyp(i))
       kt1(i) = osu_funct_kt(smc(i,1),soiltyp(i))
       fx(i) = -2. * df1(i) * (smc(i,1) - .23) / zsoil(i,1)- kt1(i)
       fx(i) = min(fx(i), ep(i)/hvap)
       fx(i) = max(fx(i),0.)
!
!  sigmaf is the fraction of area covered by vegetation
!
       edir(i) = fx(i) * (1. - sigmaf(i)) * partlnd(i)
     endif
   enddo
!
!  transpiration from all levels of the soil
!
   do i = its,ite
     if(flag(i)) then
       canfac(i) = (canopy(i) / scanop) ** cfactr
       etpfac(i) = sigmaf(i) * etp(i)* (1. - canfac(i)) / hvap
       gx(i) = (smc(i,1) - twilt(i)) / (tref(i) - twilt(i))
       gx(i) = max(gx(i),0.)
       gx(i) = min(gx(i),1.)
       et(i,1) = (zsoil(i,1) / zsoil(i,lsoil)) * gx(i) * etpfac(i)* partlnd(i)
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
       if(flag(i)) then
         gx(i) = (smc(i,k) - twilt(i)) / (tref(i) - twilt(i))
         gx(i) = max(gx(i),0.)
         gx(i) = min(gx(i),1.)
         et(i,k) =(zsoil(i,k) - zsoil(i,k-1)) / zsoil(i,lsoil)                    &
                 * gx(i) * etpfac(i) * partlnd(i)
       endif
     enddo
   enddo
!
!  canopy re-evaporation
!
   do i = its,ite
     if(flag(i)) then
       ec(i) = sigmaf(i) * canfac(i) * ep(i) / hvap
       ec(i) = ec(i) * partlnd(i)
       ec(i) = min(ec(i),canopy(i)/delt)
     endif
   enddo
!
!  sum up total evaporation
!
   do i = its,ite
     if(flag(i)) then
       evap(i) = edir(i) + ec(i)
     endif
   enddo
!
   do k = 1,lsoil
     do i = its,ite
       if(flag(i)) then
         evap(i) = evap(i) + et(i,k)
       endif
     enddo
   enddo
!
!  return evap unit from kg m-2 s-1 to watts m-2
!
   do i = its,ite
     if(flag(i)) then
       evap(i) = min(evap(i)*hvap,ep(i))
     endif
   enddo
!
!  evaporation over bare sea ice
!
   do i = its,ite
     if(slimsk(i).eq.2.) then
       evap(i) = partlnd(i) * ep(i)
     endif
   enddo
!
!  treat downward moisture flux situation
!  (evap was preset to zero so no update needed)
!  dew is converted from kg m-2 to m to conform to precip unit
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.ep(i).le.0.
   enddo
!
   do i = its,ite
     if(flag(i)) then
       dew(i) = -ep(i) * delt / (hvap * rhoh2o)
       evap(i) = ep(i)
       dm(i) = 1.
     endif
   enddo
!
!  snow covered land and sea ice
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.snowd(i).gt.0.
   enddo
!
!  change of snow depth due to evaporation or sublimation
!
!  convert evap from kg m-2 s-1 to m s-1 to determine the reduction of s
!
   do i = its,ite
     if(flag(i)) then
       bfact = snowd(i) / (delt * ep(i) / ((hvap+hfus) * rhoh2o))
       bfact = min(bfact,1.)
!
!  the evaporation of snow
!
       if(ep(i).le.0.) bfact = 1.
       if(snowd(i).le..001) then
         evap(i) = (snowd(i)/.001)*bfact*ep(i) + evap(i)
         snowev(i) = (snowd(i)/.001)*bfact*ep(i)
       else
         evap(i) = bfact * ep(i)
         snowev(i) = bfact * ep(i)
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
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.snowd(i).gt..0
   enddo
!
   do i = its,ite
     if(flag(i).and.tsurf(i).gt.t0c) then
       snowmt(i) = rch(i) * rsmall(i) * (tsurf(i) - t0c) / (rhoh2o * hfus)
       snowmt(i) = min(snowmt(i),snowd(i)/delt)
       snowd(i) = max(snowd(i) - snowmt(i) * delt, 0.)
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
       qss(i) = eps * qss(i) / (cb2pa * prsi1(i) + epsm1 * qss(i))
       evap(i) = elocp * rch(i) * (qss(i) - q0(i))
     endif
   enddo
!
!  prepare tendency terms for the soil moisture field without precipitat
!  the unit of moisture flux needs to become m s-1 for soil moisture
!   hence the factor of rhoh2o
!
   do i = its,ite
     flag(i) = slimsk(i).eq.1.
   enddo
!
   do i = its,ite
     if(flag(i)) then
       rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt
       smcz(i) = max(smc(i,1), smc(i,2))
       dmdz(i) = (smc(i,1) - smc(i,2)) / (-.5 * zsoil(i,2))
       df1(i) = osu_funct_df(smcz(i),soiltyp(i))
       kt1(i) = osu_funct_kt(smcz(i),soiltyp(i))
       rhsmc(i,1) = (df1(i) * dmdz(i) + kt1(i)                                 &
               + (edir(i) + et(i,1))) / (zsoil(i,1) * rhoh2o)
       ddz(i) = 1. / (-.5 * zsoil(i,2))
!
!  aim, bim, and cim are the elements of the tridiagonal matrix for the
!  implicit update of the soil moisture
!
       aim(i,1) = 0.
       bim(i,1) = df1(i) * ddz(i) / (-zsoil(i,1) * rhoh2o)
       cim(i,1) = -bim(i,1)
     endif
   enddo
!
   do k = 2,lsoil
     if(k.lt.lsoil) then
       do i = its,ite
         if(flag(i)) then
           dmdz2(i) = (smc(i,k) - smc(i,k+1))                                  &
                      / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
           smcz(i) = max(smc(i,k), smc(i,k+1))
           df2(i) = osu_funct_df(smcz(i),soiltyp(i))
           kt2(i) = osu_funct_kt(smcz(i),soiltyp(i))
           rhsmc(i,k) = (df2(i) * dmdz2(i) + kt2(i)                            &
                        - df1(i) * dmdz(i) - kt1(i) + et(i,k))                 &
                        / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
           ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
           cim(i,k) = -df2(i) * ddz2(i)                                        &
                   / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
         endif
        enddo
     else
       do i = its,ite
         if(flag(i)) then
           kt2(i) = osu_funct_kt(smc(i,k),soiltyp(i))
           rhsmc(i,k) = (kt2(i)                                                &
                        - df1(i) * dmdz(i) - kt1(i) + et(i,k))                 &
                        / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
           drain(i) = kt2(i)
           cim(i,k) = 0.
         endif
        enddo
     endif
     do i = its,ite
       if(flag(i)) then
         aim(i,k) = -df1(i) * ddz(i)                                           &
                   / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
         bim(i,k) = -(aim(i,k) + cim(i,k))
         df1(i) = df2(i)
         kt1(i) = kt2(i)
         dmdz(i) = dmdz2(i)
         ddz(i) = ddz2(i)
       endif
     enddo
   enddo
!
!  update soil temperature and sea ice temperature
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0.
   enddo
!
!  surface temperature is part of the update when snow is absent
!
   do i = its,ite
     if(flag(i).and..not.flagsnw(i)) then
       yy(i) = t1(i) +                                                         &
             (rcap(i)-gflux(i) - evap(i)) / (rsmall(i) * rch(i))
       zz(i) = 1. + dft0(i) / (-.5 * zsoil(i,1) * rch(i) * rsmall(i))
       xx(i) = dft0(i) * (stsoil(i,1) - yy(i)) / (.5 * zsoil(i,1) * zz(i))
     endif
     if(flag(i).and.flagsnw(i)) then
       yy(i) = stsoil(i,1)
!
!  heat flux from snow is explicit in time
!
       zz(i) = 1.
       xx(i) = dfsnow * (stsoil(i,1) - tsurf(i))                               &
               / (-factsnw(i) * max(snowd(i),.001))
     endif
   enddo
!
!  compute the forcing and the implicit matrix elements for update
!
!  ch2o is the heat capacity of water and csoil is the heat capacity of
!
   do i = its,ite
     if(flag(i)) then
       smcz(i) = max(smc(i,1), smc(i,2))
       dtdz1(i) = (stsoil(i,1) - stsoil(i,2)) / (-.5 * zsoil(i,2))
       if(slimsk(i).eq.1.) then
         dft1(i) = osu_funct_ktsoil(smcz(i),soiltyp(i))
         hcpct(i) = smc(i,1) * ch2o + (1. - smc(i,1)) * csoil
       else
         dft1(i) = dft0(i)
         hcpct(i) = cice
       endif
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
       rhstc(i,1) = (dft1(i) * dtdz1(i) - xx(i))                               &
                    / (zsoil(i,1) * hcpct(i))
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
       if(slimsk(i).eq.1.) then
         hcpct(i) = smc(i,k) * ch2o + (1. - smc(i,k)) * csoil
       elseif(slimsk(i).eq.2.) then
         hcpct(i) = cice
       endif
     enddo
     if(k.lt.lsoil) then
       do i = its,ite
         if(flag(i)) then
           dtdz2(i) = (stsoil(i,k) - stsoil(i,k+1))                            &
                      / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
           smcz(i) = max(smc(i,k), smc(i,k+1))
           if(slimsk(i).eq.1.) then
             dft2(i) = osu_funct_ktsoil(smcz(i),soiltyp(i))
           endif
           ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
           ci(i,k) = -dft2(i) * ddz2(i)                                        &
                     / ((zsoil(i,k-1) - zsoil(i,k)) * hcpct(i))
         endif
       enddo
     else
!
!  at the bottom, climatology is assumed at 2m depth for land and
!  freezing temperature is assumed for sea ice at z(i,lsoil)
!
       do i = its,ite
         if(slimsk(i).eq.1.) then
           dtdz2(i) = (stsoil(i,k) - tg3(i))                                   &
                      / (.5 * (zsoil(i,k-1) + zsoil(i,k)) - zbot)
           dft2(i) = osu_funct_ktsoil(smc(i,k),soiltyp(i))
           ci(i,k) = 0.
         endif
         if(slimsk(i).eq.2.) then
           dtdz2(i) = (stsoil(i,k) - tgice)                                    &
                      / (.5 * zsoil(i,k-1) - .5 * zsoil(i,k))
           dft2(i) = dft1(i)
           ci(i,k) = 0.
         endif
       enddo
     endif
     do i = its,ite
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
   enddo
!
!  solve the tri-diagonal matrix
!
   do k = 1,lsoil
     do i = its,ite
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
   do i = its,ite
     if(flag(i)) then
       ci(i,1) = -ci(i,1) / bi(i,1)
       rhstc(i,1) = rhstc(i,1) / bi(i,1)
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
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
   do i = its,ite
     if(flag(i)) then
       ci(i,lsoil) = rhstc(i,lsoil)
     endif
   enddo
   do k = lsoil-1,1
     do i = its,ite
       if(flag(i)) then
         ci(i,k) = ci(i,k) * ci(i,k+1) + rhstc(i,k)
       endif
     enddo
   enddo
!
!  update soil and ice temperature
!
   do k = 1,lsoil
     do i = its,ite
       if(flag(i)) then
         stsoil(i,k) = stsoil(i,k) + ci(i,k)
       endif
     enddo
   enddo
!
!  update surface temperature for snow free surfaces
!
   do i = its,ite
     if(slimsk(i).ne.0..and..not.flagsnw(i)) then
       tsurf(i) = t1(i)+(rcap(i)-evap(i))/                                     &
                 (rsmall(i)*rch(i)+dft0(i)/(-.5 * zsoil(i,1)))
     endif
     if(slimsk(i).eq.2..and..not.flagsnw(i)) then
       tsurf(i) = min(tsurf(i),t0c)
     endif
   enddo
   do k = 1,lsoil
     do i = its,ite
       if(slimsk(i).eq.2) then
         stc(i,k) = min(stsoil(i,k),t0c)
       endif
     enddo
   enddo
!
!  time filter for soil and skin temperature
!
#ifdef DBG
#ifdef MP
   if(iope) then
#endif
     call print_maxmin_six(tsurf,ite,1,1,1,'tsurf in phys_lsm_osu1')
#ifdef MP
   endif
#endif
#endif
   if(inistp.eq.0) then
     do i = its,ite
       if(slimsk(i).ne.0.) then
         tskin(i) = tsurf(i)
       endif
     enddo
#ifdef DBG
#ifdef MP
     if(iope) then
#endif
       call print_maxmin_six(tskin,ite,1,1,1,'tskin in phys_lsm_osu1')
#ifdef MP
     endif
#endif
#endif
     do k = 1,lsoil
       do i = its,ite
         if(slimsk(i).ne.0.) then
           stc(i,k) = stsoil(i,k)
         endif
       enddo
     enddo
   endif
!
!  gflux calculation
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0.                                                 &
               .and.flagsnw(i)
   enddo
!
   do i = its,ite
     if(flag(i)) then
       gflux(i) = dfsnow * (stcx(i,1) - tskin(i))                              &
                /(factsnw(i) * max(snowd(i),.001))
     endif
   enddo
!
   do i = its,ite
     if(slimsk(i).ne.0..and..not.flagsnw(i)) then
       gflux(i) = dft0(i) * (stcx(i,1) - tskin(i))/ (-.5 * zsoil(i,1))
     endif
   enddo
!
!  calculate sensible heat flux
!
   do i = its,ite
     hflx(i) = rch(i) * (tskin(i) - theta1(i))
   enddo
!
!  the rest of the output
!
   do i = its,ite
     qsurf(i) = q0(i) + evap(i) / (elocp * rch(i))
     dm(i) = 1.
!
!  convert snow depth back to mm of water equivalent
!
     snoweq(i) = snowd(i) * 1000.
   enddo
#ifdef DBG
   if(iope) then
     call print_maxmin_six(qsurf,ite,1,1,1,'qsurf in phys_lsm_osu1')
     call print_maxmin_six(evap,ite,1,1,1,'evap in phys_lsm_osu1')
     call print_maxmin_six(rch,ite,1,1,1,'rch in phys_lsm_osu1')
   endif
#endif
!
   return
   end subroutine phys_lsm_osu1
!
!-------------------------------------------------------------------------------
#include <define.h>
   subroutine phys_lsm_osu2(lsoil,t1,q1,                                       &
                         snoweq,tskin,qsurf,inistp,lat,                        &
                         smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,             &
                         dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,             &
                         gflux,zsoil,                                          &
                         cm, ch, rhscnpy,rhsmc,aim,bim,cim,                    &
                         rcl,prsi1,prsl1,prsik1,prslk1,                        &
                         zl1,slimsk,                                           &
                         drain,evap,hflx,ep,wind,                              &
                         ims,ime, its,ite)
!-------------------------------------------------------------------------------
   use paramodel, only : nvtype_,nstype_
   use constant, only : cp=>cp_,rd=>rd_,rv=>rv_,g=>g_,rhoh2o=>rhoh2o_,         &
                        rvordm1=>rvordm1_,cal_,hvap=>hvap_,hfus=>hfus_,        &
                        sbc=>sbc_,t0c=>t0c_,akapa_,elocp=>elocp_,eps=>rdorv_,  &
                        convrad=>convrad_,epsm1=>rdorvm1_,cb2pa=>cb2pa_
#ifdef DFS
   use dfsvar, only : iope
#else
   use comio, only : iope
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include "abort.h"
!
! passing array
!
   integer              ::  lsoil,inistp,lat,                                  &
                            ims,ime, its,ite
   integer              ::  vegtype(ims:ime)
   integer              ::  soiltyp(ims:ime)
!
   real                 ::  rcl,rdelx,delt
   real                 ::  prsi1(ims:ime),prsl1(ims:ime)
   real                 ::  prsik1(ims:ime),prslk1(ims:ime),zl1(ims:ime)
   real                 ::  t1(ims:ime),q1(ims:ime)
   real                 ::  snoweq(ims:ime),snowmt(ims:ime),snowev(ims:ime)
   real                 ::  cm(ims:ime),ch(ims:ime)
   real                 ::  tskin(ims:ime),qsurf(ims:ime)
   real                 ::  dm(ims:ime),slrad(ims:ime)
   real                 ::  smc(ims:ime,1:lsoil),stc(ims:ime,1:lsoil)
   real                 ::  tg3(ims:ime),canopy(ims:ime)
   real                 ::  z0cm(ims:ime),plantr(ims:ime),gflux(ims:ime)                
   real                 ::  slimsk(ims:ime),rhscnpy(ims:ime)
   real                 ::  rhsmc(ims:ime,1:lsoil)
   real                 ::  aim(ims:ime,1:lsoil),bim(ims:ime,1:lsoil)
   real                 ::  cim(ims:ime,1:lsoil)
   real                 ::  drain(ims:ime)
   real                 ::  zsoil(ims:ime,1:lsoil),sigmaf(ims:ime)
   real                 ::  evap(ims:ime),hflx(ims:ime)
   real                 ::  ep(ims:ime)
   real                 ::  ustar(ims:ime),wind(ims:ime)
   real                 ::  dlwflx(ims:ime)
!
! parameters
!
   real, parameter      ::  charnock=.014
   real, parameter      ::  ca=.4         ! ca is the von karman constant
   real, parameter      ::  dfsnow=.31,ch2o=4.2e6,csoil=1.26e6
   real, parameter      ::  alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
   real, parameter      ::  a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
   real, parameter      ::  aa1=-1.076,bb1=.7045,cc1=-.05808
   real, parameter      ::  bb2=-.1954,cc2=.009999
#ifndef HYDRO
   real, parameter      ::  scanop=2.0,cfactr=.5,zbot=-3.,tgice=271.2
#else
   real, parameter      ::  scanop=0.5,cfactr=.5,zbot=-3.,tgice=271.2
#endif
!  real, parameter      ::  scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2
   real, parameter      ::  cice=1880.*917.,topt=298.
   real, parameter      ::  ctfil1=.5,ctfil2=1.-ctfil1
   real, parameter      ::  rnu=1.51e-5,arnu=.135*rnu
#ifndef ICE
   integer, parameter   ::  nx=7501
   real                 ::  c1xpvs0,c2xpvs0
   real                 ::  tbpvs0(nx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
!
! local array
!
   integer              ::  i,k,latd,lond,delt2
   integer              ::  minvegtype,maxvegtype
   integer              ::  minsoiltyp,maxsoiltyp
   integer              ::  jx1
   real                 ::  dthv,vconv,vsgd,vcpr
   real                 ::  aa,aa0,bb,bb0,cc,ff,fms,fhs
   real                 ::  tflx,hlt,hl0,hltinf,hl0inf,rcs,rct,rcq,rss
   real                 ::  rsi,eth,bfact,xj1,sig2k,xrcl
   real                 ::  fpvs,fpvs0,fpvs01,osu_funct_df,osu_funct_kt
   real                 ::  osu_funct_ktsoil 
   real                 ::  rs(its:ite),theta1(its:ite)
   real                 ::  tv1(its:ite),tvs(its:ite)
   real                 ::  z1(its:ite),thv1(its:ite)
   real                 ::  rho(its:ite),qs1(its:ite)
   real                 ::  qss(its:ite),snowd(its:ite)
   real                 ::  etpfac(its:ite),tsurf(its:ite)
   real                 ::  q0(its:ite)
   real                 ::  stsoil(its:ite,1:lsoil),dew(its:ite)
   real                 ::  edir(its:ite),et(its:ite,1:lsoil),ec(its:ite)
   real                 ::  rcap(its:ite),rsmall(its:ite)
   real                 ::  rch(its:ite)
   real                 ::  dft0(its:ite),t12(its:ite),t14(its:ite)
   real                 ::  delta(its:ite),tref(its:ite)
   real                 ::  twilt(its:ite),df1(its:ite)
   real                 ::  kt1(its:ite),fx(its:ite)
   real                 ::  gx(its:ite),canfac(its:ite)
   real                 ::  smcz(its:ite),dmdz(its:ite)
   real                 ::  ddz(its:ite),dmdz2(its:ite)
   real                 ::  ddz2(its:ite),df2(its:ite),kt2(its:ite)
   real                 ::  xx(its:ite),yy(its:ite),zz(its:ite)
   real                 ::  dtdz2(its:ite),dft2(its:ite)
   real                 ::  dtdz1(its:ite),dft1(its:ite),hcpct(its:ite)
   real                 ::  ai(its:ite,1:lsoil),bi(its:ite,1:lsoil)
   real                 ::  ci(its:ite,1:lsoil)
   real                 ::  rhstc(its:ite,1:lsoil)
   real                 ::  factsnw(its:ite),z0(its:ite)
   real                 ::  slwd(its:ite)
   logical              ::  flagsnw(its:ite)
   logical              ::  flag(its:ite)
   real                 ::  term1(its:ite),term2(its:ite),partlnd(its:ite)
!
   real                 ::  snet(its:ite),smcdry(its:ite)
   real                 ::  rsmax(1:nvtype_),rgl(1:nvtype_)
   real                 ::  rsmin(1:nvtype_),hs(1:nvtype_)
   real                 ::  smdry(1:nstype_)
   real                 ::  smref(1:nstype_),smwlt(1:nstype_)
!
   real                 ::  stcx(its:ite,1:lsoil)
#ifdef HYDRO
   real                 ::  hydrow(its:ite,1:lsoil)
   real                 ::  ecx(its:ite)
#endif
!
! initialize local variables
!
   rs=0.     ;  theta1=0.    ;  tv1=0.    ;  tvs=0.    ;  z1=0.
   thv1=0.   ;  rho=0.       ;  qs1=0.    ;  qss=0.    ;  snowd=0.
   etpfac=0. ;  tsurf=0.     ;  q0=0.     ;  stsoil=0. ;  dew=0.
   edir=0.   ;  et=0.        ;  ec=0.     ;  rcap=0.   ;  rsmall=0.
   rch=0.    ;  dft0=0.      ;  t12=0.    ;  t14=0.    ;  delta=0.
   tref=0.   ;  twilt=0.     ;  df1=0.    ;  kt1=0.    ;  fx=0.
   gx=0.     ;  canfac=0.    ;  smcz=0.   ;  dmdz=0.   ;  ddz=0.
   dmdz2=0.  ;  ddz2=0.      ;  df2=0.    ;  kt2=0.    ;  xx=0.
   yy=0.     ;  zz=0.        ;  dtdz2=0.  ;  dft2=0.   ;  dtdz1=0.
   dft1=0.   ;  hcpct=0.     ;  ai=0.     ;  bi=0.     ;  ci=0.
   rhstc=0.  ;  factsnw=0.   ;  z0=0.     ;  slwd=0.   ;  term1=0.
   term2=0.  ;  partlnd=0.   ;  snet=0.   ;  smcdry=0. ;  rsmax=0.
   rgl=0.    ;  rsmin=0.     ;  hs=0.     ;  smdry=0.  ;  smref=0.
   smwlt=0.  ;  stcx=0.
#ifdef HYDRO
   hydrow=0. ;  ecx=0.
#endif
!
!  constants
!
!
!  the 13 vegetation types are:
!
!  1  ...  broadleave-evergreen trees (tropical forest)
!  2  ...  broadleave-deciduous trees
!  3  ...  broadleave and needle leave trees (mixed forest)
!  4  ...  needleleave-evergreen trees
!  5  ...  needleleave-deciduous trees (larch)
!  6  ...  broadleave trees with groundcover (savanna)
!  7  ...  groundcover only (perenial)
!  8  ...  broadleave shrubs with perenial groundcover
!  9  ...  broadleave shrubs with bare soil
! 10  ...  dwarf trees and shrubs with ground cover (trunda)
! 11  ...  bare soil
! 12  ...  cultivations (use parameters from type 7)
! 13  ...  glacial
!
!  the 12 vegetation type of USGS data are:
!
!  1  ... mixed farming tall grassland  (assumed to be type 12 of old veg)
!  2  ... tall/medium grassland, shrubland            (type  8 of old veg)
!  3  ... short grassland meadow and shrubland        (type  7 of old veg)
!  4  ... tundra                                      (type 10 of old veg)
!  5  ... sandy dessert                               (type 11 of old veg)
!  6  ... rocky dessert                               (type 11 of old veg)
!  7  ... tropical evergreen broadleaved forest       (type  1 of old veg)
!  8  ... evergreen forest, needleleaved forest       (type  4 of old veg)
!  9  ... medium grassland, woodland                  (type  9 of old veg)
! 10  ... deciduous forest                            (type  2 of old veg)
! 11  ... mixed deciduous and evergreen forest        (type  3 of old veg)
! 12  ... ice                                         (type 13 of old veg)
!
#ifndef USGS
   do i = 1,13
     rsmax(i)=5000.
   enddo
   rsmin=(/150.,100.,125.,150.,100.,70.,40.,300.,400.,150.,999.,40.,999./)
   rgl=(/30.,30.,30.,30.,30.,65.,100.,100.,100.,100.,999.,100.,999./)
   hs=(/41.69,54.53,51.93,47.35,47.35,54.53,36.35,42.00,42.00,42.00,999.,36.35,999./)
#else
   do i = 1,12
     rsmax(i)=5000.
   enddo
   rsmin=(/40.,300.,40.,150.,999.,999.,150.,150.,400.,100.,125.,999./)
   rgl=(/100.,100.,100.,100.,999.,999.,30.,30.,100.,30.,30.,999./)
   hs=(/36.35,42.00,36.35,42.00,999.,999.,41.69,47.35,42.00,54.53,51.93,999./)
#endif
!
#ifdef STATSGO_SOIL
   smref=(/0.236, 0.283, 0.312, 0.360, 0.360, 0.329, 0.314,                    &
         0.387, 0.382, 0.338, 0.404, 0.412, 0.329, 0.0, 0.108, 0.283/)
   smwlt=(/0.010, 0.028, 0.047, 0.084, 0.084, 0.066, 0.067,                    &
         0.120, 0.103, 0.100, 0.126, 0.138, 0.066, 0.0, 0.006, 0.028/)
#else
!
!  smmax=(/.421,.464,.468,.434,.406,.465,.404,.439,.421/)
!
   smdry=(/.07,.14,.22,.08,.18,.16,.12,.10,.07/)
   smref=(/.283,.387,.412,.312,.338,.382,.315,.329,.283/)
   smwlt=(/.029,.119,.139,.047,.010,.103,.069,.066,.029/)
#endif
!-------------------------------------------------------------------------------
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
   delt2 = delt * 2.
!
!     estimate sigma ** k at 2 m
!
   sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  1000*prsi1(i) is in pascals
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
!  surface
!  convert slrad to the civilized unit from langley minute-1 k-4
!  surface roughness length is converted to m from cm
!
   xrcl = sqrt(rcl)
   do i = its,ite
     slwd(i) = slrad(i) * convrad
!
!  dlwflx has been given a negative sign for downward longwave
!  snet is the net shortwave flux
!
     snet(i) = -slwd(i) - dlwflx(i)
     q0(i) = max(q1(i),1.e-8)
     tsurf(i) = tskin(i)
     theta1(i) = t1(i) / prslk1(i) * prsik1(i)
     tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
     thv1(i) = theta1(i) * (1. + rvordm1 * q0(i))
     tvs(i) = tsurf(i) * (1. + rvordm1 * q0(i))
     rho(i) = cb2pa * prsl1(i) / (rd * tv1(i))
#ifdef ICE
     qs1(i) = cb2pa * fpvs(t1(i))
#else
     qs1(i) = cb2pa * fpvs0(t1(i))
#endif
     qs1(i) = eps * qs1(i) / (cb2pa * prsl1(i) + epsm1 * qs1(i))
     qs1(i) = max(qs1(i), 1.e-8)
     q0(i) = min(qs1(i),q0(i))
#ifdef ICE
     qss(i) = cb2pa * fpvs(tsurf(i))
#else
     qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
     qss(i) = eps * qss(i) / (cb2pa * prsi1(i) + epsm1 * qss(i))
     if(slimsk(i).eq.1) then
       if(vegtype(i).lt.minvegtype.or.vegtype(i).gt.maxvegtype) then
         if (iope) then
           write(6,*) 'illegal vegetation type'
           write(6,*) 'vegtype(i)=',vegtype(i)
         endif
         call MPABORT
       endif
     endif
     if(vegtype(i).gt.0) rs(i) = rsmin(vegtype(i))
     z0(i) = .01 * z0cm(i)
#ifndef HYDRO
     canopy(i)= max(canopy(i),0.)
#endif
     dm(i) = 1.
     factsnw(i) = 10.
     if(slimsk(i).eq.2.) factsnw(i) = 3.
!
!  snow depth in water equivalent is converted from mm to m unit
!
     snowd(i) = snoweq(i) / 1000.
     flagsnw(i) = .false.
!
!  when snow depth is less than 1 mm, a patchy snow is assumed and
!  soil is allowed to interact with the atmosphere.
!  we should eventually move to a linear combination of soil and
!  snow under the condition of patchy snow.
!
     if(snowd(i).gt..001.or.slimsk(i).eq.2) rs(i) = 0.
     if(snowd(i).gt..001) flagsnw(i) = .true.
   enddo
!----------
   do i = its,ite
     z1(i) = zl1(i)
   enddo
!
   do k = 1,lsoil
     do i = its,ite
       rhsmc(i,k) = 0.
       aim(i,k) = 0.
       bim(i,k) = 1.
       cim(i,k) = 0.
       stsoil(i,k) = stc(i,k)
       stcx(i,k) = stc(i,k)
     enddo
   enddo
!
   do i = its,ite
     rhscnpy(i) = 0.
   enddo
#ifdef HYDRO
!
   do i = its,ite
     if(canopy(i).lt.0.) then
       ecx(i)=canopy(i)/delt
       canopy(i)= 0.
     endif
   enddo
#endif
!
!  rcp = rho cp ch v
!
   do i = its,ite
     rch(i) = rho(i) * cp * ch(i) * wind(i)
   enddo
!
!  sensible and latent heat flux over open water
!
   do i = its,ite
     if(slimsk(i).eq.0.) then
       evap(i) = elocp * rch(i) * (qss(i) - q1(i))
       dm(i) = 1.
       qsurf(i) = qss(i)
     endif
   enddo
!
!  compute soil/snow/ice heat flux in preparation for surface energy
!  balance calculation
!
   do i = its,ite
     gflux(i) = 0.
     if(slimsk(i).eq.1.) then
       smcz(i) = .5 * (smc(i,1) + .20)
       if(soiltyp(i).lt.minsoiltyp.or.soiltyp(i).gt.maxsoiltyp) then
         if (iope) then
           write(6,*) 'illegal soil type'
           write(6,*) 'soiltyp(i)=',soiltyp(i)
         endif
         call MPABORT
       endif
       dft0(i) = osu_funct_ktsoil(smcz(i),soiltyp(i))
     elseif(slimsk(i).eq.2.) then
!
!  df for ice is taken from maykut and untersteiner
!  df is in si unit of w k-1 m-1
!
       dft0(i) = 2.2
     endif
   enddo
!
 300  continue
!
   do i = its,ite
     if(slimsk(i).ne.0.) then
       if(flagsnw(i)) then
!
!  when snow covered, ground heat flux comes from snow
!
         tflx = min(t1(i), tsurf(i))
         gflux(i) = -dfsnow * (tflx - stsoil(i,1))                             &
                    / (factsnw(i) * max(snowd(i),.001))
       else
          gflux(i) = dft0(i) * (stsoil(i,1) - t1(i))                           &
                    / (-.5 * zsoil(i,1))
       endif
       gflux(i) = max(gflux(i),-200.)
       gflux(i) = min(gflux(i),+200.)
     endif
   enddo
   do i = its,ite
     flag(i) = slimsk(i).ne.0.
     partlnd(i) = 1.
     if(snowd(i).gt.0..and.snowd(i).le..001) then
       partlnd(i) = 1. - snowd(i) / .001
     endif
   enddo
   do i = its,ite
     snowev(i) = 0.
     if(snowd(i).gt..001) partlnd(i) = 0.
   enddo
!
!  compute potential evaporation for land and sea ice
!
   do i = its,ite
     if(flag(i)) then
       t12(i) = t1(i) * t1(i)
       t14(i) = t12(i) * t12(i)
!
!  rcap = fnet - sig t**4 + gflx - rho cp ch v (t1-theta1)
!
       rcap(i) = -slwd(i) - sbc * t14(i) + gflux(i)                            &
                 - rch(i) * (t1(i) - theta1(i))
!
!  rsmall = 4 sig t**3 / rch + 1
!
       rsmall(i) = 4. * sbc * t1(i) * t12(i) / rch(i) + 1.
!
!  delta = l / cp * dqs/dt
!
       delta(i) = elocp * eps * hvap * qs1(i) / (rd * t12(i))
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
!  actual evaporation over land in three parts : edir, et, and ec
!  direct evaporation from soil, the unit goes from m s-1 to kg m-2 s-1
!
   do i = its,ite
     flag(i) = slimsk(i).eq.1..and.ep(i).gt.0.
   enddo
!
   do i = its,ite
     if(flag(i)) then
!
! soiltype???
!
       df1(i) = osu_funct_df(smc(i,1),soiltyp(i))
       kt1(i) = osu_funct_kt(smc(i,1),soiltyp(i))
     endif
     if(flag(i).and.stc(i,1).lt.t0c) then
       df1(i) = 0.
       kt1(i) = 0.
     endif
     if(flag(i)) then
       tref(i) = smref(soiltyp(i))
       twilt(i) = smwlt(soiltyp(i))
#ifndef STATSGO_SOIL
       smcdry(i) = smdry(soiltyp(i))
#endif
!
!  sigmaf is the fraction of area covered by vegetation
!
#ifdef NCAR_EDIR
       fx(i) = min(max((smc(i,1) - twilt(i))                                   &
           / (tref(i) - twilt(i)),0.),1.)
       edir(i) = fx(i) * (1. - sigmaf(i)) * partlnd(i) * ep(i) / hvap
#else
       fx(i) = -2. * df1(i) * (smc(i,1) - smcdry(i)) / zsoil(i,1)- kt1(i)
       fx(i) = min(fx(i), ep(i)/hvap)
       fx(i) = max(fx(i),0.)
       edir(i) = fx(i) * (1. - sigmaf(i)) * partlnd(i)
#ifdef DBG
       if(lat.eq.latd) then
         if (iope) write(6,*)' i j edir smc sigmaf partland ',                 &
                  i,lat,edir(i), smc(i,1),sigmaf(i),partlnd(i)
       endif
#endif
#endif
     endif
   enddo
!
!  calculate stomatal resistance
!
   do i = its,ite
     if(flag(i)) then
!
!  resistance due to par. we use net solar flux as proxy at the present
!
       if(vegtype(i).lt.minvegtype.or.vegtype(i).gt.maxvegtype) then
         if (iope) then
               write(6,*) 'illegal vegetation type for use in rgl'
               write(6,*) 'vegtype(i)=',vegtype(i)
         endif
         call MPABORT
       endif
       ff = .55 * 2. * snet(i) / rgl(vegtype(i))
       rcs = (ff + rs(i)/rsmax(vegtype(i))) / (1. + ff)
       rcs = max(rcs,.0001)
       rct = 1.
       rcq = 1.
!
!  compute resistance without the effect of soil moisture
!
       rs(i) = rs(i) / (rcs * rct * rcq)
     endif
   enddo
!
!  transpiration from all levels of the soil
!
   do i = its,ite
     if(flag(i)) then
       canfac(i) = (canopy(i) / scanop) ** cfactr
       etpfac(i) = sigmaf(i)* (1. - canfac(i)) / hvap
       gx(i) = (smc(i,1) - twilt(i)) / (tref(i) - twilt(i))
       gx(i) = max(gx(i),0.)
       gx(i) = min(gx(i),1.)
!
!  resistance due to soil moisture deficit
!
       rss = gx(i) * (zsoil(i,1) / zsoil(i,lsoil))
       rss = max(rss,.0001)
       rsi = rs(i) / rss
!
!  transpiration a la monteith
!
       eth = (term1(i) + term2(i)) /                                           &
             (delta(i) + rsmall(i) * (1. + rsi * ch(i) * wind(i)))
       et(i,1) = etpfac(i) * eth* partlnd(i)
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
       if(flag(i)) then
         gx(i) = (smc(i,k) - twilt(i)) / (tref(i) - twilt(i))
         gx(i) = max(gx(i),0.)
         gx(i) = min(gx(i),1.)
!
!  resistance due to soil moisture deficit
!
         rss = gx(i) * ((zsoil(i,k) - zsoil(i,k-1)) / zsoil(i,lsoil))
         rss = max(rss,1.e-6)
         rsi = rs(i) / rss
!
!  transpiration a la monteith
!
         eth = (term1(i) + term2(i)) /                                         &
               (delta(i) + rsmall(i) * (1. + rsi * ch(i) * wind(i)))
         et(i,k) = eth* etpfac(i) * partlnd(i)
       endif
     enddo
   enddo
!
 400  continue
!
!  canopy re-evaporation
!
   do i = its,ite
     if(flag(i)) then
       ec(i) = sigmaf(i) * canfac(i) * ep(i) / hvap
       ec(i) = ec(i) * partlnd(i)
       ec(i) = min(ec(i),canopy(i)/delt)
#ifdef HYDRO
       canopy(i)=canopy(i)-ec(i)*delt
       ec(i)=ec(i)+ecx(i)
#endif
     endif
   enddo
!
!  sum up total evaporation
!
   do i = its,ite
     if(flag(i)) then
       evap(i) = edir(i) + ec(i)
#ifdef HYDRO
       ec(i)=ec(i)*hvap
#endif
     endif
   enddo
!
   do k = 1,lsoil
     do i = its,ite
       if(flag(i)) then
         evap(i) = evap(i) + et(i,k)
       endif
     enddo
   enddo
!
!  return evap unit from kg m-2 s-1 to watts m-2
!
   do i = its,ite
     if(flag(i)) then
#ifndef HYDRO
       evap(i) = min(evap(i)*hvap,ep(i))
#else
       evap(i) = evap(i)*hvap
#endif
     endif
   enddo
!
!  evaporation over bare sea ice
!
   do i = its,ite
     if(slimsk(i).eq.2.) then
       evap(i) = partlnd(i) * ep(i)
     endif
   enddo
!
!  treat downward moisture flux situation
!  (evap was preset to zero so no update needed)
!  dew is converted from kg m-2 to m to conform to precip unit
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.ep(i).le.0.
#ifndef HYDRO
     dew(i) = 0.
#endif
   enddo
!
   do i = its,ite
     if(flag(i)) then
       dew(i) = -ep(i) * delt / (hvap * rhoh2o)
       evap(i) = ep(i)
       dew(i) = dew(i) * partlnd(i)
       evap(i) = evap(i) * partlnd(i)
       dm(i) = 1.
#ifdef HYDRO
       ec(i)=evap(i)*sigmaf(i)/hvap
       canopy(i)=canopy(i)-ec(i)*delt
       edir(i)=evap(i)*(1.-sigmaf(i))/hvap
       ec(i) = ec(i)*hvap
       evap(i) = edir(i)*hvap + ec(i)
#endif
     endif
   enddo
!
!  snow covered land and sea ice
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.snowd(i).gt.0.
   enddo
!
!  change of snow depth due to evaporation or sublimation
!
!  convert evap from kg m-2 s-1 to m s-1 to determine the reduction of s
!
   do i = its,ite
     if(flag(i)) then
       bfact = snowd(i) / (delt * ep(i) / ((hvap+hfus) * rhoh2o))
       bfact = min(bfact,1.)
!
!  the evaporation of snow
!
       if(ep(i).le.0.) bfact = 1.
       if(snowd(i).le..001) then
         snowev(i) = bfact * ep(i) * (1. - partlnd(i))
#ifdef HYDRO
         snoex=snowev(i)*delt/(rhoh2o * (hvap+hfus))
         if(snoex.gt.snowd(i)) then
               snowev(i)=snowd(i)*(rhoh2o * (hvap+hfus))/delt
         endif
#endif
         evap(i) = evap(i) + snowev(i)
       else
         snowev(i) = bfact * ep(i)
#ifdef HYDRO
         snoex=snowev(i)*delt/(rhoh2o * (hvap+hfus))
         if(snoex.gt.snowd(i)) then
           snowev(i)=snowd(i)*(rhoh2o * (hvap+hfus))/delt
         endif
#endif
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
 500  continue
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.snowd(i).gt..0
   enddo
!
   do i = its,ite
     if(flag(i).and.tsurf(i).gt.t0c) then
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
#ifndef HYDRO
#ifdef ICE
       qss(i) = cb2pa * fpvs(tsurf(i))
#else
       xj1=min(max(c1xpvs0+c2xpvs0*tsurf(i),1.),float(nx))
       jx1=min(xj1,nx-1.)
       fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
       qss(i) = cb2pa * fpvs01
#endif
       qss(i) = eps * qss(i) / (cb2pa * prsi1(i) + epsm1 * qss(i))
       evap(i) = elocp * rch(i) * (qss(i) - q0(i))
#endif
     endif
   enddo
!
!  prepare tendency terms for the soil moisture field without precipitat
!  the unit of moisture flux needs to become m s-1 for soil moisture
!   hence the factor of rhoh2o
!
   do i = its,ite
     flag(i) = slimsk(i).eq.1.
   enddo
!
   do i = its,ite
     smcz(i) = max(smc(i,1), smc(i,2))
      if(flag(i)) then
        df1(i) = osu_funct_df(smcz(i),soiltyp(i))
        kt1(i) = osu_funct_kt(smcz(i),soiltyp(i))
      endif
      if(flag(i).and.stc(i,1).lt.t0c) then
        df1(i) = 0.
        kt1(i) = 0.
      endif
     if(flag(i)) then
!
#ifndef HYDRO
       rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt
#endif
       dmdz(i) = (smc(i,1) - smc(i,2)) / (-.5 * zsoil(i,2))
#ifndef HYDRO
       rhsmc(i,1) = (df1(i) * dmdz(i) + kt1(i)                                 &
               + (edir(i) + et(i,1))) / (zsoil(i,1) * rhoh2o)
       rhsmc(i,1) = rhsmc(i,1) - (1. - sigmaf(i)) * dew(i) /                   &
                    ( zsoil(i,1) * delt )
#else
       rhsmc(i,1) = -(df1(i) * dmdz(i) + kt1(i))                               &
            / (-zsoil(i,1) * rhoh2o)
       hydrow(i,1)=edir(i)+et(i,1)
#endif
       ddz(i) = 1. / (-.5 * zsoil(i,2))
!
!  aim, bim, and cim are the elements of the tridiagonal matrix for the
!  implicit update of the soil moisture
!
#ifndef HYDRO
       aim(i,1) = 0.
       bim(i,1) = df1(i) * ddz(i) / (-zsoil(i,1) * rhoh2o)
       cim(i,1) = -bim(i,1)
#else
       aim(i,1) = df1(i)*ddz(i)
       bim(i,1) = -df1(i)*ddz(i)
       cim(i,1) = -kt1(i)
#endif
     endif
   enddo
   do k = 2,lsoil
     if(k.lt.lsoil) then
       do i = its,ite
         if(flag(i)) then
           df2(i) = osu_funct_df(smcz(i),soiltyp(i))
#ifdef NOBSFLW
!
!- for turning off backgound flow, set kt2(i)=0.
!
           kt2(i) = 0.
#else
           kt2(i) = osu_funct_kt(smcz(i),soiltyp(i))
#endif
         endif
         if(flag(i).and.stc(i,k).lt.t0c) then
           df2(i) = 0.
           kt2(i) = 0.
         endif
         if(flag(i)) then
           dmdz2(i) = (smc(i,k) - smc(i,k+1))                                  &
                      / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
           smcz(i) = max(smc(i,k), smc(i,k+1))
#ifndef HYDRO
           rhsmc(i,k) = (df2(i) * dmdz2(i) + kt2(i)                            &
                        - df1(i) * dmdz(i) - kt1(i) + et(i,k))                 &
                        / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
#else
!
! acr   do something here?  should never get here in a 2-layer model...
!
           rhsmc(i,k) = (df2(i) * dmdz2(i) + kt2(i)                            &
                        - df1(i) * dmdz(i) - kt1(i) )                          &
                        / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
           hydrow(i,k)=hydrow(i,k)+et(i,k)
#endif
           ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
           cim(i,k) = -df2(i) * ddz2(i)                                        &
                     / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
         endif
       enddo
     else
       do i = its,ite
         if(flag(i)) then
#ifdef NOBSFLW
           kt2(i) = 0.
#else
!
!- for turning off backgound flow, set kt2(i)=0.
!
           kt2(i) = osu_funct_kt(smc(i,k),soiltyp(i))
#endif
         endif
         if(flag(i).and.stc(i,k).lt.t0c) kt2(i) = 0.
         if(flag(i)) then
#ifndef HYDRO
           rhsmc(i,k) = (kt2(i)                                                &
                       - df1(i) * dmdz(i) - kt1(i) + et(i,k))                  &
                      / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
#else
!
! acr  adjust rhsmc?  Not used...
!
           rhsmc(i,k) = (df1(i) * dmdz(i) + kt1(i))                            &
                        / (-rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
           hydrow(i,k)=hydrow(i,k)+kt2(i)+et(i,k)
#endif
           drain(i) = kt2(i)
           cim(i,k) = 0.
         endif
       enddo
     endif
!
     do i = its,ite
#ifndef HYDRO
       if(flag(i)) then
         aim(i,k) = -df1(i) * ddz(i)                                           &
                / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
         bim(i,k) = -(aim(i,k) + cim(i,k))
         df1(i) = df2(i)
         kt1(i) = kt2(i)
         dmdz(i) = dmdz2(i)
         ddz(i) = ddz2(i)
       endif
#else
       if(flag(i)) then
         aim(i,k) = -df1(i)*ddz(i)
         bim(i,k) =  df1(i)*ddz(i)
         cim(i,k) = kt1(i)
         df1(i) = df2(i)
         kt1(i) = kt2(i)
         dmdz(i) = dmdz2(i)
         ddz(i) = ddz2(i)
       endif
#endif
     enddo
   enddo
!
 600  continue
!
!  update soil temperature and sea ice temperature
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0.
   enddo
!
!  surface temperature is part of the update when snow is absent
!
   do i = its,ite
     if(flag(i).and..not.flagsnw(i)) then
       yy(i) = t1(i) +                                                         &
             (rcap(i)-gflux(i)- evap(i)) / (rsmall(i) * rch(i))
       zz(i) = 1. + dft0(i) / (-.5 * zsoil(i,1) * rch(i) * rsmall(i))
       xx(i) = dft0(i) * (stsoil(i,1) - yy(i)) / (.5 * zsoil(i,1) * zz(i))
     endif
     if(flag(i).and.flagsnw(i)) then
       yy(i) = stsoil(i,1)
!
!  heat flux from snow is explicit in time
!
       zz(i) = 1.
       xx(i) = dfsnow * (stsoil(i,1) - tsurf(i))                               &
               / (-factsnw(i) * max(snowd(i),.001))
     endif
   enddo
!
!  compute the forcing and the implicit matrix elements for update
!
!  ch2o is the heat capacity of water and csoil is the heat capacity of
!
   do i = its,ite
     if(flag(i)) then
       smcz(i) = max(smc(i,1), smc(i,2))
       dtdz1(i) = (stsoil(i,1) - stsoil(i,2)) / (-.5 * zsoil(i,2))
       if(slimsk(i).eq.1.) then
         dft1(i) = osu_funct_ktsoil(smcz(i),soiltyp(i))
         hcpct(i) = smc(i,1) * ch2o + (1. - smc(i,1)) * csoil
       else
         dft1(i) = dft0(i)
         hcpct(i) = cice
       endif
       dft2(i) = dft1(i)
       ddz(i) = 1. / (-.5 * zsoil(i,2))
!
!  ai, bi, and ci are the elements of the tridiagonal matrix for the
!  implicit update of the soil temperature
!
       ai(i,1) = 0.
       bi(i,1) = dft1(i) * ddz(i) / (-zsoil(i,1) * hcpct(i))
       ci(i,1) = -bi(i,1)
       bi(i,1) = bi(i,1) + dft0(i) / (.5 * zsoil(i,1) **2 * hcpct(i) * zz(i))
       rhstc(i,1) = (dft1(i) * dtdz1(i) - xx(i)) / (zsoil(i,1) * hcpct(i))
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
       if(slimsk(i).eq.1.) then
         hcpct(i) = smc(i,k) * ch2o + (1. - smc(i,k)) * csoil
       elseif(slimsk(i).eq.2.) then
         hcpct(i) = cice
       endif
     enddo
     if(k.lt.lsoil) then
       do i = its,ite
         if(flag(i)) then
           dtdz2(i) = (stsoil(i,k) - stsoil(i,k+1))                            &
                      / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
           smcz(i) = max(smc(i,k), smc(i,k+1))
           if(slimsk(i).eq.1.) then
             dft2(i) = osu_funct_ktsoil(smcz(i),soiltyp(i))
           endif
           ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
           ci(i,k) = -dft2(i) * ddz2(i)                                        &
                     / ((zsoil(i,k-1) - zsoil(i,k)) * hcpct(i))
         endif
       enddo
     else
!
!  at the bottom, climatology is assumed at 2m depth for land and
!  freezing temperature is assumed for sea ice at z(i,lsoil)
!
       do i = its,ite
         if(slimsk(i).eq.1.) then
           dtdz2(i) = (stsoil(i,k) - tg3(i))                                   &
                  / (.5 * (zsoil(i,k-1) + zsoil(i,k)) - zbot)
           dft2(i) = osu_funct_ktsoil(smc(i,k),soiltyp(i))
           ci(i,k) = 0.
         endif
         if(slimsk(i).eq.2.) then
           dtdz2(i) = (stsoil(i,k) - tgice)                                    &
                  / (.5 * zsoil(i,k-1) - .5 * zsoil(i,k))
           dft2(i) = dft1(i)
           ci(i,k) = 0.
         endif
       enddo
     endif
     do i = its,ite
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
   enddo
!
 700  continue
!
!  solve the tri-diagonal matrix
!
   do k = 1,lsoil
     do i = its,ite
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
   do i = its,ite
     if(flag(i)) then
       ci(i,1) = -ci(i,1) / bi(i,1)
       rhstc(i,1) = rhstc(i,1) / bi(i,1)
     endif
   enddo
!
   do k = 2,lsoil
     do i = its,ite
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
   do i = its,ite
    if(flag(i)) then
       ci(i,lsoil) = rhstc(i,lsoil)
     endif
   enddo
!
   do k = lsoil-1,1
     do i = its,ite
       if(flag(i)) then
         ci(i,k) = ci(i,k) * ci(i,k+1) + rhstc(i,k)
       endif
     enddo
   enddo
!
!  update soil and ice temperature
!
   do k = 1,lsoil
     do i = its,ite
       if(flag(i)) then
         stsoil(i,k) = stsoil(i,k) + ci(i,k)
       endif
     enddo
   enddo
!
!  update surface temperature for snow free surfaces
!
   do i = its,ite
     if(slimsk(i).ne.0..and..not.flagsnw(i)) then
       tsurf(i) = t1(i)+(rcap(i)-evap(i))/                                     &
                 (rsmall(i)*rch(i)+dft0(i)/(-.5 * zsoil(i,1)))
     endif
     if(slimsk(i).eq.2..and..not.flagsnw(i)) then
       tsurf(i) = min(tsurf(i),t0c)
     endif
   enddo
   do k = 1,lsoil
     do i = its,ite
       if(slimsk(i).eq.2) then
         stsoil(i,k) = min(stsoil(i,k),t0c)
       endif
     enddo
   enddo
!
!  time filter for soil and skin temperature
!
#ifdef DBG
#ifdef MP
   if(iope) then
#endif
      call print_maxmin_six(tsurf,ite,1,1,1,'tsurf in phys_lsm_osu2')
#ifdef MP
   endif
!
#endif
#endif
   if(inistp.eq.0) then
     do i = its,ite
       if(slimsk(i).ne.0.) then
         tskin(i) = ctfil1 * tsurf(i) + ctfil2 * tskin(i)
         tskin(i) = tsurf(i)
       endif
     enddo
#ifdef DBG
#ifdef MP
     if(iope) then
#endif
       call print_maxmin_six(tskin,ite,1,1,1,'tskin in phys_lsm_osu2')
#ifdef MP
     endif
#endif
#endif
     do k = 1,lsoil
       do i = its,ite
         if(slimsk(i).ne.0.) then
           stc(i,k) = stsoil(i,k)
         endif
       enddo
     enddo
   endif
!
!  gflux calculation
!
   do i = its,ite
     flag(i) = slimsk(i).ne.0..and.flagsnw(i)
   enddo
!
   do i = its,ite
     if(flag(i)) then
       gflux(i) = dfsnow * (stcx(i,1) - tskin(i))                              &
                  /(factsnw(i) * max(snowd(i),.001))
     endif
   enddo
!
   do i = its,ite
    if(slimsk(i).ne.0..and..not.flagsnw(i)) then
      gflux(i) = dft0(i) * (stcx(i,1) - tskin(i))/ (-.5 * zsoil(i,1))
     endif
   enddo
!
!  calculate sensible heat flux
!
   do i = its,ite
     hflx(i) = rch(i) * (tskin(i) - theta1(i))
   enddo
!
!  the rest of the output
!
   do i = its,ite
     qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
     dm(i) = 1.
!
!  convert snow depth back to mm of water equivalent
!
     snoweq(i) = snowd(i) * 1000.
   enddo
#ifdef DBG
!
   if(iope) then
      call print_maxmin_six(qsurf,ite,1,1,1,'qsurf in phys_lsm_osu2')
      call print_maxmin_six(evap,ite,1,1,1,'evap in phys_lsm_osu2')
      call print_maxmin_six(rch,ite,1,1,1,'rch in phys_lsm_osu2')
   endif
#endif
!
   return
   end subroutine phys_lsm_osu2
#endif /* ~OSU end */
!
