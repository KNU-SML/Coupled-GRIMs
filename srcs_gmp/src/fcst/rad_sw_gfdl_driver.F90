#include <define.h>
!-------------------------------------------------------------------------------
   subroutine rad_sw_gfdl_driver(ipts,s0,isrc,ibnd,pl,ta,wa,oa,co2,cosz,taucl, & 
               ccly,cfac,icfc,icwp,                                            &
#ifdef ICECLOUD
               cwp,cip,rew,rei,fice,                                           &
#endif
               albuvb,albuvd,albirb,albird,paer,                               &
               htrc,tupfxc,tdnflx,supfxc,sdnfxc,                               &
               tupfx0,supfx0,sdnfx0,                                           &
#ifdef CHEM
               lat,ibeg,                                                       &
#endif
#ifdef NIM_DIAG
               sdnfvb,sdnfvd,sdnfnb,sdnfnd,fnet)
#else
               sdnfvb,sdnfvd,sdnfnb,sdnfnd)
#endif
!-------------------------------------------------------------------------------
!
! subprogram: rad_sw_gfdl_driver      computes short-wave radiative heating
!
! abstract: this code is a modified version of m.d. chous sw
!   radiation code to fit nmc mrf and climate models.  it computes
!   sw atmospheric absorption and scattering effects due to o3,
!   h2o,co2,o2,clouds, and aerosols, etc.
!   it has 4 uv+vis bands and 3 nir bands (10 k-values each).
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! references: chou (1986, j. clim. appl.meteor.)
!             chou (1990, j. clim.), and chou (1992, j. atms. sci.)
!
! usage:       call rad_sw_gfdl_driver
!
! input parameters:
!   s0     : solar constant
!   isrc   : flags for selecting absorbers
!            1:aerosols, 2:o2, 3:co2, 4:h2o, 5:o3
!            =0:without it,  =1: with it.
!   pl     : model level pressure in mb
!   ta     : model layer temperature in k
!   wa     : layer specific humidity in gm/gm
!   oa     : layer ozone concentration in gm/gm
!   co2    : co2 mixing ration by volumn
!   cosz   : cosine of solar zenith angle
!   taucl  : optical depth of cloud layers
!   ccly   : layer cloud fraction
!   cfac   : fraction of clear sky view at the layer interface
!   icfc   : =0 no cloud factor to weigh clear and cloudy fluxes
!            =1 use cloud factor to weigh clear and cloudy fluxes
!   icwp   : flag indicates the method used for cloud properties
!            calculations, =0 use t-p; =1 use cwc/cic.
!   cwp    : layer cloud water path (g/m**2)
!   cip    : layer cloud ice path (g/m**2)
!   rew    : layer water cloud drop effective radius (micron)
!   rei    : layer ice cloud drop effective radius
!   fice   : fraction of cloud ice content
!   albuvb : uv+vis surf direct albedo
!   albuvd : uv+vis surf diffused albedo
!   albirb : nir surf direct albedo
!   albird : nir surf diffused albedo
!   paer   : aerosol profiles (fraction)
!
! output parameter:
!   htrc   : heating rates for cloudy sky in  k/day
!   tupfxc : upward flux at toa for cloudy sky  w/m**2
!   tdnflx : dnward flux at toa for all sky  w/m**2
!   supfxc : upward flux at sfc for cloudy sky  w/m**2
!   sdnfxc : dnward flux at sfc for cloudy sky  w/m**2
!   tupfx0 : upward flux at toa for clear sky   w/m**2
!   supfx0 : upward flux at sfc for clear sky   w/m**2
!   sdnfx0 : dnward flux at sfc for clear sky   w/m**2
!   sdnfvb : downward surface vis beam flux     w/m**2
!   sdnfnb : downward surface nir beam flux     w/m**2
!   sdnfvd : downward surface vis diff flux     w/m**2
!   sdnfnd : downward surface nir diff flux     w/m**2
!
! note:
!   for all quantities, k=1 is the top level/layer, except
!   si and sl, for which k=1 is the surface level/layer.
!
!  ::: structure :::
!
!    [rad_sw_gfdl_driver]
!      |
!      |--- A-1) [rad_ir_ice] *
!      |    A-2) [rad_ir] *
!      |--- B-1) [rad_uv_ice] *
!      |    B-2) [rad_uv] *
!      |
!      |--- [co2_flux_gfdl] *
!      |--- [rad_sw_aeros_tau_ice] *
!      |--- [rad_sw_aeros_tau] *
!      |--- [rad_sw_flux] *
!
!-------------------------------------------------------------------------------
   use paramodel
   use comio
   use rdparm
#ifdef DFS
   use dfsvar, only : iope
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! ---  input
!
   integer              ::  ipts,ibnd,icfc,icwp
   real                 ::  s0,co2
   real                 ::  pl (imbx,lp1), ta(imbx,lp1),   wa(imbx,l)
   real                 ::  oa(imbx,l)
   real                 ::  taucl(imbx,l), ccly(imbx,l), cfac(imbx,lp1)
   real                 ::  cosz(imax)
   real                 ::  albuvb(imax),  albuvd(imax), albirb(imax)
   real                 ::  albird(imax)
   real                 ::  paer(imbx,6)
   integer              ::  isrc(nsrc)
#ifdef ICECLOUD
   real                 ::  fice(imax,l)
   real                 ::  cwp(imax,l),   cip(imax,l), rew(imax,l)
   real                 ::  rei(imax,l)
#endif
#ifdef CHEM
   integer              ::  lat
   integer              ::  ibeg
#endif
!
! ---  output
!
   real                 ::  tupfxc(imax), supfxc(imax), sdnfxc(imax)
   real                 ::  tdnflx(imax)
   real                 ::  tupfx0(imax), supfx0(imax), sdnfx0(imax)
   real                 ::  htrc(imbx,l)
   real                 ::  sdnfvb(imax), sdnfvd(imax), sdnfnb(imax)
   real                 ::  sdnfnd(imax)
   real                 ::  sdn0vb(imax), sdn0vd(imax), sdn0nb(imax)
   real                 ::  sdn0nd(imax)
#ifdef NIM_DIAG
   real                 ::  fnet(imax)
#endif
!
! ---  internal array
!
   real                 ::  fnet0(imbx,lp1), fnetc(imbx,lp1), htr0 (imbx,lp1)
   real                 ::  dflx0(imbx,lp1), dflxc(imax),     dp   (imbx,l)
   real                 ::  scal (imbx,l),   swh  (imbx,lp1), so2  (imbx,lp1)
   real                 ::  wh   (imbx,l),   csm  (imax),  cf0(imax)
   real                 ::  cf1(imax)
   real                 ::  dwsfb0(imax), dwsfd0(imax), dwsfbc(imax)
   real                 ::  dwsfdc(imax)
#ifdef ICECLOUD
   real                 ::  rewi(imax,l), reii(imax,l), oh(imbx,l)
#endif
   logical              ::  daytm(imax)
   integer              ::  i,k,nday
   real                 ::  xa,to2,xx,fac,rcf1,ccc
   !
   integer,save         ::  ifpr
   real                 ::  taucrt
   data taucrt / 0.05 /, ifpr / 0 /
!
!===> ... ibnd=1:use one nir band, =2:use three nir bands
!     data ibnd / 1 /  ! define in rad_initialize
!===> ... begin here
!
   do i = 1,ipts
     dwsfb0(i) = 0.0e0
     dwsfd0(i) = 0.0e0
     dwsfbc(i) = 0.0e0
     dwsfdc(i) = 0.0e0
   enddo
!
   if (ifpr .eq. 0) then
#ifndef RMP
#ifndef NOPRINT
#ifdef OPENMP
!$omp master
#endif
     if (iope) write(6,12) (isrc(i),i=1,nsrc)
12   format(3x,'aerosol, o2, co2, h2o, o3 =',5i3)
#ifdef OPENMP
!$omp end master
#endif
#endif /* ~NOPRINT end */
#endif /* ~RMP end */
     ifpr = 1
   endif
!
   nday = 0
   do i = 1,ipts
     swh (i,1) = 0.0e0
     so2 (i,1) = 0.0e0
     tdnflx(i) = s0 * cosz(i)
     tupfxc(i) = 0.0e0
     tupfx0(i) = 0.0e0
     supfxc(i) = 0.0e0
     supfx0(i) = 0.0e0
     sdnfxc(i) = 0.0e0
     sdnfx0(i) = 0.0e0
     sdn0vb(i) = 0.0e0
     sdn0vd(i) = 0.0e0
     sdn0nb(i) = 0.0e0
     sdn0nd(i) = 0.0e0
     dflxc(i)  = 0.0e0
     cf0(i)    = cfac(i,lp1)
     cf1(i)    = 1.0e0 - cf0(i)
!
!===> ... csm is the effective secant of the solar xenith angle
!
     csm   (i) = 35.0e0/(sqrt(1224.0e0*cosz(i)*cosz(i)+1.0e0))
     daytm(i) = cosz(i) .gt. 0.0e0
     if (daytm(i)) nday = nday + 1
   enddo
!
   if (nday .eq. 0) then
     do i = 1,ipts
       sdnfvb(i) = 0.0e0
       sdnfvd(i) = 0.0e0
       sdnfnb(i) = 0.0e0
       sdnfnd(i) = 0.0e0
     enddo
     do k = 1,l
       do i = 1,ipts
         htrc(i,k) = 0.0e0
       enddo
     enddo
!
     return
   endif
!
   do k = 1,l
     do i = 1,ipts
!
!===> ... layer thickness and pressure scaling function for
!         water vapor absorption
       dp  (i,k) = pl(i,k+1) - pl(i,k)
       scal(i,k) = dp(i,k)                                                     &
               * (0.5e0*(pl(i,k)+pl(i,k+1))/300.0e0)**0.8e0
!
!===> ... scaled absorber amounts for h2o(wh,swh), unit is g/cm**2
!
       wh(i,k) = 1.02e0 * wa(i,k) * scal(i,k)                                  &
!                * exp(0.00135e0*(ta(i,k)-240.0e0))
                 * (1.0e0 + 0.00135e0*(ta(i,k)-240.0e0))
       swh(i,k+1) = swh(i,k) + wh(i,k)
     enddo
   enddo
!
!===> ... initialize fluxes
!
   do k = 1,lp1
     do i = 1,ipts
       fnet0(i,k) = 0.0e0
       fnetc(i,k) = 0.0e0
       dflx0(i,k) = 0.0e0
     enddo
   enddo
!
   if (icfc .eq. 1) then
     do i=1,ipts
       cfac(i,lp1) = 0.0
     enddo
     do k=1,l
       do i=1,ipts
         if (cf1(i) .gt. 0.0) then
           rcf1 = 1.0 / cf1(i)
           cfac(i,k) = (cfac(i,k) - cf0(i)) * rcf1
           ccly(i,k) = ccly(i,k) * rcf1
         endif
       enddo
     enddo
   endif
!
!  if (icwp.ne. 1) then
!    do k = 1,l
!      do i = 1,ipts
!        taucl(i,k) = taucl(i,k) * ccly(i,k)
!      enddo
!    enddo
!  else
!
#ifdef ICECLOUD
   if (icwp.eq. 1) then
     do k = 1,l
       do i = 1,ipts
!0799    ccc = ccly(i,k) * sqrt(ccly(i,k))
         ccc = ccly(i,k)
         cwp(i,k) = cwp(i,k) * ccc
         cip(i,k) = cip(i,k) * ccc
         rewi(i,k) = 1.0 / rew(i,k)
         reii(i,k) = 1.0 / rei(i,k)
       enddo
     enddo
   endif
#endif
!
!===> ... compute nir fluxes
!
   if (isrc(4) .eq. 1) then
     do i = 1,ipts
       dwsfb0(i) = 0.0e0
       dwsfd0(i) = 0.0e0
       dwsfbc(i) = 0.0e0
       dwsfdc(i) = 0.0e0
     enddo
#ifdef ICECLOUD
     call rad_ir_ice(ipts,wh,ta,taucl,csm,daytm,ibnd,fice,                     &
                      isrc(1),paer,albirb,albird,                              &
                      icwp,cwp,cip,ccly,rew,rei,rewi,reii,                     &
                      tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,               &
                      fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#else
     call rad_ir(ipts,wh,taucl,csm,daytm,ibnd,                                 &
                  isrc(1),paer,albirb,albird,                                  &
                  tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,                   &
#ifdef CHEM
                  lat,ibeg,                                                    &
#endif
                  fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#endif
   endif
!
!===> ... save surface nir band fluxes
!
   do i = 1,ipts
     sdnfnb(i) = cf0(i)*dwsfb0(i) + cf1(i)*dwsfbc(i)
     sdnfnd(i) = cf0(i)*dwsfd0(i) + cf1(i)*dwsfdc(i)
   enddo
!
!===> ... compute uv+visible fluxes
!         scaled amounts for o3(wh), unit is (cm-amt)stp for o3.
!
   if (isrc(5) .eq. 1) then
!
     do i = 1,ipts
       dwsfb0(i) = 0.0e0
       dwsfd0(i) = 0.0e0
       dwsfbc(i) = 0.0e0
       dwsfdc(i) = 0.0e0
     enddo
!
#ifdef ICECLOUD
!               scaled amounts for o3(wh), unit : (cm-amt)stp for o3.
!
     xa = 1.02 * 466.7
     do k = 1,l
       do i = 1,ipts
         oh(i,k) = xa * oa(i,k) * dp(i,k) + 1.0e-11
       enddo
     enddo
     call rad_uv_ice(ipts,wh,oh,ta,taucl,csm,daytm,fice,                       &
                      isrc(1),paer,albuvb,albuvd,                              &
                      icwp,cwp,cip,ccly,rew,rei,rewi,reii,                     &
                      tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,               &
                      fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#else
     xa = 1.02e0 * 466.7e0
     do k = 1,l
       do i = 1,ipts
         wh(i,k) = xa * oa(i,k) * dp(i,k)
       enddo
     enddo
     call rad_uv(ipts,wh,taucl,csm,daytm,                                      &
                  isrc(1),paer,albuvb,albuvd,                                  &
                  tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,                   &
#ifdef CHEM
                  lat,ibeg,                                                    &
#endif
                  fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#endif
   end if
!
   do i = 1,ipts
!
!===> ... save surface downward vis band fluxes
!
     sdnfvb(i) = cf0(i)*dwsfb0(i) + cf1(i)*dwsfbc(i)
     sdnfvd(i) = cf0(i)*dwsfd0(i) + cf1(i)*dwsfdc(i)
!
!===> ... compute final fluxes
!
     tupfxc(i) = cf0(i)*tupfx0(i) + cf1(i)*tupfxc(i)
     supfxc(i) = cf0(i)*supfx0(i) + cf1(i)*supfxc(i)
     sdnfxc(i) = cf0(i)*sdnfx0(i) + cf1(i)*sdnfxc(i)
   enddo
!
   do k = 1,lp1
     do i = 1,ipts
       fnetc (i,k) = cf0(i)*fnet0(i,k) + cf1(i)*fnetc(i,k)
     enddo
   enddo
!
!===> ... compute the absorption due to oxygen,chou(1990,j.climate,209-217)
!         scaled amounts for o2(o2,so2), unit is (cm-atm)stp for o2.
!
   if (isrc(2) .eq. 1) then
     do k = 1,l
       do i = 1,ipts
         so2(i,k+1) = so2(i,k) + 165.22e0 * scal(i,k)
       enddo
     enddo
!
!===> ... to2 is the broadband transmission function for oxygen
!         0.0287 is the fraction of solar flux in the o2 bands
!
     do k = 2,lp1
       do i = 1,ipts
         to2 = exp(-0.00027e0 * sqrt(so2(i,k) * csm(i)))
         dflx0(i,k) = 0.0287e0 * (1.0e0 - to2)
       enddo
     enddo
   endif
!
!===> ... table look-up for the absorption due to co2
!         compute scaled amounts for co2(wc,so2).
!
   if (isrc(3) .eq. 1) then
     xa = co2 * 789.e0
     do k = 1,l
       do i = 1,ipts
         so2(i,k+1) = so2(i,k) + xa*scal(i,k)
       enddo
     enddo
!
     call co2_flux_gfdl(ipts,so2,swh,csm,daytm,dflx0)
!
   endif
!
!===> ... adjust for the effect of o2 and co2 on clear sky net fluxe
!
   if (isrc(2).eq.1 .or. isrc(3).eq.1) then
!lear do 180 k=1,lp1
!       do 180 i=1,ipts
!         fnet0(i,k) = fnet0(i,k) - dflx0(i,k)
!180    continue
!
!===> ... adjust for the effect of o2 and co2 on cloud sky net fluxe
!
     do k = 1,l
       do i = 1,ipts
         if (ccly(i,k) .gt. 0.01e0)                                            &
                  dflxc(i) = dflxc(i) + dflx0(i,k)*cfac(i,k)*ccly(i,k)
           fnetc(i,k+1) = fnetc(i,k+1) - dflxc(i)                              &
                    - dflx0(i,k+1)*cfac(i,k+1)
       enddo
     enddo
!
!===> ... adjust for other fluxes
!
     do i = 1,ipts
       xx = dflxc(i) + cf0(i)*dflx0(i,lp1)
       sdnfx0(i) = sdnfx0(i) - dflx0(i,lp1)
       sdnfxc(i) = sdnfxc(i) - xx
       sdnfnb(i) = sdnfnb(i) - xx
     enddo
   endif
!
   if (icfc .eq. 1) then
!
!===> ... compute final fluxes at top and surface
!
     do i = 1,ipts
       sdnfvb(i) = cf0(i)*sdn0vb(i) + cf1(i)*sdnfvb(i)
       sdnfvd(i) = cf0(i)*sdn0vd(i) + cf1(i)*sdnfvd(i)
       sdnfnb(i) = cf0(i)*sdn0nb(i) + cf1(i)*sdnfnb(i)
       sdnfnd(i) = cf0(i)*sdn0nd(i) + cf1(i)*sdnfnd(i)
       tupfxc(i) = cf0(i)*tupfx0(i) + cf1(i)*tupfxc(i)
       supfxc(i) = cf0(i)*supfx0(i) + cf1(i)*supfxc(i)
       sdnfxc(i) = cf0(i)*sdnfx0(i) + cf1(i)*sdnfxc(i)
     enddo
     do k = 1,lp1
       do i = 1,ipts
         fnetc (i,k) = cf0(i)*fnet0(i,k) + cf1(i)*fnetc(i,k)
       enddo
     enddo
   endif
!
!===> ... convert flux unit to w/m**2
!
   do k = 1,lp1
     do i = 1,ipts
!lear  fnet0 (i,k) = fnet0(i,k) * tdnflx(i)
       fnetc (i,k) = fnetc(i,k) * tdnflx(i)
     enddo
   enddo
!
   do i = 1,ipts
     sdnfnb(i) = sdnfnb(i) * tdnflx(i)
     sdnfnd(i) = sdnfnd(i) * tdnflx(i)
     sdnfvb(i) = sdnfvb(i) * tdnflx(i)
     sdnfvd(i) = sdnfvd(i) * tdnflx(i)
     tupfx0(i) = tupfx0(i) * tdnflx(i)
     tupfxc(i) = tupfxc(i) * tdnflx(i)
     supfx0(i) = supfx0(i) * tdnflx(i)
     supfxc(i) = supfxc(i) * tdnflx(i)
     sdnfx0(i) = sdnfx0(i) * tdnflx(i)
     sdnfxc(i) = sdnfxc(i) * tdnflx(i)
#ifdef NIM_DIAG
     fnet(i)   = fnetc(i,1)
#endif
   enddo
!
!===> ... fac is the factor for heating rates (in k/day)
!         if use k/sec, result should be devided by 86400.
!
!  fac = 3.6*24./10.03*.98
   fac = 8.4418744e0
!
   do k = 1,l
     do i = 1,ipts
!lear  htr0(i,k) = (fnet0(i,k)-fnet0(i,k+1)) * fac / dp(i,k)
       htrc(i,k) = (fnetc(i,k)-fnetc(i,k+1)) * fac / dp(i,k)
     enddo
   enddo
!
   return
   end subroutine rad_sw_gfdl_driver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_ir_ice(ipts,wh,ta,taucl,csm,daytm,ibnd,fice,                 &
                         kaer,paer,albb,albd,                                  &
                         icwp,cwp,cip,ccly,rew,rei,rewi,reii,                  &
                         tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,            &
                         fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
!-------------------------------------------------------------------------------
!
!  compute solar flux in the nir region (3 bands, 10-k per band)
!  the nir region has three water vapor bands, ten k's for each band.
!    1.   1000-4400 (/cm)         2.27-10.0 (micron)
!    2.   4400-8200               1.22-2.27
!    3.   8200-14300              0.70-1.22
!
!  input parameters:                           units 
!    wh,ta,taucl,csm,ibnd,fice,kaer,paer,albb,albd
!    icwp,cwp,cip,cclv,rew,rei
!  fixed input data:
!    h2o absorption coefficient (xk)           cm**2/gm
!    k-distribution function    (hk)           fraction
!
!  the following parameters must specified by users:
!    cloud single scattering albedo (sacl)     n/d    
!    cloud asymmetry factor (asycl)            n/d   
!  aerosols optical parameters are obtained from calling
!    subprogram aeros
!
!  output parameters:
!    fnet0  : clear sky net flux
!    fnetc  : cloudy sky net flux
!    tupfxc : cloudy sky upward flux at toa
!    supfxc : cloudy sky upward flux at sfc
!    sdnfxc : cloudy sky downward flux at sfc
!    tupfx0 : clear sky upward flux at toa
!    supfx0 : clear sky upward flux at sfc
!    sdnfx0 : clear sky downward flux at sfc
!    dwsfb0 : clear sky sfc down dir. flux
!    dwsfd0 : clear sky sfc down dif. flux
!    dwsfbc : cloudy sky sfc down dir. flux
!    dwsfdc : cloudy sky sfc down dif. flux
!
!  program history log:
!   1994-06-12   m.d. chou, gla.
!   1995-02-09   yu-tai hou      - recode for nmc models
!   1998-08-03   yu-tai hou      - updated cloud radiative properties
!                calculation. use slingo's method (jas 1989) on water
!                cloud, ebert and curry's method (jgr 1992) on ice cloud.
!   1999-03-25   yu-tai hou      - updated cloud properties use the
!                most recent chou et al. data (j. clim 1998)
!   1999-04-27   yu-tai hou      - updated cloud radiative property
!                calculations use linear t-adjusted method.
!   1999-04-27   yu-tai hou      - updated cloud radiative property
!                calculations use linear t-adjusted method.
!   1999-09-13   yu-tai hou      - updated to chou's june,1999 version
!
!-------------------------------------------------------------------------------
   use rdparm8
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  wh(imax,l),   taucl(imax,l), csm(imax)
   real                 ::  paer(imax,nae)
   real                 ::  albb(imax),   albd(imax),    zth(imax,l)
   real                 ::  cwp(imax,l),  cip(imax,l),   rew(imax,l)
   real                 ::  rei(imax,l)
   real                 ::  ccly(imax,l), ta(imax,lp1),    fice(imax,l)
   real                 ::  rewi(imax,l), reii(imax,l)
   logical              ::  daytm(imax)
!
! --- output
!
   real                 ::  fnet0 (imax,lp1), dwsfb0(imax), dwsfd0(imax)
   real                 ::  fnetc (imax,lp1), dwsfbc(imax), dwsfdc(imax)
   real                 ::  tupfxc(imax),     supfxc(imax), sdnfxc(imax)
   real                 ::  tupfx0(imax),     supfx0(imax), sdnfx0(imax)
!
   integer              ::  ncloud
   logical              ::  cloudy(imax)
!
!     logical lprnt
!
! --- temporary array
!
   real                 ::  upflux(imax,lp1), dwflux(imax,lp1)
   real                 ::  dwsfxb(imax),     dwsfxd(imax)
   real                 ::  tauto (imax,l),   ssato (imax,l),  asyto (imax,l)
   real                 ::  taurs (l),   ssat1 (imax,l),     asyt1 (imax,l)
   real                 ::  tauaer(imax,l),   ssaaer(imax,l),  asyaer(imax,l)
   real                 ::  fffcw (imax,l), ffft1 (imax,l), fffto (imax,l)
   real                 ::  asycw(imax,l),  ssacw(imax,l)
   real,save            ::  xk  (nk0),        hk  (nk0,nrb)
!
!0499 --- t adjusted cld property method
!
   real,save            ::  ssaw0(nrb,2),   ssai0(nrb,2),   asyw0(nrb,2)
   real,save            ::  asyi0(nrb,2)
   real,save            ::  a0w(nrb,2), a1w(nrb,2), b0w(nrb,2), b1w(nrb,2)
   real,save            ::  a0i(nrb,2), a1i(nrb,2), c0w(nrb,2), c1w(nrb,2)
   real,save            ::  b0i(nrb,2), b1i(nrb,2), b2i(nrb,2)
   real,save            ::  c0i(nrb,2), c1i(nrb,2), c2i(nrb,2)
   real                 ::  facw(imax,l),faci(imax,l)
!
   data xk / 0.0010, 0.0133, 0.0422, 0.1334, 0.4217,                           &
             1.3340, 5.6230, 31.620, 177.80, 1000.0 /
   data hk / .01074, .00360, .00411, .00421, .00389,                           &
             .00326, .00499, .00465, .00245, .00145,                           &
             .08236, .01157, .01133, .01143, .01240,                           &
             .01258, .01381, .00650, .00244, .00094,                           &
             .20673, .03497, .03011, .02260, .01336,                           &
             .00696, .00441, .00115, .00026, .00000,                           &
             .29983, .05014, .04555, .03824, .02965,                           &
             .02280, .02321, .01230, .00515, .00239 /
!
   data ssaw0/.7578,.9869,.9997,.9869, .7570,.9868,.9998,.9916/                &
   ,    asyw0/.8678,.8185,.8354,.8315, .8723,.8182,.8354,.8311/                &
   ,    ssai0/.7283,.9442,.9994,.9620, .7368,.9485,.9995,.9750/                &
   ,    asyi0/.9058,.8322,.8068,.8220, .9070,.8304,.8067,.8174/
   real,save            ::  fffrs0,fpmin,fpmax
   data fffrs0 / 0.1 /
   data fpmin,fpmax /1.0e-8, 0.999999/
!
!0499 - t-adjusted cld prop coeff, water cloud
!
   data                                                                        &
       a0w / 1.466e-2, 2.276e-2, 2.654e-2, 2.494e-2                            &
   ,         1.528e-2, 2.286e-2, 2.642e-2, 2.517e-2 /                          &
   ,   a1w / 1.617e+0, 1.451e+0, 1.351e+0, 1.392e+0                            &
   ,         1.611e+0, 1.449e+0, 1.353e+0, 1.386e+0 /                          &
   ,   b0w / 1.708e-1, 5.314e-4,-4.594e-6, 6.473e-3                            &
   ,         1.674e-1, 5.427e-4,-3.306e-6, 3.218e-3 /                          &
   ,   b1w / 7.142e-3, 1.258e-3, 2.588e-5, 6.649e-4                            &
   ,         7.561e-3, 1.263e-3, 2.287e-5, 5.217e-4 /                          &
   ,   c0w / 8.266e-1, 7.507e-1, 7.925e-1, 7.811e-1                            &
   ,         8.344e-1, 7.501e-1, 7.922e-1, 7.808e-1 /                          &
   ,   c1w / 4.119e-3, 6.770e-3, 4.297e-3, 5.034e-3                            &
   ,         3.797e-3, 6.812e-3, 4.323e-3, 5.031e-3 /
!
!0499 - t-adjusted cld prop coeff, ice cloud
!
    data                                                                       &
       a0i / 2.822e-4,-3.248e-5,-3.758e-5,-1.214e-5                            &
   ,         2.712e-4,-4.308e-5,-3.917e-5,-2.456e-5 /                          &
   ,   a1i / 2.491e+0, 2.522e+0, 2.522e+0, 2.520e00                            &
   ,         2.489e+0, 2.523e+0, 2.522e+0, 2.521e00 /                          &
   ,   b0i / 1.853e-1, 2.544e-3,-7.701e-7, 1.461e-2                            &
   ,         1.738e-1, 2.461e-3,-8.979e-7, 7.083e-3 /                          &
   ,   b1i / 1.841e-3, 1.023e-3, 9.849e-6, 4.612e-4                            &
   ,         1.887e-3, 9.436e-4, 8.102e-6, 3.495e-4 /                          &
   ,   b2i /-6.671e-6,-2.266e-6,-.3988e-9,-1.202e-6                            &
   ,        -6.615e-6,-2.107e-6,-.1862e-9,-8.500e-7 /                          &
   ,   c0i / 8.388e-1, 7.572e-1, 7.519e-1, 7.600e-1                            &
   ,         8.414e-1, 7.566e-1, 7.519e-1, 7.566e-1 /                          &
   ,   c1i / 1.519e-3, 1.563e-3, 1.099e-3, 1.275e-3                            &
   ,         1.477e-3, 1.537e-3, 1.097e-3, 1.241e-3 /                          &
   ,   c2i /-6.702e-6,-5.232e-6,-3.081e-6,-4.020e-6                            &
   ,        -6.403e-6,-5.130e-6,-3.070e-6,-3.804e-6 /
!
   do k = 1,l
     do i = 1,ipts
       facw(i,k) = max(0.0, min(10.0,273.15-ta(i,k)))*0.1
       faci(i,k) = max(0.0, min(30.0,263.15-ta(i,k)))/30.0
     enddo
   enddo
!
!===> ... loop over three nir bands
!
   if (ibnd .eq. 1) then
     ibb1 = nrb
     ibb2 = nrb
   else
     ibb1 = 1
     ibb2 = nrb - 1
   end if
!
   do ib = ibb1,ibb2
!
     ib1 = nvb + ib
!
!===> ... layer optical depth due to rayleigh scattering
!
     do k = 1,l
       do i = 1,ipts
         ssaaer(i,k) = 0.0
         asyaer(i,k) = 0.0
         tauaer(i,k) = 0.0
       enddo
     enddo
     ib1 = nvb + ib
     call rad_sw_aeros_tau_ice(ipts,ib1,kaer,paer,tauaer,ssaaer,asyaer,taurs)
     cloudy(:) = .false.
!
!0898 ... get cloud properties from cwp and cip
!
     if (icwp .eq. 1) then
       do k = 1,l
         do i = 1,ipts
           if (ccly(i,k) .gt. 0.0) then
!
! --- t-adj method
!
             tau1=cwp(i,k)*(  facw(i,k) *(a0w(ib,1)+a1w(ib,1)*rewi(i,k))       &
                        +(1.0-facw(i,k))*(a0w(ib,2)+a1w(ib,2)*rewi(i,k)))
             tau2=cip(i,k)*(  faci(i,k) *(a0i(ib,1)+a1i(ib,1)*reii(i,k))       &
                        +(1.0-faci(i,k))*(a0i(ib,2)+a1i(ib,2)*reii(i,k)))
             taucl(i,k) = tau1 + tau2
             ssa1 = 1.0 - (  facw(i,k) *(b0w(ib,1)+b1w(ib,1)*rew(i,k))         &
                        + (1.0-facw(i,k))*(b0w(ib,2)+b1w(ib,2)*rew(i,k)))
             ssa2 = 1.0 - (  faci(i,k) *(b0i(ib,1)                             &
                        + (b1i(ib,1)+b2i(ib,1)*rei(i,k))*rei(i,k))             &
                        + (1.0-faci(i,k))*(b0i(ib,2)                           &
                        + (b1i(ib,2)+b2i(ib,2)*rei(i,k))*rei(i,k)) )
             ssaw1 = ssa1 * tau1
             ssaw2 = ssa2 * tau2
             ssacw(i,k) = ssaw1 + ssaw2
             asy1 =    facw(i,k) *(c0w(ib,1)+c1w(ib,1)*rew(i,k))               &
                        + (1.0-facw(i,k))*(c0w(ib,2)+c1w(ib,2)*rew(i,k))
             asy2 =    faci(i,k) *(c0i(ib,1)                                   &
                        + (c1i(ib,1)+c2i(ib,1)*rei(i,k))*rei(i,k))             &
                        + (1.0-faci(i,k))*(c0i(ib,2)                           &
                        + (c1i(ib,2)+c2i(ib,2)*rei(i,k))*rei(i,k))
             asyw1 = asy1 * ssaw1
             asyw2 = asy2 * ssaw2
             asycw(i,k) = asyw1 + asyw2
             fffcw(i,k) = asy1*asyw1 + asy2*asyw2
             cloudy(i)  = .true.
           else
             taucl(i,k) = 0.0
             ssacw(i,k) = 1.0
             asycw(i,k) = 0.0
             fffcw(i,k) = 0.0
           endif
         enddo
       enddo
     else
       do k = 1,l
         do i = 1,ipts
           if (taucl(i,k) .gt. 0.0) then
             ssa1 = (1.0-fice(i,k)) * (facw(i,k) * ssaw0(ib,1)                 &
                            + (1.0-facw(i,k))* ssaw0(ib,2))
             ssa2 =      fice(i,k)  * (faci(i,k) * ssai0(ib,1)                 &
                            + (1.0-faci(i,k))* ssai0(ib,2))
             ssaw1 = ssa1 * taucl(i,k)
             ssaw2 = ssa2 * taucl(i,k)
             ssacw(i,k) = ssaw1 + ssaw2
             asy1 = (1.0-fice(i,k)) * (facw(i,k) * asyw0(ib,1)                 &
                            + (1.0-facw(i,k))* asyw0(ib,2))
             asy2 =      fice(i,k)  * (faci(i,k) * asyi0(ib,1)                 &
                            + (1.0-faci(i,k))* asyi0(ib,2))
             asyw1 = asy1 * ssaw1
             asyw2 = asy2 * ssaw2
             asycw(i,k) = asyw1 + asyw2
             fffcw(i,k) = asy1*asyw1 + asy2*asyw2
             cloudy(i) = .true.
           else
             ssacw(i,k) = 1.0
             asycw(i,k) = 0.0
             fffcw(i,k) = 0.0
           endif
         enddo
       enddo
     endif
!
     ncloud = 0
     do i = 1,ipts
       if (cloudy(i)) ncloud = ncloud + 1
     enddo
!
!===> ... ik is the index for the k-distribution function (or the
!     absorption coefficient)
!
     do ik = 1,nk0
!
       if (hk(ik,ib) .ge. 0.00001) then
!
!===> ... compute tatal optical thickness, single scattering albedo,
!         and asymmetry factor for clear sky
!
         do k = 1,l
           do i = 1,ipts
             tauwv      = xk(ik)*wh(i,k)
             tauto(i,k) = max(fpmin, tauwv+tauaer(i,k)+taurs(k))
             ssat1(i,k) = ssaaer(i,k)*tauaer(i,k)+taurs(k)
             asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
             ffft1(i,k) = asyaer(i,k)*asyt1(i,k) + fffrs0*taurs(k)
             ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
             tem        = 1.0 / max(fpmin, ssat1(i,k))
             asyto(i,k) = asyt1(i,k) * tem
             fffto(i,k) = ffft1(i,k) * tem
           enddo
         enddo
!
!===> ... clear sky fluxes calculations
!
!          call rad_sw_flux(tauto,ssato,asyto,fffto,csm,zth,albb,albd,
!    1                 upflux,dwflux,dwsfxb,dwsfxd, l, lp1, ipts)
         call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,          &
                             2,upflux,dwflux,dwsfxb,dwsfxd)
!
         do k = 1,lp1
           do i = 1,ipts
             fnet0 (i,k) = fnet0 (i,k)                                         &
                        + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
           enddo
         enddo
         do i = 1,ipts
           tupfx0(i) = tupfx0(i) + upflux(i,1)   * hk(ik,ib)
           supfx0(i) = supfx0(i) + upflux(i,lp1) * hk(ik,ib)
           sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1) * hk(ik,ib)
           dwsfb0(i) = dwsfb0(i) + dwsfxb(i)     * hk(ik,ib)
           dwsfd0(i) = dwsfd0(i) + dwsfxd(i)     * hk(ik,ib)
         enddo
         if (ncloud .gt. 0) then
!
!===> ... compute tatal optical thickness, single scattering albedo,
!         and asymmetry factor for cloudy sky
!
           do k = 1,l
             do i = 1,ipts
               if (taucl(i,k) .ge. 0.001) then
                 tauto(i,k) = taucl(i,k) + tauto(i,k)
                 ssat1(i,k) = ssacw(i,k) + ssat1(i,k)
                 ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
                 tem        = 1.0 / max(fpmin, ssat1(i,k))
                 asyto(i,k) = (asycw(i,k) + asyt1(i,k)) * tem
                 fffto(i,k) = (fffcw(i,k) + ffft1(i,k)) * tem
               endif
             enddo
           enddo
!
!===> ... cloudy sky fluxes calculations
!
!          call rad_sw_flux(tauto,ssato,asyto,fffto,csm,zth,albb,albd,
!    1                 upflux,dwflux,dwsfxb,dwsfxd, l, lp1, ipts)
           call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,        &
                                2,upflux,dwflux,dwsfxb,dwsfxd)
!
           do k = 1,lp1
             do i = 1,ipts
               fnetc(i,k) = fnetc(i,k)                                         &
                                 + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
             enddo
           enddo
           do i = 1,ipts
             tupfxc(i) = tupfxc(i) + upflux(i,1)   * hk(ik,ib)
             supfxc(i) = supfxc(i) + upflux(i,lp1) * hk(ik,ib)
             sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1) * hk(ik,ib)
             dwsfbc(i) = dwsfbc(i) + dwsfxb(i)     * hk(ik,ib)
             dwsfdc(i) = dwsfdc(i) + dwsfxd(i)     * hk(ik,ib)
           enddo
         else
           do k = 1,lp1
             do i = 1,imax
               fnetc(i,k) = fnet0(i,k)
             enddo
           enddo
           do i = 1,ipts
             tupfxc(i) = tupfx0(i)
             supfxc(i) = supfx0(i)
             sdnfxc(i) = sdnfx0(i)
             dwsfbc(i) = dwsfb0(i)
             dwsfdc(i) = dwsfd0(i)
           enddo
         endif
!
       endif
     enddo           ! k-distribution loop ends here(IK loop)
   enddo             ! loop over nir bands ends here
!
   return
   end subroutine rad_ir_ice
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_ir(ipts,wh,taucl,csm,daytm,ibnd,                             &
                     kaer,paer,albb,albd,                                      &
                     tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,                &
#ifdef CHEM
                     lat,ibeg,                                                 &
#endif
                     fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
!-------------------------------------------------------------------------------
!
!  compute solar flux in the nir region (3 bands, 10-k per band)
!  the nir region has three water vapor bands, ten ks for each band.
!    1.   1000-4400 (/cm)         2.27-10.0 (micron)
!    2.   4400-8200               1.22-2.27
!    3.   8200-14300              0.70-1.22
!
!  input parameters:                           units
!    wh,taucl,csm,daytm,ibnd,kaer,paer,albb,albd
!  fixed input data:
!    h2o absorption coefficient (xk)           cm**2/gm
!    k-distribution function    (hk)           fraction
!
!  the following parameters must specified by users:
!    cloud single scattering albedo (sacl)     n/d
!    cloud asymmetry factor (asycl)            n/d
!  aerosols optical parameters are obtained from calling
!    subprogram aeros
!
!  output parameters:
!    fnet0  : clear sky net flux
!    fnetc  : cloudy sky net flux
!    tupfxc : cloudy sky upward flux at toa
!    supfxc : cloudy sky upward flux at sfc
!    sdnfxc : cloudy sky downward flux at sfc
!    tupfx0 : clear sky upward flux at toa
!    supfx0 : clear sky upward flux at sfc
!    sdnfx0 : clear sky downward flux at sfc
!    dwsfb0 : clear sky sfc down dir. flux
!    dwsfd0 : clear sky sfc down dif. flux
!    dwsfbc : cloudy sky sfc down dir. flux
!    dwsfdc : cloudy sky sfc down dif. flux
! 
!-------------------------------------------------------------------------------
   use rdparm
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! --- input
!
   integer              ::  ipts,ibnd,kaer
   real                 ::  wh(imbx,l),  taucl(imbx,l),  csm(imax)
   real                 ::  paer(imbx,nae)
   real                 ::  albb(imax),  albd(imax)
   logical              ::  daytm(imax)
#ifdef CHEM
   integer              ::  lat
   integer              ::  ibeg
#endif
!
! --- output
!
   real                 ::  fnet0 (imbx,lp1), dwsfb0(imax), dwsfd0(imax)
   real                 ::  fnetc (imbx,lp1), dwsfbc(imax), dwsfdc(imax)
   real                 ::  tupfxc(imax),     supfxc(imax), sdnfxc(imax)
   real                 ::  tupfx0(imax),     supfx0(imax), sdnfx0(imax)
!
! --- temporary array
!
   real                 ::  upflux(imbx,lp1), dwflux(imbx,lp1)
!
!    2, dwsfxb(imax),     dwsfxd(imax),       taurs (l)
!
   real                 ::  dwsfxb(imax),     dwsfxd(imax)
   real                 ::  tauto (imbx,l),   ssato (imbx,l),     asyto (imbx,l)
   real                 ::  taurs (l),        ssat1 (imbx,l),     asyt1 (imbx,l)
!
!    4, taut1 (imbx,l),   ssat1 (imbx,l),     asyt1 (imbx,l)
!
   real                 ::  tauaer(imbx,l),   ssaaer(imbx,l),     asyaer(imbx,l)
   real                 ::  r0  (imbx,lp1),   t0  (imbx,lp1),     tb  (imbx,lp1)
   real                 ::  rf  (imbx,lp1),   tf  (imbx,lp1)
   real,save            ::  xk  (nk0),        hk  (nk0,nrb),      sacl(nrb)
   real,save            ::  asycl,fpmin,fpmax
!
   data xk / 0.0010, 0.0133, 0.0422, 0.1334, 0.4217,                           &
             1.3340, 5.6230, 31.620, 177.80, 1000.0 /
   data hk / .01074, .00360, .00411, .00421, .00389,                           &
             .00326, .00499, .00465, .00245, .00145,                           &
             .08236, .01157, .01133, .01143, .01240,                           &
             .01258, .01381, .00650, .00244, .00094,                           &
             .20673, .03497, .03011, .02260, .01336,                           &
             .00696, .00441, .00115, .00026, .00000,                           &
             .29983, .05014, .04555, .03824, .02965,                           &
             .02280, .02321, .01230, .00515, .00239 /
!
   data sacl / 0.98, 0.994, 0.9995, 0.99 /, asycl / 0.843 /
   data fpmin,fpmax /1.0e-6, 0.999999/
!  data ssawv/0.00001/
   integer              ::  ibb1,ibb2,ib,ib1,ik,k,i
   real                 ::  ssacl,tauwv
!
!===> ... loop over three nir bands
!
   if (ibnd .eq. 1) then
     ibb1 = nrb
     ibb2 = nrb
   else
     ibb1 = 1
     ibb2 = nrb - 1
   end if
!
   do ib = ibb1,ibb2
!
!===> ... get aerosols and rayleigh scattering optical properties
!
     ib1 = nvb + ib
#ifdef CHEM
     call rad_sw_aeros_tau(ipts,ib1,kaer,paer,tauaer,ssaaer,asyaer,taurs,lat,ibeg)
#else
     call rad_sw_aeros_tau(ipts,ib1,kaer,paer,tauaer,ssaaer,asyaer,taurs)
#endif
!
     ssacl=sacl(ib)
!
!===> ... ik is the index for the k-distribution function (or the
!     absorption coefficient)
!
     do ik = 1,nk0
!
       if (hk(ik,ib) .lt. 0.00001) cycle
!
!===> ... compute tatal optical thickness, single scattering albedo,
!         and asymmetry factor for clear sky
!
       do k = 1,l
         do i = 1,ipts
           tauwv = xk(ik)*wh(i,k)
!          taut1(i,k) = tauwv + tauaer(i,k) + taurs(k)
           tauto(i,k) = amax1(fpmin,tauwv+tauaer(i,k)+taurs(k))
!          ssat1(i,k) = ssawv*tauwv+ssaaer(i,k)*tauaer(i,k)+taurs(k)
           ssat1(i,k) = ssaaer(i,k)*tauaer(i,k)+taurs(k)
           asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
!          tauto(i,k) = amax1(fpmin, taut1(i,k))
           ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
           asyto(i,k) = asyt1(i,k) / amax1(fpmin, ssat1(i,k))
         enddo
       enddo
!
!===> ... clear sky fluxes calculations
!
       call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,            &
                          2,upflux,dwflux,dwsfxb,dwsfxd)
!
       do k = 1,lp1
         do i = 1,ipts
           fnet0 (i,k) = fnet0 (i,k) + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
         enddo
       enddo
       do i = 1,ipts
         tupfx0(i) = tupfx0(i) + upflux(i,1)  *hk(ik,ib)
         supfx0(i) = supfx0(i) + upflux(i,lp1)*hk(ik,ib)
         sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1)*hk(ik,ib)
         dwsfb0(i) = dwsfb0(i) + dwsfxb(i)*hk(ik,ib)
         dwsfd0(i) = dwsfd0(i) + dwsfxd(i)*hk(ik,ib)
         enddo
!
!===> ... compute tatal optical thickness, single scattering albedo,
!         and asymmetry factor for cloudy sky
!
       do k = 1,l
         do i = 1,ipts
           if (taucl(i,k) .ge. 0.001e0) then
!            tauto(i,k) = taucl(i,k) + taut1(i,k)
             tauto(i,k) = taucl(i,k) + tauto(i,k)
             ssat1(i,k) = ssacl*taucl(i,k) + ssat1(i,k)
             ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
             asyto(i,k) = (asycl*ssacl*taucl(i,k) + asyt1(i,k)) / ssat1(i,k)
           endif
         enddo
       enddo
!
!===> ... cloudy sky fluxes calculations
!
       call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,            &
                          2,upflux,dwflux,dwsfxb,dwsfxd)
!
       do k = 1,lp1
         do i = 1,ipts
           fnetc(i,k) = fnetc(i,k) + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
         enddo
       enddo
       do i = 1,ipts
         tupfxc(i) = tupfxc(i) + upflux(i,1)  *hk(ik,ib)
         supfxc(i) = supfxc(i) + upflux(i,lp1)*hk(ik,ib)
         sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1)*hk(ik,ib)
         dwsfbc(i) = dwsfbc(i) + dwsfxb(i)*hk(ik,ib)
         dwsfdc(i) = dwsfdc(i) + dwsfxd(i)*hk(ik,ib)
       enddo
     enddo
   enddo
!
   return
   end subroutine rad_ir
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_uv_ice(ipts,wz,oz,ta,taucl,csm,daytm,fice,                   &
                         kaer,paer,albb,albd,                                  &
                         icwp,cwp,cip,ccly,rew,rei,rewi,reii,                  &
                         tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,            &
                         fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
!-------------------------------------------------------------------------------
!
!  compute solar flux in the uv+visible region
!  the uv+visible region is grouped into 8 bands:
!    uv-c     (.175-.225);(.225-.245,.260-.280);(.245-.260);
!    uv-b     (.280-.295);(.295-.310);(.310-.320);
!    uv-a     (.320-.400);
!    par      (.400-.700)
!
!  input parameters:                            units
!    wz,oz,ta,taucl,csm,fice,kaer,paer,albb,albd
!    icwp,cwp,cip,cclv,rew,rei
!
!  output parameters:
!    fnet0  : clear sky net flux
!    fnetc  : cloudy sky net flux
!    tupfxc : cloudy sky upward flux at toa
!    supfxc : cloudy sky upward flux at sfc
!    sdnfxc : cloudy sky downward flux at sfc
!    tupfx0 : clear sky upward flux at toa
!    supfx0 : clear sky upward flux at sfc
!    sdnfx0 : clear sky downward flux at sfc
!    dwsfb0 : clear sky sfc down dir. flux
!    dwsfd0 : clear sky sfc down dif. flux
!    dwsfbc : cloudy sky sfc down dir. flux
!    dwsfdc : cloudy sky sfc down dif. flux
!
!  fixed input data:
!    fraction of solar flux contained
!       in the 8 bands (ss)                     fraction
!    rayleigh optical thickness (tauray)        /mb
!    ozone absorption coefficient (ak)          /(cm-atm)stp
!
!  the following parameters must be specified by users:
!    cloud asymmetry factor (asycl)             n/d 
!  aerosol parameters are from subprogram aeros:
!
!  program history log:
!   1994-06-12   m.d. chou, gla.
!   1995-02-09   yu-tai hou      - recode for nmc models
!   1998-08-03   yu-tai hou      - updated cloud radiative properties
!              calculation. use slingo's method (jas 1989) on water
!              cloud, ebert and curry's method (jgr 1992) on ice cloud.
!   1999-03-25   yu-tai hou      - updated cloud properties use the
!              most recent chou et al. data (j. clim 1998)
!   1999-04-27   yu-tai hou      - updated cloud radiative property
!              calculations use linear t-adjusted method.
!   1999-09-13   yu-tai hou      - updated to chou's june,1999 version
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm8
!-------------------------------------------------------------------------------
#include "abort.h"
!
!
!  parameter (nvbb=4)
   integer,parameter    ::  nvbb=8
!
! --- input
!
   real                 ::  oz(imax,l),   taucl(imax,l), albb(imax), albd(imax)
   real                 ::  csm(imax),  zth(imax,l),   paer(imax,nae)
   real                 ::  ta(imax,lp1),   fice(imax,l)
   real                 ::  cwp(imax,l),  cip(imax,l),  rew(imax,l),rei(imax,l)
   real                 ::  ccly(imax,l), wz(imax,l)
   real                 ::  rewi(imax,l), reii(imax,l)
   logical              ::  daytm(imax)
!
! --- output
!
   real                 ::  fnet0 (imax,lp1), dwsfb0(imax), dwsfd0(imax)
   real                 ::  fnetc (imax,lp1), dwsfbc(imax), dwsfdc(imax)
   real                 ::  tupfxc(imax),     supfxc(imax), sdnfxc(imax)
   real                 ::  tupfx0(imax),     supfx0(imax), sdnfx0(imax)
!
! --- temporary array
!
   real                 ::  upflux(imax,lp1), dwflux(imax,lp1)
   real                 ::  dwsfxb(imax),     dwsfxd(imax)
   real                 ::  tauto (imax,l),   ssato (imax,l),   asyto (imax,l)
   real                 ::  taurs (l),   ssat1 (imax,l),   asyt1 (imax,l)
   real                 ::  tauaer(imax,l),   ssaaer(imax,l),   asyaer(imax,l)
   real                 ::  fffcw (imax,l),   ffft1 (imax,l),   fffto (imax,l)
   real                 ::  asycw (imax,l),   ssacw (imax,l)
!
! --- solar flux and absorption coefficients
!
   real,save            ::  ss(nvbb),         ak(nvbb),         wk(nvbb)
!0499
! --- t adjusted cld property method
!
   real,save            ::  a0w(2), a1w(2), b0w(2), b1w(2), b0i(2), b1i(2), b2i(2)
   real,save            ::  a0i(2), a1i(2), c0w(2), c1w(2), c0i(2), c1i(2), c2i(2)
   real,save            ::  ssaw0(2), ssai0(2), asyw0(2), asyi0(2)
   real                 ::  facw(imax,l), faci(imax,l)
!
   logical              ::  cloudy(imax)
   integer              ::  ncloud
!
!     logical lprnt
!
   data ss / 0.00057, 0.00367, 0.00083, 0.00417,                               &
             0.00600, 0.00556, 0.05913, 0.39081 /
   data ak / 30.47, 187.2, 301.9, 42.83,                                       &
             7.090, 1.250, .0345, .0572 /
   data wk / 7*0.0e0, 0.75e-3 /
   data ssaw0 /.999998,.999998/, ssai0 /.999994,.999995/                       &
       ,asyw0 / 0.853,  0.853 /, asyi0 / 0.7991, 0.7998/
   real,save            ::  fffrs0, fpmin, fpmax
   data fffrs0 / 0.1 /
   data fpmin, fpmax / 1.0e-8, 0.999999 /
!
!0898 - coeff for water cloud
!
!
! --- t adjusted water/ice cloud coeff.
!
   data                                                                        &
      a0w / 0.2807e-1,0.2798e-1 /, a1w / 0.1307e+1,0.1309e+1 /                 &
   ,  b0w / -.1176e-6,-.1810e-6 /, c0w / 0.8276e+0,0.8272e+0 /                 &
   ,  b1w / 0.1770e-6,0.1778e-6 /, c1w / 0.2541e-2,0.2565e-2 /                 &
   ,  a0i / -.3011e-4,-.5975e-5 /, a1i / 0.2519e+1,0.2517e+1 /                 &
   ,  b0i / 0.1688e-6,0.1721e-6 /, c0i / 0.7473e+0,0.7480e+0 /                 &
   ,  b1i / 0.9936e-7,0.9177e-7 /, c1i / 0.1015e-2,0.1015e-2 /                 &
   ,  b2i /-.1114e-10,-.1125e-10/, c2i / -.2524e-5,-.2531e-5 /
!
   do k = 1,l
     do i = 1,ipts
       facw(i,k) = max(0.0, min(10.0,273.15-ta(i,k)))*0.1
       faci(i,k) = max(0.0, min(30.0,263.15-ta(i,k)))/30.0
     enddo
   enddo
   cloudy(:) = .false.
!
   if (nvb .ne. nvbb) then
     write(6,*) ' nvb=',nvb,' nvbb=',nvbb,' run stopped'
     call MPABORT
   endif
!
   if (icwp .ne. 1) then
     do k = 1,l
       do i = 1,ipts
         if (taucl(i,k) .gt. 0.0) then
!
!0499 - t-adj prop from specified ssa and asy
!
           ssa1 = (1.0-fice(i,k))*(facw(i,k) *ssaw0(1)                         &
                          + (1.0-facw(i,k))*ssaw0(2) )
           ssa2 =        fice(i,k) *(faci(i,k) *ssai0(1)                       &
                          + (1.0-faci(i,k))*ssai0(2) )
           ssaw1 = ssa1 * taucl(i,k)
           ssaw2 = ssa2 * taucl(i,k)
           ssacw(i,k) = ssaw1 + ssaw2
           asy1 = (1.0-fice(i,k))*(facw(i,k) *asyw0(1)                         &
                          + (1.0-facw(i,k))*asyw0(2) )
           asy2 =        fice(i,k) *(faci(i,k) *asyi0(1)                       &
                          + (1.0-faci(i,k))*asyi0(2) )
           asyw1 = asy1 * ssaw1
           asyw2 = asy2 * ssaw2
           asycw(i,k) = asyw1 + asyw2
           fffcw(i,k) = asy1*asyw1 + asy2*asyw2
           cloudy(i) = .true.
         else
           ssacw(i,k) = 1.0
           asycw(i,k) = 0.0
           fffcw(i,k) = 0.0
         endif
       enddo
     enddo
   else
     do k = 1,l
       do i = 1,ipts
         if (ccly(i,k) .gt. 0.01) then
!
!0499 --- t-adj prop from ice/water paths
!
           tau1 = cwp(i,k)*(   facw(i,k) *(a0w(1)+a1w(1)*rewi(i,k))            &
                        +(1.-facw(i,k))*(a0w(2)+a1w(2)*rewi(i,k)))
           tau2 = cip(i,k)*(   faci(i,k) *(a0i(1)+a1i(1)*reii(i,k))            &
                        +(1.-faci(i,k))*(a0i(2)+a1i(2)*reii(i,k)))
           taucl(i,k) = tau1 + tau2
           ssa1 = 1.0 - (   facw(i,k) *(b0w(1)+b1w(1)*rew(i,k))                &
                    + (1.-facw(i,k))*(b0w(2)+b1w(2)*rew(i,k)) )
           ssa2 = 1.0 - (  faci(i,k) *(b0i(1)                                  &
                    + (b1i(1)+b2i(1)*rei(i,k))*rei(i,k))                       &
                    + (1.-faci(i,k))*(b0i(2)                                   &
                    + (b1i(2)+b2i(2)*rei(i,k))*rei(i,k)) )
           ssaw1 = ssa1 * tau1
           ssaw2 = ssa2 * tau2
           ssacw(i,k) = ssaw1 + ssaw2
           asy1 =     facw(i,k) *(c0w(1)+c1w(1)*rew(i,k))                      &
                     + (1.-facw(i,k))*(c0w(2)+c1w(2)*rew(i,k))
           asy2 =     faci(i,k) *(c0i(1)                                       &
                     + (c1i(1)+c2i(1)*rei(i,k))*rei(i,k) )                     &
                     + (1.-faci(i,k))*(c0i(2)                                  &
                     + (c1i(2)+c2i(2)*rei(i,k))*rei(i,k) )
           asyw1 = asy1 * ssaw1
           asyw2 = asy2 * ssaw2
           asycw(i,k) = asyw1 + asyw2
           fffcw(i,k) = asy1*asyw1 + asy2*asyw2
           cloudy(i)  = .true.
         else
           taucl(i,k) = 0.0
           ssacw(i,k) = 1.0
           asycw(i,k) = 0.0
           fffcw(i,k) = 0.0
         endif
       enddo
     enddo
   endif
!
   ncloud = 0
   do i = 1,ipts
     if (cloudy(i)) ncloud = ncloud + 1
   enddo
!
!===> ... integration over spectral bands
!
   do iv = 1,nvb
!
!===> ... layer optical depth due to rayleigh scattering
!
     do k = 1,l
       do i = 1,ipts
         ssaaer(i,k) = 0.0
         asyaer(i,k) = 0.0
         tauaer(i,k) = 0.0
       enddo
     enddo
     call rad_sw_aeros_tau_ice(ipts,iv,kaer,paer,tauaer,ssaaer,asyaer,taurs)
!
!===> ... compute total optical thickness, single scattering albedo,
!         and asymmetry factor for clear sky
!
     do k = 1,l
       do i = 1,ipts
         tauoz = ak(iv)*oz(i,k)
         tauwv = wk(iv)*wz(i,k)
         tauto(i,k) = max(fpmin, tauoz+tauwv+tauaer(i,k)+taurs(k))
         ssat1(i,k) = ssaaer(i,k)*tauaer(i,k) + taurs(k)
         asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
         ffft1(i,k) = asyaer(i,k)*asyt1(i,k) + fffrs0*taurs(k)
!
         ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
         tem        = 1.0 / max(fpmin, ssat1(i,k))
         asyto(i,k) = asyt1(i,k) * tem
         fffto(i,k) = ffft1(i,k) * tem
       enddo
     enddo
!
!===> ... clear sky fluxes calculations
!

!       call rad_sw_flux(tauto,ssato,asyto,fffto,csm,zth,albb,albd,
!    1              upflux,dwflux,dwsfxb,dwsfxd, l, lp1, ipts)
     call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,              &
                       1,upflux,dwflux,dwsfxb,dwsfxd)
!
     do k = 1,lp1
       do i = 1,ipts
         fnet0(i,k) = fnet0(i,k) + (dwflux(i,k) - upflux(i,k))*ss(iv)
       enddo
     enddo
     do i = 1,ipts
       tupfx0(i) = tupfx0(i) + upflux(i,1)   * ss(iv)
       supfx0(i) = supfx0(i) + upflux(i,lp1) * ss(iv)
       sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1) * ss(iv)
       dwsfb0(i) = dwsfb0(i) + dwsfxb(i)     * ss(iv)
       dwsfd0(i) = dwsfd0(i) + dwsfxd(i)     * ss(iv)
     enddo
!
     if (ncloud .gt. 0) then
!
!===> ... compute total optical thickness, single scattering albedo,
!         and asymmetry factor for cloudy sky
!
       do k = 1,l
         do i = 1,ipts
           if (taucl(i,k) .gt. 0.0) then
             tauto(i,k) = taucl(i,k) + tauto(i,k)
             ssat1(i,k) = ssacw(i,k) + ssat1(i,k)
             ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
             tem        = 1.0  / max(fpmin, ssat1(i,k))
             asyto(i,k) = (asycw(i,k) + asyt1(i,k)) * tem
             fffto(i,k) = (fffcw(i,k) + ffft1(i,k)) * tem
           endif
         enddo
       enddo
!
!===> ... cloudy sky fluxes calculations
!
!      call rad_sw_flux(tauto,ssato,asyto,fffto,csm,zth,albb,albd,
!    1              upflux,dwflux,dwsfxb,dwsfxd, l, lp1, ipts)
       call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,            &
                          1,upflux,dwflux,dwsfxb,dwsfxd)
!
       do k = 1,lp1
         do i = 1,ipts
           fnetc(i,k) = fnetc(i,k) + (dwflux(i,k) - upflux(i,k))*ss(iv)
         enddo
       enddo
       do i = 1,ipts
         tupfxc(i) = tupfxc(i) + upflux(i,1)   * ss(iv)
         supfxc(i) = supfxc(i) + upflux(i,lp1) * ss(iv)
         sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1) * ss(iv)
         dwsfbc(i) = dwsfbc(i) + dwsfxb(i)     * ss(iv)
         dwsfdc(i) = dwsfdc(i) + dwsfxd(i)     * ss(iv)
       enddo
     else
       do k = 1,lp1
         do i = 1,ipts
           fnetc(i,k) = fnet0(i,k)
         enddo
       enddo
       do i = 1,ipts
         tupfxc(i) = tupfx0(i)
         supfxc(i) = supfx0(i)
         sdnfxc(i) = sdnfx0(i)
         dwsfbc(i) = dwsfb0(i)
         dwsfdc(i) = dwsfd0(i)
       enddo
     endif
!
   enddo           !    integration over spectral bands loop end
!
   return
   end subroutine rad_uv_ice
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_uv(ipts,oz,taucl,csm,daytm,                                  &
                    kaer,paer,albb,albd,                                       &
                    tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,                 &
#ifdef CHEM
                    lat,ibeg,                                                  &
#endif
                    fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
!-------------------------------------------------------------------------------
!
!  compute solar flux in the uv+visible region
!  the uv+visible region is grouped into 4 bands:
!  (.225-.285);(.175-.225,.285-.300);(.300-.325);(.325-.690)
!
!  input parameters:                            units
!    oz,taucl,csm,daytm,kaer,paer,albb,albd
!
!  output parameters:
!    fnet0  : clear sky net flux
!    fnetc  : cloudy sky net flux
!    tupfxc : cloudy sky upward flux at toa
!    supfxc : cloudy sky upward flux at sfc
!    sdnfxc : cloudy sky downward flux at sfc
!    tupfx0 : clear sky upward flux at toa
!    supfx0 : clear sky upward flux at sfc
!    sdnfx0 : clear sky downward flux at sfc
!    dwsfb0 : clear sky sfc down dir. flux
!    dwsfd0 : clear sky sfc down dif. flux
!    dwsfbc : cloudy sky sfc down dir. flux
!    dwsfdc : cloudy sky sfc down dif. flux
!
!  fixed input data:
!    fraction of solar flux contained
!       in the 8 bands (ss)                     fraction
!    rayleigh optical thickness (tauray)        /mb
!    ozone absorption coefficient (ak)          /(cm-atm)stp
!
!  the following parameters must be specified by users:
!    cloud asymmetry factor (asycl)             n/d
!  aerosol parameters are from subprogram aeros:
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  oz(imbx,l), taucl(imbx,l), albb(imax), albd(imax)
   real                 ::  csm(imax),  paer(imbx,nae)
   logical              ::  daytm(imax)
#ifdef CHEM
   integer              ::  lat
   integer              ::  ibeg
#endif
!
! --- output
!
   real                 ::  fnet0 (imbx,lp1), dwsfb0(imax), dwsfd0(imax)
   real                 ::  fnetc (imbx,lp1), dwsfbc(imax), dwsfdc(imax)
   real                 ::  tupfxc(imax),     supfxc(imax), sdnfxc(imax)
   real                 ::  tupfx0(imax),     supfx0(imax), sdnfx0(imax)
!
! --- temporary array
!
   real                 ::  upflux(imbx,lp1), dwflux(imbx,lp1)
!
!    2, dwsfxb(imax),     dwsfxd(imax),     taurs (l)
!
   real                 ::  dwsfxb(imax),     dwsfxd(imax)
   real                 ::  tauto (imbx,l),   ssato (imbx,l),   asyto (imbx,l)
   real                 ::  taurs (l),        ssat1 (imbx,l),   asyt1 (imbx,l)
!
!    4, taut1 (imbx,l),   ssat1 (imbx,l),   asyt1 (imbx,l)
!
   real                 ::  tauaer(imbx,l),   ssaaer(imbx,l),   asyaer(imbx,l)
   real                 ::  r0  (imbx,lp1),   t0  (imbx,lp1),   tb  (imbx,lp1)
   real                 ::  rf  (imbx,lp1),   tf  (imbx,lp1)
!
! --- solar flux and absorption coefficients
!
   real                 ::  ss(nvb),          ak(nvb)
!
   data ss / 0.00530, 0.00505, 0.01109, 0.44498 /
   data ak / 0.1805e+3, 0.267e+2, 0.199e+1, 0.050 /
   data asycl / 0.843 /
   data fpmin, fpmax / 1.0e-8, 0.9999999 /
!
!===> ... integration over spectral bands
!
   do iv = 1,nvb
!
!===> ... get aerosols and rayleigh scattering optical properties
!
#ifdef CHEM
     call rad_sw_aeros_tau(ipts,iv,kaer,paer,tauaer,ssaaer,asyaer,taurs,lat,ibeg)
#else
     call rad_sw_aeros_tau(ipts,iv,kaer,paer,tauaer,ssaaer,asyaer,taurs)
#endif
!
!===> ... compute total optical thickness, single scattering albedo,
!         and asymmetry factor for clear sky
!
     do k = 1,l
       do i = 1,ipts
!        taut1(i,k) = ak(iv)*oz(i,k) + tauaer(i,k) + taurs(k)
         tauto(i,k) = amax1(fpmin,ak(iv)*oz(i,k)+tauaer(i,k)+taurs(k))
         ssat1(i,k) = ssaaer(i,k)*tauaer(i,k) + taurs(k)
         asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
!        tauto(i,k) = amax1(fpmin, taut1(i,k))
         ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
         asyto(i,k) = asyt1(i,k) / amax1(fpmin, ssat1(i,k))
       enddo
     enddo
!
!===> ... clear sky fluxes calculations
!
     call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,              &
                       1,upflux,dwflux,dwsfxb,dwsfxd)
!
     do k = 1,lp1
       do i = 1,ipts
         fnet0(i,k) = fnet0(i,k) + (dwflux(i,k) - upflux(i,k))*ss(iv)
       enddo
     enddo
     do i = 1,ipts
       tupfx0(i) = tupfx0(i) + upflux(i,1)  *ss(iv)
       supfx0(i) = supfx0(i) + upflux(i,lp1)*ss(iv)
       sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1)*ss(iv)
       dwsfb0(i) = dwsfb0(i) + dwsfxb(i)*ss(iv)
       dwsfd0(i) = dwsfd0(i) + dwsfxd(i)*ss(iv)
     enddo
!
!===> ... compute total optical thickness, single scattering albedo,
!         and asymmetry factor for cloudy sky
!
     do k = 1,l
       do i = 1,ipts
         if (taucl(i,k) .gt. 0.0e0) then
!          tauto(i,k) = taucl(i,k) + taut1(i,k)
           tauto(i,k) = taucl(i,k) + tauto(i,k)
           ssat1(i,k) = taucl(i,k) + ssat1(i,k)
           ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
           asyto(i,k) = (asycl*taucl(i,k)+asyt1(i,k)) / ssat1(i,k)
         endif
       enddo
     enddo
!
!===> ... cloudy sky fluxes calculations
!
     call rad_sw_flux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,              &
                       1,upflux,dwflux,dwsfxb,dwsfxd)
!
     do k = 1,lp1
       do i = 1,ipts
         fnetc(i,k) = fnetc(i,k) + (dwflux(i,k) - upflux(i,k))*ss(iv)
       enddo
     enddo
     do i = 1,ipts
       tupfxc(i) = tupfxc(i) + upflux(i,1)  *ss(iv)
       supfxc(i) = supfxc(i) + upflux(i,lp1)*ss(iv)
       sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1)*ss(iv)
       dwsfbc(i) = dwsfbc(i) + dwsfxb(i)*ss(iv)
       dwsfdc(i) = dwsfdc(i) + dwsfxd(i)*ss(iv)
     enddo
!
   enddo
!
   return
   end subroutine rad_uv
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine co2_flux_gfdl(ipts,swc,swh,csm,daytm,dflx) 
!-------------------------------------------------------------------------------
   use rdparm
!-------------------------------------------------------------------------------
!  compute the absorption due to co2. ref: chou (j. climate, 1990,              
!     209-217)                                                                  
!  the effect of co2 absorption below the cloud top is neglected.               
!  input variables:                                                             
!     swc,swh     : column amount of co2 and water vapor                        
!     csm         : secant of solar zenith angle                                
!     daytm       : daytime flag                                                
!  output variables:                                                            
!     dflx        : flux reduction due to co2 for clear sky                     
!                                                                               
!-------------------------------------------------------------------------------
!  real                 ::  csm(imax),     swc(imbx,l),   swh(imbx,l)
   real                 ::  csm(imax),     swc(imbx,lp1),   swh(imbx,lp1)
   real                 ::  dflx(imbx,lp1),cah(22,19)
   logical              ::  daytm(imax)
!
! ...  co2 look-up table .............                                          
!
    data      (cah(i,1),  i=1,22)    / 0.9923,0.9922,0.9921,0.9920,            &
    0.9916,0.9910,0.9899,0.9882,0.9856,0.9818,0.9761,0.9678,0.9558,            &
    0.9395,0.9188,0.8945,0.8675,0.8376,0.8029,0.7621,0.7154,0.6647 /           &
   ,          (cah(i,2),  i=1,22)    / 0.9876,0.9876,0.9875,0.9873,            &
    0.9870,0.9864,0.9854,0.9837,0.9811,0.9773,0.9718,0.9636,0.9518,            &
    0.9358,0.9153,0.8913,0.8647,0.8350,0.8005,0.7599,0.7133,0.6627 /           &
   ,          (cah(i,3),  i=1,22)    / 0.9808,0.9807,0.9806,0.9805,            &
    0.9802,0.9796,0.9786,0.9769,0.9744,0.9707,0.9653,0.9573,0.9459,            &
    0.9302,0.9102,0.8866,0.8604,0.8311,0.7969,0.7565,0.7101,0.6596 /         
!                                                                               
    data      (cah(i,4),  i=1,22)    / 0.9708,0.9708,0.9707,0.9705,            &
    0.9702,0.9697,0.9687,0.9671,0.9647,0.9612,0.9560,0.9483,0.9372,            &
    0.9221,0.9027,0.8798,0.8542,0.8253,0.7916,0.7515,0.7054,0.6551 /           &
   ,          (cah(i,5),  i=1,22)    / 0.9568,0.9568,0.9567,0.9565,            &
    0.9562,0.9557,0.9548,0.9533,0.9510,0.9477,0.9428,0.9355,0.9250,            &
    0.9106,0.8921,0.8700,0.8452,0.8171,0.7839,0.7443,0.6986,0.6486 /           &
   ,          (cah(i,6),  i=1,22)    / 0.9377,0.9377,0.9376,0.9375,            &
    0.9372,0.9367,0.9359,0.9345,0.9324,0.9294,0.9248,0.9181,0.9083,            &
    0.8948,0.8774,0.8565,0.8328,0.8055,0.7731,0.7342,0.6890,0.6395 /         
!                                                                               
    data      (cah(i,7),  i=1,22)    / 0.9126,0.9126,0.9125,0.9124,            &
    0.9121,0.9117,0.9110,0.9098,0.9079,0.9052,0.9012,0.8951,0.8862,            &
    0.8739,0.8579,0.8385,0.8161,0.7900,0.7585,0.7205,0.6760,0.6270 /           &
   ,          (cah(i,8),  i=1,22)    / 0.8809,0.8809,0.8808,0.8807,            &
    0.8805,0.8802,0.8796,0.8786,0.8770,0.8747,0.8712,0.8659,0.8582,            &
    0.8473,0.8329,0.8153,0.7945,0.7697,0.7394,0.7024,0.6588,0.6105 /           &
   ,          (cah(i,9),  i=1,22)    / 0.8427,0.8427,0.8427,0.8426,            &
    0.8424,0.8422,0.8417,0.8409,0.8397,0.8378,0.8350,0.8306,0.8241,            &
    0.8148,0.8023,0.7866,0.7676,0.7444,0.7154,0.6796,0.6370,0.5897 /         
!                                                                               
    data      (cah(i,10), i=1,22)    / 0.7990,0.7990,0.7990,0.7989,            &
    0.7988,0.7987,0.7983,0.7978,0.7969,0.7955,0.7933,0.7899,0.7846,            &
    0.7769,0.7664,0.7528,0.7357,0.7141,0.6866,0.6520,0.6108,0.5646 /           &
   ,          (cah(i,11), i=1,22)    / 0.7515,0.7515,0.7515,0.7515,            &
    0.7514,0.7513,0.7511,0.7507,0.7501,0.7491,0.7476,0.7450,0.7409,            &
    0.7347,0.7261,0.7144,0.6992,0.6793,0.6533,0.6203,0.5805,0.5357 /           &
   ,          (cah(i,12), i=1,22)    / 0.7020,0.7020,0.7020,0.7019,            &
    0.7019,0.7018,0.7017,0.7015,0.7011,0.7005,0.6993,0.6974,0.6943,            &
    0.6894,0.6823,0.6723,0.6588,0.6406,0.6161,0.5847,0.5466,0.5034 /         
!                                                                               
    data      (cah(i,13), i=1,22)    / 0.6518,0.6518,0.6518,0.6518,            &
    0.6518,0.6517,0.6517,0.6515,0.6513,0.6508,0.6500,0.6485,0.6459,            &
    0.6419,0.6359,0.6273,0.6151,0.5983,0.5755,0.5458,0.5095,0.4681 /           &
   ,          (cah(i,14), i=1,22)    / 0.6017,0.6017,0.6017,0.6017,            &
    0.6016,0.6016,0.6016,0.6015,0.6013,0.6009,0.6002,0.5989,0.5967,            &
    0.5932,0.5879,0.5801,0.5691,0.5535,0.5322,0.5043,0.4700,0.4308 /           &
   ,          (cah(i,15), i=1,22)    / 0.5518,0.5518,0.5518,0.5518,            &
    0.5518,0.5518,0.5517,0.5516,0.5514,0.5511,0.5505,0.5493,0.5473,            &
    0.5441,0.5393,0.5322,0.5220,0.5076,0.4878,0.4617,0.4297,0.3929 /         
!                                                                               
    data      (cah(i,16), i=1,22)    / 0.5031,0.5031,0.5031,0.5031,            &
    0.5031,0.5030,0.5030,0.5029,0.5028,0.5025,0.5019,0.5008,0.4990,            &
    0.4960,0.4916,0.4850,0.4757,0.4624,0.4441,0.4201,0.3904,0.3564 /           &
   ,          (cah(i,17), i=1,22)    / 0.4565,0.4565,0.4565,0.4564,            &
    0.4564,0.4564,0.4564,0.4563,0.4562,0.4559,0.4553,0.4544,0.4527,            &
    0.4500,0.4460,0.4400,0.4315,0.4194,0.4028,0.3809,0.3538,0.3227 /           &
   ,          (cah(i,18), i=1,22)    / 0.4122,0.4122,0.4122,0.4122,            &
    0.4122,0.4122,0.4122,0.4121,0.4120,0.4117,0.4112,0.4104,0.4089,            &
    0.4065,0.4029,0.3976,0.3900,0.3792,0.3643,0.3447,0.3203,0.2923 /         
!                                                                               
    data      (cah(i,19), i=1,22)    / 0.3696,0.3696,0.3696,0.3696,            &
    0.3696,0.3696,0.3695,0.3695,0.3694,0.3691,0.3687,0.3680,0.3667,            &
    0.3647,0.3615,0.3570,0.3504,0.3409,0.3279,0.3106,0.2892,0.2642 /         
!-------------------------------------------------------------------------------
!                                                                             
! ... table look-up for the absorption due to co2                               
!     0.0343 is the fraction of solar flux in the co2 bands                     
!     df is the absorption of solar radiation due to co2                        
!                                                                               
! --- flux reduction due to co2                                                 
!
   xx = 1.0e0 / 0.3e0                                                      
   do k = 2,lp1                                                             
     do i = 1,ipts                                                            
       if (daytm(i)) then                                                       
         clog = alog10(swc(i,k) / csm(i))                                        
         wlog = alog10(swh(i,k) / csm(i))                                        
         ic = int( (clog+3.15e0) * xx + 1.0e0)                                 
         iw = int( (wlog+4.15e0) * xx + 1.0e0)                                 
         ic = max(2, min(22, ic))                                                
         iw = max(2, min(19, iw))                                                
         ic1 = ic - 1                                                            
         iw1 = iw - 1                                                            
         dc = (3.0 + clog) * xx - float(ic-2)                                    
         dw = (4.0 + wlog) * xx - float(iw-2)                                    
         x1 = cah(1,  iw1) + (cah(1,  iw) - cah(1,  iw1)) * dw                   
         x2 = cah(ic1,iw1) + (cah(ic1,iw) - cah(ic1,iw1)) * dw                   
         y2 = x2 + (cah(ic,iw1) - cah(ic1,iw1)) * dc                             
         dflx(i,k) = dflx(i,k) + 0.0343 * amax1(0.0, x1-y2)                      
       endif                                                                   
     enddo
   enddo
!                                                                               
   return                                                                    
   end subroutine co2_flux_gfdl
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_sw_aeros_tau_ice(ipts,ib,kaer,paer,tau,ssa,asy,taurs)
!-------------------------------------------------------------------------------
!
!  compute aerosols optical properties of six typical profiles
!  in eight uv+vis bands and four nir bands.
!   band: 1. 0.175-0.225 (uv-c)     2. 0.225-0.245;0.260-0.280 (uv-c)
!         3. 0.245-0.260 (uv-c)     4. 0.280-0.295 (uv-b)
!         5. 0.295-0.310 (uv-b)     6. 0.310-0.320 (uv-b)
!         7. 0.320-0.400 (uv-a)     8. 0.400-0.700 (par)
!         9. 2.27 - 4.0  (nir)     10. 1.22 - 2.27 (nir)
!        11. 0.70 - 1.22 (nir)     12. 0.70 - 4.0  (nir)
!
!  references: 
!    wmo report wcp-112 (1986)
!
!  formulations:
!     dtau(k) = zk(k) * dz(k)
!     zk(k)   = zk(0) * (prss(k)/press(0))**(h/hd)   for exp type
!             = const                                for  others
!     z(k)    = -hh * (ln(press(k))-ln(press(k+1)))
!     hh      - atmospheric scale height, a fuction of press
!     hd      - aerosols scale height
!   where zk is ext. coeff.; z is height; the vertical indeces are
!     k=1 at surface for input sigma levels and k=1 at top for all
!     other quantities.
!
!  input parameters:
!     ib,kaer,paer
!
!  output parameters:
!     tau  - optical depth                         n/d
!     ssa  - single scattering albedo              n/d
!     asy  - asymmetry parameter                   n/d
!     taurs- rayleigh scattering optical depth     n/d
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm8
   use comswaer8
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  paer(imbx,nae)
!
! --- output
!
   real                 ::  tau(imbx,l),  ssa(imbx,l),  asy(imbx,l),  taurs(l)
   logical              ::  laer
   real,save            ::  crt1,crt2
   data  crt1,crt2 / 30.0, 0.03333 /
!
!===> ... layer optical depth due to rayleigh scattering
!
   do k = 1,l
     taurs(k) = taur(k,ib)
   enddo
!
   do k = 1,l
     do i = 1,ipts
       ssa(i,k) = 0.0e0
       asy(i,k) = 0.0e0
       tau(i,k) = 0.0e0
     enddo
   enddo
!
   if (kaer .lt. 1) return
!
   do iaer = 1,nae
     laer = .false.
     do i = 1,ipts
       laer = laer .or. paer(i,iaer).gt.0.0e0
     enddo
!    write(6,42) iaer
!42  format(2x,'in aeros: calc aerosol profile =',i3)
     if (.not. laer) cycle
!
!===> ... find aerosol optical depth, single scattering albedo
!         and asymmetry factor
     do k = 1,l
       kk = idm(k,iaer)
       hd = haer(kk,iaer)
       zk = zaer(kk,iaer)
       if (hd .gt. 0.0e0) then
         tau0 = zk * (sig0(k,iaer)**(hh(k)/hd)) * dz(k)
       else 
         tau0 = (zk - hd*hh(k)*alog(sig0(k,iaer))) * dz(k)
       endif
       do i = 1,ipts
         tau(i,k) = tau(i,k) + paer(i,iaer)*raer(ib)*tau0
         ssa(i,k) = ssa(i,k) + paer(i,iaer)*oaer(kk,iaer,ib)
         asy(i,k) = asy(i,k) + paer(i,iaer)*gaer(kk,iaer,ib)
       enddo
     enddo
   enddo
!
!===> ... smooth profile at domain boundaries
!
   do k = 2,l
     do i = 1,ipts
       ratio = 1.0e0
       if (tau(i,k) .gt. 0.0e0) ratio = tau(i,k-1) / tau(i,k)
       tt = tau(i,k) + tau(i,k-1)
       if (ratio .gt. crt1) then
         tau(i,k) = 0.2e0 * tt
         tau(i,k-1) = tt - tau(i,k)
       else if (ratio .lt. crt2) then
         tau(i,k) = 0.8e0 * tt
         tau(i,k-1) = tt - tau(i,k)
       endif
     enddo
   enddo
!
   return
   end subroutine rad_sw_aeros_tau_ice
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#ifdef CHEM
   subroutine rad_sw_aeros_tau(ipts,ib,kaer,paer,tau,ssa,asy,taurs,lat,ibeg)
#else
   subroutine rad_sw_aeros_tau(ipts,ib,kaer,paer,tau,ssa,asy,taurs)
#endif
!-------------------------------------------------------------------------------
!
!  compute aerosols optical properties of six typical profiles
!  in four uv+vis bands and four nir bands.
!   band: 1. 0.225-0.285 (uv)       2. 0.175-0.225;0.285-0.300 (uv)
!         3. 0.300-0.325 (uv)       4. 0.325-0.690 (par)
!         5. 2.27 - 4.0  (nir)      6. 1.22 - 2.27 (nir)
!         7. 0.70 - 1.22 (nir)      8. 0.70 - 4.0  (nir)
!
!  references: 
!    wmo report wcp-112 (1986)
!
!  formulations:
!     dtau(k) = zk(k) * dz(k)
!     zk(k)   = zk(0) * (prss(k)/press(0))**(h/hd)   for exp type
!             = const                                for  others
!     z(k)    = -hh * (ln(press(k))-ln(press(k+1)))
!     hh      - atmospheric scale height, a fuction of press
!     hd      - aerosols scale height
!   where zk is ext. coeff.; z is height; the vertical indeces are
!     k=1 at surface for input sigma levels and k=1 at top for all
!     other quantities.
!
!  input parameters:
!     ib,kaer,paer
!
!  output parameters:
!     tau  - optical depth                         n/d
!     ssa  - single scattering albedo              n/d
!     asy  - asymmetry parameter                   n/d
!     taurs- rayleigh scattering optical depth     n/d
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm
   use comswaer
#ifdef CHEM
   use flexaod, only : outod, outssa, outg
   use diag_aod_mod, only : diagaod, diagssa, diagasm, scaleaod
   use logical_mod, only : lfixedaod, linstaod
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! --- input
!
   integer              ::  ipts,ib,kaer
   real                 ::  paer(imbx,nae)
#ifdef CHEM
   integer              ::  lat
   integer              ::  ibeg
   integer              ::  ir
#endif
!
! --- output
!
   real                 ::  tau(imbx,l),  ssa(imbx,l),  asy(imbx,l),  taurs(l)
!
! local
!
   logical              ::  laer
   real                 ::  crt1,crt2
   data  crt1,crt2 / 30.0, 0.03333 /
!
   real                 ::  tau0,ratio,tt,zk,hd
   integer              ::  i,k,iaer,kk
!
!===> ... layer optical depth due to rayleigh scattering
!
   do k = 1,l
     taurs(k) = taur(k,ib)
   enddo
!
   do k = 1,l
     do i = 1,ipts
       ssa(i,k) = 0.0e0
       asy(i,k) = 0.0e0
       tau(i,k) = 0.0e0
     enddo
   enddo
!
   if (kaer .lt. 1) return
!
   do iaer = 1,nae
     laer = .false.
     do i = 1,ipts
       laer = laer .or. paer(i,iaer).gt.0.0e0
     enddo
     if (.not. laer) cycle
!
!===> ... find aerosol optical depth, single scattering albedo
!         and asymmetry factor
     do k = 1,l
       kk = idm(k,iaer)
       hd = haer(kk,iaer)
       zk = zaer(kk,iaer)
       if (hd .gt. 0.0e0) then
         tau0 = zk * (sig0(k,iaer)**(hh(k)/hd)) * dz(k)
       else
         tau0 = (zk - hd*hh(k)*alog(sig0(k,iaer))) * dz(k)
       endif
       do i = 1,ipts
         tau(i,k) = tau(i,k) + paer(i,iaer)*raer(ib)*tau0
         ssa(i,k) = ssa(i,k) + paer(i,iaer)*oaer(kk,iaer,ib)
         asy(i,k) = asy(i,k) + paer(i,iaer)*gaer(kk,iaer,ib)
       enddo
     enddo
   enddo
!
!===> ... smooth profile at domain boundaries
!
   do k = 2,l
     do i = 1,ipts
       ratio = 1.0e0
       if (tau(i,k) .gt. 0.0e0) ratio = tau(i,k-1) / tau(i,k)
       tt = tau(i,k) + tau(i,k-1)
       if (ratio .gt. crt1) then
         tau(i,k) = 0.2e0 * tt
         tau(i,k-1) = tt - tau(i,k)
       else if (ratio .lt. crt2) then
         tau(i,k) = 0.8e0 * tt
         tau(i,k-1) = tt - tau(i,k)
       endif
     enddo
   enddo
!
#ifdef CHEM
   ! use aod from chemistry module
   if (linstaod) then
     do k = 1,l
       do i = 1,ipts
         ir = i + ibeg - 1
         tau(i,k) = outod(ir,lat,k,ib)
         ssa(i,k) = outssa(ir,lat,k,ib)
         asy(i,k) = outg(ir,lat,k,ib)
       enddo
     enddo
   ! zero aod
   elseif (.not. lfixedaod) then
     do k = 1,l
       do i = 1,ipts
         ssa(i,k) = 0.0e0
         asy(i,k) = 0.0e0
         tau(i,k) = 0.0e0
       enddo
     enddo
   endif
   ! save diagnostics
   do k = 1,l
     do i = 1,ipts
       ir = i + ibeg - 1
       diagaod(ir,lat,k,ib) = diagaod(ir,lat,k,ib) + tau(i,k)
       diagssa(ir,lat,k,ib) = diagssa(ir,lat,k,ib) + ssa(i,k)
       diagasm(ir,lat,k,ib) = diagasm(ir,lat,k,ib) + asy(i,k)
     enddo
   enddo
   do i = 1,ipts
     ir = i + ibeg - 1
     scaleaod(ir,lat,ib) = scaleaod(ir,lat,ib) + 1
   enddo
#endif
!
   return
   end subroutine rad_sw_aeros_tau
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_sw_flux(ipts,tau,ssc,g0,csm,alb,ald,daytm,                   &
                          isbd,upflux,dwflux,dwsfcb,dwsfcd)
!-------------------------------------------------------------------------------
!
!  uses the delta-eddington approximation to compute the bulk
!  scattering properties of a single layer coded following
!  coakley et al.  (jas, 1982)
!
!  inputs:
!    tau: the effective optical thickness
!    ssc: the effective single scattering albedo
!    g0:  the effective asymmetry factor
!    csm: the effective secant of the zenith angle
!    alb: surface albedo for direct radiation
!    ald: surface albedo for diffused radiation
!    daytm: daytime flag
!    isbd: =1 for uv+vis spectral bands
!          =2 for nir spectral bands
!
!  outputs:
!    upflux: upward fluxes
!    dwflux: downward fluxes
!    dwsfcb: downward surface flux direct component
!    dwsfcd: downward surface flux diffused component
!
!-------------------------------------------------------------------------------
   use rdparm
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  tau(imbx,l), ssc(imbx,l), g0(imbx,l)
   real                 ::  csm(imax),   alb(imax),   ald(imax)
   logical              ::  daytm(imax)
!
! --- output
!
   real                 ::  upflux(imbx,lp1),dwflux(imbx,lp1)
   real                 ::  dwsfcb(imax),dwsfcd(imax)
!
! --- temporary
!
   real                 ::  tb (imbx,lp1),r0 (imbx,lp1),t0 (imbx,lp1)
   real                 ::  rf (imbx,lp1),tf (imbx,lp1),zth(imax)
   real                 ::  ttb(imbx,lp1),tdn(imbx,lp1),rup(imbx,lp1)
   real                 ::  tfd(imbx,lp1),rfu(imbx,lp1),rfd(imbx,lp1)
!
   do i = 1,ipts
     tb(i,lp1) = 0.0e0
     r0(i,lp1) = alb(i)
     t0(i,lp1) = 0.0e0
     rf(i,lp1) = ald(i)
     tf(i,lp1) = 0.0e0
     zth(i)    = 1.0e0 / csm(i)
   enddo
!
   do k = 1,l
     do i = 1,ipts
       if (daytm(i)) then
!
!===> ... delta-eddington scaling of single scattering albedo,
!         optical thickness, and asymmetry factor, k & h eqs(27-29)
!
         ff   = g0(i,k) * g0(i,k)
         aa   = 1.0e0 - ff*ssc(i,k)
         taup = tau(i,k) * aa
         sscp = ssc(i,k) * (1.0e0 - ff) / aa
         gp   = g0(i,k) / (1.0e0 + g0(i,k))
!
         oms1 = 1.0e0 - sscp
         ogs1 = 1.0e0 - sscp*gp
         tlam = 3.0e0 * oms1*ogs1
         slam = sqrt(tlam)
         zz   = zth(i) * zth(i)
         den1 = 1.0e0 - tlam*zz
!
!===> ... safety check
!
         den  = sscp / sign(amax1(1.0e-20, abs(den1)), den1)
!
         gama = 0.50e0 * (1.0e0 + 3.0e0*gp*oms1*zz) * den
         alfa = 0.75e0 * zth(i) * (gp + ogs1) * den
         u1   = 1.50e0 * ogs1 / slam
         up1  = u1 + 1.0e0
         um1  = u1 - 1.0e0
         amg  = alfa - gama
         apg  = alfa + gama
!
!===> ... compute layer transmissions and reflections
!         r0   :  layer reflection of the direct beam
!         t0   :  layer diffuse+direct transmission of direct beam
!         rf   :  layer reflection of the diffused radiation
!         tf   :  layer transmission of the diffused radiation
!         tb   :  layer direc transmission of the direct beam
!
         e1  = exp( -taup*slam )
         ue  = u1 * e1
         uepe= ue + e1
         ueme= ue - e1
         den = 1.0e0 / ((up1 + ueme)*(up1 - ueme))
!
         arg = amin1(30.0e0, taup*csm(i))
         tb1  = exp(-arg)
         rf1  = (up1 + uepe) * (um1 - ueme) * den
         tf1  = 4.0e0 * ue * den
         za   = amg * tb1
         r01  = za * tf1 + apg * rf1 - amg
         t01  = za * rf1 + apg * tf1 - (apg - 1.0e0)*tb1
!
         tb(i,k) = amax1(0.0e0, tb1)
         r0(i,k) = amax1(0.0e0, r01)
         t0(i,k) = amax1(0.0e0, t01)
         rf(i,k) = amax1(0.0e0, rf1)
         tf(i,k) = amax1(0.0e0, tf1)
       else
         tb(i,k) = 0.0e0
         r0(i,k) = 0.0e0
         t0(i,k) = 0.0e0
         rf(i,k) = 0.0e0
         tf(i,k) = 0.0e0
       end if
     enddo
   enddo
!
   if (isbd .eq. 2) then
     do k = 1,l
       do i = 1,ipts
         if (daytm(i).and.ssc(i,k).le.0.0001e0) then
           tb(i,k) = exp(-tau(i,k)*csm(i))
           t0(i,k) = tb(i,k)
           r0(i,k) = 0.0e0
           tf(i,k) = exp(-1.66e0*tau(i,k))
           rf(i,k) = 0.0e0
         endif
       enddo
     enddo
   endif
!
   do i = 1,ipts
     tdn(i,1) = t0(i,1)
     rfd(i,1) = rf(i,1)
     tfd(i,1) = tf(i,1)
     ttb(i,1) = tb(i,1)
     ttb(i,l) = 0.0e0
   enddo
!
!===> ... layers added downward starting from top
!
! hoon: REMARK
! in following loop can make underflow because the rf and den could be too small
!
   do k = 2,lp1
     do i = 1,ipts
       if (daytm(i)) then
         den = tf(i,k) / (1.0e0 - rfd(i,k-1) * rf(i,k))
         tdn(i,k) = ttb(i,k-1)*t0(i,k) + (tdn(i,k-1)-ttb(i,k-1)                &
                     + ttb(i,k-1)*r0(i,k)*rfd(i,k-1)) * den
         rfd(i,k) = rf(i,k) + tf(i,k)*rfd(i,k-1) * den
         tfd(i,k) = tfd(i,k-1) * den
         if(abs(tfd(i,k)).lt.1.e-20) tfd(i,k)=0.
         ttb(i,k) = ttb(i,k-1) * tb(i,k)
         if(abs(ttb(i,k)).lt.1.e-20) ttb(i,k)=0.
       endif
     enddo
   enddo
!
!===> ... layers added upward starting from surface
!
   do i = 1,ipts
     rfu(i,lp1) = rf(i,lp1)
     rup(i,lp1) = r0(i,lp1)
   enddo
!
   do k = l,1,-1
     do i = 1,ipts
       if (daytm(i)) then
         den = tf(i,k) / (1.0e0 - rfu(i,k+1) * rf(i,k))
         rup(i,k) = r0(i,k) + ((t0(i,k)-tb(i,k))*rfu(i,k+1)                    &
                     + tb(i,k)*rup(i,k+1)) * den
         rfu(i,k) = rf(i,k) + tf(i,k)*rfu(i,k+1) * den
       endif
     enddo
   enddo
!
!===> ... find upward and downward fluxes
!
   do i = 1,ipts
     if (daytm(i)) then
       upflux(i,1) = rup(i,1)
       dwflux(i,1) = 1.0e0
     else
       upflux(i,1) = 0.0e0
       dwflux(i,1) = 0.0e0
     endif
   enddo
!
   do k = 2,lp1
     do i = 1,ipts
       if (daytm(i)) then
         den = 1.0e0 / (1.0e0 - rfd(i,k-1)*rfu(i,k))
         upflux(i,k) = (ttb(i,k-1)*rup(i,k) +                                  &
                     (tdn(i,k-1)-ttb(i,k-1))*rfu(i,k)) * den
         dwflux(i,k) = ttb(i,k-1) + ((tdn(i,k-1)-ttb(i,k-1))                   &
                     + ttb(i,k-1)*rup(i,k)*rfd(i,k-1)) * den
       else
         upflux(i,k) = 0.0e0
         dwflux(i,k) = 0.0e0
       endif
     enddo
   enddo
!
!===> ... surface downward fluxes
!
   do i = 1,ipts
     dwsfcb(i) = ttb(i,l)
     dwsfcd(i) = dwflux(i,lp1)-dwsfcb(i)
   enddo
!
   return
   end subroutine rad_sw_flux
