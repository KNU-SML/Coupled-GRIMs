#include "define.h"
#ifndef SWRMDC
   subroutine rad_sw_nasa_driver
   end subroutine rad_sw_nasa_driver
#else
!-------------------------------------------------------------------------------
   subroutine rad_sw_nasa_driver(ipts,s0,isrc,ibnd,pl,ta,wa,oa,co2,cosz,taucl, & 
               ccly,cfac,icfc,icwp,cwp,cip,rew,rei,fice,                       &
               albuvb,albuvd,albirb,albird,kprf,idxc,cmix,denn,rh,             &
               htrc,tupfxc,tdnflx,supfxc,sdnfxc,                               &
               tupfx0,supfx0,sdnfx0,                                           &
               sdnfvb,sdnfvd,sdnfnb,sdnfnd)
!-------------------------------------------------------------------------------
!
! subprogram: rad_sw_nasa_driver      computes short-wave radiative heating
!
! abstract: this code is a modified version of m.d. chous sw
!   radiation code to fit nmc mrf and climate models.  it computes
!   sw atmospheric absorption and scattering effects due to o3,
!   h2o,co2,o2,clouds, and aerosols, etc.
!   it has 8 uv+vis bands and 3 nir bands (10 k-values each).
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! references: chou (1986, j. clim. appl.meteor.)
!   chou (1990, j. clim.), and chou (1992, j. atms. sci.)
!   chou and suarez (1999, nasa/tm-1999-104606,vol.15)
!
! usage: call rad_sw_nasa_driver
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
!   kprf   : tropospheric aerosol profile type index
!   idxc,cmix
!          : aerosol component index and mixing ratio
!   denn   : aerosol number densities of 1st and 2nd layers
!   rh     : relative humidity in fraction
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
!    [rad_sw_nasa_driver]
!      |
!      |--- [rad_ir_nasa] *
!      |--- [rad_uv_nasa] *
!      |--- [rad_co2_flux_nasa] *
!      |--- [rad_sw_aeros_tau_nasa] *
!      |--- [rad_sw_flux_nasa] *
!
!-------------------------------------------------------------------------------
   use paramodel
   use comio
   use rdparm99
   use aerparm
   use co2tab_sw
#ifdef DFS
   use dfsvar, only : iope
#endif
!
#include <abort.h>
! ---  input
!
   real                 ::  pl (imbx,lp1), ta(imbx,lp1)
   real                 ::  wa(imbx,l),   oa(imbx,l)
   real                 ::  taucl(imbx,l), ccly(imbx,l)
   real                 ::  cfac(imbx,lp1),cosz(imax)
   real                 ::  albuvb(imax),  albuvd(imax)
   real                 ::  albirb(imax), albird(imax)
   real                 ::  rh  (imbx,l)
   real                 ::  fice(imbx,l)
   real                 ::  cwp(imbx,l),   cip(imax,l)
   real                 ::  rew(imbx,l), rei(imbx,l)
   real                 ::  denn(ndn,imax)
   real                 ::  cmix(nxc,imax)
   integer              ::  isrc(nsrc)
   integer              ::  idxc(nxc,imax),kprf(imax)
!
! ---  output
!
   real                 ::  tupfxc(imax), supfxc(imax)
   real                 ::  sdnfxc(imax), tdnflx(imax)
   real                 ::  tupfx0(imax), supfx0(imax)
   real                 ::  sdnfx0(imax), htrc(imbx,l)
   real                 ::  sdnfvb(imax), sdnfvd(imax)
   real                 ::  sdnfnb(imax), sdnfnd(imax)
   real                 ::  sdn0vb(imax), sdn0vd(imax)
   real                 ::  sdn0nb(imax), sdn0nd(imax)
!
! ---  internal array
!
   real                 ::  fnet0(imbx,lp1), fnetc(imbx,lp1)
   real                 ::  htr0 (imbx,l)
   real                 ::  dflx0(imbx,lp1), dflxc(imbx,lp1)
   real                 ::  dp   (imbx,l)
   real                 ::  scal (imbx,l),   swh  (imbx,lp1)
   real                 ::  so2  (imbx,lp1)
   real                 ::  wh   (imbx,l),   oh   (imbx,l)
   real                 ::  swu  (imbx,lp1)
   real                 ::  cf0  (imax),     cf1  (imax)
   real                 ::  snt(imax), cnt(imax)
   real                 ::  rewi(imax,l), reii(imax,l)
   logical              ::  daytm(imax)
!
!
   data ifpr / 0 /
!
!===> ... ibnd=1:use one nir band, =2:use three nir bands
!     data ibnd / 1 /  ! define in rad_initialize
!===> ... begin here
!
   if (ifpr .eq. 0) then
#ifndef NOPRINT
#ifdef OPENMP
!$omp master
#endif
     if (iope) then
       write(6,12) (isrc(i),i=1,nsrc)
12     format(3x,'aerosol, o2, co2, h2o, o3 =',5i3)
     endif
#ifdef OPENMP
!$omp end master
#endif
#endif
     ifpr = 1
   end if
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
     dflx0(i,1)= 0.0e0
     cf0(i)    = cfac(i,lp1)
     cf1(i)    = 1.0e0 - cf0(i)
     daytm(i) = cosz(i) .gt. 0.0e0
     if (daytm(i)) then
       snt(i) = 1.0e0 / cosz(i) ! snt = secant of solar zenith angle
       nday = nday + 1
     else
       snt(i) = 1.0e3
     end if
     sdnfvb(i) = 0.0e0
     sdnfvd(i) = 0.0e0
     sdnfnb(i) = 0.0e0
     sdnfnd(i) = 0.0e0
     sdn0vb(i) = 0.0e0
     sdn0vd(i) = 0.0e0
     sdn0nb(i) = 0.0e0
     sdn0nd(i) = 0.0e0
   enddo
   if (nday .eq. 0) then
     do k = 1,l
       do i = 1,ipts
         htrc(i,k) = 0.0e0
       enddo
     enddo
!
     return
   end if
!
   tfac = 0.5e0 / 300.0e0
   do k = 1,l
     do i = 1,ipts
!
!===> ... layer thickness and pressure scaling function for
!         water vapor absorption
!
       dp  (i,k) = pl(i,k+1) - pl(i,k)
       scal(i,k) = dp(i,k) * (tfac*(pl(i,k)+pl(i,k+1)))**0.8e0
!
!===> ... scaled absorber amounts for h2o(wh,swh), unit is g/cm**2
!       tem     = 0.00135e0*(ta(i,k)-240.0e0)
!
       wh(i,k) = 1.02e0 * wa(i,k) * scal(i,k)                                  &
             * exp(0.00135e0*(ta(i,k)-240.0e0))
!    1          * (1.0e0 + tem + 0.5e0*tem*tem) + 1.0e-11
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
     do i = 1,ipts
       cfac(i,lp1) = 0.0
     end do
     do k = 1,l
       do i = 1,ipts
         if (cf1(i) .gt. 0.0) then
           rcf1 = 1.0 / cf1(i)
           cfac(i,k) = (cfac(i,k) - cf0(i)) * rcf1
           ccly(i,k) = ccly(i,k) * rcf1
         end if
       end do
     end do
   end if
!
   if (icwp.ne. 1) then
     do k = 1,l
       do i = 1,ipts
!        taucl(i,k) = taucl(i,k) * ccly(i,k)
         taucl(i,k) = taucl(i,k) * ccly(i,k)*sqrt(ccly(i,k))
       end do
     end do
   else
     do k = 1,l
       do i = 1,ipts
         ccc = ccly(i,k) * sqrt(ccly(i,k))
!        ccc = ccly(i,k)
         cwp(i,k) = cwp(i,k) * ccc
         cip(i,k) = cip(i,k) * ccc
         rewi(i,k) = 1.0 / rew(i,k)
         reii(i,k) = 1.0 / rei(i,k)
       end do
     end do
   end if
!
!===> ... compute nir fluxes
!
   if (isrc(4) .eq. 1) then
!
     call rad_ir_nasa(ipts,wh,ta,taucl,cosz,snt,daytm,ibnd,fice,               &
                      isrc(1),kprf,idxc,cmix,denn,rh,albirb,albird,            &
                      icwp,cwp,cip,ccly,rew,rei,rewi,reii,                     &
                      tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,               &
                      fnet0,fnetc,sdn0nb,sdn0nd,sdnfnb,sdnfnd)
!
   end if
!
!===> ... compute uv+visible fluxes
!         scaled amounts for o3(wh), unit is (cm-amt)stp for o3.
!
   if (isrc(5) .eq. 1) then
     xa = 1.02 * 466.7
     do k = 1,l
       do i = 1,ipts
         oh(i,k) = xa * oa(i,k) * dp(i,k) + 1.0e-11
       enddo
     enddo
     call rad_uv_nasa(ipts,wh,oh,ta,taucl,cosz,snt,daytm,fice,                 &
                      isrc(1),kprf,idxc,cmix,denn,rh,albuvb,albuvd,            &
                      icwp,cwp,cip,ccly,rew,rei,rewi,reii,                     &
                      tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,               &
                      fnet0,fnetc,sdn0vb,sdn0vd,sdnfvb,sdnfvd)
!
   end if
!
!===> ... compute the absorption due to oxygen,chou(1990,j.climate,209-217)
!         scaled amounts for o2(o2,so2), unit is (cm-atm)stp for o2.
!         the constant 165.22=(1000/980)*23.14%*(22400/32)
!
   if (isrc(2) .eq. 1) then
     do i = 1,ipts
       cnt(i) = 165.22e0 * snt(i)
     end do
     do k = 1,l
       do i = 1,ipts
         so2(i,k+1) = so2(i,k) + cnt(i) * scal(i,k)
       enddo
     enddo
!
!===> ... compute flux reduction due to oxygen, the constant 0.0633 is
!         the fraction of insolation contained in the oxygen bands.
!         to2 is the broadband transmission function for oxygen
!
     do k = 2,lp1
       do i = 1,ipts
         to2 = exp( -0.145e-3 * sqrt(so2(i,k)) )
         dflx0(i,k) = 0.0633e0 * (1.0e0 - to2)
       end do
     end do
   end if
!
!===> ... table look-up for the absorption due to co2
!         compute scaled amounts for co2(wc,so2).
!         the constant 789=(1000/980)*(44/28.97)*(22400/44)
!
   if (isrc(3) .eq. 1) then
     do i = 1,ipts
       cnt(i)   = co2 * snt(i)
       so2(i,1) = max(so2(i,1), 1.0e-11)
     end do
     do k = 1,l
       do i = 1,ipts
         so2(i,k+1) = so2(i,k) + 789.0 *cnt(i)*scal(i,k)
       enddo
     enddo
!
!===> ... for co2 absorption in spectrum 1.220-2.270 micron
!         both water vapor and co2 absorptions are moderate
!         so2 and swh are the co2 and water vapor amounts
!         integrated from the top of the atmosphere
!
     u1 = -3.0
     du = 0.15
     w1 = -4.0
     dw = 0.15
     do k = 2,lp1
       do i = 1,ipts
         swu(i,k) = log10(so2(i,k))
         swh(i,k) = log10(swh(i,k)*snt(i))
       enddo
     enddo
!
!===> ... dflx0 is the updated flux reduction
!
     call co2_flux_nasa(ipts,swu,u1,du,nu,swh,w1,dw,nw,cah,daytm,dflx0)
!
!===> ... for co2 absorption in spectrum 2.270-10.00 micron
!         where the co2 absorption has a large impact on the
!         heating of middle atmosphere
!
     u1 = 0.250e-3
     du = 0.050e-3
     w1 = -2.0
     dw = 0.05
!
!===> ... co2 mixing ratio is independent of space
!         swh is the logarithm of pressure
!
     do k = 2,lp1
       do i = 1,ipts
         swu(i,k) = cnt(i)
         swh(i,k) = log10(pl(i,k))
       enddo
     enddo
!
!===> ... dflx0 is the updated flux reduction
!
     call co2_flux_nasa(ipts,swu,u1,du,nx,swh,w1,dw,ny,coa,daytm,dflx0)
!
   end if
!
!===> ... adjust for the effect of o2 and co2 on clear sky net fluxe
!
   if (isrc(2).eq.1 .or. isrc(3).eq.1) then
     do k = 1,lp1
       do i = 1,ipts
         fnet0(i,k) = fnet0(i,k) - dflx0(i,k)
       enddo
     enddo
!
!===> ... adjust for the effect of o2 and co2 on cloud sky net fluxe
!
     do i = 1,ipts
       jtop = lp1
!
!===> ... above clouds
!
       do k = 1,lp1
         dflxc(i,k) = dflx0(i,k)
         if (cfac(i,k) .lt. 1.0) then
           jtop = k
           exit
         end if
       end do
!
!===> ... below cloud top
!
       if (jtop .lt. lp1) then
         do k = jtop+1,lp1
           dflxc(i,k) = dflx0(i,k) * (fnetc(i,k)/fnet0(i,k))
         end do
       end if
       do k = 1,lp1
         fnetc(i,k) = fnetc(i,k) - dflxc(i,k)
       end do
     enddo
!
!===> ... adjust for other fluxes
!
      do i = 1,ipts
        sdnfx0(i) = sdnfx0(i) - dflx0(i,lp1)
        sdnfxc(i) = sdnfxc(i) - dflxc(i,lp1)
        sdn0nb(i) = sdn0nb(i) - dflx0(i,lp1)
        sdnfnb(i) = sdnfnb(i) - dflxc(i,lp1)
      enddo
!
   end if
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
   end if
!
!===> ... convert flux unit to w/m**2
!
   do k = 1,lp1
     do i = 1,ipts
!lear  fnet0 (i,k) = fnet0(i,k) * tdnflx(i)
       fnetc (i,k) = fnetc(i,k) * tdnflx(i)
     enddo
   enddo
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
   enddo
!
!===> ... fac is the factor for heating rates (in k/day)
!         if use k/sec, result should be devided by 86400.
!
!     fac = 3.6*24./10.031*.98
   fac = 8.4410328e0
!
   do k = 1,l
     do i = 1,ipts
!lear  htr0(i,k) = (fnet0(i,k)-fnet0(i,k+1)) * fac / dp(i,k)
       htrc(i,k) = (fnetc(i,k)-fnetc(i,k+1)) * fac / dp(i,k)
     enddo
   enddo
!
   return
   end subroutine rad_sw_nasa_driver
!
!-------------------------------------------------------------------------------
   subroutine rad_ir_nasa(ipts,wh,ta,taucl,zth,csm,daytm,ibnd,fice,            &
                          kaer,kprf,idxc,cmix,denn,rh,albb,albd,               &
                          icwp,cwp,cip,ccly,rew,rei,rewi,reii,                 &
                          tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,           &
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
!    wh,ta,taucl,csm,ibnd,fice,kaer,kprf,idxc,cmix,denn,rh,albb,albd
!    icwp,cwp,cip,cclv,rew,rei,zth
!  fixed input data:
!    h2o absorption coefficient (xk)           cm**2/gm
!    k-distribution function    (hk)           fraction
!
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
!   1999-10-13   yu-tai hou      - changed aerosol to opac algorithm
!
!-------------------------------------------------------------------------------
   use rdparm99
   use aerparm
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  wh(imbx,l),   taucl(imbx,l), csm(imax)
   real                 ::  rh (imbx,l)
   real                 ::  albb(imax),   albd(imax),    zth(imax)
   real                 ::  cwp(imbx,l),  cip(imbx,l),   rew(imbx,l)
   real                 ::  rei(imbx,l)
   real                 ::  ccly(imbx,l), ta(imbx,lp1),  fice(imbx,l)
   real                 ::  rewi(imbx,l), reii(imbx,l)
   real                 ::  cmix(nxc,imax),denn(ndn,imax)
   integer              ::  idxc(nxc,imax)
   integer              ::  kprf(imax)
   logical              ::  daytm(imax)
!
! --- output
!
   real                 ::  fnet0 (imbx,lp1), dwsfb0(imax), dwsfd0(imax)
   real                 ::  fnetc (imbx,lp1), dwsfbc(imax), dwsfdc(imax)
   real                 ::  tupfxc(imax),     supfxc(imax), sdnfxc(imax)
   real                 ::  tupfx0(imax),     supfx0(imax), sdnfx0(imax)
!
   integer              ::  ncloud
   logical              ::  cloudy(imax)
!
! --- temporary array
!
   real                 ::  upflux(imbx,lp1), dwflux(imbx,lp1)
   real                 ::  dwsfxb(imax),     dwsfxd(imax)
   real                 ::  tauto (imbx,l),   ssato (imbx,l),  asyto (imbx,l)
   real                 ::  taurs (l),        ssat1 (imbx,l),  asyt1 (imbx,l)
   real                 ::  tauaer(imbx,l),   ssaaer(imbx,l),  asyaer(imbx,l)
   real,save            ::  xk  (nk0),        hk  (nk0,nrb)
!
! --- t adjusted cld property method
!
   real,save            ::  ssaw0(nrb,2),   ssai0(nrb,2),   asyw0(nrb,2)
   real,save            ::  asyi0(nrb,2)
   real                 ::  fffcw(imbx,l),  ffft1(imbx,l),  fffto(imbx,l)
   real                 ::  asycw(imbx,l),  ssacw(imbx,l)
   real,save            ::  a0w(nrb,2), a1w(nrb,2), b0w(nrb,2), b1w(nrb,2)
   real,save            ::  a0i(nrb,2), a1i(nrb,2), c0w(nrb,2), c1w(nrb,2)
   real,save            ::  b0i(nrb,2), b1i(nrb,2), b2i(nrb,2)
   real,save            ::  c0i(nrb,2), c1i(nrb,2), c2i(nrb,2)
   real                 ::  facw(imax,l),faci(imax,l)

   real,save            ::  fffrs0, fpmin, fpmax
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
!     data ssaw0/.7578,.9869,.9997,.9869, .7570,.9868,.9998,.9916/
   data ssaw0/.7578,.9869,.9997,.9869, .7570,.9868,.9998,.9965/                &
   ,    asyw0/.8678,.8185,.8354,.8315, .8723,.8182,.8354,.8311/                &
!    2,    ssai0/.7283,.9442,.9994,.9620, .7368,.9485,.9995,.9750/
   ,    ssai0/.7283,.9442,.9994,.9620, .7368,.9485,.9995,.9823/                &
   ,    asyi0/.9058,.8322,.8068,.8220, .9070,.8304,.8067,.8174/
   data fffrs0 / 0.1 /
   data fpmin,fpmax /1.0e-8, 0.999999/
!
   data                                                                        &
!
! ---- t-adjusted cld prop coeff, water cloud
!
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
   data                                                                        &
!
! --- t-adjusted cld prop coeff, ice cloud
!
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
!===> ... get aerosols and rayleigh scattering optical properties
!
     ib1 = nvb + ib
     call rad_sw_aeros_tau_nasa(ipts,ib1,kaer,kprf,idxc,cmix,denn,rh           &
                             ,tauaer,ssaaer,asyaer,taurs)

     cloudy(:) = .false.
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
                      + (1.0-faci(i,k))*(b0i(ib,2)                             &
                        + (b1i(ib,2)+b2i(ib,2)*rei(i,k))*rei(i,k)) )
             ssaw1 = ssa1 * tau1
             ssaw2 = ssa2 * tau2
             ssacw(i,k) = ssaw1 + ssaw2
             asy1 =    facw(i,k) *(c0w(ib,1)+c1w(ib,1)*rew(i,k))               &
                + (1.0-facw(i,k))*(c0w(ib,2)+c1w(ib,2)*rew(i,k))
             asy2 =    faci(i,k) *(c0i(ib,1)                                   &
                     + (c1i(ib,1)+c2i(ib,1)*rei(i,k))*rei(i,k))                &
                + (1.0-faci(i,k))*(c0i(ib,2)                                   &
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
           end if
         enddo
       enddo
     else
       do k = 1,l
         do i = 1,ipts
           if (ccly(i,k) .gt. 0.0) then
             tau2 = fice(i,k) * taucl(i,k)
             tau1 = taucl(i,k) - tau2
             ssa1 =      facw(i,k) * ssaw0(ib,1)                               &
                  + (1.0-facw(i,k))* ssaw0(ib,2)
             ssa2 =      faci(i,k) * ssai0(ib,1)                               &
                  + (1.0-faci(i,k))* ssai0(ib,2)
             ssaw1 = ssa1 * tau1
             ssaw2 = ssa2 * tau2
             ssacw(i,k) = ssaw1 + ssaw2
             asy1 =      facw(i,k) * asyw0(ib,1)                               &
                  + (1.0-facw(i,k))* asyw0(ib,2)
             asy2 =      faci(i,k) * asyi0(ib,1)                               &
                  + (1.0-faci(i,k))* asyi0(ib,2)
             asyw1 = asy1 * ssaw1
             asyw2 = asy2 * ssaw2
             asycw(i,k) = asyw1 + asyw2
             fffcw(i,k) = asy1*asyw1 + asy2*asyw2
             cloudy(i) = .true.
           else
             ssacw(i,k) = 1.0
             asycw(i,k) = 0.0
             fffcw(i,k) = 0.0
           end if
         enddo
       enddo
     end if
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
         call rad_sw_flux_nasa(ipts,tauto,ssato,asyto,fffto,csm,zth,albb,albd, &
                               daytm,upflux,dwflux,dwsfxb,dwsfxd)
!
         do k = 1,lp1
           do i = 1,ipts
             fnet0 (i,k) = fnet0 (i,k)                                         &
                         + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
           enddo
         enddo
!
         do i = 1,ipts
          tupfx0(i) = tupfx0(i) + upflux(i,1)   * hk(ik,ib)
          supfx0(i) = supfx0(i) + upflux(i,lp1) * hk(ik,ib)
          sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1) * hk(ik,ib)
          dwsfb0(i) = dwsfb0(i) + dwsfxb(i)     * hk(ik,ib)
          dwsfd0(i) = dwsfd0(i) + dwsfxd(i)     * hk(ik,ib)
         enddo
!
         if (ncloud .gt. 0) then
!
!===> ... compute tatal optical thickness, single scattering albedo,
!         and asymmetry factor for cloudy sky
!
           do k = 1,l
             do i = 1,ipts
               if (ccly(i,k) .gt. 0.0) then
                 tauto(i,k) = taucl(i,k) + tauto(i,k)
                 ssat1(i,k) = ssacw(i,k) + ssat1(i,k)
                 ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
                 tem        = 1.0 / max(fpmin, ssat1(i,k))
                 asyto(i,k) = (asycw(i,k) + asyt1(i,k)) * tem
                 fffto(i,k) = (fffcw(i,k) + ffft1(i,k)) * tem
               end if
             enddo
           enddo
!
!===> ... cloudy sky fluxes calculations
!
           call rad_sw_flux_nasa(ipts,tauto,ssato,asyto,fffto,csm,zth,         &
                                 albb,albd,daytm,upflux,dwflux,dwsfxb,dwsfxd)
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
       endif
     enddo           ! k-distribution loop ends here
   enddo           ! loop over nir bands ends here
!
   return
   end subroutine rad_ir_nasa
!-------------------------------------------------------------------------------
   subroutine rad_uv_nasa(ipts,wz,oz,ta,taucl,zth,csm,daytm,fice,              &
                          kaer,kprf,idxc,cmix,denn,rh,albb,albd,               &
                          icwp,cwp,cip,ccly,rew,rei,rewi,reii,                 &
                          tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,           &
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
!    wz,oz,ta,taucl,csm,fice,kaer,kprf,idxc,cmix,denn,rh,albb,albd
!    icwp,cwp,cip,cclv,rew,rei,zth
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
!   1999-10-13   yu-tai hou      - changed aerosol to opac algorithm
!
!-------------------------------------------------------------------------------
   use rdparm99
   use aerparm
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  oz(imbx,l),   taucl(imbx,l), albb(imax), albd(imax)
   real                 ::  csm(imax),    zth(imax),     rh(imbx,l)
   real                 ::  ta(imbx,lp1), fice(imbx,l)
   real                 ::  cwp(imbx,l),  cip(imbx,l),                         &
                            rew(imbx,l), rei(imbx,l)
   real                 ::  ccly(imbx,l), wz(imbx,l),                          &
                            rewi(imbx,l), reii(imbx,l)
   real                 ::  cmix(nxc,imax),denn(ndn,imax)
   integer              ::  kprf(imax),idxc(nxc,imax)
   logical              ::  daytm(imax)
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
   real                 ::  dwsfxb(imax),     dwsfxd(imax)
   real                 ::  tauto (imbx,l),   ssato (imbx,l),   asyto (imbx,l)
   real                 ::  taurs (l),        ssat1 (imbx,l),   asyt1 (imbx,l)
   real                 ::  tauaer(imbx,l),   ssaaer(imbx,l),   asyaer(imbx,l)
   real                 ::  fffcw (imbx,l),   ffft1 (imbx,l),   fffto (imbx,l)
   real                 ::  asycw (imbx,l),   ssacw (imbx,l)
!
! --- solar flux and absorption coefficients
!
   real,save            ::  ss(nvb),          ak(nvb),          wk(nvb)
!
! --- t adjusted cld property method
!
   real,save            ::  a0w(2), a1w(2), b0w(2), b1w(2), b0i(2)
   real,save            ::  b1i(2), b2i(2)
   real,save            ::  a0i(2), a1i(2), c0w(2), c1w(2), c0i(2)
   real,save            ::  c1i(2), c2i(2)
   real,save            ::  ssaw0(2), ssai0(2), asyw0(2), asyi0(2)
   real                 ::  facw(imbx,l), faci(imbx,l)
   real,save            ::  fffrs0, fpmin, fpmax
!
   logical              ::  cloudy(imax)
   integer              ::  ncloud
!
   data ss / 0.00057, 0.00367, 0.00083, 0.00417,                               &
             0.00600, 0.00556, 0.05913, 0.39081 /
   data ak / 30.47, 187.2, 301.9, 42.83,                                       &
             7.090, 1.250, .0345, .0572 /
   data wk / 7*0.0e0, 0.75e-3 /
   data ssaw0 /.999998,.999998/, ssai0 /.999994,.999995/                       &
        asyw0 / 0.853,  0.853 /, asyi0 / 0.7991, 0.7998/
   !
   data fffrs0 / 0.1 /
   data fpmin, fpmax / 1.0e-8, 0.999999 /
!
! --- t adjusted water/ice cloud coeff.
!
   data                                                                        &
      a0w / 0.2807e-1,0.2798e-1 /, a1w / 0.1307e+1,0.1309e+1 /                 &
     ,b0w / -.1176e-6,-.1810e-6 /, c0w / 0.8276e+0,0.8272e+0 /                 &
     ,b1w / 0.1770e-6,0.1778e-6 /, c1w / 0.2541e-2,0.2565e-2 /                 &
     ,a0i / -.3011e-4,-.5975e-5 /, a1i / 0.2519e+1,0.2517e+1 /                 &
     ,b0i / 0.1688e-6,0.1721e-6 /, c0i / 0.7473e+0,0.7480e+0 /                 &
     ,b1i / 0.9936e-7,0.9177e-7 /, c1i / 0.1015e-2,0.1015e-2 /                 &
     ,b2i /-.1114e-10,-.1125e-10/, c2i / -.2524e-5,-.2531e-5 /
!
   do k = 1,l
     do i = 1,ipts
       facw(i,k) = max(0.0, min(10.0,273.15-ta(i,k)))*0.1
       faci(i,k) = max(0.0, min(30.0,263.15-ta(i,k)))/30.0
     enddo
   enddo
   cloudy(:) = .false.
!
   if (icwp .ne. 1) then
     do k = 1,l
       do i = 1,ipts
         if (ccly(i,k) .gt. 0.0) then
!
! --- t-adj prop from specified ssa and asy
!
           tau2 = fice(i,k) * taucl(i,k)
           tau1 = taucl(i,k) - tau2
           ssa1 =      facw(i,k) *ssaw0(1)                                     &
                + (1.0-facw(i,k))*ssaw0(2)
           ssa2 =      faci(i,k) *ssai0(1)                                     &
                + (1.0-faci(i,k))*ssai0(2)
           ssaw1 = ssa1 * tau1
           ssaw2 = ssa2 * tau2
           ssacw(i,k) = ssaw1 + ssaw2
           asy1 =      facw(i,k) *asyw0(1)                                     &
                + (1.0-facw(i,k))*asyw0(2)
           asy2 =      faci(i,k) *asyi0(1)                                     &
                + (1.0-faci(i,k))*asyi0(2)
           asyw1 = asy1 * ssaw1
           asyw2 = asy2 * ssaw2
           asycw(i,k) = asyw1 + asyw2
           fffcw(i,k) = asy1*asyw1 + asy2*asyw2
           cloudy(i) = .true.
         else
           ssacw(i,k) = 1.0
           asycw(i,k) = 0.0
           fffcw(i,k) = 0.0
         end if
       enddo
     enddo
   else
     do k = 1,l
       do i = 1,ipts
         if (ccly(i,k) .gt. 0.0) then
!
! --- t-adj prop from ice/water paths
!
           tau1 = cwp(i,k)*(   facw(i,k) *(a0w(1)+a1w(1)*rewi(i,k))            &
                +(1.-facw(i,k))*(a0w(2)+a1w(2)*rewi(i,k)))
           tau2 = cip(i,k)*(   faci(i,k) *(a0i(1)+a1i(1)*reii(i,k))            &
                +(1.-faci(i,k))*(a0i(2)+a1i(2)*reii(i,k)))
           taucl(i,k) = tau1 + tau2
           ssa1 = 1.0 - (   facw(i,k) *(b0w(1)+b1w(1)*rew(i,k))                &
                + (1.-facw(i,k))*(b0w(2)+b1w(2)*rew(i,k)) )
           ssa2 = 1.0 - (  faci(i,k) *(b0i(1)                                  &
                + (b1i(1)+b2i(1)*rei(i,k))*rei(i,k))                           &
                + (1.-faci(i,k))*(b0i(2)                                       &
                + (b1i(2)+b2i(2)*rei(i,k))*rei(i,k)) )
           ssaw1 = ssa1 * tau1
           ssaw2 = ssa2 * tau2
           ssacw(i,k) = ssaw1 + ssaw2
           asy1 =     facw(i,k) *(c0w(1)+c1w(1)*rew(i,k))                      &
                + (1.-facw(i,k))*(c0w(2)+c1w(2)*rew(i,k))
           asy2 =     faci(i,k) *(c0i(1)                                       &
                + (c1i(1)+c2i(1)*rei(i,k))*rei(i,k) )                          &
                + (1.-faci(i,k))*(c0i(2)                                       &
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
         end if
       enddo
     enddo
   end if
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
!===> ... get aerosols and rayleigh scattering optical properties
!
     call rad_sw_aeros_tau_nasa(ipts,iv,kaer,kprf,idxc,cmix,denn,rh            &
                                ,tauaer,ssaaer,asyaer,taurs)
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
     call rad_sw_flux_nasa(ipts,tauto,ssato,asyto,fffto,csm,zth,albb,albd,     &
                 daytm,upflux,dwflux,dwsfxb,dwsfxd)
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
           if (ccly(i,k) .gt. 0.0) then
             tauto(i,k) = taucl(i,k) + tauto(i,k)
             ssat1(i,k) = ssacw(i,k) + ssat1(i,k)
             ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
             tem        = 1.0  / max(fpmin, ssat1(i,k))
             asyto(i,k) = (asycw(i,k) + asyt1(i,k)) * tem
             fffto(i,k) = (fffcw(i,k) + ffft1(i,k)) * tem
           end if
         enddo
       enddo
!
!===> ... cloudy sky fluxes calculations
!
       call rad_sw_flux_nasa(ipts,tauto,ssato,asyto,fffto,csm,zth,albb,albd,   &
                   daytm,upflux,dwflux,dwsfxb,dwsfxd)
!
       do k = 1,lp1
         do i = 1,ipts
           fnetc(i,k) = fnetc(i,k)+(dwflux(i,k)-upflux(i,k))*ss(iv)
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
   end subroutine rad_uv_nasa
!
!-------------------------------------------------------------------------------
   subroutine co2_flux_nasa(ipts,swc,u1,du,nu,swh,w1,dw,nw,tbl,daytm,dflx)
!-------------------------------------------------------------------------------
!
!  compute the absorption due to co2. ref: chou (j. climate, 1990,              
!     209-217)                                                                  
!     updated sep. 1999 based on nasa/tm-1999-104606, vol 15.
!  the effect of co2 absorption below the cloud top is neglected.               
!  input variables:                                                             
!     swc,swh     : column amount of co2 and water vapor                        
!     u1,du,w1,dw : coefficients
!     tbl         : look up co2 absorption table
!     nu,nw       : table dimensions
!     daytm       : daytime flag                                                
!  output variables:                                                            
!     dflx        : flux reduction due to co2 for clear sky                     
!                                                                               
!-------------------------------------------------------------------------------
   use rdparm99
!-------------------------------------------------------------------------------
   real                 ::  swc(imbx,lp1),   swh(imbx,lp1)
   real                 ::  dflx(imbx,lp1),  tbl(nu,nw)
   logical              ::  daytm(imax)
!
! ... table look-up for the reduction of clear-sky solar
!
   x1 = u1 - 0.5*du
   y1 = w1 - 0.5*dw
   do k = 2,lp1
     do i = 1,ipts
       if (daytm(i)) then
         clog = swc(i,k)
         wlog = swh(i,k)
         ic = int( (clog - x1)/du + 1.0e0)
         iw = int( (wlog - y1)/dw + 1.0e0)
         ic = max(2, min(nu, ic))
         iw = max(2, min(nw, iw))
         ic1 = ic - 1
         iw1 = iw - 1
         dc = clog - float(ic-2)*du - u1
         dd = wlog - float(iw-2)*dw - w1
         x2 = tbl(ic1,iw1) + (tbl(ic1,iw)-tbl(ic1,iw1))/dw * dd
         y2 = x2 + (tbl(ic,iw1) - tbl(ic1,iw1))/du * dc
!
         dflx(i,k) = dflx(i,k) + y2
       end if
     end do
   end do
!                                                                               
   return                                                                    
   end subroutine co2_flux_nasa
!-------------------------------------------------------------------------------
   subroutine rad_sw_aeros_tau_nasa(ipts,ib,kaer,kprf,idxc,cmix,denn,rh        &
                                   ,tau,ssa,asy,taurs)
!-------------------------------------------------------------------------------
!
! abstract:
!  compute aerosols optical properties in eight uv+vis bands and
!  four nir bands. there are seven different vertical prifile
!  structures. in the troposphere, aerosol distribution at each
!  grid point is composed from up to six components out of a total
!  of ten different substances.
!
! program history log:
!   1999-10-13  y.h.  updated to opac data (1998)
!
!   band: 1. 0.175-0.225 (uv-c)     2. 0.225-0.245;0.260-0.280 (uv-c)
!         3. 0.245-0.260 (uv-c)     4. 0.280-0.295 (uv-b)
!         5. 0.295-0.310 (uv-b)     6. 0.310-0.320 (uv-b)
!         7. 0.320-0.400 (uv-a)     8. 0.400-0.700 (par)
!         9. 2.27 - 4.0  (nir)     10. 1.22 - 2.27 (nir)
!        11. 0.70 - 1.22 (nir)     12. 0.70 - 4.0  (nir)
!
!  input parameters:
!     ib   - spectral band index                   -    1
!     kaer - =0 do not compute aerosols            -    1
!            =1 compute aerosol profiles
!     kprf - indecies of aerosol prof structures   -    imax
!     idxc - indecies of aerosol components        -    nxc*imax
!     cmix - mixing ratioes of aerosol components  -    ncx*imax
!     denn - aerosol number densities              -    ndn*imax
!     rh   - relative humidity in fraction         -    imax*l
!
!  output parameters:
!     tau  - optical depth                         n/d
!     ssa  - single scattering albedo              n/d
!     asy  - asymmetry parameter                   n/d
!     taurs- rayleigh scattering optical depth     n/d
!
!  variables in common block:
!     haer - scale height of aerosols              km   ndm*nae
!     idm  - aerosol domain index                  -    l*nae
!     dz   - layer thickness                       km   l
!     hz   - level height                          km   l+1
!     taur - rayleigh scattering optical depth     -    l*nbd
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm99
   use aerparm
   use comswaer99
!-------------------------------------------------------------------------------
!
! --- input
   real                 ::  cmix(nxc,imax), denn(ndn,imax), rh(imbx,l)
   integer              ::  idxc(nxc,imax)
   integer              ::  kprf(imax)
! --- output
   real                 ::  tau(imbx,l),  ssa(imbx,l),  asy(imbx,l),  taurs(l)
   real,save            ::  crt1,crt2
   data  crt1,crt2 / 30.0, 0.03333 /
!
!===> ... layer optical depth due to rayleigh scattering
!
   do k = 1,l
     taurs(k) = taur(k,ib)
   end do
!
   do k = 1,l
     do i = 1,ipts
       ssa(i,k) = 0.0e0
       asy(i,k) = 0.0e0
       tau(i,k) = 0.0e0
     end do
   end do
!
   if (kaer .lt. 1) return
!
   do i = 1,ipts
!
     kpf = kprf(i)
     do k = 1,l
       idom = idm(k,kpf)
       drh = rh(i,k) - 0.5
!
       if (idom .eq. 1) then
!
! --- 1st domain - mixing layer
!
         ext1 = 0.0
         sca1 = 0.0
         ssa1 = 0.0
         asf1 = 0.0
         do icmp = 1,nxc
           ic = idxc(icmp,i)
           if (ic .gt. ncm1) then
             ic1 = ic - ncm1
             drh1 = exp(abpw(ic1,ib)*drh)
             drh2 = drh * drh
             ex00 = aext(1,ic1,ib) + aext(2,ic1,ib)*drh                        &
                   + aext(3,ic1,ib)*drh1
             sc00 = bsca(1,ic1,ib) + bsca(2,ic1,ib)*drh                        &
                   + bsca(3,ic1,ib)*drh1
             ss00 = cssa(1,ic1,ib) + cssa(2,ic1,ib)*drh                        &
                   + cssa(3,ic1,ib)*drh2
             as00 = dasf(1,ic1,ib) + dasf(2,ic1,ib)*drh                        &
                   + dasf(3,ic1,ib)*drh2
           else if (ic .gt. 0) then
             ex00 = ext0(ic,ib)
             sc00 = sca0(ic,ib)
             ss00 = ssa0(ic,ib)
             as00 = asf0(ic,ib)
           else
             ex00 = 0.0
             sc00 = 0.0
             ss00 = 0.0
             as00 = 0.0
           end if
           ext1 = ext1 + cmix(icmp,i) * ex00
           sca1 = sca1 + cmix(icmp,i) * sc00
           ssa1 = ssa1 + cmix(icmp,i) * ss00 * ex00
           asf1 = asf1 + cmix(icmp,i) * as00 * sc00
         end do
         ext2 = ext1 * denn(1,i)
         ssa2 = ssa1 / ext1
         asf2 = asf1 / sca1
       else if (idom .eq. 2) then
!
! --- 2nd domain - mineral transport layers
!
         ext2 = ext0(6,ib) * denn(2,i)
         ssa2 = ssa0(6,ib)
         asf2 = asf0(6,ib)
       else if (idom .eq. 3) then
!
! --- 3rd domain - free tropospheric layers
!   1:inso 0.17E-3; 2:soot 0.4; 7:wasO 0.59983; n:730
!
         drh1 = exp(abpw(1,ib)*drh)
         drh2 = drh * drh
         ex01 = ext0(1,ib)
         sc01 = sca0(1,ib)
         ss01 = ssa0(1,ib)
         as01 = asf0(1,ib)
         ex02 = ext0(2,ib)
         sc02 = sca0(2,ib)
         ss02 = ssa0(2,ib)
         as02 = asf0(2,ib)
         ex03 = aext(1,1,ib) + aext(2,1,ib)*drh + aext(3,1,ib)*drh1
         sc03 = bsca(1,1,ib) + bsca(2,1,ib)*drh + bsca(3,1,ib)*drh1
         ss03 = cssa(1,1,ib) + cssa(2,1,ib)*drh + cssa(3,1,ib)*drh2
         as03 = dasf(1,1,ib) + dasf(2,1,ib)*drh + dasf(3,1,ib)*drh2
         ext1 = 0.17E-3*ex01 + 0.4*ex02 + 0.59983*ex03
         sca1 = 0.17E-3*sc01 + 0.4*sc02 + 0.59983*sc03
         ssa1 = 0.17E-3*ss01*ex01 + 0.4*ss02*ex02 + 0.59983*ss03*ex03
         asf1 = 0.17E-3*as01*sc01 + 0.4*as02*sc02 + 0.59983*as03*sc03
         ext2 = ext1 * 730.0
         ssa2 = ssa1 / ext1
         asf2 = asf1 / sca1
       else if (idom .eq. 4) then
!
! --- 4th domain - stratospheric layers
!
         ext2 = estr(ib)
         ssa2 = 0.9
         asf2 = 0.6
       else
!
! --- upper stratosphere assume no aerosols
!
         ext2 = 0.0
         ssa2 = 1.0
         asf2 = 0.0
       end if
!
       hd = haer(idom,kpf)
       if (hd .gt. 0.0e0) then
            hd1 = 1.0 / hd
            sig0u = exp(-hz(k)  *hd1)
            sig0l = exp(-hz(k+1)*hd1)
            tau(i,k) = ext2 * hd*(sig0l - sig0u)
       else
            tau(i,k) = ext2 * dz(k)
       end if
       ssa(i,k) = ssa2
       asy(i,k) = asf2
     end do
   end do
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
       end if
     end do
   end do
!
   return
   end subroutine rad_sw_aeros_tau_nasa
!-------------------------------------------------------------------------------
   subroutine rad_sw_flux_nasa(ipts,tau,ssc,g0,ff,csm,zth,alb,ald,             &
                     daytm,upflux,dwflux,dwsfcb,dwsfcd)
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
!    ff:  the effective forward scattering factor
!    csm: secant of the zenith angle
!    zth: cosine of the zenith angle
!    alb: surface albedo for direct radiation
!    ald: surface albedo for diffused radiation
!    daytm: daytime flag
!
!  outputs:
!    upflux: upward fluxes
!    dwflux: downward fluxes
!    dwsfcb: downward surface flux direct component
!    dwsfcd: downward surface flux diffused component
!
!-------------------------------------------------------------------------------
   use rdparm99
!-------------------------------------------------------------------------------
!
! --- input
!
   real                 ::  tau(imbx,l), ssc(imbx,l), g0(imbx,l), ff(imbx,l)
   real                 ::  csm(imax),   zth(imax),   alb(imax),  ald(imax)
   logical              ::  daytm(imax)
!
! --- output
!
   real                 ::  upflux(imbx,lp1),dwflux(imbx,lp1)
   real                 ::  dwsfcb(imax),dwsfcd(imax)
!
! --- temporary
!
   real                 ::  ttb(imbx,lp1),tdn(imbx,lp1),rup(imbx,lp1)
   real                 ::  rfu(imbx,lp1),rfd(imbx,lp1),tb (imbx,lp1)
   real                 ::  tt (imbx,lp1,2),rr (imbx,lp1,2)
!
!===> ... diffuse incident radiation is approximated by beam radiation
!         with an incident angle of 53 degrees. cos(53) = 0.602
!
   zthd = 0.602
   csmd = 1.0 / zthd
   epsln = 1.0e-30
!
!===> ... delta-eddington scaling of single scattering albedo,
!         optical thickness, and asymmetry factor, k & h eqs(27-29)
!
   do k = 1,l
     do i = 1,ipts
       if (daytm(i)) then
!
!===> ... delta-eddington scaling of single scattering albedo,
!         optical thickness, and asymmetry factor, k & h eqs(27-29)
!
         aa   = 1.0e0 - ff(i,k)*ssc(i,k)
         taup = tau(i,k) * aa
         sscp = ssc(i,k) * (1.0e0 - ff(i,k)) / aa
         gp   = (g0(i,k) - ff(i,k)) / (1.0e0 - ff(i,k))
!
         oms1 = 1.0e0 - sscp
         ogs1 = 1.0e0 - sscp*gp
         tlam = 3.0e0 * oms1*ogs1
         slam = sqrt(tlam)

         u1   = 1.50e0 * ogs1 / slam
         u1p1 = u1 + 1.0e0
         u1m1 = u1 - 1.0e0
         e1   = exp(max(-taup*slam, -30.0))
         u1e  = u1 * e1
         u1epe= u1e + e1
         u1eme= u1e - e1
         den  = 1.0e0 / ((u1p1 + u1eme)*(u1p1 - u1eme))
         rf1  = (u1p1 + u1epe) * (u1m1 - u1eme) * den
         tf1  = 4.0e0 * u1e * den
!
!===> ... compute layer transmissions and reflections
!         (i,k,j) j=1,2 for layer k illuminated by diffuse and
!                       direct incoming radiation
!         rr   :  layer reflection
!         tt   :  layer total transmission
!         tb   :  layer direc transmission
!
!       diffuse radiation
!       -----------------
!
         zzth = zthd
         zz   = zzth * zzth
         den1 = 1.0 - tlam*zz
!
!===> ... safety check
!
         if (abs(den1) .lt. 1.0e-8) then
           zzth = zzth + 0.001
           zz   = zzth * zzth
           den1 = 1.0 - tlam*zz
         end if
         den1 = sscp / den1
!
         gama = 0.50e0 * (1.0e0 + 3.0e0*gp*oms1*zz) * den1
         alfa = 0.75e0 * zthd * (gp + ogs1) * den1
         amg  = alfa - gama
         apg  = alfa + gama
!
         tb(i,k)  = exp( -min(30.0, taup*csmd) )
         za       = amg * tb(i,k)
         rr(i,k,1)= za*tf1 + apg*rf1 - amg
         tt(i,k,1)= za*rf1 + apg*tf1 + (1.0e0-apg)*tb(i,k)
!
!       direct radiation
!       -----------------
!
         zzth = zth(i)
         zz   = zzth * zzth
         den1 = 1.0 - tlam*zz
!
!===> ... safety check
!
         if (abs(den1) .lt. 1.0e-8) then
           zzth = zzth + 0.001
           zz   = zzth * zzth
           den1 = 1.0 - tlam*zz
         end if
         den1   = sscp / den1
!
         gama = 0.50e0 * (1.0e0 + 3.0e0*gp*oms1*zz) * den1
         alfa = 0.75e0 * zth(i) * (gp + ogs1) * den1
         amg  = alfa - gama
         apg  = alfa + gama
!
         tb(i,k)  = exp( -min(30.0, taup*csm(i)) )
         za       = amg * tb(i,k)
         rr(i,k,2)= za*tf1 + apg*rf1 - amg
         tt(i,k,2)= za*rf1 + apg*tf1 + (1.0e0-apg)*tb(i,k)
!
         tb(i,k)   = max(0.0e0, tb(i,k))
         rr(i,k,2) = max(0.0e0, rr(i,k,2))
         tt(i,k,2) = max(0.0e0, tt(i,k,2))
         rr(i,k,1) = max(0.0e0, rr(i,k,1))
         tt(i,k,1) = max(0.0e0, tt(i,k,1))
!
! --- night time condition
!
       else
         tb(i,k)   = 0.0e0
         rr(i,k,1) = 0.0e0
         rr(i,k,2) = 0.0e0
         tt(i,k,1) = 0.0e0
         tt(i,k,2) = 0.0e0
       end if
     end do
   end do
!
! --- at the surface
!
   do i = 1,ipts
     tb(i,lp1)   = 0.0e0
     rr(i,lp1,2) = alb(i)
     tt(i,lp1,2) = 0.0e0
     rr(i,lp1,1) = ald(i)
     tt(i,lp1,1) = 0.0e0
   end do
!
   do i = 1,ipts
     ttb(i,1) = tb(i,1)
     tdn(i,1) = tt(i,1,2)
     rfd(i,1) = rr(i,1,1)
   end do
!
!===> ... layers added downward starting from top
!
   do k = 2,lp1
     do i = 1,ipts
       if (daytm(i)) then
         den = tt(i,k,1) / (1.0e0 - rfd(i,k-1) * rr(i,k,1))
         ttb(i,k) = ttb(i,k-1) * tb(i,k)
         if (ttb(i,k) .lt. epsln) ttb(i,k) = 0.0
         tdn(i,k) = ttb(i,k-1)*tt(i,k,2)+(tdn(i,k-1)-ttb(i,k-1)                &
                  + ttb(i,k-1)*rr(i,k,2)*rfd(i,k-1)) * den
         rfd(i,k) = rr(i,k,1) + tt(i,k,1)*rfd(i,k-1) * den
       end if
     end do
   end do
!
!===> ... layers added upward starting from surface
!
   do i = 1,ipts
     rfu(i,lp1) = rr(i,lp1,1)
     rup(i,lp1) = rr(i,lp1,2)
   end do
   do k = l,1,-1
     kp1 = k + 1
     do i = 1,ipts
       if (daytm(i)) then
         den = tt(i,k,1) / (1.0e0 - rfu(i,kp1) * rr(i,k,1))
         rup(i,k) = rr(i,k,2) + ((tt(i,k,2)-tb(i,k))*rfu(i,kp1)                &
                  + tb(i,k)*rup(i,kp1)) * den
         rfu(i,k) = rr(i,k,1) + tt(i,k,1)*rfu(i,kp1) * den
       end if
     end do
   end do
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
     end if
   end do
   do k = 2,lp1
     km1 = k - 1
     do i = 1,ipts
       if (daytm(i)) then
         den = 1.0e0 / (1.0e0 - rfd(i,km1)*rfu(i,k))
         aa  = ttb(i,km1) * rup(i,k)
         bb  = tdn(i,km1) - ttb(i,km1)
         upflux(i,k) = (aa + bb*rfu(i,k)) * den
         dwflux(i,k) = ttb(i,km1) + (aa*rfd(i,km1) + bb) * den
       else
         upflux(i,k) = 0.0e0
         dwflux(i,k) = 0.0e0
       end if
     end do
   end do
!
!===> ... surface downward fluxes
!
   do i = 1,ipts
     dwsfcb(i) = ttb(i,l)
     dwsfcd(i) = dwflux(i,lp1)-dwsfcb(i)
   end do
!
   return
   end subroutine rad_sw_flux_nasa
#endif /* ~SWRMDC end */
