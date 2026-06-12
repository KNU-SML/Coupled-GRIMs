#include <define.h>
   subroutine lsm_vic1_driver(msl,    nsl,    msub,     md,                    &
                    dtime,  month,       i,    lat,                            &
                      prc,   pgcm,    wind,   tgcm,                            &
                     qgcm,   zgcm,   flwds,   sols,                            &
                     binf,     ds,     dsm,     ws,                            &
                     expt,    kst,     dph,    cef,                            &
                      bub,    qrt,     bkd,    sld,                            &
                      wcr,    wpw,    silz,   snwz,                            &
                      smr,    smx,   dphnd,  smxnd,                            &
                    expnd,  bubnd,   alpnd,  betnd,                            &
                    gamnd,   nveg,     nvt,     wt,                            &
                       rt,   flai,     cwt,    csn,                            &
                      smc,    sic,     tnd,    swq,                            &
                      rsn,    tsf,     tpk,    sfw,                            &
                      pkw,    fmu,   sfall,  lstsn,                            &
                       sh,     lh,      gh,     ts,                            &
                     albd,  ovflw,   bsflw, snowmt,                            &
#ifdef SMP_RA2SFC
                   snowev,    gdh,   ez0,  netsw)
#else
                   snowev,    gdh,   ez0)
#endif
!-------------------------------------------------------------------------------
!
! subroutine: lsm_vic1_driver
!
! abstract: simulate the interaction between the land surface
!           and the atmosphere over a cell
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2003-11-01  ji chen                modified from 4.0.3 uw vic c version
!   2008-08-01  kyeonghee seol         debugged, scm option
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use vic_veglib
#include <vartyp.h>
!-------------------------------------------------------------------------------
!
! ------------------- input variables ----------------------------------
!
! -- model basic parameter
!
   integer  ::  msl              ! number of soil layers
   integer  ::  nsl              ! number of soil nodes for soil temperature
   integer  ::  msub             ! maximum number of subgrid points
   integer  ::  md               ! mdist (1no/2with) precipiation distribution
   real     ::  dtime            ! vic time step (second)
   integer  ::  month            ! current month
   integer  ::  i, lat           ! long and lati indexes
!
! -- atmosphere data
!
   real     ::  prc              ! precipitation rate (m h2o/time step)
   real     ::  pgcm             ! atm bottom level pressure (pa)
   real     ::  wind             ! atm bottom level wind (m/s)
   real     ::  tgcm             ! atm bottom level temperature (kelvin)
   real     ::  qgcm             ! atm bottom level specific humidity (kg/kg)
   real     ::  zgcm             ! atm bottom level height above surface (m)
   real     ::  flwds            ! downward longwave radi onto surface (w/m2)
   real     ::  sols             ! solar rad onto srf (w/m2)
#ifdef SMP_RA2SFC
   real     ::  netsw            ! net solar rad onto srf (w/m2)
#endif
!
! -- soil parameters
!
   real     ::  binf             ! variable infiltration curve parameter (n/a)
   real     ::  ds               ! fract of dsm nonlinear baseflow begins (fract)
   real     ::  dsm              ! maximum velocity of baseflow (mm/day)
   real     ::  ws               ! fract maxi sm nonlinear baseflow occurs (fract)
   real     ::  cef              ! exponent used in infiltration curve (n/a)
   real     ::  expt(msl)        ! para the vari of ksat with sm (n/a)
   real     ::  kst(msl)         ! saturated hydrologic conductivity (mm/day)
   real     ::  dph(msl)         ! thickness of soil layer (m)
   real     ::  bub(msl)         ! bubbling pressure of soil layer (cm)
   real     ::  qrt(msl)         ! quartz content of soil layer (fraction)
   real     ::  bkd(msl)         ! bulk density of soil layer (kg/m3)
   real     ::  sld(msl)         ! soil density of soil layer (kg/m3)
   real     ::  wcr(msl)         ! fract sm content at the critical point (mm)
   real     ::  wpw(msl)         ! fractional sm content wilting point (mm)
   real     ::  silz             ! surface roughness of bare soil (m)
   real     ::  snwz             ! surface roughness of snowpack (m)
   real     ::  smr(msl)         ! soil moisture residual moisture (mm)
   real     ::  smx(msl)         ! maximum soil moisture (mm)
   real     ::  dphnd(nsl)       ! thickness of soil node (m)
   real     ::  smxnd(nsl)       ! maximum soil moisture at soil node (m3/m3)
   real     ::  expnd(nsl)       ! para the vari of ksat at soil node (n/a)
   real     ::  bubnd(nsl)       ! bubbling pressure at soil node (cm)
   real     ::  alpnd(nsl)       ! para alpha at soil node
   real     ::  betnd(nsl)       ! para beta at soil node
   real     ::  gamnd(nsl)       ! para gamma at soil node
   real     ::  vz0              ! roughness length for vegetation
!
! -- vegetation parameters
!
   integer  ::  nveg             ! number of vegetation type in a grid cell
   integer  ::  nvt(msub)        ! vegetation type number
   real     ::  wt(msub)         ! fraction of grid cell covered by veg
   real     ::  rt(msl,msub)     ! fraction of root in the soil layer
   real     ::  flai(msub)       ! leaf area index
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
!
   real     ::  cwt(md,msub)     ! canopy water (m h2o)
   real     ::  csn(md,msub)     ! canopy snow  (m h2o)
   real     ::  smc(msl,md,msub) ! soil moisture content (mm)
   real     ::  sic(msl,md,msub) ! soil ice content (mm)
   real     ::  tnd(nsl,md,msub) ! temperature at soil nodes (k)
   real     ::  swq(md,msub)     ! snow water equivalent (m h2o)
   real     ::  rsn(md,msub)     ! snow density (kg/m3)
   real     ::  tsf(md,msub)     ! snow surface temperature (k)
   real     ::  tpk(md,msub)     ! snow pack temperature (k)
   real     ::  sfw(md,msub)     ! surface snow water equivalent (m h2o)
   real     ::  pkw(md,msub)     ! snow pack snow water equivalent (m h2o)
   real     ::  fmu              ! precipitation fraction (fraction)
   integer  ::  lstsn            ! number of time step since last snowevant
                                 ! > 0: snow covered ground, < 0: no snow
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
!
   real     ::  sh               ! sensible heat flux (w/m**2) [+ to ground]
   real     ::  lh               ! latent heat flux (w/m**2) [+ to ground]
   real     ::  gh               ! ground heat flux (w/m**2) [+ to ground]
   real     ::  ts               ! surface radiative temperature (kelvin)
   real     ::  albd             ! albedo (fraction)
   real     ::  ovflw            ! overland flow (mm/time step)
   real     ::  bsflw            ! base flow (mm/time step)
   real     ::  snowmt           ! snow melt (m/s)
   real     ::  snowev           ! evaporation over snow surface (w/m2)
   real     ::  sfall            ! snow fall (m h2o/time step)
   real     ::  gdh              ! ground heat storage (w/m**2) [+ to ground]
   real     ::  ez0              ! roughness length for vegetation
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
!
   integer  ::  vtype         ! vegetation type
   integer  ::  ndist         ! number of precipitation districts
!
   real     ::  rhoair           ! air density (kg/m3)
   real     ::  vpair            ! vapor pressure (pa)
   real     ::  vpd              ! vapor pressure deficit (pa)
   real     ::  svp              ! function name to compute vapor pressure
!
   real     ::  rfall            ! rain fall (m h2o/time step)
   real     ::  ppt              ! effective precipitation (m/time step)
!
   real     ::  u(3), ra(3)      ! unified wind and resistance for vic
!
   real     ::  fflai            ! fractional subgrid lai
   real     ::  fwt              ! fractional subgrid weight
   real     ::  fswq             ! fractional subgrid swq
   real     ::  frsn             ! fractional subgrid rsn
   real     ::  ftsf             ! fractional subgrid tsf
   real     ::  ftpk             ! fractional subgrid tpk
   real     ::  fsfw             ! fractional subgrid sfw
   real     ::  fpkw             ! fractional subgrid pkw
   real     ::  fcwt             ! fractional subgrid cwt
   real     ::  fcsn             ! fractional subgrid csn
   real     ::  fsmc(msl)        ! fractional subgrid smc
   real     ::  fsic(msl)        ! fractional subgrid sic
   real     ::  frt(msl)         ! fractional subgrid rt
   real     ::  ftnd(nsl)        ! fractional subgrid tnd
   real     ::  falbd            ! fractional albedo over subgrid
   real     ::  flh              ! fractional subgrid latent heat
   real     ::  fsh              ! fractional subgrid sensible heat
   real     ::  fgh              ! fractional subgrid ground heat
!
   real     ::  fsnowm           ! fractional subgrid snow melt (m/time step)
   real     ::  snowme           ! fractional subgrid snow melt (w/m2)
   real     ::  fsnowe           ! fractional subgrid snow evaporation
   real     ::  frtf             ! fractional subgrid rtf
   real     ::  fovflw           ! fractional subgrid overland flow
   real     ::  fbsflw           ! fractional subgrid base flow
   real     ::  fgdh             ! fractional subgrid ground heat storage
   real     ::  fz0              ! fractional subgrid roughness length
!
   real     ::  sncvfr           ! snow cover fraction (fraction)
   real     ::  tmpt             ! temporary soil temperature (k)
   real     ::  sfc_max_unfwat   ! function name
!
!     t0c           ::   ice/water mix temperature (k)
!
   real, parameter  ::  rd = 2.8705e+2, t0c=273.15
!
   integer  ::  m, n, nv, idist, k     ! loop indexes
!-------------------------------------------------------------------------------
!
!----------- following is land surface scheme from vic -----------------
!
#ifdef DBG
   print *,'--- enter lsm_vic1_driver --- lat=',lat,' lon i=',i
#endif
!   
   sh = 0.
   lh = 0.
   gh = 0.
   ts = 0.
   albd = 0.
   ovflw = 0.
   bsflw = 0.
   snowmt = 0.
   snowev = 0.
   gdh = 0.
   ez0 = 0.
!
! ----------------------------------------------------------------------
! compute air density, vapor pressure
! ----------------------------------------------------------------------
!
   rhoair = pgcm / (rd * tgcm)       ! air density (kg/m3)
   vpair = pgcm*qgcm/0.622
   vpd = svp(tgcm) - vpair
   if(vpd.lt.0.0) vpd = 0.0
!
! ----------------------------------------------------------------------
! depended on the air temperature,
! the precipitation is splitted into rainfall and/or snowfall.
! ----------------------------------------------------------------------
!
#ifdef DBG
   print *,'vic check in lsm_vic1_driver before vic_rainorsnow'
#endif
!
   rfall = 0.0
   sfall = 0.0
!
   call vic_rainorsnow(msl,   nsl,  msub,    md,                               &
                  nveg,   prc,  tgcm,   fmu,                                   &
                 lstsn, rfall, sfall,   cwt,                                   &
                   csn,   swq,   rsn,   tsf,                                   &
                   tpk,   sfw,   pkw,   smc,                                   &
                   sic,   tnd)
!
! ----------------------------------------------------------------------
! start land surface process simulation
! ----------------------------------------------------------------------
!
   if(fmu.ne.1.) then
     ndist = md
   else 
     ndist = 1
   end if
!
   do nv = 1,nveg     ! nveg, number of different land covers in a cell
!
     vtype = nvt(nv)           ! vegetation type
#ifdef SMP_RA2SFC
!
! replace LAI from RA2 with LAI from vic_veglib 
!
     fflai = veg_lai(month,vtype)          ! leaf area index
#else
     fflai = flai(nv)          ! leaf area index
#endif
#ifdef DBG
!
     print *,'lsm_vic1_driver vty=', vtype, fflai
#endif
!
     do m = 1,msl
       frt(m) = rt(m,nv)      ! root content
     end do
!
     do idist = 1,ndist
       if(idist .eq. 1) then           ! in the wet district
         fwt = fmu*wt(nv)
       else                            ! in the dry district
         fwt = (1.0 - fmu)*wt(nv)
       end if
!
       sncvfr = 0.0                    ! initialize snowcoverfraction
       fswq = swq(idist,nv)
       frsn = rsn(idist,nv)
       ftsf = tsf(idist,nv)
       ftpk = tpk(idist,nv)
       fsfw = sfw(idist,nv)
       fpkw = pkw(idist,nv)
       fcwt = cwt(idist,nv)
       fcsn = csn(idist,nv)
!
       do m = 1,msl
         fsmc(m) = smc(m,idist,nv)
         fsic(m) = sic(m,idist,nv)
       end do
!
       do n = 1,nsl
         ftnd(n) = tnd(n,idist,nv)
       end do
!
! ----------------------------------------------------------------------
! compute the vic_exch resistance
! ----------------------------------------------------------------------
!
#ifdef DBG
       print *,'vic check in lsm_vic1_driver before vic_exch'
#endif
       call vic_exch(vtype,  month, silz, snwz, vz0,                           &
                      zgcm,   wind,    u,   ra)
!
! ----------------------------------------------------------------------
! update the soil ice content according to the new soil condition
! ----------------------------------------------------------------------
!
       do m = 1,msl
         if(m.eq.1) then
           tmpt = (ftnd(1)+ftnd(2))/2.0
         else
           tmpt = ftnd(m+1)
         endif
         if(tmpt.lt.t0c) then
           fsic(m) = fsmc(m)-                                                  &
                sfc_max_unfwat(tmpt,smx(m),bub(m),expt(m))
           if(fsic(m).lt.0.0) fsic(m) = 0.0
         else
           fsic(m) = 0.0
         endif
       enddo
!
! ----------------------------------------------------------------------
! solve land surface fluxes
! ----------------------------------------------------------------------
!
#ifdef DBG
       print *,'vic check in lsm_vic1_driver before surface_fluxes'
       print *,'vtype ',vtype,' month=',month,fsmc,fsic
#endif
!
       call vic_surface_flux(msl,    nsl,  dtime,  month,                      &
                   rfall,  sfall,   pgcm,   tgcm,                              &
                   flwds,   sols,  vpair,    vpd,                              &
                  rhoair,      u,     ra,   binf,                              &
                     dph,    qrt,    bkd,    sld,                              &
                     wcr,    wpw,    smr,  dphnd,                              &
                   smxnd,  expnd,  bubnd,  alpnd,                              &
                   betnd,  gamnd,  vtype,    frt,                              &
                   fflai,  lstsn,   fswq,   frsn,                              &
                    ftsf,   ftpk,   fsfw,   fpkw,                              &
                    fcwt,   fcsn,   fsmc,   fsic,                              &
                    ftnd,  falbd,    flh,    fsh,                              &
                     fgh,   frtf, fsnowm, fsnowe,                              &
#ifdef SMP_RA2SFC
                    fgdh,  netsw, sncvfr)
#else
                    fgdh, sncvfr)
#endif
!
#ifdef DBG
       print *,'in lsm_vic1_driver after vic_surface_flux',fsmc,fsic
#endif
!
! ----------------------------------------------------------------------
! update the soil ice content after computing new soil temperature
! ----------------------------------------------------------------------
!
       do m = 1,msl
         if(m.eq.1) then
           tmpt = (ftnd(1)+ftnd(2))/2.0
         else
           tmpt = ftnd(m+1)
         endif
         if(tmpt.lt.t0c) then
           fsic(m) = fsmc(m)-                                                  &
                sfc_max_unfwat(tmpt,smx(m),bub(m),expt(m))
           if(fsic(m).lt.0.0) fsic(m) = 0.0
         else
           fsic(m) = 0.0
         endif
       enddo
!
! -- update variables for continous run
!
       cwt(idist,nv) = fcwt
       csn(idist,nv) = fcsn
       swq(idist,nv) = fswq
       rsn(idist,nv) = frsn
       tsf(idist,nv) = ftsf
       tpk(idist,nv) = ftpk
       sfw(idist,nv) = fsfw
       pkw(idist,nv) = fpkw
!
       do n = 1,nsl
         tnd(n,idist,nv) = ftnd(n)
       end do
!
       if(md.ne.1) then
         if(fmu.eq.1) then
           cwt(md,nv) = fcwt
           csn(md,nv) = fcsn
           swq(md,nv) = fswq
           rsn(md,nv) = frsn
           tsf(md,nv) = ftsf
           tpk(md,nv) = ftpk
           sfw(md,nv) = fsfw
           pkw(md,nv) = fpkw
           do n = 1,nsl
             tnd(n,md,nv) = ftnd(n)
           end do
         end if
       end if
!
! -- aggragate output variables
!
       albd = albd + falbd*fwt
       lh = lh + flh*fwt
       sh = sh + fsh*fwt
       gh = gh + fgh*fwt
       snowmt = snowmt + fsnowm*fwt/dtime  ! m/time step-> m/s
       snowev = snowev + fsnowe*fwt
       gdh = gdh + fgdh*fwt
!
! roughness length
!
       if (vtype .eq. 12) then
         fz0=silz
       else
         fz0=vz0
       endif
       ez0 = ez0 + fwt*(sncvfr*snwz+(1-sncvfr)*fz0)
!
       ts = ts + fwt*(sncvfr*ftsf+(1.-sncvfr)*ftnd(1))
!
! ----------------------------------------------------------------------
! compute land surface runoff, and soil moisture movement
! ----------------------------------------------------------------------
!
#ifdef DBG
       print *,'vic check in lsm_vic1_driver before runoff', fsmc,fsic
       print *,' fsnowm, frtf, sncvfr ',fsnowm,frtf,sncvfr
       print *,' others ',msl,  dtime,    ppt,  binf,                          &
                    ws,     ds,    dsm,   smx,                                 &
                   smr,    kst,   expt,                                        &
                fovflw, fbsflw
#endif
!
       ppt = fsnowm+(1.0-sncvfr)*frtf
!
       call vic_runoff_solver(msl,  dtime,    ppt,   binf,                     &
                               ws,     ds,    dsm,    cef,                     &
                              smx,    smr,    kst,   expt,                     &
                              dph,    bub,                                     &
                             fsmc,   fsic, fovflw, fbsflw)
!
#ifdef DBG
       print *,' after runoff',fsmc,fsic,fovflw, fbsflw
#endif
!
! -- update variables for continous run
!
       do m = 1,msl
         smc(m,idist,nv) = fsmc(m)
         sic(m,idist,nv) = fsic(m)
!
         if (md.ne.1) then
           if (fmu.eq.1.0) then
             smc(m,md,nv) = fsmc(m)
             sic(m,md,nv) = fsic(m)
           end if
         end if
       end do
!
! -- aggragate output variables
!
       ovflw = ovflw + fovflw*fwt
       bsflw = bsflw + fbsflw*fwt
     end do ! idist
!
   end do ! nv
!
   return
   end subroutine lsm_vic1_driver
