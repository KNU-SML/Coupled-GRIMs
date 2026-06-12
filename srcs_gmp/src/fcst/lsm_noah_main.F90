#include <define.h>
!-------------------------------------------------------------------------------
   subroutine lsm_noah_main(ims2,imx2,kmx,                                     &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,                    &
                  dlwflx,slrad,snowmt,snowev,delt,z0cm,tg3,                    &
                  gflux,zsoil,cm,ch,rcl,                                       &
                  prsi1,prsl1,prsik1,prslk1,                                   &
                  zl1,slimsk,inistp,lat,                                       &
                  drain,evap,hflx, ep, wind,                                   &
                  snwdph, slc, snoalb, slptype, shdmin, shdmax,                &
#ifndef NOAHYDRO
                  snowfl, runoff, precip, srflag)
#else
                  snowfl, runoff, precip, srflag, ec)
#endif
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [lsm_noah_main] 
!        |
!        |--- [phys_lsm_noah] *
!                  |--- [lsm_noah_module]
!                           |--- [noah_read_parameter]  
!                           |--- [noah_fresh_snow]  
!                           |--- [noah_snow_fraction]  
!                           |--- [noah_soilt_conductivity]  
!                           |--- [noah_snow_rough]  
!                           |--- [noah_albedo]  
!                           |--- [lsm_exch_coeff]
!                           |--- [noah_penman]  
!                           |--- [noah_canopy_res]  
!                           |--- [noah_bare_soil_solver]  
!                           |--- [noah_snow_cover_solver]  
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2008-05-01  kyung-hee seol         implementation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!   2012-05-01  song-you hong          constant variables and cleanup
!
!-------------------------------------------------------------------------------
   use constant, only : cal_,cp=>cp_,g=>g_,hvap=>hvap_,rd=>rd_,rv=>rv_,        &
                        rdorv=>rdorv_,rdorvm1=>rdorvm1_,rvordm1=>rvordm1_,     &
                        sigma=>sbc_,cb2pa=>cb2pa_,rhoh2o=>rhoh2o_,             &
                        convrad=>convrad_,elocp=>elocp_
   use paramodel, only : ILOTS
#ifdef DFS
   use dfsvar, only   : iope
#else
   use comio, only    : iope
#endif
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
#ifdef CRAY_THREAD
!fpp$ noconcur r
!fpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
!
! ca is the von karman constant
!
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
! passing array
!
   real                 ::  prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real                 ::  zl1(imx2),t1(imx2),q1(imx2)
   real                 ::  snoweq(imx2),snowmt(imx2),snowev(imx2)
   real                 ::  cm(imx2),ch(imx2)
   real                 ::  tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
   real                 ::  smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
   real                 ::  z0cm(imx2),gflux(imx2)
   real                 ::  slimsk(imx2)
   real                 ::  drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
   real                 ::  evap(imx2),hflx(imx2),ep(imx2)
   real                 ::  wind(imx2)
   real                 ::  dlwflx(imx2)
   integer              ::  vegtype(imx2)
   integer              ::  soiltyp(imx2)
   integer              ::  slptype(imx2)
   real                 ::  runoff(imx2)
   real                 ::  precip(imx2), srflag(imx2)
   real                 ::  snowfl(imx2)
#ifdef NOAHYDRO
   real                 ::  sn_new(imx2)
   real                 ::  ec(imx2)
#endif
   real                 ::  snoalb(imx2),shdmin(imx2),shdmax(imx2)
   real                 ::  snwdph(imx2), slc(imx2, kmx)
!   real                 ::  ktsoil
!
! local array
!
   real                 ::  tsurf(ILOTS)
   real                 ::  q0(ILOTS)
   real                 ::  theta1(ILOTS),tv1(ILOTS), rho(ILOTS)
   real                 ::  qs1(ILOTS),qss(ILOTS)
   real                 ::  rch(ILOTS),slwd(ILOTS)
!
   integer              ::  couple
   integer              ::  ice(ILOTS)
   integer              ::  nroot(ILOTS)
   real                 ::  sldpth(ILOTS, kmx)
   real                 ::  zlvl(ILOTS), radflx(ILOTS)
   real                 ::  sfcprs(ILOTS), prcp(ILOTS), dqsdt2(ILOTS)
   real                 ::  snowh(ILOTS), sneqv(ILOTS), albedo(ILOTS)
   real                 ::  cmx(ILOTS),     chx(ILOTS),    cmc(ILOTS)
   real                 ::  runoff1(ILOTS),runoff2(ILOTS)
   real                 ::  runoff3(ILOTS)
   real                 ::  sldpth1d(kmx), smc1d(kmx)
   real                 ::  slc1d(kmx),    stc1d(kmx)
   real                 ::  et(kmx)
   real                 ::  ptu(ILOTS), alb(ILOTS)
   real                 ::  edir(ILOTS), ett(ILOTS)
#ifndef NOHYDRO
   real                 ::  ec(ILOTS)
#endif
   real                 ::  esnow(ILOTS), drip(ILOTS), dew(ILOTS)
   real                 ::  flx1(ILOTS), flx2(ILOTS), flx3(ILOTS)
   real                 ::  sncovr(ILOTS), rc(ILOTS), pc(ILOTS)
   real                 ::  rsmin(ILOTS), xlai(ILOTS), rcs(ILOTS)
   real                 ::  rct(ILOTS), rcq(ILOTS)
   real                 ::  rcsoil(ILOTS),soilw(ILOTS),soilm(ILOTS)
   real                 ::  smcwlt(ILOTS),smcdry(ILOTS)
   real                 ::  smcref(ILOTS)
   real                 ::  smcmax(ILOTS)
!
   real                 ::  delt2,delt,rcl
   integer              ::  lond,latd
   integer              ::  lat,inistp
   integer              ::  im,km,ims2,imx2,kmx
   integer              ::  i, k
#ifdef ICE
   real                 ::  fpvs
   external fpvs
#else
   real                 ::  fpvs0
   external fpvs0
#endif
!    
!  declare parameters needed for deriving noah variables
!  source: ETA's surface.f
!
   real,parameter       ::  a2=17.2693882,a3=273.16,a4=35.86,a23m4=a2*(a3-a4)
!------------------------------------------------------------------------------
!
!  surface energy/water balance over land and seaice
!
!  initialization
!
   latd = 23
   lond = 100
   delt2 = delt * 2.
   im = ims2
   km = kmx
   do i = 1,im                                           
     if(slimsk(i).eq.1.) then       
       if(srflag(i) .eq. 1.) then
         snowfl(i)=precip(i)
       else
         snowfl(i)=0.
       endif
     endif
   enddo
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  1000.*prsi is in pascals
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
!  surface
!  convert slrad to the civilized unit from langley minute-1 k-4
!  surface roughness length is converted to m from cm
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
       q0(i) = max(q1(i),1.e-8)
       tsurf(i) = tskin(i)
       theta1(i) = t1(i) / prslk1(i) * prsik1(i)
       tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
       rho(i) = (cb2pa * prsl1(i)) / (rd * tv1(i))
#ifdef ICE
       qs1(i) = cb2pa * fpvs(t1(i))
#else
       qs1(i) = cb2pa * fpvs0(t1(i))
#endif
       qs1(i) = rdorv * qs1(i) / (cb2pa * prsl1(i) +  rdorvm1 * qs1(i))
       qs1(i) = max(qs1(i), 1.e-8)
       q0(i) = min(qs1(i),q0(i))
#ifdef ICE
       qss(i) = cb2pa * fpvs(tsurf(i))
#else
       qss(i) = cb2pa * fpvs0(tsurf(i))
#endif
       qss(i) = rdorv * qss(i) / (cb2pa * prsi1(i) +  rdorvm1 * qss(i))
     endif
   enddo
!
!  rcp = rho cp ch v
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
       rch(i) = rho(i) * cp * ch(i) * wind(i)
     endif
   enddo
!
!      gmp/rmp -> noah: prepare variables to run noah lsm
!
! 1. configuration information (c):
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
!
!  1.1   ice        sea-ice flag  (=1: sea-ice, =0: land)
!
       if(slimsk(i).eq.2.) then       
         ice(i) =1 
       else 
         ice(i) =0 
       endif
!
!  1.2   dt         timestep (sec) (dt should not exceed 3600 secs)
!        dt = delt
!  1.3   zlvl       height (m) above ground of atmospheric forcing variables
!
       zlvl(i) = zl1(i)
!
!  1.4   nsoil      number of soil layers (at least 2)
!
!        nsoil = km 
!  1.5   sldpth     the thickness of each soil layer (m)
!
!  noah:zsoil(1) = -sldpth(1)
!
       sldpth(i,1) = -zsoil(i,1)
       do k = 2,km
!
!  noah: zsoil(k) = -sldpth(k)+zsoil(k-1)
!
         sldpth(i,k) = zsoil(i,k-1) - zsoil(i,k)
       enddo
     endif
   enddo
! 
! 2. forcing data (f):
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
!
!  2.1   lwdn       lw downward radiation (w m-2; positive, not net longwave)
!
!  2.2   soldn      solar downward radiation (w m-2; positive, not net solar)
!
!      lwdn = dlwflx
!
       slwd(i) = slrad(i) * convrad
       radflx(i) = -1. * slwd(i)
!
!  2.3   sfcprs     pressure at height zlvl above ground (pascals)
!
       sfcprs(i) = cb2pa * prsl1(i)
!
!  2.4   prcp       precip rate (kg m-2 s-1)
!
       prcp(i) = rhoh2o * precip(i) / delt
!
!  2.5   sfctmp     air temperature (k) at height zlvl above ground
!        sfctmp =  t1 (temp at lowest model layer, in k)
!  2.6   th2        air potential temperature (k) at height zlvl above ground
!        th2 =  theta1 (see above for theta1 calculations)
!  2.7   q2         mixing ratio at height zlvl above ground (kg kg-1)
!        q2 =  q0 
!
     endif
   enddo
! 
! 3. other forcing (input) data (i):
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
!
!  3.1   sfcspd     wind speed (m s-1) at height zlvl above ground
!        sfcspd = wind (wind speed at lowest model layer, in m/sec)
!  3.2   q2sat      sat mixing ratio at height zlvl above ground (kg kg-1)
!        q2sat =  qs1
!  3.3   dqsdt2     slope of sat specific humidity curve at t=sfctmp (kg kg-1 k-1)
!      -eta:  dqsdt(i,j)=qlms(i,j)*a23m4/(tlm(i,j)-a4)**2
!      -eta vs mrfx: qlms(i,j) = q2sat; tlm(i,j)=t1
!
       dqsdt2(i) = qs1(i) * a23m4/(t1(i)-a4)**2
     endif
   enddo
!
! 4. canopy/soil characteristics (s):
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
!
!  4.1   vegtyp     vegetation type (integer index)
!      
!  4.2   soiltyp    soil type (integer index)
!     
!  4.3   slopetyp   class of sfc slope (integer index)
!    
     if(slptype(i) .gt.  9) slptype(i) = 9
!
!  4.4   shdfac     areal fractional coverage of green vegetation (0.0-1.0)
!        shdfac = sigmaf
!  4.5   shdmin     minimum areal fractional coverage of green vegetation
!                (fraction= 0.0-1.0) <= shdfac
!  4.6   ptu        photo thermal unit (plant phenology for annuals/crops)
!                 (not yet used)
!
     ptu(i) = 0.   ! arbitrary value (for future use)
!
!  4.7   alb        backround snow-free surface albedo (fraction)
!      mrf uses 4-type albedo; albedo determined in phys_lsm_noah is not used
!      modify phys_lsm_noah ---> comment out albedo determination
!
         alb(i) = 0.2 
!  4.8   snoalb     upper bound on maximum albedo over deep snow
!    
!  4.9   tbot       bottom soil temperature (local yearly-mean sfc air temp)
!        tbot = tg3 
!
     endif
   enddo
!
! 5. history (state) variables (h):
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
!
!  5.1  cmc         canopy moisture content (m)
!
       cmc(i) = canopy(i)/rhoh2o  !.. convert from mm to m
!
!  5.2  t1          ground/canopy/snowpack) effective skin temperature (k)
!       t1 = tskin 
!  5.3  stc(nsoil)  soil temp (k)
!       stc(nsoil) = stc(km)
!  5.4  smc(nsoil)  total soil moisture content (volumetric fraction)
!       smc(nsoil) = smc(km)
!  5.5  sh2o(nsoil) unfrozen soil moisture content (volumetric fraction)
!       sh2o(nsoil) = slc(km)
!  5.6  snowh       actual snow depth (m)
!
       snowh(i) = snwdph(i) / rhoh2o     !.. convert from mm to m
!
!  5.7  sneqv       liquid water-equivalent snow depth (m)
!
       sneqv(i) = snoweq(i) / rhoh2o     !.. convert from mm to m
!
!  5.8  albedo      surface albedo including snow effect (unitless fraction)
!     here, albedo is set to 0.2
!     once passed into phys_lsm_noah, it is used to estimate soldn
!
       albedo(i) = 0.2
!
!  5.9  ch          surface exchange coefficient for heat and moisture
!                (m s-1); note: ch is technically a conductance
!     ch and cm are computed in parent model -> input variables
!
       chx(i) = ch(i) * wind(i)            ! compute conductance
!
!  5.10 cm          surface exchange coefficient for momentum (m s-1); note:
!                cm is technically a conductance
!
       cmx(i) = cm(i) * wind(i)            ! compute conductance
     endif
   enddo
!
!    call noah lsm
!
! two new arguments are added in phys_lsm_noah routine
! couple: flag indicating whether to run noah in couple or uncouple mode
! srflag: snow-rain flag (1=snow, 0=rain)
!
   couple = 1            !!<---  run noah lsm in 'couple' mode
   do i = 1,im
     if(slimsk(i) .eq. 1.) then
!
!  from global to local array
!
       do k= 1,km
         smc1d(k) = smc(i,k)
         stc1d(k) = stc(i,k)
         slc1d(k) = slc(i,k)
         sldpth1d(k) = sldpth(i,k)
       enddo
       call  phys_lsm_noah (                                                   &
            couple,srflag(i),                                                  &
            ice(i),delt,zlvl(i),km,sldpth1d,                                   &
            dlwflx(i),radflx(i),sfcprs(i),prcp(i),t1(i),                       &
            q0(i),wind(i),                                                     &
            theta1(i),qs1(i),dqsdt2(i),                                        &
            vegtype(i),soiltyp(i),slptype(i),sigmaf(i),                        &
            shdmin(i),ptu(i),alb(i),snoalb(i),tg3(i),                          &
            cmc(i),tskin(i),stc1d,smc1d,slc1d,                                 &
            snowh(i),sneqv(i),albedo(i),chx(i),cmx(i),                         &
            evap(i),hflx(i),                                                   &
            ec(i),edir(i),et,ett(i),                                           &
            esnow(i),drip(i),dew(i),                                           &
            dm(i), ep(i), gflux(i),                                            &
            flx1(i),flx2(i),flx3(i),                                           &
            snowmt(i),sncovr(i),                                               &
            runoff1(i),runoff2(i),runoff3(i),                                  &
            rc(i),pc(i),rsmin(i),xlai(i),                                      &
            rcs(i),rct(i),rcq(i),rcsoil(i),                                    &
            soilw(i),soilm(i),                                                 &
#ifndef NOAHYDRO
            smcwlt(i),smcdry(i),smcref(i),smcmax(i),nroot(i))
#else
            sn_new(i),                                                         &
            smcwlt(i),smcdry(i),smcref(i),smcmax(i),nroot(i))
!
! acr - snowfall must be the same in phys_lsm_noah as in lsm_noah_main
!
       snowfl(i) = sn_new(i)
#endif
!
!  from local to global array
!
       do k = 1,km
         smc(i,k) = smc1d(k)
         stc(i,k) = stc1d(k)
         slc(i,k) = slc1d(k)
       enddo
#ifdef DBG
       if(iope) then
         print *, 'land diag fields',                                          &
           lat, i, evap(i),hflx(i), ec(i),edir(i),et,ett(i),                   &
           esnow(i),drip(i),dew(i), dm(i), ep(i), gflux(i),                    &
           flx1(i),flx2(i),flx3(i), snowmt(i),sncovr(i),                       &
           runoff1(i),runoff2(i),runoff3(i)
       endif
#endif
     endif 
   enddo
!
!   noah -> gmp/rmp: prepare variables for return to parent model
!
! 6. output (o):
!
!   return the following output fields to parent model
!   eta        actual latent heat flux (w m-2: positive, if upward from sfc)
!   eta = evap
!   sheat      sensible heat flux (w m-2: positive, if upward from sfc)
!   sheat = hflx
!   beta       ratio of actual/potential evap (dimensionless)
!   beta = dm
!   etp        potential evaporation (w m-2)
!   etp = ep
!   ssoil      soil heat flux (w m-2: negative if downward from surface)
!   ssoil = gflux
!   runoff1    surface runoff (m s-1), not infiltrating the surface
!   runoff1 = runoff (in mm s-1) / 1000.
!   runoff2    subsurface runoff (m s-1), drainage out bottom
!   runoff2 = drain (in mm s-1) / 1000.
!
   do i = 1,im
     if(slimsk(i) .eq. 1.) then
!
!  runoff1, runfoff2 (m s-1) -> runoff, drain (mm s-1)
!
       runoff(i) = runoff1(i) * rhoh2o
       drain(i)  = runoff2(i) * rhoh2o
!
!  cmc, snowh, sneqv (m) -> canopy, snwdph, snoweq (mm)
!
       canopy(i) = cmc(i) * rhoh2o
       snwdph(i) = snowh(i) * rhoh2o
       snoweq(i) = sneqv(i) * rhoh2o
#ifdef NOAHYDRO
!
!  transfer snow sublimation (kg/m^2/s)
!
       snowev(i) = esnow(i)
#endif
     endif
   enddo
!
!  compute qsurf and dm
!
   do i = 1,im
     if(slimsk(i) .eq. 1.) then
       qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
       dm(i) = 1.
     endif
   enddo
!
   return
   end subroutine lsm_noah_main
!
!-------------------------------------------------------------------------------
   subroutine phys_lsm_noah(                                                   &
        couple,srflag,                                                         &
        ice,dt,zlvl,nsoil,sldpth,                                              &
        lwdn,radflx,sfcprs,prcp,sfctmp,q2,sfcspd,                              &
        th2,q2sat,dqsdt2,                                                      &
        vegtyp,soiltyp,slopetyp,shdfac,shdmin,ptu,alb,snoalb,tbot,             &
        cmc,t1,stc,smc,sh2o,snowh,sneqv,albedo,ch,cm,                          &
        eta,sheat,                                                             &
! ----------------------------------------------------------------------
! outputs, diagnostics, parameters below generally not necessary when
! coupled with e.g. a nwp model (such as the noaa/nws/ncep mesoscale eta
! model).  other applications may require different output variables. 
! ----------------------------------------------------------------------
        ec,edir,et,ett,esnow,drip,dew,                                         &
        beta,etp,ssoil,                                                        &
        flx1,flx2,flx3,                                                        &
        snomlt,sncovr,                                                         &
        runoff1,runoff2,runoff3,                                               &
        rc,pc,rsmin,xlai,rcs,rct,rcq,rcsoil,                                   &
        soilw,soilm,                                                           &
#ifdef NOAHYDRO
        sn_new,                                                                &
#endif
        smcwlt,smcdry,smcref,smcmax,nroot)
!-------------------------------------------------------------------------------
!
! cheng-hsuan lu (2002)
! (1) modify calling argument
!     1. add 'couple' (=1: coupled, =0: decoupled)
!     2. add 'srflag' for snow-rain detection
!     3. replace soldn by radflx
! (2) modify source code 
!     1. if couple=1, radflx=fdown 
!        if couple=0, radflx=soldn
!     2. move the code fragment where albedo is determined (noah_albedo)
!        to the section where lsm_exch_coeff is called
!     3. invoke the section of noah_albedo and lsm_exch_coeff only if couple=0
!     4. use srflag (instead of sfctmp) to determine snow/rain
!     5. change czil from 0.2 to 0.1
!     6. pass smc (instead of slc) to ddevap
!
! ----------------------------------------------------------------------
! subroutine phys_lsm_noah - version 2.5 - 18 october 2001
! ----------------------------------------------------------------------
! sub-driver for "noah/osu lsm" family of physics subroutines for a
! soil/veg/snowpack land-surface model to update soil moisture, soil
! ice, soil temperature, skin temperature, snowpack water content,
! snowdepth, and all terms of the surface energy balance and surface
! water balance (excluding input atmospheric forcings of downward
! radiation and precip)
! ----------------------------------------------------------------------
! phys_lsm_noah argument list key:
! ----------------------------------------------------------------------
!  c  configuration information
!  f  forcing data
!  i  other (input) forcing data
!  s  surface characteristics
!  h  history (state) variables
!  o  output variables
!  d  diagnostic output
! ----------------------------------------------------------------------
! 1. configuration information (c):
! ----------------------------------------------------------------------
!   couple     couple-uncouple flag  (=1: coupled, =0: decoupled)  **clu_rev*!
!   srflag     flag for snow-rain detection)        **clu_rev*!
!   ice          sea-ice flag  (=1: sea-ice, =0: land)
!   dt          timestep (sec) (dt should not exceed 3600 secs, recommend
!                1800 secs or less)
!   zlvl       height (m) above ground of atmospheric forcing variables
!   nsoil      number of soil layers (at least 2, and not greater than
!                parameter nsold set below)
!   sldpth     the thickness of each soil layer (m)
! ----------------------------------------------------------------------
! 2. forcing data (f):
! ----------------------------------------------------------------------
!   lwdn       lw downward radiation (w m-2; positive, not net longwave)
!   radflx     radiation flux (soldn or fdown)                    **clu_rev*!
!   soldn      solar downward radiation (w m-2; positive, not net solar)
!   sfcprs     pressure at height zlvl above ground (pascals)
!   prcp       precip rate (kg m-2 s-1) (note, this is a rate)
!   sfctmp     air temperature (k) at height zlvl above ground
!   th2        air potential temperature (k) at height zlvl above ground
!   q2         mixing ratio at height zlvl above ground (kg kg-1)
! ----------------------------------------------------------------------
! 3. other forcing (input) data (i):
! ----------------------------------------------------------------------
!   sfcspd     wind speed (m s-1) at height zlvl above ground
!   q2sat      sat mixing ratio at height zlvl above ground (kg kg-1)
!   dqsdt2     slope of sat specific humidity curve at t=sfctmp
!                (kg kg-1 k-1)
! ----------------------------------------------------------------------
! 4. canopy/soil characteristics (s):
! ----------------------------------------------------------------------
!   vegtyp     vegetation type (integer index)
!   soiltyp    soil type (integer index)
!   slopetyp   class of sfc slope (integer index)
!   shdfac     areal fractional coverage of green vegetation
!                (fraction= 0.0-1.0)
!   shdmin     minimum areal fractional coverage of green vegetation
!                (fraction= 0.0-1.0) <= shdfac
!   ptu        photo thermal unit (plant phenology for annuals/crops)
!                (not yet used, but passed to noah_read_parameter for future use in
!                veg parms)
!   alb        backround snow-free surface albedo (fraction), for julian
!                day of year (usually from temporal interpolation of
!                monthly mean values' calling prog may or may not
!                include diurnal sun angle effect)
!   snoalb     upper bound on maximum albedo over deep snow (e.g. from
!                robinson and kukla, 1985, j. clim. & appl. meteor.)
!   tbot       bottom soil temperature (local yearly-mean sfc air
!                temperature)
! ----------------------------------------------------------------------
! 5. history (state) variables (h):
! ----------------------------------------------------------------------
!  cmc         canopy moisture content (m)
!  t1          ground/canopy/snowpack) effective skin temperature (k)
!  stc(nsoil)  soil temp (k)
!  smc(nsoil)  total soil moisture content (volumetric fraction)
!  sh2o(nsoil) unfrozen soil moisture content (volumetric fraction)
!                note: frozen soil moisture = smc - sh2o
!  snowh       actual snow depth (m)
!  sneqv       liquid water-equivalent snow depth (m)
!                note: snow density = sneqv/snowh
!  albedo      surface albedo including snow effect (unitless fraction)
!                =snow-free albedo (alb) when sneqv=0, or
!                =fct(msnoalb,alb,vegtyp,shdfac,shdmin) when sneqv>0
!  ch          surface exchange coefficient for heat and moisture
!                (m s-1); note: ch is technically a conductance since
!                it has been multiplied by wind speed.
!  cm          surface exchange coefficient for momentum (m s-1); note:
!                cm is technically a conductance since it has been
!                multiplied by wind speed.  cm is not needed in phys_lsm_noah
! ----------------------------------------------------------------------
! 6. output (o):
! ----------------------------------------------------------------------
! output variables necessary for a coupled numerical weather prediction
! model, e.g. noaa/nws/ncep mesoscale eta model.  for this application,
! the remaining output/diagnostic/parameter blocks below are not
! necessary.  other applications may require different output variables.
!   eta        actual latent heat flux (kg m-2 s-1: negative, if up from
!            surface)
!   sheat      sensible heat flux (w m-2: negative, if upward from
!            surface)
! ----------------------------------------------------------------------
!   ec         canopy water evaporation (m s-1)
!   edir       direct soil evaporation (m s-1)
!   et(nsoil)  plant transpiration from a particular root (soil) layer
!                 (m s-1)
!   ett        total plant transpiration (m s-1)
!   esnow      sublimation from (or deposition to if <0) snowpack
!                (m s-1)
!   drip       through-fall of precip and/or dew in excess of canopy
!                water-holding capacity (m)
!   dew        dewfall (or frostfall for t<273.15) (m)
! ----------------------------------------------------------------------
!   beta       ratio of actual/potential evap (dimensionless)
!   etp        potential evaporation (kg m-2 s-1)
!   ssoil      soil heat flux (w m-2: negative if downward from surface)
! ----------------------------------------------------------------------
!   flx1       precip-snow sfc (w m-2)
!   flx2       freezing rain latent heat flux (w m-2)
!   flx3       phase-change heat flux from snowmelt (w m-2)
! ----------------------------------------------------------------------
!   snomlt     snow melt (m) (water equivalent)
!   sncovr     fractional snow cover (unitless fraction, 0-1)
! ----------------------------------------------------------------------
!   runoff1    surface runoff (m s-1), not infiltrating the surface
!   runoff2    subsurface runoff (m s-1), drainage out bottom of last
!                soil layer
!   runoff3    numerical trunctation in excess of porosity (smcmax)
!                for a given soil layer at the end of a time step
! ----------------------------------------------------------------------
!   rc         canopy resistance (s m-1)
!   pc         plant coefficient (unitless fraction, 0-1) where pc*etp
!                = actual transp
!   xlai       leaf area index (dimensionless)
!   rsmin      minimum canopy resistance (s m-1)
!   rcs        incoming solar rc factor (dimensionless)
!   rct        air temperature rc factor (dimensionless)
!   rcq        atmos vapor pressure deficit rc factor (dimensionless)
!   rcsoil     soil moisture rc factor (dimensionless)
! ----------------------------------------------------------------------
! 7. diagnostic output (d):
! ----------------------------------------------------------------------
!   soilw      available soil moisture in root zone (unitless fraction
!            between smcwlt and smcmax)
!   soilm      total soil column moisture content (frozen+unfrozen) (m) 
! ----------------------------------------------------------------------
! 8. parameters (p):
! ----------------------------------------------------------------------
!   smcwlt     wilting point (volumetric)
!   smcdry     dry soil moisture threshold where direct evap frm top
!                layer ends (volumetric)
!   smcref     soil moisture threshold where transpiration begins to
!                stress (volumetric)
!   smcmax     porosity, i.e. saturated value of soil moisture
!                (volumetric)
!   nroot      number of root layers, a function of veg type, determined
!              in subroutine noah_read_parameter.
!-------------------------------------------------------------------------------
   use lsm_noah_module, only : noah_read_parameter,   noah_fresh_snow,         &
                               noah_snow_fraction,    noah_soilt_conductivity, &
                               noah_snow_rough,       noah_albedo,             &
                               noah_penman,           noah_canopy_res,         &
                               noah_bare_soil_solver, noah_snow_cover_solver
!-------------------------------------------------------------------------------
   implicit none
! 
   integer,parameter    ::  nsold = 20
! ----------------------------------------------------------------------
! declarations - logical
! ----------------------------------------------------------------------
   logical              ::  frzgra
   logical              ::  saturated
   logical              ::  snowng
! ----------------------------------------------------------------------
! declarations - integer
! ----------------------------------------------------------------------
   integer              ::  ice
   integer              ::  couple                !...**clu_rev*!
   integer              ::  k
   integer              ::  kz
   integer              ::  nsoil
   integer              ::  nroot
   integer              ::  slopetyp
   integer              ::  soiltyp
   integer              ::  vegtyp
! ----------------------------------------------------------------------
! declarations - real
! ----------------------------------------------------------------------
   real                 ::  albedo
   real                 ::  alb
   real                 ::  bexp
   real                 ::  beta
   real                 ::  cfactr
   real                 ::  ch
   real                 ::  cm
   real                 ::  cmc
   real                 ::  cmcmax
   real                 ::  ccsnow
   real                 ::  csoil
   real                 ::  czil
   real                 ::  dew
   real                 ::  df1
   real                 ::  df1h
   real                 ::  df1a
   real                 ::  dksat
   real                 ::  dt
   real                 ::  dwsat
   real                 ::  dqsdt2
   real                 ::  dsoil
   real                 ::  dtot
   real                 ::  drip
   real                 ::  ec
   real                 ::  edir
   real                 ::  esnow
   real                 ::  et(nsoil)
   real                 ::  ett
   real                 ::  frcsno
   real                 ::  frcsoi
   real                 ::  epsca
   real                 ::  eta
   real                 ::  etp
   real                 ::  fdown
   real                 ::  f1
   real                 ::  flx1
   real                 ::  flx2
   real                 ::  flx3
   real                 ::  fxexp
   real                 ::  frzx
   real                 ::  sheat
   real                 ::  hs
   real                 ::  kdt
   real                 ::  lwdn
   real                 ::  pc
   real                 ::  prcp
   real                 ::  ptu
   real                 ::  prcp1
   real                 ::  psisat
   real                 ::  q2
   real                 ::  q2sat
   real                 ::  quartz
   real                 ::  radflx                      ! ..**clu_rev*!
   real                 ::  rch
   real                 ::  refkdt
   real                 ::  rr
   real                 ::  rtdis(nsold)
   real                 ::  runoff1
   real                 ::  runoff2
   real                 ::  rgl
   real                 ::  runoff3
   real                 ::  rsmax
   real                 ::  rc
   real                 ::  rsmin
   real                 ::  rcq
   real                 ::  rcs
   real                 ::  rcsoil
   real                 ::  rct
   real                 ::  rsnow
   real                 ::  sndens
   real                 ::  sncond 
   real                 ::  ssoil
   real                 ::  sbeta
   real                 ::  sfcprs
   real                 ::  sfcspd
   real                 ::  sfctmp
   real                 ::  srflag                      ! ..**clu_rev*!
   real                 ::  shdfac
   real                 ::  shdmin
   real                 ::  sh2o(nsoil)
   real                 ::  sldpth(nsoil)
   real                 ::  smcdry
   real                 ::  smcmax
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  smc(nsoil)
   real                 ::  sneqv
   real                 ::  sncovr
   real                 ::  snowh
   real                 ::  sn_new
   real                 ::  slope
   real                 ::  snup
   real                 ::  salp
   real                 ::  snoalb
   real                 ::  stc(nsoil)
   real                 ::  snomlt
   real                 ::  soldn
   real                 ::  soilm
   real                 ::  soilw
   real                 ::  soilwm
   real                 ::  soilww
   real                 ::  t1
   real                 ::  t1v
   real                 ::  t24
   real                 ::  t2v
   real                 ::  tbot
   real                 ::  th2
   real                 ::  th2v
   real                 ::  topt
   real                 ::  tsnow
   real                 ::  xlai
   real                 ::  zlvl
   real                 ::  zbot
   real                 ::  z0
   real                 ::  zsoil(nsold)
! ----------------------------------------------------------------------
! declarations - parameters
! ----------------------------------------------------------------------
   real,parameter       ::  tfreez = 273.15
   real,parameter       ::  hvap=2.500e+6
   real,parameter       ::  r = 287.04
   real,parameter       ::  cp = 1004.5
! ----------------------------------------------------------------------
!   initialization
! ----------------------------------------------------------------------
   runoff1 = 0.0
   runoff2 = 0.0
   runoff3 = 0.0
   snomlt = 0.0
!
   if (ice .eq. 1) then
     write(6,*) 'error:call noah over sea ice land mask '
   endif
! ----------------------------------------------------------------------
!  the variable "ice" is a flag denoting sea-ice case 
! ----------------------------------------------------------------------
   if (ice .eq. 1) then
! ----------------------------------------------------------------------
! sea-ice layers are equal thickness and sum to 3 meters
! ----------------------------------------------------------------------
     do kz = 1,nsoil
       zsoil(kz) = -3.*float(kz)/float(nsoil)
     end do
   else
! ----------------------------------------------------------------------
! calculate depth (negative) below ground from top skin sfc to bottom of
!   each soil layer.  note:  sign of zsoil is negative (denoting below
!   ground)
! ----------------------------------------------------------------------
     zsoil(1) = -sldpth(1)
     do kz = 2,nsoil
       zsoil(kz) = -sldpth(kz)+zsoil(kz-1)
     end do
   endif
! ----------------------------------------------------------------------
! next is crucial call to set the land-surface parameters, including
! soil-type and veg-type dependent parameters.
! ----------------------------------------------------------------------
   call noah_read_parameter (                                                  &
                vegtyp,soiltyp,slopetyp,                                       &
               cfactr,cmcmax,rsmax,topt,refkdt,kdt,sbeta,                      &
               shdfac,rsmin,rgl,hs,zbot,frzx,psisat,slope,                     &
               snup,salp,bexp,dksat,dwsat,smcmax,smcwlt,smcref,                &
               smcdry,f1,quartz,fxexp,rtdis,sldpth,zsoil,                      &
               nroot,nsoil,z0,czil,xlai,csoil,ptu)
! ----------------------------------------------------------------------
!  initialize precipitation logicals.
! ----------------------------------------------------------------------
   snowng = .false.
   frzgra = .false.
! ----------------------------------------------------------------------
! if sea-ice case, assign default water-equiv snow on top
! ----------------------------------------------------------------------
   if (ice .eq. 1) then
     sneqv = 0.01
     snowh = 0.05
   endif
! ----------------------------------------------------------------------
! if input snowpack is nonzero, then compute snow density "sndens" and
!   snow thermal conductivity "sncond" (note that ccsnow is a function
!   subroutine)
! ----------------------------------------------------------------------
   if (sneqv .eq. 0.0) then
     sndens = 0.0
     snowh = 0.0
     sncond = 1.0
   else
     sndens = sneqv/snowh
!       if(sndens.gt.1.) sndens=1.
     sncond = ccsnow(sndens) 
   endif
! ----------------------------------------------------------------------
! determine if it's precipitating and what kind of precip it is.
! if it's prcping and the air temp is colder than 0 c, it's snowing!
! if it's prcping and the air temp is warmer than 0 c, but the grnd
! temp is colder than 0 c, freezing rain is presumed to be falling.
! ----------------------------------------------------------------------
   if (prcp .gt. 0.0) then
!clu    if (sfctmp .le. tfreez) then   !** clu_rev *!
     if (srflag .eq. 1.) then   
       snowng = .true.
     else
       if (t1 .le. tfreez) frzgra = .true.
     endif
   endif
! ----------------------------------------------------------------------
! if either prcp flag is set, determine new snowfall (converting prcp
! rate from kg m-2 s-1 to a liquid equiv snow depth in meters) and add
! it to the existing snowpack.
! note that since all precip is added to snowpack, no precip infiltrates
! into the soil so that prcp1 is set to zero.
! ----------------------------------------------------------------------
   if ( (snowng) .or. (frzgra) ) then
     sn_new = prcp * dt * 0.001
     sneqv = sneqv + sn_new
     prcp1 = 0.0
! ----------------------------------------------------------------------
! update snow density based on new snowfall, using old and new snow.
! update snow thermal conductivity
! ----------------------------------------------------------------------
     call noah_fresh_snow (sfctmp,sn_new,snowh,sndens)  
     sncond = ccsnow (sndens) 
   else
! ----------------------------------------------------------------------
! precip is liquid (rain), hence save in the precip variable that
! later can wholely or partially infiltrate the soil (along with 
! any canopy "drip" added to this later)
#ifdef NOAHYDRO
! acr addition - zero out new snow when it is not snowing
! ----------------------------------------------------------------------
     sn_new = 0
#else
! ----------------------------------------------------------------------
#endif
     prcp1 = prcp
   endif
!
! ** clu_rev *!
! only determine snowcover (needed for thermal conductivity calculations)
! the surface albedo modification is done together with lsm_exch_coeff (optional)
!
! ----------------------------------------------------------------------
! determine snowcover over land.
! ----------------------------------------------------------------------
   if (ice .eq. 0) then
! ----------------------------------------------------------------------
! if snow depth=0, set snow fraction=0
! ----------------------------------------------------------------------
     if (sneqv .eq. 0.0) then
       sncovr = 0.0
!clu      albedo = alb
     else
! ----------------------------------------------------------------------
! determine snow fractional coverage.
! ----------------------------------------------------------------------
       call noah_snow_fraction (sneqv,snup,salp,snowh,sncovr)
!clu      call noah_albedo (alb,snoalb,shdfac,shdmin,sncovr,tsnow,albedo)
     endif
   else
! ----------------------------------------------------------------------
! set snow cover over sea-ice
! ----------------------------------------------------------------------
     sncovr = 1.0
!clu    albedo = 0.60
   endif
! ----------------------------------------------------------------------
! thermal conductivity for sea-ice case
! ----------------------------------------------------------------------
   if (ice .eq. 1) then
     df1 = 2.2
   else
! ----------------------------------------------------------------------
! next calculate the subsurface heat flux, which first requires
! calculation of the thermal diffusivity.  treatment of the
! latter follows that on pages 148-149 from "heat transfer in 
! cold climates", by v. j. lunardini (published in 1981 
! by van nostrand reinhold co.) i.e. treatment of two contiguous 
! "plane parallel" mediums (namely here the first soil layer 
! and the snowpack layer, if any). this diffusivity treatment 
! behaves well for both zero and nonzero snowpack, including the 
! limit of very thin snowpack.  this treatment also eliminates
! the need to impose an arbitrary upper bound on subsurface 
! heat flux when the snowpack becomes extremely thin.
! ----------------------------------------------------------------------
! first calculate thermal diffusivity of top soil layer, using
! both the frozen and liquid soil moisture, following the 
! soil thermal diffusivity function of peters-lidard et al.
! (1998,jas, vol 55, 1209-1224), which requires the specifying
! the quartz content of the given soil class 
! (see routine noah_read_parameter)
! ----------------------------------------------------------------------
     call noah_soilt_conductivity (df1,smc(1),quartz,smcmax,sh2o(1))
! ----------------------------------------------------------------------
! next add subsurface heat flux reduction effect from the 
! overlying green canopy, adapted from section 2.1.2 of 
! peters-lidard et al. (1997, jgr, vol 102(d4))
! ----------------------------------------------------------------------
     df1 = df1 * exp(sbeta*shdfac)
   endif
! ----------------------------------------------------------------------
! finally "plane parallel" snowpack effect following 
! v.j. linardini reference cited above. note that dtot is
! combined depth of snowdepth and thickness of first soil layer
! ----------------------------------------------------------------------
   dsoil = -(0.5 * zsoil(1))
   if (sneqv .eq. 0.) then
     ssoil = df1 * (t1 - stc(1) ) / dsoil
   else
     dtot = snowh + dsoil
     frcsno = snowh/dtot
     frcsoi = dsoil/dtot
!
! 1. harmonic mean (series flow)
!        df1 = (sncond*df1)/(frcsoi*sncond+frcsno*df1)
     df1h = (sncond*df1)/(frcsoi*sncond+frcsno*df1)
! 2. arithmetic mean (parallel flow)
!        df1 = frcsno*sncond + frcsoi*df1
     df1a = frcsno*sncond + frcsoi*df1
!
! 3. geometric mean (intermediate between harmonic and arithmetic mean)
!        df1 = (sncond**frcsno)*(df1**frcsoi)
! test - mbek, 10 jan 2002
! weigh df by snow fraction
!        df1 = df1h*sncovr + df1a*(1.0-sncovr)
!        df1 = df1h*sncovr + df1*(1.0-sncovr)
     df1 = df1a*sncovr + df1*(1.0-sncovr)
! ----------------------------------------------------------------------
! calculate subsurface heat flux, ssoil, from final thermal diffusivity
! of surface mediums, df1 above, and skin temperature and top 
! mid-layer soil temperature
! ----------------------------------------------------------------------
     ssoil = df1 * (t1 - stc(1) ) / dtot
   endif
! ----------------------------------------------------------------------
! determine surface roughness over snowpack using snow condition from
! the previous timestep.
! ----------------------------------------------------------------------
   if (sncovr .gt. 0.) then
     call noah_snow_rough (sncovr,z0)
   endif
! ----------------------------------------------------------------------
! next call routine lsm_exch_coeff to calculate the sfc exchange coef (ch) for
! heat and moisture.
!
! note !!!
! comment out call lsm_exch_coeff, 
! if lsm_exch_coeff already called in calling program
! (such as in coupled atmospheric model).
!
! note !!!
! do not call lsm_exch_coeff until after above call to noah_read_parameter, 
! in case alternative values of roughness length (z0) and zilintinkevich coef
! (czil) are set there via namelist i/o.
!
! note !!!
! routine lsm_exch_coeff returns a ch that represents the wind spd times the
! "original" nondimensional "ch" typical in literature.  hence the ch
! returned from lsm_exch_coeff has units of m/s.  the important companion
! coefficient of ch, carried here as "rch", is the ch from lsm_exch_coeff times
! air density and parameter "cp".  "rch" is computed in "call noah_penman".
! rch rather than ch is the coeff usually invoked later in eqns.
!
! note !!!
! lsm_exch_coeff also returns the surface exchange coefficient for momentum, cm,
! also known as the surface drage coefficient, but cm is not used here.
! ----------------------------------------------------------------------
! calc virtual temps and virtual potential temps needed by subroutines
! lsm_exch_coeff and noah_penman.
! ----------------------------------------------------------------------
   t2v = sfctmp * (1.0 + 0.61 * q2 )
! ----------------------------------------------------------------------
! ** clu_rev *!
! decouple mode
! determine surface albedo 
! compute surface exchange coefficients
! set soldn = radflx
! compute fdown from soldn, albedo, and lwdn
!
! couple mode
! set fdown = radflx
! compute soldn from fdown, albedo, and lwdn
!
   if (couple .eq. 0) then            !......decouple mode
! ----------------------------------------------------------------------
! determine albedo over land.
! ----------------------------------------------------------------------
     if (ice .eq. 0) then
! ----------------------------------------------------------------------
! if snow depth=0, set albedo=snow free albedo.
! ----------------------------------------------------------------------
       if (sneqv .eq. 0.0) then
         albedo = alb
       else
! ----------------------------------------------------------------------
! determine surface albedo modification due to snowdepth state.
! ----------------------------------------------------------------------
            call noah_albedo (alb,snoalb,shdfac,shdmin,sncovr,tsnow,albedo)
       endif
     else
! ----------------------------------------------------------------------
! set albedo over sea-ice
! ----------------------------------------------------------------------
       albedo = 0.60
     endif
! ----------------------------------------------------------------------
! comment out below 2 lines if call lsm_exch_coeff is commented out, 
! i.e. in the coupled model.
! ----------------------------------------------------------------------
!
     t1v = t1 * (1.0 + 0.61 * q2)
     th2v = th2 * (1.0 + 0.61 * q2)
     call lsm_exch_coeff (zlvl,z0,t1v,th2v,sfcspd,czil,cm,ch)
! ----------------------------------------------------------------------
! calculate total downward radiation (solar plus longwave) needed in
! noah_penman ep subroutine that follows
! ----------------------------------------------------------------------
     soldn = radflx                      !..downward solar radiation  
     fdown = soldn*(1.0-albedo) + lwdn   !..total downward radiation
   else                                   !.....couple mode               
     fdown = radflx                      !..total downward radiation
     soldn = (fdown-lwdn)/(1.0-albedo)   !..dn sw rad 
   endif
!----------------------------------------------------------------------
! call noah_penman subroutine to calculate potential evaporation (etp), and
! other partial products and sums save in common/rite for later
! calculations.
! ----------------------------------------------------------------------
    call noah_penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,t24,ssoil,           &
                      q2,q2sat,etp,rch,epsca,rr,snowng,frzgra,                 &
                      dqsdt2,flx2)
! ----------------------------------------------------------------------
! call noah_canopy_res to calculate the canopy resistance and convert it into pc
! if nonzero greenness fraction
! ----------------------------------------------------------------------
   if (shdfac .gt. 0.) then
! ----------------------------------------------------------------------
!  frozen ground extension: total soil water "smc" was replaced 
!  by unfrozen soil water "sh2o" in call to noah_canopy_res below
! ----------------------------------------------------------------------
     call noah_canopy_res (soldn,ch,sfctmp,q2,sfcprs,sh2o,zsoil,nsoil,         &
                            smcwlt,smcref,rsmin,rc,pc,nroot,q2sat,dqsdt2,      &
                            topt,rsmax,rgl,hs,xlai,                            &
                            rcs,rct,rcq,rcsoil)
   endif
! ----------------------------------------------------------------------
! now decide major pathway branch to take depending on whether snowpack
! exists or not:
! ----------------------------------------------------------------------
   esnow = 0.0
   if (sneqv .eq. 0.0) then
     call noah_bare_soil_solver (etp,eta,prcp,smc,smcmax,smcwlt,               &
               smcref,smcdry,cmc,cmcmax,nsoil,dt,shdfac,                       &
               sbeta,q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,                      &
               stc,epsca,bexp,pc,rch,rr,cfactr,                                &
               sh2o,slope,kdt,frzx,psisat,zsoil,                               &
               dksat,dwsat,tbot,zbot,runoff1,runoff2,                          &
               runoff3,edir,ec,et,ett,nroot,ice,rtdis,                         &
               quartz,fxexp,csoil,                                             &
               beta,drip,dew,flx1,flx2,flx3)
   else
     call noah_snow_cover_solver (etp,eta,prcp,prcp1,                          &
                  snowng,smc,smcmax,smcwlt,                                    &
                  smcref,smcdry,cmc,cmcmax,nsoil,dt,                           &
                  sbeta,df1,                                                   &
                  q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,stc,epsca,               &
                  sfcprs,bexp,pc,rch,rr,cfactr,sncovr,sneqv,sndens,            &
                  snowh,sh2o,slope,kdt,frzx,psisat,snup,                       &
                  zsoil,dwsat,dksat,tbot,zbot,shdfac,runoff1,                  &
                  runoff2,runoff3,edir,ec,et,ett,nroot,snomlt,                 &
                  ice,rtdis,quartz,fxexp,csoil,                                &
                  beta,drip,dew,flx1,flx2,flx3)
     esnow = eta
   endif
! ----------------------------------------------------------------------
!   prepare sensible heat (h) for return to parent model
! ----------------------------------------------------------------------
   sheat = -(ch * cp * sfcprs)/(r * t2v) * ( th2 - t1 )
! ----------------------------------------------------------------------
!  convert units and/or sign of total evap (eta), potential evap (etp),
!  subsurface heat flux (s), and runoffs for what parent model expects
!  convert eta from kg m-2 s-1 to w m-2
! ----------------------------------------------------------------------
   eta = eta*hvap
   etp = etp*hvap
! ----------------------------------------------------------------------
! convert the sign of soil heat flux so that:
!   ssoil>0: warm the surface  (night time)
!   ssoil<0: cool the surface  (day time)
! ----------------------------------------------------------------------
   ssoil = -1.0*ssoil      
! ----------------------------------------------------------------------
!  convert runoff3 (internal layer runoff from supersat) from m to m s-1
!  and add to subsurface runoff/drainage/baseflow
! ----------------------------------------------------------------------
   runoff3 = runoff3/dt
   runoff2 = runoff2+runoff3
! ----------------------------------------------------------------------
! total column soil moisture in meters (soilm) and root-zone 
! soil moisture availability (fraction) relative to porosity/saturation
! ----------------------------------------------------------------------
   soilm = -1.0*smc(1)*zsoil(1)
   do k = 2,nsoil
     soilm = soilm+smc(k)*(zsoil(k-1)-zsoil(k))
   end do
   soilwm = -1.0*(smcmax-smcwlt)*zsoil(1)
   soilww = -1.0*(smc(1)-smcwlt)*zsoil(1)
   do k = 2,nroot
     soilwm = soilwm+(smcmax-smcwlt)*(zsoil(k-1)-zsoil(k))
     soilww = soilww+(smc(k)-smcwlt)*(zsoil(k-1)-zsoil(k))
   end do
   soilw = soilww/soilwm
! ----------------------------------------------------------------------
! end subroutine phys_lsm_noah
! ----------------------------------------------------------------------
   return
   end subroutine phys_lsm_noah
!
