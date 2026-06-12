#include <define.h>
   subroutine funct_lsm
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!      --- [osu_funct_df]
!      --- [osu_funct_kt]
!      --- [osu_funct_ktsoil]
!      --- [osu_funct_thsat]
!      --- [osu_funct_twlt]
!      --- [func_snoweb]
!               |--- [stabilityc]
!      --- [func_soiltemp]
!      --- [ccsnow]
!      --- [ffrh2o]
!      --- [ddevap]
!      --- [ssnksrc]
!
!-------------------------------------------------------------------------------
   end subroutine funct_lsm
!
!-------------------------------------------------------------------------------
   function osu_funct_df(theta,ktype)
!-------------------------------------------------------------------------------
   use paramodel, only : ntype=>nstype_
   use comfcst, only   : dfk, tsat
!-------------------------------------------------------------------------------
!
   w = (theta / tsat(ktype)) * 20. + 1.
   kw = w
   kw = min(kw,21)
   kw = max(kw,1)
   osu_funct_df = dfk(kw,ktype) + (w - kw) * (dfk(kw+1,ktype) - dfk(kw,ktype))
!
   return
   end function osu_funct_df
!
!-------------------------------------------------------------------------------
   function osu_funct_kt(theta,ktype)
!-------------------------------------------------------------------------------
   use paramodel, only : ntype=>nstype_
   use comfcst, only   : ktk, tsat
!-------------------------------------------------------------------------------
!
   w = (theta / tsat(ktype)) * 20. + 1.
   kw = w
   kw = min(kw,21)
   kw = max(kw,1)
   osu_funct_kt = ktk(kw,ktype)                                                & 
          + (w - kw) * (ktk(kw+1,ktype) - ktk(kw,ktype))
!
   return
   end function osu_funct_kt
!
!-------------------------------------------------------------------------------
   function osu_funct_ktsoil(theta,ktype)
!-------------------------------------------------------------------------------
   use comfcst, only   : dfkt, tsat
!-------------------------------------------------------------------------------
   integer  ::  ktype, kw
   real     ::  osu_funct_ktsoil, theta, w
!
   w = (theta / tsat(ktype)) * 20. + 1.
   kw = int(w)
   kw = min(kw,21)
   kw = max(kw,1)
   osu_funct_ktsoil = dfkt(kw,ktype) + (w-kw)*(dfkt(kw+1,ktype)-dfkt(kw,ktype))
!
   return
   end function osu_funct_ktsoil
!
!-------------------------------------------------------------------------------
   real function osu_funct_thsat(ktype)
!-------------------------------------------------------------------------------
   use comfcst, only      : tsat
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer              ::  ktype
   osu_funct_thsat = tsat(ktype)
!
   return
   end function osu_funct_thsat
!
!-------------------------------------------------------------------------------
   real function osu_funct_twlt(ktype)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ktype
!-------------------------------------------------------------------------------
   osu_funct_twlt = .1
!
   return
   end function osu_funct_twlt
!
!-------------------------------------------------------------------------------
   function func_snoweb(tsf)
!-------------------------------------------------------------------------------
!
! source file:       func_snoweb.f
! purpose:           calculate snow pack energy balance
!
! program history log:     
!     03-06-01  ji chen
!     modified from 4.0.3 uw vic (snowpackenergybalance)
!
! notes:
!    in c language, the values of function parameters are passed one-way, 
!    i.e., the value is from the upper to down; but two-way in fortran
!
! reference: bras, r. a., hydrology, an introduction to hydrologic
!            science, addisson wesley, inc., reading, etc., 1990.
!
!-------------------------------------------------------------------------------
   use vic_veglib
!-------------------------------------------------------------------------------
#include <vartyp.h>
!
! input variables
!
   real  ::  tsf            ! snow surface temperature (c)
!
! output variables 
!
   real  ::  func_snoweb    ! residual snow energy (w/m2)
#include <vic_snoweb.h>
!
! local variables 
!
   real  ::  le              ! latent heat of vaporization (j/kg)
   real  ::  ls              ! latent heat of sublimation (j/kg)
   real  ::  d               ! vegetation displace (m)
   real  ::  z0              ! vegetation roughness (m)
   real  ::  rac             ! corrected resistance (s/m)
   real  ::  longradout      ! emitted longwave radiation (w/m2)
   real  ::  rad             ! net radiation (w/m2)
   real  ::  vpsnow          ! saturated vapor pressure in snow pack (pa)
   real  ::  vapormf         ! mass flux of vapor to or from snow (m/s)
   real  ::  advectede       ! energy advected by precipitation (w/m2)
   real  ::  deltcc          ! change in cold content (w/m2)
   real  ::  spkh            ! snow pack heat flux (w/m2)
   real  ::  restterm        ! rest term in surface energy balance (w/m2)
   real  ::  stabilityc, svp ! function names
!
!  cp            ::  specific heat of moistair at constant (j/kg/k)
!  lf            ::  latent heat of freezing (j/kg) at 0c
!  rhoh2o        ::  water density (kg/m^3)
!
   real, parameter  ::  cp=1013.0,lf=3.337e5,rhoh2o=1.e3
!
!  sigma         ::  stefan-boltzmann constant (w/m^2/k^-4)
!  t0c           ::  ice/water mix temperature (k)
!
   real, parameter  ::  huge_resist=1.e2, sigma=5.67e-8, t0c=273.15
!
!  ch_water      ::  volumetric heat capacity (j/(m3*c)) of water
!  ch_ice        ::  volumetric heat capacity (j/(m3*c)) of ice
!
   real, parameter  ::  ch_water=4186.8e3,ch_ice=2100.0e3
!------------------------------------------------------------------------------
!
   le = (2.501 - 0.002361 * tsf) * 1.0e6
   ls = (677. - 0.07 * tsf) * 4.1868 * 1000.0
!
!      d  = veg_d(emonth, evtype)
!      z0 = veg_rough(emonth, evtype)
! 
! correct vic_exch conductance for stable conditions
! 
   d = 0.0
   z0 = 0.03
!
   if (ewind .gt. 0.0) then
      rac = era/stabilityc(2.0, d, tsf+t0c, etair+t0c, ewind, z0)
   else
      rac = huge_resist
   end if
! 
! calculate longwave exchange and net radiation
! 
   longradout = sigma * (tsf+t0c)**4
!
   rad = enetsht + elongw - longradout
! 
! calculate the sensible heat flux
! 
   osenhs = erhoair * cp * (etair - tsf)/rac
!
#ifdef DBGVIC
   print *,'in func_snoweb osenhs',erhoair,cp,etair,tsf,era,rac
#endif
! 
! calculate the mass flux of ice to or from the surface layer
! calculate saturated vapor pressure in snow pack, (eq. 3.32, bras 90)
! 
   vpsnow  = svp(tsf+t0c)
   vapormf = erhoair * (0.622/epgcm) * (evpair - vpsnow)/rac
   vapormf = vapormf/rhoh2o
!
#ifdef DBGVIC
   print*,'in func_snoweb vapormf ',vapormf,erhoair,epgcm,evpair,              &
          vpsnow,rac,enetsht,elongw,longradout
#endif
!
   if (evpd .eq. 0.0 .and. vapormf .lt. 0.0) vapormf=0.0
  
! 
! calculate latent heat flux
! 
   if (tsf .ge. 0.0) then
! 
! melt conditions: use latent heat of vaporization 
! 
      olaths = le * vapormf * rhoh2o
   else 
! 
! accumulation: use latent heat of sublimation (eq. 3.19, bras 1990)
! 
      olaths = ls * vapormf * rhoh2o
   end if
! 
! calculate advected heat flux from rain
! 
   advectede = (ch_water * (etair-tsf) * efrtf) / edtime
! 
! calculate change in cold content
! 
   deltcc = ch_ice * esfswq * (eoldtsf-tsf)/ edtime
! 
! discussion: (ji 2/2004) the contribution from deltcc to the snow 
!    surface energy balance should be not very high and we should limit
!    the valud of deltcc. the reason is that in the energy balance equation
!    all the energy terms are facing to a virtual snow surface instead
!    of the surface snow layer.
! 
!      deltcc = deltcc/10.0
!      if(abs(deltcc).gt.100) deltcc = sign(100.0, deltcc)
!      deltcc = 0.0
!
#ifdef DBGVIC
   print *,'in func_snoweb ---',ch_ice,esfswq,eoldtsf,tsf,edtime
#endif
! 
! calculate snow pack heat flux from the ground surface
! 
   if(esnowd.gt.0.) then
      spkh = (2.9302e-6)* efrsn*efrsn * (etgrnd-tsf)/esnowd
   else 
      spkh = 0.0
   end if
!
!      if(abs(spkh).gt.100) spkh = sign(100.0, spkh)
! 
! calculate net energy exchange at the snow surface
! 
#ifdef DBGVIC
   print *,'in func_snoweb +++',tsf,rad,osenhs,olaths,advectede,               &
           spkh, deltcc, edtime, esfwater,lf,etgrnd,esnowd
#endif
!
   restterm  = rad + osenhs + olaths + advectede + spkh + deltcc
   orefreeze = (esfwater*lf*rhoh2o)/edtime
!
#ifdef DBGVIC
   print *,'in func_snoweb ===',restterm, orefreeze
#endif
!
! please check wigmosta, vail and lettenmaier (p 1670, wrr 1994)
!   for the computation of snow melt energy
!   bascially, when tsf < 0c no snow melt, 
!                 tsf = 0, allowing the heat flux avaiable for snowmelt
! 
   if (abs(tsf).lt.1.e-6 .and. restterm .ge. (-orefreeze)) then
      orefreeze = -restterm           ! available energy input over
                                      ! cold content used to melt
      restterm  = 0.0
   else
      restterm  = restterm + orefreeze  ! add this value to the pack
   end if
!
   func_snoweb  = restterm
!
#ifdef DBGVIC
   print *,'in func_snoweb ***',restterm, orefreeze
#endif
!
   return
   end function func_snoweb
!
!-------------------------------------------------------------------------------
   function stabilityc(z, d, tsurf, tair, wind, z0)
!-------------------------------------------------------------------------------
!                 
! subprogram: stabilityc
!
! abstract: this function calculates the stability correction for 
!   exchange of sensible heat between surface&atmosphere
!
! program history log:
!   2003-06/07  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (stabilitycorrection.c)
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
!
!-------------------- input variables ------------------------------------------
!
   real z          ! reference height (m)
   real d          ! displacement height (m)
   real tsurf      ! surface temperature (k)
   real tair       ! air temperature (k)
   real wind       ! wind speed (m/s)
   real z0         ! roughness length (m)
!
!-------------------- output variables -----------------------------------------
!
   real stabilityc ! multiplier for vic_exch resistance
!
!----------------------- local variables ---------------------------------------
!
   real rilimit    ! upper limit for richardsons number
   real ri         ! richardsons number 
   real ricr       ! critical richardsons number
   real g          ! gravity acceleration (m/s2)

   parameter (ricr = 0.2, g=9.81)
!-------------------------------------------------------------------------------
!
! calculate effect of atmospheric stability using richardson number approach
!
   stabilityc = 1.0
   if (tsurf.ne.tair.and.z.gt.d) then
! 
! non-neutral conditions
! 
     ri = g*(tair-tsurf)*(z-d)/((tair+tsurf)/2.0*wind*wind)
!    
     if(z0.le.0.0) z0=0.01
!
     if(((z-d)/z0).le.1) then
       rilimit = tair/((tair+tsurf)/2.0 * 5.0)
     else
       rilimit = tair/((tair+tsurf)/2.0*(log((z-d)/z0) + 5.0))
     end if
!
     if (ri .gt. rilimit) ri = rilimit
!    
     if (ri .gt. 0.0) then 
       stabilityc = (1 - ri/ricr) * (1 - ri/ricr)
     else      
       if (ri .lt. -0.5) ri = -0.5
       stabilityc = sqrt(1.0 - 16.0 * ri)
     end if
   end if
!
#ifdef DBGVIC
   print *,'vic in stability check stability ', stabilityc
#endif
!
   return
   end function stabilityc
!
!-------------------------------------------------------------------------------
   function func_soiltemp(t)
!-------------------------------------------------------------------------------
!
! source file:       func_soiltemp.f
! purpose:           soil temperature equation
!
! program history log:
!     03-06-01  ji chen
!     modified from 4.0.3 uw vic (soil_thermal_eqn.c)
!
! comments:
!      please see equation (8) in cherkauer and lettenmaier (jgr 1999)
!      for the soil thermal equation
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
!
! input variables 
!
   real t                  ! current node soil temperature (k)
!
! output variables 
!
   real func_soiltemp      ! residure of thermal balance function
!
! common blocks
!
#include <vic_soileb.h>
! 
! local variables 
!
   real sfc_max_unfwat     ! function name
   real t0c                ! ice/water mix temperature (k)
   parameter (t0c=273.15)
!-------------------------------------------------------------------------------
   if(t.lt.t0c) then
      ice = soilm - sfc_max_unfwat(t, soilmx, bubble, expt)
   else
      ice = 0.0
   end if
!
   if(ice.lt.0.) ice = 0
!
   func_soiltemp = (t-t0c)*ee - aa*(tl-tu) -                                   &
            bb*(tl+tu-2*t0c-gamma*fprime)-cc-dd*(ice-ice0)
!
#ifdef DBGVIC   
!   print *,'in func_soiltemp ',func_soiltemp,t,tl,tu,
!     &         gamma,fprime,ice,ice0,soilm,soilmx
!        print *,'   func_soiltemp ',(t-t0c)*ee,aa*(tl-tu),
!     &          bb*(tl+tu-2*t0c-gamma*fprime),cc,dd*(ice-ice0)
#endif
   return
   end function func_soiltemp
!
!-------------------------------------------------------------------------------
   function ccsnow (dsnow)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! calculate snow termal conductivity
!
   real                 ::  c
   real                 ::  dsnow
   real                 ::  ccsnow
   real,parameter       ::  unit = 0.11631
! ----------------------------------------------------------------------
! ccsnow in units of cal/(cm*hr*c), returned in w/(m*c)
! basic version is dyachkova equation (1960), for range 0.1-0.4
! ----------------------------------------------------------------------
   c=0.328*10**(2.25*dsnow)
   ccsnow=unit*c
!
! de vaux equation (1933), in range 0.1-0.6
!      ccsnow=0.0293*(1.+100.*dsnow**2)
! e. andersen from flerchinger
!      ccsnow=0.021+2.51*dsnow**2        
! end function ccsnow
!
   return                                                      
   end function ccsnow
!
!-------------------------------------------------------------------------------
   function ffrh2o (tkelv,smc,sh2o,smcmax,bexp,psis)
!-------------------------------------------------------------------------------
!
! calculate amount of supercooled liquid soil water content if
! temperature is below 273.15k (t0).  requires newton-type iteration to
! solve the nonlinear implicit equation given in eqn 17 of koren et al
! (1999, jgr, vol 104(d16), 19569-19585).
!
! new version (june 2001): much faster and more accurate newton
! iteration achieved by first taking log of eqn cited above -- less than
! 4 (typically 1 or 2) iterations achieves convergence.  also, explicit
! 1-step solution option for special case of parameter ck=0, which
! reduces the original implicit equation to a simpler explicit form,
! known as the "flerchinger eqn". improved handling of solution in the
! limit of freezing point temperature t0.
!
! input:
!
!   tkelv.........temperature (kelvin)
!   smc...........total soil moisture content (volumetric)
!   sh2o..........liquid soil moisture content (volumetric)
!   smcmax........saturation soil moisture content (from noah_read_parameter)
!   b.............soil type "b" parameter (from noah_read_parameter)
!   psis..........saturated soil matric potential (from noah_read_parameter)
!
! output:
!   ffrh2o.........supercooled liquid water content
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real, parameter   ::  ck = 8.0
   real, parameter   ::  blim = 5.5
   real, parameter   ::  error = 0.005
   real, parameter   ::  hlice = 3.335e5
   real, parameter   ::  gs = 9.81
   real, parameter   ::  dice = 920.0
   real, parameter   ::  dh2o = 1000.0
   real, parameter   ::  t0 = 273.15
!
   real              ::  bexp
   real              ::  bx
   real              ::  denom
   real              ::  df
   real              ::  dswl
   real              ::  fk
   real              ::  ffrh2o
   real              ::  psis
   real              ::  sh2o
   real              ::  smc
   real              ::  smcmax
   real              ::  swl
   real              ::  swlk
   real              ::  tkelv
   integer           ::  nlog
   integer           ::  kcount
!-------------------------------------------------------------------------------
!
! limits on parameter b: b < 5.5  (use parameter blim)
! simulations showed if b > 5.5 unfrozen water content is
! non-realistically high at very low temperatures.
!
   bx = bexp
   if (bexp .gt. blim) bx = blim
! 
! initializing iterations counter and iterative solution flag.
!
   nlog=0
   kcount=0
! 
!  if temperature not significantly below freezing (t0), sh2o = smc
!
   if (tkelv .gt. (t0 - 1.e-3)) then                                           
     ffrh2o = smc                                                              
   else
     if (ck .ne. 0.0) then
!-------------------------------------------------------------------------------
! option 1: iterated solution for nonzero ck
! in koren et al, jgr, 1999, eqn 17
!-------------------------------------------------------------------------------
! initial guess for swl (frozen content)
!
       swl = smc-sh2o
!
! keep within bounds.
!
       if (swl .gt. (smc-0.02)) swl = smc-0.02
       if (swl .lt. 0.) swl = 0.
! 
!  start of iterations
!
       do while ( (nlog .lt. 10) .and. (kcount .eq. 0) )
         nlog = nlog+1
         df = alog(( psis*gs/hlice ) * ( ( 1.+ck*swl )**2. ) *                 &
         ( smcmax/(smc-swl) )**bx) - alog(-(tkelv-t0)/tkelv)
         denom = 2. * ck / ( 1.+ck*swl ) + bx / ( smc - swl )
         swlk = swl - df/denom
! 
! bounds useful for mathematical solution.
!
         if (swlk .gt. (smc-0.02)) swlk = smc - 0.02
         if (swlk .lt. 0.) swlk = 0.
! 
! mathematical solution bounds applied.
! 
         dswl = abs(swlk-swl)
         swl = swlk
! 
! if more than 10 iterations, use explicit method (ck=0 approx.)
! when dswl less or eq. error, no more iterations required.
!
         if ( dswl .le. error )  then
           kcount = kcount+1
         endif
       end do
! 
!  end of iterations
! 
! bounds applied within do-block are valid for physical solution.
! 
       ffrh2o = smc - swl
! 
! end option 1
! 
     endif
!-------------------------------------------------------------------------------
! option 2: explicit solution for flerchinger eq. i.e. ck=0
! in koren et al., jgr, 1999, eqn 17
! apply physical bounds to flerchinger solution
!-------------------------------------------------------------------------------
     if (kcount .eq. 0) then
       fk = (((hlice/(gs*(-psis)))*                                            &
            ((tkelv-t0)/tkelv))**(-1/bx))*smcmax
       if (fk .lt. 0.02) fk = 0.02
       ffrh2o = min (fk, smc)
! 
! end option 2
! 
     endif
   endif
! 
   return
   end function ffrh2o
!
!-------------------------------------------------------------------------------
   function ddevap (etp1,smc,zsoil,shdfac,smcmax,bexp,                         &
                   dksat,dwsat,smcdry,smcref,smcwlt,fxexp)
!-------------------------------------------------------------------------------
!
! calculate direct soil evaporation
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                 ::  bexp
   real                 ::  ddevap
   real                 ::  dksat
   real                 ::  dwsat
   real                 ::  etp1
   real                 ::  fx
   real                 ::  fxexp
   real                 ::  shdfac
   real                 ::  smc
   real                 ::  smcdry
   real                 ::  smcmax
   real                 ::  zsoil
   real                 ::  smcref
   real                 ::  smcwlt
! 
! direct evap a function of relative soil moisture availability, linear
! when fxexp=1.
! 
   fx = ((smc - smcdry) / (smcmax - smcdry))**fxexp
! 
! fx > 1 represents demand control
! fx < 1 represents flux control
! 
   fx = max ( min ( fx, 1. ) ,0. )
! 
! allow for the direct-evap-reducing effect of shade
! 
   ddevap = fx * ( 1.0 - shdfac ) * etp1
! 
! end function ddevap
! 
   return
   end function ddevap
!
!-------------------------------------------------------------------------------
   function ssnksrc (tavg,smc,sh2o,zsoil,nsoil,                                &
                    smcmax,psisat,bexp,dt,k,qtot) 
!-------------------------------------------------------------------------------
!
! calculate sink/source term of the termal diffusion equation. (sh2o) is
! available liqued water.
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!      
   integer              ::  k
   integer              ::  nsoil
!      
   real                 ::  bexp
   real                 ::  df
   real                 ::  dt
   real                 ::  dz
   real                 ::  dzh
   real                 ::  free
   real                 ::  ffrh2o
   real                 ::  psisat
   real                 ::  qtot
   real                 ::  sh2o
   real                 ::  smc
   real                 ::  smcmax
   real                 ::  ssnksrc
   real                 ::  tavg
   real                 ::  tdn
   real                 ::  tm
   real                 ::  tup
   real                 ::  tz
   real                 ::  x0
   real                 ::  xdn
   real                 ::  xh2o
   real                 ::  xup
   real                 ::  zsoil (nsoil)
!
   real,parameter       ::  dh2o = 1.0000e3
   real,parameter       ::  hlice = 3.3350e5
   real,parameter       ::  t0 = 2.7315e2
!-------------------------------------------------------------------------------
   if (k .eq. 1) then
     dz = -zsoil(1)
   else
     dz = zsoil(k-1)-zsoil(k)
   endif
! 
! via function ffrh2o, compute potential or 'equilibrium' unfrozen
! supercooled free water for given soil type and soil layer temperature.
! function frh20 invokes eqn (17) from v. koren et al (1999, jgr, vol.
! 104, pg 19573).  (aside:  latter eqn in journal in centigrade units.
! routine ffrh2o use form of eqn in kelvin units.)
! 
   free = ffrh2o(tavg,smc,sh2o,smcmax,bexp,psisat)
! 
! in next block of code, invoke eqn 18 of v. koren et al (1999, jgr,
! vol. 104, pg 19573.)  that is, first estimate the new amountof liquid
! water, 'xh2o', implied by the sum of (1) the liquid water at the begin
! of current time step, and (2) the freeze of thaw change in liquid water
! implied by the heat flux 'qtot' passed in from routine noah_soilt_land.
! second, determine if xh2o needs to be bounded by 'free' (equil amt) or
! if 'free' needs to be bounded by xh2o.
! 
   xh2o = sh2o + qtot*dt/(dh2o*hlice*dz)
! 
! first, if freezing and remaining liquid less than lower bound, then
! reduce extent of freezing, thereby letting some or all of heat flux
! qtot cool the soil temp later in routine noah_soilt_land.
! 
   if ( xh2o .lt. sh2o .and. xh2o .lt. free) then 
     if ( free .gt. sh2o ) then
       xh2o = sh2o
     else
       xh2o = free
     endif
   endif
!
! second, if thawing and the increase in liquid water greater than upper
! bound, then reduce extent of thaw, thereby letting some or all of heat
! flux qtot warm the soil temp later in routine noah_soilt_land.
! 
   if ( xh2o .gt. sh2o .and. xh2o .gt. free )  then
     if ( free .lt. sh2o ) then
       xh2o = sh2o
     else
       xh2o = free
     endif
   endif 
!
   if (xh2o .lt. 0.) xh2o = 0.
   if (xh2o .gt. smc) xh2o = smc
! 
! calculate phase-change heat source/sink term for use in 
! routine noah_soilt_land and update liquid water to reflcet 
! final freeze/thaw increment.
!
   ssnksrc = -dh2o*hlice*dz*(xh2o-sh2o)/dt
   sh2o = xh2o
!
77 return
   end function ssnksrc
!
!-------------------------------------------------------------------------------
