#include <define.h>
   subroutine lsm_vic_funct
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [lsm_vic_funct]
!      |
!      |--- [vic_funct_energy] *
!      |         |--- [soil_temprfl] *
!      |         |--- [arno_evap] *
!      |--- [vic_funct_penman] *
!      |         |--- [svp_slope] *
!      |--- [canopy_evap] *
!      |         |--- [vic_transpiration] *
!      |--- [root_brent] *
!      |--- [svp] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine lsm_vic_funct
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
#ifdef VICLSM1
   function vic_funct_energy(ts)
!-------------------------------------------------------------------------------
!
! subprogram: vic_funct_energy
!
! abstract: this subroutine computes the surface energy balance for bare 
!   soil and vegetation uncovered by snow.  it computes outgoing longwave,
!   sensible heat flux, ground heat flux, and storage of heat in the thin
!   upper layer, based on the given surface temperature.
!        the energy balance equation used comes from xu liangs paper 
!  "insights of the ground heat flux in land surface 
!   parameterization schemes."
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (func_surf_energy_bal.c)
!
!-------------------------------------------------------------------------------
   use vic_veglib
#include <vartyp.h>
#include <vic_surfeb.h>
! ------------------- input variables ----------------------------------
!
   real  ::  ts              ! surface temperature (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ----------------------------------
!
   real  ::  vic_funct_energy     ! surface energy balance error (w/m2)
! ----------------------------------------------------------------------
!
! ------------------- common blocks ------------------------------------
!
   integer  ::  kdt,jdt,inistp,limlow,maxstp,numsum,nummax
   common/comfveri/ kdt,jdt,inistp,limlow,maxstp,numsum,nummax
! ----------------------------------------------------------------------

!
! ---------------------- local variables -------------------------------
!
   real     ::  le, ls          ! latent heat of vaporization/sublimation (j/kg)
   real     ::  snow_flux       ! snow heat flux (w/m2)
   real     ::  kappa_snow      ! snow thermal conductivity (w/m/k)
   real     ::  ice             ! soil layer ice content (mm)
   real     ::  ice0            ! last soil layer ice content (mm)
   real     ::  delth           ! land layer heat storage change
   real     ::  evap            ! evapotranspiration (m/s)
   real     ::  ts_old          ! last surface temperature (k)
   real     ::  t1_old          ! last first layer soil temperature (k)
   real     ::  t1              ! soil temperature (k)
   real     ::  tmp_t           ! temporary temperature (k)
   real     ::  smmax           ! maximum soil moisture (mm)
   real     ::  rad             ! net radiation (w/m2)

   real     ::  canopy_evap     ! funcation name (m/s)
   real     ::  arno_evap       ! funcation name (m/s)
   real     ::  sfc_max_unfwat  ! funcation name

   integer  ::  m               ! loop index
   !
   ! rhoice         ::  density of ice (kg/m^3)
   ! rhoh2o         ::  water density (kg/m^3)
   !
   real, parameter  ::  rhoice = 917,rhoh2o=1.e3
   !
   ! cp       ::  specific heat of moist air at constant (j/kg/k)
   ! lf       ::  latent heat of freezing (j/kg) at 0c
   ! sigma    ::  stefan-boltzmann constant (w/m^2/k^-4)
   ! t0c      ::  ice/water mix temperature (k)
   !
   real, parameter  ::  cp=1013.0,lf=3.337e5,sigma=5.6730e-8,t0c=273.15
!
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! initialize latent heat and sensible heat fluxes
! ----------------------------------------------------------------------
   olath = 0.0
   osenh = 0.0
   mfcwt = efcwt
   mfrtf = efrtf

! ----------------------------------------------------------------------
! compute surface temperature at half time step
! ----------------------------------------------------------------------
!
   if(esnowd .gt. 0.0) then
!
! ----------------------------------------------------------------------
! compute energy flux through snow pack
! ----------------------------------------------------------------------
     kappa_snow = 2.9302e-6 * efrsn * efrsn
     snow_flux  = kappa_snow * (eftsf-ts) / esnowd
   else
     snow_flux  = 0.0
   endif
!
! ----------------------------------------------------------------------
! use finite difference method to explicitly solve ground heat
! flux at soil thermal nodes (cherkauer and lettenmaier, 1999)
! ----------------------------------------------------------------------
   ts_old = eftnd(1)
   t1_old = eftnd(2)
   ice0   = efsic(1)
!
   mftnd(1) = ts             ! the first node of new t_profile is ts
!
#ifdef DBGVIC
   if(kdt.eq.341)then
   print *,'in vic_funct_energy check eftnd ', eftnd
   endif
#endif
!
   call soil_temprfl(nsl,  edtime,   ealpnd,  ebetnd,                          &
                  egamnd,   eftnd,   ekapnd,   ecsnd,                          &
                   esmnd,  eicend,   esmxnd,  ebubnd,                          &
                  eexpnd,   mftnd)
!
   t1 = mftnd(2)
!
#ifdef DBGVIC
   if(kdt.eq.341)then
     print *,'in vic_funct_energy check mftnd ', mftnd
   endif
#endif
! ----------------------------------------------------------------------
! initialize mfsmc
! ----------------------------------------------------------------------
   do m = 1,msl
     mfsmc(m) = efsmc(m)
   end do
!
! ----------------------------------------------------------------------
! compute the ground heat flux to the top soil layer
! ----------------------------------------------------------------------
   ogrndh = ekapnd(1)*(t1 - ts)/edph(1)
!      if(abs(ogrndh).gt.100) ogrndh = sign(100.0,ogrndh)
!
! ----------------------------------------------------------------------
! compute the current ice content of the top soil layer
! ----------------------------------------------------------------------
   tmp_t = (ts+ t1)/2.0
   smmax = esmxnd(1)*edph(1)*1000.0
!
   if(tmp_t .lt. 273.15)then
     ice = efsmc(1) - sfc_max_unfwat(tmp_t,smmax,ebubnd(1),eexpnd(1))
     if(ice .lt. 0.) ice=0.
   else
     ice = 0.
   end if
#ifdef DBGVIC
!
   if(kdt.eq.341)then
     print *,'in vic_funct_energy ice=',ice, efsmc(1),tmp_t,smmax,             &
                                         ebubnd(1),eexpnd(1)
   endif
#endif
!
   delth = ecsnd(1)*((ts_old+t1_old)/2.-tmp_t)*edph(1)/edtime
#ifdef DBGVIC
!
   if(kdt.eq.341)then
     print *,'delth1=',delth,ecsnd(1),(ts_old+t1_old)/2, tmp_t
   endif
#endif
!
   delth = delth - rhoice*lf*((ice0-ice)/1000.0)/edtime
!
! ----------------------------------------------------------------------
! discussion: (ji 2/2004) the contribution from delth to the surface
!    energy balance should be not very high and we should limit the
!    valud of delth. the reason is that in the energy balance equation
!    all the energy terms are facing to a virtual land surface instead
!    of the first soil layer.
! ----------------------------------------------------------------------
   delth = delth/10.0  ! 6/2004
   if(abs(delth).gt.100) delth = sign(100.0, delth)
!
   odelth = delth
!
! temporary set delth=0.0 !ji 2/2004
!      delth = 0.0
!
#ifdef DBGVIC
   if(kdt.eq.341)then
      print *,'delth2=',delth,' ice=',ice,' ice0=',ice0,efsmc(1)
   endif
#endif
!
! ----------------------------------------------------------------------
! compute evapotranspiration if not snow covered
! ----------------------------------------------------------------------
   if(esncvfr.lt.1.0) then               ! no snow cover
     rad = enetsht + elongw - sigma*ts**4 + ogrndh + delth

#ifdef DBGVIC
     if(kdt.eq.341)then
       print *,'rad=',rad,enetsht,elongw,sigma*ts**4,ogrndh,delth
     endif
#endif
!    print *,' in vic_funct_energy efflai',efflai
!
     if(efflai.gt.0.) then ! snow free canopy
!
! ----------------------------------------------------------------------
! compute net surface radiation for evaporation estimates
! ----------------------------------------------------------------------
#ifdef DBGVIC
       if(kdt.eq.341)then
         print *,'canopy_evap', msl,  evtype,   edtime,  efflai,               &
                          efcwt,   efrtf,    efsmc,   efsic,                   &
                           ewcr,    ewpw,     efrt,     rad,                   &
                           evpd,                                               &
                          epgcm,   efrtf
       endif
#endif
       evap = canopy_evap(msl,  evtype,   edtime,  efflai,                     &
                          efcwt,   efrtf,    efsmc,   efsic,                   &
                           ewcr,    ewpw,     efrt,     rad,                   &
                        erhoair,  evpair,     evpd,     era,                   &
                        enetsht,   etair,    epgcm,   mfcwt,                   &
                          mfrtf,   mfsmc)
!      print *,'in vic_funct_energy canop evap',evap
     else
! ----------------------------------------------------------------------
! compute net surface radiation for evaporation estimates
! ----------------------------------------------------------------------
#ifdef DBGVIC
       if(kdt.eq.341)then
         print *,'arno_evap ',efsmc(1), efsic(1),   smmax, esmr(1),            &
                             ebinf,  edph(1),  edtime,     rad,                &
                            evpair, evpd,                                      &
                         epgcm
       endif
#endif                     
       mfrtf = efrtf
       evap = arno_evap(        ts,                                            &
                          efsmc(1), efsic(1),   smmax,  esmr(1),               &
                             ebinf,  edph(1),  edtime,      rad,               &
                           erhoair,   evpair,    evpd,      era,               &
                       enetsht,    etair,   epgcm, mfsmc(1))
!      print *,'in vic_funct_energy arno evap',evap
     end if
!
! ----------------------------------------------------------------------
! compute the latent heat flux from the surface and covering vegetation
!   le: latent heat of vaporization (j/kg) 
!   ls: latent heat of sublimation (j/kg)
! ----------------------------------------------------------------------
!
     if(ts.lt.t0c) then
       ls = (677.-0.07*(ts-t0c))*4.1868*1000.0
       olath = -rhoh2o*ls*evap
     else
       le = (2.501-0.002361*(ts-t0c))*1.e6
       olath = -rhoh2o*le*evap
     endif
!
#ifdef DBGVIC
     if(kdt.eq.341)then
       print *,'in vic_funct_energy check',olath,evap
     endif
#endif         
! ----------------------------------------------------------------------
! compute the sensible heat flux from the surface
! ----------------------------------------------------------------------
     osenh = erhoair*cp*(etair - ts)/era
#ifdef DBGVIC
     if(kdt.eq.341)then
       print *,'in vic_funct_energy osenh ',osenh,erhoair,etair,ts,era,cp
     endif
#endif
   end if
! ----------------------------------------------------------------------
! compute surface energy balance error
! ----------------------------------------------------------------------
   vic_funct_energy = (1.0-esncvfr)*(enetsht + elongw - sigma*ts**4 +          &
                 olath + osenh) + ogrndh + delth +                             &
                 esncvfr*snow_flux
#ifdef DBGVIC
!
   if(kdt.gt.339)then
     print *,'in vic_funct_energy esnowd(m) =',esnowd,esncvfr,enetsht,         & 
          elongw,sigma*ts**4,ts,olath,osenh,ogrndh,delth,snow_flux
     print *,'in vic_funct_energy vic_funct_energy=',vic_funct_energy
     endif
#endif
!
   return
   end function vic_funct_energy
#endif
!
!-------------------------------------------------------------------------------
   subroutine soil_temprfl(nsl,  dtime,  alpnd,  betnd,                        &
                         gamnd,     t0,  kapnd,   csnd,                        &
                          smnd,  icend,  smxnd,  bubnd,                        &
                         expnd,     tn)
#ifdef VICLSM1
#include <vartyp.h>
!-------------------------------------------------------------------------------
!
! abstract: this subroutine iteratively solve the soil temperature
!           profile using a numerical difference equation. the solution
!           equation is second order in space, and first order in time.
!
! program history log:
!   2003-06,07,08,09  ji chen   modified from 4.0.3 uw vic (frozen_soil.c)
!
! comments:
!    please check eqns (3), (4), (5), (6), (7) and (8) in cherkauer &
!    lettenmaier (jgr 1999) for the computing soiltemp profile
!-------------------------------------------------------------------------------
!
! ------------------- input variables ----------------------------------
!
   integer  ::  nsl           ! number of soil thermal nodes
!
   real     ::  dtime            ! time step (s)
   real     ::  alpnd(nsl)       ! coeff for computing soiltemp profile
   real     ::  betnd(nsl)       ! coeff for computing soiltemp profile
   real     ::  gamnd(nsl)       ! coeff for computing soiltemp profile
   real     ::  t0(nsl)          ! olde soil temperature profile (k)
   real     ::  kapnd(nsl)       ! soil layer thermal conductivity (w/m/k)
   real     ::  csnd(nsl)        ! soil layer heat capacity (j/m^3/k)
   real     ::  smnd(nsl)        ! soil smndure at node (mm/mm)
   real     ::  icend(nsl)       ! soil icend content at node (mm/mm)
   real     ::  smxnd(nsl)       ! maximum soil moisture at node (mm/mm)
   real     ::  bubnd(nsl)       ! bubbling pressure of soil (cm)
   real     ::  expnd(nsl)       ! parameter for ksat with soil moisture (n/a)   
!
! ------------------- output variables ---------------------------------
!
   real     ::  tn(nsl)          ! new soil temperature profile (k)
!
! ------------------- common block -------------------------------------
!
#include <vic_soileb.h>
!
! ---------------------- local variables -------------------------------
!
   integer  ::  j, itcount
!
   real     ::  maxdiff, diff, oldt, alpnd2(nsl)
   real     ::  a(nsl), b(nsl), c(nsl), d(nsl), e(nsl)
!
   real     ::  root_brent       ! function name
!
   character(len=80)   ::  ctrfct  ! determine function
!
!  soildt     ::  used to bracket soil temperatures while
!  threshold  ::  temperature profile iteration threshold
!
   integer, parameter  ::  maxit = 10
   real, parameter     ::  soildt=0.25, threshold = 1.e-2
!
!  rhoice     ::  density of ice (kg/m^3)
!  lf         ::  latent heat of freezing (j/kg) at 0c
!  t0c        ::  ice/water mix temperature (k)
!
   real, parameter  ::  rhoice = 917.0,lf=3.337e5,t0c=273.15
! ----------------------------------------------------------------------
!
   do j = 1,nsl
     alpnd2(j) = alpnd(j)*alpnd(j)
   end do
   
   do j = 2,nsl
     if(j.eq.nsl) then
       a(j) = betnd(j)*dtime*(kapnd(j)-kapnd(j-1))
     else
       a(j) = betnd(j)*dtime*(kapnd(j+1)-kapnd(j-1))
     endif
     b(j) = 2.*alpnd2(j)*dtime*kapnd(j)
     c(j) = alpnd2(j)*betnd(j)*csnd(j)*(t0(j)-t0c)
     d(j) = alpnd2(j)*betnd(j)*rhoice*lf
     e(j) = alpnd2(j)*betnd(j)*csnd(j)                                         &
                + 4.*kapnd(j)*alpnd2(j)*dtime
   end do
!
! tn(1) = ts, initialized in vic_funct_energy
!
   do j = 2,nsl
     tn(j)=t0(j)
   end do
!
   itcount = 0
!
 5000 if(itcount.lt.maxit) then
     itcount = itcount + 1
     maxdiff = threshold
     do j = 2,nsl
       oldt = tn(j)
! ----------------------------------------------------------------------
! 2nd order variable kapnd equation 
! ----------------------------------------------------------------------
       if(j.eq.nsl) then
         tl = tn(j)
       else
         tl = tn(j+1)
       endif
       tu = tn(j-1)
!
       fprime = (tl-tu)/alpnd(j)
!         
       if(tn(j).gt.t0c) then
         tn(j) = (a(j)*(tl-tu)                                                 &
                + b(j)*(tl+tu-2*t0c-gamnd(j)*fprime)                           &
                + c(j) + d(j)*(0.-icend(j))) / e(j) + t0c
       else
         soilm  = smnd(j)
         soilmx = smxnd(j)
         bubble = bubnd(j)
         expt   = expnd(j)
         ice0   = icend(j)
         gamma  = gamnd(j)
         aa = a(j)
         bb = b(j)
         cc = c(j)
         dd = d(j)
         ee = e(j)
         ctrfct = 'func_soiltemp'
#ifdef DBGVIC
!        write(6,*) 'in soil_temprfl bef root_brent',j,aa,bb,cc,dd,ee
#endif
         tn(j) =root_brent(t0(j)+soildt, t0(j)-soildt, ctrfct)
         if(abs(t0(j)-tn(j)).gt.5) then
           write(6,*) 'in soil_temprfl aft ',j,t0(j),tn(j)
!          stop
         endif
       end if
!
       diff= abs(oldt-tn(j))
       if(diff .gt. maxdiff) maxdiff=diff
     end do
   end if ! 5000 if
!
   if(maxdiff .gt. threshold) then
     goto 5000
   end if
!
#endif /* VICLSM1 end */
   return
   end subroutine soil_temprfl
!
!-------------------------------------------------------------------------------
   function arno_evap(   ts,                                                   &
                      fsmc1,  soili,  smmax,  smres,                           &
                       binf,     d1,  dtime,    rad,                           &
                     rhoair,  vpair,    vpd,     ra,                           &
                     netsht,   tair,   pgcm, mfsmc1)
#include <vartyp.h>
!-------------------------------------------------------------------------------
!
! source file:       arno_evap.f
! org:               ecpc/crd/sio/ucsd 
! date:              june, july, august & september 2003
! prgm history:      modified from 4.0.3 uw vic (arno_evap.c)
!
! abstract: this routine compute evaporation based on the assumption that 
!    evaporation is at the potential for the area which is saturated, and at
!    some percentage of the potential for the area which is partial saturated.
!
!-------------------------------------------------------------------------------
!
! ------------------- input variables ----------------------------------
!
   real     ::  ts             ! surface temperature (k)
   real     ::  soili          ! soil ice in layer 1 (mm)
   real     ::  smmax          ! maximum soil moisture (mm)
   real     ::  smres          ! residual soil moisture (mm)
   real     ::  binf           ! vic infiltration parameter (n/a)
   real     ::  d1             ! soil layer 1 thickness (m)
   real     ::  dtime          ! time step (second)
   real     ::  rad            ! net radiation (w/m2)
   real     ::  rhoair         ! air density (kg/m3)
   real     ::  vpair          ! actual vapor pressure of air (pa)
   real     ::  vpd            ! vapor pressure deficit (pa)
   real     ::  ra             ! aerodynamical resistance (s/m)
   real     ::  netsht         ! net shortwave radiation (w/m2)
   real     ::  tair           ! air temperature (k)
   real     ::  pgcm           ! pressure (pa)
!
! ------------------- modified variables -------------------------------
!
   real     ::  fsmc1          ! soil moisture in layer 1 (mm)
   real     ::  mfsmc1         ! modified soil moisture in layer 1 (mm)
!
! ------------------- output variables ---------------------------------
!
   real     ::  arno_evap      ! evaporation over bare soil (m/s)
!
! ---------------------- local variables -------------------------------
!
   real     ::  evap          ! evaporation
   real     ::  moist         ! available soil moisture (mm)
   real     ::  epot          ! potential bare soil evaporation (m/s)
   real     ::  max_infil     ! maximum infitration capacity (mm)
   real     ::  tmp           ! infiltration capacity with moist (mm)
   real     ::  tmp_evap      ! temporary varible for evaporation (mm/time step)
   real     ::  ratio         ! for computing tmp
   real     ::  as            ! saturated area (fraction)
   real     ::  beta_asp
   real     ::  dummy         ! temporary variable
   real     ::  tmpsum

   integer  ::  num_term, i ! loop index

   real     ::  vic_funct_penman, svp ! function name
   real, parameter  ::  rhoh2o=1.e3   ! water density (kg/m^3)
!
   evap = 0
!
   moist = fsmc1 - soili
   if(moist.gt.smmax) moist= smmax
!
! ----------------------------------------------------------------------
! calculate the potential bare soil evaporation (m/s)
! ----------------------------------------------------------------------
! in the original c-vic (uw), the negative rad is also used in noah_penman
! method, which will introduce errors in computing energy balance.
! therefore, we control the negatvie rad computation. if rad is negative,
! we use aerodynamic method to compute evap and assume canopy 
! temperature equals to air temperature.
! ----------------------------------------------------------------------
!
   if(rad.gt.0.0)then
     epot = vic_funct_penman(rad,     vpd,      ra,      0.0,                  &
                              0.0,     1.0,     1.0,     tair,                 &
                           netsht,    pgcm,       0)
   else
     epot = rhoair*(0.622/pgcm)*(svp(ts)-vpair)/ra
     epot = epot / rhoh2o
     if (vpd.eq.0.0.and.epot.gt.0.0) epot=0.0
   endif
!
#ifdef DBGVIC
   print *,'in arno_evap epot=',epot
#endif
! ----------------------------------------------------------------------
! compute temporary infiltration rate based on given soil_moist.
! ----------------------------------------------------------------------
   max_infil = (1.0+binf)*smmax
   if(binf .eq. (-1.0) .or. moist.eq.smmax) then
     tmp = max_infil
   else 
     ratio = 1.0 - moist /smmax
! ----------------------------------------------------------------------
! if(ratio < small && ratio > -small) ratio = 0.
! ----------------------------------------------------------------------
     ratio = ratio**(1.0 / (binf + 1.0))
     tmp = max_infil*(1.0 - ratio)       !infiltration capacity
   endif
!
! ----------------------------------------------------------------------
! evaporation see eqs.(14) & (15) in liang et al. (jgr 1994) derivation.
! ----------------------------------------------------------------------
   if(tmp .ge. max_infil.or.epot.le.0.0) then
     evap = epot
   else
! ----------------------------------------------------------------------
! compute as. as is % area saturated, 1-as is % area that is unsaturated
! ----------------------------------------------------------------------
     ratio = tmp/max_infil 
     ratio = 1.0 - ratio
     ratio = ratio**binf
     as = 1 - ratio                 ! saturated area
! ----------------------------------------------------------------------
! compute the beta function in the arno evaporation model using
! the first 30 terms in the power expansion expression.
! ----------------------------------------------------------------------
     ratio = ratio**(1.0/binf)
     dummy = 1.0
     do num_term = 1, 30
       tmpsum = ratio
       do i = 1, num_term
         tmpsum = tmpsum*ratio
       end do
       dummy = dummy + binf*tmpsum/(binf+num_term)
     end do
         
     beta_asp = as + (1.0-as)*(1.0-ratio)*dummy
     evap = epot*beta_asp
   endif
!
#ifdef DBGVIC
   print *,'in arno_evap evap=',evap
#endif
!
! ----------------------------------------------------------------------
! evaporation cannot exceed available soil moisture.
! evaporation second soil layer = 0.0
! ----------------------------------------------------------------------
!
   tmp_evap = evap*1000.0*dtime   ! m/s -> mm/time step
!
   if((moist-smres).gt.0.0) then
     if(tmp_evap .gt. (moist - smres)) then
       evap = (moist -  smres)/1000.0/dtime
       tmp_evap=moist-smres
     end if
   else
     tmp_evap = 0.0
     evap = 0.0
   endif
#ifdef DBGVIC
!
   print *,'in arno_evap tmp_evap',tmp_evap,moist,smres,evap
#endif
!
   mfsmc1 = fsmc1 - tmp_evap
!  mfsmc1 = fsmc1
   arno_evap = evap       ! m/s
!
   return
   end function arno_evap
!
!-------------------------------------------------------------------------------
   function vic_funct_penman(rad,     vpd,      ra,      rs,                   &
                            rarc,     lai, gsm_inv,    tair,                   &
                        netshort,    pgcm,     rgl)
#include <vartyp.h>
!-------------------------------------------------------------------------------
!
! subprogram: vic_funct_penman
!
! abstract: calculate et using the combination eq
!
! comments:
!    please check equations (11), (12), (13), (14), (15), (16), and (3)
!    in wigmosta, vail and lettenmaier (wrr 1994) for computing evap.
!
! program history log:
!   2003-06/09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic c version (penman.c)
!
!-------------------------------------------------------------------------------
!
! ------------------- input variables ----------------------------------
!
   real rad        ! net radiation (w/m2)
   real vpd        ! vapor pressure deficit (pa)
   real ra         ! vic_exch resistance (s/m)
   real rs         ! minimum stomatal resistance (s/m)
   real rarc       ! architectural resistance (s/m)
   real lai        ! lai
   real gsm_inv    ! soil moisture stress factor
   real netshort   ! net shortwave radiation (w/m2)
   real tair       ! air temperature (k)
   real pgcm       ! pressure (pa)
   real rgl        ! value of solar radiation below which there 
                   ! will be no transpiration (ranges from ~30 w/m^2
                   ! for trees to ~100 w/m^2 for crops)
!
! ------------------- output variables ---------------------------------
!
   real vic_funct_penman  ! evapotranspiration (m/s)
!
! ---------------------- local variables -------------------------------
!
   real ta         ! air temperature (c)
   real f
   real tmp_value
   real slope      ! slope of saturated vapor pressure curve (pa/c)
   real rc         ! canopy resistance
   real rhoair     ! density of air in (kg/m3)
   real le         ! latent heat of vaporization (j/kg)
   real gamma      ! psychrometric constant (pa/c)
   real tfactor    ! factor for canopy resist based on temperature
   real vpdfactor  ! factor for canopy resistance based on vpd
   real dayfactor  ! factor for canopy resist based on photosynthesis
   real svp_slope  ! function name

   real cp         ! specific heat of moist air at constant (j/kg/k)
   real t0c        ! ice/water mix temperature (k)
   parameter (cp=1013.0, t0c=273.15)

   real closure    ! pa
   real rsmax      ! maximum stomatal resistance (s/m)
   real vpdminftr  ! minimum vpd factor
   parameter (closure = 4000.0,rsmax=5000,vpdminftr=0.1)
!
! ----------------------------------------------------------------------
! calculate the slope of the saturated vapor pressure curve in pa/k
! ----------------------------------------------------------------------
!
   slope = svp_slope(tair)
!
   ta = tair - t0c
!
! ----------------------------------------------------------------------
! calculate resistance factors (wigmosta et al., 1994)
! ----------------------------------------------------------------------
!
   if(rs .gt. 0.) then
     f = netshort / rgl
     dayfactor = (1. + f)/(f + rs/rsmax)
   else
     dayfactor = 1.
   end if
!
   tfactor = .08 * ta - 0.0016 * ta * ta
   if(tfactor.le.0.0) tfactor = 1.e-10
!
   vpdfactor = 1 - vpd/closure
   if(vpdfactor .lt. vpdminftr) vpdfactor = vpdminftr
!
! ----------------------------------------------------------------------
! calculate canopy resistance in s/m
! ----------------------------------------------------------------------
!
   tmp_value = lai * gsm_inv * tfactor * vpdfactor
!
   if(tmp_value.ne.0.0) then
     rc = rs/(lai * gsm_inv * tfactor * vpdfactor) * dayfactor
     if(rc .gt. rsmax) rc = rsmax
   else
     rc = rsmax
   end if
!
!      print *,'penman rc',rc,rs,lai,gsm_inv
!      print *,'penman tfactor ',tfactor,ta,vpdfactor,vpd
!      print *,'penman dayfactor',dayfactor,netshort,rgl
!
! ----------------------------------------------------------------------
! calculate latent heat of vaporization. eq. 4.2.1 in handbook of 
! hydrology, assume ts is ta.  le (j/kg)
! ----------------------------------------------------------------------
   le = (2.501 - 0.002361 * ta) * 1.0e6 
!  
! ----------------------------------------------------------------------
! calculate gamma (pa/c) eq. 4.2.28. handbook of hydrology 
! ----------------------------------------------------------------------
   gamma = 1628.6 * pgcm/le
!  
! ----------------------------------------------------------------------
! calculate air density (kg/m3), using eq. 4.2.4 handbook of hydrology 
! ----------------------------------------------------------------------
   rhoair = 0.003486 * pgcm/(275.0 + ta)
! 
! ----------------------------------------------------------------------
! calculate the evaporation in mm/s (not dividing by the density 
! of water (~1000 kg/m3)), the result ends up being in mm instead of m 
! ----------------------------------------------------------------------
!   
   vic_funct_penman = (slope * rad + rhoair * cp * vpd/ra)/                    &
         (le * (slope + gamma * (1 + 0.1*(rc + rarc)/ra)))

!      print *,'vpenman=',vic_funct_penman
!      print *,'vic_funct_penman ',slope,rad,rhoair,cp,vpd,ra,rc,rarc,gamma,le
!      print *,'penman ',slope * rad,rhoair * cp * vpd/ra,
!     &                  gamma * (1 + (rc + rarc)/ra)

   vic_funct_penman = vic_funct_penman / 1000.0         ! mm/s -> m/s
   if (vpd .gt. 0.0 .and. vic_funct_penman .lt. 0.0) vic_funct_penman = 0.0 ! no dew
!
   return
   end function vic_funct_penman
!
!-------------------------------------------------------------------------------
   function svp_slope(temp)
!-------------------------------------------------------------------------------
#include <vartyp.h>
!-------------------------------------------------------------------------------
! 
! subprogram: svp_slope
!
! abstract: the gradient of d(svp)/dt using handbook
!   of hydrology eqn 4.2.3
!
!
! program history log:
!   2003-06/07  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (svp.c)
!
!-------------------------------------------------------------------------------
!
! --------------- input variables from atmospheric model ---------------
!
   real  ::  temp        ! air temperature (k)
!
! ------------------- output variables ---------------------------------
!
   real  ::  svp_slope   ! the gradient of saturated vapor pressure (pa/c)
!
! --------------- local variables --------------------------------------
!
   real  ::  t           ! air temperature (c)
   real  ::  svp         ! function name
   real, parameter  ::  b_svp = 17.269, c_svp = 237.3, t0c=273.15
!-------------------------------------------------------------------------------
   t = temp - t0c
!   
   svp_slope = (b_svp * c_svp) / ((c_svp + t) * (c_svp + t))                   &
               * svp(temp)

   return
   end function svp_slope
!
!-------------------------------------------------------------------------------
   function canopy_evap(msl,  vtype,   dtime,  fflai,                          &
                       fcwt,  rfall,    fsmc,   fsic,                          &
                        wcr,    wpw,     frt,    rad,                          &
                     rhoair,  vpair,     vpd,     ra,                          &
                     netsht,   tair,    pgcm,  mfcwt,                          &
                      mfrtf,  mfsmc)
!-------------------------------------------------------------------------------
!
! source file:       canopy_evap.f
! org:               ecpc/crd/sio/ucsd 
! date:              june, july, august & september 2003
! prgm history:      modified from 4.0.3 uw vic (canopy_evap.c)
!
! abstract: this function computes the evaporation, traspiration and 
!           throughfall of the vegetation types for multi-layered model.
!
! comments:
!    please check equations (1), (5), (7), (8), (9), (10), and (11) 
!    in liang, lettenmaier, wood, and burges (wrr 1994) for computing et
!
!-------------------------------------------------------------------------------
   use vic_veglib
#include <vartyp.h>
! ------------------- input variables ----------------------------------
!
   integer msl       ! number of soil layer
   integer vtype     ! vegetation type
   real    ::  dtime        ! time step (second)
   real    ::  fflai        ! leaf area index
   real    ::  fsic(msl)    ! soil ice (mm)
   real    ::  wcr(msl)     ! ~70% field capacity (mm)
   real    ::  wpw(msl)     ! wilting point (mm)
   real    ::  frt(msl)     ! root content (fraction)
   real    ::  rad          ! net radiation (w/m2)
   real    ::  rhoair       ! air density (kg/m^3)
   real    ::  vpair        ! actual vapor pressure of air (pa)
   real    ::  vpd          ! vapor pressure deficit (pa)
   real    ::  ra           ! vic_exch resistance (s/m)
   real    ::  netsht       ! net shortwave radiation (w/m2)
   real    ::  tair         ! air temperature (k)
   real    ::  pgcm         ! pressure (pa) 
!
! ------------------- modified variables -------------------------------
!
   real    ::  rfall        ! precipitation during time step (m/time step)
   real    ::  fcwt         ! canopy water (m)
   real    ::  fsmc(msl)    ! soil moisture (liquid + ice) (mm)
   real    ::  mfrtf        ! rain throughfall (m/time step)
   real    ::  mfcwt        ! canopy water (m)
   real    ::  mfsmc(msl)   ! soil moisture (liquid + ice) (mm)
!
! ------------------- output variables ---------------------------------
!
   real    ::  canopy_evap  ! et (m/s)
!
! ------------------- local variables ----------------------------------
!
   integer  ::  k 

   real     ::  wdmax        ! maximum canopy water hold capacity (m)
   real     ::  evap         ! evaporation (m/s) (from vic_funct_penman)
   real     ::  tmp_wdew     ! temporary trap canopy water (m)
   real     ::  frac         ! wet leaf fraction (fraction)
   real     ::  canopyevap   ! canopy evaporation from wet leaf (m/time step)

   real     ::  rs           ! minimum stomatal resistance (s/m)
   real     ::  rarc         ! architectural resistance (s/m)
   real     ::  rgl          ! value of solar radiation below which there 
                             ! will be no transpiration (ranges from ~30 w/m^2
                             ! for trees to ~100 w/m^2 for crops)
   real     ::  layevp(msl)  ! transpiration from each soillayer (m/time step)
   real     ::  vic_funct_penman  ! function name
   real     ::  svp               ! function name

   real     ::  tmp_evap          ! temporary total evaporation (m/time step)
   real     ::  bas, f
!
!  rhoh2o        ::  water density (kg/m^3)
!  lai_wf        ::  leaf water factor for interception storage (m)
!
   real, parameter  ::  rhoh2o=1.e3
   real, parameter  ::  lai_wf=0.0002
! ----------------------------------------------------------------------
! 
! compute evaporation from canopy intercepted water
!
! definitions:
!    wdmax - max dew and rain holding capacity
!    fcwt - dew and rain trapped on vegetation
! ----------------------------------------------------------------------
!
   rs  = veg_rmin(vtype)
   rarc= veg_rarc(vtype)
   rgl = veg_rgl(vtype)
   mfrtf = 0.0
!
   do k = 1,msl
     mfsmc(k) = fsmc(k)
   enddo
!
   tmp_wdew = fcwt
   wdmax = lai_wf * fflai
!
   if (tmp_wdew.gt.wdmax) then
     mfrtf = tmp_wdew - wdmax
     tmp_wdew = wdmax
   end if
#ifdef DBGVIC
!
   print *,'in canopy_evap',fcwt,wdmax,tmp_wdew
#endif
!
! ----------------------------------------------------------------------
! in the original c-vic (uw), the negative rad is also used in noah_penman
! method, which will introduce errors in computing energy balance.
! therefore, we change the negative rad computation. if rad is negative,
! we use aerodynamic method to compute evap and assume canopy 
! temperature equals air temperature.
! ----------------------------------------------------------------------
!
! -- evap (m/s)
!
   if(rad.gt.0.0) then
     evap = vic_funct_penman(rad,    vpd,    ra,    rs,                        &
                             rarc,  fflai,   1.0,  tair,                       &
                           netsht,   pgcm,   rgl)
   else
     evap = rhoair*(0.622/pgcm)*(svp(tair)-vpair)/ra
     evap = evap / rhoh2o
     if(vpd.gt.0.0.and.evap.lt.0.0) evap = 0.0  ! no dew
     if (vpd.eq.0.0.and.evap.gt.0.0) evap=0.0
     if (evap.gt.0.) evap = 0.0       ! negative rad, no evap
   endif
!
   if(evap.ge.0.0)then
     if(tmp_wdew.gt.0) then
       bas  = tmp_wdew / wdmax
       frac = bas**(2.0/3.0) ! wet leaf fraction
     else
       frac = 0.0
     endif
     canopyevap = frac * evap * dtime            
   else
     frac = 1.0
     canopyevap = evap*dtime     ! dew
   endif
!     
! ----------------------------------------------------------------------
! evap can not exceed current storage
! ----------------------------------------------------------------------
!
   if (canopyevap .gt. 0.0) then
     f = min(1.0,(tmp_wdew / canopyevap))
   else
     f = 1.0
   endif
   canopyevap = canopyevap*f
!
   tmp_wdew = tmp_wdew + rfall - canopyevap
   if (tmp_wdew .gt. wdmax)then
     mfrtf  = mfrtf + tmp_wdew - wdmax
     tmp_wdew = wdmax
   endif
!
! ----------------------------------------------------------------------
! compute transpiration from vegetation
! ----------------------------------------------------------------------
   tmp_evap = canopyevap
!
! if evap .lt. 0, dewing occurs and the whole leaf wet and no vic_transpiration
!
#ifdef DBGVIC
   print 77,canopyevap,evap,rad,tmp_wdew,mfrtf
   77   format('canpevp=',e12.4,4(1x,e12.4))
#endif
!
! netsht radiation > 0 for activing transpiration processes (6/2004)
! ra/10.0 (6/2004)
!
   if(evap.ge.0.and.rad.ge.0.0.and.netsht.gt.0.0) then
     call vic_transpiration(msl,  dtime,    fsmc,    fsic,                     &
                    wcr,    wpw,    frt,       f,                              &
                   frac,    rad,    vpd,      ra,                              &
                     rs,   rarc,  fflai,  netsht,                              &
                   tair,   pgcm,    rgl,   mfsmc,                              &
                 layevp)

     do k = 1,msl
       tmp_evap = tmp_evap + layevp(k)
     end do
#ifdef DBGVIC
!
     print *,'vic_transpiration ', msl,  dtime,   fsmc,    fsic,               &
                    wcr,    wpw,    frt,       f,                              &
                   frac,    rad,    vpd,      ra,                              &
                     rs,   rarc,  fflai,  netsht,                              &
                   tair,   pgcm,    rgl,   mfsmc,                              &
                 layevp
#endif
   endif
!
   mfcwt = tmp_wdew
!
   canopy_evap = tmp_evap/dtime     ! m/time step -> m/s
!
   return
   end function canopy_evap
!
!-------------------------------------------------------------------------------
#include <define.h>
   subroutine vic_transpiration(msl, dtime,  fsmc,    fsic,                    &
                       wcr,   wpw,  root,       f,                             &
                      frac,   rad,   vpd,      ra,                             &
                        rs,  rarc, fflai,  netsht,                             &
                      tair,  pgcm,   rgl,   mfsmc,                             &
                    layevp)
!-------------------------------------------------------------------------------
!
! subprogram: vic_transpiration
!
! abstract: computes the traspiration
!
! comments:
!    please check equations (5), (7), (8), (9), (10), and (11) in liang, 
!    lettenmaier, wood, and burges (wrr 1994) for computing evap.
!
! program history log:
!   2003-06/07  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (canopy_evap.c)
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
!
   integer  ::  msl         ! number of soil layer
   real     ::  dtime       ! time step (second)
   real     ::  fsic(msl)   ! soil ice (mm)
   real     ::  wcr(msl)    ! ~70% field capacity (mm)
   real     ::  wpw(msl)    ! wilting point soil moisture (mm)
   real     ::  root(msl)   ! root content (fraction)
   real     ::  f           ! fraction of time step for no canopy water
   real     ::  frac        ! fraction of wet leaf
   real     ::  rad         ! net radiation (w/m2)
   real     ::  vpd         ! vapor pressure deficit (pa)
   real     ::  ra          ! vic_exch resistance (s/m)
   real     ::  rs          ! minimum stomatal resistance (s/m)
   real     ::  rarc        ! architectural resistance (s/m)
   real     ::  fflai       ! leaf area index
   real     ::  netsht      ! net shortwave radiation (w/m2)
   real     ::  tair        ! air temperature (k)
   real     ::  pgcm        ! pressure (pa)
   real     ::  rgl         ! value of solar radiation below which there 
                            ! will be no transpiration (ranges from ~30 w/m^2
                            ! for trees to ~100 w/m^2 for crops)
!
! ------------------- modified variables -------------------------------
!
   real     ::  fsmc(msl)   ! soil moisture (liquid + ice) (mm)
   real     ::  mfsmc(msl)  ! soil moisture (liquid + ice) (mm)
!
! ------------------- output variables ---------------------------------
!
   real     ::  layevp(msl) ! transpiration from each soil layer (m/time step)
!
! ---------------------- local variables -------------------------------
!
   integer  ::  k           ! loop index
   real     ::  gsm_inv     ! soil moisture stress factor
   real     ::  evap        ! tmp holding for evap total
   real     ::  moist1      ! tmp holding of moisture top
   real     ::  moist2      ! tmp holding of moisture bottom
   real     ::  wcr1        ! tmp holding of critical water top
   real     ::  wcr2        ! tmp holding of critical water bottom
   real     ::  root1       ! tmp holding of root top
   real     ::  root2       ! tmp holding of root bottom
   real     ::  rootsum     ! proportion of roots in moist>wcr zones
   real     ::  spare_evap  ! evap for 2nd distribution
   real     ::  avlsm(msl)  ! moisture available for trans (mm)

   real     ::  vic_funct_penman   ! function name
!
! ----------------------------------------------------------------------
! computes evapotranspiration for unfrozen soils. allows multiple layers
! compute moisture content in combined upper layers
! ----------------------------------------------------------------------
!
   moist1 = 0.0
   wcr1 = 0.0
   root1= 0.0
!
   do k = 1, msl-1
     if(root(k) .gt. 0.) then
       avlsm(k) = fsmc(k)-fsic(k)
       moist1= moist1 + avlsm(k)
       wcr1  = wcr1 + wcr(k)
       root1 = root1 + root(k)
     else
       avlsm(k) = 0.0
     end if
   end do
!
! ----------------------------------------------------------------------
! compute moisture content in lowest layer
! ----------------------------------------------------------------------
!
   k = msl
   moist2   = fsmc(k)-fsic(k)
   avlsm(k) = moist2
   wcr2  = wcr(k)
   root2 = root(k)
!
! ----------------------------------------------------------------------
! please refer the content (p14417-14418) of liang et al. (jgr 1994)
! for the following two cases of transpiration 
!
! case 1: moisture in both layers exceeds wcr, or moisture in
!     layer with more than half of the roots exceeds wcr.
!     potential evapotranspiration not hindered by soil dryness.  if
!     layer with less than half the roots is dryer than wcr, extra
!     evaporation is taken from the wetter layer. otherwise layers
!     contribute to evapotransipration based on root fraction.
! ----------------------------------------------------------------------
!
   if( (moist1.ge.wcr1.and.moist2.ge.wcr2.and.wcr1.gt.0.).or.                  &
       (moist1.ge.wcr1.and.root1.ge.0.5) .or.                                  &
       (moist2.ge.wcr2.and.root2.ge.0.5) ) then

     gsm_inv=1.0

! -- evap (m/s)

     evap = vic_funct_penman(rad,    vpd,      ra,      rs,                    &
                             rarc,  fflai, gsm_inv,    tair,                   &
                           netsht,   pgcm,     rgl)

     evap = evap*(1.0-f*frac)*dtime     ! m/s -> m/time step
!
! ----------------------------------------------------------------------
! divide up evap based on root distribution
! note the indexing of the roots 
! ----------------------------------------------------------------------
!
     rootsum = 1.0
     spare_evap = 0.0
!
     do k = 1,msl
       if(avlsm(k).ge.wcr(k))then
         layevp(k) = evap*root(k)
       else
         if (avlsm(k) .ge. wpw(k)) then 
           gsm_inv = (avlsm(k) - wpw(k)) /                                     &
                     (wcr(k) - wpw(k))
         else 
           gsm_inv = 0.0
         end if
         layevp(k)  = evap*gsm_inv*root(k)
         rootsum    = rootsum - root(k)
         spare_evap = evap*root(k)*(1.0-gsm_inv)
       end if
     end do
!
! ----------------------------------------------------------------------
! assign excess evaporation to wetter layer
! ----------------------------------------------------------------------
!
     if(spare_evap.gt.0.0) then
       do k = 1, msl
         if(avlsm(k).ge.wcr(k)) then
           layevp(k) = layevp(k)+root(k)*spare_evap/rootsum
         end if
       end do
     end if
#ifdef DBGVIC
!
      print 77,evap,layevp
      77  format('case 1 tra evap',e10.4,' layevp',3(e10.4,1x))
#endif
! ----------------------------------------------------------------------
! case 2: independent evapotranspirations
!    evapotranspiration is restricted by low soil moisture. evaporation
!    is computed independantly from each soil layer.
! ----------------------------------------------------------------------
   else 
     do k = 1,msl
! ----------------------------------------------------------------------
! set evaporation restriction factor 
! ----------------------------------------------------------------------
       if(avlsm(k) .ge. wcr(k)) then
         gsm_inv = 1.0
       else if(avlsm(k) .ge. wpw(k)) then
         gsm_inv = (avlsm(k) - wpw(k)) /                                       &
                      (wcr(k) - wpw(k))
       else 
         gsm_inv = 0.0
       end if
!
       if(gsm_inv .gt. 0.0 .and. root(k).gt.0.0) then
!
! ----------------------------------------------------------------------
! compute potential evapotranspiration
! ----------------------------------------------------------------------
!
         evap = vic_funct_penman(rad,   vpd,      ra,    rs,                   &
                                   rarc, fflai, gsm_inv,  tair,                &
                                 netsht,  pgcm,     rgl)

         layevp(k) = evap*root(k)*(1.0-f*frac)*dtime
#ifdef DBGVIC
         print 66,evap,k,layevp(k)
         66  format('case 2 evap',e10.4,' k',i4,' layevp',e10.4)
#endif
       else
         layevp(k) = 0.0
       end if
     end do
   end if
!
! ----------------------------------------------------------------------
! check that transpiration does not cause soil moisture to 
!       fall below wilting point.
! ----------------------------------------------------------------------
#ifdef DBGVIC
!
   print *,' bef fsmc-wpw layevp',layevp
#endif
!
   do k = 1,msl
     if(layevp(k)*1000.0 .gt. (fsmc(k)-wpw(k))) then
       layevp(k) = (fsmc(k) - wpw(k))/1000.0
     endif
     if (layevp(k) .lt. 0.0 ) then
       layevp(k) = 0.0
     end if
   end do
#ifdef DBGVIC
!
   print *,' aft fsmc-wpw layevp',layevp
#endif
!
   do k = 1, msl
     mfsmc(k) = fsmc(k) - layevp(k)*1000.0
   end do
!
   return
   end subroutine vic_transpiration
!
!-------------------------------------------------------------------------------
#ifdef VICLSM1
   function root_brent(ubound, lbound, ctrfct)
!-------------------------------------------------------------------------------
!
! abstract: this function returns the temperature for 
!  which the sum of the energy balance terms is zero in the 
!  interval [mint, maxt]. the temperature is calculated 
!  to within a tolerance (6 * macheps * |t| + 2 * t), where macheps 
!  is the relative machine precision and t is a positive tolerance, as 
!  specified.
!            the function assures that f(mintsurf) and f(maxtsurf) have 
!  opposite signs. if this is not the case the program will stop.  in 
!  addition the program will perform not more than a certain number of 
!  iterations, as specified in brent.h, and will abort if more iterations 
!  are needed.
!
! program history log:
!   2003-06    ji chen   modifed from 4.0.3 uw vic
!       (june ~ september)
!                 
!  general documentation for this module
!
!  source: brent, r. p., 1973, algorithms for minimization without
!   derivatives,prentice hall, inc., englewood cliffs, new jersey chapter 4
!
!  the method is also discussed in:
!  press, w. h., s. a. teukolsky, w. t. vetterling, b. p. flannery, 1992,
!  numerical recipes in fortran, the art of scientific computing,
!  second edition, cambridge university press
!  (be aware that this book discusses a brent method for minimization (brent), 
!  and one for root finding (zbrent).  the latter one is similar to the one 
!  implemented here and is also copied from brent [1973].)
!-------------------------------------------------------------------------------
#include <vartyp.h>
!
! ------------------- input variables ----------------------------------
!
   real ubound           ! upper bound for root
   real lbound           ! lower bound for root
!
   character*(*)  ctrfct ! determine function
!
! ---------------------- local variables -------------------------------
!
   real a, b, c, d, e, fa, fb, fc
   real m, p, q, r, s, tol
!
   real root_brent         ! function name
   real vic_funct_energy   ! function name
   real func_soiltemp      ! function name
   real func_snoweb        ! function name
!
   real tstep, tt
   parameter (tstep = 0.5, tt = 1.e-7)
!
   integer maxtries, maxiter
   parameter (maxtries=100, maxiter=1000)
   real macheps
   parameter (macheps=3.e-8)
!
   integer i, j            ! loop index
!
! ----------------------------------------------------------------------
! initialize variable argujment list
! ----------------------------------------------------------------------
!
   a = lbound
   b = ubound
!   
   if(ctrfct.eq.'vic_funct_energy')then
     fa = vic_funct_energy(a)
     fb = vic_funct_energy(b)
   elseif(ctrfct.eq.'func_soiltemp')then
     fa = func_soiltemp(a)
     fb = func_soiltemp(b)
   elseif(ctrfct.eq.'func_snoweb')then
     fa = func_snoweb(a)
     fb = func_snoweb(b)
   end if
!
! ----------------------------------------------------------------------
! if root not bracketed attempt to bracket the root
! ----------------------------------------------------------------------
!
   j = 0
!
   do while ((fa * fb).ge. 0  .and. j .lt. maxtries)
!
     a = a-tstep
     if(ctrfct.eq.'func_snoweb') then
       b = b
     else
       b = b+tstep
     endif
!
     if(ctrfct.eq.'vic_funct_energy')then
       fa = vic_funct_energy(a)
       fb = vic_funct_energy(b)
     elseif(ctrfct.eq.'func_soiltemp')then
       fa = func_soiltemp(a)
       fb = func_soiltemp(b)
     elseif(ctrfct.eq.'func_snoweb')then
       fa = func_snoweb(a)
     end if
!
      j = j + 1
   end do
!
   if ((fa * fb) .ge. 0) then
      write(6,*) 'error: first error in root_brent fa * fb >= 0'
      write(6,*) 'ctrfct ',ctrfct, ' a=',a,' b=',b
      write(6,*) ' fa=',fa,' fb=',fb
      call vic_stop_run
   end if
!  
   fc = fb
!   
   do i = 1,maxiter
!
     if(fb*fc .gt. 0) then
       c = a
       fc = fa
       d = b - a
       e = d
     end if
!    
     if (abs(fc) .lt. abs(fb)) then
       a = b
       b = c
       c = a
       fa = fb
       fb = fc
       fc = fa
     end if
!   
     tol = 2 * macheps * abs(b) + tt
     m = 0.5 * (c - b)
!    
     if (abs(m) .le. tol .or. fb.eq.0.0) then
       if(ctrfct.eq.'vic_funct_energy')then
         fb = vic_funct_energy(b)
       elseif(ctrfct.eq.'func_soiltemp')then
         fb = func_soiltemp(b)
       elseif(ctrfct.eq.'func_snoweb')then
         fb = func_snoweb(b)
       end if
       root_brent = b
       return
     else
       if (abs(e) .lt. tol .or. abs(fa) .le. abs(fb)) then
         d = m
         e = d
       else 
         s = fb/fa
         if (a.eq.c) then
! ----------------------------------------------------------------------
! linear interpolation           
! ----------------------------------------------------------------------
           p = 2.0 * m * s
           q = 1.0 - s
         else
! ----------------------------------------------------------------------
! inverse quadratic interpolation
! ----------------------------------------------------------------------
           q = fa/fc
           r = fb/fc
           p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
           q = (q - 1.0) * (r - 1.0) * (s - 1.0)
         end if
!   
         if (p .gt. 0.0) then
           q = -q
         else
           p = -p
         end if
         s = e
         e = d
         if(((2.*p).lt.(3.0*m*q-abs(tol*q))).and.                           &
           (p.lt.abs(0.5*s*q)))then
           d = p/q
         else 
           d = m
           e = d
         end if
       end if
       a = b
       fa = fb
       if(abs(d).gt.tol) then
         b = b + d
       else
         if(m.gt.0) then
           b = b + tol
         else
           b = b - tol
         end if
       end if
!
       if(ctrfct.eq.'vic_funct_energy')then
         fb = vic_funct_energy(b)
       elseif(ctrfct.eq.'func_soiltemp')then
         fb = func_soiltemp(b)
       elseif(ctrfct.eq.'func_snoweb')then
         fb = func_snoweb(b)
       end if
!
     end if
   end do
!
   write(6,*) 'error: second error in root_brent too many iterations'
   write(6,*) 'function ctrfct ',ctrfct(1:14), ' a=',a,' b=',b
   write(6,*) 'fa=',fa,' fb=',fb
   write(6,*) 'dumping input variables - check for valid values'
   call vic_stop_run
!
   return
   end function root_brent
#endif
!
!-------------------------------------------------------------------------------
   function svp(temp)
#include <vartyp.h>
!-------------------------------------------------------------------------------
!     
! subprogram: svp
!
! abstract: compute the saturated vapor pressure using 
!   handbook of hydrology eqn 4.2.2 pressure in pa
!
! program history log:
!   2003-06/07  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (svp.c)
!
!-------------------------------------------------------------------------------
!
! --------------- input variables from atmospheric model ---------------
!
   real temp        ! air temperature (k)
!
! ------------------- output variables ---------------------------------
!
   real svp         ! saturated vapor pressure (pa)
!
! --------------- local variables --------------------------------------
!
   real t           ! air temperature (c)
   real a_svp, b_svp, c_svp, t0c
   parameter (a_svp=0.61078, b_svp=17.269, c_svp=237.3, t0c=273.15)
!
!-------------------------------------------------------------------------------
!
   t = temp - t0c
!
   svp = a_svp * exp((b_svp * t)/(c_svp+t))
!
   if(t.lt.0) then
     svp = svp*(1.0 + .00972 * t + .000042 * t * t)
   end if
!
   svp = svp * 1000.0
!
   return
   end function svp
