#include <define.h>
   module lsm_noah_module
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [lsm_noah_main]
!      |
!    [lsm_noah_module]
!      |
!      |--- [noah_read_parameter] *
!      |--- [noah_fresh_snow] *
!      |--- [noah_snow_fraction] *
!      |--- [noah_soilt_conductivity] *
!      |--- [noah_snow_rough] *
!      |--- [noah_albedo] *
!      |--- [lsm_exch_coeff]  
!      |--- [noah_penman] *
!      |--- [noah_canopy_res] *
!      |--- [noah_bare_soil_solver]  --- [noah_soil_moisture_flux] *
!      |                              |     |- [noah_soilm_transpiration] *
!      |                              |     |- [noah_soilm_budget] *
!      |                              |             |- [noah_soilm_coeff] *  
!      |                              |     |- [noah_soilm_matrix] *
!      |                              |- [noah_soilt_conductivity]
!      |                              |- [noah_soil_heat_flux] *
!      |                                    |- [noah_soilt_ice] *
!      |                                    |- [noah_soilt_matrix] *  
!      |                                    |     |- [noah_tridiag_matrix] *  
!      |                                    |- [noah_soilt_land] *  
!      |                                          |- [noah_soilt_avg] *  
!      |                                          |- [noah_soilt_bdy] *  
!      |                            
!      |--- [noah_snow_cover_slover] --- [noah_soil_moisture_flux] *
!                                     |- [noah_soil_heat_flux] *
!                                     |- [noah_snowpack] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   contains
!
!-------------------------------------------------------------------------------
   subroutine noah_read_parameter (                                            &
                      vegtyp,soiltyp,slopetyp,                                 &
               cfactr,cmcmax,rsmax,topt,refkdt,kdt,sbeta,                      &
               shdfac,rsmin,rgl,hs,zbot,frzx,psisat,slope,                     &
               snup,salp,bexp,dksat,dwsat,                                     &
                      smcmax,smcwlt,smcref,                                    &
               smcdry,f1,quartz,fxexp,rtdis,sldpth,zsoil,                      &
               nroot,nsoil,z0,czil,lai,csoil,ptu)
! ----------------------------------------------------------------------
!
! subroutine noah_read_parameter
!
! internally set (default valuess), or optionally read-in via namelist
! i/o, all soil and vegetation parameters required for the execusion of
! the noah lsm.
!
! optional non-default parameters can be read in, accommodating up to 30
! soil, veg, or slope classes, if the default max number of soil, veg,
! and/or slope types is reset.
!
! future upgrades of routine noah_read_parameter must expand to 
! incorporate some of the empirical parameters of the frozen soil and 
! noah_snowpack physics (such as in routines ffrh2o, noah_snowpack, and 
! noah_fresh_snow) not yet set in this noah_read_parameter routine, but 
! rather set in lower level subroutines.
!
! set maximum number of soil-, veg-, and slopetyp in data statement.
!
! ----------------------------------------------------------------------
#include "abort.h"
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   integer,parameter    ::  max_slopetyp = 30
   integer,parameter    ::  max_soiltyp = 30
   integer,parameter    ::  max_vegtyp = 30
!
! number of defined soil-, veg-, and slopetyps used.
!
   integer              ::  defined_veg
   integer              ::  defined_soil
   integer              ::  defined_slope
!
   data defined_veg/13/
   data defined_soil/9/
   data defined_slope/9/
!
! ----------------------------------------------------------------------
!  set-up soil parameters for given soil type
!  input: soltyp: soil type (integer index)
!  output: soil parameters:
!    maxsmc: max soil moisture content (porosity)
!    refsmc: reference soil moisture (onset of soil moisture
!        stress in transpiration)
!    wltsmc: wilting pt soil moisture contents
!    drysmc: air dry soil moist content limits
!    satpsi: saturated soil potential
!    satdk:  saturated soil hydraulic conductivity
!    bb:     the 'b' parameter
!    satdw:  saturated soil diffusivity
!    f11:    used to compute soil diffusivity/conductivity
!    quartz:  soil quartz content
! ----------------------------------------------------------------------
! soil types   zobler (1986)     cosby et al (1984) (quartz cont.(1))
!  1       coarse         loamy sand    (0.82)
!  2       medium         silty clay loam    (0.10)
!  3       fine         light clay    (0.25)
!  4       coarse-medium     sandy loam    (0.60)
!  5       coarse-fine       sandy clay    (0.52)
!  6       medium-fine       clay loam     (0.35)
!  7       coarse-med-fine   sandy clay loam    (0.60)
!  8       organic         loam       (0.40)
!  9       glacial land ice  loamy sand    (na using 0.82)
! ----------------------------------------------------------------------
   real                 ::  bb(max_soiltyp)
   real                 ::  drysmc(max_soiltyp)
   real                 ::  f11(max_soiltyp)
   real                 ::  maxsmc(max_soiltyp)
   real                 ::  refsmc(max_soiltyp)
   real                 ::  satpsi(max_soiltyp)
   real                 ::  satdk(max_soiltyp)
   real                 ::  satdw(max_soiltyp)
   real                 ::  wltsmc(max_soiltyp)
   real                 ::  qtz(max_soiltyp)
!
   real                 ::  bexp
   real                 ::  dksat
   real                 ::  dwsat
   real                 ::  f1
   real                 ::  ptu
   real                 ::  quartz
   real                 ::  refsmc1
   real                 ::  smcdry
   real                 ::  smcmax
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  wltsmc1
! ----------------------------------------------------------------------
! soil texture-related arrays.
! ----------------------------------------------------------------------
   data maxsmc/0.421, 0.464, 0.468, 0.434, 0.406, 0.465,                       &
          0.404, 0.439, 0.421, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
   data satpsi/0.04, 0.62, 0.47, 0.14, 0.10, 0.26,                             &
          0.14, 0.36, 0.04, 0.00, 0.00, 0.00,                                  &
          0.00, 0.00, 0.00, 0.00, 0.00, 0.00,                                  &
          0.00, 0.00, 0.00, 0.00, 0.00, 0.00,                                  &
          0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
   data satdk /1.41e-5, 0.20e-5, 0.10e-5, 0.52e-5, 0.72e-5,                    &
          0.25e-5, 0.45e-5, 0.34e-5, 1.41e-5, 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00/
   data bb    /4.26,  8.72, 11.55, 4.74, 10.73,  8.17,                         &
          6.77,  5.25,  4.26, 0.00,  0.00,  0.00,                              &
          0.00,  0.00,  0.00, 0.00,  0.00,  0.00,                              &
          0.00,  0.00,  0.00, 0.00,  0.00,  0.00,                              &
          0.00,  0.00,  0.00, 0.00,  0.00,  0.00/
   data qtz   /0.82, 0.10, 0.25, 0.60, 0.52, 0.35,                             &
          0.60, 0.40, 0.82, 0.00, 0.00, 0.00,                                  &
          0.00, 0.00, 0.00, 0.00, 0.00, 0.00,                                  &
          0.00, 0.00, 0.00, 0.00, 0.00, 0.00,                                  &
          0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
! ----------------------------------------------------------------------
! the following 5 parameters are derived later in noah_read_parameter.f 
! from the soil data, and are just given here for reference and to force 
! static storage allocation. -dag lohmann, feb. 2001
! ----------------------------------------------------------------------
   data refsmc/0.283, 0.387, 0.412, 0.312, 0.338, 0.382,                       &
          0.315, 0.329, 0.283, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
   data wltsmc/0.029, 0.119, 0.139, 0.047, 0.100, 0.103,                       &
          0.069, 0.066, 0.029, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
   data drysmc/0.029, 0.119, 0.139, 0.047, 0.100, 0.103,                       &
          0.069, 0.066, 0.029, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                            &
          0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
   data satdw /5.71e-6, 2.33e-5, 1.16e-5, 7.95e-6, 1.90e-5,                    &
          1.14e-5, 1.06e-5, 1.46e-5, 5.71e-6, 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00,                            &
          0.00   , 0.00   , 0.00   , 0.00   , 0.00/
   data f11  /-0.999, -1.116, -2.137, -0.572, -3.201, -1.302,                  &
         -1.519, -0.329, -0.999,  0.000,  0.000,  0.000,                       &
          0.000,  0.000,  0.000,  0.000,  0.000,  0.000,                       &
          0.000,  0.000,  0.000,  0.000,  0.000,  0.000,                       &
          0.000,  0.000,  0.000,  0.000,  0.000,  0.000/
! ----------------------------------------------------------------------
! set-up vegetation parameters for a given vegetaion type:
! input: vegtyp = vegetation type (integer index)
! ouput: vegetation parameters
!   shdfac: vegetation greenness fraction
!   rsmin:  mimimum stomatal resistance
!   rgl:    parameter used in solar rad term of
!       canopy resistance function
!   hs:     parameter used in vapor pressure deficit term of
!       canopy resistance function
!   snup:   threshold snow depth (in water equivalent m) that
!          implies 100% snow cover
! ----------------------------------------------------------------------
! ssib vegetation types (dorman and sellers, 1989; jam)
!  1:  broadleaf-evergreen trees  (tropical forest)
!  2:  broadleaf-deciduous trees
!  3:  broadleaf and needleleaf trees (mixed forest)
!  4:  needleleaf-evergreen trees
!  5:  needleleaf-deciduous trees (larch)
!  6:  broadleaf trees with groundcover (savanna)
!  7:  groundcover only (perennial)
!  8:  broadleaf shrubs with perennial groundcover
!  9:  broadleaf shrubs with bare soil
! 10:  dwarf trees and shrubs with groundcover (tundra)
! 11:  bare soil
! 12:  cultivations (the same parameters as for type 7)
! 13:  glacial (the same parameters as for type 11)
! ----------------------------------------------------------------------
   integer              ::  nroot
   integer              ::  nroot_data(max_vegtyp)
!
   real                 ::  frzfact
   real                 ::  hs
   real                 ::  hstbl(max_vegtyp)
   real                 ::  lai
   real                 ::  lai_data(max_vegtyp)
   real                 ::  psisat
   real                 ::  rsmin
   real                 ::  rgl
   real                 ::  rgltbl(max_vegtyp)
   real                 ::  rsmtbl(max_vegtyp)
   real                 ::  shdfac
   real                 ::  snup
   real                 ::  snupx(max_vegtyp)
   real                 ::  z0
   real                 ::  z0_data(max_vegtyp)
! ----------------------------------------------------------------------
! vegetation class-related arrays
! ----------------------------------------------------------------------
   data nroot_data /4,4,4,4,4,4,3,3,3,2,3,3,2,0,0,                             &
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

   data rsmtbl /150.0, 100.0, 125.0, 150.0, 100.0, 70.0,                       &
            40.0, 300.0, 400.0, 150.0, 400.0, 40.0,                            &
           150.0,   0.0,   0.0,   0.0,   0.0,  0.0,                            &
             0.0,   0.0,   0.0,   0.0,   0.0,  0.0,                            &
             0.0,   0.0,   0.0,   0.0,   0.0,  0.0/
   data rgltbl /30.0,  30.0,  30.0,  30.0,  30.0,  65.0,                       &
          100.0, 100.0, 100.0, 100.0, 100.0, 100.0,                            &
          100.0,   0.0,   0.0,   0.0,   0.0,   0.0,                            &
            0.0,   0.0,   0.0,   0.0,   0.0,   0.0,                            &
            0.0,   0.0,   0.0,   0.0,   0.0,   0.0/
   data hstbl /41.69, 54.53, 51.93, 47.35,  47.35, 54.53,                      &
          36.35, 42.00, 42.00, 42.00,  42.00, 36.35,                           &
          42.00,  0.00,  0.00,  0.00,   0.00,  0.00,                           &
           0.00,  0.00,  0.00,  0.00,   0.00,  0.00,                           &
           0.00,  0.00,  0.00,  0.00,   0.00,  0.00/
   data snupx  /0.080, 0.080, 0.080, 0.080, 0.080, 0.080,                      &
           0.040, 0.040, 0.040, 0.040, 0.025, 0.040,                           &
           0.025, 0.000, 0.000, 0.000, 0.000, 0.000,                           &
           0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                           &
           0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
   data z0_data /2.653, 0.826, 0.563, 1.089, 0.854, 0.856,                     &
            0.035, 0.238, 0.065, 0.076, 0.011, 0.035,                          &
            0.011, 0.000, 0.000, 0.000, 0.000, 0.000,                          &
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,                          &
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
   data lai_data /4.0, 4.0, 4.0, 4.0, 4.0, 4.0,                                &
             4.0, 4.0, 4.0, 4.0, 4.0, 4.0,                                     &
             4.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                     &
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                     &
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
! ----------------------------------------------------------------------
! class parameter 'slopetyp' was included to estimate linear reservoir
! coefficient 'slope' to the baseflow runoff out of the bottom layer.
! lowest class (slopetyp=0) means highest slope parameter = 1.
! definition of slopetyp from 'zobler' slope type:
! slope class  percent slope
! 1          0-8
! 2          8-30
! 3          > 30
! 4          0-30
! 5          0-8 & > 30
! 6          8-30 & > 30
! 7          0-8, 8-30, > 30
! 9          glacial ice
! blank        ocean/sea
! ----------------------------------------------------------------------
! note:
! class 9 from 'zobler' file should be replaced by 8 and 'blank' 9
! ----------------------------------------------------------------------
   real slope
   real slope_data(max_slopetyp)
!
   data slope_data /0.1,  0.6, 1.0, 0.35, 0.55, 0.8,                           &
               0.63, 0.0, 0.0, 0.0,  0.0,  0.0,                                &
               0.0 , 0.0, 0.0, 0.0,  0.0,  0.0,                                &
               0.0 , 0.0, 0.0, 0.0,  0.0,  0.0,                                &
               0.0 , 0.0, 0.0, 0.0,  0.0,  0.0/
! ----------------------------------------------------------------------
! set namelist file name
! ----------------------------------------------------------------------
   character(len=50)    ::  namelist_name
! ----------------------------------------------------------------------
! set universal parameters (not dependent on soil, veg, slope type)
! ----------------------------------------------------------------------
   integer              ::  i
   integer              ::  nsoil
   integer              ::  slopetyp
   integer              ::  soiltyp
   integer              ::  vegtyp
!
   integer              ::  bare
   data bare /11/
!
   logical              ::  lparam
   data lparam /.true./
!
   logical,save         ::  lfirst
   data lfirst /.true./
!
! ** clu_rev: retain definition status to quarantee 'one-time' execution
!
! ** clu_rev: change the default value of czil from 0.2 to 0.1
! ----------------------------------------------------------------------
! parameter used to calculate roughness length of heat.
! ----------------------------------------------------------------------
   real                 ::  czil
   real                 ::  czil_data
!clu  data czil_data /0.2/
   data czil_data /0.1/
! ----------------------------------------------------------------------
! parameter used to caluculate vegetation effect on soil heat flux.
! ----------------------------------------------------------------------
   real                 ::  sbeta
   real                 ::  sbeta_data
   data sbeta_data /-2.0/
! ----------------------------------------------------------------------
! bare soil evaporation exponent used in ddevap.
! ----------------------------------------------------------------------
   real                 ::  fxexp
   real                 ::  fxexp_data
   data fxexp_data /2.0/
! ----------------------------------------------------------------------
! soil heat capacity [j m-3 k-1]
! ----------------------------------------------------------------------
   real                 ::  csoil
   real                 ::  csoil_data
!     data csoil_data /1.26e+6/
   data csoil_data /2.00e+6/
! ----------------------------------------------------------------------
! specify snow distribution shape parameter salp - shape parameter of
! distribution function of snow cover. from anderson's data (hydro-17)
! best fit is when salp = 2.6
! ----------------------------------------------------------------------
   real                 ::  salp
   real                 ::  salp_data
   data salp_data /2.6/
! ----------------------------------------------------------------------
! kdt is defined by reference refkdt and dksat; refdk=2.e-6 is the sat.
! dk. value for the soil type 2
! ----------------------------------------------------------------------
   real                 ::  refdk
   real                 ::  refdk_data
   data refdk_data /2.0e-6/
!
   real                 ::  refkdt
   real                 ::  refkdt_data
   data refkdt_data /3.0/
!
   real                 ::  frzx
   real                 ::  kdt
! ----------------------------------------------------------------------
! frozen ground parameter, frzk, definition: ice content threshold above
! which frozen soil is impermeable reference value of this parameter for
! the light clay soil (type=3) frzk = 0.15 m.
! ----------------------------------------------------------------------
   real                 ::  frzk
   real                 ::  frzk_data
   data frzk_data /0.15/
!
   real                 ::  rtdis(nsoil)
   real                 ::  sldpth(nsoil)
   real                 ::  zsoil(nsoil)
! ----------------------------------------------------------------------
! set two canopy water parameters.
! ----------------------------------------------------------------------
   real                 ::  cfactr
   real                 ::  cfactr_data
   data cfactr_data /0.5/
!
   real                 ::  cmcmax
   real                 ::  cmcmax_data
   data cmcmax_data /0.5e-3/
! ----------------------------------------------------------------------
! set max. stomatal resistance.
! ----------------------------------------------------------------------
   real                 ::  rsmax
   real                 ::  rsmax_data
   data rsmax_data /5000.0/
! ----------------------------------------------------------------------
! set optimum transpiration air temperature.
! ----------------------------------------------------------------------
   real                 ::  topt
   real                 ::  topt_data
   data topt_data /298.0/
! ----------------------------------------------------------------------
! specify depth[m] of lower boundary soil temperature.
! ----------------------------------------------------------------------
   real                 ::  zbot
   real                 ::  zbot_data
!      data zbot_data /-3.0/
   data zbot_data /-8.0/
! ----------------------------------------------------------------------
! set two soil moisture wilt, soil moisture reference parameters
! ----------------------------------------------------------------------
   real                 ::  smlow
   real                 ::  smlow_data
   data smlow_data /0.5/
!
   real                 ::  smhigh
   real                 ::  smhigh_data
   data smhigh_data /3.0/
! ----------------------------------------------------------------------
! namelist definition:
! ----------------------------------------------------------------------
   namelist /soil_veg/ slope_data, rsmtbl, rgltbl, hstbl, snupx,               &
      bb, drysmc, f11, maxsmc, refsmc, satpsi, satdk, satdw,                   &
      wltsmc, qtz, lparam, zbot_data, salp_data, cfactr_data,                  &
      cmcmax_data, sbeta_data, rsmax_data, topt_data,                          &
      refdk_data, frzk_data, bare, defined_veg, defined_soil,                  &
      defined_slope, fxexp_data, nroot_data, refkdt_data, z0_data,             &
      czil_data, lai_data, csoil_data
! ----------------------------------------------------------------------
! read namelist file to override default parameters only once.
! namelist_name must be 50 characters or less.
! ----------------------------------------------------------------------
   if (lfirst) then
! ** clu_rev: hardwire the namelist*!
!        namelist_name = 'soil_veg_namelist_ver_2.5.1'
!
     lfirst = .false.
!
     if (defined_soil .gt. max_soiltyp) then
       write(6,*) 'warning: defined_soil too large in namelist'
       call MPABORT
     endif
     if (defined_veg .gt. max_vegtyp) then
       write(6,*) 'warning: defined_veg too large in namelist'
       call MPABORT
     endif
     if (defined_slope .gt. max_slopetyp) then
       write(6,*) 'warning: defined_slope too large in namelist'
       call MPABORT
     endif
!
     if(vegtyp.eq.0) then
       write(6,*)'vegtyp.eq.0'
     endif
     if(soiltyp.eq.0) then
       write(6,*)'soiltyp.eq.0'
     endif
     if(slopetyp.eq.0) then
       write(6,*)'slopetyp.eq.0'
     endif
!
     smlow = smlow_data
     smhigh = smhigh_data
!
     do i = 1,defined_soil
       satdw(i)  = bb(i)*satdk(i)*(satpsi(i)/maxsmc(i))
       f11(i) = alog10(satpsi(i)) + bb(i)*alog10(maxsmc(i)) + 2.0
       refsmc1 = maxsmc(i)*(5.79e-9/satdk(i))                                  &
               **(1.0/(2.0*bb(i)+3.0))
!      refsmc(i) = refsmc1 + (maxsmc(i)-refsmc1) / 3.0
       refsmc(i) = refsmc1 + (maxsmc(i)-refsmc1) / smhigh
       wltsmc1 = maxsmc(i) * (200.0/satpsi(i))**(-1.0/bb(i))
!      wltsmc(i) = wltsmc1 - 0.5 * wltsmc1
       wltsmc(i) = wltsmc1 - smlow * wltsmc1
! ----------------------------------------------------------------------
! current version drysmc values that equate to wltsmc.
! future version could let drysmc be independently set via namelist.
! ----------------------------------------------------------------------
       drysmc(i) = wltsmc(i)
     end do

! ----------------------------------------------------------------------
! end lfirst block
! ----------------------------------------------------------------------
   endif
!
   if (soiltyp .gt. defined_soil) then
     write(6,*) 'warning: too many soil types'
     call MPABORT
   endif
   if (vegtyp .gt. defined_veg) then
     write(6,*) 'warning: too many veg types'
     call MPABORT
   endif
   if (slopetyp .gt. defined_slope) then
     write(6,*) 'warning: too many slope types'
     write(6,*) 'slopetyp,defined_slope=',slopetyp,defined_slope
     call MPABORT
   endif
! ----------------------------------------------------------------------
! set-up universal parameters (not dependent on soiltyp, vegtyp or
! slopetyp)
! ----------------------------------------------------------------------
   zbot = zbot_data
   salp = salp_data
   cfactr = cfactr_data
   cmcmax = cmcmax_data
   sbeta = sbeta_data
   rsmax = rsmax_data
   topt = topt_data
   refdk = refdk_data
   frzk = frzk_data
   fxexp = fxexp_data
   refkdt = refkdt_data
   czil = czil_data
   csoil = csoil_data
! ----------------------------------------------------------------------
!  set-up soil parameters
! ----------------------------------------------------------------------
   bexp = bb(soiltyp)
   dksat = satdk(soiltyp)
   dwsat = satdw(soiltyp)
   f1 = f11(soiltyp)
!     frzfact = (smcmax / smcref) * (0.412 / 0.468)
   kdt = refkdt * dksat/refdk
   psisat = satpsi(soiltyp)
   quartz = qtz(soiltyp)
   smcdry = drysmc(soiltyp)
   smcmax = maxsmc(soiltyp)
   smcref = refsmc(soiltyp)
   smcwlt = wltsmc(soiltyp)
   frzfact = (smcmax / smcref) * (0.412 / 0.468)
! ----------------------------------------------------------------------
! to adjust frzk parameter to actual soil type: frzk * frzfact
! ----------------------------------------------------------------------
   frzx = frzk * frzfact
! ----------------------------------------------------------------------
! set-up vegetation parameters
! ----------------------------------------------------------------------
   nroot = nroot_data(vegtyp)
!
!  force noah to run on 2 layers
!
   if(nroot .gt. nsoil) nroot = nsoil
!
   snup = snupx(vegtyp)
   rsmin = rsmtbl(vegtyp)
   rgl = rgltbl(vegtyp)
   hs = hstbl(vegtyp)
   z0 = z0_data(vegtyp)
   lai = lai_data(vegtyp)
   if (vegtyp .eq. bare) shdfac = 0.0
!
   if (nroot .gt. nsoil) then
     write(6,*) 'warning: too many root layers'
     call MPABORT
   endif
! ----------------------------------------------------------------------
! calculate root distribution.  present version assumes uniform
! distribution based on soil layer depths.
! ----------------------------------------------------------------------
   do i = 1,nsoil
     rtdis(i) = -sldpth(i)/zsoil(nroot)
   end do
! ----------------------------------------------------------------------
!  set-up slope parameter
! ----------------------------------------------------------------------
   slope = slope_data(slopetyp)
!
   return
   end subroutine noah_read_parameter
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_fresh_snow (temp,newsn,snowh,sndens)
!-------------------------------------------------------------------------------
! 
! calculate snow depth and densitity to account for the new snowfall.
! new values of snow depth & density returned.
!
! temp    air temperature (k)
! newsn   new snowfall (m)
! snowh   snow depth (m)
! sndens  snow density (g/cm3=dimensionless fraction of h2o density)
!
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   real sndens
   real dsnew
   real snowhc
   real hnewc
   real snowh
   real newsn
   real newsnc
   real temp 
   real tempc
!
! conversion into simulation units      
!
   snowhc = snowh*100.
   newsnc = newsn*100.
   tempc = temp-273.15
!
! calculating new snowfall density depending on temperature
! equation from gottlib l. 'a general runoff model for snowcovered
! and glacierized basin', 6th nordic hydrological conference,
! vemadolen, sweden, 1980, 172-177pp.
!
   if (tempc .le. -15.) then
     dsnew = 0.05
   else                                                      
     dsnew = 0.05+0.0017*(tempc+15.)**1.5
   endif
!
! adjustment of snow density depending on new snowfall      
!
   hnewc = newsnc/dsnew
   sndens = (snowhc*sndens+hnewc*dsnew)/(snowhc+hnewc)
   snowhc = snowhc+hnewc
   snowh = snowhc*0.01
!
   return
   end subroutine noah_fresh_snow
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_snow_fraction (sneqv,snup,salp,snowh,sncovr)
!-------------------------------------------------------------------------------
!      
! calculate snow fraction (0 -> 1)
! sneqv   snow water equivalent (m)
! snup    threshold sneqv depth above which sncovr=1
! salp    tuning parameter
! sncovr  fractional snow cover
!
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   real                 ::  sneqv, snup, salp, sncovr, rsnow, z0n, snowh
! ----------------------------------------------------------------------
! snup is veg-class dependent snowdepth threshhold (set in routine
! noah_read_parameter) above which snocvr=1.
! ----------------------------------------------------------------------
   if (sneqv .lt. snup) then
     rsnow = sneqv/snup
     sncovr = 1. - ( exp(-salp*rsnow) - rsnow*exp(-salp))
   else
     sncovr = 1.0
   endif
!
   z0n=0.035 
!     formulation of dickinson et al. 1986
!        sncovr=snowh/(snowh + 5*z0n)
!     formulation of marshall et al. 1994
!        sncovr=sneqv/(sneqv + 2*z0n)
!
   return
   end subroutine noah_snow_fraction
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilt_conductivity ( df, smc, qz,  smcmax, sh2o)
!-------------------------------------------------------------------------------
!
! calculate thermal diffusivity and conductivity of the soil for a given
! point and time.
!
! peters-lidard approach (peters-lidard et al., 1998)
! june 2001 changes: frozen soil condition.
!
!
! we now get quartz as an input argument (set in routine noah_read_parameter):
!      data quartz /0.82, 0.10, 0.25, 0.60, 0.52, 
!     &             0.35, 0.60, 0.40, 0.82/
!
! if the soil has any moisture content compute a partial sum/product
! otherwise use a constant value which works well with most soils
! ----------------------------------------------------------------------
!  thkw ......water thermal conductivity
!  thkqtz ....thermal conductivity for quartz
!  thko ......thermal conductivity for other soil components
!  thks ......thermal conductivity for the solids combined(quartz+other)
!  thkice ....ice thermal conductivity
!  smcmax ....porosity (= smcmax)
!  qz .........quartz content (soil type dependent)
! ----------------------------------------------------------------------
! use as in peters-lidard, 1998 (modif. from johansen, 1975).
!
!                                  pablo grunmann, 08/17/98
! refs.:
!      farouki, o.t.,1986: thermal properties of soils. series on rock 
!              and soil mechanics, vol. 11, trans tech, 136 pp.
!      johansen, o., 1975: thermal conductivity of soils. ph.d. thesis,
!              university of trondheim,
!      peters-lidard, c. d., et al., 1998: the effect of soil thermal 
!              conductivity parameterization on surface energy fluxes
!              and temperatures. journal of the atmospheric sciences,
!              vol. 55, pp. 1209-1224.
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
    real                ::  df
    real                ::  gammd
    real                ::  thkdry
    real                ::  ake
    real                ::  thkice
    real                ::  thko
    real                ::  thkqtz
    real                ::  thksat
    real                ::  thks
    real                ::  thkw
    real                ::  qz
    real                ::  satratio
    real                ::  sh2o
    real                ::  smc
    real                ::  smcmax
    real                ::  xu
    real                ::  xunfroz
! ----------------------------------------------------------------------
!
! needs parameters
! porosity(soil type):
!      poros = smcmax
!
! saturation ratio:
!
   satratio = smc/smcmax
!
! parameters  w/(m.k)
!
   thkice = 2.2
   thkw = 0.57
   thko = 2.0
!      if (qz .le. 0.2) thko = 3.0
   thkqtz = 7.7
!
! solids' conductivity      
!
   thks = (thkqtz**qz)*(thko**(1.- qz))
!
! unfrozen fraction (from 1., i.e., 100%liquid, to 0. (100% frozen))
!
   xunfroz = sh2o /smc
!
! unfrozen volume for saturation (porosity*xunfroz)
!
   xu=xunfroz*smcmax 
!
! saturated thermal conductivity
!
   thksat = thks**(1.-smcmax)*thkice**(smcmax-xu)*thkw**(xu)
!
! dry density in kg/m3
!
   gammd = (1. - smcmax)*2700.
!
! dry thermal conductivity in w.m-1.k-1
!
   thkdry = (0.135*gammd + 64.7)/(2700. - 0.947*gammd)
   if ( (sh2o + 0.0005) .lt. smc ) then
!
! frozen
!
     ake = satratio
   else
!
! unfrozen
! range of validity for the kersten number (ake)
!
     if ( satratio .gt. 0.1 ) then
!
! kersten number (using "fine" formula, valid for soils containing at 
! least 5% of particles with diameter less than 2.e-6 meters.)
! (for "coarse" formula, see peters-lidard et al., 1998).
!
       ake = log10(satratio) + 1.0
     else
!
! use k = kdry
!
       ake = 0.0
     endif
   endif
!
!  thermal conductivity
!
   df = ake*(thksat - thkdry) + thkdry
!
   return
   end subroutine noah_soilt_conductivity
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_snow_rough (sncovr,z0)
!-------------------------------------------------------------------------------
!
! calculate total roughness length over snow
! sncovr  fractional snow cover
! z0      roughness length (m)
! z0s     snow roughness length:=0.001 (m)
!
! current noah lsm condition - mbek, 09-oct-2001
!
!-------------------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   real                 :: sncovr, z0, z0s
! ----------------------------------------------------------------------
   z0s = z0
!
   z0 = (1-sncovr)*z0 + sncovr*z0s
!
   return
   end subroutine noah_snow_rough
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_albedo (alb,snoalb,shdfac,shdmin,sncovr,tsnow,albedo)
!-------------------------------------------------------------------------------
!      
! calculate albedo including snow effect (0 -> 1)
!   alb     snowfree albedo
!   snoalb  maximum (deep) snow albedo
!   shdfac    areal fractional coverage of green vegetation
!   shdmin    minimum areal fractional coverage of green vegetation
!   sncovr  fractional snow cover
!   albedo  surface albedo including snow effect
!   tsnow   snow surface temperature (k)
!      
! snoalb is argument representing maximum albedo over deep snow,
! as passed into phys_lsm_noah, and adapted from the satellite-based maximum 
! snow albedo fields provided by d. robinson and g. kukla 
! (1985, jcam, vol 24, 402-411)
!
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   real                 ::  alb, snoalb, shdfac, shdmin, sncovr, albedo, tsnow
! ----------------------------------------------------------------------
   albedo = alb + (1.0-(shdfac-shdmin))*sncovr*(snoalb-alb) 
   if (albedo .gt. snoalb) albedo=snoalb
!
!     base formulation (dickinson et al., 1986, cogley et al., 1990)
!          if (tsnow.le.263.16) then
!            albedo=snoalb
!          else
!            if (tsnow.lt.273.16) then
!              tm=0.1*(tsnow-263.16)
!              albedo=0.5*((0.9-0.2*(tm**3))+(0.8-0.16*(tm**3)))
!            else
!              albedo=0.67
!            endif
!          endif
!
!     isba formulation (verseghy, 1991; baker et al., 1990)
!          if (tsnow.lt.273.16) then
!            albedo=snoalb-0.008*dt/86400
!          else
!            albedo=(snoalb-0.5)*exp(-0.24*dt/86400)+0.5
!          endif
!
   return
   end subroutine noah_albedo
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,t24,ssoil,      &
                           q2,q2sat,etp,rch,epsca,rr,snowng,frzgra,            &
                           dqsdt2,flx2)
!-------------------------------------------------------------------------------
!
! calculate potential evaporation for the current point.  various
! partial sums/products are also calculated and passed back to the
! calling routine for later use.
!
!-------------------------------------------------------------------------------
   use constant, only : cp=>cp_,r=>rd_
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   logical              ::  snowng
   logical              ::  frzgra
!
   real                 ::  a
   real                 ::  beta
   real                 ::  ch
   real                 ::  delta
   real                 ::  dqsdt2
   real                 ::  epsca
   real                 ::  etp
   real                 ::  fdown
   real                 ::  flx2
   real                 ::  fnet
   real                 ::  prcp
   real                 ::  q2
   real                 ::  q2sat
   real                 ::  rad
   real                 ::  rch
   real                 ::  rho
   real                 ::  rr
   real                 ::  ssoil
   real                 ::  sfcprs
   real                 ::  sfctmp
   real                 ::  t24
   real                 ::  t2v
   real                 ::  th2
!
   real,parameter       ::  cph2o = 4.218e+3
   real,parameter       ::  cpice = 2.106e+3
   real,parameter       ::  elcp = 2.4888e+3
   real,parameter       ::  lsubf = 3.335e+5
   real,parameter       ::  lsubc = 2.501000e+6
   real,parameter       ::  sigma = 5.67e-8
!
! executable code begins here:
!
   flx2 = 0.0
!
! prepare partial quantities for penman equation.
!
   delta = elcp * dqsdt2
   t24 = sfctmp * sfctmp * sfctmp * sfctmp
   rr = t24 * 6.48e-8 /(sfcprs * ch) + 1.0
   rho = sfcprs / (r * t2v)
   rch = rho * cp * ch
!
! adjust the partial sums / products with the latent heat
! effects caused by falling precipitation.
!
   if (.not. snowng) then
     if (prcp .gt. 0.0) rr = rr + cph2o*prcp/rch
   else
     rr = rr + cpice*prcp/rch
   endif
   fnet = fdown - sigma*t24 - ssoil
!
! include the latent heat effects of frzng rain converting to ice on
! impact in the calculation of flx2 and fnet.
!
   if (frzgra) then
     flx2 = -lsubf * prcp
     fnet = fnet - flx2
   endif
!
! finish penman equation calculations.
!
   rad = fnet/rch + th2 - sfctmp
   a = elcp * (q2sat - q2)
   epsca = (a*rr + rad*delta) / (delta + rr)
   etp = epsca * rch / lsubc
!
   return
   end subroutine noah_penman
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_canopy_res (solar,ch,sfctmp,q2,sfcprs,smc,zsoil,nsoil,      &
                               smcwlt,smcref,rsmin,rc,pc,nroot,q2sat,dqsdt2,   &
                               topt,rsmax,rgl,hs,xlai,                         &
                               rcs,rct,rcq,rcsoil)
!-------------------------------------------------------------------------------
!
! calculate canopy resistance which depends on incoming solar radiation,
! air temperature, atmospheric water vapor pressure deficit at the
! lowest model level, and soil moisture (preferably unfrozen soil
! moisture rather than total)
!
! source:  jarvis (1976), noilhan and planton (1989, mwr), jacquemin and
! noilhan (1990, blm)
! see also:  chen et al (1996, jgr, vol 101(d3), 7251-7268), eqns 12-14
! and table 2 of sec. 3.1.2         
!
! input:
!   solar   incoming solar radiation
!   ch      surface exchange coefficient for heat and moisture
!   sfctmp  air temperature at 1st level above ground
!   q2      air humidity at 1st level above ground
!   q2sat   saturation air humidity at 1st level above ground
!   dqsdt2  slope of saturation humidity function wrt temp
!   sfcprs  surface pressure
!   smc     volumetric soil moisture 
!   zsoil   soil depth (negative sign, as it is below ground)
!   nsoil   no. of soil layers
!   nroot   no. of soil layers in root zone (1.le.nroot.le.nsoil)
!   xlai    leaf area index
!   smcwlt  wilting point
!   smcref  reference soil moisture (where soil water deficit stress
!             sets in)
! rsmin, rsmax, topt, rgl, hs are canopy stress parameters set in
!   surboutine noah_read_parameter
!
! output:
!   pc  plant coefficient
!   rc  canopy resistance
!
!-------------------------------------------------------------------------------
   use constant, only : cp => cp_,rd => rd_
!
   implicit none
!
   integer,parameter    ::  nsold = 20
!
   integer              ::  k
   integer              ::  nroot
   integer              ::  nsoil
!
   real                 ::  ch
   real                 ::  delta
   real                 ::  dqsdt2
   real                 ::  ff
   real                 ::  gx
   real                 ::  hs
   real                 ::  p
   real                 ::  part(nsold) 
   real                 ::  pc
   real                 ::  q2
   real                 ::  q2sat
   real                 ::  rc
   real                 ::  rsmin
   real                 ::  rcq
   real                 ::  rcs
   real                 ::  rcsoil
   real                 ::  rct
   real                 ::  rgl
   real                 ::  rr
   real                 ::  rsmax
   real                 ::  sfcprs
   real                 ::  sfctmp
   real                 ::  smc(nsoil)
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  solar
   real                 ::  topt
   real                 ::  slvcp
   real                 ::  st1
   real                 ::  tair4
   real                 ::  xlai
   real                 ::  zsoil(nsoil)
!
   real,parameter       ::  sigma = 5.67e-8
   real,parameter       ::  slv=2.501000e6
!-------------------------------------------------------------------------------
!
! initialize canopy resistance multiplier terms.
!
   rcs = 0.0
   rct = 0.0
   rcq = 0.0
   rcsoil = 0.0
   rc = 0.0
!
! contribution due to incoming solar radiation
!
   ff = 0.55*2.0*solar/(rgl*xlai)
   rcs = (ff + rsmin/rsmax) / (1.0 + ff)
   rcs = max(rcs,0.0001)
!
! contribution due to air temperature at first model level above ground
! rct expression from noilhan and planton (1989, mwr).
!
   rct = 1.0 - 0.0016*((topt-sfctmp)**2.0)
   rct = max(rct,0.0001)
!
! contribution due to vapor pressure deficit at first model level.
! rcq expression from ssib 
!
   rcq = 1.0/(1.0+hs*(q2sat-q2))
   rcq = max(rcq,0.01)
!
! contribution due to soil moisture availability.
! determine contribution from each soil layer, then add them up.
!
   gx = (smc(1) - smcwlt) / (smcref - smcwlt)
   if (gx .gt. 1.) gx = 1.
   if (gx .lt. 0.) gx = 0.
!
! use soil depth as weighting factor
!
   part(1) = (zsoil(1)/zsoil(nroot)) * gx
!
! use root distribution as weighting factor
!      part(1) = rtdis(1) * gx
!
   do k = 2,nroot
     gx = (smc(k) - smcwlt) / (smcref - smcwlt)
     if (gx .gt. 1.) gx = 1.
     if (gx .lt. 0.) gx = 0.
!
! use soil depth as weighting factor        
!
     part(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot)) * gx
!
! use root distribution as weighting factor
!        part(k) = rtdis(k) * gx 
!
   end do
   do k = 1,nroot
     rcsoil = rcsoil+part(k)
   end do
   rcsoil = max(rcsoil,0.0001)
!
! determine canopy resistance due to all factors.  convert canopy
! resistance (rc) to plant coefficient (pc) to be used with potential
! evap in determining actual evap.  pc is determined by:
!   pc * linerized penman potential evap =
!   penman-monteith actual evaporation (containing rc term).
!
   rc = rsmin/(xlai*rcs*rct*rcq*rcsoil)
!      tair4 = sfctmp**4.
!      st1 = (4.*sigma*rd)/cp
!      slvcp = slv/cp
!      rr = st1*tair4/(sfcprs*ch) + 1.0
   rr = (4.*sigma*rd/cp)*(sfctmp**4.)/(sfcprs*ch) + 1.0
!
   delta = (slv/cp)*dqsdt2
   pc = (rr+delta)/(rr*(1.+rc*ch)+delta)
!
   return
   end subroutine noah_canopy_res
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_bare_soil_solver(etp,eta,prcp,smc,smcmax,smcwlt,            &
                    smcref,smcdry,cmc,cmcmax,nsoil,dt,shdfac,                  &
                    sbeta,q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,                 &
                    stc,epsca,bexp,pc,rch,rr,cfactr,                           &
                    sh2o,slope,kdt,frzfact,psisat,zsoil,                       &
                    dksat,dwsat,tbot,zbot,runoff1,runoff2,                     &
                    runoff3,edir,ec,et,ett,nroot,ice,rtdis,                    &
                    quartz,fxexp,csoil,                                        &
                    beta,drip,dew,flx1,flx2,flx3)
!-------------------------------------------------------------------------------
!
! calculate soil moisture and heat flux values and update soil moisture
! content and soil heat content values for the case when no snow pack is
! present.
!
!-------------------------------------------------------------------------------
   use constant, only : cp => cp_
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
!
   integer              ::  ice
   integer              ::  nroot
   integer              ::  nsoil
!
   real                 ::  bexp
   real                 ::  beta
   real                 ::  cfactr
   real                 ::  cmc
   real                 ::  cmcmax
   real                 ::  csoil
   real                 ::  dew
   real                 ::  df1
   real                 ::  dksat
   real                 ::  drip
   real                 ::  dt
   real                 ::  dwsat
   real                 ::  ec
   real                 ::  edir
   real                 ::  epsca
   real                 ::  eta
   real                 ::  eta1
   real                 ::  etp
   real                 ::  etp1
   real                 ::  et(nsoil)
   real                 ::  ett
   real                 ::  fdown
   real                 ::  f1
   real                 ::  fxexp
   real                 ::  flx1
   real                 ::  flx2
   real                 ::  flx3
   real                 ::  frzfact
   real                 ::  kdt
   real                 ::  pc
   real                 ::  prcp
   real                 ::  prcp1
   real                 ::  psisat
   real                 ::  q2
   real                 ::  quartz
   real                 ::  rch
   real                 ::  rr
   real                 ::  rtdis(nsoil)
   real                 ::  runoff1
   real                 ::  runoff2
   real                 ::  runoff3
   real                 ::  ssoil
   real                 ::  sbeta
   real                 ::  sfctmp
   real                 ::  shdfac
   real                 ::  sh2o(nsoil)
   real                 ::  slope
   real                 ::  smc(nsoil)
   real                 ::  smcdry
   real                 ::  smcmax
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  stc(nsoil)
   real                 ::  t1
   real                 ::  t24
   real                 ::  tbot
   real                 ::  th2
   real                 ::  yy
   real                 ::  yynum
   real                 ::  zbot
   real                 ::  zsoil(nsoil)
   real                 ::  zz1
!
   real,parameter       ::  sigma = 5.67e-8
!-------------------------------------------------------------------------------
!
! executable code begins here:
! convert etp from kg m-2 s-1 to ms-1 and initialize dew.
!
   prcp1 = prcp * 0.001
   etp1 = etp * 0.001
   dew = 0.0
   if (etp .gt. 0.0) then
!
! convert prcp from 'kg m-2 s-1' to 'm s-1'.
!
     call noah_soil_moisture_flux (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,     &
                    sh2o,slope,kdt,frzfact,                                    &
                    smcmax,bexp,pc,smcwlt,dksat,dwsat,                         &
                    smcref,shdfac,cmcmax,                                      &
                    smcdry,cfactr,runoff1,runoff2,runoff3,                     &
                    edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,                &
                    drip)
!
! convert modeled evapotranspiration fm  m s-1  to  kg m-2 s-1
!
     eta = eta1 * 1000.0
   else
!
! if etp < 0, assume dew forms (transform etp1 into dew and reinitialize
! etp1 to zero).
!
     dew = -etp1
     etp1 = 0.0
!
! convert prcp from 'kg m-2 s-1' to 'm s-1' and add dew amount.
!
     prcp1 = prcp1 + dew
!
     call noah_soil_moisture_flux (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,     &
               sh2o,slope,kdt,frzfact,                                         &
               smcmax,bexp,pc,smcwlt,dksat,dwsat,                              &
               smcref,shdfac,cmcmax,                                           &
               smcdry,cfactr,runoff1,runoff2,runoff3,                          &
               edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,                     &
               drip)
!
! convert modeled evapotranspiration from 'm s-1' to 'kg m-2 s-1'.
!
     eta = eta1 * 1000.0
   endif
!
! based on etp and e values, determine beta
!
   if ( etp .le. 0.0 ) then
     beta = 0.0
     if ( etp .lt. 0.0 ) then
       beta = 1.0
       eta = etp
     endif
   else
     beta = eta / etp
   endif
!
! get soil thermal diffuxivity/conductivity for top soil lyr,
! calc. adjusted top lyr soil temp and adjusted soil flux, then
! call noah_soil_heat_flux to compute/update soil heat flux and soil temps.
!
   call noah_soilt_conductivity (df1,smc(1),quartz,smcmax,sh2o(1))
!
! vegetation greenness fraction reduction in subsurface heat flux 
! via reduction factor, which is convenient to apply here to 
! thermal diffusivity that is later used in noah_soilt_land 
! to compute sub sfc heat flux
! (see additional comments on veg effect sub-sfc heat flx in 
! routine phys_lsm_noah)
!
   df1 = df1 * exp(sbeta*shdfac)
!
! compute intermediate terms passed to routine noah_soilt_land (via routine 
! noah_soil_heat_flux below) for use in computing subsurface heat flux 
! in noah_soilt_land
!
   yynum = fdown - sigma * t24
   yy = sfctmp + (yynum/rch+th2-sfctmp-beta*epsca) / rr
   zz1 = df1 / ( -0.5 * zsoil(1) * rch * rr ) + 1.0
   call noah_soil_heat_flux (ssoil,stc,smc,smcmax,nsoil,t1,dt,yy,zz1,zsoil,    &
               tbot,zbot,smcwlt,psisat,sh2o,bexp,f1,df1,ice,                   &
               quartz,csoil)
!
! set flx1 and flx3 (noah_snow_cover_solverk phase change heat fluxes) 
! to zero since they are not used here in noah_snow_cover_solver.  
! flx2 (freezing rain heat flux) was similarly initialized in 
! the noah_penman routine.
!
   flx1 = 0.0
   flx3 = 0.0
!
   return
   end subroutine noah_bare_soil_solver
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_snow_cover_solver (etp,eta,prcp,prcp1,                      &
                      snowng,smc,smcmax,smcwlt,                                &
                      smcref,smcdry,cmc,cmcmax,nsoil,dt,                       &
                      sbeta,df1,                                               &
                      q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,stc,epsca,           &
                      sfcprs,bexp,pc,rch,rr,cfactr,sncovr,esd,sndens,          &
                      snowh,sh2o,slope,kdt,frzfact,psisat,snup,                &
                      zsoil,dwsat,dksat,tbot,zbot,shdfac,runoff1,              &
                      runoff2,runoff3,edir,ec,et,ett,nroot,snomlt,             &
                      ice,rtdis,quartz,fxexp,csoil,                            &
                      beta,drip,dew,flx1,flx2,flx3)
!-------------------------------------------------------------------------------
!
! abstaract: calculate soil moisture and heat flux values & update soil moisture
!            content and soil heat content values for the case when a snow pack 
!            is present.
!
! program history log:
!   2000-01-01  fei chen               initial implementation 
!   2008-08-01  kyeong-hee seol        debugged, scm option
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real,parameter       ::  cp = 1004.5
   real,parameter       ::  cph2o = 4.218e+3
   real,parameter       ::  cpice = 2.106e+3
   real,parameter       ::  esdmin = 1.e-6
   real,parameter       ::  lsubf = 3.335e+5
   real,parameter       ::  lsubc = 2.501000e+6
   real,parameter       ::  lsubs = 2.83e+6
   real,parameter       ::  sigma = 5.67e-8
   real,parameter       ::  tfreez = 273.15
!
   integer              ::  ice
   integer              ::  nroot
   integer              ::  nsoil
!
   logical              ::  snowng
!
   real                 ::  bexp
   real                 ::  beta
   real                 ::  cfactr
   real                 ::  cmc
   real                 ::  cmcmax
   real                 ::  csoil
   real                 ::  denom
   real                 ::  dew
   real                 ::  df1
   real                 ::  dksat
   real                 ::  drip
   real                 ::  dsoil
   real                 ::  dtot
   real                 ::  dt
   real                 ::  dwsat
   real                 ::  ec
   real                 ::  edir
   real                 ::  epsca
   real                 ::  esd
   real                 ::  expsno
   real                 ::  expsoi
   real                 ::  eta
   real                 ::  eta1
   real                 ::  etp
   real                 ::  etp1
   real                 ::  etp2
   real                 ::  et(nsoil)
   real                 ::  ett
   real                 ::  ex
   real                 ::  expfac
   real                 ::  fdown
   real                 ::  fxexp
   real                 ::  flx1
   real                 ::  flx2
   real                 ::  flx3
   real                 ::  f1
   real                 ::  kdt
   real                 ::  pc
   real                 ::  prcp
   real                 ::  prcp1
   real                 ::  q2
   real                 ::  rch
   real                 ::  rr
   real                 ::  rtdis(nsoil)
   real                 ::  ssoil
   real                 ::  sbeta
   real                 ::  ssoil1
   real                 ::  sfctmp
   real                 ::  shdfac
   real                 ::  smc(nsoil)
   real                 ::  sh2o(nsoil)
   real                 ::  smcdry
   real                 ::  smcmax
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  snomlt
   real                 ::  snowh
   real                 ::  stc(nsoil)
   real                 ::  t1
   real                 ::  t11
   real                 ::  t12
   real                 ::  t12a
   real                 ::  t12b
   real                 ::  t24
   real                 ::  tbot
   real                 ::  zbot
   real                 ::  th2
   real                 ::  yy
   real                 ::  zsoil(nsoil)
   real                 ::  zz1
   real                 ::  salp
   real                 ::  sfcprs
   real                 ::  slope
   real                 ::  frzfact
   real                 ::  psisat
   real                 ::  snup
   real                 ::  runoff1
   real                 ::  runoff2
   real                 ::  runoff3
   real                 ::  quartz
   real                 ::  sndens
   real                 ::  sncond
   real                 ::  rsnow
   real                 ::  sncovr
   real                 ::  qsat
   real                 ::  etp3
   real                 ::  seh
   real                 ::  t14
   real                 ::  ccsnow
!-------------------------------------------------------------------------------
! executable code begins here:
! convert potential evap (etp) from kg m-2 s-1 to m s-1 and then to an
! amount (m) given timestep (dt) and call it an effective snowpack
! reduction amount, etp2 (m).  this is the amount the snowpack would be
! reduced due to evaporation from the snow sfc during the timestep.
! evaporation will proceed at the potential rate unless the snow depth
! is less than the expected snowpack reduction.
! if seaice (ice=1), beta remains=1.
! ----------------------------------------------------------------------
   prcp1 = prcp1*0.001
!
   etp2 = etp * 0.001 * dt
   beta = 1.0
   if (ice .ne. 1) then
     if (esd .lt. etp2) then
       beta = esd / etp2
     endif
   endif
! ----------------------------------------------------------------------
! if etp<0 (downward) then dewfall (=frostfall in this case).
! ----------------------------------------------------------------------
   dew = 0.0
   if (etp .lt. 0.0) then
     dew = -etp * 0.001
   endif
! ----------------------------------------------------------------------
! if precip is falling, calculate heat flux from snow sfc to newly
! accumulating precip.  note that this reflects the flux appropriate for
! the not-yet-updated skin temperature (t1).  assumes temperature of the
! snowfall striking the gound is =sfctmp (lowest model level air temp).
! ----------------------------------------------------------------------
   flx1 = 0.0
   if (snowng) then
     flx1 = cpice * prcp * (t1 - sfctmp)
   else
     if (prcp .gt. 0.0) flx1 = cph2o * prcp * (t1 - sfctmp)
   endif
! ----------------------------------------------------------------------
! calculate an 'effective snow-grnd sfc temp' (t12) based on heat fluxes
! between the snow pack and the soil and on net radiation.
! include flx1 (precip-snow sfc) and flx2 (freezing rain latent heat)
! fluxes.  flx1 from above, flx2 brought in via commom block rite.
! flx2 reflects freezing rain latent heat flux using t1 calculated in
! noah_penman.
! ----------------------------------------------------------------------
   dsoil = -(0.5 * zsoil(1))
   dtot = snowh + dsoil
   denom = 1.0 + df1 / (dtot * rr * rch)
   t12a = ( (fdown-flx1-flx2-sigma*t24)/rch+ th2 - sfctmp - beta*epsca ) / rr
   t12b = df1 * stc(1) / (dtot * rr * rch)
   t12 = (sfctmp + t12a + t12b) / denom      
! ----------------------------------------------------------------------
! if the 'effective snow-grnd sfc temp' is at or below freezing, no snow
! melt will occur.  set the skin temp to this effective temp.  reduce
! (by sublimation ) or increase (by frost) the depth of the snowpack,
! depending on sign of etp.
! update soil heat flux (ssoil) using new skin temperature (t1)
! since no snowmelt, set accumulated snowmelt to zero, set 'effective'
! precip from snowmelt to zero, set phase-change heat flux from snowmelt
! to zero.
! ----------------------------------------------------------------------
   if (t12 .le. tfreez) then
     t1 = t12
     ssoil = df1 * (t1 - stc(1)) / dtot
     esd = max(0.0, esd-etp2)
     flx3 = 0.0
     ex = 0.0
     snomlt = 0.0
   else
! ----------------------------------------------------------------------
! if the 'effective snow-grnd sfc temp' is above freezing, snow melt
! will occur.  call the snow melt rate,ex and amt, snomlt.  revise the
! effective snow depth.  revise the skin temp because it would have chgd
! due to the latent heat released by the melting. calc the latent heat
! released, flx3. set the effective precip, prcp1 to the snow melt rate,
! ex for use in noah_soil_moisture_flux.  
! adjustment to t1 to account for snow patches.
! calculate qsat valid at freezing point.  note that esat (saturation
! vapor pressure) value of 6.11e+2 used here is that valid at frzzing
! point.  note that etp from call noah_penman in phys_lsm_noah is ignored here in
! favor of bulk etp over 'open water' at freezing temp.
! update soil heat flux (s) using new skin temperature (t1)
! ----------------------------------------------------------------------
     t1 = tfreez * sncovr + t12 * (1.0 - sncovr)
     qsat = (0.622*6.11e2)/(sfcprs-0.378*6.11e2)
     etp = rch*(qsat-q2)/cp
     etp2 = etp*0.001*dt
     beta = 1.0
     ssoil = df1 * (t1 - stc(1)) / dtot
! ----------------------------------------------------------------------
! if potential evap (sublimation) greater than depth of snowpack.
! beta<1
! snowpack has sublimated away, set depth to zero.
! ----------------------------------------------------------------------
     if (esd .le. etp2) then
       beta = esd / etp2
       esd = 0.0
       ex = 0.0
       snomlt = 0.0
     else
! ----------------------------------------------------------------------
! potential evap (sublimation) less than depth of snowpack, retain
!   beta=1.
! snowpack (esd) reduced by potential evap rate
! etp3 (convert to flux)
! ----------------------------------------------------------------------
       esd = esd-etp2
       etp3 = etp*lsubc
       seh = rch*(t1-th2)
       t14 = t1*t1
       t14 = t14*t14
       flx3 = fdown - flx1 - flx2 - sigma*t14 - ssoil - seh - etp3
       if (flx3 .le.0.0) flx3 = 0.0
       ex = flx3*0.001/lsubf
! ----------------------------------------------------------------------
! snowmelt reduction depending on snow cover
! if snow cover less than 5% no snowmelt reduction
! ***note:  does 'if' below fail to match the melt water with the melt
!           energy?
! ----------------------------------------------------------------------
       if (sncovr .gt. 0.05) ex = ex * sncovr
       snomlt = ex * dt
! ----------------------------------------------------------------------
! esdmin represents a snowpack depth threshold value below which we
! choose not to retain any snowpack, and instead include it in snowmelt.
! ----------------------------------------------------------------------
       if (esd-snomlt .ge. esdmin) then
         esd = esd - snomlt
       else
! ----------------------------------------------------------------------
! snowmelt exceeds snow depth
! ----------------------------------------------------------------------
         ex = esd/dt
         flx3 = ex*1000.0*lsubf
         snomlt = esd
         esd = 0.0
       endif
! ----------------------------------------------------------------------
! end of 'esd .le. etp2' if-block
! ----------------------------------------------------------------------
     endif
     prcp1 = prcp1 + ex
! ----------------------------------------------------------------------
! end of 't12 .le. tfreez' if-block
! ----------------------------------------------------------------------
   endif
! ----------------------------------------------------------------------
! final beta now in hand, so compute evaporation.  evap equals etp
! unless beta<1.
! ----------------------------------------------------------------------
   eta = beta*etp
! ----------------------------------------------------------------------
! set the effective potnl evapotransp (etp1) to zero since this is snow
! case, so surface evap not calculated from edir, ec, or ett 
! in noah_soil_moisture_flux (below).
! if seaice (ice=1) skip call to noah_soil_moisture_flux.
! noah_soil_moisture_flux returns updated soil moisture values.  
! in this, the snow pack case, eta1 is not used in calculation of evap.
! ----------------------------------------------------------------------
   etp1 = 0.0
   if (ice .ne. 1) then
     call noah_soil_moisture_flux (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,     &
                 sh2o,slope,kdt,frzfact,                                       &
                 smcmax,bexp,pc,smcwlt,dksat,dwsat,                            &
                 smcref,shdfac,cmcmax,                                         &
                 smcdry,cfactr,runoff1,runoff2,runoff3,                        &
                 edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,                   &
                 drip)
   endif
! ----------------------------------------------------------------------
! before call noah_soil_heat_flux in this snowpack case, set zz1 and yy 
! arguments to special values that ensure that ground heat flux calculated 
! in noah_soil_heat_flux matches that already computer for below the snowpack, 
! thus the sfc heat flux to be computed in noah_soil_heat_flux will 
! effectively be the flux at the snow top surface.  
! t11 is a dummy arguement so we will not use the skin temp value 
! as revised by noah_soil_heat_flux.
! ----------------------------------------------------------------------
   zz1 = 1.0
   yy = stc(1)-0.5*ssoil*zsoil(1)*zz1/df1
   t11 = t1
! ----------------------------------------------------------------------
! noah_soil_heat_flux will calc/update the soil temps.  
! note:  the sub-sfc heat flux (ssoil1) and the skin temp (t11) output from
! this noah_soil_heat_flux call are not used  in any subsequent calculations. 
! rather, they are dummy variables here in the noah_snow_cover_solver case, 
! since the skin temp and sub-sfc heat flux are updated instead near 
! the beginning of the call to noah_snow_cover_solver.
! ----------------------------------------------------------------------
   call noah_soil_heat_flux (ssoil1,stc,smc,smcmax,nsoil,t11,dt,yy,zz1,zsoil,  &
               tbot,zbot,smcwlt,psisat,sh2o,bexp,f1,df1,ice,                   &
               quartz,csoil)
! ----------------------------------------------------------------------
! snow depth and density adjustment based on snow compaction.  yy is
! assumed to be the soil temperture at the top of the soil column.
! ----------------------------------------------------------------------
   if (esd .gt. 0.) then
     call noah_snowpack (esd,dt,snowh,sndens,t1,yy)
   else
     esd = 0.
     snowh = 0.
     sndens = 0.
     sncond = 1.
   endif
!
   return
   end subroutine noah_snow_cover_solver
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soil_moisture_flux (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil, &
                     sh2o,slope,kdt,frzfact,                                   &
                     smcmax,bexp,pc,smcwlt,dksat,dwsat,                        &
                     smcref,shdfac,cmcmax,                                     &
                     smcdry,cfactr,runoff1,runoff2,runoff3,                    &
                     edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,               &
                     drip)
!-------------------------------------------------------------------------------
!
! calculate soil moisture flux.  the soil moisture content (smc - a per
! unit volume measurement) is a dependent variable that is updated with
! prognostic eqns. the canopy moisture content (cmc) is also updated.
! frozen ground version:  new states added: sh2o, and frozen ground
! correction factor, frzfact and parameter slope.
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer,parameter    ::  nsold = 20
!
   integer              ::  i
   integer              ::  k
   integer              ::  nsoil
   integer              ::  nroot
!
   real                 ::  ai(nsold)
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
!
   real                 ::  bexp
   real                 ::  beta
   real                 ::  cfactr
   real                 ::  cmc
   real                 ::  cmc2ms
   real                 ::  cmcmax
   real                 ::  ddevap
   real                 ::  dksat
   real                 ::  drip
   real                 ::  dt
   real                 ::  dummy
   real                 ::  dwsat
   real                 ::  ec
   real                 ::  edir
   real                 ::  et(nsoil)
   real                 ::  eta1
   real                 ::  etp1
   real                 ::  ett
   real                 ::  excess
   real                 ::  frzfact
   real                 ::  fxexp
   real                 ::  kdt
   real                 ::  pc
   real                 ::  pcpdrp
   real                 ::  prcp1
   real                 ::  q2
   real                 ::  rhsct
   real                 ::  rhstt(nsold)
   real                 ::  rtdis(nsoil)
   real                 ::  runoff1
   real                 ::  runoff2
   real                 ::  runoff3
   real                 ::  sfctmp
   real                 ::  shdfac
   real                 ::  smc(nsoil)
   real                 ::  sh2o(nsoil)
   real                 ::  sice(nsold)
   real                 ::  sh2oa(nsold)
   real                 ::  sh2ofg(nsold)
   real                 ::  slope
   real                 ::  smcdry
   real                 ::  smcmax
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  trhsct
   real                 ::  zsoil(nsoil)
!-------------------------------------------------------------------------------
!
! temperature criteria for snowfall tfreez should have same value as in
! phys_lsm_noah.f (main subroutine)
!
   real,parameter       ::  tfreez = 273.15
!
! executable code begins here if the potential evapotranspiration is
! greater than zero.
!
   dummy = 0.
   edir = 0.
   ec = 0.
   ett = 0.
   do k = 1,nsoil
      et(k) = 0.
   end do
!
   if (etp1 .gt. 0.0) then
!
! retrieve direct evaporation from soil surface.  call this function
! only if veg cover not complete.
! frozen ground version:  sh2o states replace smc states.
!
     if (shdfac .lt. 1.) then
       edir = ddevap(etp1,smc(1),zsoil(1),shdfac,smcmax,                       &
                    bexp,dksat,dwsat,smcdry,smcref,smcwlt,fxexp)
     endif
!
! initialize plant total transpiration, retrieve plant transpiration,
! and accumulate it for all soil layers.
!
     if (shdfac.gt.0.0) then
       call noah_soilm_transpiration                                           &
                   (et,nsoil,etp1,sh2o,cmc,zsoil,shdfac,smcwlt,                &
                    cmcmax,pc,cfactr,smcref,sfctmp,q2,nroot,rtdis)
       do k = 1,nsoil
         ett = ett + et ( k )
       end do
!
! calculate canopy evaporation.
! if statements to avoid tangent linear problems near cmc=0.0.
!
       if (cmc .gt. 0.0) then
            ec = shdfac * ( ( cmc / cmcmax ) ** cfactr ) * etp1
       else
            ec = 0.0
       endif
!
! ec should be limited by the total amount of available water on the
! canopy.  -f.chen, 18-oct-1994
!
       cmc2ms = cmc / dt
       ec = min ( cmc2ms, ec )
     endif
   endif
!
! total up evap and noah_soilm_transpiration types
! to obtain actual evapotransp
!
   eta1 = edir + ett + ec
!
! compute the right hand side of the canopy eqn term ( rhsct )
!
   rhsct = shdfac * prcp1 - ec
!
! convert rhsct (a rate) to trhsct (an amount) and add it to existing
! cmc.  if resulting amt exceeds max capacity, it becomes drip and will
! fall to the grnd.
!
   drip = 0.
   trhsct = dt * rhsct
   excess = cmc + trhsct
   if (excess .gt. cmcmax) drip = excess - cmcmax
!
! pcpdrp is the combined prcp1 and drip (from cmc) that goes into the
! soil
!
   pcpdrp = (1. - shdfac) * prcp1 + drip / dt
!
! store ice content at each soil layer 
! before calling noah_soilm_budget & noah_soilm_matrix
!
   do i = 1,nsoil
     sice(i) = smc(i) - sh2o(i)
   end do
!
! call subroutines noah_soilm_budget and noah_soilm_matrix 
!   to solve the soil moisture tendency equations. 
!
! if the infiltrating precip rate is nontrivial,
!   (we consider nontrivial to be a precip total over the time step 
!   exceeding one one-thousandth of the water holding capacity of 
!   the first soil layer)
! then call the noah_soilm_budget/noah_soilm_matrix subroutine pair twice 
!   in the manner of time scheme "f" (implicit state, averaged coefficient)
!   of section 2 of kalnay and kanamitsu (1988, mwr, vol 116, 
!   pages 1945-1958)to minimize 2-delta-t oscillations in the 
!   soil moisture value of the top soil layer that can arise because
!   of the extreme nonlinear dependence of the soil hydraulic 
!   diffusivity coefficient and the hydraulic conductivity on the
!   soil moisture state
! otherwise call the noah_soilm_budget/noah_soilm_matrix subroutine pair 
!   once in the manner of time scheme "d" 
!   (implicit state, explicit coefficient) 
!   of section 2 of kalnay and kanamitsu
! pcpdrp is units of kg/m**2/s or mm/s, zsoil is negative depth in m 
!
! ----------------------------------------------------------------------
   if ( (pcpdrp*dt) .gt. (0.001*1000.0*(-zsoil(1))*smcmax) ) then
!
! frozen ground version:
! smc states replaced by sh2o states in noah_soilm_budget subr. 
! sh2o & sice states included in noah_soilm_matrix subr.  
! frozen ground correction factor, frzfact
! added.  all water balance calculations using unfrozen water
!
     call noah_soilm_budget (rhstt,edir,et,sh2o,sh2o,nsoil,pcpdrp,zsoil,       &
               dwsat,dksat,smcmax,bexp,runoff1,                                &
               runoff2,dt,smcwlt,slope,kdt,frzfact,sice,ai,bi,ci)
     call noah_soilm_matrix (sh2ofg,sh2o,dummy,rhstt,rhsct,dt,nsoil,smcmax,    &
                 cmcmax,runoff3,zsoil,smc,sice,ai,bi,ci)
     do k = 1,nsoil
       sh2oa(k) = (sh2o(k) + sh2ofg(k)) * 0.5
     end do
     call noah_soilm_budget (rhstt,edir,et,sh2o,sh2oa,nsoil,pcpdrp,zsoil,      &
               dwsat,dksat,smcmax,bexp,runoff1,                                &
               runoff2,dt,smcwlt,slope,kdt,frzfact,sice,ai,bi,ci)
     call noah_soilm_matrix (sh2o,sh2o,cmc,rhstt,rhsct,dt,nsoil,smcmax,        &
                 cmcmax,runoff3,zsoil,smc,sice,ai,bi,ci)
   else
     call noah_soilm_budget (rhstt,edir,et,sh2o,sh2o,nsoil,pcpdrp,zsoil,       &
               dwsat,dksat,smcmax,bexp,runoff1,                                &
               runoff2,dt,smcwlt,slope,kdt,frzfact,sice,ai,bi,ci)
     call noah_soilm_matrix (sh2o,sh2o,cmc,rhstt,rhsct,dt,nsoil,smcmax,        &
                 cmcmax,runoff3,zsoil,smc,sice,ai,bi,ci)
   endif
!      runof = runoff
!
   return
   end subroutine noah_soil_moisture_flux
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilm_transpiration                                         &
                     (et,nsoil,etp1,smc,cmc,zsoil,shdfac,smcwlt,               &
                      cmcmax,pc,cfactr,smcref,sfctmp,q2,nroot,rtdis)
!-------------------------------------------------------------------------------
!
! calculate transpiration for the veg class.
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer              ::  i
   integer              ::  k
   integer              ::  nsoil
   integer              ::  nroot
!
   real                 ::  cfactr
   real                 ::  cmc
   real                 ::  cmcmax
   real                 ::  denom
   real                 ::  et(nsoil)
   real                 ::  etp1
   real                 ::  etp1a
   real                 ::  gx (7)
!
!.....real part(nsoil)
!
   real                 ::  pc
   real                 ::  q2
   real                 ::  rtdis(nsoil)
   real                 ::  rtx
   real                 ::  sfctmp
   real                 ::  sgx
   real                 ::  shdfac
   real                 ::  smc(nsoil)
   real                 ::  smcref
   real                 ::  smcwlt
   real                 ::  zsoil(nsoil)
!-------------------------------------------------------------------------------
!
! initialize plant noah_soilm_transpiration to zero for all soil layers.
!
   do k = 1,nsoil
     et(k) = 0.
   end do
!
! calculate an 'adjusted' potential transpiration
! if statement below to avoid tangent linear problems near zero
! note: gx and other terms below redistribute transpiration by layer,
! et(k), as a function of soil moisture availability, while preserving
! total etp1a.
!
   if (cmc .gt. 0.0) then
     etp1a = shdfac * pc * etp1 * (1.0 - (cmc /cmcmax) ** cfactr)
   else
     etp1a = shdfac * pc * etp1
   endif
   sgx = 0.0
   do i = 1,nroot
     gx(i) = ( smc(i) - smcwlt ) / ( smcref - smcwlt )
     gx(i) = max ( min ( gx(i), 1. ), 0. )
     sgx = sgx + gx (i)
   end do
   sgx = sgx / nroot
   denom = 0.
   do i = 1,nroot
     rtx = rtdis(i) + gx(i) - sgx
     gx(i) = gx(i) * max ( rtx, 0. )
     denom = denom + gx(i)
   end do
   if (denom .le. 0.0) denom = 1.
   do i = 1,nroot
     et(i) = etp1a * gx(i) / denom
   end do
! 
! above code assumes a vertically uniform root distribution
! code below tests a variable root distribution
! 
!      et(1) = ( zsoil(1) / zsoil(nroot) ) * gx * etp1a
!      et(1) = ( zsoil(1) / zsoil(nroot) ) * etp1a
!
! using root distribution as weighting factor
! 
!      et(1) = rtdis(1) * etp1a
!      et(1) = etp1a * part(1)
! 
! loop down thru the soil layers repeating the operation above,
! but using the thickness of the soil layer (rather than the
! absolute depth of each layer) in the final calculation.
! 
!      do k = 2,nroot
!        gx = ( smc(k) - smcwlt ) / ( smcref - smcwlt )
!        gx = max ( min ( gx, 1. ), 0. )
! test canopy resistance
!        gx = 1.0
!        et(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot))*gx*etp1a
!        et(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot))*etp1a
! 
! using root distribution as weighting factor
! 
!        et(k) = rtdis(k) * etp1a
!        et(k) = etp1a*part(k)
!      end do      
!
   return
   end subroutine noah_soilm_transpiration
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilm_budget (rhstt,edir,et,sh2o,sh2oa,nsoil,pcpdrp,        &
                   zsoil,dwsat,dksat,smcmax,bexp,runoff1,                      &
                   runoff2,dt,smcwlt,slope,kdt,frzx,sice,ai,bi,ci)
!-------------------------------------------------------------------------------
!
! calculate the right hand side of the time tendency term of the soil
! water diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer,parameter    ::  nsold = 20
!
   integer              ::  ialp1
   integer              ::  iohinf
   integer              ::  j
   integer              ::  jj      
   integer              ::  k
   integer              ::  ks
   integer              ::  nsoil
!
   real                 ::  acrt
   real                 ::  ai(nsold)
   real                 ::  bexp
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
   real                 ::  dd
   real                 ::  ddt
   real                 ::  ddz
   real                 ::  ddz2
   real                 ::  denom
   real                 ::  denom2
   real                 ::  dice
   real                 ::  dksat
   real                 ::  dmax(nsold)
   real                 ::  dsmdz
   real                 ::  dsmdz2
   real                 ::  dt
   real                 ::  dt1
   real                 ::  dwsat
   real                 ::  edir
   real                 ::  et(nsoil)
   real                 ::  fcr
   real                 ::  frzx
   real                 ::  infmax
   real                 ::  kdt
   real                 ::  mxsmc
   real                 ::  mxsmc2
   real                 ::  numer
   real                 ::  pcpdrp
   real                 ::  pddum
   real                 ::  px
   real                 ::  rhstt(nsoil)
   real                 ::  runoff1
   real                 ::  runoff2
   real                 ::  sh2o(nsoil)
   real                 ::  sh2oa(nsoil)
   real                 ::  sice(nsoil)
   real                 ::  sicemax
   real                 ::  slope
   real                 ::  slopx
   real                 ::  smcav
   real                 ::  smcmax
   real                 ::  smcwlt
   real                 ::  sstt
   real                 ::  sum
   real                 ::  val
   real                 ::  wcnd
   real                 ::  wcnd2
   real                 ::  wdf
   real                 ::  wdf2
   real                 ::  zsoil(nsoil)
!
! frozen ground version:
! reference frozen ground parameter, cvfrz, is a shape parameter of
! areal distribution function of soil ice content which equals 1/cv.
! cv is a coefficient of spatial variation of soil ice content.  based
! on field data cv depends on areal mean of frozen depth, and it close
! to constant = 0.6 if areal mean frozen depth is above 20 cm.  that is
! why parameter cvfrz = 3 (int{1/0.6*0.6}).
! current logic doesn't allow cvfrz be bigger than 3
!
   integer,parameter     ::  cvfrz = 3
!-------------------------------------------------------------------------------
!
! determine rainfall infiltration rate and runoff.  include the
! infiltration formule from schaake and koren model.
! modified by q duan
!
   iohinf=1
!
! let sicemax be the greatest, if any, frozen water content within soil
! layers.
!
   sicemax = 0.0
   do ks = 1,nsoil
     if (sice(ks) .gt. sicemax) sicemax = sice(ks)
   end do
!
! determine rainfall infiltration rate and runoff
!
   pddum = pcpdrp
   runoff1 = 0.0
   if (pcpdrp .ne. 0.0) then
!
     dt1 = dt/86400.
     smcav = smcmax - smcwlt
     dmax(1)=-zsoil(1)*smcav
!
! frozen ground version:
!
     dice = -zsoil(1) * sice(1)
!          
     dmax(1)=dmax(1)*(1.0 - (sh2oa(1)+sice(1)-smcwlt)/smcav)
     dd=dmax(1)
     do ks = 2,nsoil
!
! frozen ground version:
!
       dice = dice + ( zsoil(ks-1) - zsoil(ks) ) * sice(ks)
       dmax(ks) = (zsoil(ks-1)-zsoil(ks))*smcav
       dmax(ks) = dmax(ks)*(1.0 - (sh2oa(ks)+sice(ks)-smcwlt)/smcav)
       dd = dd+dmax(ks)
     end do
!
! val = (1.-exp(-kdt*sqrt(dt1)))
! in below, remove the sqrt in above
!
     val = (1.-exp(-kdt*dt1))
     ddt = dd*val
     px = pcpdrp*dt  
     if (px .lt. 0.0) px = 0.0
     infmax = (px*(ddt/(px+ddt)))/dt
!
! frozen ground version:
! reduction of infiltration based on frozen ground parameters
!
     fcr = 1. 
     if (dice .gt. 1.e-2) then 
       acrt = cvfrz * frzx / dice 
       sum = 1.
       ialp1 = cvfrz - 1 
       do j = 1,ialp1
         k = 1
         do jj = j+1,ialp1
           k = k * jj
         end do
         sum = sum + (acrt ** ( cvfrz-j)) / float (k) 
       end do
       fcr = 1. - exp(-acrt) * sum 
     endif 
     infmax = infmax * fcr
!
! correction of infiltration limitation:
! if infmax .le. hydrolic conductivity assign infmax the value of
! hydrolic conductivity
!
!    mxsmc = max ( sh2oa(1), sh2oa(2) ) 
     mxsmc = sh2oa(1)
!
     call noah_soilm_coeff (wdf,wcnd,mxsmc,smcmax,bexp,dksat,dwsat, sicemax)
!
     infmax = max(infmax,wcnd)
     infmax = min(infmax,px)
!
     if (pcpdrp .gt. infmax) then
       runoff1 = pcpdrp - infmax
       pddum = infmax
     endif
   endif
!
! to avoid spurious drainage behavior, 'upstream differencing' in line
! below replaced with new approach in 2nd line:
! 'mxsmc = max(sh2oa(1), sh2oa(2))'
!
   mxsmc = sh2oa(1)
   call noah_soilm_coeff (wdf,wcnd,mxsmc,smcmax,bexp,dksat,dwsat,sicemax)
!
! calc the matrix coefficients ai, bi, and ci for the top layer
!
   ddz = 1. / ( -.5 * zsoil(2) )
   ai(1) = 0.0
   bi(1) = wdf * ddz / ( -zsoil(1) )
   ci(1) = -bi(1)
!
! calc rhstt for the top layer after calc'ng the vertical soil moisture
! gradient btwn the top and next to top layers.
!
   dsmdz = ( sh2o(1) - sh2o(2) ) / ( -.5 * zsoil(2) )
   rhstt(1) = (wdf * dsmdz + wcnd - pddum + edir + et(1))/zsoil(1)
   sstt = wdf * dsmdz + wcnd + edir + et(1)
!
! initialize ddz2
!
   ddz2 = 0.0
!
! loop thru the remaining soil layers, repeating the abv process
!
   do k = 2,nsoil
     denom2 = (zsoil(k-1) - zsoil(k))
     if (k .ne. nsoil) then
       slopx = 1.
!
! again, to avoid spurious drainage behavior, 'upstream differencing' in
! line below replaced with new approach in 2nd line:
! 'mxsmc2 = max (sh2oa(k), sh2oa(k+1))'
!
       mxsmc2 = sh2oa(k)
       call noah_soilm_coeff (wdf2,wcnd2,mxsmc2,                               &
                                smcmax,bexp,dksat,dwsat,sicemax)
!
! calc some partial products for later use in calc'ng rhstt
!
       denom = (zsoil(k-1) - zsoil(k+1))
       dsmdz2 = (sh2o(k) - sh2o(k+1)) / (denom * 0.5)
!
! calc the matrix coef, ci, after calc'ng its partial product
!
       ddz2 = 2.0 / denom
       ci(k) = -wdf2 * ddz2 / denom2
     else
!
! slope of bottom layer is introduced
!
       slopx = slope
!
! retrieve the soil water diffusivity and hydraulic conductivity for
! this layer
!
       call noah_soilm_coeff (wdf2,wcnd2,sh2oa(nsoil),                         &
                                smcmax,bexp,dksat,dwsat,sicemax)
!
! calc a partial product for later use in calc'ng rhstt
!
       dsmdz2 = 0.0
!
! set matrix coef ci to zero
!
       ci(k) = 0.0
     endif
!
! calc rhstt for this layer after calc'ng its numerator
!
     numer = (wdf2 * dsmdz2) + slopx * wcnd2 - (wdf * dsmdz) - wcnd + et(k)
     rhstt(k) = numer / (-denom2)
!
! calc matrix coefs, ai, and bi for this layer
!
     ai(k) = -wdf * ddz / denom2
     bi(k) = -( ai(k) + ci(k) )
!
! reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr
! runoff2:  sub-surface or baseflow runoff
!
     if (k .eq. nsoil) then
       runoff2 = slopx * wcnd2
     endif
     if (k .ne. nsoil) then
       wdf = wdf2
       wcnd = wcnd2
       dsmdz = dsmdz2
       ddz = ddz2
     endif
   end do
!
   return
   end subroutine noah_soilm_budget
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilm_coeff (wdf,wcnd,smc,                                  &
                                smcmax,bexp,dksat,dwsat,sicemax)
!-------------------------------------------------------------------------------
!
! calculate soil water diffusivity and soil hydraulic conductivity.
!
!-------------------------------------------------------------------------------
   implicit none
!
   real                 ::  bexp
   real                 ::  dksat
   real                 ::  dwsat
   real                 ::  expon
   real                 ::  factr1
   real                 ::  factr2
   real                 ::  sicemax
   real                 ::  smc
   real                 ::  smcmax
   real                 ::  vkwgt
   real                 ::  wcnd
   real                 ::  wdf
!-------------------------------------------------------------------------------
!
!     calc the ratio of the actual to the max psbl soil h2o content
!
   smc = smc
   smcmax = smcmax
   factr1 = 0.2 / smcmax
   factr2 = smc / smcmax
!
! prep an expntl coef and calc the soil water diffusivity
!
   expon = bexp + 2.0
   wdf = dwsat * factr2 ** expon
!
! frozen soil hydraulic diffusivity.  very sensitive to the vertical
! gradient of unfrozen water. the latter gradient can become very
! extreme in freezing/thawing situations, and given the relatively 
! few and thick soil layers, this gradient sufferes serious 
! trunction errors yielding erroneously high vertical transports of
! unfrozen water in both directions from huge hydraulic diffusivity.  
! therefore, we found we had to arbitrarily constrain wdf 
! --
! version d_10cm: ........  factr1 = 0.2/smcmax
! weighted approach...................... pablo grunmann, 28_sep_1999.
!
   if (sicemax .gt. 0.0)  then
     vkwgt = 1./(1.+(500.*sicemax)**3.)
     wdf = vkwgt*wdf + (1.- vkwgt)*dwsat*factr1**expon
   endif
!
! reset the expntl coef and calc the hydraulic conductivity
!
   expon = (2.0 * bexp) + 3.0
   wcnd = dksat * factr2 ** expon
!
   return
   end subroutine noah_soilm_coeff
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilm_matrix (sh2oout,sh2oin,cmc,rhstt,rhsct,dt,            &
                     nsoil,smcmax,cmcmax,runoff3,zsoil,smc,sice,               &
                     ai,bi,ci)
!-------------------------------------------------------------------------------
!
! calculate/update soil moisture content values and canopy moisture
! content values.
!
!-------------------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
   integer,parameter    ::  nsold = 20
!
   integer              ::  i
   integer              ::  k 
   integer              ::  kk11
   integer              ::  nsoil
!
   real                 ::  ai(nsold)
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
   real                 ::  ciin(nsold)
   real                 ::  cmc
   real                 ::  cmcmax
   real                 ::  ddz
   real                 ::  dt
   real                 ::  rhsct
   real                 ::  rhstt(nsoil)
   real                 ::  rhsttin(nsoil)
   real                 ::  runoff3
   real                 ::  sh2oin(nsoil)
   real                 ::  sh2oout(nsoil)
   real                 ::  sice(nsoil)
   real                 ::  smc(nsoil)
   real                 ::  smcmax
   real                 ::  stot
   real                 ::  wplus
   real                 ::  zsoil(nsoil)
! ----------------------------------------------------------------------
!
! create 'amount' values of variables to be input to the
! tri-diagonal matrix routine.
!
   do k = 1,nsoil
     rhstt(k) = rhstt(k) * dt
     ai(k) = ai(k) * dt
     bi(k) = 1. + bi(k) * dt
     ci(k) = ci(k) * dt
   end do
!
! copy values for input variables before call to noah_tridiag_matrix
!
   do k = 1,nsoil
     rhsttin(k) = rhstt(k)
   end do
!
   do k = 1,nsoil
     ciin(k) = ci(k)
   end do
!
   call noah_tridiag_matrix (ci,ai,bi,ciin,rhsttin,rhstt,nsoil)
!
! sum the previous smc value and the matrix solution to get a
! new value.  min allowable value of smc will be 0.02.
! runoff3: runoff within soil layers
!
   wplus = 0.0
   runoff3 = 0.
   ddz = -zsoil(1)
!
   do k = 1,nsoil
     if (k .ne. 1) ddz = zsoil(k - 1) - zsoil(k)
     sh2oout(k) = sh2oin(k) + ci(k) + wplus / ddz
!
     stot = sh2oout(k) + sice(k)
     if (stot .gt. smcmax) then
       if (k .eq. 1) then
         ddz = -zsoil(1)
       else
         kk11 = k - 1
         ddz = -zsoil(k) + zsoil(kk11)
       endif
       wplus = (stot-smcmax) * ddz
     else
       wplus = 0.
     endif
     smc(k) = max ( min(stot,smcmax),0.02 )
     sh2oout(k) = max((smc(k)-sice(k)),0.0)
   end do
!
   runoff3 = wplus
!
! update canopy water content/interception (cmc).  convert rhsct to 
! an 'amount' value and add to previous cmc value to get new cmc.
!
   cmc = cmc + dt * rhsct
   if (cmc .lt. 1.e-20) cmc=0.0
   cmc = min(cmc,cmcmax)
!
   return
   end subroutine noah_soilm_matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soil_heat_flux (ssoil,stc,smc,                              &
                     smcmax,nsoil,t1,dt,yy,zz1,zsoil,                          &
                     tbot,zbot,smcwlt,psisat,sh2o,bexp,f1,df1,ice,             &
                     quartz,csoil)
!-------------------------------------------------------------------------------
!      
! update the temperature state of the soil column based on the thermal
! diffusion equation and update the frozen soil moisture content based
! on the temperature.
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer,parameter    ::  nsold = 20
!
   integer              ::  i
   integer              ::  ice
   integer              ::  ifrz
   integer              ::  nsoil
!
   real                 ::  ai(nsold)
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
   real                 ::  bexp
   real                 ::  csoil
   real                 ::  df1
   real                 ::  dt
   real                 ::  f1
   real                 ::  psisat
   real                 ::  quartz
   real                 ::  rhsts(nsold)
   real                 ::  ssoil
   real                 ::  sh2o(nsoil)
   real                 ::  smc(nsoil)
   real                 ::  smcmax
   real                 ::  smcwlt
   real                 ::  stc(nsoil)
   real                 ::  stcf(nsold)
   real                 ::  t1
   real                 ::  tbot
   real                 ::  yy
   real                 ::  zbot
   real                 ::  zsoil(nsoil)
   real                 ::  zz1
!
   real,parameter       ::  t0 = 273.15
!-------------------------------------------------------------------------------
!
! noah_soilt_land routine calcs the right hand side of the soil temp dif eqn
!
   if (ice.eq.1) then
!
! sea-ice case
!
     call noah_soilt_ice (rhsts,stc,nsoil,zsoil,yy,zz1,df1,ai,bi,ci)
     call noah_soilt_matrix (stcf,stc,rhsts,dt,nsoil,ai,bi,ci)
   else
!
! land-mass case
!
     call noah_soilt_land (rhsts,stc,smc,smcmax,nsoil,zsoil,yy,zz1,tbot,       &
                zbot,psisat,sh2o,dt,                                           &
                bexp,f1,df1,quartz,csoil,ai,bi,ci)
     call noah_soilt_matrix (stcf,stc,rhsts,dt,nsoil,ai,bi,ci)
   endif
   do i = 1,nsoil
     stc(i) = stcf(i)
   end do
!
! in the no snowpack case (via routine noah_bare_soil_solver branch,) 
! update the grnd (skin) temperature here in response to 
! the updated soil temperature profile above.  
! (note: inspection of routine noah_snow_cover_solver 
! shows that t1 below is a dummy variable only, as skin temperature 
! is updated differently in routine noah_snow_cover_solver) 
!
   t1 = (yy + (zz1 - 1.0) * stc(1)) / zz1
!
! calculate surface soil heat flux
!
   ssoil = df1 * (stc(1) - t1) / (0.5 * zsoil(1))
!
   return
   end subroutine noah_soil_heat_flux
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilt_ice (rhsts,stc,nsoil,zsoil,yy,zz1,df1,ai,bi,ci)
!-------------------------------------------------------------------------------
!
! calculate the right hand side of the time tendency term of the soil
! thermal diffusion equation in the case of sea-ice pack.  also to
! compute (prepare) the matrix coefficients for the tri-diagonal matrix
! of the implicit time scheme.
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer,parameter    ::  nsold = 20
!
   integer              ::  k
   integer              ::  nsoil
!
   real                 ::  ai(nsold)
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
   real                 ::  ddz
   real                 ::  ddz2
   real                 ::  denom
   real                 ::  df1
   real                 ::  dtsdz
   real                 ::  dtsdz2
   real                 ::  rhsts(nsoil)
   real                 ::  ssoil
   real                 ::  stc(nsoil)
   real                 ::  tbot
   real                 ::  yy
   real                 ::  zbot
   real                 ::  zsoil(nsoil)
   real                 ::  zz1
!
   data tbot /271.16/
!-------------------------------------------------------------------------------
! 
! set a nominal universal value of the sea-ice specific heat capacity,
! hcpct = 1880.0*917.0.
!
   real,parameter       ::  hcpct = 1.72396e+6
! 
! the input argument df1 is a universally constant value of sea-ice
! thermal diffusivity, set in routine noah_snow_cover_solver as df1 = 2.2.
! 
! set ice pack depth.  use tbot as ice pack lower boundary temperature
! (that of unfrozen sea water at bottom of sea ice pack).  assume ice
! pack is of n=nsoil layers spanning a uniform constant ice pack
! thickness as defined by zsoil(nsoil) in routine phys_lsm_noah.
! 
   zbot = zsoil(nsoil)
! 
! calc the matrix coefficients ai, bi, and ci for the top layer
! 
   ddz = 1.0 / ( -0.5 * zsoil(2) )
   ai(1) = 0.0
   ci(1) = (df1 * ddz) / (zsoil(1) * hcpct)
   bi(1) = -ci(1) + df1/(0.5 * zsoil(1) * zsoil(1) * hcpct * zz1)
! 
! calc the vertical soil temp gradient btwn the top and 2nd soil layers.
! recalc/adjust the soil heat flux.  use the gradient and flux to calc
! rhsts for the top soil layer.
! 
   dtsdz = ( stc(1) - stc(2) ) / ( -0.5 * zsoil(2) )
   ssoil = df1 * ( stc(1) - yy ) / ( 0.5 * zsoil(1) * zz1 )
   rhsts(1) = ( df1 * dtsdz - ssoil ) / ( zsoil(1) * hcpct )
! 
! initialize ddz2
! 
   ddz2 = 0.0
! 
! loop thru the remaining soil layers, repeating the above process
! 
   do k = 2,nsoil
     if (k .ne. nsoil) then
! 
! calc the vertical soil temp gradient thru this layer.
! 
       denom = 0.5 * ( zsoil(k-1) - zsoil(k+1) )
       dtsdz2 = ( stc(k) - stc(k+1) ) / denom
! 
! calc the matrix coef, ci, after calc'ng its partial product.
! 
       ddz2 = 2. / (zsoil(k-1) - zsoil(k+1))
       ci(k) = -df1 * ddz2 / ((zsoil(k-1) - zsoil(k)) * hcpct)
     else
! 
! calc the vertical soil temp gradient thru the lowest layer.
! 
       dtsdz2 = (stc(k)-tbot)/(.5 * (zsoil(k-1) + zsoil(k))-zbot)
! 
! set matrix coef, ci to zero.
! 
       ci(k) = 0.
     endif
! 
! calc rhsts for this layer after calc'ng a partial product.
! 
     denom = ( zsoil(k) - zsoil(k-1) ) * hcpct
     rhsts(k) = ( df1 * dtsdz2 - df1 * dtsdz ) / denom
! 
! calc matrix coefs, ai, and bi for this layer.
! 
     ai(k) = - df1 * ddz / ((zsoil(k-1) - zsoil(k)) * hcpct)
     bi(k) = -(ai(k) + ci(k))
! 
! reset values of dtsdz and ddz for loop to next soil lyr.
! 
     dtsdz = dtsdz2
     ddz   = ddz2
   end do
! 
   return
   end subroutine noah_soilt_ice
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilt_matrix (stcout,stcin,rhsts,dt,nsoil,ai,bi,ci)
!-------------------------------------------------------------------------------
!
! calculate/update the soil temperature field.
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer,parameter    ::  nsold = 20
!
   integer              ::  k
   integer              ::  nsoil
!
   real                 ::  ai(nsold)
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
   real                 ::  ciin(nsold)
   real                 ::  dt
   real                 ::  rhsts(nsoil)
   real                 ::  rhstsin(nsoil)
   real                 ::  stcin(nsoil)
   real                 ::  stcout(nsoil)
!-------------------------------------------------------------------------------
! 
! create finite difference values for use in noah_tridiag_matrix routine
! 
   do k = 1,nsoil
     rhsts(k) = rhsts(k) * dt
     ai(k) = ai(k) * dt
     bi(k) = 1. + bi(k) * dt
     ci(k) = ci(k) * dt
   end do
! 
! copy values for input variables before call to noah_tridiag_matrix
! 
   do k = 1,nsoil
     rhstsin(k) = rhsts(k)
   end do
   do k = 1,nsoil
     ciin(k) = ci(k)
   end do
! 
! solve the tri-diagonal matrix equation
! 
   call noah_tridiag_matrix(ci,ai,bi,ciin,rhstsin,rhsts,nsoil)
! 
! calc/update the soil temps using matrix solution
! 
   do k = 1,nsoil
     stcout(k) = stcin(k) + ci(k)
   end do
! 
   return
   end subroutine noah_soilt_matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_tridiag_matrix (p,a,b,c,d,delta,nsoil)
!-------------------------------------------------------------------------------
!
! invert (solve) the tri-diagonal matrix problem shown below:
!
! ###                                            ### ###  ###   ###  ###
! #b(1), c(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! #a(2), b(2), c(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! # 0  , a(3), b(3), c(3),  0  ,   . . .  ,    0   # #      #   # d(3) #
! # 0  ,  0  , a(4), b(4), c(4),   . . .  ,    0   # # p(4) #   # d(4) #
! # 0  ,  0  ,  0  , a(5), b(5),   . . .  ,    0   # # p(5) #   # d(5) #
! # .                                          .   # #  .   # = #   .  #
! # .                                          .   # #  .   #   #   .  #
! # .                                          .   # #  .   #   #   .  #
! # 0  , . . . , 0 , a(m-2), b(m-2), c(m-2),   0   # #p(m-2)#   #d(m-2)#
! # 0  , . . . , 0 ,   0   , a(m-1), b(m-1), c(m-1)# #p(m-1)#   #d(m-1)#
! # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b(m) # # p(m) #   # d(m) #
! ###                                            ### ###  ###   ###  ###
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer              ::  k
   integer              ::  kk
   integer              ::  nsoil
!      
   real                 ::  a(nsoil)
   real                 ::  b(nsoil)
   real                 ::  c(nsoil)
   real                 ::  d(nsoil)
   real                 ::  delta(nsoil)
   real                 ::  p(nsoil)
!-------------------------------------------------------------------------------
! 
! initialize eqn coef c for the lowest soil layer
! 
   c(nsoil) = 0.0
! 
! solve the coefs for the 1st soil layer
! 
   p(1) = -c(1) / b(1)
   delta(1) = d(1) / b(1)
! 
! solve the coefs for soil layers 2 thru nsoil
! 
   do k = 2,nsoil
     p(k) = -c(k) * ( 1.0 / (b(k) + a (k) * p(k-1)) )
     delta(k) = (d(k)-a(k)*delta(k-1))*(1.0/(b(k)+a(k)*p(k-1)))
   end do
! 
! set p to delta for lowest soil layer
! 
   p(nsoil) = delta(nsoil)
! 
! adjust p for soil layers 2 thru nsoil
! 
   do k = 2,nsoil
     kk = nsoil - k + 1
     p(kk) = p(kk) * p(kk+1) + delta(kk)
   end do
! 
   return
   end subroutine noah_tridiag_matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilt_land (rhsts,stc,smc,smcmax,nsoil,zsoil,yy,zz1,        &
                   tbot,zbot,psisat,sh2o,dt,bexp,                              &
                   f1,df1,quartz,csoil,ai,bi,ci)
!-------------------------------------------------------------------------------
!
! calculate the right hand side of the time tendency term of the soil
! thermal diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer,parameter    ::  nsold = 20
!
   logical              ::  itavg
!
   integer              ::  i
   integer              ::  k
   integer              ::  nsoil
! 
! declare work arrays needed in tri-diagonal implicit solver
! 
   real                 ::  ai(nsold)
   real                 ::  bi(nsold)
   real                 ::  ci(nsold)
! 
! declarations
! 
   real                 ::  bexp
   real                 ::  csoil
   real                 ::  ddz
   real                 ::  ddz2
   real                 ::  denom
   real                 ::  df1
   real                 ::  df1n
   real                 ::  df1k
   real                 ::  dt
   real                 ::  dtsdz
   real                 ::  dtsdz2
   real                 ::  f1
   real                 ::  hcpct
   real                 ::  psisat
   real                 ::  quartz
   real                 ::  qtot
   real                 ::  rhsts(nsoil)
   real                 ::  ssoil
   real                 ::  sice
   real                 ::  smc(nsoil)
   real                 ::  sh2o(nsoil)
   real                 ::  smcmax
   real                 ::  ssnksrc
   real                 ::  stc(nsoil)
   real                 ::  tavg
   real                 ::  tbk
   real                 ::  tbk1
   real                 ::  tbot
   real                 ::  zbot
   real                 ::  tsnsr
   real                 ::  tsurf
   real                 ::  yy
   real                 ::  zsoil(nsoil)
   real                 ::  zz1
!
   real,parameter       ::  t0 = 273.15
!-------------------------------------------------------------------------------
!
! set specific heat capacities of air, water, ice, soil mineral       
! 
   real,parameter       ::  cair = 1004.0
   real,parameter       ::  ch2o = 4.2e6
   real,parameter       ::  cice = 2.106e6
!
! note: csoil now set in routine noah_read_parameter and passed in
!      parameter(csoil = 1.26e6)
! 
! initialize logical for soil layer temperature averaging.
! 
   itavg = .true.
!      itavg = .false.
! 
! begin section for top soil layer
! 
! calc the heat capacity of the top soil layer
! 
   hcpct = sh2o(1)*ch2o + (1.0-smcmax)*csoil + (smcmax-smc(1))*cair            &
           + ( smc(1) - sh2o(1) )*cice
! 
! calc the matrix coefficients ai, bi, and ci for the top layer
! 
   ddz = 1.0 / ( -0.5 * zsoil(2) )
   ai(1) = 0.0
   ci(1) = (df1 * ddz) / (zsoil(1) * hcpct)
   bi(1) = -ci(1) + df1 / (0.5 * zsoil(1) * zsoil(1)*hcpct*zz1)
! 
! calculate the vertical soil temp gradient btwn the 1st and 2nd soil
! layers.  then calculate the subsurface heat flux. use the temp
! gradient and subsfc heat flux to calc "right-hand side tendency
! terms", or "rhsts", for top soil layer.
! 
   dtsdz = (stc(1) - stc(2)) / (-0.5 * zsoil(2))
   ssoil = df1 * (stc(1) - yy) / (0.5 * zsoil(1) * zz1)
   rhsts(1) = (df1 * dtsdz - ssoil) / (zsoil(1) * hcpct)
! 
! next capture the vertical difference of the heat flux at top and
! bottom of first soil layer for use in heat flux constraint applied to
! potential soil freezing/thawing in routine ssnksrc.
! 
   qtot = ssoil - df1*dtsdz
! 
! if temperature averaging invoked (itavg=true; else skip):
! set temp "tsurf" at top of soil column (for use in freezing soil
! physics later in function subroutine ssnksrc).  if snowpack content is
! zero, then tsurf expression below gives tsurf = skin temp.  if
! snowpack is nonzero (hence argument zz1=1), then tsurf expression
! below yields soil column top temperature under snowpack.  then
! calculate temperature at bottom interface of 1st soil layer for use
! later in function subroutine ssnksrc
! 
   if (itavg) then 
     tsurf = (yy + (zz1-1) * stc(1)) / zz1
     call noah_soilt_bdy (stc(1),stc(2),zsoil,zbot,1,nsoil,tbk)
   endif
! 
! calculate frozen water content in 1st soil layer. 
! 
   sice = smc(1) - sh2o(1)
! 
! if frozen water present or any of layer-1 mid-point or bounding
! interface temperatures below freezing, then call ssnksrc to
! compute heat source/sink (and change in frozen water content)
! due to possible soil water phase change
! 
   if ( (sice   .gt. 0.) .or. (tsurf .lt. t0) .or.                             &
        (stc(1) .lt. t0) .or. (tbk   .lt. t0) ) then

     if (itavg) then 
       call noah_soilt_avg(tavg,tsurf,stc(1),tbk,zsoil,nsoil,1)
     else
       tavg = stc(1)
     endif
     tsnsr = ssnksrc (tavg,smc(1),sh2o(1),                                     &
         zsoil,nsoil,smcmax,psisat,bexp,dt,1,qtot)
!
     rhsts(1) = rhsts(1) - tsnsr / ( zsoil(1) * hcpct )
   endif
! 
! this ends section for top soil layer.
! 
! initialize ddz2
! 
   ddz2 = 0.0
! 
! loop thru the remaining soil layers, repeating the above process
! (except subsfc or "ground" heat flux not repeated in lower layers)
! 
   df1k = df1
   do k = 2,nsoil
! 
! calculate heat capacity for this soil layer.
! 
     hcpct = sh2o(k)*ch2o +(1.0-smcmax)*csoil +(smcmax-smc(k))*cair            &
           + ( smc(k) - sh2o(k) )*cice
     if (k .ne. nsoil) then
! 
! this section for layer 2 or greater, but not last layer.
! 
! calculate thermal diffusivity for this layer.
! 
       call noah_soilt_conductivity (df1n,smc(k),quartz,smcmax,sh2o(k))
! 
! calc the vertical soil temp gradient thru this layer
! 
       denom = 0.5 * ( zsoil(k-1) - zsoil(k+1) )
       dtsdz2 = ( stc(k) - stc(k+1) ) / denom
! 
! calc the matrix coef, ci, after calc'ng its partial product
! 
       ddz2 = 2. / (zsoil(k-1) - zsoil(k+1))
       ci(k) = -df1n * ddz2 / ((zsoil(k-1) - zsoil(k)) * hcpct)
! 
! if temperature averaging invoked (itavg=true; else skip):  calculate
! temp at bottom of layer.
! 
       if (itavg) then 
         call noah_soilt_bdy (stc(k),stc(k+1),zsoil,zbot,k,nsoil,tbk1)
       endif
     else
! 
! special case of bottom soil layer:  calculate thermal diffusivity for
! bottom layer.
! 
       call noah_soilt_conductivity (df1n,smc(k),quartz,smcmax,sh2o(k))
! 
! calc the vertical soil temp gradient thru bottom layer.
! 
       denom = .5 * (zsoil(k-1) + zsoil(k)) - zbot
       dtsdz2 = (stc(k)-tbot) / denom
! 
! set matrix coef, ci to zero if bottom layer.
! 
       ci(k) = 0.
! 
! if temperature averaging invoked (itavg=true; else skip):  calculate
! temp at bottom of last layer.
! 
       if (itavg) then 
         call noah_soilt_bdy (stc(k),tbot,zsoil,zbot,k,nsoil,tbk1)
       endif 
     endif
! 
! this ends special loop for bottom layer.
! 
! calculate rhsts for this layer after calc'ng a partial product.
! 
     denom = ( zsoil(k) - zsoil(k-1) ) * hcpct
     rhsts(k) = ( df1n * dtsdz2 - df1k * dtsdz ) / denom
     qtot = -1.0*denom*rhsts(k)
     sice = smc(k) - sh2o(k)
     if ( (sice .gt. 0.) .or. (tbk .lt. t0) .or.                               &
        (stc(k) .lt. t0) .or. (tbk1 .lt. t0) ) then
       if (itavg) then 
         call noah_soilt_avg(tavg,tbk,stc(k),tbk1,zsoil,nsoil,k)
       else
         tavg = stc(k)
       endif
       tsnsr = ssnksrc(tavg,smc(k),sh2o(k),zsoil,nsoil,                        &
                      smcmax,psisat,bexp,dt,k,qtot)
       rhsts(k) = rhsts(k) - tsnsr / denom
     endif 
! 
! calc matrix coefs, ai, and bi for this layer.
! 
     ai(k) = - df1 * ddz / ((zsoil(k-1) - zsoil(k)) * hcpct)
     bi(k) = -(ai(k) + ci(k))
! 
! reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.
! 
     tbk   = tbk1
     df1k  = df1n
     dtsdz = dtsdz2
     ddz   = ddz2
   end do
! 
   return
   end subroutine noah_soilt_land
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilt_avg (tavg,tup,tm,tdn,zsoil,nsoil,k) 
!-------------------------------------------------------------------------------
!      
! calculate soil layer average temperature (tavg) in freezing/thawing
! layer using up, down, and middle layer temperatures (tup, tdn, tm),
! where tup is at top boundary of layer, tdn is at bottom boundary of
! layer.  tm is layer prognostic state temperature.
!      
!-------------------------------------------------------------------------------
   implicit none
! 
   integer              ::  k
   integer              ::  nsoil
!
   real                 ::  dz
   real                 ::  dzh
   real                 ::  tavg
   real                 ::  tdn
   real                 ::  tm
   real                 ::  tup
   real                 ::  x0
   real                 ::  xdn
   real                 ::  xup
   real                 ::  zsoil (nsoil)
!
   real,parameter       ::  t0 = 2.7315e2
!-------------------------------------------------------------------------------
   if (k .eq. 1) then
     dz = -zsoil(1)
   else
     dz = zsoil(k-1)-zsoil(k)
   endif
   dzh=dz*0.5
   if (tup .lt. t0) then
     if (tm .lt. t0) then
       if (tdn .lt. t0) then
! 
! tup, tm, tdn < t0
! 
         tavg = (tup + 2.0*tm + tdn)/ 4.0            
       else
! 
! tup & tm < t0,  tdn >= t0
! 
         x0 = (t0 - tm) * dzh / (tdn - tm)
         tavg = 0.5 * (tup*dzh+tm*(dzh+x0)+t0*(2.*dzh-x0)) / dz
       endif      
     else
       if (tdn .lt. t0) then
! 
! tup < t0, tm >= t0, tdn < t0
! 
         xup  = (t0-tup) * dzh / (tm-tup)
         xdn  = dzh - (t0-tm) * dzh / (tdn-tm)
         tavg = 0.5 * (tup*xup+t0*(2.*dz-xup-xdn)+tdn*xdn) / dz
       else
! 
! tup < t0, tm >= t0, tdn >= t0
! 
         xup  = (t0-tup) * dzh / (tm-tup)
         tavg = 0.5 * (tup*xup+t0*(2.*dz-xup)) / dz
       endif   
     endif
   else
     if (tm .lt. t0) then
       if (tdn .lt. t0) then
! 
! tup >= t0, tm < t0, tdn < t0
! 
         xup  = dzh - (t0-tup) * dzh / (tm-tup)
         tavg = 0.5 * (t0*(dz-xup)+tm*(dzh+xup)+tdn*dzh) / dz
       else
! 
! tup >= t0, tm < t0, tdn >= t0
! 
         xup  = dzh - (t0-tup) * dzh / (tm-tup)
         xdn  = (t0-tm) * dzh / (tdn-tm)
         tavg = 0.5 * (t0*(2.*dz-xup-xdn)+tm*(xup+xdn)) / dz
       endif   
     else
       if (tdn .lt. t0) then
! 
! tup >= t0, tm >= t0, tdn < t0
! 
         xdn  = dzh - (t0-tm) * dzh / (tdn-tm)
         tavg = (t0*(dz-xdn)+0.5*(t0+tdn)*xdn) / dz                 
       else
! 
! tup >= t0, tm >= t0, tdn >= t0
! 
         tavg = (tup + 2.0*tm + tdn) / 4.0
       endif
     endif
   endif
!
   return
   end subroutine noah_soilt_avg
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_soilt_bdy (tu,tb,zsoil,zbot,k,nsoil,tbnd1)
!-------------------------------------------------------------------------------
!
! calculate temperature on the boundary of the layer by interpolation of
! the middle layer temperatures
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer              ::  nsoil
   integer              ::  k
!
   real                 ::  tbnd1
   real                 ::  tu
   real                 ::  tb
   real                 ::  zb
   real                 ::  zbot
   real                 ::  zup
   real                 ::  zsoil (nsoil)
!
   real,parameter       ::  t0 = 273.15
!-------------------------------------------------------------------------------
!
! use surface temperature on the top of the first layer
! 
   if (k .eq. 1) then
     zup = 0.
   else
     zup = zsoil(k-1)
   endif
! 
! use depth of the constant bottom temperature when interpolate
! temperature into the last layer boundary
! 
   if (k .eq. nsoil) then
     zb = 2.*zbot-zsoil(k)
   else
     zb = zsoil(k+1)
   endif
! 
! linear interpolation between the average layer temperatures
! 
   tbnd1 = tu+(tb-tu)*(zup-zsoil(k))/(zup-zb)
! 
   return
   end subroutine noah_soilt_bdy
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine noah_snowpack (esd,dtsec,snowh,sndens,tsnow,tsoil)
!-------------------------------------------------------------------------------
!
! calculate compaction of snowpack under conditions of increasing snow
! density, as obtained from an approximate solution of e. anderson's
! differential equation (3.29), noaa technical report nws 19, by victor
! koren, 03/25/95.
!
! esd     water equivalent of snow (m)
! dtsec   time step (sec)
! snowh   snow depth (m)
! sndens  snow density (g/cm3=dimensionless fraction of h2o density)
! tsnow   snow surface temperature (k)
! tsoil   soil surface temperature (k)
!
! subroutine will return new values of snowh and sndens
!
!-------------------------------------------------------------------------------
   implicit none
! 
   integer              ::  ipol, j
!
   real                 ::                                                     &
         bfac,sndens,dsx,dthr,dtsec,dw,snowhc,snowh,pexp,tavgc,                &
         tsnow,tsnowc,tsoil,tsoilc,esd,esdc,esdcx
!
   real,parameter       ::  c1 = 0.01, c2=21.0, g=9.81, kn=4000.0
! ----------------------------------------------------------------------
! conversion into simulation units
! ----------------------------------------------------------------------
   snowhc = snowh*100.
   esdc = esd*100.
   dthr = dtsec/3600.
   tsnowc = tsnow-273.15
   tsoilc = tsoil-273.15
! ----------------------------------------------------------------------
! calculating of average temperature of snow pack
! ----------------------------------------------------------------------
   tavgc = 0.5*(tsnowc+tsoilc)                                    
! ----------------------------------------------------------------------
! calculating of snow depth and density as a result of compaction
!  sndens=ds0*(exp(bfac*esd)-1.)/(bfac*esd)
!  bfac=dthr*c1*exp(0.08*tavgc-c2*ds0)
! note: bfac*esd in sndens eqn above has to be carefully treated
! numerically below:
!   c1 is the fractional increase in density (1/(cm*hr)) 
!   c2 is a constant (cm3/g) kojima estimated as 21 cms/g
! ----------------------------------------------------------------------
   if (esdc .gt. 1.e-2) then
     esdcx = esdc
   else
     esdcx = 1.e-2
   endif
   bfac = dthr*c1*exp(0.08*tavgc-c2*sndens)
!     dsx = sndens*((dexp(bfac*esdc)-1.)/(bfac*esdc))
! ----------------------------------------------------------------------
! the function of the form (e**x-1)/x imbedded in above expression
! for dsx was causing numerical difficulties when the denominator "x"
! (i.e. bfac*esdc) became zero or approached zero (despite the fact that
! the analytical function (e**x-1)/x has a well defined limit as 
! "x" approaches zero), hence below we replace the (e**x-1)/x 
! expression with an equivalent, numerically well-behaved 
! polynomial expansion.
!
! number of terms of polynomial expansion, and hence its accuracy, 
! is governed by iteration limit "ipol".
!      ipol greater than 9 only makes a difference on double
!            precision (relative errors given in percent %).
!       ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
!       ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
!       ipol=7, for rel.error <~ 1.8 e-4 % ...
! ----------------------------------------------------------------------
   ipol = 4
   pexp = 0.
   do j = ipol,1,-1
!    pexp = (1. + pexp)*bfac*esdc/real(j+1) 
     pexp = (1. + pexp)*bfac*esdcx/real(j+1) 
   end do
   pexp = pexp + 1.
   dsx = sndens*(pexp)
! ----------------------------------------------------------------------
! above line ends polynomial substitution
! ----------------------------------------------------------------------
!     end of korean formulation
!
!     base formulation (cogley et al., 1990)
!     convert density from g/cm3 to kg/m3
!       dsm=sndens*1000.0
! 
!       dsx=dsm+dtsec*0.5*dsm*g*esd/
!    &      (1e7*exp(-0.02*dsm+kn/(tavgc+273.16)-14.643))
! 
!     convert density from kg/m3 to g/cm3
!       dsx=dsx/1000.0
!
!     end of cogley et al. formulation 
!
! ----------------------------------------------------------------------
! set upper/lower limit on snow density
! ----------------------------------------------------------------------
   if (dsx .gt. 0.40) dsx = 0.40
   if (dsx .lt. 0.05) dsx = 0.05
   sndens = dsx
! ----------------------------------------------------------------------
! update of snow depth and density depending on liquid water during
! snowmelt.  assumed that 13% of liquid water can be stored in snow per
! day during snowmelt till snow density 0.40.
! ----------------------------------------------------------------------
   if (tsnowc .ge. 0.) then
     dw = 0.13*dthr/24.
     sndens = sndens*(1.-dw)+dw
     if (sndens .gt. 0.40) sndens = 0.40
   endif
! ----------------------------------------------------------------------
! calculate snow depth (cm) from snow water equivalent and snow density.
! change snow depth units to meters
! ----------------------------------------------------------------------
   snowhc = esdc/sndens
   snowh = snowhc*0.01
!
   return
   end subroutine noah_snowpack
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   end module lsm_noah_module
