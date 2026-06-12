#include <define.h>   
   module phys_ccm_module
#ifdef CCMCNV
!------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [phys_cps_ccm3]
!      |---[phys_ccm_module]
!               |
!               |--- [ccm3_initailize] ----- [ccm3_constant_setup] *
!               |                       |--- [esinti] *
!               |                               |- [gestbl] *
!               |--- [ccm3_reset] *             |- [gffgch] *
!               |--- [ccm3_virtual_temp] *
!               |--- [ccm3_height8level] *
!               |--- [ccm3_height8interface] *
!               |--- [ccm3_cloud_main] ----- [ccm3_cloud_buoyancy] *
!                                       |--- [ccm3_cloud_trigger] *
!                                       |--- [ccm3_cloud_property] *
!                                       |--- [ccm3_cloud_closure] *
!    [phys_scv_ccm2]                    |--- [ccm3_cloud_trans] *
!      |---[phys_ccm_module]            |--- [ccm3_cloud_q1q2] *
!               |
!               |--- [ccm3_initailize] *
!               |--- [ccm3_reset] *  
!               |--- [ccm3_virtual_temp] *
!               |--- [ccm3_height8level] *
!               |--- [ccm2_saturation_driver] *
!                         |- [ccm2_saturatoin_table1] *
!                         |- [ccm3_cloud_trigger] *
!                         |- [ccm2_saturatoin_table2] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!------------------------------------------------------------------------------
   contains
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
   subroutine ccm3_initialize
!------------------------------------------------------------------------------
   use constant, only  :  g_,hfus_,hvap_,rd_,rhoh2o_,rv_
   use comadj
!------------------------------------------------------------------------------
   rair   = rd_
   cpair  = hvap_
   cappa  = rair/cpair
   epsilo = rd_/rv_
   gravit = g_
   latvap = hvap_
   latice = hfus_
   rhoh2o = rhoh2o_
   clrh2o = latvap/rhoh2o
   cldcp  = latvap/cpair
   rh2o   = rv_
   zvir   = rh2o/rair - 1.

   call ccm3_constant_setup (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o  )
   call esinti (epsilo  ,latvap  ,latice  ,rh2o  ,cpair   )
!
   return
   end subroutine ccm3_initialize
!
!-------------------------------------------------------------------------------
   subroutine ccm3_constant_setup(rair  ,cpair ,gravit ,latvap ,rhowtr  )
!-------------------------------------------------------------------------------
   use comcmf
#include <implicit.h>
!-------------------------------------------------------------------------------
!
! Initialize moist convective mass flux procedure common block, 
!  ccm2_saturation_driver 
! 
!----------------------------Code History-------------------------------
!
! Original version:  J. Hack
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real                 ::  rair              ! gas constant for dry air
   real                 ::  cpair             ! specific heat of dry air
   real                 ::  gravit            ! acceleration due to gravity
   real                 ::  latvap            ! latent heat of vaporization
   real                 ::  rhowtr            ! density of liquid water (STP)
   integer              ::  k              ! vertical level index
!
!-----------------------------------------------------------------------
!
! Initialize physical constants for moist convective mass flux procedure
!
   cp     = cpair         ! specific heat of dry air     
   hlat   = latvap        ! latent heat of vaporization  
   grav   = gravit        ! gravitational constant       
   rgas   = rair          ! gas constant for dry air
   rhoh2o = rhowtr        ! density of liquid water (STP)
!
! Initialize free parameters for moist convective mass flux procedure
!
   c0     = 1.0e-4        ! rain water autoconversion coeff (1/m)
   dzmin  = 0.0           ! minimum cloud depth to precipitate (m)
   betamn = 0.10          ! minimum overshoot parameter
   cmftau = 3600.         ! characteristic adjustment time scale
!    
!     Limit convection to regions below 40 mb
!    
!         if (hypi(1) .ge. 4.e3) then
!           limcnv = 1
!         else
!           do k = 1,levs_
!             if (hypi(k).lt.4.e3 .and. hypi(k+1).ge.4.e3) then
!               limcnv = k
!               goto 10
!             end if
!           end do
!           limcnv = levp1_
!         end if
!      10 if (masterproc) then
!           write(6,*)'ccm3_constant_setup: Convection will be capped at intfc ',limcnv,
!        $            ' which is ',hypi(limcnv),' pascals'
!         end if
   tpmax  = 1.50          ! maximum acceptable t perturbation (deg C)
   shpmax = 1.50e-3       ! maximum acceptable q perturbation (g/g)
   rlxclm = .true.        ! logical variable to specify that relaxation
!                                time scale should applied to column as 
!                                opposed to triplets individually 
!
! Initialize miscellaneous (frequently used) constants
!
   rhlat  = 1.0/hlat      ! reciprocal latent heat of vaporization  
   rcp    = 1.0/cp        ! reciprocal specific heat of dry air     
   rgrav  = 1.0/grav      ! reciprocal gravitational constant       
!
! Initialize diagnostic location information for moist convection scheme
!
   iloc   = 1             ! longitude point for diagnostic info
   jloc   = 1             ! latitude  point for diagnostic info
   nsloc  = 1             ! nstep value at which to begin diagnostics
!
! Initialize other miscellaneous parameters
!
   tiny   = 1.0e-36       ! arbitrary small number (scalar transport)
   eps    = 1.0e-13       ! convergence criteria (machine dependent) 
!
   return
   end subroutine ccm3_constant_setup
!
!-------------------------------------------------------------------------------
   subroutine esinti(epslon  ,latvap  ,latice  ,rh2o    ,cpair   )
!-------------------------------------------------------------------------------
!
! Initialize es lookup tables 
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Hack
! Standardized:      L. Buja, Jun 1992, Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Feb 1996
!
!-------------------------------------------------------------------------------
!
#include <implicit.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real                 ::  epslon          ! Ratio of h2o to dry air molecular weights 
   real                 ::  latvap          ! Latent heat of vaporization
   real                 ::  latice          ! Latent heat of fusion
   real                 ::  rh2o            ! Gas constant for water vapor
   real                 ::  cpair           ! Specific heat of dry air
!
!---------------------------Local workspace-----------------------------
!
   real                 ::  tmn             ! Minimum temperature entry in table
   real                 ::  tmx             ! Maximum temperature entry in table
   real                 ::  trice           ! Trans range from es over h2o to es over ice
   logical              ::  ip           ! Ice phase (true or false)
!
!-----------------------------------------------------------------------
!
! Specify control parameters first
!
   tmn   = 173.16
   tmx   = 375.16
   trice =  20.00
   ip    = .true.
!
! Call gestbl to build saturation vapor pressure table.
!
   call gestbl(tmn     ,tmx     ,trice   ,ip      ,epslon  ,                   &
               latvap  ,latice  ,rh2o    ,cpair   )
!
   return
   end subroutine esinti
!
!-------------------------------------------------------------------------------
   subroutine gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   ,             &
                     latvap  ,latice  ,rh2o    ,cpair   )
!-------------------------------------------------------------------------------
!
! Builds saturation vapor pressure table for later lookup procedure.
! Uses Goff & Gratch (1946) relationships to generate the table
! according to a set of free parameters defined below.  Auxiliary
! routines are also included for making rapid estimates (well with 1%)
! of both es and d(es)/dt for the particular table configuration.
!
! Code history
!
! Original version:  J. Hack
! Standardized:      L. Buja, Jun 1992,  Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Aug 1992
!
!-------------------------------------------------------------------------------
#include "abort.h"
#include <implicit.h>
!
! Input arguments
!
   real                 ::  tmn           ! Minimum temperature entry in es lookup table
   real                 ::  tmx           ! Maximum temperature entry in es lookup table
   real                 ::  epsil         ! Ratio of h2o to dry air molecular weights
   real                 ::  trice         ! Transition range from es over range to es over ice
   real                 ::  latvap        ! Latent heat of vaporization
   real                 ::  latice        ! Latent heat of fusion
   real                 ::  rh2o          ! Gas constant for water vapor
   real                 ::  cpair         ! Specific heat of dry air
!
! Local variables
!
   real                 ::  t             ! Temperature
   integer              ::  n          ! Increment counter
   integer              ::  lentbl     ! Calculated length of lookup table
   integer              ::  itype      ! Ice phase: 0 -> no ice phase
                      !            1 -> ice phase, no transition
                      !           -x -> ice phase, x degree transition
   logical              ::  ip         ! Ice phase logical flag
!
! Statement function
!
! Common block and statement functions for saturation vapor pressure
! look-up procedure, J. J. Hack, February 1990
!
   integer plenest  ! length of saturation vapor pressure table
   parameter (plenest=250)
!
! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
   common/comes/estbl(plenest) ,tmin  ,tmax  ,ttrice ,pcf(6) ,                 &
                epsqs          ,rgasv ,hlatf ,hlatv  ,cp     ,                 &
                icephs
!
   real     ::  estbl      ! table values of saturation vapor pressure
   real     ::  tmin       ! min temperature (K) for table
   real     ::  tmax       ! max temperature (K) for table
   real     ::  ttrice     ! transition range from es over H2O to es over ice
   real     ::  pcf        ! polynomial coeffs -> es transition water to ice
   real     ::  epsqs      ! Ratio of h2o to dry air molecular weights
   real     ::  rgasv      ! Gas constant for water vapor
   real     ::  hlatf      ! Latent heat of vaporization
   real     ::  hlatv      ! Latent heat of fusion
   real     ::  cp         ! specific heat of dry air
   logical  ::  icephs  ! false => saturation vapor press over water only
!
! Dummy variables for statement functions
!
   real     ::  td         ! dummy variable for function evaluation
!-------------------------------------------------------------------------------
!
! Set es table parameters
!
   tmin   = tmn       ! Minimum temperature entry in table
   tmax   = tmx       ! Maximum temperature entry in table
   ttrice = trice     ! Trans. range from es over h2o to es over ice
   icephs = ip        ! Ice phase (true or false)
!
! Set physical constants required for es calculation
!
   epsqs  = epsil
   hlatv  = latvap
   hlatf  = latice
   rgasv  = rh2o
   cp     = cpair
!
   lentbl = ifix(tmax-tmin+2.000001)
   if (lentbl .gt. plenest) then
      write(6,9000) tmax, tmin, plenest
      call MPABORT
   end if
!
! Begin building es table.
! Check whether ice phase requested.
! If so, set appropriate transition range for temperature
!
   if (icephs) then
     if(ttrice.ne.0.0) then
       itype = -ttrice
     else
       itype = 1
     end if
   else
     itype = 0
   end if
!
   t = tmin - 1.0
   do n = 1,lentbl
     t = t + 1.0
     call gffgch(t,estbl(n),itype)
   end do
!
   do n = lentbl+1,plenest
     estbl(n) = -99999.0
   end do
!
! Table complete -- Set coefficients for polynomial approximation of
! difference between saturation vapor press over water and saturation
! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
! is valid in the range -40 < t < 0 (degrees C).
!
!                  --- Degree 5 approximation ---
!
   pcf(1) =  5.04469588506e-01
   pcf(2) = -5.47288442819e+00
   pcf(3) = -3.67471858735e-01
   pcf(4) = -8.95963532403e-03
   pcf(5) = -7.78053686625e-05
!
!                  --- Degree 6 approximation ---
!
!-----pcf(1) =  7.63285250063e-02
!-----pcf(2) = -5.86048427932e+00
!-----pcf(3) = -4.38660831780e-01
!-----pcf(4) = -1.37898276415e-02
!-----pcf(5) = -2.14444472424e-04
!-----pcf(6) = -1.36639103771e-06
!
!
 9000 format('GESTBL: FATAL ERROR *********************************',/,        &
        ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH',             &
        ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/,            &
        ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
   ! 
   return
   end subroutine gestbl
!
!-------------------------------------------------------------------------------
   subroutine gffgch(t       ,es      ,itype   )
!-------------------------------------------------------------------------------
!
! Computes saturation vapor pressure over water and/or over ice using
! Goff & Gratch (1946) relationships.  T (temperature), and itype are
! input parameters, while es (saturation vapor pressure) is an output
! parameter.  The input parameter itype serves two purposes: a value of
! zero indicates that saturation vapor pressures over water are to be
! returned (regardless of temperature), while a value of one indicates
! that saturation vapor pressures over ice should be returned when t is
! less than 273.16 degrees k.  If itype is negative, its absolute value
! is interpreted to define a temperature transition region below 273.16
! degrees k in which the returned saturation vapor pressure is a
! weighted average of the respective ice and water value.  That is, in
! the temperature range 0 => -itype degrees c, the saturation vapor
! pressures are assumed to be a weighted average of the vapor pressure
! over supercooled water and ice (all water at 0 c; all ice at -itype
! c).  Maximum transition range => 40 c
!
!--- Code history
!
! Original version:  J. Hack
! Standardized:      L. Buja, Jun 1992,  Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Feb 1996 
!
!-------------------------------------------------------------------------------
#include "abort.h"
#include <implicit.h>
!
! Arguments
!
! Input arguments
!
   real                 ::  t          ! Temperature
   integer              ::  itype   ! Flag for ice phase and associated transition
!
! Output arguments
!
   real                 ::  es         ! Saturation vapor pressure
!
! Local variables
!
   real                 ::  e1         ! Intermediate scratch variable for es over water
   real                 ::  e2         ! Intermediate scratch variable for es over water
   real                 ::  eswtr      ! Saturation vapor pressure over water
   real                 ::  f          ! Intermediate scratch variable for es over water
   real                 ::  f1         ! Intermediate scratch variable for es over water
   real                 ::  f2         ! Intermediate scratch variable for es over water
   real                 ::  f3         ! Intermediate scratch variable for es over water
   real                 ::  f4         ! Intermediate scratch variable for es over water
   real                 ::  f5         ! Intermediate scratch variable for es over water
   real                 ::  ps         ! Reference pressure (mb)
   real                 ::  t0         ! Reference temperature (freezing point of water)
   real                 ::  term1      ! Intermediate scratch variable for es over ice
   real                 ::  term2      ! Intermediate scratch variable for es over ice
   real                 ::  term3      ! Intermediate scratch variable for es over ice
   real                 ::  tr         ! Transition range for es over water to es over ice
   real                 ::  ts         ! Reference temperature (boiling point of water)
   real                 ::  weight     ! Intermediate scratch variable for es transition
   integer              ::  itypo   ! Intermediate scratch variable for holding itype
!-------------------------------------------------------------------------------
!      
! Check on whether there is to be a transition region for es
!
   if (itype.lt.0) then
     tr    = abs(float(itype))
     itypo = itype
     itype = 1
   else
     tr    = 0.0
     itypo = itype
   end if
   if (tr .gt. 40.0) then
     write(6,900) tr
     call MPABORT
   end if
!
   if(t .lt. (273.16 - tr) .and. itype.eq.1) go to 10
!
! Water
!
   ps = 1013.246
   ts = 373.16
   e1 = 11.344*(1.0 - t/ts)
   e2 = -3.49149*(ts/t - 1.0)
   f1 = -7.90298*(ts/t - 1.0)
   f2 = 5.02808*log10(ts/t)
   f3 = -1.3816*(10.0**e1 - 1.0)/10000000.0
   f4 = 8.1328*(10.0**e2 - 1.0)/1000.0
   f5 = log10(ps)
   f  = f1 + f2 + f3 + f4 + f5
   es = (10.0**f)*100.0
   eswtr = es
!
   if(t.ge.273.16 .or. itype.eq.0) go to 20
!
! Ice
!
10 continue
!
   t0    = 273.16
   term1 = 2.01889049/(t0/t)
   term2 = 3.56654*log(t0/t)
   term3 = 20.947031*(t0/t)
   es    = 575.185606e10*exp(-(term1 + term2 + term3))
!
   if (t.lt.(273.16 - tr)) go to 20
!
! Weighted transition between water and ice
!
   weight = min((273.16 - t)/tr,1.0)
   es = weight*es + (1.0 - weight)*eswtr
!
20 continue
!
   itype = itypo
   return
!
900 format('GFFGCH: FATAL ERROR ******************************',/,             &
          'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR',                &
          ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF',                 &
          ' 40.0 DEGREES C',/, ' TR = ',f7.2)
!
   end subroutine gffgch
!
!-------------------------------------------------------------------------------
   subroutine ccm3_reset(pa      ,kdim    ,pvalue  )
!-------------------------------------------------------------------------------
!
! Reset array pa(kdim) to pvalue
!
!--- Code history
!
! Original version:  CCM1
! Standardized:      L. Bath, Jun 1992
!                    L. Buja, Feb 1996
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer              :: kdim         ! Dimension of array pa
   real                 :: pvalue       ! Value to store in pa
!
! Output arguments
!
   real                 :: pa(kdim)     ! Array to reset
!
!---------------------------Local variable------------------------------
!
   integer              :: j         ! Loop index
!
!-------------------------------------------------------------------------------
!
   do j=1,kdim
     pa(j) = pvalue
   end do
!
   return
   end subroutine ccm3_reset
!
!-------------------------------------------------------------------------------
   subroutine ccm3_virtual_temp(t       ,q       ,zvir    ,tv      )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
!-------------------------------------------------------------------------------
!
! Compute the virtual temperature.
!
!--- Code history
!
! Original version:  B. Boville
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, J. Hack, August 1992
! Reviewed:          D. Williamson, March 1996
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!------------------------------Parameters-------------------------------
!
! Input arguments
!
   real                 ::  t(ILOTS,levs_)       ! temperature
   real                 ::  q(ILOTS,levs_)       ! specific humidity
   real                 ::  zvir                ! virtual temperature constant
!
! Output arguments
!
   real                 ::  tv(ILOTS,levs_)      ! virtual temperature
!
!---------------------------Local storage-------------------------------
!
   integer              ::  i,k              ! longitude and level indexes
!
   do k = 1,levs_
     do i = 1,ILOTS
       tv(i,k) = t(i,k)*(1.0 + zvir*q(i,k))
     end do
   end do
!
   return
   end subroutine ccm3_virtual_temp
!
!-------------------------------------------------------------------------------
   subroutine ccm3_height8level(pstarln ,pmln    ,rair    ,                    &
                                gravit  ,tv      ,z       )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
!-------------------------------------------------------------------------------
!
! Compute the geopotential height *ABOVE THE SURFACE* at layer
! midpoints from the virtual temperatures and pressures.
!
!--- Code history
!
! Original version:  B. Boville, D. Williamson, Jan 1990
! Standardized:      L. Buja, Jun 1992, Feb 1996
! Reviewed:          J. Hack, B. Boville, Aug 1992, Apr 1996
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!------------------------------Parameters-------------------------------
!
! Input arguments
!
   real                 ::  pstarln(ILOTS)        ! Log surface pressures
   real                 ::  pmln(ILOTS,levs_)      ! Log midpoint pressures
   real                 ::  rair                  ! Gas constant for dry air
   real                 ::  gravit                ! Acceleration of gravity
   real                 ::  tv(ILOTS,levs_)        ! Virtual temperature
!
! Output arguments
!
   real                 ::  z(ILOTS,levs_)         ! Height above surface at midpoints
!
!---------------------------Local variables-----------------------------
!
   integer              ::  i,k,l              ! Lon, level, level indices
   real                 ::  rog                   ! Rair / gravit
!
!-------------------------------------------------------------------------------
!
! Diagonal term of hydrostatic equation
!
   rog = rair/gravit
   do k = 1,levs_-1
     do i = 1,ILOTS
       z(i,k) = rog*tv(i,k)*0.5*(pmln(i,k+1) - pmln(i,k))
     end do
   end do
   do i = 1,ILOTS
     z(i,levs_) = rog*tv(i,levs_)*(pstarln(i) - pmln(i,levs_))
   end do
!
! Bottom level term of hydrostatic equation
!
   do  k=1,levs_-1
     do i = 1,ILOTS
       z(i,k) = z(i,k) + rog*tv(i,levs_)*(pstarln(i) -                         &
                         0.5*(pmln(i,levs_-1) + pmln(i,levs_)))
     end do
   end do
!
! Interior terms of hydrostatic equation
!
   do k = 1,levs_-2
     do l = k+1,levs_-1
       do i = 1,ILOTS
         z(i,k) = z(i,k) + rog*(tv(i,l)) *                                     &
                           0.5*(pmln(i,l+1) - pmln(i,l-1))
       end do
     end do
   end do
!
   return
   end subroutine ccm3_height8level
!
!-------------------------------------------------------------------------------
   subroutine ccm3_height8interface(piln    ,pmln    ,rair    ,                &
                                    gravit  ,tv      ,zm      ,zi      )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
#include <implicit.h>
!-------------------------------------------------------------------------------
!
! Compute the geopotential height at the interface points from *HEIGHT
! ABOVE THE SURFACE* at midlayer midpoints, using the supplied virtual
! temperatures and pressures.
!
!--- Code history
!
! Original version:  D. Williamson, J. Hack
! Standardized:      L. Buja, Jun 1992, Feb 1996
! Reviewed:          J. Hack, B. Boville, Aug 1992, Apr 1996
!
!------------------------------Commons----------------------------------
!
! Input arguments
!
   real                 ::  piln(ILOTS,levs_+1)    ! Log interface pressures
   real                 ::  pmln(ILOTS,levs_)     ! Log midpoint pressures
   real                 ::  rair                 ! Gas constant for dry air
   real                 ::  gravit               ! Acceleration of gravity
   real                 ::  tv(ILOTS,levs_)       ! Virtual temperature
   real                 ::  zm(ILOTS,levs_)       ! Height above surface at midpoints
!
! Output arguments
!
   real                 ::  zi(ILOTS,levs_+1)      ! Height above surface at interfaces
!
!---------------------------Local variables-----------------------------
!
   integer              ::  i,k               ! Lon, level indices
   real                 ::  rog                  ! Rair / gravit
!
!-----------------------------------------------------------------------
!
! Add increment to midlayer height
!
   rog = rair/gravit
   do i = 1,ILOTS
     zi(i,1) = zm(i,1) + rog*(pmln(i,1) - piln(i,1))*tv(i,1)
   end do
!
   do k = 2,levs_
     do i = 1,ILOTS
       zi(i,k) = zm(i,k) + rog*(pmln(i,k) - piln(i,k))*0.5*                    &
                 (2.0*tv(i,k) - (tv(i,k) - tv(i,k-1))/                         &
                  (pmln(i,k) - pmln(i,k-1))*                                   &
                  (pmln(i,k) - piln(i,k)))
     end do
   end do
!
! The surface height is zero by definition.
!
   do i = 1,ILOTS
     zi(i,levs_+1) = 0.0
   end do
!
   return
   end subroutine ccm3_height8interface
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_main(t       ,qh      ,pcpc    ,jctop   ,jcbot   ,    &
                       pblh    ,zm      ,geos    ,zi      ,qtg     ,           &
                       ttg     ,pap     ,paph    ,dpp     ,ts      ,           &
                       delt    ,mcon    ,cme     ,nstep   ,lat     ,           &
                       tpert   ,qpert   ,dlf     ,pflx    ,cape)
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use guang
   use comcmf
!-------------------------------------------------------------------------------
!
! This is contributed code not fully standardized by the CCM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version 
! of the CCM where it will include a more straightforward formulation 
! and will make use of the standard CCM nomenclature
!
! same as conv.up except saturation vapor pressure is calculated
! in a different way.
!
! jul 17/92 - guang jun zhang, m.lazare. 
!             calls new ccm3_cloud_buoyancy, ccm3_cloud_q1q2
!             and moment (several new work fields added for later).
!
! nov 21/91 - m.lazare. like previous conv except calls new
!                       clpdprp.
! feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             previous version conv.
! performs deep convective adjustment based on mass-flux closure
! algorithm.
!
! ************************ index of variables *********************!
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  wg * betad    downward mass flux at cloud base.
!  wg * betau    upward   mass flux at cloud base.
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      local sigma half-level thickness (i.e. dshj).
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * levs_     number of model levels.
!  ic  * ilg      lon+2 = size of grid slice.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlat     gaussian latitude index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  ic  * kount    current model timestep number.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  wg * long     gathered values of lon.
!  ic  * lev      levs_+1.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified 
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  c  * kups     number of points undergoing deep convection.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  i/o * pcp      row of precipitable water in metres.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  i  * pressg   row of surface pressure in pa.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qdb      row of qd at cloud base.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sdb      row of sd at cloud base.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  i  * shbj     grid slice of local bottom interface sigma values.
!  i  * shj      grid slice of local half-level sigma values.
!  i  * shtj     row of local top interfaces of first level.
!  wg * su       grid slice of dry static energy in updraft.
!  wg * sumde    row of vertically-integrated moist static energy 
!                change.
!  wg * sumdq    row of vertically-integrated scaled mixing ratio 
!                change.
!  wg * sumdt    row of vertically-integrated dry static energy change.
!  wg * sumq     row of vertically-integrated mixing ratio change.
!  i/o * t        grid slice of temperature at mid-layer.
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  o  * jcbot    row of base of cloud indices passed out.
!  w  * tf       grid slice of temperature at interface.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!
! multi-level i/o fields:
!
! input/output arguments:
!
   real     ::  t(ILOTS,levs_) 
   real     ::  qh(ILOTS,levs_) 
!      real utg(ILOTS,levs_) 
!      real vtg(ILOTS,levs_) 
   real     ::  qtg(ILOTS,levs_) 
   real     ::  ttg(ILOTS,levs_)
!
!       input arguments
!
   real     ::  pap(ILOTS,levs_) 
   real     ::  paph(ILOTS,levs_+1) 
   real     ::  dpp(ILOTS,levs_) 
   real     ::  zm(ILOTS,levs_) 
   real     ::  geos(ILOTS) 
   real     ::  zi(ILOTS,levs_+1)
   real     ::  pblh(ILOTS) 
   real     ::  zs(ILOTS) 
   real     ::  tpert(ILOTS) 
   real     ::  qpert(ILOTS)
!
!       output arguments
!
   real     ::  pcpck(ILOTS,levs_)
   real     ::  mup(ILOTS,levs_) 
   real     ::  mdn(ILOTS,levs_) 
   real     ::  mcon(ILOTS,levs_) 
   real     ::  dlg(ILOTS,levs_)    ! gathered version of the detraining cld h2o tend
   real     ::  dsdt2(ILOTS,levs_)  ! scattered version of the temp tend
   real     ::  dqdt2(ILOTS,levs_)  ! scattered version of the q tend
   real     ::  dlg2(ILOTS,levs_)   ! gathered version of the detraining cld h2o tend
   real     ::  dlf(ILOTS,levs_)    ! scattered version of the detraining cld h2o tend
   real     ::  pflx(ILOTS,levs_)   ! scattered precip flux at each level
   real     ::  pflxg(ILOTS,levs_)  ! gather precip flux at each level
   real     ::  cu(ILOTS,levs_)     ! scattered condensation rate
   real     ::  cug(ILOTS,levs_)    ! gathered condensation rate 
   real     ::  evpg(ILOTS,levs_)   ! gathered evap rate of rain in downdraft
   real     ::  evp(ILOTS,levs_)    ! scattered evap rate of rain in downdraft
   real     ::  mumax(ILOTS) 
   real     ::  cme(ILOTS,levs_)
   real     ::  small, msn, ms(ILOTS)
   integer  ::  km1, kp1
!
! single-level i/o fields:
!       input arguments
!
   real     ::  ts(ILOTS) 
   real     ::  pblt(ILOTS)
!
!       input/output arguments:
!
   real     ::  paprc(ILOTS) 
   real     ::  paprs(ILOTS) 
   real     ::  zctopm(ILOTS)
!
!       output arguments:
!
   real     ::  jctop(ILOTS) 
   real     ::  jcbot(ILOTS)
   real     ::  pcpr(ILOTS) 
   real     ::  pcps(ILOTS) 
   real     ::  pcpc(ILOTS)
!
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
   real     ::  u(ILOTS,levs_) 
   real     ::  v(ILOTS,levs_) 
   real     ::  q(ILOTS,levs_) 
   real     ::  p(ILOTS,levs_) 
   real     ::  z(ILOTS,levs_) 
   real     ::  s(ILOTS,levs_) 
   real     ::  qh0(ILOTS,levs_) 
   real     ::  tp(ILOTS,levs_) 
   real     ::  zf(ILOTS,levs_+1) 
   real     ::  pf(ILOTS,levs_+1) 
   real     ::  qstp(ILOTS,levs_) 

   real     ::  cape(ILOTS) 
   real     ::  tl(ILOTS) 
   real     ::  sumq(ILOTS) 
   real     ::  pcpdh(ILOTS)
!      real sumdt(ILOTS) 
!      real sumdq(ILOTS) 
!      real sumde(ILOTS)

   integer  ::  lcl(ILOTS) 
   integer  ::  lel(ILOTS) 
   integer  ::  lon(ILOTS) 
   integer  ::  maxi(ILOTS) 
   integer  ::  ideep(ILOTS) 
   integer  ::  index(ILOTS)
   real     ::  precip
!
! gathered work fields:
!
   real     ::  qg(ILOTS,levs_) 
   real     ::  tg(ILOTS,levs_) 
   real     ::  pg(ILOTS,levs_) 
   real     ::  zg(ILOTS,levs_) 
   real     ::  sg(ILOTS,levs_) 
   real     ::  tpg(ILOTS,levs_) 
   real     ::  zfg(ILOTS,levs_+1) 
   real     ::  qstpg(ILOTS,levs_) 
   real     ::  ug(ILOTS,levs_) 
   real     ::  vg(ILOTS,levs_) 
   real     ::  cmeg(ILOTS,levs_)

   real     ::  capeg(ILOTS) 
   real     ::  tlg(ILOTS)

   integer  ::  lclg(ILOTS) 
   integer  ::  lelg(ILOTS) 
   integer  ::  maxg(ILOTS)
!
! work fields arising from gathered calculations.
!
   real     ::  mu(ILOTS,levs_) 
   real     ::  eu(ILOTS,levs_) 
   real     ::  dqdt(ILOTS,levs_) 
   real     ::  dsdt(ILOTS,levs_) 
   real     ::  du(ILOTS,levs_) 
   real     ::  md(ILOTS,levs_) 
   real     ::  ed(ILOTS,levs_) 
   real     ::  alpha(ILOTS,levs_) 
   real     ::  sd(ILOTS,levs_) 
   real     ::  qd(ILOTS,levs_) 
   real     ::  mc(ILOTS,levs_) 
   real     ::  qhat(ILOTS,levs_) 
   real     ::  qu(ILOTS,levs_) 
   real     ::  su(ILOTS,levs_) 
   real     ::  qs(ILOTS,levs_) 
   real     ::  shat(ILOTS,levs_) 
   real     ::  dp(ILOTS,levs_) 
   real     ::  hmn(ILOTS,levs_) 
   real     ::  hsat(ILOTS,levs_) 
   real     ::  ql(ILOTS,levs_) 
   real     ::  dudt(ILOTS,levs_) 
   real     ::  dvdt(ILOTS,levs_) 
   real     ::  ud(ILOTS,levs_) 
   real     ::  vd(ILOTS,levs_)

!      real deltat(ILOTS,levs_) 
!      real deltaq(ILOTS,levs_)

   real     ::  betau(ILOTS) 
   real     ::  betad(ILOTS) 
   real     ::  qdb(ILOTS) 
   real     ::  sdb(ILOTS) 
   real     ::  dsubcld(ILOTS) 
   real     ::  mb(ILOTS) 
   real     ::  totpcp(ILOTS) 
   real     ::  totevp(ILOTS)
   real     ::  mu2(ILOTS,levs_) 
   real     ::  eu2(ILOTS,levs_) 
   real     ::  du2(ILOTS,levs_) 
   real     ::  md2(ILOTS,levs_) 
   real     ::  ed2(ILOTS,levs_)
   real     ::  denom
   integer  ::  jt(ILOTS) 
   integer  ::  jlcl(ILOTS) 
   integer  ::  j0(ILOTS) 
   integer  ::  jd(ILOTS)
   real     ::  capelmt
   real     ::  cpres
   real     ::  delt
   real     ::  fixer
   integer  ::  i
   integer  ::  ii
   integer  ::  k
   integer  ::  lat
   integer  ::  lengath
   integer  ::  msg
   integer  ::  nstep
   real     ::  psdiss
   real     ::  psevap
   real     ::  psheat
   real     ::  psrain
   real     ::  qdifr
   real     ::  qeff
   real     ::  qmin
   real     ::  rl
   real     ::  sdifr
   real     ::  diffc(ILOTS),diffl(ILOTS)
!
!--------------------------Data statements------------------------------
!
   logical  ::  momentm
   data capelmt/70./
   data momentm/.FALSE./
!
! set up convection cap, this was moved from ccm3_constant_setup to here to
! get rid of dependency on hybrid pressure coefficients for the
! Suarez repos code
!
   if (pap(1,1) .ge. 4.e3) then 
      limcnv =1
   else
     do k = 1,levs_-1
       if (pap(1,k).lt.4.e3 .and. pap(1,k+1).ge.4.e3) then
         limcnv = k
!
         goto 1010
!
       end if
     end do
     limcnv = levs_+1
   end if
!
1010 continue
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1
!
! initialize necessary arrays.
! zero out variables not used in ccm
!
   do i = 1,ILOTS
     paprc(i) = 0.
     paprs(i) = 0.
     zctopm(i) = 0.
   end do
   psdiss = 0.
   psheat = 0.
   psevap = 0.
   psrain = 0.
!      jlatpr = 32
   cpres = 1004.64
   a = 21.656
   b = 5418.
   c1 = 6.112
   c2 = 17.67
   c3 = 243.5
   eps1 = 0.622
   tfreez = 273.16
   qmin = 1.E-20
   rl = 2.5104E6
!
! initialize convective tendencies
!
   do k = 1,levs_
     do i = 1,ILOTS
       dqdt(i,k) = 0.
       dsdt(i,k) = 0.
       dudt(i,k) = 0.
       dvdt(i,k) = 0.
!            deltaq(i,k) = qh(i,k)
!            deltat(i,k) = t(i,k)
       dqdt2(i,k) = 0.
       dsdt2(i,k) = 0.
       pcpck(i,k) = 0.
       pflx(i,k) = 0.
       pflxg(i,k) = 0.
     end do
   end do
   if (.not.momentm) then
     do k = msg + 1,levs_
       do i = 1,ILOTS
         u(i,k) = 0.
         v(i,k) = 0.
       end do
     end do
   end if
!
   do i = 1,ILOTS
     pblt(i) = levs_
     pcpr(i) = 0.
     pcps(i) = 0.
     dsubcld(i) = 0.
     sumq(i) = 0.
!         sumdt(i) = 0.
!         sumdq(i) = 0.
     pcpdh(i) = rgrav
     jctop(i) = levs_
     jcbot(i) = 1
   end do
!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ILOTS
     zs(i) = geos(i)*rgrav
     pf(i,levs_+1) = paph(i,levs_+1)*0.01
     zf(i,levs_+1) = zi(i,levs_+1) + zs(i)
   end do
   do k = 1,levs_
     do i = 1,ILOTS
       p(i,k) = pap(i,k)*0.01
       pf(i,k) = paph(i,k)*0.01
       z(i,k) = zm(i,k) + zs(i)
       zf(i,k) = zi(i,k) + zs(i)
     end do
   end do
!
   diffl(:)=pblh(:)
   do k = levs_ - 1,msg + 1,-1
     do i = 1,ILOTS
       diffc(i)=z(i,k)-zs(i)
       if(diffc(i).ge.pblh(i).and.diffl(i).le.pblh(i)) pblt(i)=k
       diffl(i)=diffc(i)
!            if (abs(z(i,k)-zs(i)-pblh(i)).lt.
!     $          (zf(i,k)-zf(i,k+1))*0.5) pblt(i) = k
     end do
   end do
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! convert from specific humidity (bounded by qmin) to mixing ratio.
! define dry static energy (normalized by cp).
!
   do k = 1,levs_
     do i = 1,ILOTS
       qh0(i,k) = qh(i,k)
       qeff = max(qh(i,k),qmin)
       q(i,k) = qeff
       s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
       tp(i,k)=0.0
       shat(i,k) = s(i,k)
       qhat(i,k) = q(i,k)
       dp(i,k) = dpp(i,k)*0.01
       qg(i,k) = q(i,k)
       tg(i,k) = t(i,k)
       pg(i,k) = p(i,k)
       zg(i,k) = z(i,k)
       sg(i,k) = s(i,k)
       tpg(i,k) = tp(i,k)
       zfg(i,k) = zf(i,k)
       qstpg(i,k) = q(i,k)
       ug(i,k) = u(i,k)
       vg(i,k) = v(i,k)
       dlg(i,k) = 0
       dlg2(i,k) = 0
       dlf(i,k) = 0
       cu(i,k) = 0
       evp(i,k) = 0
     end do
   end do
   do i = 1,ILOTS
     zfg(i,levs_+1) = zf(i,levs_+1)
     capeg(i) = 0.
     lclg(i) = 1
     lelg(i) = levs_
     maxg(i) = 1
     tlg(i) = 400.
     dsubcld(i) = 0.
     qdb(i) = 0.
     sdb(i) = 0.
     betau(i) = 0.
     betad(i) = 0.
   end do
!
! evaluate covective available potential energy (cape).
!
   call ccm3_cloud_buoyancy     (q       ,t       ,p        ,                  &
               z       ,pf      ,tp      ,qstp    ,tl       ,                  &
               rl      ,cape    ,pblt    ,lcl     ,lel      ,                  &
               lon     ,maxi    ,rgas    ,grav    ,cpres    ,                  &
               msg     ,nstep   ,lat     ,tpert   ,qpert)
!
! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).
!
   call ccm3_cloud_trigger(ILOTS,cape,1,capelmt,index,lengath)
   if (lengath.eq.0) return
!DIR$ IVDEP
   do ii=1,lengath
     i=index(ii)
     ideep(ii)=i
   end do
!c        jyes = 0
!c        jno = ILOTS - 1 + 2
!c        do il = 1,ILOTS
!c           if (cape(il).gt.capelmt) then
!c              jyes = jyes + 1
!c              ideep(jyes) = il
!c           else
!c              jno = jno - 1
!c              ideep(jno) = il
!c           end if
!c        end do
!c        lengath = jyes
!c        if (lengath.eq.0) return
!
! obtain gathered arrays necessary for ensuing calculations.
!
   do k = 1,levs_
     do i = 1,lengath
       dp(i,k) = 0.01*dpp(ideep(i),k)
       qg(i,k) = q(ideep(i),k)
       tg(i,k) = t(ideep(i),k)
       pg(i,k) = p(ideep(i),k)
       zg(i,k) = z(ideep(i),k)
       sg(i,k) = s(ideep(i),k)
       tpg(i,k) = tp(ideep(i),k)
       zfg(i,k) = zf(ideep(i),k)
       qstpg(i,k) = qstp(ideep(i),k)
       ug(i,k) = u(ideep(i),k)
       vg(i,k) = v(ideep(i),k)
     end do
   end do
!
   do i = 1,lengath
     zfg(i,levs_+1) = zf(ideep(i),levs_+1)
   end do
   do i = 1,lengath
     capeg(i) = cape(ideep(i))
     lclg(i) = lcl(ideep(i))
     lelg(i) = lel(ideep(i))
     maxg(i) = maxi(ideep(i))
     tlg(i) = tl(ideep(i))
   end do
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = msg + 1,levs_
     do i = 1,lengath
       if (k.ge.maxg(i)) then
         dsubcld(i) = dsubcld(i) + dp(i,k)
       end if
     end do
   end do
!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,levs_
     do i = 1,lengath
       alpha(i,k) = 0.5
       sdifr = 0.
       qdifr = 0.
       if (sg(i,k).gt.0. .or. sg(i,k-1).gt. 0.)                                &
             sdifr = abs((sg(i,k)-sg(i,k-1))/                                  &
             max(sg(i,k-1),sg(i,k)))
       if (qg(i,k).gt.0. .or. qg(i,k-1).gt.0.)                                 &
             qdifr = abs((qg(i,k)-qg(i,k-1))/                                  &
             max(qg(i,k-1),qg(i,k)))
       if (sdifr.gt.1.E-6) then
         shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/                 &
                        (sg(i,k-1)-sg(i,k))
       else
         shat(i,k) = 0.5* (sg(i,k)+sg(i,k-1))
       end if
       if (qdifr.gt.1.E-6) then
         qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/                 &
                        (qg(i,k-1)-qg(i,k))
       else
         qhat(i,k) = 0.5* (qg(i,k)+qg(i,k-1))
       end if
     end do
   end do
!
! obtain cloud properties.
!
   call ccm3_cloud_property (qg     ,tg     ,ug     ,vg      ,pg     ,         &
                             zg     ,sg     ,mu     ,eu      ,du     ,         &
                             md     ,ed     ,sd     ,qd      ,ud     ,         &
                             vd     ,mc     ,qu     ,su      ,zfg    ,         &
                             qs     ,hmn    ,hsat   ,alpha   ,shat   ,         &
                             ql     ,totpcp ,totevp ,cmeg    ,maxg   ,         &
                             lelg   ,jt     ,jlcl   ,maxg    ,j0     ,         &
                             jd     ,rl     ,1      ,lengath ,rgas   ,         &
                             grav   ,cpres  ,msg    ,nstep   ,lat    ,         &
                             pflxg  ,evpg   ,cug    ,mu2     ,eu2    ,         &
                             du2    ,md2    ,ed2    ) 
!
! determine cloud base mass flux.
!
   do i = 1,lengath
     qdb(i) = qd(i,maxg(i))
     sdb(i) = sd(i,maxg(i))
     betad(i) = md(i,maxg(i))
     betau(i) = mu(i,maxg(i))
   end do
!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,levs_
     do i = 1,lengath
       du(i,k) = du(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       eu(i,k) = eu(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       ed(i,k) = ed(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       cug(i,k) = cug(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       evpg(i,k) = evpg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       du2(i,k) = du2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       eu2(i,k) = eu2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       ed2(i,k) = ed2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
     end do
   end do

   call ccm3_cloud_closure   (qg      ,tg      ,pg      ,zg      ,             &
                     sg      ,tpg     ,qs      ,qu      ,su      ,             &
                     mc      ,du      ,mu      ,md      ,qd      ,             &
                     sd      ,alpha   ,qhat    ,shat    ,dp      ,             &
                     qstpg   ,zfg     ,ql      ,dsubcld ,mb      ,             &
                     capeg   ,tlg     ,lclg    ,lelg    ,jt      ,             &
                     maxg    ,1       ,lengath ,rgas    ,grav    ,             &
                     cpres   ,rl      ,msg     ,capelmt ,nstep   ,             &
                     lat     )
!
! limit cloud base mass flux to theoretical upper bound.
!
   if (.true.) then
     do k = msg + 2,levs_
       do i = 1,lengath
         if (mu(i,k).gt.0.) then
           mb(i) = min(mb(i),dp(i,k)/ (2.*delt*mu(i,k)))
         end if
       end do
     end do
   else
!
     do i = 1,lengath
       mumax(i) = 0
     end do
!
     do k = msg + 2,levs_
       do i = 1,lengath
         mumax(i)  = max(mumax(i), mu(i,k)/dp(i,k))
       end do
     end do
!
     do i = 1,lengath
       if (mumax(i).gt.0.) then
         mb(i) = min(mb(i),0.5/(delt*mumax(i)))
       else
         mb(i) = 0.
       endif
     end do
   endif
!
   do k = msg + 1,levs_
     do i = 1,lengath
       mu(i,k) = mu(i,k)*mb(i)
       md(i,k) = md(i,k)*mb(i)
       mc(i,k) = mc(i,k)*mb(i)
       du(i,k) = du(i,k)*mb(i)
       eu(i,k) = eu(i,k)*mb(i)
       ed(i,k) = ed(i,k)*mb(i)
       cmeg(i,k) = cmeg(i,k)*mb(i)
       cug(i,k) = cug(i,k)*mb(i)
       evpg(i,k) = evpg(i,k)*mb(i)
       pflxg(i,k) = pflxg(i,k)*mb(i)*100./grav
       mu2(i,k) = mu2(i,k)*mb(i)
       md2(i,k) = md2(i,k)*mb(i)
       du2(i,k) = du2(i,k)*mb(i)
       eu2(i,k) = eu2(i,k)*mb(i)
       ed2(i,k) = ed2(i,k)*mb(i)
     end do
   end do
#ifdef DBG
! perform some consistency checks on updraft quantities
   small = 1.e-18
   do i = 1,lengath
     ms(i) = 0
   end do
   do k = levs_,1,-1
     kp1 = min(k+1,levs_)
     do i = 1,lengath
       msn = ms(i) + (eu2(i,k)-du2(i,k))*dp(i,k)
       if (abs(msn-mu2(i,k))/(abs(msn)+abs(mu2(i,k))+small)                    &
              .gt.1.e-10.and.mu2(i,k).ne.0.) then
         write (6,*) ' updraft consistency ', i, k, lat, msn,                  &
                 ms(i), mu2(i,k), mu2(i,kp1), eu2(i,k), du2(i,k),              &
                 (zfg(i,k)-zfg(i,k+1)), dp(i,k)
         write (6,*) ' dmu ', msn-ms(i), mu2(i,k)-mu2(i,kp1)
         write (6,*) ' jb or maxg ', maxg(i)
       endif
       ms(i) = msn
     end do
   end do
! perform some consistency checks on downdraft quantities
   do i = 1,lengath
     ms(i) = 0
   end do
   do k = 1,levs_
     kp1 = min(k+1,levs_)
     do i = 1,lengath
       msn = ms(i) - (ed2(i,k))*dp(i,k)
       denom = max(abs(msn),abs(md2(i,kp1)),abs(ed2(i,k)*dp(i,k)))
       if (     abs(msn-md2(i,kp1)).gt.1.e-10*denom                            &
             .and.denom.gt.1.e-15                                              &
            ) then
         write (6,*) ' inconsistent downdraft ', i, k, lat, msn,               &
                 ms(i), md2(i,kp1), md2(i,k), ed2(i,k)
         write (6,*) ' jt, jd, jb ', jt(i), jd(i), maxg(i)
         write (6,*) ' rdiff is ',                                             &
         abs(msn-md2(i,kp1))/(abs(msn)+abs(md2(i,kp1))+1.e-16)
       endif
       ms(i) = msn
     end do
   end do
!      if (lat.eq.32) then
!         stop
!      endif
!
#endif
   do i = 1,lengath
     betau(i) = betau(i)*mb(i)
     betad(i) = betad(i)*mb(i)
!
! totpcp from rad_cloud_comp has the dimension of kg/kg, here it is 
! converted to kg/(m^2*s), the precipitation rate
!
     totpcp(i) = totpcp(i)*mb(i)*100./grav
     totevp(i) = totevp(i)*mb(i)*100./grav
   end do
!
! compute temperature and moisture changes due to convection.
!
   call ccm3_cloud_q1q2(dqdt    ,dsdt    ,qg      ,sg      ,qs      ,          &
                        qu      ,su      ,mc      ,du      ,alpha   ,          &
                        qhat    ,shat    ,dp      ,mu      ,md      ,          &
                        sd      ,qd      ,ql      ,dsubcld ,qdb     ,          &
                        sdb     ,betau   ,betad   ,mb      ,lclg    ,          &
                        jt      ,maxg    ,2.*delt ,1       ,lengath ,          &
                        cpres   ,rl      ,msg     ,nstep   ,lat     ,          &
                        dlg     ,cug     )
!
! compute momentum changes due to convection, if desired (i.e
! if logical switch set).
!
!
!      if(momentm)                                                   
!       then
!        call moment(dudt,dvdt,du,alpha,dp,ed,eu,mc,md,mu,
!     1             pg,qd,qu,qhat,sd,su,shat,ud,vd,tg,ug,vg,zg,zfg,
!     2             dsubcld,maxg,jd,jt,rl,
!     3             msg,2.*delt,grav,cpres,rgas,levs_,1,lengath,ILOTS,lat)
!      endif
!
!      if (.false.) then
   if (1.gt.1) then
     call ccm3_cloud_trans                                                     &
                  (qh      ,mu2     ,md2     ,                                 &
                   du2     ,eu2     ,ed2     ,dp      ,                        &
                   dsubcld ,jt      ,maxg    ,ideep   ,                        &
                   1       ,lengath ,nstep   ,lat     ,                        &
                   delt    )
   endif
!
! gather back temperature and mixing ratio.
!
   do k = msg + 1,levs_
     do i = 1,lengath
       psdiss = psdiss + (dudt(i,k)*u(ideep(i),k)+                             &
                  dvdt(i,k)*v(ideep(i),k))*dpp(ideep(i),k)/grav
!
! q is updated to compute net precip, and then reset to old value.
! the last line is overwritten. so the input basic variables, i.e.
! q, t, u and v are updated to include convective increments. 
! (5/2/95)
!
       q(ideep(i),k) = q(ideep(i),k) + 2.*delt*dqdt(i,k)
       t(ideep(i),k) = t(ideep(i),k) + 2.*delt*dsdt(i,k)
       u(ideep(i),k) = u(ideep(i),k) + 2.*delt*dudt(i,k)
       v(ideep(i),k) = v(ideep(i),k) + 2.*delt*dvdt(i,k)
       cme(ideep(i),k) = cmeg(i,k)
       mup(ideep(i),k) = mu(i,k)
       mdn(ideep(i),k) = md(i,k)
       mcon(ideep(i),k) = mc(i,k)
       qtg(ideep(i),k) = dqdt(i,k)
       ttg(ideep(i),k) = dsdt(i,k)
!            utg(ideep(i),k) = dudt(i,k)
!            vtg(ideep(i),k) = dvdt(i,k)
       dlf(ideep(i),k) = dlg(i,k)
       pflx(ideep(i),k) = pflxg(i,k)
       cu(ideep(i),k) = cug(i,k)
       evp(ideep(i),k) = evpg(i,k)
     end do
   end do
!
   do i = 1,lengath
     jctop(ideep(i)) = jt(i)
     jcbot(ideep(i)) = jlcl(i)
     psevap = psevap + totevp(i)
     psrain = psrain + totpcp(i)
   end do
!
! convert back to specific humidity from mixing ratio.
! take into account any moisture added to ensure positiveness
! of specific humidity at start of routine.
!
   do k = msg + 1,levs_
     do i = 1,ILOTS
       qh(i,k) = q(i,k)
       qh(i,k) = qh(i,k) - max((qmin-qh0(i,k)),0.)
     end do
   end do
!
! determine change in column storage due to deep convection.
! enforce water vapor conservation if negative precip diagnosed
!
   do k = levs_,msg + 1,-1
     do i = 1,ILOTS
       sumq(i) = sumq(i) - dpp(i,k)* (qh(i,k)-qh0(i,k))
       if (sumq(i).lt.0.0) then
         fixer = -sumq(i)/dpp(i,k)
         qh(i,k) = qh(i,k) - fixer
         t(i,k) = t(i,k) + (rl/cpres)*fixer
         sumq(i) = 0.0
!
         qtg(i,k) = qtg(i,k) - fixer/ (2.0*delt)
         ttg(i,k) = ttg(i,k) + (rl/cpres)*fixer/ (2.0*delt)
       end if
       pcpck(i,k) = max(0.,sumq(i))
     end do
   end do
!
! obtain final precipitation rate.
!
   do i = 1,ILOTS
!         llo1 = ts(i) .ge. tfreez
!
! here pcpr and pcps are in units of kg/m^2, ie. precip per
! time step
!
!         pcpr(i) = cvmgt(pcpdh(i)*max(sumq(i),0.),0.,llo1)
!         pcps(i) = cvmgt(0.,pcpdh(i)*max(sumq(i),0.),llo1)
     precip = pcpdh(i)*max(sumq(i),0.)
     if (ts(i) .ge. tfreez) then
       pcpr(i) = precip
       pcps(i) = 0.
     else
       pcpr(i) = 0.
       pcps(i) = precip
     end if
   end do

!
! accumulate precipitation, the 1000. is the density of water, so
! paprc and paprs are now in units of meters.
!
   do i = 1,ILOTS
     paprc(i) = paprc(i) + (pcpr(i)+pcps(i))/1000.
     paprs(i) = paprs(i) + pcps(i)/1000.
   end do
!
! convert precipitation to m/s, ie, precip rate.
!
   do i = 1,ILOTS
     pcpr(i) = pcpr(i)/ (2.*delt)/1000.
     pcps(i) = pcps(i)/ (2.*delt)/1000.
     pcpc(i) = pcpr(i) + pcps(i)
     psheat = psheat + (pcps(i)+pcpr(i))*rl
   end do
   do k = msg + 1,levs_
     do i = 1,ILOTS
       pcpck(i,k) = pcpdh(i)*pcpck(i,k)/ (2.*delt)
     end do
   end do
!
! calculate conservation of quantities.
!

!       if(lat.eq.jlatpr)then
!        do l = msg+1,levs_
!        do i = 1,lengath
!          sumdq(i) = sumdq(i) + 2.*delt*(rl/cpres)*dpp(ideep(i),l)!
!     1                          dqdt(i,l)
!          sumdt(i) = sumdt(i) + 2.*delt*dpp(ideep(i),l)*dsdt(i,l)
!        end do
!        end do
!
!        write(6,*)'sumdq,sumdt,sumde in convection subroutine########'
!        do i = 1,lengath
!          sumde(i) = sumdt(i) + sumdq(i)
!          write(6, 901) sumdq(i), sumdt(i),sumde(i), i, ideep(i)
!        end do
!c
!        write(6,*)'sumdq,sumdt,sumde ... all points'
!      do i = 1,ILOTS
!         sumdq(i) = 0.0
!         sumdt(i) = 0.0
!      end do
!c
!      do l = msg+1,levs_
!      do i = 1,ILOTS
!        deltaq(i,l) = qh(i,l) - deltaq(i,l)
!        deltat(i,l) = t (i,l) - deltat(i,l)
!        sumdq(i) = sumdq(i) + (rl/cpres)*dpp(i,l)*deltaq(i,l)
!        sumdt(i) = sumdt(i) + dpp(i,l)*deltat(i,l)
!      end do
!      end do
!      do i = 1,ILOTS
!        sumde(i) = sumdt(i) + sumdq(i)
!      end do
!        write(6, 902) (i,sumdq(i),sumdt(i),sumde(i),i=1,ILOTS)
!
!  901   format(1x,3e20.12, i10, i10)
!  902   format(1x,i10, 3e20.12)
!
!      endif

   return
   end subroutine ccm3_cloud_main
!
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_buoyancy     (q       ,t       ,p       ,             &
                     z       ,pf      ,tp      ,qstp    ,tl      ,             &
                     rl      ,cape    ,pblt    ,lcl     ,lel     ,             &
                     lon     ,mx      ,rd      ,grav    ,cp      ,             &
                     msg     ,nstep   ,lat     ,tpert   ,qpert   )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use guang
!-------------------------------------------------------------------------------
! This is contributed code not fully standardized by the CCM core group.
!
! the documentation has been enhanced to the degree that we are able
!
! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996
!
! jul 14/92 - guang jun zhang, m.lazare, n.mcfarlane.  as in
!             previous version ccm3_cloud_buoyancy except remove pathalogical
!             cases of "zig-zags" in profiles where lel defined
!             far too high (by use of lelten array, which assumes
!             a maximum of five such crossing points).
! feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.  previous
!             version ccm3_cloud_buoyancy.
!
!-------------------------------------------------------------------------------
#include <implicit.h>
! input arguments
!
   real     ::  q(ILOTS,levs_)        ! spec. humidity
   real     ::  t(ILOTS,levs_)        ! temperature
   real     ::  p(ILOTS,levs_)        ! pressure
   real     ::  z(ILOTS,levs_)        ! height
   real     ::  pf(ILOTS,levs_+1)     ! pressure at interfaces
   real     ::  pblt(ILOTS)          ! index of pbl depth
   real     ::  tpert(ILOTS)         ! perturbation temperature by pbl processes
   real     ::  qpert(ILOTS)         ! perturbation moisture by pbl processes
!
! output arguments
!
   real     ::  tp(ILOTS,levs_)       ! parcel temperature
   real     ::  qstp(ILOTS,levs_)     ! saturation mixing ratio of parcel
   real     ::  tl(ILOTS)            ! parcel temperature at lcl
   real     ::  cape(ILOTS)          ! convective aval. pot. energy.
   integer  ::  lcl(ILOTS)        ! 
   integer  ::  lel(ILOTS)        ! 
   integer  ::  lon(ILOTS)        ! level of onset of deep convection
   integer  ::  mx(ILOTS)         ! level of max moist static energy
!
!--------------------------Local Variables------------------------------
!
   real     ::  capeten(ILOTS,5)     ! provisional value of cape
   real     ::  tv(ILOTS,levs_)       ! 
   real     ::  tpv(ILOTS,levs_)      ! 
   real     ::  buoy(ILOTS,levs_)

   real     ::  a1(ILOTS) 
   real     ::  a2(ILOTS) 
   real     ::  estp(ILOTS) 
   real     ::  pl(ILOTS) 
   real     ::  plexp(ILOTS) 
   real     ::  hmax(ILOTS) 
   real     ::  hmn(ILOTS) 
   real     ::  y(ILOTS)

   logical  ::  plge600(ILOTS) 
   integer  ::  knt(ILOTS) 
   integer  ::  lelten(ILOTS,5)

   real     ::  cp
   real     ::  e
   real     ::  grav

   integer  ::  i
   integer  ::  k
   integer  ::  lat
   integer  ::  msg
   integer  ::  n
   integer  ::  nstep

   real     ::  rd
   real     ::  rl
!
!-------------------------------------------------------------------------------
!
   do n = 1,5
     do i = 1,ILOTS
       lelten(i,n) = levs_
       capeten(i,n) = 0.
     end do
   end do
!
   do i = 1,ILOTS
     lon(i) = levs_
     knt(i) = 0
     lel(i) = levs_
     mx(i) = lon(i)
     cape(i) = 0.
     hmax(i) = 0.
   end do
!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
   do k = levs_,msg + 1,-1
     do i = 1,ILOTS
       hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
       if (k.ge.nint(pblt(i)) .and. k.le.lon(i) .and.                          &
            hmn(i).gt.hmax(i)) then
         hmax(i) = hmn(i)
         mx(i) = k
       end if
     end do
   end do
!
   do i = 1,ILOTS
     lcl(i) = mx(i)
     e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
     tl(i) = 2840./ (3.5*log(t(i,mx(i)))-log(e)-4.805) + 55.
     if (tl(i).lt.t(i,mx(i))) then
       plexp(i) = (1./ (0.2854* (1.-0.28*q(i,mx(i)))))
       pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)

     else
       tl(i) = t(i,mx(i))
       pl(i) = p(i,mx(i))
     end if
   end do
!
! calculate lifting condensation level (lcl).
!
   do k = levs_,msg + 2,-1
     do i = 1,ILOTS
       if (k.le.mx(i) .and. (p(i,k).gt.pl(i).and.                              &
            p(i,k-1).le.pl(i))) then
         lcl(i) = k - 1
       end if
     end do
   end do
!
! if lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ILOTS
     plge600(i) = pl(i).ge.600.
   end do
!
! initialize parcel properties in sub-cloud layer below lcl.
!
   do k = levs_,msg + 1,-1
     do i = 1,ILOTS
       if (k.gt.lcl(i) .and. k.le.mx(i) .and. plge600(i)) then
         tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
         qstp(i,k) = q(i,mx(i))
         tp(i,k) = t(i,mx(i))* (p(i,k)/p(i,mx(i)))**                           &
                   (0.2854* (1.-0.28*q(i,mx(i))))
!
! buoyancy is increased by 0.5 k as in tiedtke
!
!              tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))/
!         1                     (1.+q(i,mx(i)))
         tpv(i,k) = (tp(i,k)+tpert(i))*                                        &
                    (1.+1.608*q(i,mx(i)))/ (1.+q(i,mx(i)))
         buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
       end if
     end do
   end do
!
! define parcel properties at lcl (i.e. level immediately above pl).
!
   do k = levs_,msg + 1,-1
     do i = 1,ILOTS
       if (k.eq.lcl(i) .and. plge600(i)) then
         tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
         qstp(i,k) = q(i,mx(i))
         tp(i,k) = tl(i)* (p(i,k)/pl(i))**                                     &
                   (0.2854* (1.-0.28*qstp(i,k)))
!              estp(i)  =exp(a-b/tp(i,k))
! use of different formulas for est has about 1 g/kg difference
! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
! above giving larger qs.
!
         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/                              &
                   ((tp(i,k)-tfreez)+c3))

         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
         a1(i) = cp/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)*rl*                    &
                 eps1/ (rd*tp(i,k)**2)
         a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))*                       &
                 (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/                            &
                 (rd**2*tp(i,k)**4)-qstp(i,k)*                                 &
                 (1.+qstp(i,k)/eps1)*2.*eps1*rl/                               &
                 (rd*tp(i,k)**3))
         a1(i) = 1./a1(i)
         a2(i) = -a2(i)*a1(i)**3
         y(i) = q(i,mx(i)) - qstp(i,k)
         tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/                              &
                 ((tp(i,k)-tfreez)+c3))

         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
!
! buoyancy is increased by 0.5 k in cape calculation.
! dec. 9, 1994
!              tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))
!
         tpv(i,k) = (tp(i,k)+tpert(i))* (1.+1.608*qstp(i,k))/                  &
                    (1.+q(i,mx(i)))
         buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
       end if
     end do
   end do
!
! main buoyancy calculation.
!
   do k = levs_ - 1,msg + 1,-1
     do i = 1,ILOTS
       if (k.lt.lcl(i) .and. plge600(i)) then
         tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
         qstp(i,k) = qstp(i,k+1)
         tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**                              &
                   (0.2854* (1.-0.28*qstp(i,k)))
!          estp(i) = exp(a-b/tp(i,k))
         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/                              &
                   ((tp(i,k)-tfreez)+c3))

         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
         a1(i) = cp/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)*rl*                    &
                 eps1/ (rd*tp(i,k)**2)
         a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))*                       &
                 (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/                            &
                 (rd**2*tp(i,k)**4)-qstp(i,k)*                                 &
                 (1.+qstp(i,k)/eps1)*2.*eps1*rl/                               &
                 (rd*tp(i,k)**3))
         a1(i) = 1./a1(i)
         a2(i) = -a2(i)*a1(i)**3
         y(i) = qstp(i,k+1) - qstp(i,k)
         tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/                              &
                   ((tp(i,k)-tfreez)+c3))

         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
!              tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
!              (1.+q(i,mx(i)))
         tpv(i,k) = (tp(i,k)+tpert(i))* (1.+1.608*qstp(i,k))/                  &
                    (1.+q(i,mx(i)))
         buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
       end if
     end do
   end do
!
   do k = msg + 2,levs_
     do i = 1,ILOTS
       if (k.lt.lcl(i) .and. plge600(i)) then
         if (buoy(i,k+1).gt.0. .and. buoy(i,k).le.0.) then
           knt(i) = min(5,knt(i) + 1)
           lelten(i,knt(i)) = k
         end if
       end if
     end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
     do k = msg + 1,levs_
       do i = 1,ILOTS
         if (plge600(i) .and.k.le.mx(i) .and.k.gt.lelten(i,n)) then
            capeten(i,n) = capeten(i,n) +                                      &
                          rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
         end if
       end do
     end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
     do i = 1,ILOTS
       if (capeten(i,n).gt.cape(i)) then
         cape(i) = capeten(i,n)
         lel(i) = lelten(i,n)
       end if
     end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ILOTS
     cape(i) = max(cape(i), 0.)
   end do
!
   return
   end subroutine ccm3_cloud_buoyancy
!
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_trigger(n,array,inc,target,index,nval)
!-------------------------------------------------------------------------------
   real     ::  array(*)
   integer  ::  index(*)
!-------------------------------------------------------------------------------
   ina=1
   nval=0
   if(inc .lt. 0) ina=(-inc)*(n-1)+1
   do i = 1,n
     if(array(ina) .gt. target) then
       nval=nval+1
       index(nval)=i
     endif
     ina=ina+inc
   enddo
!
   return
   end subroutine ccm3_cloud_trigger         
!
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_property (q      ,t      ,u      ,v      ,p      ,    &
                                   z      ,s      ,mu     ,eu     ,du     ,    &
                                   md     ,ed     ,sd     ,qd     ,ud     ,    &
                                   vd     ,mc     ,qu     ,su     ,zf     ,    &
                                   qst    ,hmn    ,hsat   ,alpha  ,shat   ,    &
                                   ql     ,totpcp ,totevp ,cmeg   ,jb     ,    &
                                   lel    ,jt     ,jlcl   ,mx     ,j0     ,    &
                                   jd     ,rl     ,il1g   ,il2g   ,rd     ,    &
                                   grav   ,cp     ,msg    ,nstep  ,lat    ,    &
                                   pflx   ,evp    ,cu     ,mu2    ,eu2    ,    &
                                   du2    ,md2    ,ed2     )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use guang
!-------------------------------------------------------------------------------
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996
!
!**** PLEASE NOTE ***!
!
! we are aware of a specific problem in this code 
! (identified by the string ---> PROBLEM ONE)
! during the calculation of the updraft cloud properties,
! rather than adding a perturbation to the updraft temperature of 
! half a degree, (there was an inadvertant addition of cp*0.5) degrees
! or about 500 degrees. (This problem was in the code prior to its 
! contribution to the NCAR effort)

! Fortunately, the erroneous values
! are overwritten later in the code. The problem is quite subtle.
! The erroneous values would persist between cloud base and the lifting 
! condensation level. The addition of the very high perturbation to the updraft
! temperature causes the saturation mixing ratio to be set to zero, 
! and later the lcl to be set to one level above cloud base.
! There are therefore no levels between cloud base and the lcl. Therefore
! all erroneous values are overwritten.

! The only manifestation we are aware of with respect to this problem
! is that the lifting condensation level is constrained to be one level above
! cloud base.

! We discovered the problem after too much had been invested in
! very long integrations (in terms of computer time)
! to allow for a modification and model retuning. It is our expectation that
! this problem will be fixed in the next release of the model.
!
! nov 20/92 - guang jun zhang,m.lazare. now has deeper (more
!             realistic) downdrafts.
! jul 14/92 - guang jun zhang,m.lazare. add shallow mixing
!             formulation.
! nov 21/91 - m.lazare. like previous cldprop except minimum "f"
!                       now 0.0004 instead of 0.001 (more
!                       realistic with more deep).
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!------------------------------------------------------------------------------
!
! Input arguments
!
   real     ::  q(ILOTS,levs_)        ! spec. humidity of env
   real     ::  t(ILOTS,levs_)        ! temp of env
   real     ::  p(ILOTS,levs_)        ! pressure of env
   real     ::  z(ILOTS,levs_)        ! height of env
   real     ::  s(ILOTS,levs_)        ! normalized dry static energy of env
   real     ::  zf(ILOTS,levs_+1)      ! height of interfaces
   real     ::  u(ILOTS,levs_)        ! zonal velocity of env
   real     ::  v(ILOTS,levs_)        ! merid. velocity of env

   integer  ::  jb(ILOTS)         ! updraft base level
   integer  ::  lel(ILOTS)        ! updraft launch level
   integer  ::  jt(ILOTS)         ! updraft plume top
   integer  ::  jlcl(ILOTS)       ! updraft lifting cond level
   integer  ::  mx(ILOTS)         ! updraft base level (same is jb)
   integer  ::  j0(ILOTS)         ! level where updraft begins detraining
   integer  ::  jd(ILOTS)         ! level of downdraft
!
! output
!
   real     ::  alpha(ILOTS,levs_)    !
   real     ::  cmfdqr(ILOTS,levs_)   ! rate of production of precip at that layer
   real     ::  du(ILOTS,levs_)       ! detrainement rate of updraft
   real     ::  ed(ILOTS,levs_)       ! entrainment rate of downdraft
   real     ::  eu(ILOTS,levs_)       ! entrainment rate of updraft
   real     ::  hmn(ILOTS,levs_)      ! moist stat energy of env
   real     ::  hsat(ILOTS,levs_)     ! sat moist stat energy of env
   real     ::  mc(ILOTS,levs_)       ! net mass flux
   real     ::  md(ILOTS,levs_)       ! downdraft mass flux
   real     ::  mu(ILOTS,levs_)       ! updraft mass flux
   real     ::  pflx(ILOTS,levs_)     ! precipitation flux thru layer
   real     ::  qd(ILOTS,levs_)       ! spec humidity of downdraft
   real     ::  ql(ILOTS,levs_)       ! liq water of updraft
   real     ::  qst(ILOTS,levs_)      ! saturation spec humidity of env.
   real     ::  qu(ILOTS,levs_)       ! spec hum of updraft
   real     ::  sd(ILOTS,levs_)       ! normalized dry stat energy of downdraft
   real     ::  shat(ILOTS,levs_)     ! interface values of dry stat energy
   real     ::  su(ILOTS,levs_)       ! normalized dry stat energy of updraft
   real     ::  ud(ILOTS,levs_)       ! downdraft u
   real     ::  vd(ILOTS,levs_)       ! downdraft v
!
!     these version of the mass fluxes conserve mass 
!
   real     ::  mu2(ILOTS,levs_)      ! updraft mass flux
   real     ::  eu2(ILOTS,levs_)      ! updraft entrainment
   real     ::  du2(ILOTS,levs_)      ! updraft detrainment
   real     ::  md2(ILOTS,levs_)      ! downdraft mass flux
   real     ::  ed2(ILOTS,levs_)      ! downdraft entrainment
   real     ::  rl                   ! latent heat of vap

   integer  ::  il1g              !CORE GROUP REMOVE
   integer  ::  il2g              !CORE GROUP REMOVE

   real     ::  rd                   ! gas constant for dry air
   real     ::  grav                 ! gravity
   real     ::  cp                   ! heat capacity of dry air

   integer  ::  msg               ! missing moisture vals (always 0)
   integer  ::  nstep             ! time step index
   integer  ::  lat               ! lat index
!
! Local workspace
!
   real     ::  gamma(ILOTS,levs_)  
   real     ::  dz(ILOTS,levs_)  
   real     ::  iprm(ILOTS,levs_)  
   real     ::  hu(ILOTS,levs_)  
   real     ::  hd(ILOTS,levs_)  
   real     ::  eps(ILOTS,levs_)  
   real     ::  f(ILOTS,levs_)  
   real     ::  k1(ILOTS,levs_)  
   real     ::  i2(ILOTS,levs_)  
   real     ::  ihat(ILOTS,levs_)  
   real     ::  i3(ILOTS,levs_)  
   real     ::  idag(ILOTS,levs_)  
   real     ::  i4(ILOTS,levs_)  
   real     ::  qsthat(ILOTS,levs_)  
   real     ::  hsthat(ILOTS,levs_)  
   real     ::  gamhat(ILOTS,levs_)  
   real     ::  cu(ILOTS,levs_)  
   real     ::  evp(ILOTS,levs_)  
   real     ::  cmeg(ILOTS,levs_)  
   real     ::  qds(ILOTS,levs_) 
   real     ::  hmin(ILOTS)  
   real     ::  expdif(ILOTS)  
   real     ::  expnum(ILOTS)  
   real     ::  ftemp(ILOTS)  
   real     ::  eps0(ILOTS)  
   real     ::  rmue(ILOTS)  
   real     ::  zuef(ILOTS)  
   real     ::  zdef(ILOTS)  
   real     ::  epsm(ILOTS)  
   real     ::  ratmjb(ILOTS)  
   real     ::  est(ILOTS)  
   real     ::  totpcp(ILOTS)  
   real     ::  totevp(ILOTS)  
   real     ::  alfa(ILOTS) 
   real     ::  beta
   real     ::  c0
   real     ::  ql1
   real     ::  weight
   real     ::  tu
   real     ::  estu
   real     ::  qstu

   real     ::  small
   real     ::  mdt  
   real     ::  cu2

   integer  ::  khighest
   integer  ::  klowest  
   integer  ::  kount 
   integer  ::  i,k

   logical  ::  doit(ILOTS)
   logical  ::  done(ILOTS)
!
!------------------------------------------------------------------------------
!
   do i = 1,il2g
     ftemp(i) = 0.
     expnum(i) = 0.
     expdif(i) = 0.
   end do
!
!   Change from msg+1 to 1 to prevent blowup
!
   do k = 1,levs_
     do i = 1,il2g
       dz(i,k) = zf(i,k) - zf(i,k+1)
     end do
   end do

!
! initialize many output and work variables to zero
!
   do k = msg + 1,levs_
     do i = 1,il2g
       k1(i,k) = 0.
       i2(i,k) = 0.
       i3(i,k) = 0.
       i4(i,k) = 0.
       mu(i,k) = 0.
       f(i,k) = 0.
       eps(i,k) = 0.
       eu(i,k) = 0.
       du(i,k) = 0.
       ql(i,k) = 0.
       cu(i,k) = 0.
       evp(i,k) = 0.
       cmeg(i,k) = 0.
       qds(i,k) = q(i,k)
       md(i,k) = 0.
       ed(i,k) = 0.
       sd(i,k) = s(i,k)
       qd(i,k) = q(i,k)
       ud(i,k) = u(i,k)
       vd(i,k) = v(i,k)
       mc(i,k) = 0.
       qu(i,k) = q(i,k)
       su(i,k) = s(i,k)
!        est(i)=exp(a-b/t(i,k))
       est(i) = c1*exp((c2* (t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
       qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
       gamma(i,k) = qst(i,k)*(1. + qst(i,k)/eps1)*eps1*rl/                     &
                    (rd*t(i,k)**2)*rl/cp
       hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
       hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
       hu(i,k) = hmn(i,k)
       hd(i,k) = hmn(i,k)
       mu2(i,k) = 0.
       eu2(i,k) = 0.
       du2(i,k) = 0.
       md2(i,k) = 0.
       ed2(i,k) = 0.
       pflx(i,k) = 0.
       cmfdqr(i,k) = 0.
     end do
   end do
!
!   Set to zero things which make this routine blow up
!
   do k = 1,msg
     do i = 1,il2g
       cmfdqr(i,k) = 0.
       mu2(i,k) = 0.
       eu2(i,k) = 0.
       du2(i,k) = 0.
       md2(i,k) = 0.
       ed2(i,k) = 0.
     end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do i = 1,il2g
     hsthat(i,msg+1) = hsat(i,msg+1)
     qsthat(i,msg+1) = qst(i,msg+1)
     gamhat(i,msg+1) = gamma(i,msg+1)
     totpcp(i) = 0.
     totevp(i) = 0.
   end do
   do k = msg + 2,levs_
     do i = 1,il2g
       if (abs(qst(i,k-1)-qst(i,k)).gt.1.E-6) then
         qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*                    &
                       qst(i,k)/ (qst(i,k-1)-qst(i,k))
       else
         qsthat(i,k) = qst(i,k)
       end if
       hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
       if (abs(gamma(i,k-1)-gamma(i,k)).gt.1.E-6) then
         gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*                           &
                       gamma(i,k-1)*gamma(i,k)/                                &
                       (gamma(i,k-1)-gamma(i,k))
       else
         gamhat(i,k) = gamma(i,k)
       end if
     end do
   end do
!
! initialize cloud top to highest plume top.
!
   do i = 1,il2g
     jt(i) = max(lel(i),4)
     jd(i) = levs_
     jlcl(i) = lel(i)
     hmin(i) = 1.E6
   end do
!
! find the level of minimum hsat, where detrainment starts
!
   do k = msg + 1,levs_
     do i = 1,il2g
       if (hsat(i,k).le.hmin(i) .and. k.ge.jt(i).and.k.le.jb(i)) then
         hmin(i) = hsat(i,k)
         j0(i) = k
       end if
     end do
   end do
   do i = 1,il2g
     j0(i) = min(j0(i),jb(i)-2)
     j0(i) = max(j0(i),jt(i)+2)
!
! Fix from Guang Zhang to address out of bounds array reference
!
     j0(i) = min(j0(i),levs_)
   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,levs_
     do i = 1,il2g
       if (k.ge.jt(i) .and. k.le.jb(i)) then
         hu(i,k) = hmn(i,mx(i)) + cp*0.5
         su(i,k) = s(i,mx(i)) + 0.5
!*** PROBLEM ONE **!
         su(i,k) = s(i,mx(i)) + cp*0.5
       end if
     end do
   end do
!
! ********************************************************!
! compute taylor series for approximate eps(z) below
! ********************************************************!
!
   do k = levs_ - 1,msg + 1,-1
     do i = 1,il2g
       if (k.lt.jb(i) .and. k.ge.jt(i)) then
         k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
         ihat(i,k) = 0.5* (k1(i,k+1)+k1(i,k))
         i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
         idag(i,k) = 0.5* (i2(i,k+1)+i2(i,k))
         i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
         iprm(i,k) = 0.5* (i3(i,k+1)+i3(i,k))
         i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
       end if
     end do
   end do
!
! re-initialize hmin array for ensuing calculation.
!
   do i = 1,il2g
     hmin(i) = 1.E6
   end do
   do k = msg + 1,levs_
     do i = 1,il2g
       if (k.ge.j0(i).and.k.le.jb(i) .and. hmn(i,k).le.hmin(i)) then
         hmin(i) = hmn(i,k)
         expdif(i) = hmn(i,mx(i)) - hmin(i)
       end if
     end do
   end do
!
! ********************************************************!
! compute approximate eps(z) using above taylor series
! ********************************************************!
!
   do k = msg + 2,levs_
     do i = 1,il2g
       expnum(i) = 0.
       ftemp(i) = 0.
       if (k.lt.jt(i) .or. k.ge.jb(i)) then
         k1(i,k) = 0.
         expnum(i) = 0.
       else
         expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) +            &
                     hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
       end if
       if ((expdif(i).gt.100..and.expnum(i).gt.0.) .and.                       &
           k1(i,k).gt.expnum(i)*dz(i,k)) then
         ftemp(i) = expnum(i)/k1(i,k)
         f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 +                     &
                   (2.*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2*                 &
                   ftemp(i)**3 + (-5.*k1(i,k)*i2(i,k)*i3(i,k)+                 &
                   5.*i2(i,k)**3+k1(i,k)**2*i4(i,k))/                          &
                   k1(i,k)**3*ftemp(i)**4
         f(i,k) = max(f(i,k),0.)
         f(i,k) = min(f(i,k),0.0002)
       end if
     end do
   end do
   do i = 1,il2g
     if (j0(i).lt.jb(i)) then
       if (f(i,j0(i)).lt.1.E-6 .and. f(i,j0(i)+1).gt.f(i,j0(i)))               &
         j0(i) = j0(i) + 1
     end if
   end do
   do k = msg + 2,levs_
     do i = 1,il2g
       if (k.ge.jt(i) .and. k.le.j0(i)) then
         f(i,k) = max(f(i,k),f(i,k-1))
       end if
     end do
   end do
   do i = 1,il2g
     eps0(i) = f(i,j0(i))
     eps(i,jb(i)) = eps0(i)
   end do
   do k = levs_,msg+1,-1
     do i = 1,il2g
       if (k.ge.j0(i)) then
         if (k.le.jb(i)) eps(i,k) = f(i,j0(i))
       else
         if (k.ge.jt(i)) eps(i,k) = f(i,k)
       end if
     end do
   end do
!
! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
!
   do i = 1,il2g
     if (eps0(i).gt.0.) then
       mu(i,jb(i)) = 1.
       eu(i,jb(i)) = eps0(i)/2.
       mu2(i,jb(i)) = 1.
       eu2(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
     end if
   end do
   do k = levs_,msg + 1,-1
     do i = 1,il2g
       if (eps0(i).gt.0. .and. k.ge.jt(i) .and. k.lt.jb(i)) then
         zuef(i) = zf(i,k) - zf(i,jb(i))
         rmue(i) = (1./eps0(i))* (exp(eps(i,k+1)*zuef(i))-1.)/zuef(i)
         mu(i,k) = (1./eps0(i))* (exp(eps(i,k)*zuef(i))-1.)/zuef(i)
         eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
         du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
         mu2(i,k) = mu(i,k)
         eu2(i,k) = eu(i,k)
         du2(i,k) = du(i,k)
       end if
     end do
   end do
!
   khighest = levs_+1
   klowest = 1
   do i = 1,il2g
     khighest = min(khighest,lel(i))
     klowest = max(klowest,jb(i))
   end do
   do k = klowest-1,khighest,-1
!dir$ ivdep
     do i = 1,il2g
       if (k.le.jb(i)-1 .and. k.ge.lel(i) .and. eps0(i).gt.0.) then
         if (mu(i,k).lt.0.01) then
           hu(i,k) = hu(i,jb(i))
           mu(i,k) = 0.
           mu2(i,k) = mu(i,k)
           eu2(i,k) = 0.
           du2(i,k) = mu2(i,k+1)/dz(i,k)
         else
           hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) +                             &
                    dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)-                        &
                    du(i,k)*hsat(i,k))
         end if
       end if
     end do
   end do
!
! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
!
   do i = 1,il2g
     doit(i) = .true.
   end do
   do k = klowest-2,khighest-1,-1
     do i = 1,il2g
       if (doit(i) .and. k.le.jb(i)-2 .and. k.ge.lel(i)-1) then
         if (hu(i,k  ).le.hsthat(i,k) .and.                                    &
             hu(i,k+1).gt.hsthat(i,k+1) .and. mu(i,k).ge.0.02) then
           if (hu(i,k)-hsthat(i,k).lt.-2000.) then
             jt(i) = k + 1
             doit(i) = .false.
           else
             jt(i) = k
             if (eps0(i).le.0.) doit(i) = .false.
           end if
         else if (hu(i,k).gt.hu(i,jb(i)) .or. mu(i,k).lt.0.01) then
           jt(i) = k + 1
           doit(i) = .false.
         end if
       end if
     end do
   end do
   do k = levs_,msg + 1,-1
!dir$ ivdep
     do i = 1,il2g
       if (k.ge.lel(i) .and. k.le.jt(i) .and. eps0(i).gt.0.) then
         mu(i,k) = 0.
         eu(i,k) = 0.
         du(i,k) = 0.
         mu2(i,k) = 0.
         eu2(i,k) = 0.
         du2(i,k) = 0.
         hu(i,k) = hu(i,jb(i))
       end if
       if (k.eq.jt(i) .and. eps0(i).gt.0.) then
         du(i,k) = mu(i,k+1)/dz(i,k)
         du2(i,k) = mu2(i,k+1)/dz(i,k)
         eu2(i,k) = 0.
         mu2(i,k) = 0.
       end if
     end do
   end do
!
! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
   do i = 1,il2g
!
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
!
     alfa(i) = 0.1
     jt(i) = min(jt(i),jb(i)-1)
     jd(i) = max(j0(i),jt(i)+1)
     jd(i) = min(jd(i),jb(i))
     hd(i,jd(i)) = hmn(i,jd(i)-1)
     ud(i,jd(i)) = u(i,jd(i)-1)
     vd(i,jd(i)) = v(i,jd(i)-1)
     if (jd(i).lt.jb(i) .and. eps0(i).gt.0.) then
       epsm(i) = eps0(i)
!          alfa(i)=2.*epsm(i)*( zf(i,jd(i))-zf(i,jb(i)) )/
!     1         (  exp(2.*epsm(i)*( zf(i,jd(i))-
!               zf(i,jb(i)) ))-1.  )
       md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
       md2(i,jd(i)) = md(i,jd(i))
     end if
   end do
   do k = msg + 1,levs_
     do i = 1,il2g
       if ((k.gt.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
         zdef(i) = zf(i,jd(i)) - zf(i,k)
         md(i,k) = -alfa(i)/ (2.*eps0(i))*                                     &
                  (exp(2.*epsm(i)*zdef(i))-1.)/zdef(i)
         md2(i,k) = md(i,k)
       end if
     end do
   end do
#ifdef DBG
   k = levs_
   do i = 1,il2g
     res = mu2(i,k)/dz(i,k) - eu2(i,k) + du2(i,k)
     if (res.gt.1.e-10*max(eu2(i,k),du2(i,k))) then
       write (6,*) 'inconsistent mass fluxes ', i,k,lat
       write (6,*) 'mt, mb, eu*dz, du*dz ', mu2(i,k),                          &
            0., eu2(i,k), du2(i,k)
       stop
     endif
   end do
   do k = 1,levs_-1
     do i = 1,il2g
       res = (mu2(i,k)-mu2(i,k+1))/dz(i,k) - eu2(i,k) + du2(i,k)
       if (res.gt.1.e-10*max(eu2(i,k),du2(i,k))) then
         write (6,*) 'inconsistent mass fluxes ', i,k,lat
         write (6,*) 'mt, mb, eu*dz, du*dz ', mu2(i,k),                        &
              mu2(i,k+1), eu2(i,k), du2(i,k)
         stop
       endif
     end do
   end do
#endif
   do k = msg + 1,levs_
!dir$ ivdep
     do i = 1,il2g
       if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.               &
            jd(i).lt.jb(i)) then
         ratmjb(i) = min(abs(mu2(i,jb(i))/md2(i,jb(i))),1.)
         md2(i,k) = md2(i,k)*ratmjb(i)
         ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1.)
         md(i,k) = md(i,k)*ratmjb(i)
       end if
     end do
   end do
   do k = msg + 1,levs_
      do i = 1,il2g
         if ((k.gt.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            hd(i,k) = md(i,k-1)/md(i,k)*hd(i,k-1) -                            &
                       dz(i,k-1)/md(i,k)*ed(i,k-1)*hmn(i,k-1)
            ud(i,k) = md(i,k-1)/md(i,k)*ud(i,k-1) -                            &
                       dz(i,k-1)/md(i,k)*ed(i,k-1)*u(i,k-1)
            vd(i,k) = md(i,k-1)/md(i,k)*vd(i,k-1) -                            &
                       dz(i,k-1)/md(i,k)*ed(i,k-1)*v(i,k-1)
         end if
      end do
   end do
   small = 1.e-20
   do k = msg + 1,levs_
     do i = 1,il2g
       if ((k.ge.jt(i).and.k.le.levs_) .and. eps0(i).gt.0.) then
!         if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
!         if ((k.gt.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
         ed2(i,k-1) = (md2(i,k-1)-md2(i,k))/dz(i,k-1)
!            mdt = min(md2(i,k),-small)
!            hd(i,k) = (md(i,k-1)*hd(i,k-1) -
!     $                 dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
       end if
     end do
   end do
!
! calculate updraft and downdraft properties.
!
   do k = msg + 2,levs_
     do i = 1,il2g
       if ((k.ge.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.               &
           jd(i).lt.jb(i)) then
!         sd(i,k) = shat(i,k)
!    1             +              (hd(i,k)-hsthat(i,k))/
!    2               (cp    *(1.+gamhat(i,k)))
         qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/           &
                    (rl*(1. + gamhat(i,k)))
       end if
     end do
   end do
!
   do i = 1,il2g
      done(i) = .false.
   end do
   kount = 0
   do k = levs_,msg + 2,-1
     do i = 1,il2g
       if ((.not.done(i) .and. k.gt.jt(i) .and. k.lt.jb(i)) .and.              &
            eps0(i).gt.0.) then
         su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) +                               &
                   dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
         qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) +                               &
                   dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)-                           &
                   du(i,k)*qst(i,k))
         tu = su(i,k) - grav/cp*zf(i,k)
         estu = c1*exp((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
         qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2.-estu)
         if (qu(i,k).ge.qstu) then
           jlcl(i) = k
           kount = kount + 1
           done(i) = .true.
         end if
       end if
     end do
     if (kount.ge.il2g) goto 690
   end do
 690  continue
   do k = msg + 2,levs_
     do i = 1,il2g
       if (k.eq.jb(i) .and. eps0(i).gt.0.) then
         qu(i,k) = q(i,mx(i))
         su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
       end if
       if ((k.gt.jt(i).and.k.le.jlcl(i)) .and. eps0(i).gt.0.) then
         su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/                          &
                  (cp* (1.+gamhat(i,k)))
         qu(i,k) = qsthat(i,k) + gamhat(i,k)*                                  &
                  (hu(i,k)-hsthat(i,k))/                                       &
                  (rl* (1.+gamhat(i,k)))
       end if
     end do
   end do
!
   do k = levs_,msg + 2,-1
     do i = 1,il2g
       if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0.) then
         cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/                     &
                  dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/                          &
                  (rl/cp)
         if (k.eq.jt(i)) cu(i,k) = 0.
!               cu(i,k) = max(0.,cu(i,k))
!               cu2     = max(0.,
!               cu2     = max(-1.e99,
!     $                   +(eu(i,k)*q(i,k) - du(i,k)*qst(i,k))
!     $                   -(mu(i,k)*qu(i,k)-mu(i,k+1)*qu(i,k+1))/dz(i,k)
!     $                   )
!                  
!               if (abs(cu(i,k)-cu2)/(abs(cu(i,k))+abs(cu2)+1.e-50)
!     $              .gt.0.0000001) then
!                  write (6,*) ' inconsistent condensation rates ', 
!     $                 i, k, lat,
!     $                 cu(i,k), cu2, jt(i), jb(i), jlcl(i), lel(i)
!     $                 ,mu(i,k)
!               endif
       end if
     end do
   end do
!
   beta = 0.
   c0 = 2.E-3
   do k = levs_,msg + 2,-1
     do i = 1,il2g
       cmfdqr(i,k) = 0.
! this modification is for test3 run, modified on 6/20/1995
!          if(t(i,jt(i) ).gt.tfreez)    c0=0.
!          if(t(i,jt(i) ).le.tfreez   )    c0=2.e-3
       if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0. .and.               &
            mu(i,k).ge.0.0) then
         if (mu(i,k).gt.0.) then
           ql1 = 1./mu(i,k)* (mu(i,k+1)*ql(i,k+1)-                             &
                dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
           ql(i,k) = ql1/ (1.+dz(i,k)*c0)
         else
           ql(i,k) = 0.
         end if
         totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*                     &
                     (beta*ql(i,k) + (1. - beta)*ql(i,k+1)))
         cmfdqr(i,k) = c0*mu(i,k)*ql(i,k)
       end if
     end do
   end do
!
   do i = 1,il2g
     qd(i,jd(i)) = qds(i,jd(i))
     sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
#ifdef NEWVER
!
!v     mod to make my downdraft calculaton match guangs assumption
!      but still conserve energy. For consistency we should probably 
!      also set sd(i,jd(i) = s(i,jd(i)-1)
!      it would also make more sense to me to use qhat and shat
!      but I am trying to minimize changes to guangs algorithm
!
!         qd(i,jd(i)) = q(i,jd(i)-1)
!
#endif
#undef NEWVER
   end do
!
#undef GUANGSWAY3
#define GUANGSWAY3
#ifdef GUANGSWAY3
   do k = msg + 2,levs_
     do i = 1,il2g
       if ((k.ge.jd(i).and.k.lt.jb(i)) .and. eps0(i).gt.0. .and.               &
            jd(i).lt.jb(i)) then
         qd(i,k+1) = qds(i,k+1)
         sd(i,k+1) = (hd(i,k+1)-rl*qd(i,k+1))/cp
         evp(i,k) = -ed(i,k)*q(i,k) +                                          &
                   (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
         if (k.eq.jd(i)) then
           evp(i,k) = -ed(i,k)*q(i,k) +                                        &
                     (md(i,k)*q(i,k-1)-md(i,k+1)*qd(i,k+1))/dz(i,k)
         end if
         evp(i,k) = max(evp(i,k),0.)
         totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   end do
!
   do i = 1,il2g
     totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -                        &
                 md(i,jb(i))*qd(i,jb(i))
   end do
!
! no evaporation below cloud base is considered
!
#else
   do k = msg + 2,levs_
     do i = 1,il2g
       if (k.ge.jd(i).and.k.lt.jb(i) .and. eps0(i).gt.0.) then
         qd(i,k+1) = qds(i,k+1)
         evp(i,k) = -ed(i,k)*q(i,k) +                                          &
                   (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
         evp(i,k) = max(evp(i,k),0.)
         mdt = min(md(i,k+1),-small)
         sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) +                &
                       md(i,k)*sd(i,k))/mdt
         totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
       end if
     end do
   end do
   do i = 1,il2g
!    totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
     totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) -                         &
                 md(i,jb(i))*qd(i,jb(i))
   end do
   if (.true.) then
     do i = 1,il2g
       k = jb(i)
       if (eps0(i).gt.0.) then
         evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
         evp(i,k) = max(evp(i,k),0.)
         totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
       end if
     end do
   endif
#endif
   do i = 1,il2g
     totpcp(i) = max(totpcp(i),0.)
     totevp(i) = max(totevp(i),0.)
   end do
!
   weight = 1.0
   do k = msg + 2,levs_
     do i = 1,il2g
       if (totevp(i).gt.0. .and. totpcp(i).gt.0. .and.                         &
            k.ge.jd(i) .and. k.le.jb(i)) then
         md(i,k) = md(i,k)*min(1.,weight*totpcp(i)/                            &
                  (totevp(i)+weight*totpcp(i)))
         ed(i,k) = ed(i,k)*min(1.,weight*totpcp(i)/                            &
                  (totevp(i)+weight*totpcp(i)))
         evp(i,k) = evp(i,k)*min(1.,                                           &
                   weight*totpcp(i)/ (totevp(i)+                               &
                   weight*totpcp(i)))
       else
         md(i,k) = 0.
         ed(i,k) = 0.
         evp(i,k) = 0.
       end if
!
! cmeg is the cloud water condensed - rain water evaporated
! cmfdqr  is the cloud water converted to rain - (rain evaporated)
!
       cmeg(i,k) = cu(i,k) - evp(i,k)
       cmfdqr(i,k) = cmfdqr(i,k)-evp(i,k)
     end do
   end do
!
   do k = msg + 2,levs_
     do i = 1,il2g
       if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
         md2(i,k) = md2(i,k)*min(1.,weight*totpcp(i)/                          &
                   (totevp(i)+weight*totpcp(i)))
         ed2(i,k) = ed2(i,k)*min(1.,weight*totpcp(i)/                          &
                   (totevp(i)+weight*totpcp(i)))
       else
         md2(i,k) = 0.
         ed2(i,k) = 0.
       end if
     end do
   end do
   do k = 2,levs_
     do i = 1,il2g
       pflx(i,k) = pflx(i,k-1) + cmfdqr(i,k)*dz(i,k)
     end do
   end do
   do i = 1,il2g
     if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
       totevp(i) = totevp(i)*min(1.,                                           &
                  weight*totpcp(i)/(totevp(i) + weight*totpcp(i)))
     else
       totevp(i) = 0.
     end if
   end do
!
   do k = msg + 1,levs_
     do i = 1,il2g
       if (k.ge.jt(i) .and. k.le.jb(i)) then
         mc(i,k) = mu(i,k) + md(i,k)
       end if
     end do
   end do
#undef GUANGSWAY3
!
   return
   end subroutine ccm3_cloud_property
!
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_closure  (q       ,t       ,p       ,                 &
                          z       ,s       ,tp      ,qs      ,qu      ,        &
                          su      ,mc      ,du      ,mu      ,md      ,        &
                          qd      ,sd      ,alpha   ,qhat    ,shat    ,        &
                          dp      ,qstp    ,zf      ,ql      ,dsubcld ,        &
                          mb      ,cape    ,tl      ,lcl     ,lel     ,        &
                          jt      ,mx      ,il1g    ,il2g    ,rd      ,        &
                          grav    ,cp      ,rl      ,msg     ,capelmt ,        &
                          nstep   ,lat     )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use guang
!-------------------------------------------------------------------------------
!
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996
!
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!
!-----------------------------Arguments---------------------------------
!
   real                 ::  q(ILOTS,levs_)        ! spec humidity
   real                 ::  t(ILOTS,levs_)        ! temperature
   real                 ::  p(ILOTS,levs_)        ! pressure (mb)
   real                 ::  z(ILOTS,levs_)        ! height (m)
   real                 ::  s(ILOTS,levs_)        ! normalized dry static energy 
   real                 ::  tp(ILOTS,levs_)       ! parcel temp
   real                 ::  qs(ILOTS,levs_)       ! sat spec humidity
   real                 ::  qu(ILOTS,levs_)       ! updraft spec. humidity
   real                 ::  su(ILOTS,levs_)       ! normalized dry stat energy of updraft
   real                 ::  mc(ILOTS,levs_)       ! net convective mass flux 
   real                 ::  du(ILOTS,levs_)       ! detrainment from updraft
   real                 ::  mu(ILOTS,levs_)       ! mass flux of updraft
   real                 ::  md(ILOTS,levs_)       ! mass flux of downdraft
   real                 ::  qd(ILOTS,levs_)       ! spec. humidity of downdraft
   real                 ::  sd(ILOTS,levs_)       ! dry static energy of downdraft
   real                 ::  alpha(ILOTS,levs_)
   real                 ::  qhat(ILOTS,levs_)     ! environment spec humidity at interfaces
   real                 ::  shat(ILOTS,levs_)     ! env. normalized dry static energy at intrfcs
   real                 ::  dp(ILOTS,levs_)       ! pressure thickness of layers
   real                 ::  qstp(ILOTS,levs_)     ! spec humidity of parcel
   real                 ::  zf(ILOTS,levs_+1)     ! height of interface levels
   real                 ::  ql(ILOTS,levs_)       ! liquid water mixing ratio

   real                 ::  mb(ILOTS)            ! cloud base mass flux
   real                 ::  cape(ILOTS)          ! available pot. energy of column
   real                 ::  tl(ILOTS)
   real                 ::  dsubcld(ILOTS)       ! thickness of subcloud layer

   integer              ::  lcl(ILOTS)        ! index of lcl
   integer              ::  lel(ILOTS)        ! index of launch leve
   integer              ::  jt(ILOTS)         ! top of updraft
   integer              ::  mx(ILOTS)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
   real                 ::  dtpdt(ILOTS,levs_)
   real                 ::  dqsdtp(ILOTS,levs_)
   real                 ::  dtmdt(ILOTS,levs_)
   real                 ::  dqmdt(ILOTS,levs_)
   real                 ::  dboydt(ILOTS,levs_)
   real                 ::  thetavp(ILOTS,levs_)
   real                 ::  thetavm(ILOTS,levs_)

   real                 ::  dtbdt(ILOTS),dqbdt(ILOTS),dtldt(ILOTS)
   real                 ::  beta
   real                 ::  capelmt
   real                 ::  cp
   real                 ::  dadt
   real                 ::  debdt
   real                 ::  dltaa
   real                 ::  eb
   real                 ::  grav

   integer              ::  i
   integer              ::  il1g
   integer              ::  il2g
   integer              ::  k
   integer              ::  lat
   integer              ::  msg
   integer              ::  nstep

   real                 ::  rd
   real                 ::  rl
   real,save            ::  tau
!
! tau=4800. were used in canadian climate center. however, when it
! is used here in echam3 t42, convection is too weak, thus 
! adjusted to 2400. i.e the e-folding time is 1 hour now.
!
   data tau/7200./
!-------------------------------------------------------------------------------
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
   do i = il1g,il2g
     mb(i) = 0.
     eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
     dtbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*                                 &
                  (shat(i,mx(i))-su(i,mx(i)))+                                 &
                  md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
     dqbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*                                 &
                  (qhat(i,mx(i))-qu(i,mx(i)))+                                 &
                  md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
     debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
     dtldt(i) = -2840.* (3.5/t(i,mx(i))*dtbdt(i)-debdt/eb)/                    &
                  (3.5*log(t(i,mx(i)))-log(eb)-4.805)**2
   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,levs_
     do i = il1g,il2g
       dtmdt(i,k) = 0.
       dqmdt(i,k) = 0.
     end do
   end do
!
   do k = msg + 1,levs_ - 1
     do i = il1g,il2g
       if (k.eq.jt(i)) then
         dtmdt(i,k) = (1./dp(i,k))*                                            &
                          (mu(i,k+1)* (su(i,k+1)-shat(i,k+1)-                  &
                          rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-              &
                          shat(i,k+1)))
         dqmdt(i,k) = (1./dp(i,k))*(mu(i,k+1)* (qu(i,k+1)-                     &
                          qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*                    &
                          (qd(i,k+1)-qhat(i,k+1)))
        end if
     end do
   end do
!
   beta = 0.
   do k = msg + 1,levs_ - 1
     do i = il1g,il2g
       if (k.gt.jt(i) .and. k.lt.mx(i)) then
         dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+                            &
                          mc(i,k+1)* (s(i,k)-shat(i,k+1)))/                    &
                          dp(i,k) - rl/cp*du(i,k)*                             &
                          (beta*ql(i,k)+ (1-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
!     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     2                +du(i,k)*(qs(i,k)-q(i,k))
!     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

         dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+                      &
                          cp/rl* (su(i,k+1)-s(i,k)))-                          &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*                   &
                          (su(i,k)-s(i,k)))+md(i,k+1)*                         &
                          (qd(i,k+1)-qhat(i,k+1)+cp/rl*                        &
                          (sd(i,k+1)-s(i,k)))-md(i,k)*                         &
                          (qd(i,k)-qhat(i,k)+cp/rl*                            &
                          (sd(i,k)-s(i,k))))/dp(i,k) +                         &
                          du(i,k)* (beta*ql(i,k)+                              &
                          (1-beta)*ql(i,k+1))
       end if
     end do
   end do
!
   do k = msg + 1,levs_
     do i = il1g,il2g
       if (k.ge.lel(i) .and. k.le.lcl(i)) then
         thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)*                     &
                            (1.+1.608*qstp(i,k)-q(i,mx(i)))
         thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)*                      &
                            (1.+0.608*q(i,k))
         dqsdtp(i,k) = qstp(i,k)* (1.+qstp(i,k)/eps1)*eps1*rl/                 &
                            (rd*tp(i,k)**2)
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
         dtpdt(i,k) = tp(i,k)/ (1.+                                            &
                          rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))*             &
                           (dtbdt(i)/t(i,mx(i))+                               &
                          rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/                   &
                          tl(i)**2*dtldt(i)))
!
! dboydt is the integrand of cape change.
!
         dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1./                                &
                           (1.+1.608*qstp(i,k)-q(i,mx(i)))*                    &
                           (1.608 * dqsdtp(i,k) * dtpdt(i,k) -                 &
                           dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608/              &
                           (1.+0.608*q(i,k))*dqmdt(i,k)))*grav*                &
                           thetavp(i,k)/thetavm(i,k)
       end if
     end do
   end do
!
   do k = msg + 1,levs_
     do i = il1g,il2g
       if (k.gt.lcl(i) .and. k.lt.mx(i)) then
         thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)*                     &
                            (1.+0.608*q(i,mx(i)))
         thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)*                      &
                            (1.+0.608*q(i,k))
!
! dboydt is the integrand of cape change.
!
         dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+                                   &
                           0.608/ (1.+0.608*q(i,mx(i)))*dqbdt(i)-              &
                           dtmdt(i,k)/t(i,k)-                                  &
                           0.608/ (1.+0.608*q(i,k))*dqmdt(i,k))*               &
                           grav*thetavp(i,k)/thetavm(i,k)
       end if
     end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
   do i = il1g,il2g
     dadt = 0.
     do k = lel(i),mx(i) - 1
       dadt = dadt + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
     end do
!
     dltaa = -1.* (cape(i)-capelmt)
     if (dadt.ne.0.) mb(i) = max(dltaa/tau/dadt,0.)
   end do
!
   return
   end subroutine ccm3_cloud_closure
!
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_trans(q     ,mu     ,md      ,du     ,eu     ,        &
                               ed    ,dp     ,dsubcld ,jt     ,mx     ,        &
                               ideep ,il1g   ,il2g    ,nstep  ,lat    ,        &
                               delt    )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
!-------------------------------------------------------------------------------
!
! Convective transport of water species  : only water vapor here --- hong
!
! Note that we are still assuming that the water are in a moist mixing ratio
! this will change soon
!
!-------------------------Code History----------------------------------
!
! Original version:  P. Rasch, Jan 1996 
! Standardized:      L. Buja,  Feb 1996
! Reviewed:          P. Rasch, Feb 1996      
! 
!-------------------------------------------------------------------------------
#include <implicit.h>
!--------------------------Commons--------------------------------------
! 
! Input
!
   real                 ::  mu(ILOTS,levs_)       ! Mass flux up
   real                 ::  md(ILOTS,levs_)       ! Mass flux down
   real                 ::  du(ILOTS,levs_)       ! Mass detraining from updraft
   real                 ::  eu(ILOTS,levs_)       ! Mass entraining from updraft
   real                 ::  ed(ILOTS,levs_)       ! Mass entraining from downdraft
   real                 ::  dp(ILOTS,levs_)       ! Delta pressure between interfaces
   real                 ::  dsubcld(ILOTS)       ! Delta pressure from cloud base to sfc

   integer              ::  jt(ILOTS)         ! Index of cloud top for each column
   integer              ::  mx(ILOTS)         ! Index of cloud top for each column
   integer              ::  ideep(ILOTS)      ! Gathering array
   integer              ::  il1g              ! Gathered min lon indices over which to operate
   integer              ::  il2g              ! Gathered max lon indices over which to operate
   integer              ::  lat               ! Latitude
   integer              ::  nstep             ! Time step index

   real                 ::  delt                 ! Time step

! input/output

   real                 ::  q(ILOTS,levs_,1)  ! water array including moisture

!--------------------------Local Variables------------------------------

   integer              ::  i                 ! Work index
   integer              ::  k                 ! Work index
   integer              ::  kbm               ! Highest altitude index of cloud base
   integer              ::  kk                ! Work index
   integer              ::  kkp1              ! Work index
   integer              ::  km1               ! Work index
   integer              ::  kp1               ! Work index
   integer              ::  ktm               ! Highest altitude index of cloud top
   integer              ::  m                 ! Work index

   real                 ::  cabv                 ! Mix ratio of constituent above
   real                 ::  cbel                 ! Mix ratio of constituent below
   real                 ::  cdifr                ! Normalized diff between cabv and cbel
   real                 ::  chat(ILOTS,levs_)     ! Mix ratio in env at interfaces
   real                 ::  cond(ILOTS,levs_)     ! Mix ratio in downdraft at interfaces
   real                 ::  const(ILOTS,levs_)    ! Gathered gases array 
   real                 ::  conu(ILOTS,levs_)     ! Mix ratio in updraft at interfaces
   real                 ::  dcondt(ILOTS,levs_)   ! Gathered tend array 
   real                 ::  small                ! A small number
   real                 ::  mbsth                ! Threshold for mass fluxes
   real                 ::  mupdudp              ! A work variable
   real                 ::  minc                 ! A work variable
   real                 ::  maxc                 ! A work variable
   real                 ::  qn                   ! A work variable
   real                 ::  fluxin               ! A work variable
   real                 ::  fluxout              ! A work variable
   real                 ::  netflux              ! A work variable

!-------------------------------------------------------------------------------
   small = 1.e-36
!
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
!
   mbsth = 1.e-15
!
! Find the highest level top and bottom levels of convection
!
   ktm = levs_
   kbm = levs_
   do i = il1g, il2g
     ktm = min(ktm,jt(i))
     kbm = min(kbm,mx(i))
   end do
!
! Loop ever each constituent
!
   do m = 1,1
!
! Gather up the constituent and set tend to zero
!
     do k = 1,levs_
       do i =il1g,il2g
         const(i,k) = q(ideep(i),k,m)
       end do
     end do
!
! From now on work only with gathered data
!
! Interpolate environment gases values to interfaces
!
     do k = 1,levs_
       km1 = max(1,k-1)
       do i = il1g, il2g
         minc = min(const(i,km1),const(i,k))
         maxc = max(const(i,km1),const(i,k))
         if (minc.lt.0) then
           cdifr = 0.
         else
           cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
         endif
!
! If the two layers differ significantly use a geometric averaging
! procedure
!
         if (cdifr.gt.1.E-6) then
           cabv = max(const(i,km1),maxc*1.e-12)
           cbel = max(const(i,k),maxc*1.e-12)
           chat(i,k) = log(cabv/cbel)                                          &
                             /(cabv-cbel)                                      &
                             *cabv*cbel
         else             ! Small diff, so just arithmetic mean
           chat(i,k) = 0.5* (const(i,k)+const(i,km1))
         end if
!
! Provisional up and down draft values
!
         conu(i,k) = chat(i,k)
         cond(i,k) = chat(i,k)
!
!              provisional tends
!
         dcondt(i,k) = 0.

       end do
     end do

#ifdef DBG
     do k = 1,levs_
       km1 = max(1,k-1)
       do i = il1g, il2g
         if (chat(i,k).lt.0.) then
           write (6,*) ' negative chat ', i, k, lat, chat(i,k),                &
                    const(i,km1), const(i,k)
           stop
         endif
       end do
     end do
#endif
!
! Do levels adjacent to top and bottom
!
     k = 2
     km1 = 1
     kk = levs_ 
     do i = il1g,il2g
       mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
       if (mupdudp.gt.mbsth) then
         conu(i,kk) = (                                                        &
                          +eu(i,kk)*const(i,kk)*dp(i,kk)                       &
                         )/mupdudp
       endif
       if (md(i,k).lt.-mbsth) then
         cond(i,k) =  (                                                        &
                          -ed(i,km1)*const(i,km1)*dp(i,km1)                    &
                         )/md(i,k)
       endif
     end do
!
! Updraft from bottom to top
!
     do kk = levs_-1,1,-1
       kkp1 = min(levs_,kk+1)
       do i = il1g,il2g
         mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
         if (mupdudp.gt.mbsth) then
           conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)                             &
                              +eu(i,kk)*const(i,kk)*dp(i,kk)                   &
                            )/mupdudp
         endif
       end do
     end do
!
! Downdraft from top to bottom
!
     do k = 3,levs_
       km1 = max(1,k-1)
       do i = il1g,il2g
         if (md(i,k).lt.-mbsth) then
           cond(i,k) =  (  md(i,km1)*cond(i,km1)                               &
                              -ed(i,km1)*const(i,km1)*dp(i,km1)                &
                            )/md(i,k)
         endif
       end do
     end do
!
#ifdef DBG
     do k = ktm,levs_
       do i = il1g,il2g
!               if (conu(i,k)*mu(i,k).lt.0.) then
         if (conu(i,k).lt.0.) then
           write (6,*) ' warning negativue cu ',                               &
                    i, k, lat, m, conu(i,k), mu(i,k)
         endif
!               if (cond(i,k)*md(i,k).gt.0.) then
         if (cond(i,k).lt.0.) then
           write (6,*) ' warning negativue cd ',                               &
                    i, k, lat, m, cond(i,k)
           km1 = max(1,k-1)
           write (6,*) 'mda, cda, eda, cone, mdb ',                            &
                    md(i,km1),cond(i,km1) ,                                    &
                            -ed(i,km1),const(i,km1),md(i,k) 
         endif
       end do
     end do
#endif
!
     do k = ktm,levs_
       km1 = max(1,k-1)
       kp1 = min(levs_,k+1)
       do i = il1g,il2g
!
! version 1 hard to check for roundoff errors
!               dcondt(i,k) = 
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)
!
! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k) 
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))
!
! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
!
         fluxin =  mu(i,kp1)*conu(i,kp1)                                       &
                    + mu(i,k)*min(chat(i,k),const(i,km1))                      &
                    -(md(i,k)  *cond(i,k)                                      &
                      + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
         fluxout = mu(i,k)*conu(i,k)                                           &
                     +mu(i,kp1)*min(chat(i,kp1),const(i,k))                    &
                    -(md(i,kp1)*cond(i,kp1)                                    &
                      + md(i,k)*min(chat(i,k),const(i,k)))

         netflux = fluxin - fluxout
         if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
           netflux = 0.
         endif
         dcondt(i,k) = netflux/dp(i,k)
       end do
     end do
!
     do k = kbm,levs_             
       km1 = max(1,k-1)
       do i = il1g,il2g
         if (k.eq.mx(i)) then
!
! version 1
!
!                  dcondt(i,k) = (1./dsubcld(i))!
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )
!
! version 2
!
           fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
           fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
!
! version 3
!
           fluxin =  mu(i,k)*min(chat(i,k),const(i,km1))                       &
                       - md(i,k)*cond(i,k)
           fluxout = mu(i,k)*conu(i,k)                                         &
                    - md(i,k)*min(chat(i,k),const(i,k))

           netflux = fluxin - fluxout
           if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
             netflux = 0.
           endif
!                  dcondt(i,k) = netflux/dsubcld(i)
           dcondt(i,k) = netflux/dp(i,k)
         else if (k.gt.mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
           dcondt(i,k) = 0.
         end if
       end do
     end do
!
! Update and scatter data back to full arrays
!
     do k = 1,levs_
       kp1 = min(levs_,k+1)
       do i = il1g,il2g
         qn = const(i,k)+dcondt(i,k)*2.*delt
#ifdef DBG
         if (qn.lt.0) then
           write (6,*) ' qn less than zero ', i, k, lat, m, qn
           write (6,*) ' qo, dcondt*2*dt ', q(ideep(i),k,m),                   &
                    dcondt(i,k)*2.*delt
           write (6,*) ' jt, jb ', jt(i), mx(i)
           write (6,*) ' mu md bel ', mu(i,kp1), md(i,kp1)
           write (6,*) ' mu md abv ', mu(i,k), md(i,k)
           write (6,*) ' conu, cond abv ', conu(i,k), cond(i,k)
           write (6,*) ' conu, cond bel ', conu(i,kp1),                        &
                                               cond(i,kp1)
           write (6,*) ' chat abv, bel ', chat(i,k), chat(i,kp1)
           write (6,*) ' dp ', dp(i,k)
           fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)                 &
                    -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
           fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)             &
                    -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))
           netflux = fluxin - fluxout
           write (6,*) ' fluxin, fluxout, netflux ', fluxin,                   &
                 fluxout, netflux
           write (6,*) ' term1 ',                                              &
                       +mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))/dp(i,k)
           write (6,*) ' term2 ',                                              &
                       -mu(i,k)*   (conu(i,k)-chat(i,k))/dp(i,k)
           write (6,*) ' term3 ',                                              &
                       +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))/dp(i,k)
           write (6,*) ' term4 ',                                              &
                       -md(i,k)*   (cond(i,k)-chat(i,k))/dp(i,k)
           stop
         endif
#endif
         q(ideep(i),k,m) = qn
       end do
     end do
   end do                    ! m = 1, 1
!
   return
   end subroutine ccm3_cloud_trans
!
!-------------------------------------------------------------------------------
   subroutine ccm3_cloud_q1q2(dqdt    ,dsdt    ,q       ,s       ,qs      ,    &
                              qu      ,su      ,mc      ,du      ,alpha   ,    &
                              qhat    ,shat    ,dp      ,mu      ,md      ,    &
                              sd      ,qd      ,ql      ,dsubcld ,qdb     ,    &
                              sdb     ,betau   ,betad   ,mb      ,lcl     ,    &
                              jt      ,mx      ,dt      ,il1g    ,il2g    ,    &
                              cp      ,rl      ,msg     ,nstep   ,lat     ,    &
                              dl      ,cu      )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
!-------------------------------------------------------------------------------
!
! This is contributed code which is not fully standardized.  Both the
! internal documentation and the coding style itself are somewhat
! different from the rest of the standard CCM3 model code.  It is hoped
! that the consistency will be improved in future versions
!
! Reference: Zhang and McFarland, 1995, Atmos. Ocean
!
! PLEASE NOTE:
!
! The form of the equations used in this subroutine differ from that of
! the above reference. The forms are equivalent in the continuous form
! but differ in the descrete form. We have preserved the form as
! supplied by Zhang in order to prevent changes in the simulation that
! would result from a change to a more familiar form.
!
! This code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments
! disabled.  These are to facilitate future development.  We expect to
! release a cleaner code in a future release
!
!---------------------------Code history--------------------------------
!
! Original version:  G. Zhang
! Modified:          G. Zhang, P. Rasch and J. Hack, 1994 and 1995
! Standardized:      CCM Core Group, 1995
! Reviewed:          P. Rasch, April 1996
!
! jul 14/92 - guang jun zhang, m.lazare.  like previous version
!             ccm3_cloud_q1q2 except modify definitions of qhat and shat to
!             avoid pathological cases where very small gradients
!             exist between adjacent levels .
! feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.  previous
!             version ccm3_cloud_q1q2.
!
!-------------------------------------------------------------------------------
!
!-----------------------------Arguments---------------------------------
!
   real                 ::  dqdt(ILOTS,levs_)     ! moisture tendency (kg/kg/s)
   real                 ::  dsdt(ILOTS,levs_)     ! heating rate (K/s)
!
   real                 ::  q(ILOTS,levs_)        ! specific humidity (kg/kg)
   real                 ::  s(ILOTS,levs_)        ! dry static energy divided by Cp (K)
   real                 ::  qs(ILOTS,levs_)       ! saturation specific humidity (kg/kg)
   real                 ::  qu(ILOTS,levs_)       ! updraft spec. humidity (kg/kg)
   real                 ::  su(ILOTS,levs_)       ! updraft dry static energy divided by Cp (K)
   real                 ::  mc(ILOTS,levs_)       ! sum of mu and md
   real                 ::  du(ILOTS,levs_)       ! detraining mass from updraft
   real                 ::  alpha(ILOTS,levs_)    !
   real                 ::  qhat(ILOTS,levs_)     ! env. specific humidity at interfaces
   ! env. dry static energy divided by Cp (K) at interfaces
   real                 ::  shat(ILOTS,levs_)     
   real                 ::  dp(ILOTS,levs_)       ! pressure thickness
   real                 ::  mu(ILOTS,levs_)       ! updraft mass flux
   real                 ::  md(ILOTS,levs_)       ! downdraft mass flux
   real                 ::  sd(ILOTS,levs_)       ! downdraft normalized dry stat energy
   real                 ::  qd(ILOTS,levs_)       ! downdraft mass flux
   real                 ::  ql(ILOTS,levs_)       ! liquid water in updraft
   real                 ::  dl(ILOTS,levs_)       ! detrained liquid water
   real                 ::  cu(ILOTS,levs_)       ! condensation rate in updraft

   real                 ::  dsubcld(ILOTS)       ! pressure thickness below cloud base
   real                 ::  qdb(ILOTS)           ! 
   real                 ::  sdb(ILOTS)           ! 
   real                 ::  betau(ILOTS)         ! 
   real                 ::  betad(ILOTS)         ! 
   real                 ::  mb(ILOTS)            ! mass flux at cloud base

   integer              ::  lcl(ILOTS)        !
   integer              ::  jt(ILOTS)         ! index of top of updraft
   integer              ::  mx(ILOTS)         ! index of base of updraft

! work fields:

   real                 ::  facq(ILOTS,levs_)     ! 
   real                 ::  facs(ILOTS,levs_)     ! 
   real                 ::  dsds(ILOTS,levs_)     ! 
   real                 ::  dqdq(ILOTS,levs_)     ! 

   real                 ::  rl                   ! latent heat of vap.
!-------------------------------------------------------------------
   rnu = 0.
   do k = msg + 1,levs_
     do i = il1g,il2g
       dsdt(i,k) = 0.
       dqdt(i,k) = 0.
       if (k.ge.jt(i) .and. k.lt.levs_) then
         sdif = 0.
         if (s(i,k).gt.0. .or. s(i,k+1).gt.                                    &
                0.) sdif = abs((s(i,k)-s(i,k+1))/                              &
                max(s(i,k),s(i,k+1)))
         if (sdif.gt.1.E-6) then
           dsds(i,k) = s(i,k+1)/ (s(i,k)-s(i,k+1)) -                           &
                            log(s(i,k)/s(i,k+1))*                              &
                            (s(i,k+1)/ (s(i,k)-s(i,k+1)))**2
         else
           dsds(i,k) = 0.5
         end if
         qdif = 0.
         if (q(i,k).gt.0. .or. q(i,k+1).gt.                                    &
                0.) qdif = abs((q(i,k)-q(i,k+1))/                              &
                max(q(i,k),q(i,k+1)))
         if (qdif.gt.1.E-6) then
           dqdq(i,k) = q(i,k+1)/ (q(i,k)-q(i,k+1)) -                           &
                            log(q(i,k)/q(i,k+1))*                              &
                            (q(i,k+1)/ (q(i,k)-q(i,k+1)))**2
         else
           dqdq(i,k) = 0.5
         end if
         facs(i,k) = 1./ (1.+mc(i,k+1)*dt/dp(i,k)*rnu*                         &
                         dsds(i,k))
         facq(i,k) = 1./ (1.+mc(i,k+1)*dt/dp(i,k)*rnu*                         &
                         dqdq(i,k))
       end if
     end do
   end do
!
! this effect is included in test4
! this feature is included on dec. 8, 1994
! following tiedtke, if rh > 80% detrained liquid water has no
! cooling or moistening effect, thus effectively becomes 
! convective precipitation.
!
   do i = il1g,il2g
     if (mb(i).gt.0.) then
       fact = 1.
       if (q(i,jt(i)).gt.0.8*qs(i,jt(i)) .and.                                 &
             jt(i).lt.levs_-3) fact = 0.
       dsdt(i,jt(i)) = facs(i,jt(i))/dp(i,jt(i))*                              &
                           (mu(i,jt(i)+1)* (su(i,                              &
                           jt(i)+1)-shat(i,jt(i)+1)-                           &
                           rl/cp*fact*ql(i,jt(i)+1))+                          &
                           md(i,jt(i)+1)* (sd(i,jt(i)+1)-shat(i,               &
                           jt(i)+1)))
       dqdt(i,jt(i)) = facq(i,jt(i))/dp(i,jt(i))*                              &
                           (mu(i,jt(i)+1)* (qu(i,                              &
                           jt(i)+1)-qhat(i,jt(i)+1)+fact*ql(i,                 &
                           jt(i)+1))+md(i,jt(i)+1)*                            &
                           (qd(i,jt(i)+1)-qhat(i,jt(i)+1)))
       dl(i,jt(i)) = facq(i,jt(i))/dp(i,jt(i))*                                &
                           (mu(i,jt(i)+1)*ql(i,jt(i)+1))
     end if
   end do
!
   beta = 0.
   do k = msg + 1,levs_
     do i = il1g,il2g
       if (k.gt.jt(i) .and. k.lt.mx(i) .and. mb(i).gt.0.) then
         fact = 1.
         if (q(i,k).gt.0.8*qs(i,k) .and. k.lt.levs_-3) fact = 0.
         shat(i,k) = shat(i,k) + dsdt(i,k-1)*dt*rnu*                           &
                         dsds(i,k-1)
         qhat(i,k) = qhat(i,k) + dqdt(i,k-1)*dt*rnu*                           &
                         dqdq(i,k-1)
         dsdt(i,k) = facs(i,k)* ((mc(i,k)* (shat(i,k)-                         &
                         s(i,k))+mc(i,k+1)* (s(i,k)-                           &
                         shat(i,k+1)))/dp(i,k)-rl/cp*fact*du(i,k)*             &
                         (beta*ql(i,k)+ (1-beta)*ql(i,k+1)))
!          dqdt(i,k)=facq(i,k)
!     1               *((mc(i,k)*(qhat(i,k)-q(i,k))
!     2               +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     3               +du(i,k)*(qs(i,k)-q(i,k))
!     4         +fact*du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1)) )

         dqdt(i,k) = facq(i,k)* ((mu(i,k+1)* (qu(i,k+1)-                       &
                         qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))-               &
                         mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*                    &
                         (su(i,k)-s(i,k)))+md(i,k+1)* (qd(i,k+1)-              &
                         qhat(i,k+1)+cp/rl* (sd(i,k+1)-s(i,k)))-               &
                         md(i,k)* (qd(i,k)-qhat(i,k)+cp/rl*                    &
                         (sd(i,k)-s(i,k))))/dp(i,k)+fact*du(i,k)*              &
                         (beta*ql(i,k)+ (1-beta)*ql(i,k+1)))
         dl(i,k) = facq(i,k)*(du(i,k)*                                         &
                         (beta*ql(i,k)+ (1-beta)*ql(i,k+1)))
       end if
     end do
   end do
!
   do k = msg + 1,levs_             
     do i = il1g,il2g
       if (k.eq.mx(i) .and. mb(i).gt.0.) then
         shat(i,k) = shat(i,k) + dsdt(i,k-1)*dt*rnu*                           &
                         dsds(i,k-1)
         qhat(i,k) = qhat(i,k) + dqdt(i,k-1)*dt*rnu*                           &
                         dqdq(i,k-1)
         dsdt(i,k) = (1./dsubcld(i))*                                          &
                         (mu(i,k)* (shat(i,k)-                                 &
                         su(i,k))+md(i,k)*                                     &
                         (shat(i,k)-sd(i,k)))
         dqdt(i,k) = (1./dsubcld(i))*                                          &
                         (mu(i,k)* (qhat(i,k)-                                 &
                         qu(i,k))+md(i,k)*                                     &
                         (qhat(i,k)-qd(i,k)))
       else if (k.gt.mx(i) .and. mb(i).gt.0.) then
         dsdt(i,k) = dsdt(i,k-1)
         dqdt(i,k) = dqdt(i,k-1)
       end if
     end do
   end do
!
   return
   end subroutine ccm3_cloud_q1q2
!
!-------------------------------------------------------------------------------
   subroutine ccm2_saturation_driver                                           &
                    (lat     ,nstep   ,tdt     ,pmid    ,pdel    ,             &
                     rpdel   ,zm      ,tpert   ,qpert   ,phis    ,             &
                     pblht   ,t       ,q       ,cmfdt   ,cmfdq   ,             &
                     cmfmc   ,cmfdqr  ,cmfsl   ,cmflq   ,precc   ,             &
                     qc      ,cnt     ,cnb     )
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS, levs_
   use comcmf
!-------------------------------------------------------------------------------
!
! Moist convective mass flux procedure:
! If stratification is unstable to nonentraining parcel ascent,
! complete an adjustment making successive use of a simple cloud model
! consisting of three layers (sometimes referred to as a triplet)
! 
! Code generalized to allow specification of parcel ("updraft")
! properties, as well as convective transport of an arbitrary
! number of passive constituents (see q array).  The code
! is written so the water vapor field is passed independently
! in the calling list from the block of other transported
! constituents, even though as currently designed, it is the
! first component in the constituents field. 
!
!----------------------------Code History-------------------------------
!
! Original version:  J. J. Hack, March 22, 1990
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992
! Implemented to mrf song-you hong, september 1999
!
!-------------------------------------------------------------------------------
!
#include <implicit.h>
!------------------------------Parameters-------------------------------
   real,parameter       ::  ssfac = 1.001
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer              ::  lat              ! latitude index (S->N)
   integer              ::  nstep            ! current time step index

   real                 ::  tdt                 ! 2 delta-t (seconds)
   real                 ::  pmid(ILOTS,levs_)    ! pressure
   real                 ::  pdel(ILOTS,levs_)    ! delta-p
   real                 ::  rpdel(ILOTS,levs_)   ! 1./pdel
   real                 ::  zm(ILOTS,levs_)      ! height abv sfc at midpoints
   real                 ::  tpert(ILOTS)        ! PBL perturbation theta
   real                 ::  qpert(ILOTS,1)      ! PBL perturbation specific humidity 
   real                 ::  phis(ILOTS)         ! surface geopotential
   real                 ::  pblht(ILOTS)        ! PBL height (provided by PBL routine)
!
! Input/output arguments
!
   real                 ::  t(ILOTS,levs_)       ! temperature (t bar)
   real                 ::  q(ILOTS,levs_,1)     ! specific humidity (sh bar)
!
! Output arguments
!
   real                 ::  cmfdt(ILOTS,levs_)   ! dt/dt due to moist convection
   real                 ::  cmfdq(ILOTS,levs_)   ! dq/dt due to moist convection
   real                 ::  cmfmc(ILOTS,levs_ )  ! moist convection cloud mass flux
   real                 ::  cmfdqr(ILOTS,levs_)  ! dq/dt due to convective rainout 
   real                 ::  cmfsl(ILOTS,levs_ )  ! convective lw static energy flux
   real                 ::  cmflq(ILOTS,levs_ )  ! convective total water flux
   real                 ::  precc(ILOTS)        ! convective precipitation rate
   real                 ::  qc(ILOTS,levs_)      ! dq/dt due to rainout terms
   real                 ::  cnt(ILOTS)          ! top level of convective activity   
   real                 ::  cnb(ILOTS)          ! bottom level of convective activity
!
!---------------------------Local workspace-----------------------------
!
   real                 ::  gam(ILOTS,levs_)     ! 1/cp (d(qsat)/dT)
   real                 ::  sb(ILOTS,levs_)      ! dry static energy (s bar)
   real                 ::  hb(ILOTS,levs_)      ! moist static energy (h bar)
   real                 ::  shbs(ILOTS,levs_)    ! sat. specific humidity (sh bar star)
   real                 ::  hbs(ILOTS,levs_)     ! sat. moist static energy (h bar star)
   real                 ::  shbh(ILOTS,levs_+1)   ! specific humidity on interfaces
   real                 ::  sbh(ILOTS,levs_+1)    ! s bar on interfaces
   real                 ::  hbh(ILOTS,levs_+1)    ! h bar on interfaces
   real                 ::  cmrh(ILOTS,levs_+1)   ! interface constituent mixing ratio 
   real                 ::  prec(ILOTS)         ! instantaneous total precipitation
   real                 ::  dzcld(ILOTS)        ! depth of convective layer (m)
   real                 ::  beta(ILOTS)         ! overshoot parameter (fraction)
   real                 ::  betamx(ILOTS)       ! local maximum on overshoot
   real                 ::  eta(ILOTS)          ! convective mass flux (kg/m^2 s)
   real                 ::  etagdt(ILOTS)       ! eta*grav*dt
   real                 ::  cldwtr(ILOTS)       ! cloud water (mass)
   real                 ::  rnwtr(ILOTS)        ! rain water  (mass)
   real                 ::  sc  (ILOTS)         ! dry static energy   ("in-cloud")
   real                 ::  shc (ILOTS)         ! specific humidity   ("in-cloud")
   real                 ::  hc  (ILOTS)         ! moist static energy ("in-cloud")
   real                 ::  cmrc(ILOTS)         ! constituent mix rat ("in-cloud")
   real                 ::  dq1(ILOTS)          ! shb  convective change (lower lvl)
   real                 ::  dq2(ILOTS)          ! shb  convective change (mid level)
   real                 ::  dq3(ILOTS)          ! shb  convective change (upper lvl)
   real                 ::  ds1(ILOTS)          ! sb   convective change (lower lvl)
   real                 ::  ds2(ILOTS)          ! sb   convective change (mid level)
   real                 ::  ds3(ILOTS)          ! sb   convective change (upper lvl)
   real                 ::  dcmr1(ILOTS)        ! q convective change (lower lvl)
   real                 ::  dcmr2(ILOTS)        ! q convective change (mid level)
   real                 ::  dcmr3(ILOTS)        ! q convective change (upper lvl)
   real                 ::  estemp(ILOTS,levs_)  ! saturation vapor pressure (scratch)
   real                 ::  vtemp1(2*ILOTS)     ! intermediate scratch vector
   real                 ::  vtemp2(2*ILOTS)     ! intermediate scratch vector
   real                 ::  vtemp3(2*ILOTS)     ! intermediate scratch vector
   real                 ::  vtemp4(2*ILOTS)     ! intermediate scratch vector
   integer              ::  indx1(ILOTS)     ! longitude indices for condition true
   logical              ::  etagt0           ! true if eta > 0.0
   real                 ::  sh1                 ! dummy arg in qhalf statement func.
   real                 ::  sh2                 ! dummy arg in qhalf statement func.
   real                 ::  shbs1               ! dummy arg in qhalf statement func.
   real                 ::  shbs2               ! dummy arg in qhalf statement func.
   real                 ::  cats                ! modified characteristic adj. time
   real                 ::  rtdt                ! 1./tdt
   real                 ::  qprime              ! modified specific humidity pert.
   real                 ::  tprime              ! modified thermal perturbation
   real                 ::  pblhgt              ! bounded pbl height (max[pblh,1m])
   real                 ::  fac1                ! intermediate scratch variable
   real                 ::  shprme              ! intermediate specific humidity pert.
   real                 ::  qsattp              ! sat mix rat for thermally pert PBL parcels 
   real                 ::  dz                  ! local layer depth
   real                 ::  temp1               ! intermediate scratch variable
   real                 ::  b1                  ! bouyancy measure in detrainment lvl
   real                 ::  b2                  ! bouyancy measure in condensation lvl
   real                 ::  temp2               ! intermediate scratch variable
   real                 ::  temp3               ! intermediate scratch variable
   real                 ::  g                   ! bounded vertical gradient of hb
   real                 ::  tmass               ! total mass available for convective exch
   real                 ::  denom               ! intermediate scratch variable
   real                 ::  qtest1              ! used in negative q test (middle lvl) 
   real                 ::  qtest2              ! used in negative q test (lower lvl) 
   real                 ::  fslkp               ! flux lw static energy (bot interface)
   real                 ::  fslkm               ! flux lw static energy (top interface)
   real                 ::  fqlkp               ! flux total water (bottom interface)
   real                 ::  fqlkm               ! flux total water (top interface)
   real                 ::  botflx              ! bottom constituent mixing ratio flux
   real                 ::  topflx              ! top constituent mixing ratio flux
   real                 ::  efac1               ! ratio q to convectively induced chg (btm lvl)
   real                 ::  efac2               ! ratio q to convectively induced chg (mid lvl)
   real                 ::  efac3               ! ratio q to convectively induced chg (top lvl)
   real                 ::  tb(ILOTS,levs_)      ! working storage for temp (t bar)
   real                 ::  shb(ILOTS,levs_)     ! working storage for spec hum (sh bar)
   real                 ::  adjfac              ! adjustment factor (relaxation related)
   real                 ::  qmin
#ifdef DBG
!
!  Following 7 real variables are used in diagnostics calculations   
!
   real                 ::  rh                  ! relative humidity 
   real                 ::  es                  ! sat vapor pressure 
   real                 ::  hsum1               ! moist static energy integral
   real                 ::  qsum1               ! total water integral 
   real                 ::  hsum2               ! final moist static energy integral
   real                 ::  qsum2               ! final total water integral
   real                 ::  fac                 ! intermediate scratch variable
#endif
   integer              ::  i,k              ! longitude, level indices
   integer              ::  ii               ! index on "gathered" vectors
   integer              ::  len1             ! vector length of "gathered" vectors
   integer              ::  m                ! constituent index
   integer              ::  ktp              ! tmp indx used to track top of convective layer
#ifdef DBG
   integer              ::  n                ! vertical index     (diagnostics)
   integer              ::  kp               ! vertical index     (diagnostics)
   integer              ::  kpp              ! index offset, kp+1 (diagnostics)
   integer              ::  kpm1             ! index offset, kp-1 (diagnostics)
#endif
!
!---------------------------Statement functions-------------------------
   real                 ::  qhalf
!
   qhalf(sh1,sh2,shbs1,shbs2) = min(max(sh1,sh2),                              &
                               (shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
!
!-------------------------------------------------------------------------------
!
! set up convection cap, this was moved from ccm3_constant_setup to here to
! get rid of dependency on hybrid pressure coefficients for the
! Suarez repos code
!
   if (pmid(1,1) .ge. 4.e3) then 
     limcnv =1
   else
     do k = 1,levs_-1
       if (pmid(1,k).lt.4.e3 .and. pmid(1,k+1).ge.4.e3) then
         limcnv = k
!
         goto 1010
!
       end if
     end do
     limcnv = levs_+1
   end if
!
 1010 continue

!
! Ensure that characteristic adjustment time scale (cmftau) assumed
! in estimate of eta isn't smaller than model time scale (tdt)
! The time over which the convection is assumed to act (the adjustment
! time scale) can be applied with each application of the three-level
! cloud model, or applied to the column tendencies after a "hard"
! adjustment (i.e., on a 2-delta t time scale) is evaluated
!
   if (rlxclm) then
     cats   = tdt               ! relaxation applied to column
     adjfac = tdt/(max(tdt,cmftau))
   else
     cats   = max(tdt,cmftau) ! relaxation applied to triplet
     adjfac = 1.0
   endif
   rtdt = 1.0/tdt
!
! Move temperature and moisture into working storage
!
   qmin = 1.E-20          ! blow up near tropause with negative values
   do k = limcnv,levs_
     do i = 1,ILOTS
       tb (i,k) = t(i,k)
       shb(i,k) = max(q(i,k,1),qmin)
     end do
   end do
!
! Compute sb,hb,shbs,hbs
!
   call ccm2_saturation_table1                                                 &
              (tb      ,pmid    ,estemp ,shbs    ,gam     ,                    &
               ILOTS   ,ILOTS    ,levs_   ,limcnv  ,levs_    )
!
   do k = limcnv,levs_
     do i = 1,ILOTS
       sb (i,k) = cp*tb(i,k) + zm(i,k)*grav + phis(i)
       hb (i,k) = sb(i,k) + hlat*shb(i,k)
       hbs(i,k) = sb(i,k) + hlat*shbs(i,k)
     end do
   end do
!
! Compute sbh, shbh
!
   do k = limcnv+1,levs_
     do i = 1,ILOTS
       sbh (i,k) = 0.5*(sb(i,k-1) + sb(i,k))
       shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
       hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
     end do
   end do
!
! Specify properties at top of model (not used, but filling anyway)
!
   do i = 1,ILOTS
     sbh (i,limcnv) = sb(i,limcnv)
     shbh(i,limcnv) = shb(i,limcnv)
     hbh (i,limcnv) = hb(i,limcnv)
   end do
!
! Zero vertically independent control, tendency & diagnostic arrays
!
   do i = 1,ILOTS
     prec(i)  = 0.0
     dzcld(i) = 0.0
     cnb(i)   = 0.0
     cnt(i)   = float(levs_+1)
   end do
#ifdef DBG
!#######################################################################
!#                                                                     #
!#    output initial thermodynamic profile if debug diagnostics        #
!#                                                                     #
   if (lat.eq.jloc .and. nstep.ge.nsloc) then
     i = iloc
!#                                                                     #
!#       approximate vertical integral of moist static energy          #
!#       and total preciptable water                                   #
!#                                                                     #
     hsum1 = 0.0
     qsum1 = 0.0
     do k = limcnv,levs_
       hsum1 = hsum1 + pdel(i,k)*rgrav*hb(i,k)
       qsum1 = qsum1 + pdel(i,k)*rgrav*shb(i,k)
     end do
!#                                                                     #
     write (6,8010)
     fac = grav*864.
     do k = limcnv,levs_
       rh = shb(i,k)/shbs(i,k)
       write(6,8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc(i,k),               &
            cmfsl(i,k), cmflq(i,k)
       write(6,8040) tb(i,k),shb(i,k),rh,sb(i,k),hb(i,k),hbs(i,k),             &
            tdt*cmfdt(i,k),tdt*cmfdq(i,k),tdt*cmfdqr(i,k)
     end do
     write(6, 8000) prec(i)
   end if
#endif
!#                                                                     #
!#                                                                     #
!#######################################################################
!
! Begin moist convective mass flux adjustment procedure.
! Formalism ensures that negative cloud liquid water can never occur
!
   do 70 k = levs_-1,limcnv+1,-1
     do 10 i = 1,ILOTS
       etagdt(i) = 0.0
       eta   (i) = 0.0
       beta  (i) = 0.0
       ds1   (i) = 0.0
       ds2   (i) = 0.0
       ds3   (i) = 0.0
       dq1   (i) = 0.0
       dq2   (i) = 0.0
       dq3   (i) = 0.0
!
! Specification of "cloud base" conditions
!
       qprime    = 0.0
       tprime    = 0.0
!
! Assign tprime within the PBL to be proportional to the quantity
! tpert (which will be bounded by tpmax), passed to this routine by 
! the PBL routine.  Don't allow perturbation to produce a dry 
! adiabatically unstable parcel.  Assign qprime within the PBL to be 
! an appropriately modified value of the quantity qpert (which will be 
! bounded by shpmax) passed to this routine by the PBL routine.  The 
! quantity qprime should be less than the local saturation value 
! (qsattp=qsat[t+tprime,p]).  In both cases, tpert and qpert are
! linearly reduced toward zero as the PBL top is approached.
!
       pblhgt = max(pblht(i),1.0)
       if( (zm(i,k+1) .le. pblhgt) .and.dzcld(i).eq.0.0 ) then
         fac1   = max(0.0,1.0-zm(i,k+1)/pblhgt)
         tprime = min(tpert(i),tpmax)*fac1
         qsattp = shbs(i,k+1) + cp*rhlat*gam(i,k+1)*tprime
         shprme = min(min(qpert(i,1),shpmax)*fac1,max(qsattp-shb(i,k+1),0.0))
         qprime = max(qprime,shprme)
       else
         tprime = 0.0
         qprime = 0.0
       end if
!
! Specify "updraft" (in-cloud) thermodynamic properties
!
       sc (i)    = sb (i,k+1) + cp*tprime
       shc(i)    = shb(i,k+1) + qprime
       hc (i)    = sc (i    ) + hlat*shc(i)
       vtemp4(i) = hc(i) - hbs(i,k)
       dz        = pdel(i,k)*rgas*tb(i,k)*rgrav/pmid(i,k)
       if (vtemp4(i).gt.0.0) then
         dzcld(i) = dzcld(i) + dz
       else
         dzcld(i) = 0.0
       end if
   10   continue
#ifdef DBG
!#######################################################################
!#                                                                     #
!#    output thermodynamic perturbation information                    #
!#                                                                     #
     if (lat.eq.jloc .and. nstep.ge.nsloc) then
       write (6,8090) k+1,sc(iloc),shc(iloc),hc(iloc)
     end if
!#                                                                     #
!#######################################################################
#endif
!
! Check on moist convective instability
! Build index vector of points where instability exists
!
     call ccm3_cloud_trigger(ILOTS,vtemp4,1,0.0,indx1,len1)
     if (len1.le.0) go to 70
!
! Current level just below top level => no overshoot
!
     if (k.le.limcnv+1) then
       do ii = 1,len1
         i = indx1(ii)
         temp1     = vtemp4(i)/(1.0 + gam(i,k))
         cldwtr(i) = max(0.0,(sb(i,k) - sc(i) + temp1))
         beta(i)   = 0.0
         vtemp3(i) = (1.0 + gam(i,k))*(sc(i) - sbh(i,k))
       end do
     else
!
! First guess at overshoot parameter using crude buoyancy closure
! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
! If pre-existing supersaturation in detrainment layer, beta=0
! cldwtr is temporarily equal to hlat*l (l=> liquid water)
!
!DIR$ IVDEP
       do ii = 1,len1
         i = indx1(ii)
         temp1     = vtemp4(i)/(1.0 + gam(i,k))
         cldwtr(i) = max(0.0,(sb(i,k)-sc(i)+temp1))
         betamx(i) = 1.0 - c0*max(0.0,(dzcld(i)-dzmin))
         b1        = (hc(i) - hbs(i,k-1))*pdel(i,k-1)
         b2        = (hc(i) - hbs(i,k  ))*pdel(i,k  )
         beta(i)   = max(betamn,min(betamx(i), 1.0 + b1/b2))
         if (hbs(i,k-1).le.hb(i,k-1)) beta(i) = 0.0
!
! Bound maximum beta to ensure physically realistic solutions
!
! First check constrains beta so that eta remains positive
! (assuming that eta is already positive for beta equal zero)
!
         vtemp1(i) = -(hbh(i,k+1) - hc(i))*pdel(i,k)*rpdel(i,k+1)+             &
                    (1.0 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i))
         vtemp2(i) = (1.0 + gam(i,k))*(sc(i) - sbh(i,k))
         vtemp3(i) = vtemp2(i)
         if ((beta(i)*vtemp2(i) - vtemp1(i)).gt.0.) then
           betamx(i) = 0.99*(vtemp1(i)/vtemp2(i))
           beta(i)   = max(0.0,min(betamx(i),beta(i)))
         end if
       end do
!
! Second check involves supersaturation of "detrainment layer"
! small amount of supersaturation acceptable (by ssfac factor)
!
!DIR$ IVDEP
       do ii = 1,len1
         i = indx1(ii)
         if (hb(i,k-1).lt.hbs(i,k-1)) then
           vtemp1(i) = vtemp1(i)*rpdel(i,k)
           temp2 = gam(i,k-1)*(sbh(i,k) - sc(i) + cldwtr(i)) -                 &
                      hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
           temp3 = vtemp3(i)*rpdel(i,k)
           vtemp2(i) = (tdt/cats)*(hc(i) - hbs(i,k))*temp2/                    &
                      (pdel(i,k-1)*(hbs(i,k-1) - hb(i,k-1))) + temp3
           if ((beta(i)*vtemp2(i) - vtemp1(i)).gt.0.) then
             betamx(i) = ssfac*(vtemp1(i)/vtemp2(i))
             beta(i)   = max(0.0,min(betamx(i),beta(i)))
           end if
         else 
           beta(i) = 0.0
         end if
       end do
!
! Third check to avoid introducing 2 delta x thermodynamic
! noise in the vertical ... constrain adjusted h (or theta e)
! so that the adjustment doesn't contribute to "kinks" in h
!
!DIR$ IVDEP
       do ii = 1,len1
         i = indx1(ii)
         g = min(0.0,hb(i,k) - hb(i,k-1))
         temp1 = (hb(i,k) - hb(i,k-1) - g)*(cats/tdt)/                         &
                    (hc(i) - hbs(i,k))
         vtemp1(i) = temp1*vtemp1(i) + (hc(i) - hbh(i,k+1))*                   &
                        rpdel(i,k)
         vtemp2(i) = temp1*vtemp3(i)*rpdel(i,k) +                              &
                        (hc(i) - hbh(i,k) - cldwtr(i))*                        &
                        (rpdel(i,k) + rpdel(i,k+1))
         if ((beta(i)*vtemp2(i) - vtemp1(i)).gt.0.) then
           if (vtemp2(i).ne.0.0) then
             betamx(i) = vtemp1(i)/vtemp2(i)
           else
             betamx(i) = 0.0
           end if
           beta(i) = max(0.0,min(betamx(i),beta(i)))
         end if
       end do
     end if
!
! Calculate mass flux required for stabilization.
!
! Ensure that the convective mass flux, eta, is positive by
! setting negative values of eta to zero..
! Ensure that estimated mass flux cannot move more than the
! minimum of total mass contained in either layer k or layer k+1.
! Also test for other pathological cases that result in non-
! physical states and adjust eta accordingly.
!
!DIR$ IVDEP
     do ii = 1,len1
       i = indx1(ii)
       beta(i) = max(0.0,beta(i))
       temp1 = hc(i) - hbs(i,k)
       temp2 = ((1.0 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i)) -            &
                  beta(i)*vtemp3(i))*rpdel(i,k) -                              &
                 (hbh(i,k+1) - hc(i))*rpdel(i,k+1)
       eta(i) = temp1/(temp2*grav*cats)
       tmass = min(pdel(i,k),pdel(i,k+1))*rgrav
       if (eta(i).gt.tmass*rtdt .or. eta(i).le.0.0) eta(i) = 0.0
!
! Check on negative q in top layer (bound beta)
!
       if (shc(i)-shbh(i,k).lt.0.0 .and. beta(i)*eta(i).ne.0.0) then
         denom = eta(i)*grav*tdt*(shc(i) - shbh(i,k))*rpdel(i,k-1)
         beta(i) = max(0.0,min(-0.999*shb(i,k-1)/denom,                        &
                                      beta(i)))
       end if
!
! Check on negative q in middle layer (zero eta)
!
       qtest1 = shb(i,k) + eta(i)*grav*tdt*((shc(i) - shbh(i,k+1)) -           &
                  (1.0 - beta(i))*cldwtr(i)*rhlat -                            &
                  beta(i)*(shc(i) - shbh(i,k)))*rpdel(i,k)
       if (qtest1.le.0.0) eta(i) = 0.0
!
! Check on negative q in lower layer (bound eta)
!
       fac1 = -(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
       qtest2 = shb(i,k+1) - eta(i)*grav*tdt*fac1
       if (qtest2 .lt. 0.0) then
         eta(i) = 0.99*shb(i,k+1)/(grav*tdt*fac1)
       end if
       etagdt(i) = eta(i)*grav*tdt
     end do
!
#ifdef DBG
!#######################################################################
!#                                                                     #
     if (lat.eq.jloc .and. nstep.ge.nsloc) then
       write(6,8080) beta(iloc), eta(iloc)
     end if
!#                                                                     #
!#######################################################################
#endif
!
! Calculate cloud water, rain water, and thermodynamic changes
!
!DIR$ IVDEP
     do 30 ii = 1,len1
       i = indx1(ii)
       cldwtr(i) = etagdt(i)*cldwtr(i)*rhlat*rgrav
       rnwtr(i) = (1.0 - beta(i))*cldwtr(i)
       ds1(i) = etagdt(i)*(sbh(i,k+1) - sc(i))*rpdel(i,k+1)
       dq1(i) = etagdt(i)*(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
       ds2(i) = (etagdt(i)*(sc(i) - sbh(i,k+1)) +                              &
                   hlat*grav*cldwtr(i) - beta(i)*etagdt(i)*                    &
                   (sc(i) - sbh(i,k)))*rpdel(i,k)
       dq2(i) = (etagdt(i)*(shc(i) - shbh(i,k+1)) -                            &
                   grav*rnwtr(i) - beta(i)*etagdt(i)*                          &
                   (shc(i) - shbh(i,k)))*rpdel(i,k)
       ds3(i) = beta(i)*(etagdt(i)*(sc(i) - sbh(i,k)) -                        &
                  hlat*grav*cldwtr(i))*rpdel(i,k-1)
       dq3(i) = beta(i)*etagdt(i)*(shc(i) - shbh(i,k))*rpdel(i,k-1)
!
! Isolate convective fluxes for later diagnostics
!
       fslkp = eta(i)*(sc(i) - sbh(i,k+1))
       fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) -                            &
                          hlat*cldwtr(i)*rtdt)
       fqlkp = eta(i)*(shc(i) - shbh(i,k+1))
       fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))
!
! Update thermodynamic profile (update sb, hb, & hbs later)
!
       tb (i,k+1) = tb(i,k+1)  + ds1(i)*rcp
       tb (i,k  ) = tb(i,k  )  + ds2(i)*rcp
       tb (i,k-1) = tb(i,k-1)  + ds3(i)*rcp
       shb(i,k+1) = shb(i,k+1) + dq1(i)
       shb(i,k  ) = shb(i,k  ) + dq2(i)
       shb(i,k-1) = shb(i,k-1) + dq3(i)
!
! ** Update diagnostic information for final budget *!
! Tracking precipitation, temperature & specific humidity tendencies,
! rainout term, convective mass flux, convective liquid
! water static energy flux, and convective total water flux
! The variable afac makes the necessary adjustment to the
! diagnostic fluxes to account for adjustment time scale based on
! how relaxation time scale is to be applied (column vs. triplet)
!
       prec(i)    = prec(i) + (rnwtr(i)/rhoh2o)*adjfac
!
! The following variables have units of "units"/second
!
       cmfdt (i,k+1) = cmfdt (i,k+1) + ds1(i)*rcp*rtdt*adjfac
       cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)*rcp*rtdt*adjfac
       cmfdt (i,k-1) = cmfdt (i,k-1) + ds3(i)*rcp*rtdt*adjfac
       cmfdq (i,k+1) = cmfdq (i,k+1) + dq1(i)*rtdt*adjfac
       cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rtdt*adjfac
       cmfdq (i,k-1) = cmfdq (i,k-1) + dq3(i)*rtdt*adjfac
       qc    (i,k  ) = (grav*rnwtr(i)*rpdel(i,k))*rtdt*adjfac
       cmfdqr(i,k  ) = cmfdqr(i,k  ) + qc(i,k)
       cmfmc (i,k+1) = cmfmc (i,k+1) + eta(i)*adjfac
       cmfmc (i,k  ) = cmfmc (i,k  ) + beta(i)*eta(i)*adjfac
!
! The following variables have units of w/m**2
!
       cmfsl (i,k+1) = cmfsl (i,k+1) + fslkp*adjfac
       cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm*adjfac
       cmflq (i,k+1) = cmflq (i,k+1) + hlat*fqlkp*adjfac
       cmflq (i,k  ) = cmflq (i,k  ) + hlat*fqlkm*adjfac
   30   continue
!
! Next, convectively modify passive constituents
! For now, when applying relaxation time scale to thermal fields after 
! entire column has undergone convective overturning, constituents will 
! be mixed using a "relaxed" value of the mass flux determined above
! Although this will be inconsistant with the treatment of the thermal
! fields, it's computationally much cheaper, no more-or-less justifiable,
! and consistent with how the history tape mass fluxes would be used in
! an off-line mode (i.e., using an off-line transport model)
!
   if(1.ge.1) then
     do 50 m = 1,1         ! note: indexing assumes water is first field
!DIR$ IVDEP
       do 40 ii = 1,len1
         i = indx1(ii)
!
! If any of the reported values of the constituent is negative in
! the three adjacent levels, nothing will be done to the profile
!
         if ((q(i,k+1,m).lt.0.0) .or. (q(i,k,m).lt.0.0) .or.                   &
              (q(i,k-1,m).lt.0.0)) go to 40
!
! Specify constituent interface values (linear interpolation)
!
         cmrh(i,k  ) = 0.5*(q(i,k-1,m) + q(i,k  ,m))
         cmrh(i,k+1) = 0.5*(q(i,k  ,m) + q(i,k+1,m))
!
! Specify perturbation properties of constituents in PBL
!
         pblhgt = max(pblht(i),1.0)
         if( (zm(i,k+1) .le. pblhgt).and.dzcld(i).eq.0.0 ) then
           fac1 = max(0.0,1.0-zm(i,k+1)/pblhgt)
           cmrc(i) = q(i,k+1,m) + qpert(i,m)*fac1
         else
           cmrc(i) = q(i,k+1,m)
         end if
!
! Determine fluxes, flux divergence => changes due to convection
! Logic must be included to avoid producing negative values. A bit
! messy since there are no a priori assumptions about profiles.
! Tendency is modified (reduced) when pending disaster detected.
!
         botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
         topflx   = beta(i)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
         dcmr1(i) = -botflx*rpdel(i,k+1)
         efac1    = 1.0
         efac2    = 1.0
         efac3    = 1.0
!
         if (q(i,k+1,m)+dcmr1(i) .lt. 0.0) then
           efac1 = max(tiny,abs(q(i,k+1,m)/dcmr1(i)) - eps)
         end if
!
         if (efac1.eq.tiny .or. efac1.gt.1.0) efac1 = 0.0
         dcmr1(i) = -efac1*botflx*rpdel(i,k+1)
         dcmr2(i) = (efac1*botflx - topflx)*rpdel(i,k)
!  
         if (q(i,k,m)+dcmr2(i) .lt. 0.0) then
           efac2 = max(tiny,abs(q(i,k  ,m)/dcmr2(i)) - eps)
         end if
!
         if (efac2.eq.tiny .or. efac2.gt.1.0) efac2 = 0.0
         dcmr2(i) = (efac1*botflx - efac2*topflx)*rpdel(i,k)
         dcmr3(i) = efac2*topflx*rpdel(i,k-1)
!
         if (q(i,k-1,m)+dcmr3(i) .lt. 0.0) then
           efac3 = max(tiny,abs(q(i,k-1,m)/dcmr3(i)) - eps)
         end if
!
         if (efac3.eq.tiny .or. efac3.gt.1.0) efac3 = 0.0
         efac3    = min(efac2,efac3)
         dcmr2(i) = (efac1*botflx - efac3*topflx)*rpdel(i,k)
         dcmr3(i) = efac3*topflx*rpdel(i,k-1)
!
         q(i,k+1,m) = q(i,k+1,m) + dcmr1(i)
         q(i,k  ,m) = q(i,k  ,m) + dcmr2(i)
         q(i,k-1,m) = q(i,k-1,m) + dcmr3(i)
   40     continue
   50   continue                ! end of m=1,1 loop
   endif
!
! Constituent modifications complete
!
     if (k.eq.limcnv+1) go to 60
!
! Complete update of thermodynamic structure at integer levels
! gather/scatter points that need new values of shbs and gamma
!
     do ii = 1,len1
       i = indx1(ii)
       vtemp1(ii     ) = tb(i,k)
       vtemp1(ii+len1) = tb(i,k-1)
       vtemp2(ii     ) = pmid(i,k)
       vtemp2(ii+len1) = pmid(i,k-1)
     end do
     call ccm2_saturation_table2                                               &
                 (vtemp1  ,vtemp2  ,estemp  ,vtemp3  , vtemp4  ,               &
                  2*len1   )    ! using estemp as extra long vector
!DIR$ IVDEP
     do ii = 1,len1
       i = indx1(ii)
       shbs(i,k  ) = vtemp3(ii     )
       shbs(i,k-1) = vtemp3(ii+len1)
       gam(i,k  ) = vtemp4(ii     )
       gam(i,k-1) = vtemp4(ii+len1)
       sb (i,k  ) = sb(i,k  ) + ds2(i)
       sb (i,k-1) = sb(i,k-1) + ds3(i)
       hb (i,k  ) = sb(i,k  ) + hlat*shb(i,k  )
       hb (i,k-1) = sb(i,k-1) + hlat*shb(i,k-1)
       hbs(i,k  ) = sb(i,k  ) + hlat*shbs(i,k  )
       hbs(i,k-1) = sb(i,k-1) + hlat*shbs(i,k-1)
     end do
!
! Update thermodynamic information at half (i.e., interface) levels
!
!DIR$ IVDEP
     do ii = 1,len1
       i = indx1(ii)
       sbh (i,k) = 0.5*(sb(i,k) + sb(i,k-1))
       shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
       hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
       sbh (i,k-1) = 0.5*(sb(i,k-1) + sb(i,k-2))
       shbh(i,k-1) = qhalf(shb(i,k-2),shb(i,k-1),                              &
                           shbs(i,k-2),shbs(i,k-1))
       hbh (i,k-1) = sbh(i,k-1) + hlat*shbh(i,k-1)
     end do
!
#ifdef DBG
!#######################################################################
!#                                                                     #
!#    this update necessary, only if debugging diagnostics requested   #
!#                                                                     #
     if (lat.eq.jloc .and. nstep.ge.nsloc) then
       do i = 1,ILOTS
         call ccm2_saturation_table2                                           &
                  (tb(i,k+1),pmid(i,k+1),estemp,shbs(i,k+1),gam(i,k+1),1)
         sb (i,k+1) = sb(i,k+1) + ds1(i)
         hb (i,k+1) = sb(i,k+1) + hlat*shb(i,k+1)
         hbs(i,k+1) = sb(i,k+1) + hlat*shbs(i,k+1)
         kpp = k + 2
         if(k+1.eq.levs_) kpp = k + 1
         do kp=k+1,kpp
           kpm1 = kp-1
           sbh(i,kp)  = 0.5*(sb(i,kpm1) + sb(i,kp))
           shbh(i,kp) = qhalf(shb(i,kpm1),shb(i,kp),shbs(i,kpm1),              &
                           shbs(i,kp))
           hbh(i,kp)  = sbh(i,kp) + hlat*shbh(i,kp)
         end do
       end do
!#                                                                     #
!#          diagnostic output                                          #
!#                                                                     #
       i = iloc
       write(6, 8060) k
       fac = grav*864.
       do n = limcnv,levs_
         rh  = shb(i,n)/shbs(i,n)
         write(6,8020)shbh(i,n),sbh(i,n),hbh(i,n),fac*cmfmc(i,n),              &
                 cmfsl(i,n), cmflq(i,n)
!--------------write(6, 8050)
!--------------write(6, 8030) fac*cmfmc(i,n),cmfsl(i,n), cmflq(i,n)
         write(6, 8040) tb(i,n),shb(i,n),rh,sb(i,n),hb(i,n),                   &
                 hbs(i,n), tdt*cmfdt(i,n),tdt*cmfdq(i,n),                      &
              tdt*cmfdqr(i,n)
       end do
       write(6, 8000) prec(i)
     end if
!#                                                                     #
!#                                                                     #
!#######################################################################
#endif
!
! Ensure that dzcld is reset if convective mass flux zero
! specify the current vertical extent of the convective activity
! top of convective layer determined by size of overshoot param.
!
   60   do i = 1,ILOTS
       etagt0 = eta(i).gt.0.0
       if (.not.etagt0) dzcld(i) = 0.0
       if (etagt0 .and. beta(i).gt.betamn) then
         ktp = k-1
       else
         ktp = k
       end if
       if (etagt0) then
         cnt(i) = min(cnt(i),float(ktp))
         cnb(i) = max(cnb(i),float(k))
       end if
     end do
   70 continue                  ! end of k loop
!
! ** apply final thermodynamic tendencies *!
!
   do k = limcnv,levs_
     do i = 1,ILOTS
       t (i,k) = t (i,k) + cmfdt(i,k)*tdt
       q(i,k,1) = q(i,k,1) + cmfdq(i,k)*tdt
     end do
   end do
!
! Kludge to prevent cnb-cnt from being zero (in the event
! someone decides that they want to divide by this quantity)
!
   do i = 1,ILOTS
     if (cnb(i).ne.0.0 .and. cnb(i).eq.cnt(i)) then
       cnt(i) = cnt(i) - 1.0
     end if
   end do
!
   do i = 1,ILOTS
     precc(i) = prec(i)*rtdt
   end do
!
#ifdef DBG
!#######################################################################
!#                                                                     #
!#    we're done ... show final result if debug diagnostics requested  #
!#                                                                     #
   if (lat.eq.jloc .and. nstep.ge.nsloc) then
     i=iloc
     fac = grav*864.
     write(6, 8010)
     do k = limcnv,levs_
       rh = shb(i,k)/shbs(i,k)
       write(6, 8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc(i,k),              &
              cmfsl(i,k), cmflq(i,k)
       write(6, 8040) tb(i,k),shb(i,k),rh   ,sb(i,k),hb(i,k),                  &
              hbs(i,k), tdt*cmfdt(i,k),tdt*cmfdq(i,k),tdt*cmfdqr(i,k)
     end do
     write(6, 8000) prec(i)
!#                                                                     #
!#       approximate vertical integral of moist static energy and      #
!#       total preciptable water after adjustment and output changes   #
!#                                                                     #
     hsum2 = 0.0
     qsum2 = 0.0
     do k = limcnv,levs_
       hsum2 = hsum2 + pdel(i,k)*rgrav*hb(i,k)
       qsum2 = qsum2 + pdel(i,k)*rgrav*shb(i,k)
     end do
!#                                                                     #
     write (6,8070) hsum1, hsum2, abs(hsum2-hsum1)/hsum2,                      &
           qsum1, qsum2, abs(qsum2-qsum1)/qsum2
   end if
!#                                                                     #
!#                                                                     #
!#######################################################################
#endif
   return                 ! we're all done ... return to calling procedure
#ifdef DBG
!
! Formats
!
 8000 format(///,10x,'PREC = ',3pf12.6,/)
 8010 format('1**        TB      SHB      RH       SB',                        &
        '       HB      HBS      CAH      CAM       PRECC ',                   &
        '     ETA      FSL       FLQ     **', /)
 8020 format(' ----- ',     9x,3p,f7.3,2x,2p,     9x,-3p,f7.0,2x,              &
        f7.0, 37x, 0p,2x,f8.2,  0p,2x,f8.2,2x,f8.2, ' ----- ')
 8030 format(' ----- ',  0p,82x,f8.2,  0p,2x,f8.2,2x,f8.2,                     &
          ' ----- ')
 8040 format(' - - - ',f7.3,2x,3p,f7.3,2x,2p,f7.3,2x,-3p,f7.0,2x,              &
        f7.0, 2x,f8.0,2x,0p,f7.3,3p,2x,f7.3,2x,f7.3,30x,                       &
         ' - - - ')
 8050 format(' ----- ',110x,' ----- ')
 8060 format('1 K =>',  i4,/,                                                  &
           '           TB      SHB      RH       SB',                          &
           '       HB      HBS      CAH      CAM       PREC ',                 &
           '     ETA      FSL       FLQ', /)
 8070 format(' VERTICALLY INTEGRATED MOIST STATIC ENERGY BEFORE, AFTER',       &
        ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/,               &
        ' VERTICALLY INTEGRATED MOISTURE            BEFORE, AFTER',            &
        ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/)
 8080       format(' BETA, ETA => ', 1p,2e12.3)
 8090 format (' k+1, sc, shc, hc => ', 1x, i2, 1p, 3e12.4)
#endif
!
   end subroutine ccm2_saturation_driver
!
!-------------------------------------------------------------------------------
   subroutine ccm2_saturation_table1                                           &
           (t       ,p       ,es      ,qs      ,gam     ,                      &
                     ii      ,ilen    ,kk      ,kstart  ,kend    )
!-------------------------------------------------------------------------------
!
! Utility procedure to look up and return saturation vapor pressure from 
! precomputed table, calculate and return saturation specific humidity 
! (g/g).   Differs from aqsat by also calculating and returning
!  gamma (l/cp)*(d(qsat)/dT)
! Input arrays temperature and pressure (dimensioned ii,kk).
!
!---------------------------Code history--------------------------------
!
! Original version: J. Hack, Feb 1990
! Standardized:     L. Buja, Feb 1996
! Reviewed:         J. Hack, Aug 1992 
!                   J. Hack, P. Rasch, Apr 1996 
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer  ::  ii            ! I dimension of arrays t, p, es, qs
   integer  ::  kk            ! K dimension of arrays t, p, es, qs 
   real     ::  t(ii,kk)      ! Temperature
   real     ::  p(ii,kk)      ! Pressure
   integer  ::  ilen          ! Vector length in I direction
   integer  ::  kstart        ! Starting location in K direction
   integer  ::  kend          ! Ending location in K direction
! 
! Output arguments
!
   real     ::  es(ii,kk)        ! Saturation vapor pressure
   real     ::  qs(ii,kk)        ! Saturation specific humidity
   real     ::  gam(ii,kk)       ! (l/cp)*(d(qs)/dt)
!
!---------------------------Local workspace-----------------------------
!
   logical  ::  lflg          ! True if in temperature transition region
   integer  ::  i             ! i index for vector calculations
   integer  ::  k             ! k index 
   real     ::  omeps            ! 1. - 0.622
   real     ::  trinv            ! Reciprocal of ttrice (transition range)
   real     ::  tc               ! Temperature (in degrees C)
   real     ::  weight           ! Weight for es transition from water to ice
   real     ::  hltalt           ! Appropriately modified hlat for T derivatives
   real     ::  hlatsb           ! hlat weighted in transition region
   real     ::  hlatvp           ! hlat modified for t changes above 273.16
   real     ::  tterm            ! Account for d(es)/dT in transition region
   real     ::  desdt            ! d(es)/dT
!
!--------------------------Statement functions--------------------------
#include <eslookup.h>
!-------------------------------------------------------------------------------
!
   omeps = 1.0 - epsqs
   do k = kstart,kend
     do i = 1,ilen
       es(i,k) = estblf(t(i,k))
!
! Saturation specific humidity
!
       qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
!
! The following check is to avoid the generation of negative qs
! values which can occur in the upper stratosphere and mesosphere
!
       qs(i,k) = min(1.0,qs(i,k))
!
       if (qs(i,k) .lt. 0.0) then
         qs(i,k) = 1.0
         es(i,k) = p(i,k)
       end if
     end do
   end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
   trinv = 0.0
   if ((.not. icephs) .or. (ttrice.eq.0.0)) go to 10
   trinv = 1.0/ttrice
!
   do k = kstart,kend
     do i = 1,ilen
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above 273.16 where constant slope is given by -2369 j/(kg c) =cpv - cw
!
       tc     = t(i,k) - 273.16
       lflg   = (tc.ge.-ttrice .and. tc.lt.0.0)
       weight = min(-tc*trinv,1.0)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0*tc
       if (t(i,k).lt.273.16) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       if (lflg) then
         tterm = pcf(1) +                                                      &
                 tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       else
         tterm = 0.0
       end if
       desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
       gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/                                 &
                    (cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
       if(qs(i,k).eq.1.0) gam(i,k) = 0.0
     end do
   end do
!
   go to 20
!
! No icephs or water to ice transition
!
   10 do k = kstart,kend
     do i = 1,ilen
!
! Account for change of hlatv with t above 273.16 where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0*(t(i,k)-273.16)
       if (icephs) then
         hlatsb = hlatv + hlatf
       else
         hlatsb = hlatv
       end if
       if (t(i,k).lt.273.16) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
       gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/                                 &
                     (cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
       if (qs(i,k) .eq. 1.0) gam(i,k) = 0.0
     end do
   end do
!
   20 return
   end subroutine ccm2_saturation_table1
!-------------------------------------------------------------------------------
   subroutine ccm2_saturation_table2                                           &
                   (t       ,p       ,es      ,qs      ,gam      , len     )
!-------------------------------------------------------------------------------
!
! Utility procedure to look up and return saturation vapor pressure from 
! precomputed table, calculate and return saturation specific humidity 
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure
!
!----------------------------Code History-------------------------------
!
! Original version:  J. Hack
! Standardized:      J. Rosinski, June 1992
!                    T. Acker, March 1996
! Reviewed:          J. Hack, August 1992
!
!-------------------------------------------------------------------------------
#include <implicit.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer              ::  len       ! vector length
   real                 ::  t(len)       ! temperature
   real                 ::  p(len)       ! pressure
! 
! Output arguments
!
   real                 ::  es(len)   ! saturation vapor pressure
   real                 ::  qs(len)   ! saturation specific humidity
   real                 ::  gam(len)  ! (l/cp)*(d(qs)/dt)
!
!--------------------------Local Variables------------------------------
!
   logical              ::  lflg   ! true if in temperature transition region
!
   integer              ::  i      ! index for vector calculations
!
   real                 ::  omeps     ! 1. - 0.622
   real                 ::  trinv     ! reciprocal of ttrice (transition range)
   real                 ::  tc        ! temperature (in degrees C)
   real                 ::  weight    ! weight for es transition from water to ice
   real                 ::  hltalt    ! appropriately modified hlat for T derivatives  
!
   real                 ::  hlatsb    ! hlat weighted in transition region
   real                 ::  hlatvp    ! hlat modified for t changes above 273.16
   real                 ::  tterm     ! account for d(es)/dT in transition region
   real                 ::  desdt     ! d(es)/dT
!
!-------------------------------------------------------------------------------
#include <eslookup.h>
!-------------------------------------------------------------------------------
   omeps = 1.0 - epsqs
   do i = 1,len
     es(i) = estblf(t(i))
!
! Saturation specific humidity
!
     qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
     qs(i) = min(1.0,qs(i))
!
     if (qs(i) .lt. 0.0) then
       qs(i) = 1.0
       es(i) = p(i)
     end if
   end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
   trinv = 0.0
   if ((.not. icephs) .or. (ttrice.eq.0.0)) go to 10
   trinv = 1.0/ttrice
   do i = 1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition 
! range from ice to water also accounting for change of hlatv with t 
! above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
!
     tc     = t(i) - 273.16
     lflg   = (tc.ge.-ttrice .and. tc.lt.0.0)
     weight = min(-tc*trinv,1.0)
     hlatsb = hlatv + weight*hlatf
     hlatvp = hlatv - 2369.0*tc
     if (t(i).lt.273.16) then
       hltalt = hlatsb
     else
       hltalt = hlatvp
     end if
     if (lflg) then
       tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) +                 &
                 tc*pcf(5))))
     else
       tterm = 0.0
     end if
     desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
     gam(i) = hltalt*qs(i)*p(i)*desdt/                                         &
              (cp*es(i)*(p(i) - omeps*es(i)))
     if(qs(i).eq.1.0) gam(i) = 0.0
   end do
   return
!
! No icephs or water to ice transition
!
10 do i = 1,len
!
! Account for change of hlatv with t above 273.16 where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
     hlatvp = hlatv - 2369.0*(t(i)-273.16)
     if (icephs) then
       hlatsb = hlatv + hlatf
     else
       hlatsb = hlatv
     end if
     if (t(i).lt.273.16) then
       hltalt = hlatsb
     else
       hltalt = hlatvp
     end if
     desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
     gam(i) = hltalt*qs(i)*p(i)*desdt/                                         &
               (cp*es(i)*(p(i) - omeps*es(i)))
     if (qs(i) .eq. 1.0) gam(i) = 0.0
   end do
!
   return
   end subroutine ccm2_saturation_table2
!
!-------------------------------------------------------------------------------
#endif /* CCMCNV end */
   end module phys_ccm_module
