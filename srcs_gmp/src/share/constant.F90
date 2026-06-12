#include "define.h"
module constant
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real, parameter :: rerth_ =6.3712e+6
   real, parameter :: rrerth_=1./rerth_
   real, parameter :: g_     =9.8060e+0
   real, parameter :: omega_ =7.2921e-5
   real, parameter :: rd_    =2.8705e+2
   real, parameter :: rv_    =4.6150e+2
   real, parameter :: rdorv_ =rd_/rv_
   real, parameter :: rdorvm1_ =rd_/rv_-1.
   real, parameter :: rdog_  =rd_/g_
   real, parameter :: rvrdm1_=0.6077338
   real, parameter :: rvordm1_=rv_/rd_-1.
   real, parameter :: fv_    =rv_/rd_-1.
   real, parameter :: cp_    =1.0046e+3
   real, parameter :: cv_    =7.1760e+2
   real, parameter :: cvap_  =1.8460e+3
   real, parameter :: cliq_  =4.1855e+3
   real, parameter :: cice_  =2.1060E+3
   real, parameter :: hvap_  =2.5000e+6
   real, parameter :: hfus_  =3.3358e+5
   real, parameter :: hsub_  =2.8340E+6
   real, parameter :: psat_  =6.1078e+2
   real, parameter :: sbc_   =5.6730e-8
   real, parameter :: solr_  =1.3533e+3
   real, parameter :: t0c_   =2.7315e+2        ! ice/water mix temperature (k)
   real, parameter :: ttp_   =2.7316e+2
   real, parameter :: cal_   =4.1855e+0
   real, parameter :: convrad_ =cal_*1.e4/60.
   real, parameter :: rhoh2o_=1000.            ! water density (kg/m^3)
   real, parameter :: pi_    =3.14159265358979
   real, parameter :: qmin_  =1.0e-30
   real, parameter :: qmin8_ =1.0e-8
   real, parameter :: qmin9_ =1.0e-9
   real, parameter :: qmin30_=1.0e-30
   real, parameter :: aday_  =86400.
   real, parameter :: akapa_ =rd_/cp_
   real, parameter :: elocp_ =hvap_/cp_
   real, parameter :: karman_=0.4
   real, parameter :: rhoair0_  = 1.28
   real, parameter :: rhosnow_  = 100.
   real, parameter :: cb2pa_    = 1000.
   real, parameter :: degrad_=180.0e0/pi_
   real, parameter :: daysec_=1.1574e-5
!
end module constant
!-------------------------------------------------------------------------------
