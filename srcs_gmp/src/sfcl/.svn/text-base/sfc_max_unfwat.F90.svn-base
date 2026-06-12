#include <define.h>
   function sfc_max_unfwat(t, soilmx, bubble, expt)
!-------------------------------------------------------------------------------
!
! subprogram:  sfc_max_unfwat
! abstract: compute the maximum amount of unfrozen 
!   water that can exist at the current temperature.
!
! program history log:
!   2003-03-09  ji chen                cvs verion setup
!   2008-03-09  kyeng hee seol         debugging
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist in
!
! notes:
!    please check equnation (14) in cherkauer and lettenmaier 
!    (jgr 1999) for the equation of computing unfrozen soil moisture
!
!-------------------------------------------------------------------------------
   use constant, only : g=>g_
#include <vartyp.h>
!-------------------------------------------------------------------------------
!
!  input variables
!
   real                 ::  t      ! soil temperature (k)
   real                 ::  soilmx ! maximum soil moisture (mm or m3/m3)
   real                 ::  bubble ! bubbling pressure of soil (cm)
   real                 ::  expt   ! parameter for ksat with soil moisture (n/a)
!
!  output variables 
!
   real                 ::  sfc_max_unfwat ! unfrozen soil moisture (mm or m3/m3)
!
!  local variables 
!
   real, parameter      ::  lf = 3.337e5   ! latent heat of freezing (j/kg) at 0c
   real                 ::  unfrozen, bas
!
   bas = (-lf * (t-273.15)) / t  / (g * bubble / 100.)
   unfrozen = soilmx * bas **( -2.0/(expt-3.0))
!
   if(unfrozen .gt. soilmx) unfrozen = soilmx
   if(unfrozen .lt. 0) unfrozen = 0
!   
   sfc_max_unfwat = unfrozen
!
   return
   end function sfc_max_unfwat
!-------------------------------------------------------------------------------
