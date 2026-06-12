! ------------------------ code history --------------------------------
! source file:       vic_surfeb.h
! purpose:           variable common block for surface energy balance
! prgmmr:            ji chen
! org:               ecpc/crd/sio/ucsd 
! date:              november 2003
! 
! notes:
!    upcase "e" before variables means computing energy balance
! ----------------------------------------------------------------------

! ------------------- input variables ----------------------------------

! -- model basic parameter
!
#include <define.h>
!
      integer, parameter :: msl=_lsoil_, nsl=_nsoil_ 
                           ! msl : number of soil layers
                           ! nsl : number of soil nodes for soil temperature
      real :: edtime       ! time step (second)

! -- atmosphere data

      real :: epgcm        ! pressure (pa)
      real :: etair        ! air temperature (k)
      real :: erhoair      ! air density (kg/m3)
      real :: evpair       ! actual vapor pressure of air (pa)
      real :: evpd         ! vapor pressure deficit (pa)
      real :: era          ! aerodynamical resistance (s/m)
      real :: ewind        ! wind speed (m/s)
      real :: enetsht      ! net shortwave radiation (w/m2)
      real :: elongw       ! longwave radiation (w/m2)

! -- soil parameters

      real :: ebinf        ! vic infiltration parameter (n/a)
      real :: edph(msl)    ! soil moisture layer thickness (m)
      real :: ewcr(msl)    ! ~70% of field capacity (mm)
      real :: ewpw(msl)    ! wilting point soil moisture (mm)
      real :: esmr(msl)    ! residual moisture (mm)
      real :: esld(msl)    ! soil particle density (kg/m^3)
      real :: ebkd(msl)    ! bulk density of soil layer (kg/m^3)
      real :: eqrt(msl)    ! quartz content of soil (fraction)

      real :: edphnd(nsl)  ! soil thermal node thicknesses (m)
      real :: ealpnd(nsl)  ! coef for computing soiltemp profile
      real :: ebetnd(nsl)  ! coef for computing soiltemp profile
      real :: egamnd(nsl)  ! coef for computing soiltemp profile
      real :: esmxnd(nsl)  ! maximum soil moisture (m3/m3)
      real :: ebubnd(nsl)  ! bubbling pressure of soil (cm)
      real :: eexpnd(nsl)  ! parameter for ksat with soil moisture (n/a)

      real :: ekapnd(nsl)  ! soil layer thermal conductivity (w/m/k)
      real :: ecsnd(nsl)   ! soil layer heat capacity (j/m^3/k)
      real :: esmnd(nsl)   ! soil moisture at node (m3/m3)
      real :: eicend(nsl)  ! soil ice at node (m3/m3)

! -- vegetation parameters

      integer :: evtype    ! vegetation type
      real :: efflai       ! leaf area index (fraction)
      real :: efrt(msl)    ! root content (fraction)

! -- land surface variables

      real :: esncvfr      ! snow cover fraction (fraction)
      real :: efcsn        ! snow on vegetation (m)
      real :: eftsf        ! snow surface temperature (k)
      real :: efswq        ! snow water equivalent (m)
      real :: efrsn        ! snow density (kg/m3)
      real :: esnowd       ! snow depth (m)

      real :: efsic(msl)   ! soil ice (mm)

! ----------------------------------------------------------------------

      common /vicsurfeb_inp/edtime,   epgcm,   etair, erhoair,                 &
                            evpair,    evpd,     era,   ewind, enetsht,        &
                            elongw,   ebinf,    edph,    ewcr,    ewpw,        &
                              esmr,    esld,    ebkd,    eqrt,  edphnd,        &
                            ealpnd,  ebetnd,  egamnd,  esmxnd,  ebubnd,        &
                            eexpnd,  ekapnd,   ecsnd,   esmnd,  eicend,        &
                            evtype,  efflai,    efrt, esncvfr,   efcsn,        &
                             eftsf,   efswq,   efrsn,  esnowd,   efsic

! ------------------- modified variables -------------------------------
! the following for variables will be modified
! the inputs start with e..., and the modified ones start with m...
!
      real :: efcwt        ! dew and rain trapped on vegetation (m)
      real :: efrtf        ! rain throughfall (m/time step)
      real :: efsmc(msl)   ! soil moisture (liquid + ice) (mm)
      real :: eftnd(nsl)   ! soil temperature profiles (k)
      real :: mfcwt        ! modified dew and rain trapped on vegetation (m)
      real :: mfrtf        ! modified rain throughfall (m/time step)
      real :: mfsmc(msl)   ! modified soil moisture (liquid + ice) (mm)
      real :: mftnd(nsl)   ! modified soil temperature profiles (k)
! ----------------------------------------------------------------------

      common /vicsurfeb_mod/efcwt, efrtf, efsmc, eftnd,                        &
                            mfcwt, mfrtf, mfsmc, mftnd

! ------------------- output variables ---------------------------------
      real :: olath        ! latent heat flux (w/m2)
      real :: osenh        ! sensible heat flux (w/m2)
      real :: ogrndh       ! ground heat flux (w/m2)
      real :: odelth       ! ground heat storage (w/m2)
!-----------------------------------------------------------------------

      common /vicsurfeb_out/olath, osenh, ogrndh, odelth
