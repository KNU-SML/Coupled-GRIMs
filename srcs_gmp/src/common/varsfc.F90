#include "define.h"
   module varsfc
!-------------------------------------------------------------------------------
!
!  varsfc.h
!
!  define sfc file variables and their properties.
!
!  this subroutine should list sfc file variables for
!  each land surface program
!  
!  if you are introducing new land surface model, you need to 
!  introduce new numsfcs and numsfcv and add subsequent 3 parameters 
!  for each new field.
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer,parameter    ::  lsoil_=_lsoil_
   integer,parameter    ::  nsoil_=_nsoil_
   integer,parameter    ::  msub_=_msub_
   integer,parameter    ::  lalbd_=_lalbd_
   integer,parameter    ::  kslmb_=lsoil_*msub_
   integer,parameter    ::  nslmb_=nsoil_*msub_
   integer,parameter    ::  naer_=5
   integer,parameter    ::  nden_=2
!
!  1. pointer to /comsfc/
!
#ifdef OSULSM1
!  osulsm1
!
   integer, parameter   ::  numsfcs=12+lsoil_*2+lalbd_+nden_+3*naer_
   integer, parameter   ::  numsfcv=19
   integer, parameter   ::  jtsf=     1       ! sst and land surface temperature
   integer, parameter   ::  jsmc=jtsf+1       ! soil moisture
   integer, parameter   ::  jsno=jsmc+lsoil_  ! snow depth
   integer, parameter   ::  jstc=jsno+1       ! soil temperature
   integer, parameter   ::  jtg3=jstc+lsoil_  ! deep soil temperature
   integer, parameter   ::  jzor=jtg3+1       ! surface roughness
   integer, parameter   ::  jcv =jzor+1       ! convective cloud cover
   integer, parameter   ::  jcvb=jcv +1       ! convective cloud base
   integer, parameter   ::  jcvt=jcvb+1       ! convective cloud top
   integer, parameter   ::  jalb=jcvt+1       ! albedo
   integer, parameter   ::  jsli=jalb+lalbd_  ! land/sea ice/snow mask
   integer, parameter   ::  jplr=jsli+1       ! stomato resistance
   integer, parameter   ::  jcpy=jplr+1       ! canopy
   integer, parameter   ::  jf10=jcpy+1       ! 10m and sig1 conv factor
   integer, parameter   ::  jvet=jf10+1       ! vegetation type
   integer, parameter   ::  jsot=jvet+1       ! soil type
   integer, parameter   ::  jalf=jsot+1       ! albedo fraction
   integer, parameter   ::  just=jalf+2       ! fricational velocity
   integer, parameter   ::  jffm=just+1       ! coefficient for momentun exchange
   integer, parameter   ::  jffh=jffm+1       ! coefficient for heat exchange
   integer, parameter   ::  joml=jffh+1       ! ocean mixed layer depth
   integer, parameter   ::  jaer=joml+1       ! aerosol distribution
   integer, parameter   ::  jkpr=jaer+naer_   ! aerosol distribution
   integer, parameter   ::  jden=jkpr+1       ! aerosol distribution
   integer, parameter   ::  jdxc=jden+nden_   ! aerosol distribution
   integer, parameter   ::  jmix=jdxc+naer_   ! aerosol distribution
#endif
#ifdef OSULSM2
!
!  osulsm2
!
   integer, parameter   ::  numsfcs=18+lsoil_*2+lalbd_+2+3*naer_+nden_
   integer, parameter   ::  numsfcv=26
   integer, parameter   ::  jtsf=     1       ! sst and land surface temperature
   integer, parameter   ::  jsmc=jtsf+1       ! soil moisture
   integer, parameter   ::  jsno=jsmc+lsoil_  ! snow depth
   integer, parameter   ::  jstc=jsno+1       ! soil temperature
   integer, parameter   ::  jtg3=jstc+lsoil_  ! deep soil temperature
   integer, parameter   ::  jzor=jtg3+1       ! surface roughness
   integer, parameter   ::  jcv =jzor+1       ! convective cloud cover
   integer, parameter   ::  jcvb=jcv +1       ! convective cloud base
   integer, parameter   ::  jcvt=jcvb+1       ! convective cloud top
   integer, parameter   ::  jalb=jcvt+1       ! albedo
   integer, parameter   ::  jsli=jalb+lalbd_  ! land/sea ice/snow mask
   integer, parameter   ::  jveg=jsli+1       ! vegetation cover
   integer, parameter   ::  jcpy=jveg+1       ! canopy water content
   integer, parameter   ::  jf10=jcpy+1       ! 10m factor
   integer, parameter   ::  jvet=jf10+1       ! vegetation type
   integer, parameter   ::  jsot=jvet+1       ! soil type
   integer, parameter   ::  jalf=jsot+1       ! albedo fraction
   integer, parameter   ::  just=jalf+2       ! fricational velocity
   integer, parameter   ::  jffm=just+1       ! coefficient for momentun exchange
   integer, parameter   ::  jffh=jffm+1       ! coefficient for heat exchange
   integer, parameter   ::  joml=jffh+1       ! ocean mixed layer depth
   integer, parameter   ::  jaer=joml+1       ! aerosol distribution
   integer, parameter   ::  jkpr=jaer+naer_   ! aerosol distribution
   integer, parameter   ::  jden=jkpr+1       ! aerosol distribution
   integer, parameter   ::  jdxc=jden+nden_   ! aerosol distribution
   integer, parameter   ::  jmix=jdxc+naer_   ! aerosol distribution
#endif
#ifdef NOALSM1
!
!  noalsm1
!
   integer, parameter   ::  numsfcs=25+lsoil_*3+lalbd_+2+3*naer_+nden_
   integer, parameter   ::  numsfcv=34
   integer, parameter   ::  jtsf=     1       ! sst and land surface temperature
   integer, parameter   ::  jsmc=jtsf+1       ! soil moisture
   integer, parameter   ::  jsno=jsmc+lsoil_  ! snow depth
   integer, parameter   ::  jstc=jsno+1       ! soil temperature
   integer, parameter   ::  jtg3=jstc+lsoil_  ! deep soil temperature
   integer, parameter   ::  jzor=jtg3+1       ! surface roughness
   integer, parameter   ::  jcv =jzor+1       ! convective cloud cover
   integer, parameter   ::  jcvb=jcv +1       ! convective cloud base
   integer, parameter   ::  jcvt=jcvb+1       ! convective cloud top
   integer, parameter   ::  jalb=jcvt+1       ! albedo
   integer, parameter   ::  jsli=jalb+lalbd_  ! land/sea ice/snow mask
   integer, parameter   ::  jveg=jsli+1       ! vegetation cover
   integer, parameter   ::  jcpy=jveg+1       ! canopy water content
   integer, parameter   ::  jf10=jcpy+1       ! 10m and sig1 conv factor
   integer, parameter   ::  jvet=jf10+1       ! vegetation type
   integer, parameter   ::  jsot=jvet+1       ! soil type
   integer, parameter   ::  jalf=jsot+1       ! albedo fraction
   integer, parameter   ::  just=jalf+2       ! fricational velocity
   integer, parameter   ::  jffm=just+1       ! coefficient for momentun exchange
   integer, parameter   ::  jffh=jffm+1       ! coefficient for heat exchange
   integer, parameter   ::  jprc=jffh+1       ! precip for noah
   integer, parameter   ::  jsrf=jprc+1       ! sr flag for noah
   integer, parameter   ::  jsnd=jsrf+1       ! snow depth for noah
   integer, parameter   ::  jslc=jsnd+1       ! slc for noah
   integer, parameter   ::  jsmn=jslc+lsoil_  ! min vegetation cover for noah
   integer, parameter   ::  jsmx=jsmn+1       ! max vegetation cover for noah
   integer, parameter   ::  jslo=jsmx+1       ! slope type for noah
   integer, parameter   ::  jsna=jslo+1       ! snow albedo for noah
   integer, parameter   ::  joml=jsna+1       ! ocean mixed layer depth
   integer, parameter   ::  jaer=joml+1       ! aerosol distribution
   integer, parameter   ::  jkpr=jaer+naer_   ! aerosol distribution
   integer, parameter   ::  jden=jkpr+1       ! aerosol distribution
   integer, parameter   ::  jdxc=jden+nden_   ! aerosol distribution
   integer, parameter   ::  jmix=jdxc+naer_   ! aerosol distribution
#endif
#ifdef VICLSM1
!
!  viclsm1
!
   integer, parameter   ::  numsfcs=34+lsoil_*14+nsoil_*8+lalbd_+2+3*naer_+nden_
   integer, parameter   ::  numsfcv=62
   integer, parameter   ::  jtsf=     1       ! sst and land surface temperature
   integer, parameter   ::  jsmc=jtsf+1       ! soil moisture
   integer, parameter   ::  jsno=jsmc+lsoil_  ! snow depth
   integer, parameter   ::  jstc=jsno+1       ! soil temperature
   integer, parameter   ::  jtg3=jstc+nsoil_  ! deep soil temperature
   integer, parameter   ::  jzor=jtg3+1       ! surface roughness
   integer, parameter   ::  jcv =jzor+1       ! convective cloud cover
   integer, parameter   ::  jcvb=jcv +1       ! convective cloud base
   integer, parameter   ::  jcvt=jcvb+1       ! convective cloud top
   integer, parameter   ::  jalb=jcvt+1       ! albedo
   integer, parameter   ::  jsli=jalb+lalbd_  ! land/sea ice/snow mask
   integer, parameter   ::  jveg=jsli+1       ! vegetation cover
   integer, parameter   ::  jcpy=jveg+1       ! canopy water content
   integer, parameter   ::  jf10=jcpy+1       ! 10m and sig1 conv factor
   integer, parameter   ::  jvet=jf10+1       ! vegetation type
   integer, parameter   ::  jrot=jvet+1       ! root content
   integer, parameter   ::  jalf=jrot+lsoil_  ! albedo fraction
   integer, parameter   ::  just=jalf+2       ! fricational velocity
   integer, parameter   ::  jffm=just+1       ! coefficient for momentun exchange
   integer, parameter   ::  jffh=jffm+1       ! coefficient for heat exchange
   integer, parameter   ::  jprc=jffh+1       ! precip for vic
   integer, parameter   ::  jsrf=jprc+1       ! sr flag for vic
   integer, parameter   ::  jbif=jsrf+1       ! Variable infil curve parameter (N/A)
   integer, parameter   ::  jds =jbif+1       ! Fract of Dsm nonlinear baseflow begins
   integer, parameter   ::  jdsm=jds +1       ! Maximum velocity of baseflow (mm/day)
   integer, parameter   ::  jws =jdsm+1       ! Fract maxi sm nonlinear baseflow occurs
   integer, parameter   ::  jcef=jws +1       ! c
   integer, parameter   ::  jexp=jcef+1       ! Para the vari of Ksat with sm
   integer, parameter   ::  jkst=jexp+lsoil_  ! Saturated hydrologic conductivity (mm/day)
   integer, parameter   ::  jdph=jkst+lsoil_  ! thickness of soil layer (m)
   integer, parameter   ::  jbub=jdph+lsoil_  ! Bubbling pressure of soil layer (cm)
   integer, parameter   ::  jqrt=jbub+lsoil_  ! Quartz content of soil layer (fraction)
   integer, parameter   ::  jbkd=jqrt+lsoil_  ! Bulk density of soil layer (kg/m3)
   integer, parameter   ::  jsld=jbkd+lsoil_  ! Soil density of soil layer (kg/m3)
   integer, parameter   ::  jwcr=jsld+lsoil_  ! sm content at the critical point (mm)
   integer, parameter   ::  jwpw=jwcr+lsoil_  ! sm content wilting point (mm)
   integer, parameter   ::  jsmr=jwpw+lsoil_  ! Soil moisture residual moisture (mm)
   integer, parameter   ::  jsmx=jsmr+lsoil_  ! maximum soil moisture (mm)
   integer, parameter   ::  jdpn=jsmx+lsoil_  ! thickness of soil node (m)
   integer, parameter   ::  jsxn=jdpn+nsoil_  ! maximum sm at soil node (m3/m3)
   integer, parameter   ::  jepn=jsxn+nsoil_  ! Para the vari of Ksat at soil node (N/A)
   integer, parameter   ::  jbbn=jepn+nsoil_  ! bubbling pressure at soil node (cm)
   integer, parameter   ::  japn=jbbn+nsoil_  ! para alpha at soil node
   integer, parameter   ::  jbtn=japn+nsoil_  ! para beta at soil node
   integer, parameter   ::  jgmn=jbtn+nsoil_  ! para gamma at soil node
   integer, parameter   ::  jlai=jgmn+nsoil_  ! leaf area index
   integer, parameter   ::  jslz=jlai+1       ! surface roughness of bare soil (m)
   integer, parameter   ::  jsnz=jslz+1       ! surface roughness of snow pack (m)
   integer, parameter   ::  jsic=jsnz+1       ! soil ice content (mm)
   integer, parameter   ::  jcsn=jsic+lsoil_  ! canopy snow  (mm h2o)
   integer, parameter   ::  jrsn=jcsn+1       ! snow density (kg/m^3)
   integer, parameter   ::  jtsn=jrsn+1       ! snow surface temperature (K)
   integer, parameter   ::  jtpk=jtsn+1       ! snow pack temperature (K)
   integer, parameter   ::  jsfw=jtpk+1       ! surface snow water equivalent (mm h2o)
   integer, parameter   ::  jpkw=jsfw+1       ! snow pack snow water equivalent (mm h2o)
   integer, parameter   ::  jlst=jpkw+1       ! time step since last snow fall
   integer, parameter   ::  joml=jlst+1       ! ocean mixed layer depth
   integer, parameter   ::  jaer=joml+1       ! aerosol distribution
   integer, parameter   ::  jkpr=jaer+naer_   ! aerosol distribution
   integer, parameter   ::  jden=jkpr+1       ! aerosol distribution
   integer, parameter   ::  jdxc=jden+nden_   ! aerosol distribution
   integer, parameter   ::  jmix=jdxc+naer_   ! aerosol distribution
#endif
#ifdef VICLSM2
!
!  viclsm2
!
   integer, parameter   ::  kslmb=lsoil_*msub_, nslmb=nsoil_*msub_
   integer, parameter   ::  numsfcs=27+lsoil_*12+nsoil_*8+msub_*11+            & 
                            kslmb*4+nslmb+lalbd_+2+3*naer_+nden_
   integer, parameter   ::  numsfcv=69
   integer, parameter   ::  jtsf=     1       ! sst and land surface temperature
   integer, parameter   ::  jsmc=jtsf+1       ! soil moisture
   integer, parameter   ::  jsno=jsmc+lsoil_  ! snow depth
   integer, parameter   ::  jstc=jsno+1       ! soil temperature
   integer, parameter   ::  jtg3=jstc+nsoil_  ! deep soil temperature
   integer, parameter   ::  jzor=jtg3+1       ! surface roughness
   integer, parameter   ::  jcv =jzor+1       ! convective cloud cover
   integer, parameter   ::  jcvb=jcv +1       ! convective cloud base
   integer, parameter   ::  jcvt=jcvb+1       ! convective cloud top
   integer, parameter   ::  jalb=jcvt+1       ! albedo
   integer, parameter   ::  jsli=jalb+lalbd_  ! land/sea ice/snow mask
   integer, parameter   ::  jveg=jsli+1       ! fraction of vegetation cover
   integer, parameter   ::  jcpy=jveg+1       ! canopy water content
   integer, parameter   ::  jf10=jcpy+1       ! 10m and sig1 conv factor
   integer, parameter   ::  jvet=jf10+1       ! vegetation type
   integer, parameter   ::  jrot=jvet+1       ! root content
   integer, parameter   ::  jalf=jrot+kslmb   ! albedo fraction
   integer, parameter   ::  just=jalf+2       ! fricational velocity
   integer, parameter   ::  jffm=just+1       ! coefficient for momentun exchange
   integer, parameter   ::  jffh=jffm+1       ! coefficient for heat exchange
   integer, parameter   ::  jprc=jffh+1       ! precip for vic
   integer, parameter   ::  jsrf=jprc+1       ! sr flag for vic
   integer, parameter   ::  jbif=jsrf+1       ! Variable infil curve parameter (N/A)
   integer, parameter   ::  jds =jbif+1       ! Fract of Dsm nonlinear baseflow begins
   integer, parameter   ::  jdsm=jds +1       ! Maximum velocity of baseflow (mm/day)
   integer, parameter   ::  jws =jdsm+1       ! Fract maxi sm nonlinear baseflow occurs
   integer, parameter   ::  jcef=jws +1       ! c
   integer, parameter   ::  jexp=jcef+1       ! Para the vari of Ksat with sm
   integer, parameter   ::  jkst=jexp+lsoil_  ! Saturated hydrologic conductivity (mm/day)
   integer, parameter   ::  jdph=jkst+lsoil_  ! thickness of soil layer (m)
   integer, parameter   ::  jbub=jdph+lsoil_  ! Bubbling pressure of soil layer (cm)
   integer, parameter   ::  jqrt=jbub+lsoil_  ! Quartz content of soil layer (fraction)
   integer, parameter   ::  jbkd=jqrt+lsoil_  ! Bulk density of soil layer (kg/m3)
   integer, parameter   ::  jsld=jbkd+lsoil_  ! Soil density of soil layer (kg/m3)
   integer, parameter   ::  jwcr=jsld+lsoil_  ! sm content at the critical point (mm)
   integer, parameter   ::  jwpw=jwcr+lsoil_  ! sm content wilting point (mm)
   integer, parameter   ::  jsmr=jwpw+lsoil_  ! Soil moisture residual moisture (mm)
   integer, parameter   ::  jsmx=jsmr+lsoil_  ! maximum soil moisture (mm)
   integer, parameter   ::  jdpn=jsmx+lsoil_  ! thickness of soil node (m)
   integer, parameter   ::  jsxn=jdpn+nsoil_  ! maximum sm at soil node (m3/m3)
   integer, parameter   ::  jepn=jsxn+nsoil_  ! Para the vari of Ksat at soil node (N/A)
   integer, parameter   ::  jbbn=jepn+nsoil_  ! bubbling pressure at soil node (cm)
   integer, parameter   ::  japn=jbbn+nsoil_  ! para alpha at soil node
   integer, parameter   ::  jbtn=japn+nsoil_  ! para beta at soil node
   integer, parameter   ::  jgmn=jbtn+nsoil_  ! para gamma at soil node
   integer, parameter   ::  jslz=jgmn+nsoil_  ! surface roughness of bare soil (m)
   integer, parameter   ::  jsnz=jslz+1       ! surface roughness of snow pack (m)
   integer, parameter   ::  jnve=jsnz+1       ! number of tiles in a grid (N/A)
   integer, parameter   ::  jmfr=jnve+1       ! vegetation fraction for tiles
   integer, parameter   ::  jmcp=jmfr+msub_   ! canopy water for tiles
   integer, parameter   ::  jmvt=jmcp+msub_   ! vegetation type for tiles
   integer, parameter   ::  jlai=jmvt+msub_   ! leaf area index
   integer, parameter   ::  jmsn=jlai+msub_   ! snow depth for tiles
   integer, parameter   ::  jmsm=jmsn+msub_   ! snow moisture for tiles (mm)
   integer, parameter   ::  jmsi=jmsm+kslmb   ! soil ice content for tiles (mm)
   integer, parameter   ::  jmst=jmsi+kslmb   ! soil temperature for tiles (K)
   integer, parameter   ::  jcsn=jmst+nslmb   ! canopy snow  (mm h2o)
   integer, parameter   ::  jrsn=jcsn+msub_   ! snow density (kg/m^3)
   integer, parameter   ::  jtsn=jrsn+msub_   ! snow surface temperature (K)
   integer, parameter   ::  jtpk=jtsn+msub_   ! snow pack temperature (K)
   integer, parameter   ::  jsfw=jtpk+msub_   ! surface snow water equivalent (mm h2o)
   integer, parameter   ::  jpkw=jsfw+msub_   ! snow pack snow water equivalent (mm h2o)
   integer, parameter   ::  jlst=jpkw+msub_   ! time step since last snow fall
   integer, parameter   ::  joml=jlst+1       ! ocean mixed layer depth
   integer, parameter   ::  jaer=joml+1       ! aerosol distribution
   integer, parameter   ::  jkpr=jaer+naer_   ! aerosol distribution
   integer, parameter   ::  jden=jkpr+1       ! aerosol distribution
   integer, parameter   ::  jdxc=jden+nden_   ! aerosol distribution
   integer, parameter   ::  jmix=jdxc+naer_   ! aerosol distribution
#endif
!
!  2. number of layers
!
!  define sfc file variables and their properties.
   integer              ::  ksfc(numsfcs)
!
#ifdef OSULSM1
! osulsm1
!
   data ksfc(jtsf)/1/
   data ksfc(jsmc)/lsoil_/
   data ksfc(jsno)/1/
   data ksfc(jstc)/lsoil_/
   data ksfc(jtg3)/1/
   data ksfc(jzor)/1/
   data ksfc(jcv )/1/
   data ksfc(jcvb)/1/
   data ksfc(jcvt)/1/
   data ksfc(jalb)/lalbd_/
   data ksfc(jsli)/1/
   data ksfc(jplr)/1/
   data ksfc(jcpy)/1/
   data ksfc(jf10)/1/
#endif
#ifdef OSULSM2
!
! osulsm2
!
   data ksfc(jtsf)/1/
   data ksfc(jsmc)/lsoil_/
   data ksfc(jsno)/1/
   data ksfc(jstc)/lsoil_/
   data ksfc(jtg3)/1/
   data ksfc(jzor)/1/
   data ksfc(jcv )/1/
   data ksfc(jcvb)/1/
   data ksfc(jcvt)/1/
   data ksfc(jalb)/lalbd_/
   data ksfc(jsli)/1/
   data ksfc(jveg)/1/
   data ksfc(jcpy)/1/
   data ksfc(jf10)/1/
   data ksfc(jvet)/1/
   data ksfc(jsot)/1/
   data ksfc(jalf)/2/
   data ksfc(just)/1/
   data ksfc(jffm)/1/
   data ksfc(jffh)/1/
   data ksfc(joml)/1/
   data ksfc(jaer)/naer_/
   data ksfc(jkpr)/1/
   data ksfc(jden)/nden_/
   data ksfc(jdxc)/naer_/
   data ksfc(jmix)/naer_/
#endif
#ifdef NOALSM1
!
! noalsm1
!
   data ksfc(jtsf)/1/
   data ksfc(jsmc)/lsoil_/
   data ksfc(jsno)/1/
   data ksfc(jstc)/lsoil_/
   data ksfc(jtg3)/1/
   data ksfc(jzor)/1/
   data ksfc(jcv )/1/
   data ksfc(jcvb)/1/
   data ksfc(jcvt)/1/
   data ksfc(jalb)/lalbd_/
   data ksfc(jsli)/1/
   data ksfc(jveg)/1/
   data ksfc(jcpy)/1/
   data ksfc(jf10)/1/
   data ksfc(jvet)/1/
   data ksfc(jsot)/1/
   data ksfc(jalf)/2/
   data ksfc(just)/1/
   data ksfc(jffm)/1/
   data ksfc(jffh)/1/
   data ksfc(jprc)/1/
   data ksfc(jsrf)/1/
   data ksfc(jsnd)/1/
   data ksfc(jslc)/lsoil_/
   data ksfc(jsmn)/1/
   data ksfc(jsmx)/1/
   data ksfc(jslo)/1/
   data ksfc(jsna)/1/
   data ksfc(joml)/1/
   data ksfc(jaer)/naer_/
   data ksfc(jkpr)/1/
   data ksfc(jden)/nden_/
   data ksfc(jdxc)/naer_/
   data ksfc(jmix)/naer_/
#endif
#ifdef VICLSM1
!
! viclsm1
!
   data ksfc(jtsf)/1/
   data ksfc(jsmc)/lsoil_/
   data ksfc(jsno)/1/
   data ksfc(jstc)/nsoil_/
   data ksfc(jtg3)/1/
   data ksfc(jzor)/1/
   data ksfc(jcv )/1/
   data ksfc(jcvb)/1/
   data ksfc(jcvt)/1/
   data ksfc(jalb)/lalbd_/
   data ksfc(jsli)/1/
   data ksfc(jveg)/1/
   data ksfc(jcpy)/1/
   data ksfc(jf10)/1/
   data ksfc(jvet)/1/
   data ksfc(jrot)/lsoil_/
   data ksfc(jalf)/2/
   data ksfc(just)/1/
   data ksfc(jffm)/1/
   data ksfc(jffh)/1/
   data ksfc(jprc)/1/
   data ksfc(jsrf)/1/
   data ksfc(jbif)/1/
   data ksfc(jds )/1/
   data ksfc(jdsm)/1/
   data ksfc(jws )/1/
   data ksfc(jcef)/1/
   data ksfc(jexp)/lsoil_/
   data ksfc(jkst)/lsoil_/
   data ksfc(jdph)/lsoil_/
   data ksfc(jbub)/lsoil_/
   data ksfc(jqrt)/lsoil_/
   data ksfc(jbkd)/lsoil_/
   data ksfc(jsld)/lsoil_/
   data ksfc(jwcr)/lsoil_/
   data ksfc(jwpw)/lsoil_/
   data ksfc(jsmr)/lsoil_/
   data ksfc(jsmx)/lsoil_/
   data ksfc(jdpn)/nsoil_/
   data ksfc(jsxn)/nsoil_/
   data ksfc(jepn)/nsoil_/
   data ksfc(jbbn)/nsoil_/
   data ksfc(japn)/nsoil_/
   data ksfc(jbtn)/nsoil_/
   data ksfc(jgmn)/nsoil_/
   data ksfc(jlai)/1/
   data ksfc(jslz)/1/
   data ksfc(jsnz)/1/
   data ksfc(jsic)/lsoil_/
   data ksfc(jcsn)/1/
   data ksfc(jrsn)/1/
   data ksfc(jtsn)/1/
   data ksfc(jtpk)/1/
   data ksfc(jsfw)/1/
   data ksfc(jpkw)/1/
   data ksfc(jlst)/1/
   data ksfc(joml)/1/
   data ksfc(jaer)/naer_/
   data ksfc(jkpr)/1/
   data ksfc(jden)/nden_/
   data ksfc(jdxc)/naer_/
   data ksfc(jmix)/naer_/
#endif
#ifdef VICLSM2
!
! viclsm2
!
   data ksfc(jtsf)/1/
   data ksfc(jsmc)/lsoil_/
   data ksfc(jsno)/1/
   data ksfc(jstc)/nsoil_/
   data ksfc(jtg3)/1/
   data ksfc(jzor)/1/
   data ksfc(jcv )/1/
   data ksfc(jcvb)/1/
   data ksfc(jcvt)/1/
   data ksfc(jalb)/lalbd_/
   data ksfc(jsli)/1/
   data ksfc(jveg)/1/
   data ksfc(jcpy)/1/
   data ksfc(jf10)/1/
   data ksfc(jvet)/1/
   data ksfc(jrot)/kslmb /
   data ksfc(jalf)/2/
   data ksfc(just)/1/
   data ksfc(jffm)/1/
   data ksfc(jffh)/1/
   data ksfc(jprc)/1/
   data ksfc(jsrf)/1/
   data ksfc(jbif)/1/
   data ksfc(jds )/1/
   data ksfc(jdsm)/1/
   data ksfc(jws )/1/
   data ksfc(jcef)/1/
   data ksfc(jexp)/lsoil_/
   data ksfc(jkst)/lsoil_/
   data ksfc(jdph)/lsoil_/
   data ksfc(jbub)/lsoil_/
   data ksfc(jqrt)/lsoil_/
   data ksfc(jbkd)/lsoil_/
   data ksfc(jsld)/lsoil_/
   data ksfc(jwcr)/lsoil_/
   data ksfc(jwpw)/lsoil_/
   data ksfc(jsmr)/lsoil_/
   data ksfc(jsmx)/lsoil_/
   data ksfc(jdpn)/nsoil_/
   data ksfc(jsxn)/nsoil_/
   data ksfc(jepn)/nsoil_/
   data ksfc(jbbn)/nsoil_/
   data ksfc(japn)/nsoil_/
   data ksfc(jbtn)/nsoil_/
   data ksfc(jgmn)/nsoil_/
   data ksfc(jslz)/1/
   data ksfc(jsnz)/1/
   data ksfc(jnve)/1/
   data ksfc(jmfr)/msub_/
   data ksfc(jmcp)/msub_/
   data ksfc(jmvt)/msub_/
   data ksfc(jlai)/msub_/
   data ksfc(jmsn)/msub_/
   data ksfc(jmsm)/kslmb/
   data ksfc(jmsi)/kslmb/
   data ksfc(jmst)/nslmb/
   data ksfc(jcsn)/msub_/
   data ksfc(jrsn)/msub_/
   data ksfc(jtsn)/msub_/
   data ksfc(jtpk)/msub_/
   data ksfc(jsfw)/msub_/
   data ksfc(jpkw)/msub_/
   data ksfc(jlst)/msub_/
   data ksfc(joml)/1/
   data ksfc(jaer)/naer_/
   data ksfc(jkpr)/1/
   data ksfc(jden)/nden_/
   data ksfc(jdxc)/naer_/
   data ksfc(jmix)/naer_/
#endif
!
!  define sfc file variables and their properties.
!
!  3.  qc parameters (vmaxmin)
!    max over ocean   without sea ice, min over ocean   without sea ice
!    max over  land   without snow   , min over land    without snow
!    max over sea ice without snow   , min over sea ice without snow
!    max over  land   with    snow   , min over land    with    snow
!    max over sea ice with    snow   , min over sea ice with    snow
!       does not check if vmaxmin(1,1,:).lt.vmaxmin(2,1,:)
!       these checks are applied to the value after scaling
!
#ifdef OSULSM1
! osulsm1
!
   real                 ::  vmaxmin(2,5,numsfcs)
   integer              ::  isf,jsf
   data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)     & ! sst and land temp
         / 313.00, 271.21,                                                     &
           353.00, 173.00,                                                     & 
           271.21, 173.00,                                                     &
           273.16, 173.00,                                                     &
           273.16, 173.00/
   data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)     & ! soilm 
         /   0.55,    0.0,                                                     &
             0.55,    0.0,                                                     &
             0.55,    0.0,                                                     &
             0.55,    0.0,                                                     &
             0.55,    0.0/
   data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)     & ! snow depth
         /     0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
           55000.,     0.01,                                                   &
           10000.,     0.01/
   data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)     & ! soilt 
         / 313.00, 200.00,                                                     &
           353.00, 173.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)     & ! deep soil temp
         / 310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)     & ! surface roughness
         /   1.00, 1.e-05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05/
   data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)     & ! convective cloud cover
         /    -1.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               1.,     1.,                                                     &
               1.,     1./
   data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)     & ! 1-type albedo
         /   0.06,   0.06,                                                     &
             0.80,   0.06,                                                     &
             0.80,   0.80,                                                     &
             0.80,   0.06,                                                     &
             0.80,   0.80/
   data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)     & ! land-sea-seaice mask
         /     0.,     0.,                                                     &
               1.,     1.,                                                     &
               2.,     2.,                                                     &
               1.,     1.,                                                     &
               2.,     2./
   data ((vmaxmin(isf,jsf,jplr),isf=1,2),jsf=1,5)     & ! plant resis
         / 1000.0,    0.0,                                                     &
           1000.0,    0.0,                                                     &
           1000.0,    0.0,                                                     &
           1000.0,    0.0,                                                     &
           1000.0,    0.0/
   data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)     & ! canopy water
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)     & ! 10m and sig1 conv factor
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
#endif
#ifdef OSULSM2
!
! osulsm2
!
   real                 ::  vmaxmin(2,5,numsfcs)
   integer              ::  isf,jsf
   data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)     & ! sst and land temp
         / 313.00, 271.21,                                                     &
           353.00, 173.00,                                                     &
           271.21, 173.00,                                                     &
           273.16, 173.00,                                                     &
           273.16, 173.00/
   data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)     & ! soilm 
         /   0.55,   0.00,                                                     &
             0.55,   0.00,                                                     &
             0.55,   0.00,                                                     &
             0.55,   0.00,                                                     &
             0.55,   0.00/
   data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)     & ! snow depth
         /     0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
           55000.,     0.01,                                                   &
           10000.,     0.01/
   data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)     & ! soilt 
         / 313.00, 200.00,                                                     &
           353.00, 173.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)     & ! deep soil temp
         / 310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)     & ! surface roughness
         /   1.00, 1.e-05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05/
   data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)     & ! convective cloud cover
         /     0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               1.,     1.,                                                     &
               1.,     1./
   data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)     & ! 4-type albedo
         /   0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01/
   data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)     & ! land-sea-seaice mask
         /     0.,     0.,                                                     &
               1.,     1.,                                                     &
               2.,     2.,                                                     &
               1.,     1.,                                                     &
               2.,     2./
   data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)     & ! vegetation cover
         /     0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)     & ! canopy water
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)     & ! 10m and sig1 conv factor
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
#ifdef USGS
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            12.00,  12.00,                                                     &
            12.00,   1.00,                                                     &
            12.00,  12.00/
#else
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00/
#endif
#ifdef USGS
   data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)     & ! soil type
         /  14.00,  14.00,                                                     &
            15.00,   0.00,                                                     &
             0.00,   0.00,                                                     &
            15.00,   0.00,                                                     &
             0.00,   0.00/
#else
   data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)     & ! soil type
         /   0.00,   0.00,                                                     &
             8.00,   1.00,                                                     &
             9.00,   9.00,                                                     &
             8.00,   1.00,                                                     &
             9.00,   0.00/
#endif
   data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)     & ! albedo frac
         /   0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00/
   data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)     & ! ustar
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)     & ! drag coeff for m
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)     & ! drag coeff for h
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,joml),isf=1,2),jsf=1,5)     & ! ocena mld
         / 5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0/
   data ((vmaxmin(isf,jsf,jaer),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
   data ((vmaxmin(isf,jsf,jkpr),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jden),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0/
   data ((vmaxmin(isf,jsf,jdxc),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jmix),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
#endif
#ifdef NOALSM1
!
! noalsm1
!
   real                 ::  vmaxmin(2,5,numsfcs)
   integer              ::  isf,jsf
!
   data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)     & ! sst and land temp
         / 313.00, 271.21,                                                     &
           400.00, 100.00,                                                     &
           271.21, 100.00,                                                     &
          273.16, 100.00,                                                      &
           273.16, 100.00/
   data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)     & ! soilm 
         /   0.55,   0.55,                                                     &
             0.55,   0.05,                                                     &
             0.55,   0.55,                                                     &
             0.55,   0.05,                                                     &
             0.55,   0.55/
   data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)     & ! snow depth
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)     & ! soilt 
         / 313.00, 200.00,                                                     &
           353.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)     & ! deep soil temp
         / 310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)     & ! surface roughness
         /   1.00, 1.e-05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05/
   data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)     & ! convective cloud cover
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)     & ! 4-type albedo
         /   0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01/
   data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)     & ! land-sea-seaice mask
         /     0.,     0.,                                                     &
               1.,     1.,                                                     &
               2.,     2.,                                                     &
               1.,     1.,                                                     &
               2.,     2./
   data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)     & ! vegetation cover
         /     0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)     & ! canopy water
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)     & ! 10m and sig1 conv factor
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
#ifdef USGS
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            12.00,  12.00,                                                     &
            12.00,   1.00,                                                     &
            12.00,  12.00/
#else
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00/
#endif
#ifdef USGS
   data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)     & ! soil type
         /  14.00,  14.00,                                                     &
            15.00,   1.00,                                                     &
             0.00,   0.00,                                                     &
            15.00,   1.00,                                                     &
             0.00,   0.00/
#else
   data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)     & ! soil type
         /   0.00,   0.00,                                                     &
             8.00,   1.00,                                                     &
             9.00,   9.00,                                                     &
             8.00,   1.00,                                                     &
             9.00,   9.00/
#endif
   data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)     & ! albedo frac
         /   0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00/
   data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)     & ! ustar
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)     & ! drag coeff for m
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)     & ! drag coeff for h
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jprc),isf=1,2),jsf=1,5)     & ! precip
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jsrf),isf=1,2),jsf=1,5)     & ! sr flag
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jsnd),isf=1,2),jsf=1,5)     & ! snow depth
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jslc),isf=1,2),jsf=1,5)     & ! slc
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jsmn),isf=1,2),jsf=1,5)     & ! min veg cov
         /   1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00/
   data ((vmaxmin(isf,jsf,jsmx),isf=1,2),jsf=1,5)     & ! max veg cov
         /   1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00/
   data ((vmaxmin(isf,jsf,jslo),isf=1,2),jsf=1,5)     & ! slope type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00/
   data ((vmaxmin(isf,jsf,jsna),isf=1,2),jsf=1,5)     & ! snow alb
         /   1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00/
   data ((vmaxmin(isf,jsf,joml),isf=1,2),jsf=1,5)     & ! ocena mld
         / 5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0/
   data ((vmaxmin(isf,jsf,jaer),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
   data ((vmaxmin(isf,jsf,jkpr),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jden),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0/
   data ((vmaxmin(isf,jsf,jdxc),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jmix),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
#endif
#ifdef VICLSM1
!
! viclsm1
!
   real                 ::  vmaxmin(2,5,numsfcs)
   integer              ::  isf,jsf
!
   data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)     & ! sst and land temp
         / 313.00, 271.21,                                                     &
           400.00, 100.00,                                                     &
           271.21, 100.00,                                                     &
           273.16, 100.00,                                                     &
           273.16, 100.00/
   data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)     & ! soilm 
         /   0.55,   0.55,                                                     &
             0.55,   0.05,                                                     &
             0.55,   0.55,                                                     &
             0.55,   0.05,                                                     &
             0.55,   0.55/
   data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)     & ! snow depth
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)     & ! soilt 
         / 313.00, 200.00,                                                     &
           353.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)     & ! deep soil temp
         / 310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)     & ! surface roughness
         /   1.00, 1.e-05,                                                     &
           800.00,   0.05,                                                     &
             1.00,   0.05,                                                     &
           800.00,   0.05,                                                     &
             1.00,   0.05/
   data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)     & ! convective cloud cover
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)     & ! convective cloud base
          /    -1.,    1.,                                                     &  
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)     & ! 4-type albedo
         /   0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01/
   data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)     & ! land-sea-seaice mask
         /     0.,     0.,                                                     &
               1.,     1.,                                                     &
               2.,     2.,                                                     &
               1.,     1.,                                                     &
               2.,     2./
   data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)     & ! vegetation cover
         /     0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)     & ! canopy water
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)     & ! 10m and sig1 conv factor
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
#ifdef USGS
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            12.00,  12.00,                                                     &
            12.00,   1.00,                                                     &
            12.00,  12.00/
#else
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00/
#endif
   data ((vmaxmin(isf,jsf,jrot),isf=1,2),jsf=1,5)     & ! root content
         /   1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00/
   data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)     & ! albedo frac
         /   0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00/
   data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)     & ! ustar
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)     & ! drag coeff for m
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)     & ! drag coeff for h
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jprc),isf=1,2),jsf=1,5)     & ! precip
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jsrf),isf=1,2),jsf=1,5)     & ! sr flag
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jbif),isf=1,2),jsf=1,5)     & ! bif
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jds ),isf=1,2),jsf=1,5)     & ! ds
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jdsm),isf=1,2),jsf=1,5)     & ! dsm
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jws ),isf=1,2),jsf=1,5)     & ! ws
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jcef),isf=1,2),jsf=1,5)     & ! cef
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jexp),isf=1,2),jsf=1,5)     & ! expt
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jkst),isf=1,2),jsf=1,5)     & ! kst
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jdph),isf=1,2),jsf=1,5)     & ! dph
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbub),isf=1,2),jsf=1,5)     & ! bub
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jqrt),isf=1,2),jsf=1,5)     & ! qrt
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbkd),isf=1,2),jsf=1,5)     & ! bkd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsld),isf=1,2),jsf=1,5)     & ! sld
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jwcr),isf=1,2),jsf=1,5)     & ! wcr
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jwpw),isf=1,2),jsf=1,5)     & ! wpw
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsmr),isf=1,2),jsf=1,5)     & ! smr
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsmx),isf=1,2),jsf=1,5)     & ! smx
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jdpn),isf=1,2),jsf=1,5)     & ! dphnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsxn),isf=1,2),jsf=1,5)     & ! smxnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jepn),isf=1,2),jsf=1,5)     & ! exptnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbbn),isf=1,2),jsf=1,5)     & ! bubnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,japn),isf=1,2),jsf=1,5)     & ! alpnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbtn),isf=1,2),jsf=1,5)     & ! betnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jgmn),isf=1,2),jsf=1,5)     & ! gamnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jlai),isf=1,2),jsf=1,5)     & ! flai
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jslz),isf=1,2),jsf=1,5)     & ! silz
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsnz),isf=1,2),jsf=1,5)     & ! snwz
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsic),isf=1,2),jsf=1,5)     & ! sic
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jcsn),isf=1,2),jsf=1,5)     & ! csn
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jrsn),isf=1,2),jsf=1,5)     & ! rsn
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jtsn),isf=1,2),jsf=1,5)     & ! tsf
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jtpk),isf=1,2),jsf=1,5)     & ! tpk
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsfw),isf=1,2),jsf=1,5)     & ! sfw
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jpkw),isf=1,2),jsf=1,5)     & ! pkw
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jlst),isf=1,2),jsf=1,5)     & ! last snowfall
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,joml),isf=1,2),jsf=1,5)     & ! ocena mld
         / 5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0/
   data ((vmaxmin(isf,jsf,jaer),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
   data ((vmaxmin(isf,jsf,jkpr),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jden),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0/
   data ((vmaxmin(isf,jsf,jdxc),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jmix),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
#endif
#ifdef VICLSM2
!
! viclsm2
!
   real                 ::  vmaxmin(2,5,numsfcs)
   integer              ::  isf,jsf
!
   data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)     & ! sst and land temp
         / 313.00, 271.21,                                                     &
           400.00, 100.00,                                                     &
           271.21, 100.00,                                                     &
           273.16, 100.00,                                                     &
           273.16, 100.00/
   data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)     & ! soilm 
         /   0.55,   0.55,                                                     &
             0.55,   0.05,                                                     &
             0.55,   0.55,                                                     &
             0.55,   0.05,                                                     &
             0.55,   0.55/
   data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)     & ! snow depth
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)     & ! soilt 
         / 313.00, 200.00,                                                     &
           353.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)     & ! deep soil temp
         / 310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00,                                                     &
           310.00, 200.00/
   data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)     & ! surface roughness
         /   1.00, 1.e-05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05,                                                     &
           300.00,   0.05,                                                     &
             1.00,   0.05/
   data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)     & ! convective cloud cover
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)     & ! convective cloud base
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)     & ! 4-type albedo
         /   0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01,                                                     &
             0.80,   0.01,                                                     &
             0.01,   0.01/
   data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)     & ! land-sea-seaice mask
         /     0.,     0.,                                                     &
               1.,     1.,                                                     &
               2.,     2.,                                                     &
               1.,     1.,                                                     &
               2.,     2./
   data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)     & ! vegetation cover
         /     0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0.,                                                     &
               1.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)     & ! canopy water
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)     & ! 10m and sig1 conv factor
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
#ifdef USGS
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   0.00,                                                     &
            12.00,  12.00,                                                     &
            12.00,  12.00,                                                     &
            12.00,  12.00/
#else
   data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)     & ! veg type
         /   0.00,   0.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00,                                                     &
            12.00,   1.00,                                                     &
            13.00,  13.00/
#endif
   data ((vmaxmin(isf,jsf,jrot),isf=1,2),jsf=1,5)     & ! root content
         /   1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             1.00,   0.00/
   data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)     & ! albedo frac
         /   0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00,                                                     &
             1.00,   0.00,                                                     &
             0.00,   0.00/
   data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)     & ! ustar
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)     & ! drag coeff for m
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)     & ! drag coeff for h
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jprc),isf=1,2),jsf=1,5)     & ! precip
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jsrf),isf=1,2),jsf=1,5)     & ! sr flag
         /    -1.,     1.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0.,                                                     &
               0.,     0./
   data ((vmaxmin(isf,jsf,jbif),isf=1,2),jsf=1,5)     & ! bif
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jds ),isf=1,2),jsf=1,5)     & ! ds
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jdsm),isf=1,2),jsf=1,5)     & ! dsm
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jws ),isf=1,2),jsf=1,5)     & ! ws
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jcef),isf=1,2),jsf=1,5)     & ! cef
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jexp),isf=1,2),jsf=1,5)     & ! expt
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jkst),isf=1,2),jsf=1,5)     & ! kst
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jdph),isf=1,2),jsf=1,5)     & ! dph
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbub),isf=1,2),jsf=1,5)     & ! bub
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jqrt),isf=1,2),jsf=1,5)     & ! qrt
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbkd),isf=1,2),jsf=1,5)     & ! bkd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsld),isf=1,2),jsf=1,5)     & ! sld
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jwcr),isf=1,2),jsf=1,5)     & ! wcr
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jwpw),isf=1,2),jsf=1,5)     & ! wpw
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsmr),isf=1,2),jsf=1,5)     & ! smr
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsmx),isf=1,2),jsf=1,5)     & ! smx
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jdpn),isf=1,2),jsf=1,5)     & ! dphnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsxn),isf=1,2),jsf=1,5)     & ! smxnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jepn),isf=1,2),jsf=1,5)     & ! exptnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbbn),isf=1,2),jsf=1,5)     & ! bubnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,japn),isf=1,2),jsf=1,5)     & ! alpnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jbtn),isf=1,2),jsf=1,5)     & ! betnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jgmn),isf=1,2),jsf=1,5)     & ! gamnd
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jslz),isf=1,2),jsf=1,5)     & ! silz
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsnz),isf=1,2),jsf=1,5)     & ! snwz
         /  0.0000,  1.00,                                                     &
            0.0005,  1.00,                                                     &
            0.0005,  1.00,                                                     &
            0.0005,  1.00,                                                     &
            0.0005,  1.00/
   data ((vmaxmin(isf,jsf,jnve),isf=1,2),jsf=1,5)     & ! nve
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmfr),isf=1,2),jsf=1,5)     & ! mfr
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmcp),isf=1,2),jsf=1,5)     & ! mcp
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmvt),isf=1,2),jsf=1,5)     & ! mvt
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jlai),isf=1,2),jsf=1,5)     & ! flai
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmsn),isf=1,2),jsf=1,5)     & ! msn
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmsm),isf=1,2),jsf=1,5)     & ! msm
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmsi),isf=1,2),jsf=1,5)     & ! msic
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jmst),isf=1,2),jsf=1,5)     & ! mst
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jcsn),isf=1,2),jsf=1,5)     & ! csn
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jrsn),isf=1,2),jsf=1,5)     & ! rsn
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jtsn),isf=1,2),jsf=1,5)     & ! tsf
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jtpk),isf=1,2),jsf=1,5)     & ! tpk
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jsfw),isf=1,2),jsf=1,5)     & ! sfw
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jpkw),isf=1,2),jsf=1,5)     & ! pkw
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,jlst),isf=1,2),jsf=1,5)     & ! last snowfall
         /  -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00,                                                      &
            -1.00,  1.00/
   data ((vmaxmin(isf,jsf,joml),isf=1,2),jsf=1,5)     & ! ocena mld
         / 5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0,                                                     &
           5000.0,    0.0/
   data ((vmaxmin(isf,jsf,jaer),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
   data ((vmaxmin(isf,jsf,jkpr),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jden),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0,                                                     &
          40000.0,    0.0/
   data ((vmaxmin(isf,jsf,jdxc),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /   10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0,                                                     &
             10.0,    0.0/
   data ((vmaxmin(isf,jsf,jmix),isf=1,2),jsf=1,5)     & ! aerosol distribution
         /    1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0,                                                     &
              1.0,    0.0/
#endif
!
   end module varsfc
!
! end of varsfc.h
!
