#include <define.h>
   subroutine rad_diurnal_cycle(ims2,imx2,kmx,                                 &
                    solhr,slag,sinlab,coslab,sdec,cdec,                        &
                    xlon,czmn,                                                 &
                    sfcdlw,sfcnsw,                                             &
#ifdef VIC
                    sfcdsw,dswsfc,                                             &
#endif
                    tf,tsea,tsflw,swh,hlw,                                     &
#ifdef NIM_DIAG
                    xmu,                                                       &
#endif
                    dlwsfc,ulwsfc,slrad,tau)
!-------------------------------------------------------------------------------
!
! subroutine: rad_diurnal_cycle     
!
! abstract: a diurnal cycle approximation is applied to previously
!   computed radiative fluxes and heating rates. first,the current
!   local-time value (for this particular model time step) of the
!   cosine solar zenith angle (cosz) is computed for all gaussian grid
!   points. shortwave (sw) heating rates which were computed with
!   latitudinal mean cosz in the separate radiation calculation
!   are weighted by the ratio of actual to mean cosz (see mrf model
!   documentation,1988,chapter 3,'radiative processes',authored by
!   k. campana,..). surface sw fluxes are also cosz weighted. surface
!   longwave (lw) flux from the atmosphere is altered each timestep
!   to account for diurnal changes of model temperature in the lower
!   atmosphere. lw heating rates from the separate radiation
!   computation are untouched.
!
! program history log:
!   1988-05-06  kenneth campana        development
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call rad_diurnal_cycle(lat,ssdec,solhr,colrad,czmn,sfcdlw,sfcnsw,
!                      tf,tov,slrad,swh,hlw,tau)
!   input argument list:
!     lat      - row number of gaussian latitude(n.h.).
!     ssdec    - sine of the solar declination for todays date-
!                         part of the output from radiation codes.
!     solhr    - time in hours after 00 hr greenwich.
!     colrad   - co-latitudes of gaussian grid in radians(n.h.).
!     czmn     - mean cosine solar zenith angle for all gaussian lats-
!                         part of the output from radiation codes.
!     sfcdlw   - downward lw flux at earth sfc(from radiation code)
!                in cal cm-2 min-1.
!     sfcnsw   - net sw flux at earth sfc(from radiation code using
!                czmn) in cal cm-2 min-1.
!     tf       - current value of model temperatures(w/o basic state)
!                in deg k.
!     tov      - basic state temperature for all model layers(deg k).
!     swh      - model lyr sw heating rates(from radiation code,using
!                czmn) - deg/sec.
!     hlw      - model lyr lw heating rates(from radiation code),
!                in deg/sec.
!
!   output argument list:
!     slrad    - surface net radiative flux (except lw upward flux
!                from sfc,which is added in land sfc model) -
!                units are cal cm-2 min-1 .
!     tau      - layer values of temperature tendency after adding
!                altered radiative heating rates-units are deg/sec.
!
!-------------------------------------------------------------------------------
   use paramodel, only  :  ILOTS,levs_
   use constant, only   :  cal_,pi_,hsigma=>sbc_
!-------------------------------------------------------------------------------
   real,parameter       ::  cnwatt=-cal_*1.e4/60.
   integer              ::  ims2,imx2,kmx
   real                 ::  xlon(imx2),czmn(imx2),sfcdlw(imx2),sfcnsw(imx2)
   real                 ::  tf(imx2),tsea(imx2),tsflw(imx2),sinlab(imx2)
   real                 ::  swh(imx2,kmx),hlw(imx2,kmx),coslab(imx2)
   real                 ::  dlwsfc(imx2),ulwsfc(imx2),slrad(imx2),tau(imx2,kmx)
#ifdef VIC
   real                 ::  dswsfc(imx2),sfcdsw(imx2)
#endif
!
! local dimension
!
   real                 ::  xmu(ILOTS)
#ifdef NO_DIURNAL
   real                 ::  xmum(ILOTS)
#endif
!-------------------------------------------------------------------------------
   lon2=ims2
   levs=kmx
!
!  compute cosine of solar zenith angle for both hemispheres.
!
#ifdef NO_DIURNAL
#ifdef DBG
   print *, 'sdec cdec in rad_diurnal_cycle', sdec, cdec
!
#endif
#endif
   cns=pi_*(solhr-12.)/12.+slag
!
   do i = 1,lon2
     ss=sinlab(i)*sdec
     cc=coslab(i)*cdec
     ch=cc*cos(xlon(i)+cns)
     xmu(i)=ch+ss
!
!    xmu(i)=(sinlab(i)*sdec)+(coslab(i)*cdec)*cos(xlon(i)+cns)
!
#ifdef NO_DIURNAL
!
! daily mean cos zenith angle : performing the simple integration
!                               from sunrise to sunset
!
! -> xmu(i) = (sinlab(i)*sdec*H + coslab(i)*cdec*sin(H))/pi_
!   H [radians] ; a half-day (from sunrise or sunset to solar noon)
!                 H = acos(-1*tanlab(i)*tandec) = acos(-1.*ss/cc)
!
     hday=acos(-1.*ss/cc)
     xmum(i)=(ss*hday+cc*sin(hday))/pi_
#endif
   enddo
!
   do i = 1,lon2
!
! normalize by average value over radiation period for daytime.
!
#ifdef NO_DIURNAL
!
! xmu is the same during the day and night
!
     xmu(i)=xmum(i)/czmn(i)
#else
     if(xmu(i).gt.0.01.and.czmn(i).gt.0.01) then
       xmu(i)=xmu(i)/czmn(i)
     else
       xmu(i)=0.
     endif
#endif
!
! adjust longwave flux at surface to account for t changes in layer 1.
!
     sdlw=sfcdlw(i)*(tf(i)/tsflw(i))**4
#ifdef SMP_NUDGING_RAD
     sdlw=sfcdlw(i)
     print*,'rad_diurnal_cycle obs sdlw',sdlw
#endif
#ifdef VIC
#ifndef SMP_NUDGING_RAD
     dswsfc(i)=sfcdsw(i)*xmu(i)
#else
     dswsfc(i)=sfcdsw(i)
#endif
#endif
!
! return net surface radiative flux.
!
#ifndef SMP_NUDGING_RAD
     slrad(i)=sfcnsw(i)*xmu(i)+sdlw
#else
     slrad(i)=sfcnsw(i)+sdlw
#endif
!
! return downward and upward longwave flux at ground, respectively.
!
     dlwsfc(i)=sdlw*cnwatt
     ulwsfc(i)=hsigma*tsea(i)**4
   enddo
#ifdef NO_DIURNAL
#ifdef DBG
   print*, 'xmu,sfcnsw(75) in rad_diurnal_cycle=',xmu(75), sfcnsw(75)*xmu(75)
#endif
#endif
!
! add radiative heating to temperature tendency
!
   do k = 1,levs
     do i = 1,lon2
       tau(i,k)=tau(i,k)+swh(i,k)*xmu(i)+hlw(i,k)
     enddo
   enddo
!
   return
   end subroutine rad_diurnal_cycle
