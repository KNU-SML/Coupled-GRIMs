#include "define.h"
   subroutine rad_ozone_time_interp(fhour,idate,ko3,o3out,pstr,jerr)
!------------------------------------------------------------------------------- 
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!     *******************************************************!
!     *  computes o3 climo from 12 month dataset, linearly   !
!     *   interpolated to day,mon of the fcst.  then create  !
!     *   a 5 deg array from the 10 deg climatology...for    !
!     *   ease when doing a latitudinal interpolation        !
!     *  thanks to s moorthi for new o3 climo...kac  dec 1996!
!     * input:                                               !
!     *   idate=nmc date-time                                !
!     *   fhour=forecast hour                                !
!     *   ko3=unit number of o3 climatology                  !
!     * output :                                             !
!     *   o3out=5-deg o3 climo for forecast date(np->spole)  !
!     *   pstr=pressure (mb) for the climo lyrs (k=1 is top) !
!     *   jerr=0 if o3 file exists, =1 if not (gfdl=default) !
!     *******************************************************!
!
!------------------------------------------------------------------------------- 
#if defined(DFS)
   use dfsvar, only : iope
#endif
   use comio
#ifdef MP
   use commpi
#endif
!------------------------------------------------------------------------------- 
   implicit none
!------------------------------------------------------------------------------- 
!
! geos ozone data
!
   real,parameter     ::  blte=-85.0, dlte=10.0
   integer,parameter  ::  jmr=18
   integer,parameter  ::  loz=17,jmout=37
!
   integer            ::  imon(12),days(12),idate(4)
   integer            ::  ilat(jmr,12)
   integer            ::  ida, imo, imo1, jday, ken
   integer            ::  numdyz, nmdtot, ndayr, nm, mday, ko3
   integer            ::  monl, monc, monr, midl, midc, midr
   integer            ::  j, j1, j2, jmr1, jerr, l
   real               ::  fhour
   real               ::  o3r(jmr,loz,12),o3tmp(jmr,loz),o3tmp2(jmr,loz)
   real               ::  o3out(jmout,loz),pstr(loz)
   real               ::  coef, difl, difr, delday 
   data  days/31,28,31,30,31,30,31,31,30,31,30,31/
!
   ida=idate(3)
   imo=idate(2)
!
!   find current day and month, initial values in ida,imo!
!       will not worry about leap year, since it will take a
!       120-year (what?) forecast to be off by 1 month.  if this
!       is deemed a problem, need to redo this calculation.
!
   if (fhour.ge.24.) then
!
!  number of days into the forecast
!
     numdyz=int(fhour/24.0 + 0.01)
!
!  get day-of-year, remember climate runs are for years
!
     imo1=imo-1
     jday = ida
!
     if (imo1.gt.0) then
       jday=0
       do ken = 1,imo1
         jday=jday+days(ken)
       enddo
       jday=jday+ida
     endif
!
     nmdtot = jday+numdyz
     ndayr = mod(nmdtot,365)
     if (ndayr.eq.0) ndayr=365
!
!   now get month from day-of-year
!
     mday=0
     do ken = 1,11
       mday=mday+days(ken)
       imo=ken
       if (ndayr.le.mday) then
         ida=ndayr-(mday-days(imo))
         !
         go to 9
         !
       endif
     enddo
!
     imo=12
9    continue
!
#ifndef NOPRINT
     if(iope) then
!      write(6,66)fhour,numdyz,jday,nmdtot,ndayr
       print*, 'sbuvo3 climo hr,numdyz,jday,nmdtot,ndayr=',                    &
                fhour,numdyz,jday,nmdtot,ndayr
     endif
!
66   format(' sbuvo3 climo hr=',f10.1,                                         &
            ' numdyz,jday,nmdtot,ndayr=',4i8)
#endif
   endif
!
   coef=1.655e-6
   jerr=0
!
   if( iope ) then
     rewind ko3
!
     do l = 1,loz
       read (ko3,15,err=998,end=999) pstr(l)
15     format(f10.3)
     enddo
!
16   format(1h ,' o3 pressures=',8f10.3)
!
     do nm = 1,12
       do j = 1,jmr
         read (ko3,19,err=998,end=999) imon(nm),ilat(j,nm),                    &
                                       (o3r(j,l,nm), l=1,10)
         read (ko3,20,err=998,end=999) (o3r(j,l,nm), l=11,loz)
       enddo
       do j = 1,jmr
         do l = 1,loz
           o3r(j,l,nm) = o3r(j,l,nm) * coef
         enddo
       enddo
     enddo
19   format(i2,i4,10f6.2)
20   format(6x,7f6.2)
!
   endif
#ifdef MP
#ifdef RMP
   call rmpbcastr(pstr,loz)
   call rmpbcastr(o3r,loz*jmr*12)
   call rmpbcasti(imon,12)
   call rmpbcasti(ilat,jmr*12)
#else
   call mpbcastr(pstr,loz)
   call mpbcastr(o3r,loz*jmr*12)
   call mpbcasti(imon,12)
   call mpbcasti(ilat,jmr*12)
#endif
#endif
!
!      do a linear interpolation in time, where we assume that
!       the ozone data is valid for mid-month
!      monl is the preceeding month, monc for current mo, and
!      monr is the future month..
!
   monl=imo-1
   monc=imo
   monr=imo+1
   if (monl.lt.1) monl=12
   if (monr.gt.12) monr=1
!
!  difl=number of days beteen mid-months of the current and
!            preceeding mo, difr=same for current and future mo..
!  dell=number of days between current day and mon,
!       delr=same for current day and next month.
!       sign convention as if we were using day of year calculations.
!
   midl=days(monl)/2
   midc=days(monc)/2
   midr=days(monr)/2
   difl=-(days(monl)-midl+midc)
   difr= (days(monc)-midc+midr)
   delday=ida-midc
!
   if (ida.gt.midc) then
     do j = 1,jmr
       do l = 1,loz
         o3tmp(j,l)=o3r(j,l,monc)+(o3r(j,l,monr)-o3r(j,l,monc))                &
                     * delday/difr
       enddo
     enddo
   else if (ida.lt.midc) then
     do j = 1,jmr
       do l = 1,loz
         o3tmp(j,l)=o3r(j,l,monc)+(o3r(j,l,monl)-o3r(j,l,monc))                &
                     * delday/difl
       enddo
     enddo
   else if (ida.eq.midc) then
     do j = 1,jmr
       do l = 1,loz
         o3tmp(j,l)=o3r(j,l,monc)
       enddo  
     enddo
   end if
!
#ifndef NOPRINT
   if(iope) then
     write(6,200)imo,ida
   endif
200 format(1x,'ozone climatology for month,day=',2i4)
#endif
!
!  flip y-direction for ngases = 1 for time being
!
   do j = 1,jmr
     do l = 1,loz
       o3tmp2(j,l)=o3tmp(jmr-j+1,l)
     enddo
   enddo
!
!  linearly interpolate to 5 deg zonal means
!
   jmr1=jmr-1
   do j = 1,jmr1
     j1=j*2
     j2=j1+1
     do l = 1,loz
       o3out(j1,l)=O3TMP2(j,l)
       o3out(j2,l)=0.5*(O3TMP2(j,l)+O3TMP2(j+1,l))
     enddo
   enddo
!
   do l = 1,loz
     o3out(1,l)=O3TMP2(1,l)
     o3out(jmout-1,l)=O3TMP2(jmr,l)
     o3out(jmout,l)=O3TMP2(jmr,l)
   enddo
!
   return
#ifdef MP
#ifdef RMP
#define FAIL rmpabort
#else
#define FAIL mpabort
#endif
#else
#define FAIL abort
#endif
998 write(6,988)ko3
   call FAIL
999 write(6,989)ko3
   call FAIL
  988 format(1h ,'error reading nasa ozone, unit=',i4)
  989 format(1h ,'e.o.f reading nasa ozone, unit=',i4)
#undef FAIL
!
   end subroutine rad_ozone_time_interp
