#include <define.h>
   subroutine setgrd2(flat,flon,delx,dely,dlamda0)
!-------------------------------------------------------------------------------
!
! subprogram:  setgrd2
!
! abstract:  setup grid latitude, longitude, and map factor etc for
!            regional grid.
!
! program history log:
!   1993-02-06  henry juang            development
!   2000-03-09  songyou hong           cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call setgrd2(flat,flon,delx,dely,dlamda0)
!   output argument list:
!     flat   - latitude
!     flon   - logitude
!     delx   - grid spacing in x direction
!     dely   - grid spacing in y direction
!     dlamda0   - angle between 0 longitude and x positive axies
!
!-------------------------------------------------------------------------------
   use constant, only : pi_,rerth_
#ifdef RMP
   use paramodel, only : lngrd_
#ifdef MP
   use commpi
#endif
   use rscomloc
!-------------------------------------------------------------------------------
   real                 ::  flat(lngrd_),flon(lngrd_)
   real,parameter       ::  pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi
   real,parameter       ::  rad=pi/180.
!-------------------------------------------------------------------------------
#ifdef DBG
#ifdef MP
   if (mype.eq.master) then
#endif
     print 1234
 1234 format(' ==== in routine setgrd2 === ')
#ifdef MP
   endif
#endif
#endif
!
   delx = rdelx
   dely = rdely
!
! --------- setup regional lat/lon and map factor -----
!
! if proj=0  do mercater projection
! if proj=1  do north polar projection
! if proj=-1 do south polar projection
!
   nproj = rproj
!
   if( nproj.eq.1 .or. nproj.eq.-1 ) then
   ! ++++++++++++++++++++++++++++++++++++++
   ! polar projection
   ! ++++++++++++++++++++++++++++++++++++++
     truth  = rtruth * rad
     truth  = nproj * truth
     orient  = rorient * rad
     dlamda0 = orient + hfpi
     a2 =  rerth_ * ( 1.0 + sin(truth) )
     radlat = rcenlat * rad
     radlon = rcenlon * rad - dlamda0
     radlat = nproj * radlat
     radlon = nproj * radlon
     yyy = a2 * cos( radlat )/(1. + sin( radlat ) )
     cenlat = rcenlat
     if( abs(cenlat) .eq. 90. ) then yyy = 0.0
     y00 = yyy * sin( radlon ) - ( rbtmgrd -1.) * dely
     x00 = yyy * cos( radlon ) - ( rlftgrd -1.) * delx
#ifdef DBG
#ifdef MP
     if (mype.eq.master) then
#endif
       print *,' delx x00 y00 ',delx,x00,y00
#ifdef MP
     endif
#endif
#endif
     !
     ! ========= lat loop =======
     !
     do j = 1,jgrd12_
       lats = j
       ijlats = (lats-1)*igrd12_
       ys = y00 + (lats-1)*dely
       latn = jgrd1_+1-j
       ijlatn = ijlats+igrd1_
       yn = y00 + (latn-1)*dely
       !
       do i = 1,igrd1_
         x = x00 + (i-1)*delx
         if( x .gt. 0. e 0 ) then
           flons = atan(ys/x)
           flonn = atan(yn/x)
         else if ( x .lt. 0. e 0 ) then
           flons = pi + atan(ys/x)
           flonn = pi + atan(yn/x)
         else
           flons = hfpi
           if( ys .lt. 0. e 0 ) flons = flons * 3.0
           flonn = hfpi
           if( yn .lt. 0. e 0 ) flonn = flonn * 3.0
         endif
         flons = nproj * flons + dlamda0
         flons = mod(flons,twopi)
         if(flons.lt.0. e 0) flons = twopi + flons
         flonn = nproj * flonn + dlamda0
         flonn = mod(flonn,twopi)
         if(flonn.lt.0. e 0) flonn = twopi + flonn
         !
         rsoa2 = sqrt( x*x + ys*ys )/a2
         flats = hfpi - 2.0 e 0 * atan(rsoa2)
         flat(ijlats+i) = nproj * flats
         flon(ijlats+i) = flons
         !
         rsoa2 = sqrt( x*x + yn*yn )/a2
         flatn = hfpi - 2.0 e 0 * atan(rsoa2)
         flat(ijlatn+i) = nproj * flatn
         flon(ijlatn+i) = flonn
         !
       enddo
     enddo
     !
   else if ( nproj.eq.0 ) then
   !
   ! ++++++++++++++++++++++++++++
   ! do mercater
   ! ++++++++++++++++++++++++++++
     truth  = rtruth * rad
     cenlat = rcenlat * rad
     cenlon = rcenlon * rad
     cenlon = mod(cenlon,twopi)
     if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
     a2 =  rerth_ * cos( truth )
     x0 = 0.0
     y0 = a2 * log( abs( tan( qtpi + 0.5 * cenlat ) ) )
     x00 = -( rlftgrd - 1. ) * delx - 0.5 * delx
     y00 = -( rbtmgrd - 1. ) * dely - 0.5 * dely
     dlamda0 = 0.0
     !
     do j = 1,jgrd12_
       lats = j
       ijlats = (lats-1)*igrd12_
       ys = y00 + (lats-1)*dely + y0
       latn = jgrd1_+1-j
       ijlatn = ijlats+igrd1_
       yn = y00 + (latn-1)*dely + y0
       !
       do i = 1,igrd1_
         x = x00 + (i-1)*delx + x0
         flons = x / a2 + cenlon
         flons = mod(flons,twopi)
         if(flons.lt.0. e 0) flons = twopi + flons
         flonn = flons
         !
         flats = 2.0 *( atan( exp( ys/a2 ) ) - qtpi )
         flat(ijlats+i) = flats
         flon(ijlats+i) = flons
         !
         flatn = 2.0 *( atan( exp( yn/a2 ) ) - qtpi )
         flat(ijlatn+i) = flatn
         flon(ijlatn+i) = flonn
         !
       enddo
     enddo
     !
   endif
!
#endif
!
   return
   end
