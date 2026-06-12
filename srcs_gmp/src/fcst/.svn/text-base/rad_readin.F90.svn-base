#include <define.h>
   subroutine rad_readin
!-------------------------------------------------------------------------------
!
!  ::: structure ::: 
!
!    [rad_readin]
!        |-------  For NEW_READ_IN in GRIMs
!            |---  [rad_read_init]
!            |---  [rad_read_prepare]
!
!
!-------------------------------------------------------------------------------
  end subroutine rad_readin
!-------------------------------------------------------------------------------
   subroutine rad_read_init
!-------------------------------------------------------------------------------
#ifndef RMP
   use comgrad 
#else
   use rscomgrad
#endif
   use comfcst, only            : rco2
#ifdef DFS
   use dfsvar,  only            : iope
#else
   use comio,   only            : iope
#endif
   use comfver, only            : thour,kdt
   use comreado3
#include "abort.h"
#ifdef MP
   use commpi
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer, parameter           ::  nfile=15, nozon=48
#ifdef MP
   integer, parameter           ::  cld=3, seal=2
   integer, parameter           ::  bin=100, lon=2, lat=4
#endif
   integer                      ::  nall
   real                         ::  fhour
!-------------------------------------------------------------------------------
!
!  from rad_initialize
!
   if (iope) then
     call rad_cloud_read(rhcl,ier)
   endif
#ifdef MP
#ifndef RMP
   call mpbcasti(ier,1)
#else
   call rmpbcasti(ier,1)
#endif
#endif
   if (ier.lt.0) then
     istrat = 0 
     if( iope ) then
       write(6,*)'===>tuning tables not available..abort'
     endif
#ifdef MP
     call MPABORT
#endif
   end if
#ifdef MP
   nall=bin*lon*lat*cld*seal
#ifndef RMP
   call mpbcastr(rhcl,nall)
#else
   call rmpbcastr(rhcl,nall)
#endif
#endif
!
#ifdef DBG
   if(iope) write(6,*)'..for diagnosed clds....istrat = ',istrat
#endif
!
   call rad_constant_init
   call rad_co2_read(nfile,rco2)
!
   return
   end subroutine rad_read_init
!-------------------------------------------------------------------------------
   subroutine rad_read_prepare(fhour,idate)
!-------------------------------------------------------------------------------
#ifndef RMP
   use comgrad
#else
   use rscomgrad
#endif
   use paramodel, only          : ngases_
   use comfver,   only          : thour,kdt
   use comreado3
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer, parameter           ::  nozon=48
   integer,parameter            ::  no3p=28, no3l=29
   integer                      ::  idate(4)
   real                         ::  fhour
!-------------------------------------------------------------------------------
!
!   from rad_prepare
!
   if (thour.eq.0 .or. ngases_.lt.1) then
     call rad_ozone_time_interp(fhour,idate,nozon,o3out,pstr,jerr)
   endif
!
!   from phys_main_driver
!
!   if (ngases_ .ge. 1) then
!      call rad_ozone_setup(idate,fhour,no3p,no3l)
!   endif
!
   return
   end subroutine rad_read_prepare
!-------------------------------------------------------------------------------
