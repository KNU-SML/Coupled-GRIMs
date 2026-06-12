#include <define.h>
   subroutine mpdimset(jcap,levs,lonf,latg)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpdimset
!            
! abstract: preset all starting point and length for 
!           all pe for global spectral model.
!
! program history log:
!    99-06-27  henry juang 	finish entire test for gsm
!
! usage:	call mpdimset(jcap,levs,lonf,latg)
!
!    input argument lists:
!	jcap	- integer spectral wavenumber
!	levs	- integer vertical layer number
!	lonf	- integer gaussian grid for longitude
!	latg	- integer gaussian grid for latitude
!
!    output argument list:
! 
! subprograms called:
!   dimset	- to compute all dimension, starting point for all pe
!
!-------------------------------------------------------------------------------
   use paramodel, only : levsp_, lonfp_, latgp_, jcapp_, lntp_, llnp_
   use commpi       ! npes, ncol, nrow, levstr,levlen, lerstr,lerlen,          &
                    ! lonstr,lonlen, latstr,latlen, lwvstr,lwvlen,             &
                    ! lnpstr,lnplen, lntstr,lntlen, lwvdef,latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer jcap,levs,lonf,latg                                                 &
          ,levmax,lonmax,latmax,lntmax,lnpmax,lwvmax                           &
          ,lermax
!
   call dimset(jcap,levs,lonf,latg,                                            &
               npes,ncol,nrow,                                                 &
               levstr,levlen,levmax,                                           &
               lerstr,lerlen,lermax,                                           &
               lonstr,lonlen,lonmax,                                           &
               latstr,latlen,latmax,                                           &
               lwvstr,lwvlen,lwvmax,                                           &
               lntstr,lntlen,lntmax,                                           &
               lnpstr,lnplen,lnpmax,                                           &
               lwvdef,latdef)
!
   if( levmax .ne. levsp_ ) then
     print *,' levmax levsp_ ',levmax,levsp_
     print *,' Error in mpdimset for levmax '
     call mpabort
   endif
!
   if( lonmax .ne. lonfp_ ) then
     print *,' lonmax lonfp ',lonmax,lonfp_
     print *,' Error in mpdimset for lonmax '
     call mpabort
   endif
!
   if( latmax .ne. latgp_ ) then
     print *,' latmax latgp_ ',latmax,latgp_
     print *,' Error in mpdimset for latmax '
     call mpabort
   endif
!
   if( lwvmax .ne. jcapp_ ) then
     print *,' lwvmax jcapp_ ',lwvmax,jcapp_
     print *,' Error in mpdimset for lwvmax '
     call mpabort
   endif
!
   if( lntmax .ne. lntp_ ) then
     print *,' lntmax lntp_ ',lntmax,lntp_
     print *,' Error in mpdimset for lntmax '
     call mpabort
   endif
!
   if( lnpmax .ne. llnp_ ) then
     print *,' lnpmax llnp_ ',lnpmax,llnp_
     print *,' Error in mpdimset for lnpmax '
     call mpabort
   endif
!
   return
   end subroutine mpdimset
!-------------------------------------------------------------------------------
