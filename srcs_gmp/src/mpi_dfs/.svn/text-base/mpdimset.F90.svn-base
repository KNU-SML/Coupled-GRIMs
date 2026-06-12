#include "define.h"
   subroutine mpdimset(jcap,lcap,levs,lonf,latg,igs,jgs,mgs,ngs,latdef)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpdimset
!            
! abstract: preset all starting point and length for 
!           all pe for global spectral model.
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
   use commpi, only : npes,ncol,nrow,mype,                                     &
               levstr,levlen,lonstr,lonlen,                                    &
               lerstr,lerlen,                                                  &
               latstr,latlen,lmstr,lmlen,                                      &
               lnstr,lnlen
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer jcap,lcap,levs,lonf,latg
   integer igs,jgs,mgs,ngs
!
   integer levmax,lonmax,latmax,lntmax,lnpmax,lwvmax
   integer lmmax,lnmax,lnmax1,latmax1,lermax
   integer latdef(latg/2)
#ifdef DBG
   integer i
   character*100 fmt
#endif
!
   if (mod(latg,ncol)/=0) then
     print*,'Invalid ncol numver in mpdimset',ncol
     call mpabort
   endif
!
   call dimset(jcap,lcap,levs,lonf,latg,                                       &
               npes,ncol,nrow,                                                 &
               levstr,levlen,levmax,                                           &
               lerstr,lerlen,lermax,                                           &
               lonstr,lonlen,lonmax,                                           &
               latstr,latlen,latmax,                                           &
#ifdef DFS
               lmstr,lmlen,lmmax,                                              &
               lnstr,lnlen,lnmax,latdef,                                       &
#else
               lwvstr,lwvlen,lwvmax,                                           &
               lntstr,lntlen,lntmax,                                           &
               lnpstr,lnplen,lnpmax,                                           &
               lwvdef,latdef,                                                  &
#endif
               igs,jgs,mgs,ngs)
#ifdef DBG
   if (mype.eq.0) then
      write(fmt,'(A3,I5,A3)')'(A,',npes,'I5)'
      write(6,trim(fmt))'latstr=',latstr(:)
      write(6,trim(fmt))'latlen=',latlen(:)
      write(6,trim(fmt))'lonstr=',lonstr(:)
      write(6,trim(fmt))'lonlen=',lonlen(:)
      write(6,trim(fmt))'lmstr=',lmstr(:)
      write(6,trim(fmt))'lmlen=',lmlen(:)
      write(6,trim(fmt))'lnstr=',lnstr(:)
      write(6,trim(fmt))'lnlen=',lnlen(:)
   endif
#endif
#ifdef OLD_1
   if( levmax .ne. levsp_ ) then
     print *,' levmax levsp_ ',levmax,levsp_
     stop ' Error in mpdimset for levmax '
   endif
!
   if( lonmax .ne. lonfp_ ) then
     print *,' lonmax lonfp ',lonmax,lonfp_
     stop ' Error in mpdimset for lonmax '
   endif
!
   if( latmax .ne. latgp_ ) then
     print *,' latmax latgp_ ',latmax,latgp_
     stop ' Error in mpdimset for latmax '
   endif
   if( lwvmax .ne. jcapp_ ) then
     print *,' lwvmax jcapp_ ',lwvmax,jcapp_
     stop ' Error in mpdimset for lwvmax '
   endif
!
   if( lntmax .ne. lntp_ ) then
     print *,' lntmax lntp_ ',lntmax,lntp_
     stop ' Error in mpdimset for lntmax '
   endif
!
   if( lnpmax .ne. llnp_ ) then
     print *,' lnpmax llnp_ ',lnpmax,llnp_
     stop ' Error in mpdimset for lnpmax '
   endif
#endif
!
   return
   end subroutine mpdimset
!-------------------------------------------------------------------------------
