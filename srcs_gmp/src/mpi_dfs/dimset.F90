#include "define.h"
   subroutine dimset(jcap,lcap,levs,lonf,latg,                                 &
                     npes,ncol,nrow,                                           &
                     levstr,levlen,levmax,                                     &
                     lerstr,lerlen,lermax,                                     &
                     lonstr,lonlen,lonmax,                                     &
                     latstr,latlen,latmax,                                     &
#ifdef DFS
                     lmstr,lmlen,lmmax,                                        &
                     lnstr,lnlen,lnmax,                                        &
                     latdef,                                                   &
#else
                     lwvstr,lwvlen,lwvmax,                                     &
                     lntstr,lntlen,lntmax,                                     &
                     lnpstr,lnplen,lnpmax,                                     &
                     lwvdef,latdef,                                            &
#endif
                     igs,jgs,mgs,ngs)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	dimset
!            
! abstract: preset all starting point and length for 
!           all pe for global spectral model.
!
! usage:	call dimset(jcap,levs,lonf,latg,
!    *                  npes,ncol,nrow,
!    *                  levstr,levlen,levmax,
!    *                  lerstr,lerlen,lermax,
!    *                  lonstr,lonlen,lonmax,
!    *                  latstr,latlen,latmax,
!    *                  lwvstr,lwvlen,lwvmax,
!    *                  lntstr,lntlen,lntmax,
!    *                  lnpstr,lnplen,lnpmax,
!    *                  lwvdef,latdef,
!    *                  igs,jgs,mgs,ngs)
!
!    input argument lists:
!	jcap	- integer spectral wavenumber
!	levs	- integer vertical layer number
!	lonf	- integer gaussian grid for longitude
!	latg	- integer gaussian grid for latitude
!	npes	- integer number of pe used: npes=ncol*nrow
!	ncol	- integer number of column
!	nrow	- integer number of nrow
!
!    output argument list:
!	lev*	- integer (npes) related to layers for each pe
!	lon*	- integer (npes) related to longitude for each pe
!	lat*	- integer (npes) related to latitude for each pe
!	lnt*	- integer (npes) related to npes cut of spectral
!	lnp*	- integer (npes) related to nrow cut of spectral
!	lwv*	- integer (npes) related to group of spectral in l
!	*str	- integer (npes) related to each starting point
!	*len	- integer (npes) related to each length
!	*max	- integer related to maximal length of the kind
!	lwvdef	- integer (jcap+1) index of l wave after distribution
!	latdef	- integer (latg/2) index of latitude after distribution
! 
! subprograms called:
!   equdiv	- to compute about equal number of subgroup by division
!   equdis	- to compute about equal number of subgroup by distribution
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer jcap,lcap,levs,lonf,latg,npes,ncol,nrow                             &
          ,latg2,lonf2,nr,nc,nn,jcaprm,nremain                                 &
          ,lnp,lx,lh,n,lerpnt                                                  &
#ifdef DFS
          ,levmax,lonmax,latmax,lmmax,lnmax                                    &
          ,levpnt,lonpnt,latpnt,lmpnt,lnpnt,lermax                             &
#else
          ,levmax,lonmax,latmax,lwvmax,lntmax,lnpmax                           &
          ,levpnt,lonpnt,latpnt,lwvpnt,lntpnt,lnppnt                           &
#endif
          ,igs,jgs,mgs,ngs
   integer                                                                     &
           levstr(0:npes-1),levlen(0:npes-1)                                   &
          ,lerstr(0:npes-1),lerlen(0:npes-1)                                   &
          ,lonstr(0:npes-1),lonlen(0:npes-1)                                   &
          ,latstr(0:npes-1),latlen(0:npes-1)                                   &
#ifdef DFS
          ,lmstr(0:npes-1),lmlen(0:npes-1)                                     &
          ,lnstr(0:npes-1),lnlen(0:npes-1)                                     &
          ,latdef(latg/2)
#else
          ,lwvstr(0:npes-1),lwvlen(0:npes-1)                                   &
          ,lntstr(0:npes-1),lntlen(0:npes-1)                                   &
          ,lnpstr(0:npes-1),lnplen(0:npes-1)                                   &
          ,lwvdef(jcap+1),latdef(latg/2)
#endif
   integer, allocatable ::                                                     &
             levdis(:),londis(:),latdis(:),lerdis(:)
#ifdef DFS
   integer, allocatable ::                                                     &
             lmdis(:),lndis(:),lntdis(:)
#else
   integer, allocatable ::                                                     &
             lwvdis(:),lnpdis(:),lntdis(:)
#endif
!
   allocate (londis(ncol))
   allocate (latdis(nrow))
   allocate (levdis(ncol))
   allocate (lerdis(nrow))
   allocate (lntdis(ncol))
#ifdef DFS
   allocate (lndis(ncol))
   allocate (lmdis(nrow))
#else
   allocate (lwvdis(nrow))
   allocate (lnpdis(nrow))
#endif
!
   call equdiv(levs ,ncol,levdis)
   call equdiv(levs ,nrow,lerdis)
#ifdef DFS
   call equdiv(lcap+1,ncol,lndis)    ! 0:JLA
   call equdiv(jcap+1,nrow,lmdis)    ! 0:MTA
   if (mod(lonf,ncol).ne.0) then
     write(6,*)'in dimset:lon & ncol error',lonf ,ncol
     call mpabort
   endif
   call equdiv(lonf ,ncol,londis)
#else
   call equdis( 1,1,jcap+1,nrow,lwvdis,lwvdef)
!
   lh=0
   do nr = 1,nrow
     lnpdis(nr)=0
     lx=lh+1
     lh=lh+lwvdis(nr)
     do n = lx,lh
       lnpdis(nr)=lnpdis(nr)+jcap+2-lwvdef(n)
     enddo
   enddo
   !
   call equdiv(lonf ,ncol,londis)
#endif
   latg2=latg/2
   if (nrow.gt.latg2) then
     write(6,*)'ERRor in latdiv:latg2,nrow=',latg2,nrow
     call mpabort
   endif
   call equdis(1,jgs,latg2+jgs-1,nrow,latdis,latdef)
!
   levmax=0
   lonmax=0
   latmax=0
   lermax=0
   lerpnt=1
#ifdef DFS
   latpnt=jgs
   lmpnt=mgs
   lnmax=0
   lmmax=0
!
#else
   lntpnt=1
   latpnt=1
   lwvpnt=1
   lnppnt=1
   lntmax=0
   latmax=0
   lwvmax=0
   lnpmax=0
#endif
   n=0
!
   do nr = 1,nrow
!
     levpnt=1
#ifdef DFS
     lonpnt=igs
     lnpnt=ngs
#else
     lonpnt=1
     call equdiv(lnpdis(nr),ncol,lntdis)
#endif
!
     do nc = 1,ncol
!
       levstr(n)=levpnt
       levlen(n)=levdis(nc)
       levpnt=levpnt+levdis(nc)
       levmax=max(levmax,levlen(n))
!
       lerstr(n)=lerpnt
       lerlen(n)=lerdis(nr)
       lermax=max(lermax,lerlen(n))
#ifdef DFS
       lonstr(n)=lonpnt
       lonlen(n)=londis(nc)
       lonpnt=lonpnt+londis(nc)
       lonmax=max(lonmax,lonlen(n))
!
       lnstr(n)=lnpnt
       lnlen(n)=lndis(nc)
       lnpnt=lnpnt+lndis(nc)
       lnmax=max(lnmax,lnlen(n))
!
       latstr(n)=latpnt
       latlen(n)=latdis(nr)
       latmax=max(latmax,latlen(n))
!
       lmstr(n)=lmpnt
       lmlen(n)=lmdis(nr)
       lmmax=max(lmmax,lmlen(n))
#else
       lonstr(n)=lonpnt
       lonlen(n)=londis(nc)
       lonpnt=lonpnt+londis(nc)
       lonmax=max(lonmax,lonlen(n))
!
       lntstr(n)=lntpnt
       lntlen(n)=lntdis(nc)
       lntpnt=lntpnt+lntdis(nc)
       lntmax=max(lntmax,lntlen(n))
!
       latstr(n)=latpnt
       latlen(n)=latdis(nr)
       latmax=max(latmax,latlen(n))
!
       lwvstr(n)=lwvpnt
       lwvlen(n)=lwvdis(nr)
       lwvmax=max(lwvmax,lwvlen(n))
!
       lnpstr(n)=lnppnt
       lnplen(n)=lnpdis(nr)
       lnpmax=max(lnpmax,lnplen(n))
#endif
!
       n=n+1
!
     enddo
!
     if( nr.lt.nrow ) then
       lerpnt=lerpnt+lerdis(nr)
#ifdef DFS
       latpnt=latpnt+latdis(nr)
       lmpnt=lmpnt+lmdis(nr)
#else
       latpnt=latpnt+latdis(nr)
       lwvpnt=lwvpnt+lwvdis(nr)
       lnppnt=lnppnt+lnpdis(nr)
#endif
     endif
!
   enddo
!
#ifdef DFS
   do n = 0,npes-1
     lmlen(n)=lmlen(n)*2                       ! -mle:-mls,mls:mle
     if (lmstr(n)==0) lmlen(n)=lmlen(n)-1      ! -mle:mle, with 0
   enddo
   lmmax=maxval(lmlen)
#else
   do n = 0,npes-1
     lwvstr(n)=lwvstr(n)-1
     lnpstr(n)=lnpstr(n)-1
     lntstr(n)=lntstr(n)-1
   enddo
   do n = 1,jcap+1
     lwvdef(n)=lwvdef(n)-1
   enddo
#endif
   latmax=latmax*2
!
   deallocate (levdis)
   deallocate (lerdis)
   deallocate (londis)
   deallocate (latdis)
#ifdef DFS
   deallocate (lmdis,lndis)
#else
   deallocate (lntdis)
   deallocate (lwvdis)
   deallocate (lnpdis)           ! SCC 12/20/05
#endif
!
   return
   end subroutine dimset
!-------------------------------------------------------------------------------
