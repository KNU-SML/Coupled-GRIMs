!
   subroutine rmpiset_sub(levr,igrd,jgrd,npes,ncol,nrow,                       &
                         levrp,igrd1p,jgrd1p,iwav1p,lnwavp,llwavp)
!-------------------------------------------------------------------------------
!
!--- This subroutine is based on rmpiset.F90 in library
!
   implicit none
   integer              ::  npes,ncol,nrow,n,nc,igrd1,jgrd1,iwav1,jwav1
!
   integer              ::  levr,igrd,jgrd
   integer              ::  levmax,lonmax,latmax,lntmax,lnpmax,lwvmax
   integer              ::  levrp,igrd1p,jgrd1p,iwav1p,lnwavp,llwavp
!
   real                 ::  pesx
!
   igrd1=igrd+1
   jgrd1=jgrd+1
   iwav1=int((igrd-12)/3)*2+1
   jwav1=int((jgrd-12)/3)*2+1
!
   if( ncol.eq.0 ) then
     pesx=npes
     nc=sqrt(pesx)
     do n = nc,2,-1
       if( mod(npes,n).eq.0 ) then
         ncol=n
         go to 1234
       endif
     enddo
     ncol=1
   else
     if( mod(npes,ncol).ne.0 ) then
       print *,' user provided npes=',npes
       print *,' with an invalid ncol=',ncol
       call abort
       stop
     endif
   endif
1234  nrow=npes/ncol   
!
   call rdimset_sub(iwav1,jwav1,levr,igrd1,jgrd1,                              &
               npes,ncol,nrow,                                                 &
               levmax,                                                         &
               lonmax,                                                         &
               latmax,                                                         &
               lwvmax,                                                         &
               lntmax,                                                         &
               lnpmax)
!
   levrp=levmax
   igrd1p=lonmax
   jgrd1p=latmax
   iwav1p=lwvmax
   lnwavp=lntmax
   llwavp=lnpmax
!
   return
   end subroutine rmpiset_sub
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rdimset_sub(iwav1,jwav1,levr,igrd1,jgrd1,                        &
                     npes,ncol,nrow,                                           &
                     levmax,                                                   &
                     lonmax,                                                   &
                     latmax,                                                   &
                     lwvmax,                                                   &
                     lntmax,                                                   &
                     lnpmax)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    rdimset
!            
! abstract: preset all starting point and length for 
!           all pe for global spectral model.
!
! program history log:
!   1999-06-27  henry juang    finish entire test for gsm
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:   call rdimset(levr,igrd1,jgrd1,
!    *                  npes,ncol,nrow,
!    *                  levstr,levlen,levmax,
!    *                  lonstr,lonlen,lonmax,
!    *                  latstr,latlen,latmax,
!    *                  lwvstr,lwvlen,lwvmax,
!    *                  lntstr,lntlen,lntmax,
!    *                  lnpstr,lnplen,lnpmax)
!
!    input argument lists:
!   iwav1   - integer spectral wavenumber
!   levr   - integer vertical layer number
!   igrd1   - integer gaussian grid for longitude
!   jgrd1   - integer gaussian grid for latitude
!   npes   - integer number of pe used: npes=ncol*nrow
!   ncol   - integer number of column
!   nrow   - integer number of nrow
!
!    output argument list:
!   lev*   - integer (npes) related to layers for each pe
!   lon*   - integer (npes) related to longitude for each pe
!   lat*   - integer (npes) related to latitude for each pe
!   lnt*   - integer (npes) related to npes cut of spectral
!   lnp*   - integer (npes) related to nrow cut of spectral
!   lwv*   - integer (npes) related to group of spectral in l
!   *str   - integer (npes) related to each starting point
!   *len   - integer (npes) related to each length
!   *max   - integer related to maximal length of the kind
! 
! subprograms called:
!   equdiv_sub   - to compute about equal number of subgroup by division
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  levr,igrd1,jgrd1,npes,ncol,nrow,iwav1,jwav1
   integer              ::  latg2,lonf2,nr,nc,nn,jcaprm,nremain
   integer              ::  lnp,n2a,lx,ll,lh,n
   integer              ::  levmax,lonmax,latmax,lwvmax,lntmax,lnpmax
   integer              ::  levpnt,lonpnt,latpnt,lwvpnt,lntpnt,lnppnt
   integer,parameter    ::  npesx=10000
   integer              ::  levstr(0:npesx-1),levlen(0:npesx-1)
   integer              ::  lonstr(0:npesx-1),lonlen(0:npesx-1)
   integer              ::  latstr(0:npesx-1),latlen(0:npesx-1)
   integer              ::  lwvstr(0:npesx-1),lwvlen(0:npesx-1)
   integer              ::  lntstr(0:npesx-1),lntlen(0:npesx-1)
   integer              ::  lnpstr(0:npesx-1),lnplen(0:npesx-1)
   integer              ::  levdis(npesx),londis(npesx),lntdis(npesx)
   integer              ::  lwvdis(npesx),latdis(npesx),lnpdis(npesx)
   integer              ::  lwvdef(npesx)
!
   call equdiv_sub(levr ,ncol,levdis)
   call equdiv_sub(iwav1 ,nrow,lwvdis)
!
   do nr = 1,nrow
     lnpdis(nr)=lwvdis(nr)*jwav1
   enddo                  
!
   latg2=jgrd1/2
   call equdiv_sub(igrd1 ,ncol,londis)
   call equdiv_sub(latg2 ,nrow,latdis)
!
   levmax=0
   lonmax=0
   latmax=0
   lwvmax=0
   lntmax=0
   lnpmax=0
!
   latpnt=1
   lwvpnt=1
   lntpnt=1
   lnppnt=1
   n=0
!
   do nr = 1,nrow
     levpnt=1
     lonpnt=1
     call equdiv_sub(lnpdis(nr),ncol,lntdis)
!
     do nc = 1,ncol
!
       levstr(n)=levpnt
       levlen(n)=levdis(nc)
       levpnt=levpnt+levdis(nc)
       levmax=max(levmax,levlen(n))
!
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
!
       n=n+1
!
     enddo
!
     if( nr.lt.nrow ) then
       latpnt=latpnt+latdis(nr)
       lwvpnt=lwvpnt+lwvdis(nr)
       lnppnt=lnppnt+lnpdis(nr)
     endif
!
   enddo
!
   do n = 0,npes-1
     lwvstr(n)=lwvstr(n)-1
     lnpstr(n)=lnpstr(n)-1
     lntstr(n)=lntstr(n)-1
   enddo
!
   levmax=levmax
   lonmax=lonmax
   latmax=latmax*2
   lwvmax=lwvmax   ! consider as iwv1p
   lntmax=lntmax
   lnpmax=lnpmax
!
   return
   end subroutine rdimset_sub
!-------------------------------------------------------------------------------
