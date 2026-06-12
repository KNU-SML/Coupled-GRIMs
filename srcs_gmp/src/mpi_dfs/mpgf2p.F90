#include "define.h"
   subroutine mpgf2p(a,lonf2,latg2,b,lonf2p,latg2p,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgf2p
!            
! abstract: transpose (if,jf,kf) to (ip,jp,kf)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpgf2p(a,lonf2,latg2,b,lonf2p,latg2p,ntotal)
!
!    input argument lists:
!   a   - real (lonf2,latg2,ntotal) full field 
!   lonf2   - integer total grid in longitude time 2
!   latg2   - integer total grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (lonf2p,latg2p,ntotal) partial field for each pe 
!   lonf2p   - integer partial grid in longitude time 2
!   latg2p   - integer partial grid in latitude divide 2
! 
! subprograms called:
!   mpi_scatterv   - scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use commpi, only : real_type,comm_world,lonlen,lonstr,latlen,latstr,        &
                      npes,mype,master
   use dfsvar, only : latdef,iba,jgs
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonf2,latg2,lonf2p,latg2p,ntotal
   integer            ::  ii,jj,n,mk,i,j,k,ierr,ilen,ist
   real               ::  a(lonf2,latg2,ntotal),b(lonf2p,latg2p,ntotal)
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
   allocate(tmpsnd(lonf2 *latg2 *ntotal))
   allocate(tmprcv(lonf2p*latg2p*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       ilen=lonlen(n)
       ist=lonstr(n)
       do k = 1,ntotal
         do j = 1,latlen(n)
           jj=j+latstr(n)-jgs
           jj=latdef(jj)
           tmpsnd(mk+1:mk+ilen)=a(ist:ist+ilen-1,jj,k)
           mk=mk+ilen
           tmpsnd(mk+1:mk+ilen)=a(iba+ist:iba+ist+ilen-1,jj,k)
           mk=mk+ilen
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
     len(mype)=ntotal*lonlen(mype)*2*latlen(mype)
   endif
!
   call mpi_scatterv(tmpsnd,len(0),loc(0),real_type,                           &
                     tmprcv,len(mype),real_type,0,comm_world,ierr)
!
   mk=0
   ilen=lonlen(mype)
   do k = 1,ntotal
     do j = 1,latlen(mype)
       b(1:ilen,j,k)=tmprcv(mk+1:mk+ilen)
       mk=mk+ilen
       b(ilen+1:2*ilen,j,k)=tmprcv(mk+1:mk+ilen)
       mk=mk+ilen
!
! fill extra space
!
        b(ilen*2+1:lonf2p,j,k)=tmprcv(mk)
     enddo
!
! fill extra space
!
     do j = latlen(mype)+1,latg2p
       b(1:lonf2p,j,k)=tmprcv(mk)
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpgf2p
!-------------------------------------------------------------------------------
